#!/usr/bin/env python3
"""
aggregate_sppIDer_report.py

Agrega _MQsummary.txt + _speciesAvgDepth-d.txt de sppIDer y genera:
  1. species_assignment_summary.tsv
  2. sppIDer_report.html  (pestanas por warning, cobertura media)

Uso:
    python3 aggregate_sppIDer_report.py /path/to/sppider_out
"""

import sys
import re
import csv
from pathlib import Path
from collections import OrderedDict

# ==============================================================================
# UMBRALES — ajusta a tu criterio
# ==============================================================================
MIN_READS        = 500_000
MAX_UNMAPPED_PCT = 20.0
PURE_THRESHOLD   = 70.0
HYBRID_BEST_MIN  = 30.0
HYBRID_SEC_MIN   = 20.0
CONTAM_THRESHOLD = 10.0
LOW_COV          = 5.0

# ==============================================================================
# PARSEO
# ==============================================================================

def _to_float(val, default=None):
    """Convierte a float ignorando 'NA' y valores no numéricos."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def parse_mqsummary(filepath):
    info = {"file": str(filepath), "species_pct": OrderedDict()}
    with open(filepath) as fh:
        lines = fh.readlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if "Num reads =" in line:
            info["total_reads"] = int(line.split("=")[1].strip())
        elif "Num mapped reads =" in line:
            info["mapped_reads"] = int(line.split("=")[1].strip())
        elif "Unmapped reads =" in line:
            m = re.search(r"([\d.]+)\s*%", line)
            info["unmapped_pct"] = float(m.group(1)) if m else 0.0
        elif "Average MQ =" in line:
            info["avg_mq"] = _to_float(line.split("=")[1].strip())
        elif "Median MQ" in line and "=" in line:
            info["median_mq"] = _to_float(line.split("=")[1].strip())
        elif line.startswith("Species"):
            continue
        else:
            parts = line.split()
            if len(parts) >= 2:
                sp_name = parts[0]
                pct_matches = re.findall(r"([\d.]+)\s*%", line)
                pct_all = float(pct_matches[0]) if pct_matches else 0.0
                info["species_pct"][sp_name] = pct_all
    return info


def parse_species_depth(filepath):
    """Parsea _speciesAvgDepth-d.txt o _speciesAvgDepth-g.txt"""
    depths = {}
    if not filepath.exists():
        return depths
    with open(filepath) as f:
        header = True
        for line in f:
            line = line.strip()
            if not line:
                continue
            if header:
                header = False
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            sp = parts[1].strip()
            if sp == "all":
                depths["__global__"] = {
                    "genome_len": int(parts[2]),
                    "mean_depth": float(parts[3]),
                    "max": int(float(parts[5])) if len(parts) > 5 else 0,
                    "median": float(parts[6]) if len(parts) > 6 else 0,
                }
                continue
            try:
                depths[sp] = {
                    "genome_len": int(parts[2]),
                    "mean_depth": float(parts[3]),
                    "log2": float(parts[4]),
                    "max": int(float(parts[5])) if len(parts) > 5 else 0,
                    "median": float(parts[6]) if len(parts) > 6 else 0,
                }
            except (ValueError, IndexError):
                pass
    return depths


def classify_sample(info):
    sp_pct = info.get("species_pct", {})
    depths = info.get("depths", {})
    flags = []
    flag_keys = set()

    sp_mapped = {k: v for k, v in sp_pct.items() if k != "Unmapped"}
    sorted_sp = sorted(sp_mapped.items(), key=lambda x: x[1], reverse=True)
    best_sp, best_pct = sorted_sp[0] if sorted_sp else ("unknown", 0)
    second_sp, second_pct = sorted_sp[1] if len(sorted_sp) > 1 else ("none", 0)

    info["best_species"] = best_sp
    info["best_pct"] = best_pct
    info["second_species"] = second_sp
    info["second_pct"] = second_pct
    info["best_cov"] = depths.get(best_sp, {}).get("mean_depth")
    info["global_cov"] = depths.get("__global__", {}).get("mean_depth")

    if best_pct >= PURE_THRESHOLD:
        info["call"] = f"pure_{best_sp}"
    elif best_pct >= HYBRID_BEST_MIN and second_pct >= HYBRID_SEC_MIN:
        info["call"] = f"hybrid_{best_sp}_x_{second_sp}"
    else:
        info["call"] = "ambiguous"
        flags.append("AMBIGUOUS: asignacion no clara")
        flag_keys.add("ambiguous")

    if info.get("total_reads", 0) < MIN_READS:
        flags.append(f"LOW_READS: {info.get('total_reads',0):,}")
        flag_keys.add("low_reads")

    if info.get("unmapped_pct", 0) > MAX_UNMAPPED_PCT:
        flags.append(f"HIGH_UNMAPPED: {info['unmapped_pct']:.1f}%")
        flag_keys.add("high_unmapped")

    if info["best_cov"] is not None and info["best_cov"] < LOW_COV:
        flags.append(f"LOW_COV: {info['best_cov']:.1f}x on {best_sp}")
        flag_keys.add("low_coverage")

    expected = {best_sp}
    if "hybrid" in info["call"]:
        expected.add(second_sp)
    for sp, pct in sp_mapped.items():
        if sp not in expected and pct >= CONTAM_THRESHOLD:
            flags.append(f"CONTAMINATION: {sp}={pct:.1f}%")
            flag_keys.add("contamination")

    info["flags"] = flags
    info["flag_keys"] = flag_keys
    return info


# ==============================================================================
# TSV
# ==============================================================================

def write_tsv(samples, outpath):
    if not samples:
        return
    sp_names = [s for s in samples[0]["species_pct"].keys()]
    sp_no_unmap = [s for s in sp_names if s != "Unmapped"]
    with open(outpath, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        header = (["Sample", "Total_reads", "Mapped_reads", "Unmapped_pct",
                   "Avg_MQ", "Global_cov"] +
                  [f"{sp}_pct" for sp in sp_names] +
                  [f"{sp}_cov" for sp in sp_no_unmap] +
                  ["Best_species", "Best_pct", "Best_cov",
                   "Second_species", "Second_pct", "Call", "Flags"])
        w.writerow(header)
        for s in samples:
            gc = s.get("global_cov")
            avg_mq = s.get("avg_mq")
            row = [s["sample"], s.get("total_reads",""), s.get("mapped_reads",""),
                   f"{s.get('unmapped_pct',0):.1f}",
                   f"{avg_mq:.1f}" if avg_mq is not None else "NA",
                   f"{gc:.1f}" if gc is not None else "NA"]
            for sp in sp_names:
                row.append(f"{s['species_pct'].get(sp,0):.2f}")
            depths = s.get("depths", {})
            for sp in sp_no_unmap:
                d = depths.get(sp, {}).get("mean_depth")
                row.append(f"{d:.2f}" if d is not None else "NA")
            bc = s.get("best_cov")
            row += [s["best_species"], f"{s['best_pct']:.2f}",
                    f"{bc:.1f}" if bc is not None else "NA",
                    s["second_species"], f"{s['second_pct']:.2f}",
                    s["call"],
                    "; ".join(s["flags"]) if s["flags"] else "OK"]
            w.writerow(row)


# ==============================================================================
# HTML
# ==============================================================================

def _esc(t):
    return t.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")


def generate_html(samples, outpath):
    if not samples:
        return
    all_sp = list(samples[0]["species_pct"].keys())
    sp_names = [s for s in all_sp if s != "Unmapped"]
    total = len(samples)
    n_pure = sum(1 for s in samples if s["call"].startswith("pure_"))
    n_hybrid = sum(1 for s in samples if s["call"].startswith("hybrid_"))
    n_ambiguous = sum(1 for s in samples if s["call"] == "ambiguous")
    n_flagged = sum(1 for s in samples if s["flags"])

    call_counts = {}
    for s in samples:
        call_counts[s["call"]] = call_counts.get(s["call"], 0) + 1

    ft = {"low_reads":0, "high_unmapped":0, "ambiguous":0,
          "contamination":0, "low_coverage":0}
    for s in samples:
        for fk in s["flag_keys"]:
            if fk in ft:
                ft[fk] += 1

    has_cov = any(s.get("best_cov") is not None for s in samples)
    max_bcov = max((s.get("best_cov") or 0) for s in samples)
    max_bcov_slider = int(max_bcov) + 10 if max_bcov > 0 else 100
    low_cov_default = int(LOW_COV)

    # --- Table rows ---
    trows = []
    for s in samples:
        hf = len(s["flags"]) > 0
        rc = "flagged" if hf else "clean"
        fd = " ".join(s["flag_keys"]) if s["flag_keys"] else "ok"
        if hf:
            fhtml = "".join(f"<span class='fb'>{_esc(f)}</span>" for f in s["flags"])
        else:
            fhtml = "<span class='ob'>OK</span>"
        spcells = []
        for sp in sp_names:
            p = s["species_pct"].get(sp, 0)
            cc = ("ph" if p >= PURE_THRESHOLD else
                  ("pm" if p >= HYBRID_SEC_MIN else
                   ("pw" if p >= CONTAM_THRESHOLD else "")))
            spcells.append(f"<td class='{cc}'>{p:.1f}</td>")
        clc = ("cp" if s["call"].startswith("pure_") else
               ("ch" if s["call"].startswith("hybrid_") else "ca"))
        uc = "pw" if s.get("unmapped_pct",0) > MAX_UNMAPPED_PCT else ""

        # Coverage cells
        gc = s.get("global_cov")
        gc_cell = f"<td>{gc:.1f}x</td>" if gc is not None else "<td class='na'>-</td>"
        bc = s.get("best_cov")
        if bc is not None:
            bcls = "pw" if bc < LOW_COV else ("ph" if bc >= 30 else "")
            bc_cell = f"<td class='bcc {bcls}'>{bc:.1f}x</td>"
        else:
            bc_cell = "<td class='bcc na'>-</td>"

        # avg_mq cell
        avg_mq = s.get("avg_mq")
        mq_cell = f"<td>{avg_mq:.1f}</td>" if avg_mq is not None else "<td class='na'>NA</td>"

        bcov_data = f'{bc:.2f}' if bc is not None else '-1'
        trows.append(
            f'<tr class="{rc}" data-f="{fd}" data-c="{s["call"]}" data-bcov="{bcov_data}">'
            f'<td class="sn">{s["sample"]}</td>'
            f'<td>{s.get("total_reads",0):,}</td>'
            f'<td>{s.get("mapped_reads",0):,}</td>'
            f'<td class="{uc}">{s.get("unmapped_pct",0):.1f}%</td>'
            f'{mq_cell}'
            f'{"".join(spcells)}'
            f'<td class="{clc}">{s["call"]}</td>'
            f'{gc_cell}{bc_cell}'
            f'<td class="fc">{fhtml}</td></tr>')
    trows_html = "\n".join(trows)

    # --- Distribution ---
    drows = []
    for call, count in sorted(call_counts.items(), key=lambda x: -x[1]):
        pct = count / total * 100
        drows.append(
            f'<tr><td>{call}</td><td>{count}</td><td>'
            f'<div class="bc"><div class="b" style="width:{pct}%"></div>'
            f'<span class="bl">{pct:.1f}%</span></div></td></tr>')
    drows_html = "\n".join(drows)

    sph = "".join(f"<th class='sc'>{sp}</th>" for sp in sp_names)
    nsp = len(sp_names)

    ci_call = 5 + nsp
    ci_gcov = 6 + nsp
    ci_bcov = 7 + nsp

    cov_tab_btn = (f'<button class="tt wt" onclick="st(\'low_coverage\',this)" '
                   f'data-t="low_coverage">Low coverage <span class="bg">{ft["low_coverage"]}</span></button>') if has_cov else ""

    cov_card = ""
    if has_cov:
        cov_card = (f'<div class="pc wo" onclick="go(\'low_coverage\')">'
                    f'<div class="pt">Baja cobertura (&lt;{LOW_COV}x)</div>'
                    f'<div class="pn">{ft["low_coverage"]}</div>'
                    f'<div class="ph">Click para ver</div></div>')

    cov_th = (f'<th onclick="so({ci_gcov})">Global cov &#8597;</th>'
              f'<th onclick="so({ci_bcov})">Best cov &#8597;</th>') if has_cov else ""

    html = f'''<!DOCTYPE html>
<html lang="es"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>sppIDer Report</title>
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Source+Sans+3:wght@300;400;600;700&display=swap');
:root{{--bg:#0f1117;--sf:#1a1d27;--s2:#232733;--bd:#2e3348;--tx:#e2e4ed;--td:#8b8fa3;--ac:#6c8aff;--ag:rgba(108,138,255,.15);--gn:#34d399;--gd:rgba(52,211,153,.12);--yw:#fbbf24;--yd:rgba(251,191,36,.12);--rd:#f87171;--rdd:rgba(248,113,113,.12);--or:#fb923c;--pr:#a78bfa;--pd:rgba(167,139,250,.12);--cy:#22d3ee}}
*{{margin:0;padding:0;box-sizing:border-box}}
body{{font-family:'Source Sans 3',sans-serif;background:var(--bg);color:var(--tx);line-height:1.6;padding:2rem}}
.ct{{max-width:1800px;margin:0 auto}}
.hd{{border-bottom:1px solid var(--bd);padding-bottom:2rem;margin-bottom:2rem}}
.hd h1{{font-size:1.8rem;font-weight:700;letter-spacing:-.02em;margin-bottom:.25rem}}
.hd p{{color:var(--td);font-size:.9rem}}
.sg{{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:.65rem;margin-bottom:2rem}}
.sk{{background:var(--sf);border:1px solid var(--bd);border-radius:10px;padding:.9rem}}
.sk .lb{{font-size:.68rem;text-transform:uppercase;letter-spacing:.08em;color:var(--td);margin-bottom:.3rem}}
.sk .vl{{font-family:'JetBrains Mono',monospace;font-size:1.35rem;font-weight:600}}
.vg{{color:var(--gn)}}.vy{{color:var(--yw)}}.vr{{color:var(--rd)}}.va{{color:var(--ac)}}.vo{{color:var(--or)}}.vp{{color:var(--pr)}}.vc{{color:var(--cy)}}
.se{{background:var(--sf);border:1px solid var(--bd);border-radius:10px;padding:1.5rem;margin-bottom:1.5rem}}
.se h2{{font-size:1.05rem;font-weight:600;margin-bottom:1rem;display:flex;align-items:center;gap:.5rem}}
.se h2 .ic{{font-size:1.2rem}}.se h2 .ht{{font-size:.75rem;color:var(--td);font-weight:400;margin-left:.5rem}}
.pg{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:.7rem}}
.pc{{background:var(--s2);border-radius:8px;padding:.9rem;border-left:3px solid var(--rd);cursor:pointer;transition:all .15s}}
.pc:hover{{transform:translateY(-1px);border-left-width:4px}}
.pc .pt{{font-size:.68rem;text-transform:uppercase;letter-spacing:.06em;color:var(--rd);margin-bottom:.15rem}}
.pc .pn{{font-family:'JetBrains Mono',monospace;font-size:1.25rem;font-weight:600}}
.pc .ph{{font-size:.7rem;color:var(--td);margin-top:.25rem}}
.pc.wo{{border-left-color:var(--or)}}.pc.wo .pt{{color:var(--or)}}
.pc.wp{{border-left-color:var(--pr)}}.pc.wp .pt{{color:var(--pr)}}
.pc.wy{{border-left-color:var(--yw)}}.pc.wy .pt{{color:var(--yw)}}
.tb{{display:flex;gap:0;margin-bottom:1.1rem;border-bottom:2px solid var(--bd);overflow-x:auto}}
.tt{{background:none;border:none;border-bottom:2px solid transparent;padding:.55rem .9rem;margin-bottom:-2px;color:var(--td);font-family:'Source Sans 3',sans-serif;font-size:.8rem;font-weight:600;cursor:pointer;white-space:nowrap;transition:all .15s;display:flex;align-items:center;gap:.3rem}}
.tt:hover{{color:var(--tx)}}.tt.ac{{color:var(--ac);border-bottom-color:var(--ac)}}
.tt .bg{{font-family:'JetBrains Mono',monospace;font-size:.65rem;font-weight:600;padding:.06rem .35rem;border-radius:8px;background:var(--s2);color:var(--td)}}
.tt.ac .bg{{background:var(--ag);color:var(--ac)}}
.tt.ex{{opacity:.5;text-decoration:line-through;border-bottom-color:var(--rd)}}
.tt.ex .bg{{background:var(--rdd);color:var(--rd)}}
.tt.wt .bg{{background:var(--rdd);color:var(--rd)}}.tt.wt.ac{{color:var(--rd);border-bottom-color:var(--rd)}}.tt.wt.ac .bg{{background:var(--rdd);color:var(--rd)}}
.tt.gt .bg{{background:var(--gd);color:var(--gn)}}.tt.gt.ac{{color:var(--gn);border-bottom-color:var(--gn)}}.tt.gt.ac .bg{{background:var(--gd);color:var(--gn)}}
.tt.ht .bg{{background:var(--yd);color:var(--yw)}}.tt.ht.ac{{color:var(--yw);border-bottom-color:var(--yw)}}.tt.ht.ac .bg{{background:var(--yd);color:var(--yw)}}
.tt.at .bg{{background:var(--pd);color:var(--pr)}}.tt.at.ac{{color:var(--pr);border-bottom-color:var(--pr)}}.tt.at.ac .bg{{background:var(--pd);color:var(--pr)}}
.ts{{width:1px;background:var(--bd);margin:.25rem .4rem;flex-shrink:0}}
.tl{{font-size:.72rem;color:var(--td);margin-bottom:.6rem;font-style:italic}}
.sr{{display:flex;gap:.75rem;margin-bottom:1rem;align-items:center}}
.sb{{background:var(--s2);border:1px solid var(--bd);border-radius:6px;padding:.45rem .7rem;color:var(--tx);font-family:'JetBrains Mono',monospace;font-size:.83rem;width:260px;outline:none;transition:border-color .2s}}
.sb:focus{{border-color:var(--ac)}}
.rc{{font-family:'JetBrains Mono',monospace;font-size:.78rem;color:var(--td)}}
.eb{{background:var(--s2);border:1px solid var(--bd);border-radius:6px;padding:.4rem .8rem;color:var(--td);font-size:.8rem;cursor:pointer;transition:all .15s;margin-left:auto;font-family:'Source Sans 3',sans-serif}}
.eb:hover{{border-color:var(--ac);color:var(--tx)}}
.tw{{overflow-x:auto}}table{{width:100%;border-collapse:collapse;font-size:.8rem}}
th{{background:var(--s2);padding:.55rem .45rem;text-align:left;font-weight:600;font-size:.7rem;text-transform:uppercase;letter-spacing:.06em;color:var(--td);border-bottom:2px solid var(--bd);position:sticky;top:0;cursor:pointer;white-space:nowrap;user-select:none;z-index:2}}
th:hover{{color:var(--ac)}}th.sc{{text-align:center}}
td{{padding:.45rem;border-bottom:1px solid var(--bd);white-space:nowrap;font-family:'JetBrains Mono',monospace;font-size:.78rem}}
tr:hover{{background:var(--ag)}}tr.flagged{{background:var(--rdd)}}tr.flagged:hover{{background:rgba(248,113,113,.18)}}
.sn{{font-weight:600}}.ph{{color:var(--gn);font-weight:600}}.pm{{color:var(--yw);font-weight:600}}.pw{{color:var(--or)}}
.cp{{color:var(--gn);font-weight:600}}.ch{{color:var(--yw);font-weight:600}}.ca{{color:var(--rd);font-weight:600}}
.na{{color:var(--td);font-style:italic}}
.fc{{white-space:normal;max-width:360px}}
.fb{{display:inline-block;background:var(--rdd);color:var(--rd);border:1px solid rgba(248,113,113,.3);border-radius:4px;padding:.08rem .32rem;font-size:.65rem;margin:.06rem;font-family:'JetBrains Mono',monospace}}
.ob{{display:inline-block;background:var(--gd);color:var(--gn);border:1px solid rgba(52,211,153,.3);border-radius:4px;padding:.08rem .32rem;font-size:.65rem;font-family:'JetBrains Mono',monospace}}
.bc{{display:flex;align-items:center;gap:.5rem}}.b{{height:16px;background:var(--ac);border-radius:3px;min-width:2px;opacity:.7}}
.bl{{font-family:'JetBrains Mono',monospace;font-size:.75rem;color:var(--td)}}
.dt{{width:auto;min-width:450px}}.dt td{{padding:.35rem .7rem}}
tr.hd{{display:none}}
.cv{{display:flex;align-items:center;gap:.8rem;padding:.75rem 1rem;background:var(--s2);border:1px solid var(--bd);border-radius:8px;margin-bottom:1rem;flex-wrap:wrap}}
.cv label{{font-size:.8rem;font-weight:600;white-space:nowrap}}
.cs{{-webkit-appearance:none;appearance:none;flex:1;max-width:320px;height:6px;border-radius:3px;background:var(--bd);outline:none;cursor:pointer}}
.cs::-webkit-slider-thumb{{-webkit-appearance:none;width:16px;height:16px;border-radius:50%;background:var(--or);cursor:pointer;border:2px solid var(--bg)}}
.cs::-moz-range-thumb{{width:16px;height:16px;border-radius:50%;background:var(--or);cursor:pointer;border:2px solid var(--bg)}}
.cvl{{font-family:'JetBrains Mono',monospace;font-size:.85rem;font-weight:600;color:var(--or);min-width:45px}}
.cvn{{font-family:'JetBrains Mono',monospace;font-size:.78rem;color:var(--td)}}
tr.clw td.bcc{{background:rgba(251,146,60,.18)!important;color:var(--or)!important;font-weight:600}}
@media print{{body{{background:#fff;color:#111;padding:1rem}}.tb,.sr,.eb{{display:none}}.se{{border:1px solid #ccc;break-inside:avoid}}tr.flagged{{background:#fff0f0}}tr.hd{{display:table-row!important}}}}
</style></head><body><div class="ct">

<div class="hd"><h1>sppIDer &mdash; Species Assignment Report</h1>
<p>{total} muestras &middot; pura &ge;{PURE_THRESHOLD}% | hibrido &ge;{HYBRID_BEST_MIN}%+{HYBRID_SEC_MIN}% | contam. &ge;{CONTAM_THRESHOLD}% | unmapped &gt;{MAX_UNMAPPED_PCT}% | min reads {MIN_READS:,} | low cov &lt;{LOW_COV}x</p></div>

<div class="sg">
<div class="sk"><div class="lb">Total</div><div class="vl va">{total}</div></div>
<div class="sk"><div class="lb">Puras</div><div class="vl vg">{n_pure}</div></div>
<div class="sk"><div class="lb">Hibridos</div><div class="vl vy">{n_hybrid}</div></div>
<div class="sk"><div class="lb">Ambiguas</div><div class="vl vp">{n_ambiguous}</div></div>
<div class="sk"><div class="lb">Con alertas</div><div class="vl vr">{n_flagged}</div></div>
<div class="sk"><div class="lb">Low reads</div><div class="vl vo">{ft["low_reads"]}</div></div>
<div class="sk"><div class="lb">High unmap</div><div class="vl vr">{ft["high_unmapped"]}</div></div>
<div class="sk"><div class="lb">Contam.</div><div class="vl vc">{ft["contamination"]}</div></div>
<div class="sk"><div class="lb">Low cov</div><div class="vl vo">{ft["low_coverage"]}</div></div>
</div>

<div class="se"><h2><span class="ic">&#9888;</span> Alertas<span class="ht">(click filtra tabla)</span></h2>
<div class="pg">
<div class="pc wo" onclick="go('low_reads')"><div class="pt">Low reads (&lt;{MIN_READS:,})</div><div class="pn">{ft["low_reads"]}</div><div class="ph">Click para ver</div></div>
<div class="pc" onclick="go('high_unmapped')"><div class="pt">High unmapped (&gt;{MAX_UNMAPPED_PCT}%)</div><div class="pn">{ft["high_unmapped"]}</div><div class="ph">Click para ver</div></div>
<div class="pc wp" onclick="go('ambiguous')"><div class="pt">Ambigua</div><div class="pn">{ft["ambiguous"]}</div><div class="ph">Click para ver</div></div>
<div class="pc wy" onclick="go('contamination')"><div class="pt">Contaminacion (&ge;{CONTAM_THRESHOLD}%)</div><div class="pn">{ft["contamination"]}</div><div class="ph">Click para ver</div></div>
{cov_card}
</div></div>

<div class="se"><h2><span class="ic">&#128202;</span> Distribucion</h2>
<div class="tw"><table class="dt"><tr><th>Asignacion</th><th>N</th><th style="min-width:280px">Proporcion</th></tr>
{drows_html}</table></div></div>

<div class="se"><h2><span class="ic">&#128300;</span> Detalle por muestra</h2>
<div class="tb" id="tB">
<button class="tt ac" onclick="st('all',this)" data-t="all">Todas <span class="bg">{total}</span></button>
<div class="ts"></div>
<button class="tt gt" onclick="st('pure',this)" data-t="pure">Puras <span class="bg">{n_pure}</span></button>
<button class="tt ht" onclick="st('hybrid',this)" data-t="hybrid">Hibridos <span class="bg">{n_hybrid}</span></button>
<button class="tt at" onclick="st('ambiguous_call',this)" data-t="ambiguous_call">Ambiguas <span class="bg">{n_ambiguous}</span></button>
<div class="ts"></div>
<button class="tt wt" onclick="st('flagged',this)" data-t="flagged">Alertas <span class="bg">{n_flagged}</span></button>
<button class="tt wt" onclick="st('low_reads',this)" data-t="low_reads">Low reads <span class="bg">{ft["low_reads"]}</span></button>
<button class="tt wt" onclick="st('high_unmapped',this)" data-t="high_unmapped">High unmap <span class="bg">{ft["high_unmapped"]}</span></button>
<button class="tt wt" onclick="st('contamination',this)" data-t="contamination">Contam. <span class="bg">{ft["contamination"]}</span></button>
<button class="tt wt" onclick="st('ambiguous',this)" data-t="ambiguous">Ambigua (w) <span class="bg">{ft["ambiguous"]}</span></button>
{cov_tab_btn}
</div>
<div class="tl">1 click = incluir | 2 clicks = excluir | 3 clicks = reset</div>
<div class="cv" id="cvC" style="display:none">
<label>Umbral cobertura min:</label>
<input type="range" class="cs" id="cvS" min="0" max="{max_bcov_slider}" step="1" value="0" oninput="uc()">
<span class="cvl" id="cvV">0x</span>
<span class="cvn" id="cvN"></span>
</div>
<div class="sr">
<input type="text" class="sb" id="sB" placeholder="Buscar muestra..." oninput="af()">
<span class="rc" id="rC"></span>
<button class="eb" onclick="ev()">Exportar TSV</button>
</div>
<div class="tw"><table id="mT"><thead><tr>
<th onclick="so(0)">Sample &#8597;</th>
<th onclick="so(1)">Reads &#8597;</th>
<th onclick="so(2)">Mapped &#8597;</th>
<th onclick="so(3)">Unmap% &#8597;</th>
<th onclick="so(4)">AvgMQ &#8597;</th>
{sph}
<th onclick="so({ci_call})">Asignacion &#8597;</th>
{cov_th}
<th>Alertas</th>
</tr></thead><tbody id="tBody">
{trows_html}
</tbody></table></div></div>

</div>
<script>
var inc=new Set(),exc=new Set();
function st(t,b){{if(t==='all'){{inc.clear();exc.clear();document.querySelectorAll('.tt').forEach(x=>x.classList.remove('ac','ex'));b.classList.add('ac')}}else{{document.querySelector('.tt[data-t="all"]').classList.remove('ac');if(inc.has(t)){{inc.delete(t);b.classList.remove('ac');exc.add(t);b.classList.add('ex')}}else if(exc.has(t)){{exc.delete(t);b.classList.remove('ex')}}else{{inc.add(t);b.classList.add('ac');b.classList.remove('ex')}}}}af()}}
function go(t){{inc.clear();exc.clear();document.querySelectorAll('.tt').forEach(x=>x.classList.remove('ac','ex'));inc.add(t);var btn=document.querySelector('.tt[data-t="'+t+'"]');if(btn)btn.classList.add('ac');document.getElementById('sB').value='';af();document.getElementById('tB').scrollIntoView({{behavior:'smooth',block:'start'}})}}
function mt(r,t){{var f=r.dataset.f||'';var c=r.dataset.c||'';switch(t){{case'pure':return c.startsWith('pure_');case'hybrid':return c.startsWith('hybrid_');case'ambiguous_call':return c==='ambiguous';case'flagged':return r.classList.contains('flagged');case'low_reads':return f.includes('low_reads');case'high_unmapped':return f.includes('high_unmapped');case'contamination':return f.includes('contamination');case'ambiguous':return f.includes('ambiguous');case'low_coverage':return f.includes('low_coverage');default:return false}}}}
function af(){{var s=document.getElementById('sB').value.toLowerCase();var R=document.querySelectorAll('#tBody tr');var v=0;var noFilters=inc.size===0&&exc.size===0;var sl=document.getElementById('cvS');var covThr=sl?parseFloat(sl.value):-1;var below=0;R.forEach(function(r){{var n=r.cells[0].textContent.toLowerCase();var sh=n.includes(s);if(sh&&!noFilters){{if(inc.size>0){{var any=false;inc.forEach(function(t){{if(mt(r,t))any=true}});sh=any}}if(sh&&exc.size>0){{exc.forEach(function(t){{if(mt(r,t))sh=false}})}}}}var bv=parseFloat(r.dataset.bcov);var covOk=bv>=0;var isLow=covOk&&bv<covThr;r.classList.toggle('clw',isLow);if(sh&&covThr>0&&covOk&&isLow)sh=false;r.classList.toggle('hd',!sh);if(sh)v++;if(isLow)below++}});document.getElementById('rC').textContent=v+' / '+R.length;if(sl){{document.getElementById('cvV').textContent=covThr+'x';document.getElementById('cvN').textContent=below+' muestras < '+covThr+'x'}}}}
let sC=-1,sA=true;
function so(c){{const t=document.getElementById('tBody');const R=Array.from(t.rows);if(sC===c)sA=!sA;else{{sC=c;sA=true}};R.sort((a,b)=>{{let x=a.cells[c].textContent.replace(/[,%x]/g,'').trim();let y=b.cells[c].textContent.replace(/[,%x]/g,'').trim();let na=parseFloat(x),nb=parseFloat(y);if(!isNaN(na)&&!isNaN(nb))return sA?na-nb:nb-na;return sA?x.localeCompare(y):y.localeCompare(x)}});R.forEach(r=>t.appendChild(r))}}
function ev(){{var TAB=String.fromCharCode(9),NL=String.fromCharCode(10);var R=document.querySelectorAll('#tBody tr:not(.hd)');var H=Array.from(document.querySelectorAll('#mT thead th')).map(function(h){{return h.textContent.replace(/\u2195/g,'').trim()}});var t=H.join(TAB)+NL;R.forEach(function(r){{t+=Array.from(r.cells).map(function(c){{return c.textContent.trim()}}).join(TAB)+NL}});var bl=new Blob([t],{{type:'text/tab-separated-values'}});var u=URL.createObjectURL(bl);var a=document.createElement('a');var nm='sppIDer';if(inc.size>0)nm+='_inc-'+Array.from(inc).join('-');if(exc.size>0)nm+='_exc-'+Array.from(exc).join('-');var sl=document.getElementById('cvS');if(sl&&parseFloat(sl.value)>0)nm+='_mincov'+sl.value+'x';a.href=u;a.download=nm+'.tsv';a.click();URL.revokeObjectURL(u)}}
function uc(){{af()}}
(function(){{var cc=document.getElementById('cvC');if(cc){{var rows=document.querySelectorAll('#tBody tr[data-bcov]');var any=false;rows.forEach(function(r){{if(parseFloat(r.dataset.bcov)>=0)any=true}});if(any)cc.style.display='flex'}}}})();
af();
</script></body></html>'''

    with open(outpath, "w") as fh:
        fh.write(html)


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    if len(sys.argv) < 2:
        print("Uso: python3 aggregate_sppIDer_report.py /path/to/sppider_out",
              file=sys.stderr)
        sys.exit(1)

    base_dir = Path(sys.argv[1]).resolve()
    if not base_dir.is_dir():
        print(f"ERROR: {base_dir} no es un directorio", file=sys.stderr)
        sys.exit(1)

    mq_files = sorted(base_dir.glob("*/*_MQsummary.txt"))
    if not mq_files:
        print(f"ERROR: No se encontraron *_MQsummary.txt en {base_dir}/*/",
              file=sys.stderr)
        sys.exit(1)

    print(f"Encontrados {len(mq_files)} archivos MQsummary")

    samples = []
    for f in mq_files:
        sample_name = f.parent.name
        try:
            info = parse_mqsummary(f)
            info["sample"] = sample_name

            depth_d = f.parent / f"{sample_name}_speciesAvgDepth-d.txt"
            depth_g = f.parent / f"{sample_name}_speciesAvgDepth-g.txt"
            if depth_d.exists():
                info["depths"] = parse_species_depth(depth_d)
            elif depth_g.exists():
                info["depths"] = parse_species_depth(depth_g)
            else:
                info["depths"] = {}

            info = classify_sample(info)
            samples.append(info)
        except Exception as e:
            print(f"WARN: Error parseando {f}: {e}", file=sys.stderr)

    samples.sort(key=lambda x: x["sample"])

    tsv_path = base_dir / "species_assignment_summary.tsv"
    write_tsv(samples, tsv_path)
    print(f"TSV: {tsv_path}")

    html_path = base_dir / "sppIDer_report.html"
    generate_html(samples, html_path)
    print(f"HTML: {html_path}")

    n_flagged = sum(1 for s in samples if s["flags"])
    print(f"\nTotal: {len(samples)} muestras | {n_flagged} con alertas")
    for s in samples:
        if s["flags"]:
            print(f"  ! {s['sample']}: {'; '.join(s['flags'])}")


if __name__ == "__main__":
    main()
