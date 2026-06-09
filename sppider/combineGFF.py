#!/usr/bin/env python3
"""
combineGFF.py — Combina GFFs de regiones CDS para mitoSppIDer.

Los GFFs individuales pueden estar en cualquier directorio.

Uso:
    python3 combineGFF.py \
        --key    /lustre/genomes/annotations/mitoGFFKey.txt \
        --out    SaccharomycesMitoCombo.gff \
        --outdir /lustre/genomes/refs/combined

KEY.txt (TSV — rutas absolutas o relativas al directorio del key):
    Saccharomyces_cerevisiae    /lustre/genomes/annotations/S288c_mito.gff
    Saccharomyces_paradoxus     /lustre/genomes/annotations/CBS432_mito.gff
"""

import argparse
from pathlib import Path

from utils import make_outdir, resolve_input, setup_logger


def parse_args():
    p = argparse.ArgumentParser(
        description="Combina GFFs de regiones CDS para mitoSppIDer.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--key",    required=True,
                   help="Ruta al archivo TSV nombre_unico<TAB>ruta_gff")
    p.add_argument("--out",    required=True,
                   help="Nombre del GFF combinado de salida (solo nombre, sin ruta)")
    p.add_argument("--outdir", required=True,
                   help="Directorio de salida para el GFF combinado")
    return p.parse_args()


def main():
    args    = parse_args()
    outdir  = make_outdir(args.outdir)
    key     = resolve_input(args.key, "key GFF")
    key_dir = key.parent

    out_name = Path(args.out).name
    out_gff  = outdir / out_name

    log_file = outdir / "combineGFF.log"
    logger   = setup_logger(log_file)

    logger.info("=" * 62)
    logger.info("combineGFF iniciado")
    logger.info("  key     : %s", key)
    logger.info("  outdir  : %s", outdir)
    logger.info("  out     : %s", out_gff)
    logger.info("=" * 62)

    cds_count = 0

    with open(out_gff, "w") as out:
        out.write("Species\tStart\tEnd\tMidpoint\tName\n")

        for raw_line in key.read_text().splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 2:
                logger.warning("Línea ignorada: %s", line)
                continue

            uni_id, gff_path_str = parts[0].strip(), parts[1].strip()

            gp = Path(gff_path_str)
            if not gp.is_absolute():
                gp = key_dir / gp
            gp = gp.resolve()
            if not gp.exists():
                raise FileNotFoundError(
                    f"GFF no encontrado para '{uni_id}': {gp}"
                )

            n_cds = 0
            with open(gp) as gff_in:
                for raw in gff_in:
                    raw = raw.strip()
                    if not raw or raw.startswith("#"):
                        continue
                    fields = raw.split("\t")
                    if len(fields) < 9 or fields[2] != "CDS":
                        continue

                    start    = int(fields[3])
                    end      = int(fields[4])
                    midpoint = start + (end - start) // 2

                    name = "unknown"
                    for attr in fields[8].split(";"):
                        attr = attr.strip()
                        if attr.upper().startswith(("ID=", "NAME=")):
                            name = attr.split("=", 1)[1].split()[0]
                            break

                    out.write(f"{uni_id}\t{start}\t{end}\t{midpoint}\t{name}\n")
                    n_cds     += 1
                    cds_count += 1

            logger.info("  %-40s %d regiones CDS", uni_id, n_cds)

    logger.info("=" * 62)
    logger.info("GFF combinado: %d regiones CDS → %s", cds_count, out_gff)
    logger.info("=" * 62)


if __name__ == "__main__":
    main()
