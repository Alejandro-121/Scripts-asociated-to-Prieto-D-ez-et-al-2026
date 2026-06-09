#!/usr/bin/env python3
"""
combineRefGenomes.py — Construye el genoma de referencia combinado.

Los fastas individuales y el key pueden estar en cualquier directorio.
La salida (fasta combinado + índices) va a --outdir.

Uso:
    python3 combineRefGenomes.py \
        --key    /lustre/genomes/refs/SaccharomycesGenomesKey.txt \
        --out    SaccharomycesCombo.fasta \
        --outdir /lustre/genomes/refs/combined \
        [--trim  1000]

KEY.txt (TSV — sin guiones en el nombre único, rutas absolutas o relativas
         al directorio donde está el key):
    Saccharomyces_cerevisiae    /lustre/genomes/S288c.fasta
    Saccharomyces_paradoxus     /lustre/genomes/GCA_002079055.1.fasta

Si las rutas del key son relativas, se resuelven desde el directorio
donde está el propio key file.
"""

import argparse
from pathlib import Path

from Bio import SeqIO

from utils import make_outdir, require_tool, resolve_input, run, setup_logger


def parse_args():
    p = argparse.ArgumentParser(
        description="Combina genomas de referencia en uno solo.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--key",    required=True,
                   help="Ruta al archivo TSV nombre_unico<TAB>ruta_fasta")
    p.add_argument("--out",    required=True,
                   help="Nombre del fasta combinado de salida (solo nombre, sin ruta)")
    p.add_argument("--outdir", required=True,
                   help="Directorio donde se escriben el fasta combinado y sus índices")
    p.add_argument("--trim",   type=int, default=0,
                   help="Excluir contigs más cortos que este valor en pb")
    return p.parse_args()


def format_bp(n: int) -> str:
    if n >= 1_000_000_000:
        gb, r = divmod(n, 1_000_000_000)
        mb, r = divmod(r, 1_000_000)
        kb, bp = divmod(r, 1_000)
        return f"{gb} Gb {mb} Mb {kb} Kb {bp} bp"
    elif n >= 1_000_000:
        mb, r = divmod(n, 1_000_000)
        kb, bp = divmod(r, 1_000)
        return f"{mb} Mb {kb} Kb {bp} bp"
    elif n >= 1_000:
        kb, bp = divmod(n, 1_000)
        return f"{kb} Kb {bp} bp"
    return f"{n} bp"


def main():
    args   = parse_args()
    outdir = make_outdir(args.outdir)
    key    = resolve_input(args.key, "key")
    key_dir = key.parent          # directorio del key, base para rutas relativas en él

    out_name  = Path(args.out).name
    out_fasta = outdir / out_name
    len_file  = outdir / f"comboLength_{out_name}.txt"

    log_file = outdir / "combineRefGenomes.log"
    logger   = setup_logger(log_file)

    for tool in ["bwa", "samtools"]:
        require_tool(tool)

    logger.info("=" * 62)
    logger.info("combineRefGenomes iniciado")
    logger.info("  key     : %s", key)
    logger.info("  outdir  : %s", outdir)
    logger.info("  out     : %s", out_fasta)
    logger.info("  trim    : %d bp", args.trim)
    logger.info("=" * 62)

    combo_total = 0

    with open(out_fasta, "w") as out_fa, open(len_file, "w") as lf:
        lf.write(f"{out_name}\tContigs con longitud >= {args.trim} bp\n")

        for raw_line in key.read_text().splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 2:
                logger.warning("Línea ignorada (formato incorrecto): %s", line)
                continue

            uni_id, genome_path_str = parts[0].strip(), parts[1].strip()

            if "-" in uni_id:
                raise ValueError(
                    f"El nombre único '{uni_id}' contiene guiones (-).\n"
                    f"Usa guiones bajos u otro separador."
                )

            # Ruta absoluta → se usa directamente
            # Ruta relativa → se resuelve desde el directorio del key
            gp = Path(genome_path_str)
            if not gp.is_absolute():
                gp = key_dir / gp
            gp = gp.resolve()
            if not gp.exists():
                raise FileNotFoundError(
                    f"Fasta no encontrado para '{uni_id}': {gp}"
                )

            sum_len = 0
            counter = 0
            with open(gp) as fasta_in:
                for rec in SeqIO.parse(fasta_in, "fasta"):
                    if len(rec.seq) < args.trim:
                        continue
                    counter += 1
                    header = f"{uni_id}-{counter}"
                    out_fa.write(f">{header}\n{rec.seq}\n")
                    lf.write(f"{header}\t{len(rec.seq)}\n")
                    sum_len += len(rec.seq)

            lf.write(f"{uni_id}-totalGenome\t{format_bp(sum_len)}\n")
            combo_total += sum_len
            logger.info("  %-40s %d contigs  %s", uni_id, counter, format_bp(sum_len))

        lf.write(f"Combo-totalGenome\t{format_bp(combo_total)}\n")

    logger.info("Total combinado: %s", format_bp(combo_total))

    logger.info("Indexando con BWA...")
    run(["bwa", "index", str(out_fasta)], logger=logger, cwd=outdir)

    logger.info("Indexando con samtools faidx...")
    run(["samtools", "faidx", str(out_fasta)], logger=logger, cwd=outdir)

    logger.info("=" * 62)
    logger.info("Genoma combinado listo: %s", out_fasta)
    logger.info("=" * 62)


if __name__ == "__main__":
    main()
