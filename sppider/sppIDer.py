#!/usr/bin/env python3
"""
sppIDer.py — Pipeline principal para análisis de genomas híbridos.

El fasta de referencia y los fastq pueden estar en directorios completamente
distintos. Cada archivo se indica con su propia ruta (absoluta o relativa).

Uso:
    python3 sppIDer.py \
        --sample  muestra1 \
        --ref     /lustre/genomes/refs/SaccharomycesCombo.fasta \
        --r1      /scratch/reads/muestra1_R1.fastq.gz \
        --r2      /scratch/reads/muestra1_R2.fastq.gz \
        --outdir  /scratch/results/muestra1 \
        [--byGroup] [--threads N] [--keep-intermediates]
"""

import argparse
import multiprocessing
import time
from pathlib import Path

from utils import (
    SCRIPT_DIR,
    format_elapsed,
    make_outdir,
    remove_files,
    require_tool,
    resolve_input,
    run,
    run_r,
    setup_logger,
)


def parse_args():
    p = argparse.ArgumentParser(
        description="sppIDer: análisis de composición genómica en híbridos.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sample",  required=True,
                   help="Nombre de la muestra (prefijo de los archivos de salida)")
    p.add_argument("--ref",     required=True,
                   help="Ruta al genoma de referencia combinado (.fasta)")
    p.add_argument("--r1",      required=True,
                   help="Ruta a las lecturas R1 (.fastq / .fastq.gz)")
    p.add_argument("--r2",      default=None,
                   help="Ruta a las lecturas R2, opcional (.fastq / .fastq.gz)")
    p.add_argument("--outdir",  required=True,
                   help="Directorio de salida (se crea si no existe)")
    p.add_argument("--byGroup", action="store_true",
                   help="Cobertura por grupos -bga (más rápido en genomas grandes)")
    p.add_argument("--threads", type=int, default=multiprocessing.cpu_count(),
                   help="Hilos para BWA y samtools")
    p.add_argument("--keep-intermediates", action="store_true",
                   help="No eliminar SAM/BAM/bedgraph intermedios")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = make_outdir(args.outdir)
    sample = args.sample

    # ── Resolver archivos de entrada (cada uno su propia ruta) ───────────────
    ref = resolve_input(args.ref, "ref")
    r1  = resolve_input(args.r1,  "r1")
    r2  = resolve_input(args.r2,  "r2") if args.r2 else None

    log_file = outdir / f"{sample}_sppIDer.log"
    logger   = setup_logger(log_file)

    for tool in ["bwa", "samtools", "genomeCoverageBed", "Rscript"]:
        require_tool(tool)

    t0      = time.time()
    threads = str(args.threads)

    logger.info("=" * 62)
    logger.info("sppIDer iniciado")
    logger.info("  sample  : %s", sample)
    logger.info("  ref     : %s", ref)
    logger.info("  r1      : %s", r1)
    logger.info("  r2      : %s", r2 or "—")
    logger.info("  outdir  : %s", outdir)
    logger.info("  modo    : %s", "byGroup (-bga)" if args.byGroup else "byBP (-d)")
    logger.info("  threads : %s", threads)
    logger.info("=" * 62)

    # ── Archivos intermedios (todos en outdir) ────────────────────────────────
    sam      = outdir / f"{sample}.sam"
    bam_view = outdir / f"{sample}.view.bam"
    bam_sort = outdir / f"{sample}.sort.bam"
    bedgraph = outdir / (f"{sample}.bedgraph" if args.byGroup
                         else f"{sample}-d.bedgraph")

    # ── 1. BWA mem ────────────────────────────────────────────────────────────
    logger.info("── 1/7  BWA mem ─────────────────────────────────────────")
    bwa_cmd = ["bwa", "mem", "-t", threads, str(ref), str(r1)]
    if r2:
        bwa_cmd.append(str(r2))
    run(bwa_cmd, logger=logger, stdout_file=sam, cwd=outdir)
    logger.info("Acumulado: %s", format_elapsed(time.time() - t0))

    # ── 2. parseSamFile ───────────────────────────────────────────────────────
    logger.info("── 2/7  parseSamFile ────────────────────────────────────")
    run(["python3", str(SCRIPT_DIR / "parseSamFile.py"), str(outdir), sample],
        logger=logger, cwd=outdir)
    logger.info("Acumulado: %s", format_elapsed(time.time() - t0))

    # ── 3. MQ scores plot ─────────────────────────────────────────────────────
    logger.info("── 3/7  MQ scores plot ──────────────────────────────────")
    run_r("MQscores_sumPlot.R", [str(outdir), sample], logger=logger)
    logger.info("Acumulado: %s", format_elapsed(time.time() - t0))

    # ── 4. samtools view → sort ───────────────────────────────────────────────
    logger.info("── 4/7  samtools view + sort ────────────────────────────")
    run(["samtools", "view", "-@", threads, "-q", "3", "-bhSu",
         str(sam), "-o", str(bam_view)],
        logger=logger, cwd=outdir)
    if not args.keep_intermediates:
        remove_files(sam, logger=logger)

    run(["samtools", "sort", "-@", threads, "-m", "8G", "-T", f"/tmp/{sample}_sort", str(bam_view),
         "-o", str(bam_sort)],
        logger=logger, cwd=outdir)
    if not args.keep_intermediates:
        remove_files(bam_view, logger=logger)
    logger.info("Acumulado: %s", format_elapsed(time.time() - t0))

    # ── 5. bedtools coverage ──────────────────────────────────────────────────
    logger.info("── 5/7  bedtools coverage ───────────────────────────────")
    flag = "-bga" if args.byGroup else "-d"
    run(["genomeCoverageBed", flag, "-ibam", str(bam_sort)],
        logger=logger, stdout_file=bedgraph, cwd=outdir)
    if not args.keep_intermediates:
        remove_files(bam_sort, logger=logger)
    logger.info("Acumulado: %s", format_elapsed(time.time() - t0))

    # ── 6. meanDepth ──────────────────────────────────────────────────────────
    logger.info("── 6/7  meanDepth ───────────────────────────────────────")
    r_script = "meanDepth_sppIDer-bga.R" if args.byGroup else "meanDepth_sppIDer-d.R"
    run_r(r_script, [str(outdir), sample], logger=logger)
    if not args.keep_intermediates:
        remove_files(bedgraph, logger=logger)
    logger.info("Acumulado: %s", format_elapsed(time.time() - t0))

    # ── 7. Plots ──────────────────────────────────────────────────────────────
    logger.info("── 7/7  gráficas ────────────────────────────────────────")
    run_r("sppIDer_depthPlot_forSpc.R", [str(outdir), sample], logger=logger)
    run_r("sppIDer_depthPlot.R",        [str(outdir), sample], logger=logger)

    logger.info("=" * 62)
    logger.info("sppIDer completado en %s", format_elapsed(time.time() - t0))
    logger.info("Resultados en: %s", outdir)
    logger.info("=" * 62)


if __name__ == "__main__":
    main()
