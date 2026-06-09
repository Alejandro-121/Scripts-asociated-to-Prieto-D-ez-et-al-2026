#!/usr/bin/env python3
"""
mitoSppIDer.py — Pipeline para genomas mitocondriales.

Uso:
    python3 mitoSppIDer.py \
        --sample  muestra1_mito \
        --ref     /lustre/genomes/refs/SaccharomycesMitoCombo.fasta \
        --r1      /scratch/reads/muestra1_R1.fastq.gz \
        --r2      /scratch/reads/muestra1_R2.fastq.gz \
        --outdir  /scratch/results/mito/muestra1 \
        [--gff    /lustre/genomes/refs/SaccharomycesMitoCombo.gff] \
        [--threads N] [--keep-intermediates]
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
        description="mitoSppIDer: análisis de genomas mitocondriales híbridos.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sample",  required=True,
                   help="Nombre de la muestra")
    p.add_argument("--ref",     required=True,
                   help="Ruta al genoma mitocondrial combinado (.fasta)")
    p.add_argument("--r1",      required=True,
                   help="Ruta a las lecturas R1 (.fastq / .fastq.gz)")
    p.add_argument("--r2",      default=None,
                   help="Ruta a las lecturas R2, opcional")
    p.add_argument("--gff",     default=None,
                   help="Ruta al GFF combinado de regiones codificantes, opcional")
    p.add_argument("--outdir",  required=True,
                   help="Directorio de salida (se crea si no existe)")
    p.add_argument("--threads", type=int, default=multiprocessing.cpu_count(),
                   help="Hilos para BWA y samtools")
    p.add_argument("--keep-intermediates", action="store_true",
                   help="No eliminar SAM/BAM/bedgraph intermedios")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = make_outdir(args.outdir)
    sample = args.sample

    ref = resolve_input(args.ref, "ref")
    r1  = resolve_input(args.r1,  "r1")
    r2  = resolve_input(args.r2,  "r2")  if args.r2  else None
    gff = resolve_input(args.gff, "gff") if args.gff else None

    log_file = outdir / f"{sample}_mitoSppIDer.log"
    logger   = setup_logger(log_file)

    for tool in ["bwa", "samtools", "genomeCoverageBed", "Rscript"]:
        require_tool(tool)

    t0      = time.time()
    threads = str(args.threads)

    logger.info("=" * 62)
    logger.info("mitoSppIDer iniciado")
    logger.info("  sample  : %s", sample)
    logger.info("  ref     : %s", ref)
    logger.info("  r1      : %s", r1)
    logger.info("  r2      : %s", r2  or "—")
    logger.info("  gff     : %s", gff or "—")
    logger.info("  outdir  : %s", outdir)
    logger.info("  threads : %s", threads)
    logger.info("=" * 62)

    sam      = outdir / f"{sample}.sam"
    bam_view = outdir / f"{sample}.view.bam"
    bam_sort = outdir / f"{sample}.sort.bam"
    bedgraph = outdir / f"{sample}-d.bedgraph"

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
    run(["samtools", "view", "-@", threads, "-q", "1", "-bhSu",
         str(sam), "-o", str(bam_view)],
        logger=logger, cwd=outdir)
    if not args.keep_intermediates:
        remove_files(sam, logger=logger)

    run(["samtools", "sort", "-@", threads, str(bam_view),
         "-o", str(bam_sort)],
        logger=logger, cwd=outdir)
    if not args.keep_intermediates:
        remove_files(bam_view, logger=logger)
    logger.info("Acumulado: %s", format_elapsed(time.time() - t0))

    # ── 5. bedtools coverage (-d) ─────────────────────────────────────────────
    logger.info("── 5/7  bedtools coverage (-d) ──────────────────────────")
    run(["genomeCoverageBed", "-d", "-ibam", str(bam_sort)],
        logger=logger, stdout_file=bedgraph, cwd=outdir)
    if not args.keep_intermediates:
        remove_files(bam_sort, logger=logger)
    logger.info("Acumulado: %s", format_elapsed(time.time() - t0))

    # ── 6. meanDepth ──────────────────────────────────────────────────────────
    logger.info("── 6/7  meanDepth ───────────────────────────────────────")
    run_r("meanDepth_sppIDer-d.R", [str(outdir), sample], logger=logger)
    if not args.keep_intermediates:
        remove_files(bedgraph, logger=logger)
    logger.info("Acumulado: %s", format_elapsed(time.time() - t0))

    # ── 7. Plot ───────────────────────────────────────────────────────────────
    logger.info("── 7/7  gráficas ────────────────────────────────────────")
    r_args = [str(outdir), sample]
    if gff:
        r_args.append(str(gff))
    run_r("mitoSppIDer_depthPlot-d.R", r_args, logger=logger)

    logger.info("=" * 62)
    logger.info("mitoSppIDer completado en %s", format_elapsed(time.time() - t0))
    logger.info("Resultados en: %s", outdir)
    logger.info("=" * 62)


if __name__ == "__main__":
    main()
