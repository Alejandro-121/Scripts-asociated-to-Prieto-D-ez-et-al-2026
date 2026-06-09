#!/usr/bin/env python3
"""
parseSamFile.py — Parsea el SAM de BWA y tabula lecturas por especie y MQ.

Uso:
    python3 parseSamFile.py OUTDIR SAMPLE
"""

import sys
import time
from pathlib import Path


def main():
    if len(sys.argv) != 3:
        print("Uso: parseSamFile.py OUTDIR SAMPLE", file=sys.stderr)
        sys.exit(1)

    outdir = Path(sys.argv[1]).resolve()
    sample = sys.argv[2]

    sam_file    = outdir / f"{sample}.sam"
    out_mq      = outdir / f"{sample}_MQ.txt"
    out_chr_len = outdir / f"{sample}_chrLens.txt"

    if not sam_file.exists():
        print(f"ERROR: SAM no encontrado: {sam_file}", file=sys.stderr)
        sys.exit(1)

    t0 = time.time()

    species_dict  = {"*": {i: 0 for i in range(61)}}
    species_order = ["*"]

    with open(out_chr_len, "w") as chr_out, open(sam_file) as sam_in:
        for raw in sam_in:
            if raw.startswith("@SQ"):
                fields  = raw.split("\t")
                sn      = fields[1].split(":")[1]
                ln      = fields[2].split(":")[1].strip()
                chr_out.write(f"{sn}\t{ln}\n")
                sp_name = sn.rsplit("-", 1)[0]
                if sp_name not in species_dict:
                    species_dict[sp_name]  = {i: 0 for i in range(61)}
                    species_order.append(sp_name)
            elif not raw.startswith("@"):
                fields  = raw.split("\t")
                chr_ref = fields[2]
                mq      = min(int(fields[4]), 60)
                sp_name = chr_ref.rsplit("-", 1)[0] if chr_ref != "*" else "*"
                if sp_name not in species_dict:
                    species_dict[sp_name]  = {i: 0 for i in range(61)}
                    species_order.append(sp_name)
                species_dict[sp_name][mq] += 1

    with open(out_mq, "w") as mq_out:
        mq_out.write("Species\tMQscore\tcount\n")
        for sp in species_order:
            for score, count in species_dict[sp].items():
                mq_out.write(f"{sp}\t{score}\t{count}\n")

    print(f"parseSamFile: {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
