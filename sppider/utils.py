"""
utils.py — utilidades compartidas para sppIDer

Convención de paths
───────────────────
Cada archivo de entrada se indica con su propia ruta independiente.
No existe un directorio de entrada global: el fasta de referencia y
los fastq pueden (y normalmente estarán) en directorios distintos.

  --ref     /lustre/genomes/refs/SaccharomycesCombo.fasta
  --r1      /scratch/project/reads/muestra1_R1.fastq.gz
  --r2      /scratch/project/reads/muestra1_R2.fastq.gz
  --outdir  /scratch/project/results/muestra1
  --sample  muestra1

Todas las rutas se resuelven a absolutas en el momento de parsear los
argumentos, antes de que el pipeline empiece. Si un archivo no existe
se obtiene un error inmediato con la ruta completa que se intentó.
"""

import logging
import subprocess
import sys
import time
from pathlib import Path

# Directorio donde viven los scripts R (junto a este utils.py)
SCRIPT_DIR = Path(__file__).resolve().parent


# ── Paths ─────────────────────────────────────────────────────────────────────

def make_outdir(outdir: str) -> Path:
    """Crea y devuelve --outdir como Path absoluto."""
    p = Path(outdir).resolve()
    p.mkdir(parents=True, exist_ok=True)
    return p


def resolve_input(path_str: str, label: str = "") -> Path:
    """
    Resuelve un archivo de entrada a su Path absoluto.
    Acepta rutas absolutas o relativas desde CWD.
    Lanza FileNotFoundError con mensaje claro si no existe.
    """
    p = Path(path_str).resolve()
    if not p.exists():
        tag = f" ({label})" if label else ""
        raise FileNotFoundError(
            f"Archivo de entrada{tag} no encontrado: {p}"
        )
    return p


# ── Logging ───────────────────────────────────────────────────────────────────

def setup_logger(log_file: Path) -> logging.Logger:
    """Logger con salida a consola (INFO) y a fichero (DEBUG)."""
    logger = logging.getLogger("sppIDer")
    if logger.handlers:
        return logger
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger


# ── Subprocesos ───────────────────────────────────────────────────────────────

def run(cmd: list, logger: logging.Logger,
        stdout_file: Path = None,
        cwd: Path = None) -> None:
    """
    Ejecuta un comando externo.
    stdout_file: si se indica, redirige stdout a ese fichero.
    Lanza RuntimeError si el código de retorno no es 0.
    """
    cmd_str = " ".join(str(c) for c in cmd)
    logger.debug("CMD: %s", cmd_str)
    t0 = time.time()

    kwargs = dict(
        stderr=subprocess.PIPE,
        cwd=str(cwd) if cwd else None,
    )
    if stdout_file:
        with open(stdout_file, "wb") as fout:
            result = subprocess.run([str(c) for c in cmd],
                                    stdout=fout, **kwargs)
    else:
        result = subprocess.run([str(c) for c in cmd],
                                stdout=subprocess.PIPE, **kwargs)

    elapsed = time.time() - t0
    logger.info("OK (%.1f s): %s", elapsed, Path(str(cmd[0])).name)

    if result.returncode != 0:
        err = result.stderr.decode(errors="replace")
        logger.error("STDERR de '%s':\n%s", cmd[0], err)
        raise RuntimeError(
            f"Falló (código {result.returncode}): {cmd_str}"
        )


def run_r(script_name: str, r_args: list, logger: logging.Logger,
          cwd: Path = None) -> None:
    """Ejecuta un script R ubicado junto a utils.py."""
    script = SCRIPT_DIR / script_name
    if not script.exists():
        raise FileNotFoundError(f"Script R no encontrado: {script}")
    run(["Rscript", str(script)] + [str(a) for a in r_args],
        logger=logger, cwd=cwd)


# ── Limpieza ──────────────────────────────────────────────────────────────────

def remove_files(*paths: Path, logger: logging.Logger = None) -> None:
    """Elimina archivos intermedios y registra el espacio liberado."""
    for p in paths:
        p = Path(p)
        if p.exists():
            size_mb = p.stat().st_size / 1_048_576
            p.unlink()
            if logger:
                logger.info("Eliminado: %s (%.1f MB)", p.name, size_mb)


# ── Tiempo ────────────────────────────────────────────────────────────────────

def format_elapsed(seconds: float) -> str:
    s = int(seconds)
    if s < 60:
        return f"{s} s"
    elif s < 3600:
        m, s = divmod(s, 60)
        return f"{m} min {s} s"
    elif s < 86400:
        h, r = divmod(s, 3600)
        m, s = divmod(r, 60)
        return f"{h} h {m} min {s} s"
    else:
        d, r = divmod(s, 86400)
        h, r = divmod(r, 3600)
        m, s = divmod(r, 60)
        return f"{d} d {h} h {m} min {s} s"


# ── Herramientas externas ─────────────────────────────────────────────────────

def require_tool(name: str) -> None:
    """Verifica que una herramienta está en PATH."""
    from shutil import which
    if which(name) is None:
        raise EnvironmentError(
            f"Herramienta no encontrada en PATH: '{name}'. "
            f"Cárgala con 'module load' o instálala."
        )
