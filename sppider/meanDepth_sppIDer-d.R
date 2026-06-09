#!/usr/bin/env Rscript
# meanDepth_sppIDer-d.R
# Uso: Rscript meanDepth_sppIDer-d.R WORKDIR PREFIX
#
# Calcula profundidad media por especie, cromosoma y ventana (modo -d, por pb).

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# ── Cargar splitBy local (reemplazo de doBy) ──────────────────────────────────
local({
  cmd_args <- commandArgs(FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(sub("--file=", "", file_arg[1]))
  } else {
    script_dir <- "."
  }
  source(file.path(script_dir, "splitBy_local.R"), local = globalenv())
})

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Uso: Rscript meanDepth_sppIDer-d.R WORKDIR PREFIX")

workingDir <- file.path(args[1], "")
strainName <- args[2]

# ── Archivos ──────────────────────────────────────────────────────────────────
dataFile   <- file.path(workingDir, paste0(strainName, "-d.bedgraph"))
chrLenFile <- file.path(workingDir, paste0(strainName, "_chrLens.txt"))
spcAvgFile <- file.path(workingDir, paste0(strainName, "_speciesAvgDepth-d.txt"))
chrAvgFile <- file.path(workingDir, paste0(strainName, "_chrAvgDepth-d.txt"))
winFile    <- file.path(workingDir, paste0(strainName, "_winAvgDepth-d.txt"))

for (f in c(dataFile, chrLenFile))
  if (!file.exists(f)) stop("No se encuentra: ", f)

# ── Longitudes de cromosomas ──────────────────────────────────────────────────
chrLens    <- read.table(chrLenFile, header = FALSE,
                         col.names = c("chrom", "length"))
genomeEnd  <- sum(as.numeric(chrLens$length))
stepSize   <- signif(genomeEnd %/% 10000, digits = 2)
message("Step size: ", stepSize, " bp")

# Estadísticas globales (lectura columna ligera)
valCol    <- fread(dataFile, select = 3, col.names = "value")
totalMean <- mean(valCol$value)
maxValue  <- max(valCol$value)
medValue  <- median(valCol$value)
rm(valCol); gc()

# Info de especies en chrLens
sp_split        <- strsplit(as.character(chrLens$chrom), "-")
chrLens$species <- vapply(sp_split, `[`, character(1), 1)
chrLens$chr     <- as.integer(vapply(sp_split, `[`, character(1), 2))
uniSpecies      <- unique(chrLens$species)

# ── Cabeceras de salida ───────────────────────────────────────────────────────
hdr_spc <- data.frame(Genome_Pos = "wholeGenome", species = "all",
                      genomeLen  = genomeEnd,      meanValue = totalMean,
                      log2mean   = 1,              max = maxValue,
                      median     = medValue)
write.table(format(hdr_spc, scientific = FALSE),
            spcAvgFile, row.names = FALSE, sep = "\t", quote = FALSE)

hdr_chr <- data.frame(Genome_Pos = "wholeGenome", chrom = "all",
                      chrLen     = genomeEnd,      meanValue = totalMean,
                      log2mean   = 1,              max = maxValue,
                      median     = medValue)
write.table(format(hdr_chr, scientific = FALSE),
            chrAvgFile, row.names = FALSE, sep = "\t", quote = FALSE)

hdr_win <- data.frame(Genome_Pos = "wholeGenome", species = "all",
                      chrom  = "all", winStart = 0, winEnd = genomeEnd,
                      meanValue = totalMean, log2mean = 1,
                      max = maxValue, median = medValue)
write.table(format(hdr_win, scientific = FALSE),
            winFile, row.names = FALSE, sep = "\t", quote = FALSE)

# ── Función auxiliar: log2 seguro ─────────────────────────────────────────────
safe_log2 <- function(val, ref) {
  l <- log2(val / ref)
  if (is.na(l) || is.infinite(l) || l < 0) 0 else l
}

spcCumLen <- 0
chrCumLen <- 0

for (species in uniSpecies) {
  message("Procesando especie: ", species)
  spcChrLens <- chrLens[chrLens$species == species, ]

  # Leer solo las filas de esta especie (grep vía pipe)
  spData <- tryCatch(
    read.table(
      pipe(paste("awk '{print $1,$2,$3}'", dataFile, "| grep", species)),
      col.names = c("chrom", "chromPos", "value")
    ),
    error = function(e) data.frame(chrom = character(), chromPos = integer(),
                                   value = numeric())
  )

  # Añadir cromosomas sin cobertura
  missingChr <- setdiff(spcChrLens$chrom, unique(spData$chrom))
  if (length(missingChr) > 0) {
    for (chrName in missingChr) {
      chrLen <- spcChrLens$length[spcChrLens$chrom == chrName]
      spData <- rbind(spData,
                      data.frame(chrom    = rep(chrName, chrLen),
                                 chromPos = seq_len(chrLen),
                                 value    = rep(0, chrLen)))
    }
  }

  sp_split2      <- strsplit(as.character(spData$chrom), "-")
  spData$chr     <- as.integer(vapply(sp_split2, `[`, character(1), 2))

  spcLen <- sum(spcChrLens$length)
  spcRow <- data.frame(
    Genome_Pos = spcCumLen, species = species, genomeLen = spcLen,
    meanValue  = mean(spData$value),
    log2mean   = safe_log2(mean(spData$value), totalMean),
    max        = max(spData$value, na.rm = TRUE),
    median     = median(spData$value)
  )
  write.table(format(spcRow, scientific = FALSE),
              spcAvgFile, col.names = FALSE, row.names = FALSE,
              sep = "\t", quote = FALSE, append = TRUE)
  spcCumLen <- spcCumLen + spcLen

  # Por cromosoma y ventana
  spOrdered <- spData[order(spData$chr), ]
  rm(spData); gc()
  chrData <- splitBy("chr", spOrdered)

  for (k in seq_along(chrData)) {
    cd    <- chrData[[k]]
    chrLen <- nrow(cd)

    chrRow <- data.frame(
      Genome_Pos = chrCumLen,
      chrom      = as.character(cd$chrom[1]),
      chrLen     = chrLen,
      meanValue  = mean(cd$value),
      log2mean   = safe_log2(mean(cd$value), totalMean),
      max        = max(cd$value, na.rm = TRUE),
      median     = median(cd$value)
    )
    write.table(format(chrRow, scientific = FALSE),
                chrAvgFile, col.names = FALSE, row.names = FALSE,
                sep = "\t", quote = FALSE, append = TRUE)

    cd$Genome_Pos <- cd$chromPos + chrCumLen
    wBounds <- unique(c(seq(0, cd$chromPos[nrow(cd)], stepSize),
                        max(cd$chromPos)))
    cd$bin  <- cut(cd$chromPos, breaks = wBounds, include.lowest = TRUE)
    winData <- splitBy("bin", cd)

    for (win in winData) {
      winRow <- data.frame(
        Genome_Pos = win$Genome_Pos[1],
        species    = species,
        chrom      = as.character(win$chrom[1]),
        winStart   = win$chromPos[1],
        winEnd     = win$chromPos[nrow(win)],
        meanValue  = mean(win$value),
        log2mean   = safe_log2(mean(win$value), totalMean),
        max        = max(win$value, na.rm = TRUE),
        median     = median(win$value)
      )
      write.table(format(winRow, scientific = FALSE),
                  winFile, col.names = FALSE, row.names = FALSE,
                  sep = "\t", quote = FALSE, append = TRUE)
    }
    chrCumLen <- chrCumLen + chrLen
  }
  rm(spOrdered, chrData); gc()
}

message("meanDepth_sppIDer-d.R completado.")