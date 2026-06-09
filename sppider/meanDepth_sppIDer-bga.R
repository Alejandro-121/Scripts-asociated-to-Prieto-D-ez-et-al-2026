#!/usr/bin/env Rscript
# meanDepth_sppIDer-bga.R
# Uso: Rscript meanDepth_sppIDer-bga.R WORKDIR PREFIX
#
# Calcula profundidad media por especie, cromosoma y ventana (modo -bga, por grupos).

suppressPackageStartupMessages({
  library(dplyr)
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
if (length(args) < 2) stop("Uso: Rscript meanDepth_sppIDer-bga.R WORKDIR PREFIX")

workingDir <- file.path(args[1], "")
strainName <- args[2]

# ── Archivos ──────────────────────────────────────────────────────────────────
dataFile   <- file.path(workingDir, paste0(strainName, ".bedgraph"))
chrLenFile <- file.path(workingDir, paste0(strainName, "_chrLens.txt"))
spcAvgFile <- file.path(workingDir, paste0(strainName, "_speciesAvgDepth-g.txt"))
chrAvgFile <- file.path(workingDir, paste0(strainName, "_chrAvgDepth-g.txt"))
winFile    <- file.path(workingDir, paste0(strainName, "_winAvgDepth-g.txt"))

for (f in c(dataFile, chrLenFile))
  if (!file.exists(f)) stop("No se encuentra: ", f)

# ── Datos ─────────────────────────────────────────────────────────────────────
strain <- read.table(dataFile, header = FALSE,
                     col.names = c("chrom", "regionStart", "regionEnd", "value"))
strain$regionLen <- strain$regionEnd - strain$regionStart

chrLens   <- read.table(chrLenFile, header = FALSE,
                        col.names = c("chrom", "length"))
genomeEnd <- sum(as.numeric(chrLens$length))
stepSize  <- signif(genomeEnd %/% 10000, digits = 2)
message("Step size: ", stepSize, " bp")

# Añadir cromosomas sin cobertura
dataUniChr <- unique(strain$chrom)
missingChr <- setdiff(chrLens$chrom, dataUniChr)
for (chrName in missingChr) {
  chrLen <- chrLens$length[chrLens$chrom == chrName]
  strain <- rbind(strain,
                  data.frame(chrom = chrName, regionStart = 0,
                             regionEnd = chrLen, value = 0,
                             regionLen = chrLen))
}

strain$wtValue <- strain$value * strain$regionLen
totalMean  <- sum(as.numeric(strain$wtValue)) / genomeEnd
medianVal  <- median(rep(strain$value, strain$regionLen))
maxVal     <- max(strain$value)

# Función log2 segura
safe_log2 <- function(val, ref) {
  l <- log2(val / ref)
  if (is.na(l) || is.infinite(l) || l < 0) 0 else l
}

# Info especie / cromosoma
sp_split       <- strsplit(as.character(strain$chrom), "-")
strain$species <- vapply(sp_split, `[`, character(1), 1)
strain$chr     <- suppressWarnings(as.integer(vapply(sp_split, `[`, character(1), 2)))

# ── Por especie ───────────────────────────────────────────────────────────────
speciesData    <- splitBy("species", strain)
speciesSummary <- data.frame()
spcCumLen      <- 0

spcHeader <- data.frame(Genome_Pos = "wholeGenome", species = "all",
                        genomeLen = genomeEnd, meanValue = totalMean,
                        log2mean  = 1, max = maxVal, median = medianVal)
write.table(format(spcHeader, scientific = FALSE),
            spcAvgFile, row.names = FALSE, sep = "\t", quote = FALSE)

speciesDataOrdered <- vector("list", length(speciesData))

for (i in seq_along(speciesData)) {
  sd      <- speciesData[[i]]
  spcLen  <- sum(as.numeric(sd$regionLen))
  spcMean <- sum(as.numeric(sd$wtValue)) / spcLen
  spcRow  <- data.frame(
    Genome_Pos = spcCumLen,
    species    = names(speciesData)[i],
    genomeLen  = spcLen,
    meanValue  = spcMean,
    log2mean   = safe_log2(spcMean, totalMean),
    max        = max(sd$value, na.rm = TRUE),
    median     = median(rep(sd$value, sd$regionLen))
  )
  write.table(format(spcRow, scientific = FALSE),
              spcAvgFile, col.names = FALSE, row.names = FALSE,
              sep = "\t", quote = FALSE, append = TRUE)
  spcCumLen <- spcCumLen + spcLen
  speciesDataOrdered[[i]] <- sd[order(sd$chr), ]
}

# ── Por cromosoma y ventana ───────────────────────────────────────────────────
chrHeader <- data.frame(Genome_Pos = "wholeGenome", chrom = "all",
                        chrLen = genomeEnd, meanValue = totalMean,
                        log2mean = 1, max = maxVal, median = medianVal)
write.table(format(chrHeader, scientific = FALSE),
            chrAvgFile, row.names = FALSE, sep = "\t", quote = FALSE)

winHeader <- data.frame(Genome_Pos = "wholeGenome", species = "all",
                        chrom = "all", winStart = 0, winEnd = genomeEnd,
                        meanValue = totalMean, log2mean = 1,
                        max = maxVal, median = medianVal)
write.table(format(winHeader, scientific = FALSE),
            winFile, row.names = FALSE, sep = "\t", quote = FALSE)

chrCumLen <- 0

for (i in seq_along(speciesDataOrdered)) {
  chrData <- splitBy("chr", speciesDataOrdered[[i]])

  for (k in seq_along(chrData)) {
    cd     <- chrData[[k]]
    chrLen <- sum(as.numeric(cd$regionLen))
    chrMean <- sum(as.numeric(cd$wtValue)) / chrLen

    chrRow <- data.frame(
      Genome_Pos = chrCumLen,
      chrom      = as.character(cd$chrom[1]),
      chrLen     = chrLen,
      meanValue  = chrMean,
      log2mean   = safe_log2(chrMean, totalMean),
      max        = max(cd$value, na.rm = TRUE),
      median     = median(rep(cd$value, cd$regionLen))
    )
    write.table(format(chrRow, scientific = FALSE),
                chrAvgFile, col.names = FALSE, row.names = FALSE,
                sep = "\t", quote = FALSE, append = TRUE)

    cd$Genome_Pos <- cd$regionStart + chrCumLen
    wBounds <- unique(c(seq(0, cd$regionEnd[nrow(cd)], stepSize),
                        max(cd$regionEnd)))
    cd$bin  <- cut(cd$regionStart, breaks = wBounds, include.lowest = TRUE)
    winData <- splitBy("bin", cd)

    for (win in winData) {
      wLen    <- sum(as.numeric(win$regionLen))
      wMean   <- sum(as.numeric(win$wtValue)) / wLen
      winRow  <- data.frame(
        Genome_Pos = win$Genome_Pos[1],
        species    = names(speciesData)[i],
        chrom      = as.character(win$chrom[1]),
        winStart   = win$regionStart[1],
        winEnd     = win$regionEnd[nrow(win)],
        meanValue  = wMean,
        log2mean   = safe_log2(wMean, totalMean),
        max        = max(win$value, na.rm = TRUE),
        median     = median(rep(win$value, win$regionLen))
      )
      write.table(format(winRow, scientific = FALSE),
                  winFile, col.names = FALSE, row.names = FALSE,
                  sep = "\t", quote = FALSE, append = TRUE)
    }
    chrCumLen <- chrCumLen + chrLen
  }
}

message("meanDepth_sppIDer-bga.R completado.")