#!/usr/bin/env Rscript
# sppIDer_depthPlot.R
# Uso: Rscript sppIDer_depthPlot.R WORKDIR PREFIX

suppressPackageStartupMessages({
  library(ggplot2)
  library(modes)
})

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Uso: Rscript sppIDer_depthPlot.R WORKDIR PREFIX")

workingDir <- file.path(args[1], "")
prefix     <- args[2]

# ── Archivos ──────────────────────────────────────────────────────────────────
dFile <- file.path(workingDir, paste0(prefix, "_winAvgDepth-d.txt"))
gFile <- file.path(workingDir, paste0(prefix, "_winAvgDepth-g.txt"))
dataFile <- if (file.exists(dFile)) dFile else if (file.exists(gFile)) gFile else
  stop("No se encuentra winAvgDepth para: ", prefix)

plotCovDist <- file.path(workingDir, paste0(prefix, "_covDistPlots.pdf"))
plotCovBins <- file.path(workingDir, paste0(prefix, "_sppIDerDepthPlot-covBins.pdf"))
binSumFile  <- file.path(workingDir, paste0(prefix, "_covBinsSummary.txt"))
binWinFile  <- file.path(workingDir, paste0(prefix, "_winAvgDepth_wCovBins.txt"))
globBinFile <- file.path(workingDir, paste0(prefix, "_globalBinsIncluded.txt"))

bedData <- read.table(dataFile, header = FALSE, skip = 2,
                      col.names = c("Genome_Pos", "species", "chrom",
                                    "start", "end", "meanValue",
                                    "log2rawCov", "max", "median"))

checkMean  <- mean(bedData$meanValue)
bedData$log2 <- pmax(0, log2(bedData$meanValue / checkMean))
logQuant     <- quantile(bedData$log2, 0.99, na.rm = TRUE)

pdf(NULL)
uniSpecies <- factor(unique(bedData$species), levels = unique(bedData$species))
colorList  <- if (length(uniSpecies) > 1) rainbow(length(uniSpecies)) else "red"
dev.off()

colors   <- setNames(colorList, levels(uniSpecies))
fillScale  <- scale_fill_manual(name  = "Species", values = colors,
                                breaks = levels(uniSpecies))
colorScale <- scale_color_manual(name = "Species", values = colors,
                                 breaks = levels(uniSpecies))
thm <- theme_bw() + theme(legend.text = element_text(face = "italic"))

# ── Bins globales (antimodos) ─────────────────────────────────────────────────
introgressCutoff <- 0.01
ampOut    <- amps(bedData$log2)
antimodes <- ampOut$Antimode[, 1]

binInfo <- data.frame(binName = character(), log2BinStart = numeric(),
                      log2BinEnd = numeric(), numWin = numeric(),
                      perData = numeric(), stringsAsFactors = FALSE)

binStart <- 0
for (i in seq_len(length(antimodes) + 1)) {
  binEnd  <- if (i > length(antimodes)) max(bedData$log2) else antimodes[i]
  numWin  <- sum(bedData$log2 >= binStart & bedData$log2 <= binEnd)
  binInfo <- rbind(binInfo, data.frame(
    binName      = paste0("bin", i - 1),
    log2BinStart = binStart,
    log2BinEnd   = binEnd,
    numWin       = numWin,
    perData      = round(100 * numWin / nrow(bedData), 2)
  ))
  binStart <- binEnd
}
binInfo$meanCovUpper <- (2 ^ binInfo$log2BinEnd) * checkMean

numBinsWanted <- max(which(binInfo$numWin >= 500), 1) + 1
binsWanted    <- binInfo[seq_len(numBinsWanted), ]
binsToWrite   <- binInfo[seq_len(numBinsWanted + 1), ]
binsToWrite$binName[1]                <- "belowThreshold"
binsToWrite$binName[numBinsWanted + 1]<- "aboveUpperBin"
binsToWrite$log2BinEnd[numBinsWanted + 1] <- binInfo$log2BinEnd[nrow(binInfo)]
binsToWrite$numWin[numBinsWanted + 1] <-
  sum(binInfo$numWin[seq(numBinsWanted + 1, nrow(binInfo))])
binsToWrite$perData[numBinsWanted + 1] <-
  sum(binInfo$perData[seq(numBinsWanted + 1, nrow(binInfo))])
binsToWrite$meanCovUpper[numBinsWanted + 1] <- binInfo$meanCovUpper[nrow(binInfo)]
write.table(binsToWrite, globBinFile, row.names = FALSE)

binOut <- bedData
binOut$covBin <- "aboveUpperBin"
for (i in seq_len(nrow(binsWanted)))
  binOut$covBin[bedData$log2 >= binsWanted$log2BinStart[i] &
                bedData$log2 <= binsWanted$log2BinEnd[i]] <- binsWanted$binName[i]
write.table(binOut, binWinFile, row.names = FALSE)

limitValue  <- binInfo$meanCovUpper[nrow(binsWanted) + 1]
lowerCutoff <- binsWanted$log2BinEnd[1]

# ── Resumen por especie ───────────────────────────────────────────────────────
spSum <- data.frame(species = character(), perAboveCutoff = character(),
                    contributes = logical(), stringsAsFactors = FALSE)

for (spName in levels(uniSpecies)) {
  spData    <- bedData[bedData$species == spName, ]
  nWinMin   <- nrow(spData) * introgressCutoff
  cutoffData<- spData[spData$log2 >= lowerCutoff, ]
  pct       <- paste0(round(nrow(cutoffData) / nrow(spData) * 100, 2), "%")
  spSum     <- rbind(spSum,
                     data.frame(species        = spName,
                                perAboveCutoff = pct,
                                contributes    = nrow(cutoffData) >= nWinMin))
}
write.table(spSum, binSumFile, row.names = FALSE, quote = FALSE)

sigSpecies <- factor(spSum$species[spSum$contributes],
                     levels = spSum$species[spSum$contributes])

# ── Preparar datos de especies significativas ─────────────────────────────────
bedData$meanValueLimited <- pmin(bedData$meanValue, limitValue)

stepSize <- bedData$Genome_Pos[2]
sigTable <- NULL
prevEnd  <- 0
colors_sig <- c()
speciesBreaks <- numeric()
speciesLabelPos <- numeric()
chrBreaks <- numeric()
prevChrEnd <- 0

for (spName in levels(uniSpecies)) {
  spData    <- bedData[bedData$species == spName, ]
  spStart   <- spData$Genome_Pos[1]
  spEnd     <- spData$Genome_Pos[nrow(spData)]
  speciesBreaks   <- c(speciesBreaks, spStart)
  speciesLabelPos <- c(speciesLabelPos, spStart + (spEnd - spStart) / 2)

  if (spName %in% levels(sigSpecies)) {
    spTable <- spData
    spTable$speciesPos <- spTable$Genome_Pos - prevEnd
    sigTable <- rbind(sigTable, spTable)
    colors_sig[gsub("_", " ", spName)] <- colors[spName]
  }

  for (chr in unique(spData$chrom)) {
    cd     <- spData[spData$chrom == chr, ]
    chrEnd <- cd$Genome_Pos[nrow(cd)]
    if ((chrEnd - prevChrEnd) >= 10 * stepSize)
      chrBreaks <- c(chrBreaks, chrEnd)
    prevChrEnd <- chrEnd
  }
  prevEnd <- spEnd
}

spcLabels <- gsub("_", "\n", levels(uniSpecies))
xaxGlobal <- if (length(levels(uniSpecies)) < 11) {
  scale_x_continuous(breaks = speciesLabelPos, labels = spcLabels,
                     name = "Genome Position", limits = c(0, NA))
} else {
  scale_x_continuous(breaks = NULL, labels = NULL,
                     name = "Genome Position", limits = c(0, NA))
}

fillLeg  <- scale_fill_manual(name = "Species", values = colors,
                              breaks = levels(uniSpecies),
                              labels = gsub("_", " ", levels(uniSpecies)))
colorLeg <- scale_color_manual(name = "Species", values = colors,
                               breaks = levels(uniSpecies),
                               labels = gsub("_", " ", levels(uniSpecies)))
vertLines <- geom_vline(xintercept = speciesBreaks)
hBinAvg   <- geom_hline(yintercept = binsWanted$meanCovUpper,
                         color = "gray30", linetype = "dotdash")
hBinLog   <- geom_hline(yintercept = binsWanted$log2BinEnd,
                         color = "gray30", linetype = "dotdash")

# ── PDF distribuciones ────────────────────────────────────────────────────────
pdf(plotCovDist)
print(ggplot(bedData, aes(x = log2)) +
  geom_density(aes(y = after_stat(scaled))) +
  ggtitle(paste(prefix, "Windowed log2 Coverage Density")) +
  coord_flip() + thm +
  geom_vline(xintercept = binsWanted$log2BinEnd,
             color = "gray30", linetype = "dotdash"))

print(ggplot(bedData, aes(x = log2, fill = species)) +
  geom_density(alpha = .75) + fillScale +
  ggtitle(paste(prefix, "Windowed log2 Coverage by Species")) +
  coord_flip() + thm +
  geom_vline(xintercept = binsWanted$log2BinEnd,
             color = "gray30", linetype = "dotdash"))

print(ggplot(bedData, aes(x = species, y = log2, fill = species)) +
  geom_violin(scale = "width") + fillScale +
  ggtitle(paste(prefix, "Windowed log2 Coverage Violin")) +
  scale_x_discrete(limits = levels(uniSpecies)) + thm +
  geom_hline(yintercept = binsWanted$log2BinEnd,
             color = "gray30", linetype = "dotdash"))

if (!is.null(sigTable) && nrow(sigTable) > 0) {
  sigData <- bedData[bedData$species %in% levels(sigSpecies), ]
  print(ggplot(sigData, aes(x = log2, fill = species)) +
    geom_density(alpha = .75) + fillScale + colorScale +
    ggtitle(paste(prefix, "log2 Coverage Density (contributing species)")) +
    coord_flip() + thm +
    geom_vline(xintercept = binsWanted$log2BinEnd,
               color = "gray30", linetype = "dotdash"))
}
dev.off()

# ── PDF plots de cobertura con bins ───────────────────────────────────────────
if (!is.null(sigTable) && nrow(sigTable) > 0) {
  sigTable$species <- gsub("_", " ", sigTable$species)
  sigLevels <- gsub("_", " ", levels(sigSpecies))

  endPos <- max(sigTable$speciesPos)
  breakPos <- seq(0, endPos, by = endPos / 10)
  breakLabels <- vapply(breakPos, function(lv) {
    if      (lv >= 1e9) paste0(round(lv / 1e9, 3), " Gb")
    else if (lv >= 1e6) paste0(round(lv / 1e6, 3), " Mb")
    else if (lv >= 1e3) paste0(round(lv / 1e3, 3), " Kb")
    else                as.character(round(lv))
  }, character(1))

  xaxSig  <- scale_x_continuous(breaks = breakPos, labels = breakLabels,
                                 name = "Genome Position", limits = c(0, endPos))
  yaxMean <- scale_y_continuous(name = "Average Depth", limits = c(0, NA))
  yaxLog  <- scale_y_continuous(name = "log2(Average Depth)", limits = c(0, NA))
  yaxLim  <- scale_y_continuous(name = "Average Depth (limited)",
                                limits = c(0, limitValue))
  fillSig <- scale_fill_manual(name = "Species", values = colors_sig,
                               breaks = sigLevels)
  colorSig<- scale_color_manual(name = "Species", values = colors_sig,
                                breaks = sigLevels)
  hBinAvgSig <- geom_hline(yintercept = binsWanted$meanCovUpper,
                            color = "gray30", linetype = "dotdash")
  hBinLogSig <- geom_hline(yintercept = binsWanted$log2BinEnd,
                            color = "gray30", linetype = "dotdash")

  # Líneas de cromosoma
  facetDF <- data.frame()
  endsDF  <- data.frame()
  prevChr2 <- 0
  for (spName in sigLevels) {
    spd <- sigTable[sigTable$species == spName, ]
    endsDF <- rbind(endsDF, data.frame(end = max(spd$speciesPos), species = spName))
    for (chr in unique(spd$chrom)) {
      cd  <- spd[spd$chrom == chr, ]
      cEnd <- cd$speciesPos[1]
      if ((cEnd - prevChr2) >= 20 * stepSize)
        facetDF <- rbind(facetDF,
                         data.frame(Breaks = cEnd, species = spName))
      prevChr2 <- cEnd
    }
  }
  plotVertEnd <- geom_vline(aes(xintercept = end), endsDF)
  plotVert    <- if (nrow(facetDF) > 0) {
    geom_vline(aes(xintercept = Breaks), facetDF, linetype = 2)
  } else {
    geom_blank()
  }

  pdf(plotCovBins, width = 14)
  for (yax in list(
    list(col = "meanValue",        scale = yaxMean, hbin = hBinAvgSig),
    list(col = "log2",             scale = yaxLog,  hbin = hBinLogSig),
    list(col = "meanValueLimited", scale = yaxLim,  hbin = hBinAvgSig)
  )) {
    col <- yax$col
    print(
      ggplot(transform(sigTable, species = factor(species, levels = sigLevels)),
             aes_string("speciesPos", col, colour = "species")) +
        yax$hbin + geom_point(size = 0.5) +
        facet_grid(species ~ .) + theme_classic() +
        colorSig + xaxSig + yax$scale +
        ggtitle(paste(prefix, "Avg depth of coverage")) +
        theme(legend.text = element_text(face = "italic"),
              strip.text  = element_text(face = "italic")) +
        plotVert + geom_vline(xintercept = 0) + plotVertEnd
    )
    print(
      ggplot(transform(sigTable, species = factor(species, levels = sigLevels)),
             aes_string("speciesPos")) +
        yax$hbin + geom_ribbon(aes_string(ymin = "0", ymax = col,
                                           fill = "species")) +
        facet_grid(species ~ .) + theme_classic() +
        fillSig + xaxSig + yax$scale +
        ggtitle(paste(prefix, "Avg depth of coverage")) +
        theme(legend.text = element_text(face = "italic"),
              strip.text  = element_text(face = "italic")) +
        plotVert + geom_vline(xintercept = 0) + plotVertEnd
    )
  }
  dev.off()
}

message("sppIDer_depthPlot.R completado.")
