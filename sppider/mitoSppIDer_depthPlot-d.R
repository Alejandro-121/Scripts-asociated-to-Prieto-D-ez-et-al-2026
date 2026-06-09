#!/usr/bin/env Rscript
# mitoSppIDer_depthPlot-d.R
# Uso: Rscript mitoSppIDer_depthPlot-d.R OUTDIR SAMPLE [GFF]
#
# Genera gráficas de profundidad de cobertura para genomas mitocondriales.
# Si se proporciona un GFF combinado, sombrea las regiones CDS.

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Uso: Rscript mitoSppIDer_depthPlot-d.R OUTDIR SAMPLE [GFF]")

workingDir <- file.path(args[1], "")
prefix     <- args[2]
gffFile    <- if (length(args) >= 3) args[3] else NULL

# ── Archivos de entrada ───────────────────────────────────────────────────────
winFile <- file.path(workingDir, paste0(prefix, "_winAvgDepth-d.txt"))
if (!file.exists(winFile)) stop("No se encuentra: ", winFile)

plotOut <- file.path(workingDir, paste0(prefix, "_mitoDepthPlot.pdf"))

# ── Leer datos de cobertura por ventana ───────────────────────────────────────
bedData <- read.table(winFile, header = FALSE, skip = 2,
                      col.names = c("Genome_Pos", "species", "chrom",
                                    "start", "end", "meanValue",
                                    "log2rawCov", "max", "median"))

checkMean   <- mean(bedData$meanValue)
bedData$log2 <- pmax(0, log2(bedData$meanValue / checkMean))

pdf(NULL)
uniSpecies <- factor(unique(bedData$species), levels = unique(bedData$species))
colorList  <- if (length(uniSpecies) > 1) rainbow(length(uniSpecies)) else "red"
dev.off()

colors    <- setNames(colorList, levels(uniSpecies))
fillScale <- scale_fill_manual(name  = "Species", values = colors,
                               breaks = levels(uniSpecies),
                               labels = gsub("_", " ", levels(uniSpecies)))
colScale  <- scale_color_manual(name = "Species", values = colors,
                                breaks = levels(uniSpecies),
                                labels = gsub("_", " ", levels(uniSpecies)))
thm <- theme_classic() +
  theme(legend.text = element_text(face = "italic"),
        axis.text.x = element_text(face = "italic"))

# ── GFF (opcional) ────────────────────────────────────────────────────────────
gffData <- NULL
if (!is.null(gffFile) && file.exists(gffFile)) {
  gffData <- read.table(gffFile, header = TRUE, stringsAsFactors = FALSE)
  message("GFF cargado: ", nrow(gffData), " regiones CDS")
} else if (!is.null(gffFile)) {
  warning("GFF no encontrado, se ignorará: ", gffFile)
}

# ── Acumulación de posiciones genómicas ───────────────────────────────────────
speciesBreaks   <- numeric()
speciesLabelPos <- numeric()
chrBreaks       <- numeric()
stepSize        <- if (nrow(bedData) > 1) bedData$Genome_Pos[2] else 1

prevChrEnd <- 0
for (sp in levels(uniSpecies)) {
  spData <- bedData[bedData$species == sp, ]
  spStart <- spData$Genome_Pos[1]
  spEnd   <- spData$Genome_Pos[nrow(spData)]
  speciesBreaks   <- c(speciesBreaks, spStart)
  speciesLabelPos <- c(speciesLabelPos, spStart + (spEnd - spStart) / 2)

  for (chr in unique(spData$chrom)) {
    cd     <- spData[spData$chrom == chr, ]
    chrEnd <- cd$Genome_Pos[nrow(cd)]
    if ((chrEnd - prevChrEnd) >= 10 * stepSize)
      chrBreaks <- c(chrBreaks, chrEnd)
    prevChrEnd <- chrEnd
  }
}

spcLabels  <- gsub("_", "\n", levels(uniSpecies))
maxPos     <- max(bedData$Genome_Pos)
breakPos   <- seq(0, maxPos, length.out = 6)
breakLabels <- vapply(breakPos, function(lv) {
  if      (lv >= 1e6) paste0(round(lv / 1e6, 2), " Mb")
  else if (lv >= 1e3) paste0(round(lv / 1e3, 2), " Kb")
  else                as.character(round(lv))
}, character(1))

xax       <- scale_x_continuous(breaks = breakPos, labels = breakLabels,
                                 name = "Posición genómica", limits = c(0, NA))
vertSpc   <- geom_vline(xintercept = speciesBreaks, color = "black")
vertChr   <- if (length(chrBreaks) > 0) {
  geom_vline(xintercept = chrBreaks, linetype = "dashed", color = "gray50")
} else {
  geom_blank()
}

# Sombrado CDS (si hay GFF)
cds_rects <- NULL
if (!is.null(gffData) && nrow(gffData) > 0) {
  # Mapear posición genómica de cada CDS usando Genome_Pos acumulado
  # Construir tabla chrom → offset genómico
  chrOffset <- data.frame()
  for (sp in levels(uniSpecies)) {
    spData <- bedData[bedData$species == sp, ]
    for (chr in unique(spData$chrom)) {
      cd <- spData[spData$chrom == chr, ]
      chrOffset <- rbind(chrOffset,
                         data.frame(chrom  = chr,
                                    species = sp,
                                    offset = cd$Genome_Pos[1] - cd$start[1]))
    }
  }

  gffMapped <- merge(gffData, chrOffset,
                     by.x = "Species", by.y = "species", all.x = FALSE)
  if (nrow(gffMapped) > 0) {
    gffMapped$xmin <- gffMapped$Start  + gffMapped$offset
    gffMapped$xmax <- gffMapped$End    + gffMapped$offset
    cds_rects <- geom_rect(data = gffMapped,
                            aes(xmin = xmin, xmax = xmax,
                                ymin = -Inf, ymax = Inf),
                            fill = "gray80", alpha = 0.4,
                            inherit.aes = FALSE)
  }
}

# ── PDF ───────────────────────────────────────────────────────────────────────
pdf(plotOut, width = 14)

for (yvar in c("meanValue", "log2")) {
  ylab <- if (yvar == "meanValue") "Profundidad media" else "log2(profundidad / media global)"

  # Puntos
  g <- ggplot(bedData, aes_string("Genome_Pos", yvar, colour = "species")) +
    xax + scale_y_continuous(name = ylab) +
    colScale + vertSpc + vertChr +
    ggtitle(paste(prefix, "—", ylab)) + thm
  if (!is.null(cds_rects)) g <- g + cds_rects
  g <- g + geom_point(size = 0.4)
  print(g)

  # Ribbons por especie
  g2 <- ggplot(bedData, aes_string("Genome_Pos")) +
    xax + scale_y_continuous(name = ylab) +
    fillScale + vertSpc + vertChr +
    ggtitle(paste(prefix, "—", ylab, "(área)")) + thm
  if (!is.null(cds_rects)) g2 <- g2 + cds_rects
  g2 <- g2 + geom_ribbon(aes_string(ymin = "0", ymax = yvar, fill = "species"))
  print(g2)
}

dev.off()
message("mitoSppIDer_depthPlot-d.R completado: ", plotOut)
