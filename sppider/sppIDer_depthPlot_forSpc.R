#!/usr/bin/env Rscript
# sppIDer_depthPlot_forSpc.R
# Uso: Rscript sppIDer_depthPlot_forSpc.R WORKDIR PREFIX

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Uso: Rscript sppIDer_depthPlot_forSpc.R WORKDIR PREFIX")

workingDir <- file.path(args[1], "")
prefix     <- args[2]

# ── Archivos ──────────────────────────────────────────────────────────────────
dFile <- file.path(workingDir, paste0(prefix, "_speciesAvgDepth-d.txt"))
gFile <- file.path(workingDir, paste0(prefix, "_speciesAvgDepth-g.txt"))
dataFile <- if (file.exists(dFile)) dFile else if (file.exists(gFile)) gFile else
  stop("No se encuentra ningún archivo speciesAvgDepth para: ", prefix)

plotOut <- file.path(workingDir, paste0(prefix, "_speciesDepth.pdf"))

bedData <- read.table(dataFile, header = FALSE, skip = 2,
                      col.names = c("Genome_Pos", "species", "end",
                                    "meanValue", "relativeMean", "max", "median"))

# ── Especies ──────────────────────────────────────────────────────────────────
uniSpecies <- unique(vapply(strsplit(as.character(bedData$species), "-"),
                            `[`, character(1), 1))

checkMean <- mean(bedData$meanValue)
maxMean   <- max(bedData$meanValue)
bedData$log2 <- pmax(0, log2(bedData$meanValue / checkMean))

pdf(NULL)
colorList <- if (length(uniSpecies) > 1) rainbow(length(uniSpecies)) else "red"
dev.off()

colors <- setNames(colorList, uniSpecies)
spcLabels <- gsub("_", "\n", uniSpecies)

speciesBreaks <- numeric()
labelPos      <- numeric()
spcLabeled    <- data.frame()

for (k in seq_along(uniSpecies)) {
  spName <- uniSpecies[k]
  spc    <- bedData[grepl(spName, bedData$species, fixed = TRUE), ]
  spc    <- cbind(species_name = spName, spc)
  spcLabeled <- rbind(spcLabeled, spc)
  speciesBreaks <- c(speciesBreaks, spc$Genome_Pos[1])
  labelPos <- c(labelPos, spc$Genome_Pos[1] + spc$end[nrow(spc)] / 2)
}

# Duplicar para que geom_ribbon dibaje escalones correctos
GenomeEnds <- c(spcLabeled$Genome_Pos[-1] - 1,
                tail(spcLabeled$Genome_Pos, 1) + tail(spcLabeled$end, 1))
spcDup <- spcLabeled
spcDup$Genome_Pos <- GenomeEnds
spcBoth <- rbind(spcLabeled, spcDup)

fillLeg  <- scale_fill_manual(name = "Species", values = colors,
                              breaks = uniSpecies,
                              labels = gsub("_", " ", uniSpecies))
xax <- if (length(uniSpecies) < 11) {
  scale_x_continuous(breaks = labelPos, labels = spcLabels,
                     name = "Genome Position", limits = c(0, NA))
} else {
  scale_x_continuous(breaks = NULL, labels = NULL,
                     name = "Genome Position", limits = c(0, NA))
}

vertLines <- geom_vline(xintercept = speciesBreaks)
thm <- theme_classic() +
  theme(axis.text.x  = element_text(face = "italic"),
        legend.text  = element_text(face = "italic"))

pdf(plotOut, width = 14)

ggplot(spcBoth, aes(x = Genome_Pos)) +
  geom_ribbon(aes(ymin = 0, ymax = meanValue, fill = species_name)) +
  fillLeg + xax +
  scale_y_continuous(name = "Average Depth",
                     limits = c(0, maxMean * 1.1)) +
  ggtitle(paste(prefix, "Avg depth of coverage")) +
  vertLines + geom_abline(intercept = 0, slope = 0) + thm

ggplot(spcBoth, aes(x = Genome_Pos)) +
  geom_ribbon(aes(ymin = 0, ymax = log2, fill = species_name)) +
  fillLeg + xax +
  scale_y_continuous(name = "log2(avg/whole genome avg)",
                     limits = c(0, max(spcLabeled$log2) * 1.1)) +
  ggtitle(paste(prefix, "log2 Mean Avg depth of coverage")) +
  vertLines + geom_abline(intercept = 0, slope = 0) + thm

dev.off()
message("sppIDer_depthPlot_forSpc.R completado: ", plotOut)
