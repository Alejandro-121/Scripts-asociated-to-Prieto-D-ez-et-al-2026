#!/usr/bin/env Rscript
# MQscores_sumPlot.R
# Uso: Rscript MQscores_sumPlot.R WORKDIR PREFIX
#
# Genera plots de porcentaje de lecturas mapeadas y calidad de mapeo (MQ).
# Salidas: PREFIX_MQsummary.txt, PREFIX_MQ_chiSquared.txt, PREFIX_plotMQ.pdf

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript MQscores_sumPlot.R WORKDIR PREFIX")
}

workingDir <- file.path(args[1], "")   # asegura trailing slash
strainName <- args[2]

# ── Archivos de entrada/salida ────────────────────────────────────────────────
mq_file    <- file.path(workingDir, paste0(strainName, "_MQ.txt"))
chiSqOut   <- file.path(workingDir, paste0(strainName, "_MQ_chiSquared.txt"))
summaryOut <- file.path(workingDir, paste0(strainName, "_MQsummary.txt"))
plotOut    <- file.path(workingDir, paste0(strainName, "_plotMQ.pdf"))

if (!file.exists(mq_file)) stop("No se encuentra: ", mq_file)

# ── Leer datos ────────────────────────────────────────────────────────────────
MQdf <- read.table(mq_file, header = TRUE, stringsAsFactors = FALSE)
levels(factor(MQdf$Species))

MQdf$Species[MQdf$Species == "*"] <- "Unmapped"
MQdf$Species <- factor(MQdf$Species)
uniSpecies   <- levels(MQdf$Species)

# Separador de etiquetas: salto de línea si hay pocas especies
spRe <- if (length(uniSpecies) <= 11) "\n" else " "

# ── Subconjuntos ──────────────────────────────────────────────────────────────
numReads       <- sum(MQdf$count)
maps           <- subset(MQdf, Species != "Unmapped")
mapsNonZero    <- subset(maps, MQscore > 0)
MQscoresZero   <- subset(MQdf, MQscore == 0)
MQscoresNonZero<- subset(MQdf, MQscore > 0)
propZero       <- (sum(MQscoresZero$count) / numReads) * 100

# ── Resumen escrito en texto ──────────────────────────────────────────────────
sink(summaryOut)
cat(strainName, "Num reads =", numReads, "\n")
cat(strainName, "Num mapped reads =", sum(maps$count), "\n")
cat(strainName, "Unmapped reads =", round(propZero, 2), "%\n")
cat(strainName, "Average MQ =", weighted.mean(MQdf$MQscore, MQdf$count), "\n")
cat(strainName, "Median MQ  =", median(rep(MQdf$MQscore, times = MQdf$count)), "\n")
cat("Species\tTotal mapped reads\t% of all reads\t% Nonzero mapped reads",
    "\tAll average MQ\tNonZero average MQ\tAll Median MQ\tNonzero Median MQ\n")

MQprop    <- data.frame(Species = character(), Percentage = numeric(),
                        stringsAsFactors = FALSE)
MQpropPos <- data.frame(Species = character(), Percentage = numeric(),
                        stringsAsFactors = FALSE)
mapsExisting <- maps
violin_df    <- data.frame()

for (species in uniSpecies) {
  spMQ     <- subset(MQdf,          Species == species)
  spNonZero<- subset(MQscoresNonZero, Species == species)

  totalCnt <- sum(spMQ$count)
  nzCnt    <- sum(spNonZero$count)
  propSp   <- (totalCnt / numReads) * 100
  propNZ   <- ifelse(sum(MQscoresNonZero$count) > 0,
                    (nzCnt / sum(MQscoresNonZero$count)) * 100, 0)

  meanAll  <- weighted.mean(spMQ$MQscore,      spMQ$count)
  meanNZ   <- weighted.mean(spNonZero$MQscore, spNonZero$count)
  medAll   <- median(rep(spMQ$MQscore,      times = spMQ$count))
  medNZ    <- median(rep(spNonZero$MQscore, times = spNonZero$count))

  for (v in c("meanAll","meanNZ","medAll","medNZ"))
    if (is.na(get(v))) assign(v, 0)

  cat(species, "\t", totalCnt, "\t", round(propSp,2), "%\t",
      round(propNZ,2), "%\t", round(meanAll,2), "\t", round(meanNZ,2),
      "\t", round(medAll,2), "\t", round(medNZ,2), "\n")

  MQprop    <- rbind(MQprop,    data.frame(Species = species, Percentage = propSp))
  MQpropPos <- rbind(MQpropPos, data.frame(Species = species, Percentage = propNZ))
  if (nzCnt > 0) violin_df <- rbind(violin_df, spMQ)
  if (totalCnt == 0)
    mapsExisting <- mapsExisting[mapsExisting$Species != species, ]
}
sink()

# Quitar "Unmapped" del gráfico sin unmapped
MQpropPos <- MQpropPos[MQpropPos$Species != "Unmapped", ]

# ── Chi-cuadrado ──────────────────────────────────────────────────────────────
# Construir matriz de conteos [MQ 0-60] x [especie]
countMat <- matrix(0, nrow = 61, ncol = length(uniSpecies),
                   dimnames = list(0:60, uniSpecies))
for (sp in uniSpecies) {
  spd <- subset(MQdf, Species == sp)
  for (i in seq_len(nrow(spd)))
    countMat[spd$MQscore[i] + 1, sp] <- spd$count[i]
}

unmapCnt <- countMat[1, "Unmapped"]
posMat   <- countMat[-1, colnames(countMat) != "Unmapped", drop = FALSE]
posSum   <- c(colSums(posMat), Unmapped = unmapCnt)
cqAll    <- chisq.test(posSum)

mq60row  <- posMat[nrow(posMat), , drop = FALSE]
cq60     <- chisq.test(mq60row)

sink(chiSqOut)
cat(strainName, "count of reads mapped - Chi Squared P-value", cqAll$p.value, "\n")
cat(names(cqAll$residuals), "\n")
cat(cqAll$residuals, "\n")
cat(strainName, "MQ60 - Chi Squared P-value", cq60$p.value, "\n")
cat(names(cq60$residuals), "\n")
cat(cq60$residuals, "\n")
sink()

# ── Colores ───────────────────────────────────────────────────────────────────
pdf(NULL)   # abrir dispositivo nulo para que palette() funcione
nSpc <- length(uniSpecies)
colorList <- if (nSpc > 2) rainbow(nSpc - 1) else c("red")
colors <- setNames(
  c(colorList[seq_along(uniSpecies[-which(uniSpecies == "Unmapped")])],
    "darkgrey"),
  c(uniSpecies[uniSpecies != "Unmapped"], "Unmapped")
)
dev.off()

# Sustituir _ por separador de línea en etiquetas
relabel <- function(x) gsub("_", spRe, x)
MQdf$Species       <- relabel(MQdf$Species)
MQprop$Species     <- relabel(MQprop$Species)
MQpropPos$Species  <- relabel(MQpropPos$Species)
mapsExisting$Species <- relabel(mapsExisting$Species)
spcLabels          <- relabel(uniSpecies)
names(colors)      <- relabel(names(colors))

fillLeg  <- scale_fill_manual(name = "Species",   values = colors, breaks = spcLabels)
colorLeg <- scale_colour_manual(name = "Species", values = colors, breaks = spcLabels)
theme_it <- theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        axis.text.x = element_text(face = "italic"))
themeUnmap <- if (length(uniSpecies) >= 11) {
  theme_bw() + theme(axis.text.x = element_blank(),
                     legend.text = element_text(face = "italic"))
} else {
  theme_it
}

# ── PDF ───────────────────────────────────────────────────────────────────────
pdf(plotOut, compress = TRUE, width = 10)

ggplot(MQprop, aes(x = factor(Species, levels = spcLabels))) +
  geom_bar(aes(fill = Species, weight = Percentage)) +
  ggtitle(paste(strainName, "Mapping bar plot w/ unmapped")) +
  labs(x = "Species", y = "Percentage") +
  fillLeg + themeUnmap + scale_y_continuous(limits = c(0, 100))

ggplot(MQpropPos, aes(x = factor(Species, levels = spcLabels[spcLabels != relabel("Unmapped")]))) +
  geom_bar(aes(fill = Species, weight = Percentage)) +
  ggtitle(paste(strainName, "Mapping bar plot w/out unmapped")) +
  labs(x = "Species", y = "Percentage") +
  fillLeg + theme_it + scale_y_continuous(limits = c(0, 100))

if (nrow(mapsExisting) > 0) {
  suppressWarnings(
    print(
      ggplot(mapsExisting,
             aes(factor(Species, levels = spcLabels), MQscore, weight = count)) +
        geom_violin(bw = 1, scale = "count", draw_quantiles = c(.25, .5, .75),
                    aes(fill = Species)) +
        fillLeg + colorLeg +
        ggtitle(paste(strainName, "Mapping quality of species with mapped reads")) +
        labs(x = "Species") + theme_it
    )
  )
}

dev.off()
message("MQscores_sumPlot.R completado: ", plotOut)
