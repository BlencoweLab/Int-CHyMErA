### Identify hits with increased cell growth among regions targeted in Int-CHyMErA screen
###
### To run after scoring hits with decreased growth (because we tune the hit threshold based on that)
### - For each time point/treatment vs. T0, find exon-deletion guide pairs with a dropout < 5th percentile of intergenic-intergenic guides. 
### - Consider only guides where none of the intronic-intergenic controls drop out.
### - Identify ‘fitness’ and ‘non-fitness’ genes as the bottom 20% and top 60% of genes, respectively, of the mean of exonic Cas9 targeting guides.
### - Based on exon-deletion pairs, find the threshold of the fraction of overlapping guide pairs that best separates 
###   positive controls (exons in fitness genes) from negative controls (exons in non-fitness genes). Require >= 2 changing guide pairs among >= 2 available.
### - Then Apply same theshold and criteria to intronic regions to identify regions with dropout phenotype.

library(parallel)
source("./R/CHyInt_regionScoringFunctions.R")


### Parameters

inDir  <- "./input"
outDir <- "./output/increase"

performance.HAP1.T12 <- "./output/dropout/tables_HAP1.exonDel/HAP1.exonDel_T12_NT:T0_Performance.csv" # we get the threshold of effective hgRNAs from here
performance.HAP1.T18 <- "./output/dropout/tables_HAP1.exonDel/HAP1.exonDel_T18_NT:T0_Performance.csv" 
performance.RPE1.T21 <- "./output/dropout/tables_RPE1.exonDel/RPE1.exonDel_T21_NT:T0_Performance.csv"
performance.RPE1.T27 <- "./output/dropout/tables_RPE1.exonDel/RPE1.exonDel_T27_NT:T0_Performance.csv"

cores <- 4  # CPUs
minInfPairs <- 2          # Minimum number of informative guide pairs to score a region, where informative means not NA and
                          # either gRNA within the hgRNA has not significant LFC when paired with an intergenic partner
minChangePairs <- 2       # Minimum number of changing guide pairs to score a region as a hit, where changing means 
                          # that the LFC is < the 5th percentile of the intergenic-intergenic guides for dropout or 
                          # higher than the 95th percentile for enrichment.
ii.quant  <- 0.95         # Threshold, in LFC of this quantil of the LFC of intergenic-intergenic guides, to call an hgRNA
                          # significantly changing. In the case of scoring increase, this is the lower bound.
fit.quant <- c(0.9, 0.8)  # Bounds (in quantile of the LFC distribution of genes knockouts) on fitness and non-fitness genes.
                          # In the case of scoring increase, the first is the lower bound and the second is the upper bound.


### Function defs

getThreshold <- function(x) {
    ### Read the performance table from a screen scoring and extract the thresold with the best
    ### separation of hits in fitness vs. non-fitness genes
    dat <- read.csv(x, check.names = F)
    maxDelta <- which(
        dat$"DoubleX fitness" - dat$"DoubleX nonfitness" == 
        max(dat$"DoubleX fitness" - dat$"DoubleX nonfitness", na.rm=T)
    )
    dat$FractionCutoff[maxDelta]
}                     


if (!dir.exists(outDir))  {dir.create(outDir, recursive = TRUE)}
if (!dir.exists(file.path(outDir, "tables_HAP1.exonDel")))   {dir.create(file.path(outDir, "tables_HAP1.exonDel"), recursive = TRUE)}
if (!dir.exists(file.path(outDir, "tables_RPE1.exonDel")))   {dir.create(file.path(outDir, "tables_RPE1.exonDel"), recursive = TRUE)}
if (!dir.exists(file.path(outDir, "tables_HAP1.intronDel"))) {dir.create(file.path(outDir, "tables_HAP1.intronDel"), recursive = TRUE)}
if (!dir.exists(file.path(outDir, "tables_RPE1.intronDel"))) {dir.create(file.path(outDir, "tables_RPE1.intronDel"), recursive = TRUE)}


### Read LFC

deH <- read.delim(file.path(outDir, "../edgeR/edgeR_HAP1.tab.gz"), check.names = F)
deR <- read.delim(file.path(outDir, "../edgeR/edgeR_RPE1.tab.gz"), check.names = F)
info <- read.delim(file.path(outDir, "../edgeR/GuideInfo.tab.gz"))

contrH <- list(
    c("T0","T12_NT"),
    c("T0","T18_NT")
)
contrR <- list(
    c("T0","T21_NT"),
    c("T0","T27_NT")
)


### Generate gene LFC based on exon-targeting Cas9 guides (TKOv3 and additional guides)

gFitnessH <- lapply(contrH, getGeneFitness, de = deH)
gFitnessR <- lapply(contrR, getGeneFitness, de = deR)


### Get percentage of hgRNAs with dropout per deleted exon region from the hits with decreased fitness

incrPairFracH <- incrPairFracR <- numeric(length(contrH))
names(incrPairFracH) <- sapply(contrH, FUN = function(x) {paste(rev(x), collapse = ":")})
names(incrPairFracR) <- sapply(contrR, FUN = function(x) {paste(rev(x), collapse = ":")})

incrPairFracH[1] <- getThreshold(performance.HAP1.T12)
incrPairFracH[2] <- getThreshold(performance.HAP1.T18)

incrPairFracR[1] <- getThreshold(performance.RPE1.T21)
incrPairFracR[2] <- getThreshold(performance.RPE1.T27)



### Use the treshold from exon deletions/decrease to call hits on intron deletions/increase

intronIncrH <- mcmapply(
    FUN = getRegionDropFrac, contrH, gFitness = gFitnessH, SIMPLIFY = FALSE, 
    MoreArgs = list(
        de = deH, info = info, region = "introns",
        direction="increase", ii.quant = ii.quant, fit.quant = fit.quant,
        goodPairsOnly = TRUE, minInfPairs = minInfPairs, minChangePairs = minChangePairs
    ), 
    mc.cores = cores
)

intronIncrR <- mcmapply(
    FUN = getRegionDropFrac, contrR, gFitness = gFitnessR, SIMPLIFY = FALSE, 
    MoreArgs = list(
        de = deR, info = info, region = "introns",
        direction="increase", ii.quant = ii.quant, fit.quant = fit.quant,
        goodPairsOnly = TRUE, minInfPairs = minInfPairs, minChangePairs = minChangePairs
    ), 
    mc.cores = cores
)


### Plot performance on intron deletion (note that this is not necessarily expected to yield in a separation)

pdf(file.path(outDir, "ScreenPerformance_HAP1.intronDel.pdf"), width = 9 , height = 5)
for (i in 1:length(contrH)) {
    plotRegionDropFrac(intronIncrH[[i]], fracOverride = incrPairFracH[i])
    legend("topright", paste(rev(contrH[[i]]), collapse=":"), bty = "n", text.col = "grey60")
}
dev.off()

pdf(file.path(outDir, "ScreenPerformance_RPE1.intronDel.pdf"), width = 9 , height = 5)
for (i in 1:length(contrR)) {
    plotRegionDropFrac(intronIncrR[[i]], fracOverride = incrPairFracR[i])
    legend("topright", paste(rev(contrR[[i]]), collapse=":"), bty="n", text.col="grey60")
}
dev.off()


### Save intron deletion results

for (i in 1:length(contrH)) {
    saveTables(
        intronIncrH[[i]], pairFracCO = incrPairFracH[i], gFitness = gFitnessH[[i]], 
        outDir = file.path(outDir, "tables_HAP1.intronDel"), 
        outName = paste0("HAP1.intronDel_", paste(rev(contrH[[i]]), collapse = ":"))
    )
    saveGuidePairs(
        contrH[[i]], de = deH, gFitness = gFitnessH[[i]], region="introns",
        direction="increase", ii.quant = ii.quant, fit.quant = fit.quant,
        outDir = file.path(outDir, "tables_HAP1.intronDel"), 
        outName = paste0("HAP1.intronDel_", paste(rev(contrH[[i]]), collapse = ":"))
    )
}

for (i in 1:length(contrR)) {
    saveTables(
        intronIncrR[[i]], pairFracCO = incrPairFracR[i], gFitness = gFitnessR[[i]],
        outDir = file.path(outDir, "tables_RPE1.intronDel"), 
        outName = paste0("RPE1.intronDel_", paste(rev(contrR[[i]]), collapse=":"))
    )
    saveGuidePairs(
        contrR[[i]], de = deR, gFitness = gFitnessR[[i]], region = "introns",
        direction="increase", ii.quant = ii.quant, fit.quant = fit.quant,
        outDir = file.path(outDir, "tables_RPE1.intronDel"),
        outName = paste0("RPE1.intronDel_", paste(rev(contrR[[i]]), collapse = ":")))
}

