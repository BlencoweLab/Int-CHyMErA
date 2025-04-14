### Identify hits with decreasing cell growth among regions targeted in Int-CHyMErA screen
###
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
outDir <- "./output/dropout"

cores <- 4                # CPUs
minInfPairs <- 2          # Minimum number of informative guide pairs to score a region, where informative means not NA and
                          # either gRNA within the hgRNA has not significant LFC when paired with an intergenic partner
minChangePairs <- 2       # Minimum number of changing guide pairs to score a region as a hit, where changing means 
                          # that the LFC is < the 5th percentile of the intergenic-intergenic guides for dropout or 
                          # higher than the 95th percentile for enrichment.
ii.quant  <- 0.05         # Threshold, in LFC of this quantile of the LFC of intergenic-intergenic guides, to call an hgRNA
                          # significantly changing. In the case of scoring decrease, this is the upper bound.
fit.quant <- c(0.2, 0.4)  # Bounds (in quantile of the LFC distribution of genes knockouts) on fitness and non-fitness genes.
                          # In the case of scoring decrease, the first is the upper bound and the second is the upper lower bound.
goodPairsOnly <- TRUE     # If TRUE, only hgRNAs where neither of the component gRNAs results in a significant
                          # fitness changes when paired with a distal intergenic control, presumably resulting in a single cut
                          # that is repaired leaving small InDels.


if (!dir.exists(outDir))  {dir.create(outDir, recursive = TRUE)}
if (!dir.exists(file.path(outDir, "tables_HAP1.exonDel")))   {dir.create(file.path(outDir, "tables_HAP1.exonDel"), recursive = TRUE)}
if (!dir.exists(file.path(outDir, "tables_RPE1.exonDel")))   {dir.create(file.path(outDir, "tables_RPE1.exonDel"), recursive = TRUE)}
if (!dir.exists(file.path(outDir, "tables_HAP1.intronDel"))) {dir.create(file.path(outDir, "tables_HAP1.intronDel"), recursive = TRUE)}
if (!dir.exists(file.path(outDir, "tables_RPE1.intronDel"))) {dir.create(file.path(outDir, "tables_RPE1.intronDel"), recursive = TRUE)}


### Read LFC

deH <- read.delim(file.path(outDir, "../edgeR/edgeR_HAP1.tab.gz"), check.names = F)
deR <- read.delim(file.path(outDir, "../edgeR/edgeR_RPE1.tab.gz"), check.names = F)
info <- read.delim(file.path(inDir, "IntCHyMErA_library.design.tab.gz"))
info <- info[match(deH$Guide.ID, info$Guide.ID),]

contrH <- list(
    c("T0","T12_NT"),
    c("T0","T18_NT")
)
contrR <- list(
    c("T0","T21_NT"),
    c("T0","T27_NT")
)


### Generate gene LFC based on exon-targeting Cas9 guides (TKOv3 and additional guides)

gFitnessH <- lapply(contrH, getGeneFitness, de = deH, info)
gFitnessR <- lapply(contrR, getGeneFitness, de = deR, info)

for (i in 1:length(contrH)) {
    write.csv(
        gFitnessH[[i]], row.names = F,
        file = file.path(outDir, paste0("HAP1.geneDropout.Cas9exonic_", paste(rev(contrH[[i]]), collapse=":"), ".csv"))
    )
}
for (i in 1:length(contrR)) {
    write.csv(
        gFitnessR[[i]], row.names = F,
        file = file.path(outDir, paste0("RPE1.geneDropout.Cas9exonic_", paste(rev(contrR[[i]]), collapse=":"), ".csv"))
    )
}


### Get percentage of hgRNAs with dropout per deleted exon region 

exonDropH <- mcmapply(
    FUN = getRegionDropFrac, contrH, gFitness=gFitnessH, SIMPLIFY = FALSE,
    MoreArgs = list(
        de = deH, info = info, 
        direction="decrease", ii.quant = ii.quant, fit.quant = fit.quant,
        goodPairsOnly = goodPairsOnly, minInfPairs = minInfPairs, minChangePairs = minChangePairs
    ),
    mc.cores = cores
)
exonDropR <- mcmapply(
    FUN = getRegionDropFrac, contrR, gFitness=gFitnessR, SIMPLIFY = FALSE,
    MoreArgs = list(
        de = deR, info = info, 
        direction="decrease", ii.quant = ii.quant, fit.quant = fit.quant,
        goodPairsOnly = goodPairsOnly, minInfPairs = minInfPairs, minChangePairs = minChangePairs
    ), 
    mc.cores = cores
)


### Plot performance on exon deletion and determine optimal threshold

dropPairFracH <- dropPairFracR <- numeric(length(contrH))
names(dropPairFracH) <- sapply(contrH, FUN = function(x) {paste(rev(x), collapse = ":")})
names(dropPairFracR) <- sapply(contrR, FUN = function(x) {paste(rev(x), collapse = ":")})


pdf(file.path(outDir, "ScreenPerformance_HAP1.exonDel.pdf"), width = 9, height = 5)
for (i in 1:length(contrH)) {
    dropPairFracH[i] <- plotRegionDropFrac(exonDropH[[i]])
    legend("topright", paste(rev(contrH[[i]]), collapse=":"), bty = "n", text.col = "grey60")
}
dev.off()

pdf(file.path(outDir, "ScreenPerformance_RPE1.exonDel.pdf"), width = 9, height = 5)
for (i in 1:length(contrR)) {
    dropPairFracR[i] <- plotRegionDropFrac(exonDropR[[i]])
    legend("topright", paste(rev(contrR[[i]]), collapse=":"), bty = "n", text.col = "grey60")
}
dev.off()


### Save exon deletion results
for (i in 1:length(contrH)) {
    saveTables(
        exonDropH[[i]], pairFracCO = dropPairFracH[i], gFitness = gFitnessH[[i]],
        outDir = file.path(outDir, "tables_HAP1.exonDel"), 
        outName = paste0("HAP1.exonDel_", paste(rev(contrH[[i]]), collapse=":"))
    )
    saveGuidePairs(
        contrH[[i]], de = deH, gFitness = gFitnessH[[i]], region = "exons",
        outDir = file.path(outDir, "tables_HAP1.exonDel"), 
        outName = paste0("HAP1.exonDel_", paste(rev(contrH[[i]]), collapse=":"))
    )
}

for (i in 1:length(contrR)) {
    saveTables(
        exonDropR[[i]], pairFracCO=dropPairFracR[i], gFitness = gFitnessR[[i]],
        outDir = file.path(outDir, "tables_RPE1.exonDel"), 
        outName = paste0("RPE1.exonDel_", paste(rev(contrR[[i]]), collapse=":"))
    )
    saveGuidePairs(
        contrR[[i]], de = deR, gFitness = gFitnessR[[i]], region="exons",
        outDir = file.path(outDir, "tables_RPE1.exonDel"), 
        outName = paste0("RPE1.exonDel_", paste(rev(contrR[[i]]), collapse=":"))
    )

}


### Use the treshold from exon deletions to call hits on intron deletions

intronDropH <- mcmapply(
    FUN = getRegionDropFrac, contrH, gFitness = gFitnessH, SIMPLIFY = FALSE, 
    MoreArgs = list(
        de = deH, info = info, region = "introns", direction = "drop",
        direction="decrease", ii.quant = ii.quant, fit.quant = fit.quant,
        goodPairsOnly = goodPairsOnly, minInfPairs = minInfPairs, minChangePairs = minChangePairs
    ), 
    mc.cores = cores
)

intronDropR <- mcmapply(
    FUN = getRegionDropFrac, contrR, gFitness = gFitnessR, SIMPLIFY = FALSE, 
    MoreArgs = list(
        de = deR, info = info, region = "introns", direction = "drop",
        direction="decrease", ii.quant = ii.quant, fit.quant = fit.quant,
        goodPairsOnly = goodPairsOnly, minInfPairs = minInfPairs, minChangePairs = minChangePairs
    ), 
    mc.cores = cores
)


### Plot performance on intron deletion

pdf(file.path(outDir, "ScreenPerformance_HAP1.intronDel.pdf"), width = 9 , height = 5)
for (i in 1:length(contrH)) {
    plotRegionDropFrac(intronDropH[[i]], fracOverride = dropPairFracH[i])
    legend("topright", paste(rev(contrH[[i]]), collapse=":"), bty = "n", text.col = "grey60")
}
dev.off()

pdf(file.path(outDir, "ScreenPerformance_RPE1.intronDel.pdf"), width = 9 , height = 5)
for (i in 1:length(contrR)) {
    plotRegionDropFrac(intronDropR[[i]], fracOverride = dropPairFracR[i])
    legend("topright", paste(rev(contrR[[i]]), collapse=":"), bty="n", text.col="grey60")
}
dev.off()


### Save intron deletion results

for (i in 1:length(contrH)) {
    saveTables(
        intronDropH[[i]], pairFracCO = dropPairFracH[i], gFitness = gFitnessH[[i]], 
        outDir = file.path(outDir, "tables_HAP1.intronDel"), 
        outName = paste0("HAP1.intronDel_", paste(rev(contrH[[i]]), collapse = ":"))
    )
    saveGuidePairs(
        contrH[[i]], de = deH, gFitness = gFitnessH[[i]], region="introns",
        outDir = file.path(outDir, "tables_HAP1.intronDel"), 
        outName = paste0("HAP1.intronDel_", paste(rev(contrH[[i]]), collapse = ":"))
    )
}

for (i in 1:length(contrR)) {
    saveTables(
        intronDropR[[i]], pairFracCO = dropPairFracR[i], gFitness = gFitnessR[[i]],
        outDir = file.path(outDir, "tables_RPE1.intronDel"), 
        outName = paste0("RPE1.intronDel_", paste(rev(contrR[[i]]), collapse=":"))
    )
    saveGuidePairs(
        contrR[[i]], de = deR, gFitness = gFitnessR[[i]], region = "introns",
        outDir = file.path(outDir, "tables_RPE1.intronDel"),
        outName = paste0("RPE1.intronDel_", paste(rev(contrR[[i]]), collapse = ":")))
}

