### Obtain hgRNA log2-fold change (LFC) using edgeR, based on raw read and library size normalized counts

library(edgeR)



### Set parameters

inDir  <- "./input"
outDir <- "./output/edgeR"

minCounts <- 20


### Function definitions

mergeToAll <- function(x, allgenes=cts) {
### Merge a reduced table x to a larger one based on the row names in allgenes
    tmp <- merge(data.frame(ID = rownames(allgenes), geneInd = 1:nrow(allgenes)),
                 x,
                 by.x=1, by.y=0, all.x=T)
    tmp <- tmp[order(tmp$geneInd), -2]
    names(tmp)[2:3] <- c("log2FC","log2CPM")
    tmp
}

de.table <- function(de, contrasts) {
### Generate a table with DE for all contrasts
    out <- data.frame(Guide.ID = de[[1]]$ID)
    for (i in 1:length(contrasts)) {
        tmp <- signif(de[[i]][, c("log2FC","PValue","FDR")], 4)
        names(tmp) <- paste(contrasts[[i]][2], ":", contrasts[[i]][1], ".", names(tmp), sep="")
        out <- cbind(out, tmp)
    }
    out
}



### Read input

if (!dir.exists(outDir))  {dir.create(outDir, recursive = TRUE)}


raw   <- read.delim(file.path(inDir, "chyint2_screen2_merged_counts_11Oct22.txt"),
                     stringsAsFactors = T)

norm   <- read.delim(file.path(inDir, "chyint2_screen2_merged_norm.counts_11Oct22.txt"),
                     stringsAsFactors = T)


### Remove data we did not use
raw  <- raw[,c(1:16, 23:25, 32:37, 44:46, 53:55)]
norm <- norm[,c(1:16, 23:25, 32:37, 44:46, 53:55)]

info <- norm[,1:15]
raw  <- as.matrix(raw[,16:ncol(raw)])
norm <- as.matrix(norm[,16:ncol(norm)])
rownames(raw) <- rownames(norm) <- info$Guide.ID


### Filter out low abundance hgRNAs and perform edgeR analysis to get LFC

samples <- data.frame(Sample = colnames(raw))
samples$Type <- sub("_[ABC]$", "", samples$Sample)
samples$Treatment <- sub("(HAP1|RPE1)_", "", samples$Type)
sampH <- samples[c(8,9,1:7),]
sampR <- samples[c(8,9,10:16),]
groupsH <- relevel(as.factor(sampH$Treatment), ref = "T0")
groupsR <- relevel(as.factor(sampR$Treatment), ref = "T0")

rawH  <- raw[,sampH$Sample]
rawR  <- raw[,sampR$Sample]
normH <- norm[,sampH$Sample]
normR <- norm[,sampR$Sample]

# Determine what RPM corresponds to a minium count of minCounts in the smalles sample
minRPM.H <- 1000000 * minCounts / min(colSums(raw)[c("HAP1_T0","Library_pool_Lig2")])
minRPM.R <- 1000000 * minCounts / min(colSums(raw)[c("RPE1_T0","Library_pool_Lig2")])
cat("Minimum RPM in HAP1:", minRPM.H, "\n")
cat("Minimum RPM in RPE1:", minRPM.R, "\n")

# Filter out hgRNAs with low RPM in ligation or T0 samples
keepH <- normH[,"Library_pool_Lig2"] >= minRPM.H & normH[,"HAP1_T0"] >= minRPM.H  # Ligation 2 is smaller library
keepR <- normR[,"Library_pool_Lig2"] >= minRPM.R & normR[,"RPE1_T0"] >= minRPM.R
cat("Discarding", length(which(!keepH)), "of", nrow(rawH), "guide pairs in HAP1\n")
cat("Discarding", length(which(!keepR)), "of", nrow(rawR), "guide pairs in RPE1\n")

yH <- DGEList(counts = rawH, group=groupsH)
yR <- DGEList(counts = rawR, group=groupsR)

yH <- calcNormFactors(yH, method = "TMMwsp")
yR <- calcNormFactors(yR, method = "TMMwsp")

yH.all <- yH
yR.all <- yR

yH <- yH[keepH,]
yR <- yR[keepR,]

# Estimate dispersion and test for dropout/increase
yH <- estimateDisp(yH)
yR <- estimateDisp(yR)

pdf(file.path(outDir, "edgeR_dispersionEstimation.pdf"), width=5, height=8)
par(mfrow=c(2,1))
plotBCV(yH, main="edgeR dispersion estimation, HAP1", ylim=c(0.19, 2.9), xlim=c(-2, 18))
plotBCV(yR, main="edgeR dispersion estimation, RPE1", ylim=c(0.19, 2.9), xlim=c(-2, 18))
dev.off()


contrH <- list(
    c("T0","T12_NT"),
    c("T0","T18_NT")
)
contrR <- list(
    c("T0","T21_NT"),
    c("T0","T27_NT")
)

# edgeR exact test
etH <- lapply(contrH, FUN=function(x) {exactTest(yH, pair=x)})
etR <- lapply(contrR, FUN=function(x) {exactTest(yR, pair=x)})

deH <- lapply(etH, topTags, n=length(which(keepH)))
deR <- lapply(etR, topTags, n=length(which(keepR)))

deH <- lapply(deH, mergeToAll, allgenes=raw)
deR <- lapply(deR, mergeToAll, allgenes=raw)

# edgeR CPM
#rpmH <- cpmByGroup(yH.all)
#rpmR <- cpmByGroup(yR.all)

#colnames(rpmH) <- paste0(colnames(rpmH), ".RPM")
#colnames(rpmR) <- paste0(colnames(rpmR), ".RPM")
    
#rpmH.sing <- cpm(yH.all)
#rpmR.sing <- cpm(yR.all)
#colnames(rpmH.sing) <- paste0(colnames(rpmH.sing), ".RPM")
#colnames(rpmR.sing) <- paste0(colnames(rpmR.sing), ".RPM")

deH.all <- de.table(deH, contrH)
deR.all <- de.table(deR, contrR)
deH.all<- data.frame(Guide.ID = deH.all$Guide.ID, round(rpmH, 2), deH.all[,-1], check.names=F)
deR.all<- data.frame(Guide.ID = deR.all$Guide.ID, round(rpmR, 2), deR.all[,-1], check.names=F)

#gz <- gzfile(file.path(outDir, "RPM.samples.HAP1_UB221019.tab.gz"), "w")
#write.table(data.frame(Guide.ID = rownames(rpmH.sing), round(rpmH.sing, 2)),
#            row.names=F, col.names=T, quote=F, sep="\t", file = gz)
#close(gz)

#gz <- gzfile(file.path(outDir, "RPM.samples.RPE1_UB221019.tab.gz"), "w")
#write.table(data.frame(Guide.ID = rownames(rpmR.sing), round(rpmR.sing, 2)),
#            row.names=F, col.names=T, quote=F, sep="\t", file = gz)
#close(gz)

gz <- gzfile(file.path(outDir, "edgeR_HAP1.tab.gz"), "w")
write.table(deH.all, row.names=F, col.names=T, quote=F, sep="\t", file = gz)
close(gz)

gz <- gzfile(file.path(outDir, "edgeR_RPE1.tab.gz"), "w")
write.table(deR.all, row.names=F, col.names=T, quote=F, sep="\t", file = gz)
close(gz)

gz <- gzfile(file.path(outDir, "GuideInfo.tab.gz"), "w")
write.table(info, row.names=F, col.names=T, quote=F, sep="\t", file = gz)
close(gz)
