### Functions for scoring growth phenotypes (dropout or increase) of targeted regions
### in Int-CHyMErA screen.


getGeneFitness <- function(x, de) {
    ### Get mean gene LFC from exon-targeting Cas9 guides (TKOv3 and additional guides)
    ### x:             Length 2 vector specifying contrast
    ### de:            Data frame produced during edgeR analysis, with sample RPM, LFC and FDR
    i <- which(names(de) == paste0(paste(rev(x), collapse=":"), ".log2FC"))
    geneDrop <- aggregate(de[info$Data.Subset %in% c("TKOv3","Cas9 Exonic Ctrl"), i], 
                          by=list(gene = info$Target.gene[info$Data.Subset %in% c("TKOv3","Cas9 Exonic Ctrl")]), 
                          FUN=mean, na.rm=T)
    names(geneDrop)[2] <- "meanLFC"
    geneDrop
}

getRegionDropFrac <- function(x, de, gFitness, info, ii.quant = 0.05, fit.quant = c(0.2, 0.4),
                              region = c("exons","introns")[1], direction = c("drop","increase")[1],
                              goodPairsOnly = TRUE, minInfPairs = 2, minChangePairs = 2, co = seq(0, 1, 0.01)) {
    ### Calculate the fraction of changing hgRNAs among the informative hgRNAs for each target region
    ### Value:         List with the fraction of changing hgRNAs per region and other information
    ### x:             Slot of contrast list
    ### de:            deH/deR
    ### ii.quant:      Quantile of the distribution of log2FC of the intergenic-intergenic controls below which a guide pair
    ###                is considered to drop out
    ### fit.quant:     Quantiles of the gene fitness distribution obtained from exonic Cas9 guides that delineate
    ###                fitness and non-fitness genes (numeric vector length 2)
    ### region:        "exons" or "introns"; used to select guides and controls
    ### goodPairsOnly: Only consider guide pairs for dropout where none of the controls paired with intergenic guides
    ###                drop out.
    ### minInfPairs:   Minimum number of informative deletion guide pairs per region, i.e. non-NA if all are considered,
    ###                or else 'good' pairs if goodPairsOnly=TRUE.
    ### minChangePairs:  Minimum number of deletion guide pairs with dropout for the region to be scored a hit.
    ### co:            Vector of thresholds for the fraction of guide pairs with dropout at which to calculate.
    
    if (!(region %in% c("exons","introns")))    {stop("region must be either 'exons' or 'introns")}
    if (!(direction %in% c("drop","increase"))) {stop("direction must be either 'drop' or 'increase'")}
    if (!all(de$Guide.ID == info$Guide.ID))     {stop("Guide pairs in de and info files do not match")}
    
    ## To facilitate obtaining the 'effective deletion', extract guide coordinates
    gcoord <- data.frame(chrom        = sub("([^:]+):.*", "\\1", info$Cas9.Target.Site),
                         Cas9.start   = as.integer(sub("[^:]+:([0-9]+)-.*", "\\1", info$Cas9.Target.Site)),
                         Cas9.end     = as.integer(sub("[^:]+:[0-9]+-([0-9]+)_.*", "\\1", info$Cas9.Target.Site)),
                         Cas12a.start = as.integer(sub("[^:]+:([0-9]+)-.*", "\\1", info$Cpf1.Target.Site)),
                         Cas12a.end   = as.integer(sub("[^:]+:[0-9]+-([0-9]+)_.*", "\\1", info$Cpf1.Target.Site))
    )
    
    ## Decide which genes have fitness effect
    i <- which(names(de) == paste0(paste(rev(x), collapse=":"), ".log2FC"))
    thres <- quantile(de[info$Data.Subset == "II_Ctrl", i], p=ii.quant, na.rm=T)
    hit.thres <- quantile(gFitness$meanLFC, p=fit.quant, na.rm=T)
    if (direction == "drop") {
        genes.fit    <- unique(gFitness$gene[gFitness$meanLFC < hit.thres[1]])
        genes.nonfit <- unique(gFitness$gene[gFitness$meanLFC > hit.thres[2]])
    }
    if (direction == "increase") {
        genes.fit    <- unique(gFitness$gene[gFitness$meanLFC > hit.thres[1]])
        genes.nonfit <- unique(gFitness$gene[gFitness$meanLFC < hit.thres[2]])
    }
    
    
    ## Check which double-deletion guides have controls that change
    if (region == "exons") {
        uqEvent <- unique(info$Event[info$Data.Subset %in% c("Exon Deletion", "ExonDeletion")])
    }
    if (region == "introns") {
        uqEvent <- unique(info$Event[info$Data.Subset %in% c("Intron Deletion")])
    }
    uqEvent <- data.frame(Event   = uqEvent,
                          Gene    = info$Target.gene[match(uqEvent, info$Event)]
    )
    if (direction == "drop") {
        badCas9  <- unique(info$Cas9.Guide[which(info$Cpf1.Guide.Type == "Intergenic" & de[,i] < thres)])
        badCas12 <- unique(info$Cpf1.Guide[which(info$Cas9.Guide.Type == "Intergenic" & de[,i] < thres)])
    }
    if (direction == "increase") {
        badCas9  <- unique(info$Cas9.Guide[which(info$Cpf1.Guide.Type == "Intergenic" & de[,i] > thres)])
        badCas12 <- unique(info$Cpf1.Guide[which(info$Cas9.Guide.Type == "Intergenic" & de[,i] > thres)])
    }
    notBadPairs <- 
        !(info$Cas9.Guide %in% badCas9) &
        !(info$Cpf1.Guide %in% badCas12)
    if (direction == "drop") {
        pairChange <- de[,i] < thres
    }
    if (direction == "increase") {
        pairChange <- de[,i] > thres
    }
    
    ## Get numbers of double deletion pairs for each unique locus
    if (region == "exons") {
        regionDelPairs <- info$Data.Subset %in% c("Exon Deletion","ExonDeletion")
    }
    if (region == "introns") {
        regionDelPairs <- info$Data.Subset %in% c("Intron Deletion")
    }
    pairCounts.double <- t(sapply(uqEvent$Event, FUN=function(y) {
        c(totPairs      = length(which(regionDelPairs & 
                                           info$Event == y &
                                           !is.na(de[,i]))),
          okPairs       = length(which(regionDelPairs & 
                                           info$Event == y &
                                           notBadPairs &
                                           !is.na(de[,i]))),
          totPairsDrop  = length(which(regionDelPairs &
                                           info$Event == y &
                                           pairChange)),
          okPairsDrop   = length(which(regionDelPairs &
                                           info$Event == y &
                                           notBadPairs  &
                                           pairChange))
        )
    }))
    
    ## Get effective deletion for each locus, i.e. the maximum region spanned by sign. guides
    sel <- regionDelPairs & pairChange 
    if (goodPairsOnly) {sel <- sel & notBadPairs}
    sel <- which(sel)
    eff.min <- aggregate(c(gcoord$Cas9.start[sel], gcoord$Cas12a.start[sel]),
                         by=list(Event = rep(info$Event[sel], 2)), min)
    eff.max <- aggregate(c(gcoord$Cas9.start[sel], gcoord$Cas12a.start[sel]),
                         by=list(Event = rep(info$Event[sel], 2)), max)
    eff.coord <- paste0(gcoord$chrom[match(eff.min$Event, info$Event)], ":",
                        eff.min$x, "-", eff.max$x)
    eff.coord <- data.frame(Event = uqEvent$Event,
                            effectiveDeletion = eff.coord[match(uqEvent$Event, eff.min$Event)]
                            )
    
    ## Get numbers of single cut control pairs for each unique locus
    if (region == "exons") {
        regionCtlPairs <- info$Data.Subset %in% c("Exon Deletion - Cas9Ctrl","Exon Deletion - Cas12aCtrl","ExonDeletion_Ctrl")
    }
    if (region == "introns") {
        regionCtlPairs <- info$Data.Subset %in% c("Intron Deletion - Cas9Ctrl","Intron Deletion - Cas12aCtrl")
    }
    pairCounts.single <-t(sapply(uqEvent$Event, FUN=function(y) {
        c(totPairs      = length(which(regionCtlPairs & 
                                           info$Event == y &
                                           !is.na(de[,i]))),
          okPairs       = length(which(regionCtlPairs & 
                                           info$Event == y &
                                           notBadPairs &
                                           !is.na(de[,i]))),
          totPairsDrop  = length(which(regionCtlPairs &
                                           info$Event == y &
                                           pairChange)),
          okPairsDrop   = length(which(regionCtlPairs &
                                           info$Event == y &
                                           notBadPairs  &
                                           pairChange))
        )
    }))
    perc.double <- pairCounts.double[,3:4] / pairCounts.double[,1:2]
    perc.single <- pairCounts.single[,3:4] / pairCounts.single[,1:2]
    
    ## Compare loci in fitness vs. non-fitness genes
    use <- pairCounts.double[,1 + goodPairsOnly] >= minInfPairs
    regionDel.fit    <- use & uqEvent$Gene %in% genes.fit
    regionDel.nonfit <- use & uqEvent$Gene %in% genes.nonfit
    
    percRegions <- matrix(nrow=length(co), ncol=4,
                          dimnames=list(co, c("DoubleX fitness","DoubleX nonfitness",
                                              "SingleX fitness","SingleX nonfitness"))
    )
    for (j in 1:length(co)) {
        percRegions[j,"DoubleX fitness"]    <- length(which(regionDel.fit & 
                                                                perc.double[,1 + goodPairsOnly] > co[j] &
                                                                pairCounts.double[,3 + goodPairsOnly] >= minChangePairs))
        percRegions[j,"DoubleX nonfitness"] <- length(which(regionDel.nonfit &
                                                                perc.double[,1 + goodPairsOnly] > co[j] &
                                                                pairCounts.double[,3 + goodPairsOnly] >= minChangePairs))
        percRegions[j,"SingleX fitness"]    <- length(which(regionDel.fit &
                                                                perc.single > co[j] &
                                                                pairCounts.single[,3 + goodPairsOnly] >= minChangePairs))
        percRegions[j,"SingleX nonfitness"] <- length(which(regionDel.nonfit &
                                                                perc.single > co[j] &
                                                                pairCounts.single[,3 + goodPairsOnly] >= minChangePairs))
    }
    percRegions[,c(1,3)] <- 100 * percRegions[,c(1,3)]/length(which(regionDel.fit))
    percRegions[,c(2,4)] <- 100 * percRegions[,c(2,4)]/length(which(regionDel.nonfit))
    if (length(which(regionDel.fit))    == 0) {percRegions[,c(1,3)] <- 0}
    if (length(which(regionDel.nonfit)) == 0) {percRegions[,c(2,4)] <- 0}
    
    list(region            = region,
         direction         = direction,
         contrast          = paste(rev(x), collapse=":"),
         ii.quant          = ii.quant,           # quantile of interngenic-intergenic controls for threshold
         ii.thres          = thres,              # value of this threshold
         fit.quant         = fit.quant,          # quantiles for fitness genes and non-fitness genes
         goodPairsOnly     = goodPairsOnly,
         minInfPairs       = minInfPairs,
         minChangePairs    = minChangePairs,
         pairCounts.double = pairCounts.double,
         pairCounts.single = pairCounts.single,
         effDeletion       = eff.coord,          
         regionDel.fit     = regionDel.fit,
         regionDel.nonfit  = regionDel.nonfit,
         perc.double       = perc.double,
         perc.single       = perc.single,
         percRegions       = percRegions
    )   
}


plotRegionDropFrac <- function(x, xlim = NA, ylim = NA, cols = c("brown1","dodgerblue1"), fracOverride = NA) {
    ### Plot the fraction of hit regions across a range of cutoffs on the fraction of changing hgRNAs
    ### and the same plot for hits if the intronic-intergenic controls were used instead.
    ### Value:        Threshold with the best separation between regions in fitness- and non-fitness genes
    ### x:            Contrast (vector of length 2), a slot of contrH/contrR
    ### fracOverride: Rather than find the optimal cutoff, use this threshold for the fraction of pairs
    ###               with dropouts to call a region a hit
    
    maxDelta <- which(x$percRegions[,"DoubleX fitness"] - x$percRegions[,"DoubleX nonfitness"] == 
                          max(x$percRegions[,"DoubleX fitness"] - x$percRegions[,"DoubleX nonfitness"], na.rm=T))
    
    co <- as.numeric(rownames(x$percRegions))
    if (any(is.na(xlim))) {xlim <- range(co)}
    if (any(is.na(ylim))) {ylim <- c(0, max(x$percRegions))}
    
    if (is.na(fracOverride)) {
        hitCO <- co[min(maxDelta)]
        hitInd <- min(maxDelta)
    } else {
        hitCO <- fracOverride
        hitInd <- max(which(co <= fracOverride))
    }
    
    if (x$region == "exons") {
        label.lower <- "exon"
        label.upper <- "Exon"
    }
    if (x$region == "introns") {
        label.lower <- "intron region"
        label.upper <- "Intron region"
    }
    
    if (x$direction == "drop") {
        dir.lower <- "dropout"
        dir.upper <- "Dropout"
    }
    if (x$direction == "increase") {
        dir.lower <- "increase"
        dir.upper <- "Increase"
    }
    
    par(mfrow=c(1,2))
    plot(0, 0, type="n", xlim=xlim, ylim=ylim,
         xlab=paste0("Cutoff for calling ", label.lower, " phenotype\n(fraction of significant pairs)"),
         ylab=paste0("% ", label.upper, "s with phenotype"),
         main=paste0(dir.upper, " upon ", label.lower, " deletion\n(intronic-intronic pairs)"))
    lines(c(rep(co[-length(co)], each=2), co[length(co)]), 
          c(x$percRegions[1,"DoubleX nonfitness"], rep(x$percRegions[-1,"DoubleX nonfitness"], each=2)),
          lwd=2, col=cols[2]
    )
    lines(c(rep(co[-length(co)], each=2), co[length(co)]), 
          c(x$percRegions[1,"DoubleX fitness"], rep(x$percRegions[-1,"DoubleX fitness"], each=2)),
          lwd=2, col=cols[1]
    )
    abline(h=0)
    abline(v=hitCO, lty=2)
    text(hitCO, 0.98 * ylim[2], round(hitCO, 2), pos=4, cex=0.8)
    text(hitCO, x$percRegions[hitInd,"DoubleX nonfitness"], 
         paste0(round(x$percRegions[hitInd,"DoubleX nonfitness"], 1), "%"), pos=4, cex=0.8, col=cols[2])
    text(hitCO, x$percRegions[hitInd,"DoubleX fitness"], 
         paste0(round(x$percRegions[hitInd,"DoubleX fitness"], 1), "%"), pos=4, cex=0.8, col=cols[1])
    legend("topright", bty="n", lwd=2, col=cols, cex=0.8,
           paste0(label.upper, c("s in fitness genes","s in non-fitness genes")))
    legend ("right", cex=0.6, bty="n",
            c(ifelse(x$goodPairsOnly, "Only 'clean' guide pairs", "All guide pairs"),
              paste0("Only ", label.lower, "s with min. ", x$minInfPairs, " informative guide pairs"),
              paste0("Requiring min. ", x$minChangePairs, " pairs with ", dir.lower)
            ))
    
    plot(0, 0, type="n", xlim=xlim, ylim=ylim,
         xlab=paste0("Cutoff for calling ", label.lower, " phenotype\n(fraction of significant pairs)"),
         ylab=paste0("% ", label.upper, "s with phenotype"),
         main=paste0(dir.upper, " upon single cut\n(intronic-intergenic pairs)"))
    lines(c(rep(co[-length(co)], each=2), co[length(co)]), 
          c(x$percRegions[1,"SingleX nonfitness"], rep(x$percRegions[-1,"SingleX nonfitness"], each=2)),
          lwd=2, col=cols[2]
    )
    lines(c(rep(co[-length(co)], each=2), co[length(co)]), 
          c(x$percRegions[1,"SingleX fitness"], rep(x$percRegions[-1,"SingleX fitness"], each=2)),
          lwd=2, col=cols[1]
    )
    abline(h=0)
    abline(v=hitCO, lty=2)
    text(hitCO, 0.98 * ylim[2], round(hitCO, 2), pos=4, cex=0.8)
    text(hitCO, x$percRegions[hitInd,"SingleX nonfitness"], 
         paste0(round(x$percRegions[hitInd,"SingleX nonfitness"], 1), "%"), pos=4, cex=0.8, col=cols[2])
    text(hitCO, x$percRegions[hitInd,"SingleX fitness"], 
         paste0(round(x$percRegions[hitInd,"SingleX fitness"], 1), "%"), pos=4, cex=0.8, col=cols[1])
    legend ("right", cex=0.6, bty="n",
            c(paste0("Only ", label.lower, "s with min. ", x$minInfPairs, " informative guide pairs"),
              paste0("Requiring min. ", x$minChangePairs, " pairs with ", dir.lower)
            ))
    
    hitCO
}

saveTables <- function(x, pairFracCO, gFitness, outDir, outName = "Output") {
    ### Save target region result tables and parameter settings
    ### x:          A slot of the list structure produced by getRegionDropFrac()
    ### pairFracCO: Fraction of significant pairs for calling a region a hit
    ### gFitness:   Data frame with gene fitness information (LFC)
    ### outDir:     Directory to save the tables
    ### outName:    Prefix for the output files

    dirDrop <- x$direction == "drop"
    
    tmp1 <- x$pairCounts.double
    tmp2 <- x$pairCounts.single
    colnames(tmp1) <- paste0(colnames(tmp1), ".doubleX")
    colnames(tmp2) <- paste0(colnames(tmp2), ".singleX")

    tmp3 <- round(x$perc.double, 4)
    tmp4 <- round(x$perc.single, 4)
    colnames(tmp3) <- c("fractionTotPairs.doubleX","fractionOkPairs.doubleX")
    colnames(tmp4) <- c("fractionTotPairs.singleX","fractionOkPairs.singleX")
    tmp <- data.frame(info[match(rownames(tmp1), info$Event),c(5,4,2,6,7)], 
                      hit                = NA,
                      effectiveDeletion  = NA,
                      geneFitness        = NA,
                      geneLFC.exonicCas9 = NA,
                      tmp3, tmp4, tmp1, tmp2)
    tmp$hit <- x$perc.double[, 1 + x$goodPairsOnly] >= pairFracCO & 
        x$pairCounts.double[,3 + x$goodPairsOnly] >= x$minChangePairs
    tmp$hit[x$pairCounts.double[,1 + x$goodPairsOnly] < x$minInfPairs] <- NA
    tmp$effectiveDeletion <- x$effDeletion$effectiveDeletion
    tmp$effectiveDeletion[!tmp$hit | is.na(tmp$hit)] <- NA
    tmp$geneFitness <- ifelse(x$regionDel.fit, "fitness", "")
    tmp$geneFitness[x$regionDel.nonfit] <- "non-fitness"
    tmp$geneLFC.exonicCas9 <- round(gFitness$meanLFC[match(tmp$Target.gene, gFitness$gene)], 4)

    write.csv(tmp, row.names=F,
              file=file.path(outDir, paste0(outName, "_", 
                                            ifelse(dirDrop, "Dropout", "Increase"), 
                                            "_Counts.csv")))
    
    write.csv(data.frame(FractionCutoff = rownames(x$percRegions), round(x$percRegions, 2), check.names=F), 
              row.names=F,
              file=file.path(outDir, paste0(outName, "_Performance.csv")))
    
    tmp <- data.frame(Parameter = names(x)[1:9], 
                      Value     = sapply(x[1:9], "[[", 1))
    tmp$Value[6] <- paste0("fitness:", x[[6]][1], ", non-fitness:", x[[6]][2])
    tmp <- rbind(tmp, data.frame(Parameter = "fit.thres",
                                 Value     = paste0("fitness",
                                                    ifelse(dirDrop, "<", ">"),
                                                    round(quantile(gFitness$meanLFC, p=x[[6]][1], na.rm=T),4), 
                                                    ", non-fitness", 
                                                    ifelse(dirDrop, ">", "<"),
                                                    round(quantile(gFitness$meanLFC, p=x[[6]][2], na.rm=T),4))
    )
    )
    tmp <- rbind(tmp, data.frame(Parameter = "changeFracForHit",
                                 Value     = pairFracCO
    )
    )
    
    write.csv(tmp[c(2,1,3:5,9,10,6:8),], row.names=F,
              file=file.path(outDir, paste0(outName, "_Parameters.csv")))
    
}

saveGuidePairs <- function(x, de, gFitness, ii.quant = 0.05, fit.quant = c(0.2, 0.4), region = c("exons","introns")[1],
                           direction = c("drop","increase")[1], outDir, outName = "Output") {
    ### Save tables of deletion guide pairs, intergenic-intergenic controls, and exonic Cas9 deletions
    ### x:          A slot of the list structure produced by getRegionDropFrac()
    ### de:         Data frame produced during edgeR analysis, with sample RPM, LFC and FDR
    ### gFitness:   Data frame with gene fitness information (LFC)
    ### ii.quant:   Quantile of the distribution of log2FC of the intergenic-intergenic controls below which (if dropout) 
    ###             a guide pair is considered to change significantly
    ### fit.quant:  Quantiles of the gene fitness distribution obtained from exonic Cas9 guides that delineate
    ###             fitness and non-fitness genes (numeric vector length 2)
    ### region:     "exons" or "introns"; used to select guides and controls
    ### direction:  "drop" or "increase"; depending on this ii.quant is a lower or upper bound
    ###             and fit.quant is c(lower bound of fitness genes, upper bound of non-fitness genes)
    ###             c(upper bound of fitness genes, lower bound of non-fitness genes)
    ### outDir:     Directory to save the tables
    ### outName:    Prefix for the output files

    
    ## Get threshold for dropout
    i <- which(names(de) == paste0(paste(rev(x), collapse=":"), ".log2FC"))
    thres <- quantile(de[info$Data.Subset == "II_Ctrl", i], p=ii.quant, na.rm=T)
    
    ### Deletion guide pair table
    ## Check which double-deletion guides have controls that change
    if (region == "exons") {
        regionDelPairs <- info$Data.Subset %in% c("Exon Deletion","ExonDeletion")
        regionCtlPairs <- info$Data.Subset %in% c("Exon Deletion - Cas9Ctrl","Exon Deletion - Cas12aCtrl","ExonDeletion_Ctrl")
    }
    if (region == "introns") {
        regionDelPairs <- info$Data.Subset %in% c("Intron Deletion")
        regionCtlPairs <- info$Data.Subset %in% c("Intron Deletion - Cas9Ctrl","Intron Deletion - Cas12aCtrl")
    }
    
    if (direction == "drop") {
        badCas9  <- unique(info$Cas9.Guide[which(info$Cpf1.Guide.Type == "Intergenic" & de[,i] < thres)])
        badCas12 <- unique(info$Cpf1.Guide[which(info$Cas9.Guide.Type == "Intergenic" & de[,i] < thres)])
    }
    if (direction == "increase") {
        badCas9  <- unique(info$Cas9.Guide[which(info$Cpf1.Guide.Type == "Intergenic" & de[,i] > thres)])
        badCas12 <- unique(info$Cpf1.Guide[which(info$Cas9.Guide.Type == "Intergenic" & de[,i] > thres)])
    }
    notBadPairs <- 
        !(info$Cas9.Guide %in% badCas9) &
        !(info$Cpf1.Guide %in% badCas12)
    
    ## Make list of double-deletion guide pairs
    if (direction == "drop") {
        minLFC.cas9  <- suppressWarnings(aggregate(de[regionCtlPairs & info$Cas9.Guide.Type == "Intronic",i], min, na.rm=T, 
                                                   by=list(Cas9.Guide=info$Cas9.Guide[regionCtlPairs & 
                                                                                          info$Cas9.Guide.Type == "Intronic"])))
        minLFC.cas12 <- suppressWarnings(aggregate(de[regionCtlPairs & info$Cpf1.Guide.Type == "Intronic",i], min, na.rm=T, 
                                                   by=list(Cpf1.Guide=info$Cpf1.Guide[regionCtlPairs & 
                                                                                          info$Cpf1.Guide.Type == "Intronic"])))
        minLFC.cas9$x[is.infinite(minLFC.cas9$x)] <- NA
        minLFC.cas12$x[is.infinite(minLFC.cas12$x)] <- NA
        
        delPairTab <- info[regionDelPairs,]
        delPairTab$drop <- de[regionDelPairs, i] < thres
        delPairTab$log2FC <- de[regionDelPairs, i]
        delPairTab$okPair <- notBadPairs[regionDelPairs]
        delPairTab$worstCas9Ctrl   <- minLFC.cas9$x[match(delPairTab$Cas9.Guide, minLFC.cas9$Cas9.Guide)]
        delPairTab$worstCas12aCtrl <- minLFC.cas12$x[match(delPairTab$Cpf1.Guide, minLFC.cas12$Cpf1.Guide)]
    }
    if (direction == "increase") {
        maxLFC.cas9  <- suppressWarnings(aggregate(de[regionCtlPairs & info$Cas9.Guide.Type == "Intronic",i], max, na.rm=T, 
                                                   by=list(Cas9.Guide=info$Cas9.Guide[regionCtlPairs & 
                                                                                          info$Cas9.Guide.Type == "Intronic"])))
        maxLFC.cas12 <- suppressWarnings(aggregate(de[regionCtlPairs & info$Cpf1.Guide.Type == "Intronic",i], max, na.rm=T, 
                                                   by=list(Cpf1.Guide=info$Cpf1.Guide[regionCtlPairs & 
                                                                                          info$Cpf1.Guide.Type == "Intronic"])))
        maxLFC.cas9$x[is.infinite(maxLFC.cas9$x)] <- NA
        maxLFC.cas12$x[is.infinite(maxLFC.cas12$x)] <- NA
        
        delPairTab <- info[regionDelPairs,]
        delPairTab$increase <- de[regionDelPairs, i] > thres
        delPairTab$log2FC <- de[regionDelPairs, i]
        delPairTab$okPair <- notBadPairs[regionDelPairs]
        delPairTab$worstCas9Ctrl   <- maxLFC.cas9$x[match(delPairTab$Cas9.Guide, maxLFC.cas9$Cas9.Guide)]
        delPairTab$worstCas12aCtrl <- maxLFC.cas12$x[match(delPairTab$Cpf1.Guide, maxLFC.cas12$Cpf1.Guide)]
    }
    
    write.csv(delPairTab, row.names=F,
              file=file.path(outDir, paste0(outName, "_deletionGuidePairs.csv")))
    
    
    ### Cas9 exonic guide table
    sel <- info$Data.Subset %in% c("TKOv3", "Cas9 Exonic Ctrl")
    hit.thres <- quantile(gFitness$meanLFC, p=fit.quant, na.rm=T)
    if (direction == "drop") {
        genes.fit    <- unique(gFitness$gene[gFitness$meanLFC < hit.thres[1]])
        genes.nonfit <- unique(gFitness$gene[gFitness$meanLFC > hit.thres[2]])
    }
    if (direction == "increase") {
        genes.fit    <- unique(gFitness$gene[gFitness$meanLFC > hit.thres[1]])
        genes.nonfit <- unique(gFitness$gene[gFitness$meanLFC < hit.thres[2]])
    }
    exonicCas9tab <- data.frame(info[sel,],
                                log2FC = de[sel, i],
                                geneFitness = NA
                                )
    exonicCas9tab$geneFitness[exonicCas9tab$Target.gene %in% genes.fit]    <- "fitness"
    exonicCas9tab$geneFitness[exonicCas9tab$Target.gene %in% genes.nonfit] <- "non-fitness"
    exonicCas9tab <- exonicCas9tab[order(exonicCas9tab$Target.gene,
                                         exonicCas9tab$log2FC * ifelse(direction == "drop", 1, -1)),]

    write.csv(exonicCas9tab, row.names=F,
              file=file.path(outDir, paste0(outName, "_exonicCas9.csv")))
    
    
    ### Cas9 exonic guide table
    sel <- info$Data.Subset == "II_Ctrl"
    iiCtrls <- data.frame(info[sel,],
                          log2FC = de[sel, i]
    )
    iiCtrls <- iiCtrls[order(iiCtrls$log2FC * ifelse(direction == "drop", -1, 1)), -c(2:5,7)]
    
    write.csv(iiCtrls, row.names=F,
              file=file.path(outDir, paste0(outName, "_IIcontrols.csv")))
}
