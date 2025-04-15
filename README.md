# A genome-wide functional analysis of conserved intronic regions reveals essential roles for speckle-associated detained introns
Farhangmehr, Braunschweig et al. (2025)

## Int-CHyMErA screen of intronic deletions

### Overview
2,600 intronic deletions were created using the ‘Cas Hybrid for Multiplexed Editing and screening Applications’ (CHyMErA) system [(Gonatopoulos-Pournatzis, Aregger et al., 2020)](https://pubmed.ncbi.nlm.nih.gov/32249828/) in HAP1 and RPE-1 cells, and fitness phenotypes were measured by assessing changes in the abundance of cells with a given deletion at time point T0 and T12/T18 (HAP1) or T21/T27 (RPE-1) through sequencing of hybrid guide RNAs (hgRNAs). Each hgRNA consists of a Cas9 and a Cas12a gRNA connected with a linker, that are processed by Cas12a into the mature gRNAs and mediate targeting of a Cas9 and a Cas12a endonuclease complex to sites flanking a region targeted for deletion.

This repository hosts the library design file and raw/normalized hgRNA counts, as well as the methods to calculate abundance changes from hgRNA counts and identify regions that confer decreased or increased cell fitness when deleted, i.e., screen hits. Please see the manuscript for details. Applying the steps below to the input read counts yields the files in the output. Raw sequencing data can be found on GEO ###########.

Time point T18 for HAP1 and T21 for RPE-1 were used for further analysis.


### Types of hgRNAs
- **Intronic deletion hgRNAs**: These are the main part of the screen, designed to crease intronic deletions either (i) deleting most of an entire intron, (ii) a conserved block within an intron, or (iii) a conserved region overlapping an alternative exon mediating nonsense-mediated mRNA decay. Up to 16 overlapping hgRNAs target each region.
- **Intergenic-intergenic controls:** Negative controls that should not result in deletion with a fitness effect. The distribution of log2-fold changes (log2FC) of these controls is used to call significant effects: Deletions with a log2FC lower than the 5th percentile are assumed to mediate a significant dropout, and above the 95th percentile are assumed to mediate a significant fitness increase.
- **Deletion specificity controls:** Each Cas9 and Cas12a part of an hgRNA desiged to crease a deletion is also present in combination with a partner targeting a distal intergenic regions. This control is expected to cause no or a weaker effect. Deletion hgRNAs where one of these controls does cause an effect are disregarded.
- **Exon-targeting controls**: hgRNAs where either the Cas9 or Cas12a gRNA targets a coding exon as in a classic CRISPR deletion screen. The large majority of exon-targeting Cas9 gRNAs are derived from the [Toronto CRISPR human knockout library](https://www.addgene.org/pooled-library/moffat-crispr-knockout-tkov3/) (TKOv3). The Cas9 gRNAs are used to measure fitness effect of gene depletion. 
- **Exon-deletion controls:** hgRNAs designed to delete a coding exon with the Cas9 and Cas12a gRNA targeted to intronic regions on either side of the exon. These controls are used to assess the efficacy of the screen, and to infer the optimal fraction of overlapping deletion hgRNAs associated with a fitness phenotyope to call a region a hit by maximizing the differential between exonic deletions in fitness- vs. non-fitness genes.


### Steps
**1. Normalize raw counts and derive log2FC using edgeR**
```
Rscript GetLFC.R
```
**2. Call dropout hits**
Generate lists of fitness (i.e., causing decreased cell fitness when deleted) and non-fitness genes using exon-targeting controls. Identify hgRNAs with a dropout; remove those where a deletion specificity control also has a fitness effect. Use exon-deletion controls to derive the optimal fraction of overlapping hgRNAs that must have a dropout to call a region a hit. Apply criteria to intron deletions and call dropout hits.
```
Rscript RegionScoring_dropout.R
```
**3. Call increase hits**
Similar to calling dropout hits, but re-use the optimized fraction of required overlapping hgRNAs with a fitness effect from the analysis of dropout-causing exon deletion controls, since they should have a similar deletion efficacy.
```
Rscript RegionScoring_increase.R
```


### Output files
**Folder _edgeR_**  
**edgeR_dispersionEstimation.pdf**  
Plots of common, trended and tag-wise dispersion from edgeR.

**edgeR_[cell].tab.gz**  
edgeR output including RPM, log2FC and FDR for every hgRNA.  

**Folders _dropout_ and _increase_**  
Results from scoring deletions causing decreased fitness  
**[cell].geneDropout.Cas9exonic_[time point].csv**  
Average log2FC per gene of coding exon-targeting Cas9 guides. Used to segment genes into fitness and non-fitness genes.

**ScreenPerformance_[cell].exonDel.pdf**  
Left, percentage of deletions of exons within fitness genes or non-fitness genes that cause dropout (y-axis), out of those scored as hits at a given threshold of the fraction of overlapping hgRNAs with a phenotype (x-axis). Dashed line indicates threshold with best separation, which is then applied to scoring intronic deletions. Right, same but for intronic-intergenic hgRNAs. Note that due to the settings in this analysis, none of these controls were allowed to have a fitness effect. Additional criteria are stated in the figure.

**ScreenPerformance_[cell].intronDel.pdf**  
Same as above, but for intronic deletions. Dashed line indicates threshold derived from exon deletions and applied to intronic deletions.

Within subfolders:  
**[cell].exonDel_[time point]_deletionGuidePairs.csv**  
**[cell].intronDel_[time point]_deletionGuidePairs.csv**  
hgRNA-level information. *Guide.ID*, unique hgRNA ID; *Event*, unique ID for deletion target; *Data.Subset*, type of hgRNA; *Event.Type, type of target; *Cas9.Guide.Type*, type of target of Cas9 part of hgRNA; *Doench.Score*, Cas9 score based on Doench et al., 2016; *Cpf1.Guide.Type*, type of target of Cas12a part of hgRNA; *Cpf1.CNN.Score*, [CHyMErA-Net](https://github.com/BlencoweLab/CHyMErA-Net) Cas12a guide score; *drop*/*increase*: hgRNA had significant dropout/increase; *log2FC* log1-fold change of hgRNA vs T0; *okPair*, TRUE if neither of the two component gRNAs had a fitness effect when combined with intergenic control; *worstCas9Ctrl*, lowest (for dropout) log2FC of Cas9 gRNA when combined with intergenic control Cas12a gRNA; *worstCas12aCtrl*  lowest (for dropout) log2FC of Cas12a gRNA when combined with intergenic control Cas9 gRNA.

**[cell].intronDel_[time point]_Dropout_Counts.csv**  
**[cell].intronDel_[time point]_Increase_Counts.csv**  
Event-level information for intronic deletions. *Event.Type*, type of target deletion (cons.intron, deletion of most of an intron; cons.region, deletion of conserved intronic region; AS, deletion of regions containing PTC exon; non-cons, deletion of non-conserved intron; Cas12a.exonic, positive Cas12a control); *Span*, width of target region; *hit*, target deletion was scored as a hit; *effectiveDeletion*, largest deletion in common across overlapping hgRNAs with fitness effect; *geneFitness*, fitness effect of host gene depletion from Cas9 exonic targeting; *geneLFC.exonicCas9*, mean log2FC of Cas9 gRNAs targeting coding exons; *totPairs.doubleX*, total number of deletion hgRNAs; *okPairs.doubleX*, number of deletion hgRNAs of which neither intronic-intergenic control showed a phenoytpe; *totPairsDrop.doubleX*, number of deletion hgRNAs with dropout; *okPairsDrop.doubleX*, number of deletion hgRNAs with dropout and of which neither intronic-intergenic control showed a fitness effect; *totPairs.singleX*, total number of intronic-intergenic hgRNAs; *okPairs.singleX*, number of intronic-intergenic hgRNAs that have no fitness effect; *totPairsDrop.singleX*, total number of intronic-intergenic hgRNAs with dropout; *okPairsDrop.singleX*, by definition 0 using current parameters; *fractionTotPairs.doubleX*, fraction of totPairs.doubleX with dropout; *fractionOKPairs.doubleX*, fraction of okPairs.doubleX with dropout; *fractionTotPairs.singleX*, fraction of totPairs.singleX with dropout; *fractionOKPairs.singleX*, fraction of okPairs.singleX with dropout (0 by definition using current parameters). Analogous for increase.

**[cell].intronDel_[time point]_exonicCas9.csv**  
**[cell].exonDel_[time point]_exonicCas9.csv**  
hgRNA-level information for coding exon-targeting Cas9 gRNAs. See above for column definitions. *geneFitness*, aggregate gene label based on mean of all coding exon-targeting Cas9 gRNAs per gene.

**[cell].intronDel_[time point]_IIcontrols.csv**  
**[cell].intronDel_[time point]_IIcontrols.csv**  
hgRNA-level information for intronic-intronic negative control hgRNAs.

**[cell].intronDel_[time point]_Parameters.csv**  
**[cell].exonDel_[time point]_Parameters.csv**  
User options of the analysis run and empirical thresholds. *direction*, drop or increase; *region*, exons or introns, *contrast*, time points compared; *ii.quant*, quantile of intergenic-intergenic controls below (dropout) or above (increase) which log2FC are considered a significant fitness effect; *ii.thres*, log2FC corresponding to ii.quant; *minChangePairs*, minimum hgRNAs with significant effect required for hit call; *fit.thres*, log2FC threshold corresponding to fit.quant; *fit.quant*, quantiles of gene fitness scores beyond which to consider genes fitness or noon-fitness genes; *goodPairsOnly*, whether to require that hgRNAs considered for hit calling must exclude those where intronic-intergenic controls result in a fitness effect; *minInfPairs*, minmum number of hgRNAs passing filtering for a targeted deletion to be evaluated.

**[cell].intronDel_[time point]_Performance.csv**  
**[cell].exonDel_[time point]_Performance.csv**  
Data underlying screen performance plots.



### Reference
Shaghayegh Farhangmehr*, Ulrich Braunschweig*, Mingkun Wu, Syed Nabeel-Shah1, Kevin Brown, Jason Moffat, and Benjamin J. Blencowe. A genome-wide functional analysis of conserved intronic regions reveals essential roles for speckle-associated detained introns. *Manuscript submitted* (2025)

### Contact
For questions or comments, please raise an issue or email u -dot- braunschweig -at- utoronto -dot- ca.