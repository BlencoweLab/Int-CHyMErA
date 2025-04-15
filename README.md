# A genome-wide functional analysis of conserved intronic regions reveals essential roles for speckle-associated detained introns
Farhangmehr, Braunschweig et al., 2025

## Int-CHyMErA screen of intronic deletions

### Overview
2,600 intronic deletions were created using the ‘Cas Hybrid for Multiplexed Editing and screening Applications’ (CHyMErA) system [(Gonatopoulos-Pournatzis, Aregger et al., 2020)](https://pubmed.ncbi.nlm.nih.gov/32249828/) in HAP1 and RPE-1 cells, and fitness phenotypes were measured by assessing changes in the abundance of cells with a given deletion at time point T0 and T18 (HAP1) or T21 (RPE-1) through sequencing of hybrid guide RNAs (hgRNAs). Each hgRNA consists of a Cas9 and a Cas12a gRNA connected with a linker, that are processed by Cas12a into the mature gRNAs and mediate targeting of a Cas9 and a Cas12a endonuclease complex to sites flanking a region targeted for deletion.

This repository hosts the library design file and raw/normalized hgRNA counts, as well as the methods to calculate abundance changes from hgRNA counts and identify regions that confer decreased or increased cell fitness when deleted, i.e., screen hits. Please see the manuscript for details. Applying the steps below to the input (raw and library-size normalized read counts) yields the files in the output. Raw sequencing data can be found on GEO ###########.

### Types of hgRNAs
- **Intronic deletion hgRNAs**: These are the main part of the screen: They are designed to crease intronic deletions, either (i) deleting most of an entire intron, (ii) a conserved block within an intron, or (iii) a conserved region overlapping an alternative exon mediating nonsense-mediated mRNA decay. Up to 16 overlapping hgRNAs target each region.
- **Intergenic-intergenic controls:** Negative controls that should not result in deletion with a phenotype. The distribution of log2-fold changes (log2FC) of these controls is used to call significant phenotypes: Deletions with a log2FC lower than the 5th percentile are assumed to mediate a significant dropout, and above the 95th percentile are assumed to mediate a significant fitness increase.
- **Deletion specificity controls:** Each Cas9 and Cas12a of an hgRNA desiged to crease a deletion is also present in combination with a partner targeting a distal intergenic regions. This control is expected to cause no or a weaker phenotype. hgRNAs where one of these controls does cause a phenotype are disregarded.
- **Exon-targeting controls**: hgRNAs where either the Cas9 or Cas12a gRNA targets a coding exon as in a classic CRISPR deletion screen. The large majority of exon-targeting Cas9 gRNAs are derived from the [Toronto CRISPR human knockout library](https://www.addgene.org/pooled-library/moffat-crispr-knockout-tkov3/) (TKOv3). The Cas9 gRNAs are used to measure fitness phenotypes of gene depletion. 
- **Exon-deletion controls:** hgRNAs designed to delete a coding exon with the Cas9 and Cas12a gRNA targeted to intronic regions on either side of the exon. These controls are used to assess the efficacy of the screen, and to infer the optimal fraction of overlapping deletion hgRNAs associated with a fitness phenotyope to call a region a hit by maximizing the differential between exonic deletions in fitness- vs. non-fitness genes.


### Steps
**1. Normalize raw counts and derive log2FC using edgeR**
```
Rscript GetLFC.R
```
**2. Call dropout hits**
Generate lists of fitness (i.e., causing decreased cell fitness when deleted) and non-fitness genes using exon-targeting controls. Identify hgRNAs with a dropout; remove those where a deletion specificity control also has a phenotype. Use exon-deletion controls to derive the optimal fraction of overlapping hgRNAs that must have a dropout to call a region a hit. Apply criteria to intron deletions and call dropout hits.
```
Rscript RegionScoring_dropout.R
```
**3. Call increase hits**
Similar to calling dropout hits, but re-use the optimized fraction of required overlapping hgRNAs with a phenotype from the analysis of dropout-causing exon deletion controls, since they should have a similar deletion efficacy.
```
Rscript RegionScoring_increase.R
```

### Reference
Shaghayegh Farhangmehr*, Ulrich Braunschweig*, Mingkun Wu, Syed Nabeel-Shah1, Kevin Brown, Jason Moffat, and Benjamin J. Blencowe. A genome-wide functional analysis of conserved intronic regions reveals essential roles for speckle-associated detained introns. *Manuscript submitted* (2025)

### Contact
For questions or comments, please raise an issue or email u -dot- braunschweig -at- utoronto -dot- ca.