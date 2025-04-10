# Whole proteome analysis

## Data processing

Whole proteome data consists of six biological replicates.

- 041521_MOBN_20_WP_2024

- 051121_BN_WP_rep3_2024

- 070521_EV1-46B_10plex_2024

- 091521_EV1-55B_WP_2024

- 20201024_ev_vi_expt7_WP_2024

- 20240408_EV2_68B_2024 

Each experiment was processed individually using following filters:

- Peptides with more than 2 internal missed cleavage sites were removed.

- Half-tryptic peptides were removed.

- Contaminant peptides were removed

- Peptides with low average of reporter ion intensities (<5,000 for all technical replicates) were removed ("floating control")

- Peptides with high variation between all technical replicate channels (CV >0.5) were removed ("floating control")

- Proteins were required to have at least two unique quantified peptides to pass into the final list

- Proteins were required to be quantified in at least two unique biological replicates for analysis

- Proteins detected in two replicates with high variability (29 proteins where the protein was not changing (FC < 1.5) in one replicate and the ratio of the two replicates was greater than two) were also removed ("two replicate variability filter")

## Principal component analysis

Prior to PCA, protein ratio values were log2 transformed and only proteins quantified in all experiments were used

PCA was performed with scikit-learn version 1.4.2 for Python version 3.1.1.

## Volcano plots

p-values were calculated with T-test for the means of two independent samples and volcano plots show uncorrected −log10 p-values on y axis and fold change of KO/AB on x axis. Significant peptides (< 0.05, -log10(0.05) = 1.3) which show a > 1.5 fold change are highlighted.

## Reference gene list sources

**Metabolic genes** Metabolic genes were identified by a CRISPR screen performed by the Sabatini lab. ([Article](https://pubmed.ncbi.nlm.nih.gov/26232224/) [Data](https://www.addgene.org/pooled-library/sabatini-human-crispr-metabolic-knockout/))

**Peroxisomal genes** Peroxisome genes were identified through a text search of Gene Ontology (GO) terms in the basic gene ontology from the GO-consortium ([go-basic.obo] (https://geneontology.org/docs/download-ontology/)) containing "peroxisome". [GOATOOLS] (https://github.com/tanghaibao/goatools) was used to perform the search.
Gene names were then mapped to protein IDs through [gProfiler](https://biit.cs.ut.ee/gprofiler/convert)

**Integrated stress response** ISR signature genes were pulled from Table S1 of [Chandel 2023](https://doi.org/10.1038/s41586-023-06423-8) This gene list comes from a second study. 99 were identified from Bulk RNA-Seq using Clustering by Inferred Co-expression analysis. 35 additional genes were added through manual curation. Mouse genes were mapped to human orthologs using Ensembl v113. 2 genes were manually mapped using NCBI evidence.

**Cell cycle**

**Nucleotide metabolism**

**Redox-related**

Multiple reference lists of genes of interest were created by searching data from the Gene Ontology Consortium. 

The basic gene ontology, [go-basic.obo](https://geneontology.org/docs/download-ontology/) downloaded on June 20, 2023, was queried using [GOATOOLS v1.2.3](https://github.com/tanghaibao/goatools).

Titles, definitions, and synonyms of GO terms were searched for substrings. The reference gene lists are comprised of genes mapping to GO-terms or children terms that matched the query.

**01_go_terms_matching_query.csv** contains the gene ontology terms and associated genes that matched the query.

**02_genes_matching_query.csv** contains a list of unique genes matching to the query.


### redox related genes

- **substrings used:** ["oxid","redox","perox","oxygen"]
- **domains searched:** Biological Process

### protein import

- **substrings used:** ["protein import", "mitochondrial import","protein transport","TIM23"]
- **domains searched:** Biological Process, Molecular Function, Cellular Component

### mitochondrial fusion fisson

- **substrings used:** ["mitochondrial fusion", "mitochondrial fission"]
- **domains searched:** Biological Process, Molecular Function

### protein homeostasis folding

- **substrings used:** ["protein fold", "folded protein", "protein stabil", "protein quality"]
- **domains searched:** Biological Process, Molecular Function

### Mitochondrial genes
Mitochondiral genes were identified by the [MitoCarta3.0](https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways). This list is curated through multiple approaches, including GFP tagging, mass spectrometry, and literature review. This resource already included UniProt identifiers so no ID mapping  was performed.

### Redox related genes

Redox related genes were identified through a text search of Gene Ontology (GO) terms in the basic gene ontology from the GO-consortium ([go-basic.obo] (https://geneontology.org/docs/download-ontology/)) containing "oxid","oxygen",  or "redox". Search results were manually curated to 15 GO terms based on relevance. [GOATOOLS] (https://github.com/tanghaibao/goatools) was used to perform the search.
Gene names were then mapped to protein IDs through [gProfiler](https://biit.cs.ut.ee/gprofiler/convert)

## GO-term enrichment analysis

Proteins that passed significance (adjusted p-value < 0.05) and fold change (FC > 1.5) cutoffs were used for analysis. The basic gene ontology from the GO consortium ([go-basic.obo](https://geneontology.org/docs/download-ontology/) downloaded on June 20, 2023) was used. An experiment-specific background gene list was used. Both up and down regulated proteins were used for analysis. Uniprot IDs were mapped to gene IDs using gProfiler. GO term analysis was performed with [GOATOOLS v1.2.3](https://github.com/tanghaibao/goatools). Plotting was performed in RStudio, Sankey GO-term plots were created in SRPlot.

## Whole proteome vs RNA correlation analysis

R represents Pearson’s R. Mean was used to combine replicates before calculating fold change. A fold change cutoff of 1.5 was used to color proteins by quadrant. RNA differential expression analysis was performed by Jahan Rahman in DESeq2. Fold change values were not filtered by significance.


References

1.	Klopfenstein, D. V., Zhang, L., Pedersen, B. S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., ... & Tang, H. (2018). GOATOOLS: A Python library for Gene Ontology analyses. Scientific reports, 8(1), 1-17.

2.	Tang, D., Chen, M., Huang, X., Zhang, G., Zeng, L., Zhang, G., ... & Wang, Y. (2023). SRplot: A free online platform for data visualization and graphing. PLoS One, 18(11), e0294236.
