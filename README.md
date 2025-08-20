# FOXC2 RNA-seq Differential Expression Analysis

## Introduction

Understanding how transcription factors regulate gene expression is essential in cancer research. FOXC2 is a key transcription factor implicated in processes such as tumor progression, epithelial-to-mesenchymal transition (EMT), and immune modulation. 

In this study, we analyzed bulk RNA-seq data from mouse melanoma cells to investigate the impact of FOXC2 knockout on the transcriptome. Using public datasets from the GEO database (accession: [GSE134296](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134296)), we compared gene expression profiles between wild-type and FOXC2-deficient samples to identify differentially expressed genes and enriched biological pathways. 

The goal was to determine whether FOXC2 deletion alters key processes such as interferon signaling and immune response, contributing to the tumor’s ability to evade immune surveillance.


## Project Structure

├── 1_download_and_trim.sh # Download, decompress, and trim FASTQ files
├── 2_kallisto_quant.sh # Transcript quantification with Kallisto
├── 3_rnaseq_analysis.R # R script for normalization, DE analysis, PCA, clustering, GO enrichment
├── samples.csv # Sample-to-condition mapping
├── data/ # Raw, trimmed, and quantified sequencing data
├── last_sec_results/ # Analysis outputs: plots, DEG tables, enrichment results

## Dataset Summary

- Organism: Mus musculus
- Cell line: B16-F1 melanoma
- Conditions: WT vs. FOXC2 knockout
- Design: 2 replicates per condition

## Software Requirements

### Shell Tools

- [SRA Toolkit](https://github.com/ncbi/sra-tools) (`prefetch`, `fasterq-dump`)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [Kallisto](https://pachterlab.github.io/kallisto/)

### R Packages

- `tximport`
- `biomaRt`
- `edgeR`
- `limma`
- `clusterProfiler`
- `org.Mm.eg.db`
- `tidyverse`
- `ggplot2`

## Running the Pipeline

1. Download, decompress, and trim reads:

1_download_and_trim.sh

2. Run Kallisto quantification:

2_kallisto_quant.sh

3. Perform differential expression and enrichment analysis in R:

3_rnaseq_analysis.R

## Output Files

+ - `voom_plot.png`: Mean-variance plot from voom transformation
+ - `volcano_plot.png`: Volcano plot of log fold change vs. significance
+ - `pca_plot.png`: Principal component analysis plot
+ - `hierarchical_clustering.png`: Clustering dendrogram
+ - `GO_upregulated.csv`, `GO_downregulated.csv`: GO enrichment results
+ - `go_enrichment_upregulated.png`, `go_enrichment_downregulated.png`: Barplots

## Results

### Quality Control and Trimming

Raw sequencing reads were trimmed using **Trimmomatic**, resulting in high-quality reads with Phred scores above 30. While FastQC reports showed minor adapter presence, pseudoalignment via **Kallisto** achieved high mapping rates (>80%), confirming the data was suitable for downstream analysis.

### Differential Expression Analysis

Differential gene expression between wild-type (WT) and FOXC2-knockout melanoma cells was assessed using the `limma` and `voom` pipeline. Significant differentially expressed genes (DEGs) were defined as:

- Adjusted p-value (FDR) < 0.05  
- |log₂ fold change| > 1  

The volcano plot illustrates DEGs, with several genes upregulated in FOXC2-knockout cells linked to immune and stress response pathways.

### Exploratory Data Analysis

- **Principal Component Analysis (PCA)** showed clear separation between WT and knockout samples along PC1, indicating FOXC2 status is the dominant source of variation.
- **Hierarchical clustering** grouped samples by condition, supporting consistency between replicates and minimal batch effects.

### Enrichment Analysis

Gene Ontology (GO) enrichment analysis (via `clusterProfiler`) revealed:

- **Upregulated genes** were significantly enriched in:
  - Type I interferon signaling
  - Response to virus
  - Regulation of viral life cycle

These findings support the hypothesis that **FOXC2 may function as an immune suppressor**, helping melanoma cells evade immune detection. Enrichment among downregulated genes was limited, potentially reflecting sample size constraints.


## Reference

- Hargadon, K.M., & Williams, C.J. (2020).  
  *FOXC2 represses interferon response genes to promote melanoma progression.*  
  [PubMed Central Article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7056877/)  

- GEO Dataset: [GSE134296](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134296)

