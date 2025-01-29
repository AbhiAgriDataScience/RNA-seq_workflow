# Transcriptomic Response of Barley to *Pseudomonas syringae* Infection

## Project Overview

This repository contains an **RNA-seq analysis workflow** for studying the transcriptomic response of wild-type barley (*Hordeum vulgare*) to *Pseudomonas syringae* infection. The analysis involves **quality control, read alignment, gene quantification, differential expression analysis (DEA), and functional enrichment analysis**.

The workflow was implemented in two phases:
- **Phase 1 (Python, Jupyter Notebooks)**: Preprocessing of RNA-seq data, including quality control, read trimming, genome alignment using HISAT2, and gene quantification.
- **Phase 2 (R, DESeq2, clusterProfiler)**: Differential expression analysis (DEA), functional annotation, and Gene Ontology (GO)/KEGG enrichment analysis.

---

## **Dataset Information**
The RNA-seq dataset used in this project was sourced from:

**Paper**:  
> Gao J, Bi W, Li H, Wu J, Yu X, Liu D, and Wang X (2018)  
> *WRKY Transcription Factors Associated With NPR1-Mediated Acquired Resistance in Barley Are Potential Resources to Improve Wheat Resistance to Puccinia triticina*.  
> *Front. Plant Sci. 9:1486. doi: [10.3389/fpls.2018.01486](https://doi.org/10.3389/fpls.2018.01486)*  

This project analyzed a **subset of the dataset**, focusing only on **wild-type barley plants** due to computational and storage constraints.

---

## **Pipeline Workflow**

### ** Phase 1: RNA-seq Preprocessing (Python)**
Performed in **Jupyter Notebooks**:
1. **Download Annotation & FASTQ Files**
2. **Quality Control** using `FastQC`
3. **Read Trimming** with `Trimmomatic`
4. **Index Reference Genome** using `HISAT2`
5. **Read Alignment** to Reference Genome using `HISAT2`
6. **Convert & Sort SAM to BAM Files** using `Samtools`
7. **Gene Quantification** using `featureCounts`

### ** Phase 2: Differential Expression & Enrichment (R)**
Conducted using **R scripts**:
1. **Import Gene Counts** into `DESeq2`
2. **Data Quality Assessment** (PCA, Boxplots)
3. **Differential Expression Analysis (DEA)**
4. **Gene Annotation** with `biomaRt`
5. **GO/KEGG Enrichment Analysis** using `clusterProfiler`
6. **Generate Plots** (Volcano, MA, PCA, Heatmaps)

---


## **Results & Outputs**
All processed results, including normalized gene counts, differentially expressed genes, GO/KEGG enrichment are saved in the results/ directory.

Key outputs:

Differential Expression R_results: results/differential_expression_results.txt
Normalized Gene Counts: R_results/normalized_counts.txt
GO Enrichment Results: R_results/GO_enrichment_results.txt
KEGG Pathway Analysis: R_results/KEGG_pathway_results.txt
Plots (PCA, Volcano, MA, Heatmap)


## **Acknowledgements**
The dataset is from Gao et al., (2018), available in Frontiers in Plant Science.

This workflow was inspired by high-throughput transcriptomics pipelines for plant-pathogen interactions.

## **License**
This project is licensed under the MIT License. See the LICENSE file for details.
