# Systematic evaluation of single-cell RNA-seq analyses performance based on long-read sequencing platforms

## Abstract

The rapid development of next-generation sequencing (NGS)-based single-cell RNA sequencing (scRNA-seq) allows for detecting and quantifying gene expression in a high-throughput manner, providing a powerful tool for comprehensively understanding cellular function in various biological processes. However, the NGS-based scRNA-seq only quantifies gene expression and cannot reveal the exact transcript structures (isoforms) of each gene due to the limited read length. On the other hand, the long read length of third-generation sequencing (TGS) technologies, including Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio), enable direct reading of intact cDNA molecules. Both ONT and PacBio have been used in conjunction with scRNA-seq, but their performance in single-cell analyses has not been systematically evaluated. To address this, we generated ONT and PacBio data from the same single-cell cDNA libraries containing different amount of cells. Using NGS as a control, we assessed the performance of each platform in cell type identification. Additionally, the reliability in identifying novel isoforms and allele-specific gene/isoform expression by both platforms was verified, providing a systematic evaluation to design the sequencing strategies in single-cell transcriptome studies.

## Authors

Enze Deng<sup>1,2,\*</sup>, Qingmei Shen<sup>2,\*</sup>, Jingna Zhang<sup>2</sup>, Yaowei Fang<sup>2</sup>, Lei Chang<sup>2</sup>, Guanzheng Luo<sup>1</sup>, Xiaoying Fan<sup>2,†</sup>

1. MOE Key Laboratory of Gene Function and Regulation, Guangdong Province Key Laboratory of Pharmaceutical Functional Genes, State Key Laboratory of Biocontrol, School of Life Sciences, Sun Yat-sen University, Guangzhou 510275, China.
2. GMU-GIBH Joint School of Life Sciences, The Fifth Affiliated Hospital of Guangzhou Medical University, Guangzhou National Laboratory, Guangzhou Medical University, Guangzhou 510005, China, 

\* - These authors contributed equally  
† - Corresponding authors

## Data

- All sequencing data used for this study is publicly available on [National Genomics Data Center](https://ngdc.cncb.ac.cn/) under the accession token [PRJCA019416](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA019416).

## Code

This repository contains automated pipelines for analyzing raw sequencing data used in the paper. The pipelines are powered by [Snakemake](https://snakemake.github.io/), a workflow management system.

These Snakemake files comprise a comprehensive workflow that automates the analysis of sequencing data. The workflow is composed of numerous sub-tasks, each defined in separate Smk files. These sub-tasks are seamlessly integrated into the main Snakemake workflow, allowing for a streamlined and efficient analysis process.

## Environment

- python=3.10.6
- bcftools=1.17
- samtools=1.17
- sambamba=1.0
- pychopper=2.7.2
- fastp=0.23.3
- starcode=1.4
- bwa=0.7.17
- gatk=4.4.0.0
- longbow=0.6.14
- minimap2=2.24
- cellranger=7.0.1
- isoquant=3.3.0
- singularity=3.6.4
- deepvariant=1.5.0
- anndata=0.8.0
- scanpy=1.9.3
- biopython=1.80
- pysam=0.20.0
- gffutils=0.11.1
- pandas=2.0.2
- scipy=1.10.1
- numpy=1.23.5
- matplotlib=3.7.1
- seaborn=0.12.2

## Usage

1. Clone this repository to your local machine:
    <code>git clone https://github.com/fanxylab/scRNAseq-Pipeline.git</code>
2. Customize the workflow by editing the configuration files and adjusting parameters as needed.
3. Run Snakemake to execute the automated analysis:
   <code>snakemake --cores <number_of_cores> </code>
   Replace <number_of_cores> with the desired number of CPU cores for parallel execution.

## Contact

For any questions or issues, please open an issue in this repository, and we will be happy to assist you.
