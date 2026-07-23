# 🍎 Apple RNA-Seq Analysis

<p align="center">

![GitHub last commit](https://img.shields.io/github/last-commit/jeanssresende/apple-rnaseq-analysis)
![GitHub repo size](https://img.shields.io/github/repo-size/jeanssresende/apple-rnaseq-analysis)
![GitHub](https://img.shields.io/github/license/jeanssresende/apple-rnaseq-analysis)
![R](https://img.shields.io/badge/R-4.5+-276DC3?logo=r)
![Linux](https://img.shields.io/badge/Linux-Ubuntu-E95420?logo=ubuntu)
![Bioconductor](https://img.shields.io/badge/Bioconductor-3.21-green)
![RNA--seq](https://img.shields.io/badge/RNA--seq-Analysis-blue)
![Status](https://img.shields.io/badge/Status-In%20Progress-orange)

</p>

---

## 📖 About

This repository contains a complete and reproducible RNA-seq analysis workflow for **apple (Malus domestica)** transcriptomic data.

The pipeline starts from raw FASTQ files and covers all major preprocessing and downstream analysis steps using Linux, R, and Bioconductor.

---

## 🚀 Workflow

```text
Raw FASTQ
    │
    ▼
Quality Control (FastQC)
    │
    ▼
MultiQC
    │
    ▼
Read Trimming
    │
    ▼
Quality Control
    │
    ▼
Alignment / Quantification
    │
    ▼
Gene Count Matrix
    │
    ▼
Differential Expression
    │
    ▼
Functional Enrichment
    │
    ▼
Visualization
```

---

## 📂 Repository Structure

```
apple-rnaseq-analysis/
│
├── scripts/
├── metadata/
├── config/
├── docs/
├── workflow/
├── figures/
├── results/
├── README.md
└── LICENSE
```

---

## 🛠 Tools

| Tool | Purpose |
|------|---------|
| FastQC | Raw read quality assessment |
| MultiQC | QC report aggregation |
| fastp / Trimmomatic | Read trimming |
| STAR / HISAT2 / Salmon | Alignment / Quantification |
| R | Statistical analysis |
| Bioconductor | RNA-seq packages |
| DESeq2 | Differential expression |
| ggplot2 | Data visualization |

---

## 📋 Pipeline Steps

- [ ] Project setup
- [ ] Import FASTQ files
- [ ] Quality Control
- [ ] Trimming
- [ ] Alignment
- [ ] Gene Quantification
- [ ] Count Matrix
- [ ] Differential Expression
- [ ] PCA
- [ ] Heatmap
- [ ] Volcano Plot
- [ ] GO Enrichment
- [ ] KEGG Enrichment
- [ ] Final Report

---

## 💻 Operating System

Developed primarily on

- Ubuntu Linux
- Bash
- R
- RStudio

---

## 📚 Requirements

- R ≥ 4.5
- Bioconductor
- FastQC
- MultiQC
- fastp
- Salmon

---

## 📄 License

MIT License

---

## 👨‍🔬 Author

**Jean Resende**

PhD Student in Biotechnology

Bioinformatics • Transcriptomics • Immunogenomics
