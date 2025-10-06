# FAP-paper

Code accompanying the manuscript:  
**"A Low-Inflammatory Diet Rewires Metaflammation and Trained Immunity to Reprogram Immune Tone in Early Colorectal Tumorigenesis"**  
(Daveri, Lalli et al., 2025)

---

## 📌 Overview
This repository contains the R code used for statistical and multi-omics analyses performed in the study on the effects of a low-inflammatory diet (LID) in patients with Familial Adenomatous Polyposis (FAP).

The script reproduces the analyses described in the manuscript, including:
- Flow cytometry immunophenotyping  
- Lipidomic, proteomic, and transcriptomic integration  
- Statistical modeling and visualization of clinical outcomes

---

## ⚙️ Requirements
- R (≥ 4.2.0)  
- The following R packages (install if missing):  
  - dplyr, ggplot2, tidyr, ComplexHeatmap, DESeq2  
  - SNFtool, circlize, BloodGen3Module  
  - ggpubr, survival, pheatmap  
  - and other standard R packages.

> The complete list of required packages is reported at the top of the script.

---

## ▶️ Usage
1. Clone or download this repository:
   ```bash
   git clone https://github.com/USERNAME/FAP-paper.git
   ```
2. Open the R script:
   ```r
   source("FAP_LID_analysis.R")
   ```
3. Run the analysis following the instructions inside the script.

---

## 📊 Data Availability & Privacy
The data analyzed in this study are **sensitive clinical datasets** and cannot be shared publicly in this repository.

- RNA sequencing data are deposited in GEO:  
  - PBMCs: **GSE296403**  
  - Tissue: **GSE222298**  
- Other raw or clinical datasets are available **upon reasonable request** to the corresponding author of the manuscript, subject to ethical approval and privacy regulations.

If you wish to reproduce the analysis, please contact the corresponding author for data access:  
📧 elena.daveri@istitutotumori.mi.it

For testing purposes, you may substitute your own dataset with the same structure as described in the Methods section of the paper.

---

## 📜 License
This project is released under the **MIT License** (see LICENSE file).

---

## ✒️ Citation
If you use this code, please cite the associated publication:

> Daveri E*, Lalli L*, Ferrero G, et al.  
> *A Low-Inflammatory Diet Rewires Metaflammation and Trained Immunity to Reprogram Immune Tone in Early Colorectal Tumorigenesis*. 2025.

---
