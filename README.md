<img width="2809" height="1870" alt="GBM-SiaNet_logo" src="https://github.com/user-attachments/assets/1cab93c4-f74c-46a3-b988-814f55556646" />

<p align="center">
  <img src="https://img.shields.io/badge/R-4.3.2+-blue" />
  <img src="https://img.shields.io/badge/Python-3.10.14+-blue" />
  <img src="https://img.shields.io/badge/License-MIT-green" />
  <img src="https://img.shields.io/badge/Status-Beta-orange" />
  <img src="https://img.shields.io/github/last-commit/wangbo17/GBM-SiaNet" />
</p>

## 📖 Overview

**GBM-SiaNet** is a siamese neural network–based interpretable machine learning classifier designed to distinguish between primary and recurrent glioblastoma (GBM) using longitudinal gene expression data, while identifying key genes associated with therapy resistance.

## 🧠 Background

Glioblastoma (GBM) is an incurable brain cancer, with fatal recurrence after standard treatment. Standard of care includes maximal surgical resection followed by radiotherapy with concurrent and adjuvant temozolomide chemotherapy. However, progression-free survival remains poor, with recurrence typically occurring 6 to 9 months post-therapy. Current evidence suggests that treatment resistance may not be a consequence of somatic genomic alterations, shifting the research focus towards understanding the role of transcriptional heterogeneity.

## ⚙️ Methods and Objectives

This study developed a classifier based on longitudinal gene expression data to accurately differentiate between primary and recurrent GBM samples. iML methods were applied to the best-performing model to identify the genes most critical for classification. The biological functions of these key genes were subsequently investigated using bioinformatics approaches. 

**Hypotheses**

- Recurrent GBMs show distinct gene expression patterns compared to primary GBMs.
- Gene expression differences between recurrent and primary GBMs are linked to therapy resistance.
- Interpretable machine learning (iML) can capture and interpret these gene expression patterns.

**Objectives**

- Develop ML classifiers using gene expression data to differentiate primary from recurrent GBMs.
- Apply iML methods to the best-performing model to identify key genes involved in this classification.
- Investigate the roles of these genes in GBM recurrence and drug resistance through bioinformatics.

## 📦 Data Preprocessing

**Data Collection**

- RNA-seq data was obtained from longitudinally matched GBM tumour samples, collected from two surgical performed on the same patients.
- Patients underwent initial resection of primary tumours followed by standard treatment; recurrent tumours were resected upon recurrence.

**Selection Criteria**

- IDH-wildtype GBM  
- Received standard treatment  
- Local recurrence  
- ≥40% tumour purity  
- Total RNA-Seq library  

**Data Normalization**

- Transcripts Per Million (TPM) values were calculated and log-transformed as log₂(TPM + 1).

**Batch Correction**

- Batch correction was avoided to preserve biological variability and prevent extra assumptions.

**Data Splitting**

- Training and Test sets were selected through stratified sampling based on data sources.

**Table 1: Summary of Sample Source and Distribution**

| **Data Source**                 | **Training Set** | **Test Set** | **Total** |
| ------------------------------- | ---------------- | ------------ | --------- |
| Stead (Tanner et al., 2024)     | 48               | 16           | 64        |
| DFKZ (Körber et al., 2019)      | 16               | 6            | 22        |
| Rabadan (Wang et al., 2016)     | 4                | 2            | 6         |
| EORTC (Hoogstrate et al., 2023) | 90               | 26           | 116       |
| **Total**                       | **158**          | **50**       | **208**   |

*Test Set was held out from model training for unbiased evaluation.*

## Feature Selection


## Model Development


![GBM_SNN](https://github.com/user-attachments/assets/0dc7456a-8cf1-48be-9fb0-b0d6d21ba1a7)
