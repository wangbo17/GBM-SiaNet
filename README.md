<p align="center"><img height="800" alt="GBM-SiaNet_logo" src="https://github.com/user-attachments/assets/1cab93c4-f74c-46a3-b988-814f55556646" /></p>

<p align="center">
  <img src="https://img.shields.io/badge/R-4.3.2+-blue" />
  <img src="https://img.shields.io/badge/Python-3.10.14+-blue" />
  <img src="https://img.shields.io/badge/License-MIT-green" />
  <img src="https://img.shields.io/badge/Status-Beta-orange" />
  <img src="https://img.shields.io/github/last-commit/wangbo17/GBM-SiaNet" />
</p>

## 📖 Overview

**GBM-SiaNet** is an interpretable machine learning classifier based on a Siamese neural network, designed to identify key genes associated with therapy resistance in glioblastoma (GBM) by distinguishing between primary and recurrent tumours using longitudinal gene expression data. This MSc Research Project in precision medicine and bioinformatics was independently designed and implemented.

## 🧠 Background

Glioblastoma (GBM) is an incurable brain cancer, with fatal recurrence after standard treatment. Recurrence typically occurs within 6–9 months, with a median post-recurrence survival time of only 3–9 months. Current evidence suggests that treatment resistance may not be a consequence of somatic genomic alterations, shifting the research focus towards understanding the role of transcriptional heterogeneity.

## ⚙️ Methods and Objectives

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

## 📋 Feature Selection

A combined feature selection strategy was implemented by integrating two filter-based methods. Each algorithm was independently applied to rank features on the training set. The top 100 to 500 features (in increments of 100) from each were merged to generate multiple integrated feature sets.

**Minimum Redundancy Maximum Relevance (mRMR) Algorithm** (Ding and Peng, 2005)

- Effectively reduces feature redundancy by selecting features highly relevant to the target class.
- Limitation: Does not consider interactions between features.

**Relief Algorithm** (Urbanowicz et al., 2018)

- Captures feature interactions by ranking features by differences between neighbouring samples.
- Limitation: Less efficient at reducing redundancy.

<p align="center">
  <img height="300" alt="F0" src="https://github.com/user-attachments/assets/84de1fff-7174-4211-a8bf-6c54c21fd420" />
  <br>
  <strong>Figure 1. Combined Feature Sets Selected by mRMR and Relief Algorithms.</strong>
</p>

## 📈 Model Development

### Traditional Machine Learning Models

To establish baseline performance, three traditional machine learning models were implemented: Logistic Regression (LR), Random Forest (RF), and Support Vector Machine (SVM). Hyperparameters were tuned via 5-fold cross-validation on the training set. 

**Table 2: Average Cross-Validation Accuracies of ML Models.**

| **Model** | **Top 100** | **Top 200** | **Top 300** | **Top 400** | **Top 500** |
| --------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| **LR**    | 0.9046      | 0.9429      | 0.9554      | **0.9746**  | 0.9683      |
| **RF**    | 0.8283      | 0.8162      | 0.8283      | **0.8404**  | 0.8354      |
| **SVM**   | 0.9046      | 0.9054      | **0.9113**  | 0.9050      | 0.9108      |

The models were then evaluated on the test set using feature sets of varying sizes. Among them, SVM consistently achieved the highest accuracy, particularly with the top 300 combined features.

<p align="center">
  <img height="350" alt="F2" src="https://github.com/user-attachments/assets/267d3e32-41b4-46ac-bdaf-52df38332d51" />
  <br>
  <strong>Figure 2. Performance of ML Models with Combined Feature Sets on the Test Set.</strong>
</p>

⚠️ *These models do not account for the paired structure or batch effects, and their failure to capture key patterns limits interpretability.*

### Siamese Neural Network-Based Hybrid Models

To address these limitations, a Siamese Neural Network (SNN) was developed to distinguish between primary and recurrent GBM samples by mapping them into a common feature space. The SNN consists of  three identical neural networks to process three inputs sample at the same time: anchor, positive, and negative. In this study, the anchor could be either a primary or recurrent tumour, with the negative sample being the other tumour class as the anchor from the same patient. The positive sample should be the same tumour class as the anchor but from a different patient within the same centre. The training objective was to maximize the distance between samples of different classes from the same patients and minimize the distance between samples of the same class from different patients within the same centre. This design enables the model to capture the key differences and similarities between matched tumour samples, while accounting for centre-specific batch effects by comparing samples within each centre.

<p align="center">
  <img height="400" alt="F3" src="https://github.com/user-attachments/assets/16c234bb-6422-4caa-900d-71a8e6490c0c" />
  <br>
  <strong>Figure 3. Objective of the Siamese Neural Network (SNN).</strong>
</p>

Given the one-dimensional nature of bulk RNA-seq data, a Fully Connected Network (FCN) was selected as the sub-network architecture. The FCN includes three fully connected layers, each followed by batch normalization, ReLU activation, and dropout layers. He initialization was used for each dense layer, with L2 regularization applied to reduce the risk of overfitting.

<p align="center">
  <img height="400" alt="F4" src="https://github.com/user-attachments/assets/25f016ef-4807-4153-aada-e27a76f878d6" />
  <br>
  <strong>Figure 4. Architecture of the Siamese Neural Network (SNN).</strong>
</p>

The core of the learning process is the triplet loss function, which compares the anchor with positive and negative samples based on their Euclidean distances. The function ensures that the negative is farther from the anchor than the positive by a certain margin, allowing same-type samples to cluster together while different ones are separated. However, it has a limitation: once the negative distance is far enough from the positive distance, the loss becomes zero, and the model stops learning. To optimise this, the study used a modified triplet loss function, replacing the Max function with the SoftPlus function. SoftPlus is a smooth approximation of Max function and continues to contribute to the loss even when the margin is satisfied. This enables smoothly continuous separation between primary and recurrent GBM, making it highly effective for the binary classification task.

<p align="center">
  <img height="300" alt="F5" src="https://github.com/user-attachments/assets/d5456f34-9c2b-47a0-8a53-cb4c76c72d11" />
  <br>
  <strong>Figure 5. Comparison Between Standard and SoftPlus-Modified Triplet Loss.</strong>
</p>

SNN was developed using the top 300 combined feature set, selected for their strong performance in traditional ML models. Due to time limitations, only this feature set was used, although ideally, all combined sets would have been explored. Optimal hyperparameters were identified through 5-fold cross-validation, achieving an average accuracy of 0.895 in distinguishing between positive and negative samples in each triplet from the validation subsets. The final model was trained on 3,368 triplets from the training subset, with 168 triplets used for validation. t-SNE showed a clear separation of primary and recurrent samples in the shared feature space, significantly improving compared to the raw input data, but there is still room for improvement.

<p align="center">
  <img height="600" alt="F6" src="https://github.com/user-attachments/assets/c0d2710c-a98f-458a-ac45-3f57858dbf46" />
  <br>
  <strong>Figure 6. t-SNE Visualization of Feature Spaces Before and After SNN Transformation.</strong>
</p>

Building on the feature space created by the SNN, traditional ML models were used to create hybrid models for further classification. In these hybrid models, the SNN acted as a feature extractor, producing representations that were then fed into traditional ML models for final classification. After optimizing the models through cross-validation, the final models were trained on the full training set and evaluated on the test set. These hybrid models consistently performed better than traditional ML models, with the SNN-LR achieving the best performance, showing an accuracy of 0.90.

<p align="center">
  <img height="325" alt="F7" src="https://github.com/user-attachments/assets/f0787aec-27ab-4009-8fcf-72c0c3f17f94" />
  <br>
  <strong>Figure 7. Performance Comparison of SNN-Based Hybrid and Traditional ML Models on the Test Set.</strong>
</p>

## 🔎 Model Interpretation

### Shapley Additive exPlanations (SHAP) Analysis

To identify the most critical genes involved in classification, SHAP analysis was used to interpret the decision-making process of the best-performing model, SNN-LR. SHAP values quantify each feature's contribution to the model's predictions. Positive SHAP values suggest a push towards recurrent GBM, while negative values suggest a push towards primary GBM. The colour bar represents gene expression levels, with red indicating high expression and blue indicating low expression. For example, lower expression levels of genes like MAG are linked with recurrent samples, while higher expression levels of genes like HOTAIR are linked with recurrent samples.

<p align="center">
  <img width="650" height="1595" alt="F8" src="https://github.com/user-attachments/assets/2cbdb0a0-a418-40e2-ae36-d403f278f94e" />
  <br>
  <strong>Figure 8. SHAP Summary Plot of the 20 Most Influential Features.</strong>
</p>

A literature search was conducted focusing on these 20 genes in the context of GBM recurrence and therapy resistance. The search confirmed that several of these genes are known to be associated with therapy resistance in GBM. Also, these genes are associated with drug resistance in other cancers, although their roles in GBM are not yet described. While these genes are typically involved in neural integrity and signalling, their roles in GBM are less clear. However, studies suggest that GBM's interaction with neural circuits may worsen tumour growth and affect survival. The study also identified novel genes which could lead to further research on GBM recurrence and therapy resistance.

**Table 3: Key Genes Identified by SHAP with Supporting Literature**

| **Gene**    | **Category**                | **Description / Mechanism**                                  |
| ----------- | --------------------------- | ------------------------------------------------------------ |
| HOTAIR      | GBM Treatment Resistance    | Epigenetic regulation, miRNA sponging, enhances proliferation, invasion, chemoresistance, angiogenesis |
| HOXC10      | GBM Treatment Resistance    | Activates PI3K/AKT pathway, promotes cell growth, inhibits apoptosis |
| L1CAM       | GBM Treatment Resistance    | Activates DNA damage repair, promotes vasculogenic mimicry, contributes to drug resistance |
| CNTN2       | GBM Treatment Resistance    | Interacts with APP, activates RTK/Ras/MAPK pathway, promotes proliferation and inhibits differentiation |
| ENPP2       | GBM Treatment Resistance    | Increases LPA production, drives migration and invasion, reduces treatment efficacy |
| GPNMB       | GBM Treatment Resistance    | Facilitates mesenchymal transition, immune evasion, enhances invasion and resistance |
| TPPP        | GBM Treatment Resistance    | Promotes epithelial-mesenchymal transition (EMT), increases invasiveness, reduces apoptosis |
| DCAF4L2     | Resistance in Other Cancers | Liver cancer, Colorectal cancer                              |
| TCEAL6      | Resistance in Other Cancers | Cervical cancer                                              |
| RNASE1      | Resistance in Other Cancers | Breast cancer                                                |
| MAG         | Neuronal Signalling         | Myelin associated glycoprotein                               |
| MYRF        | Neuronal Signalling         | Myelin regulatory factor                                     |
| KCNH8       | Neuronal Signalling         | Potassium voltage-gated channel subfamily H member 8         |
| SYT4        | Neuronal Signalling         | Synaptotagmin 4                                              |
| ERMN        | Neuronal Signalling         | Cytoskeleton-related oligodendroglial protein                |
| KHDRBS2-OT1 | Novel Candidate Gene        | RNA, long non-coding                                         |
| SNX18P3     | Novel Candidate Gene        | Sorting Nexin 18 Pseudogene 3                                |
| RN7SL683P   | Novel Candidate Gene        | RNA, 7SL, Pseudogene 683                                     |
| FAM162B     | Novel Candidate Gene        | Family With Sequence Similarity 162 Member B                 |



## 🧬 Biological Validation

### Differential Gene Expression Analysis 

Differential Gene Expression (DGE) analysis was conducted using DESeq2 with a paired design to identify genes with significantly different expression levels between primary and matched recurrent GBM samples. Genes were filtered to include only those that were expressed in at least 50% of samples from either primary or recurrent. The analysis revealed 839 Differentially Expressed Genes (DEGs), with 601 upregulated and 238 downregulated. The DGE analysis results were then compared with SHAP results, leading to some interesting findings. Although most genes showed strong agreement between these two methods, some differences were noted. For example, MAG and MYRF had higher expression in recurrent samples according to DGE, but SHAP results linked their increased expression to primary GBM. Therefore, the predicted sample types based on the expression patterns of certain key genes sometimes did not match the actual types.

<p align="center">
  <img height="400" alt="F9" src="https://github.com/user-attachments/assets/efc191a5-b7c4-4490-ba6d-50a1667a822f" />
  <br>
  <strong>Figure 9. Volcano Plot of Differential Gene Expression (DGE) Analysis Between Primary and Recurrent GBM Samples.</strong>
</p>

Further per-sample analysis showed variability in SHAP values for these two genes across individual samples. In these two plots, each point represents a sample, with quadrants indicating different predictive tendencies: the top-right area with positive SHAP values suggests a bias towards predicting recurrent GBM, while the bottom-left quadrant with negative SHAP values suggests a bias towards predicting primary GBM. Notably, for primary GBM samples, SHAP values often incorrectly suggested a recurrent classification, while in recurrent GBM samples, about 50% correctly suggested recurrence GBM, with the rest incorrectly suggesting primary GBM. These findings suggest that the roles of these genes in GBM recurrence and therapy resistance may differ among individuals, indicating that patient stratification could be necessary. However, it is important to note that while SHAP effectively captures feature interactions, but it assumes that features are independent. When features are highly correlated, this assumption may lead to misleading interpretations.

<p align="center">
  <img height="325" alt="F10" src="https://github.com/user-attachments/assets/669b2281-f9c8-42e4-9001-0031f64df357" />
  <br>
  <strong>Figure 10. SHAP Value Scatter Plot for MAG and MYRF Genes.</strong>
</p>

<p align="center">
  <img height="575" alt="F11" src="https://github.com/user-attachments/assets/c88f3c54-f7a7-4e55-bf62-9685f78aa34a" />
  <br>
  <strong>Figure 11. Correlation Plot for the 20 Most Influential Features Identified by SHAP Analysis.</strong>
</p>

### Pathway and Functional Enrichment Analysis

The resulting DEGs were analysed with Gene Set Enrichment Analysis (GSEA) to identify enriched pathways associated with GBM recurrence and therapy resistance. This analysis identified significant upregulation in pathways related to neurotransmitter signalling and synaptic transmission, consistent with the top genes identified through SHAP analysis. This further supports the idea that GBM cells might integrate into neural circuits, promoting tumour growth and therapy resistance. In contrast, pathways associated with cell cycle regulation were downregulated, potentially allowing tumour cells to survive and grow despite treatment.

<p align="center">
  <img height="600" alt="F12" src="https://github.com/user-attachments/assets/d67c70b6-5377-4449-b63c-e511cfdc643a" />
  <br>
  <strong>Figure 12. GSEA of DEGs Between Primary and Recurrent GBM Samples.</strong>
</p>

Over-Representation Analysis (ORA) was conducted to determine whether the key genes identified by SHAP were significantly enriched in specific gene sets, using an adjusted p-value threshold of < 0.01. This analysis identified significant enrichment in the 'Glial Cell Differentiation' pathway and pathways related to oligodendrocytes (OLs), which are cells in the central nervous system that form the myelin sheath. This suggests that OLs and their precursor cells may interact with the tumour, affecting cancer progression. Additionally, a gene set associated with quiescent neural stem cells (NSCs) was enriched, indicating that GBM cells might use these quiescent mechanisms to survive therapy and promote recurrence. These findings highlight the complex interactions between GBM and its microenvironment, particularly involving OLs and quiescent NSCs, which may help the tumour escape therapies and contribute to recurrence.

**Table 4: Significantly Enriched Gene Sets Identified by Over-Representation Analysis (ORA)**

| **Gene Set**                                                 | **Size** | **P.adjust** **(GSEA)** | **P.adjust** **(ORA)** | **Overlap Genes**                      | **Source**   |
| ------------------------------------------------------------ | -------- | ----------------------- | ---------------------- | -------------------------------------- | ------------ |
| Allen 116 OligoL4 6OPALIN (OLs)                              | 53       | 6.72e-14                | 2.28e-06               | *MYRF,  KCNH8, CNTN2, ENPP2, RNASE1*   | Literature   |
| Darmanis-Barres Oligodendrocytes  (OLs)                      | 20       | 1.58e-09                | 2.95e-06               | *MAG,  ERMN, ENPP2, RNASE1*            | Literature   |
| Allen 115 OligoL4 6MOBPCOL18A1 (OLs)                         | 57       | 3.37e-12                | 1.57e-04               | *KCNH8,  CNTN2, ENPP2, RNASE1*         | Literature   |
| Llorens-Bobadilla AdultSvcNSCs cluster1 (quiescent neural stem  cells ) | 334      | 8.31e-32                | 1.98e-04               | *MAG,  MYRF, ERMN, CNTN2, ENPP2, TPPP* | Literature   |
| Zhang-Barres oligodendrocytes (OLs)                          | 49       | 5.65e-18                | 6.27e-03               | *ERMN,  CNTN2, ENPP2*                  | Literature   |
| Glial Cell Differentiation  (WP2276)                         | 7        | 1.19e-02                | 8.17e-03               | *MAG,  TPPP*                           | WikiPathways |

## 🔬 Conclusion

**High-Performance Model**

- **Siamese Neural Network (SNN) combined with Logistic Regression (LR)** achieved **90%** accuracy, outperforming other models in distinguishing primary and recurrent GBM samples.

**Key Genes Identification**

- Critical genes associated with recurrence and resistance, including **HOTAIR**, **HOXC10**, **L1CAM**, and novel candidates **KHDRBS2-OT1**, **SNX18P3**, **RN7SL683P**, and **FAM162B**.

**Gene-Specific Mechanisms**

- Evidence of **gene-specific** **resistance** **mechanisms**, varying by individual patient, highlighting the necessity for patient stratification in GBM treatment.

## 💾 Reference

Korotkevich, G., Sukhov, V. and Sergushichev, A. 2016. Fast Gene Set Enrichment Analysis. *bioRxiv*.

Krishna, S., Choudhury, A., Keough, M.B., Seo, K., Ni, L., Kakaizada, S., Lee, A., Aabedi, A., Popova, G., Lipkin, B., Cao, C., Nava Gonzales, C., Sudharshan, R., Egladyous, A., Almeida, N., Zhang, Y., Molinaro, A.M., Venkatesh, H.S., Daniel, A.G.S. and Shamardani, K. 2023. Glioblastoma Remodelling of Human Neural Circuits Decreases Survival. *Nature*. **617**(7961), pp.599–607.

Lee, H.-H., Wang, Y.-N., Yang, W.-H., Xia, W., Wei, Y., Chan, L.-C., Wang, Y.-H., Jiang, Z., Xu, S., Yao, J., Qiu, Y., Hsu, Y.-H., Hwang, W.-L., Yan, M., Cha, J.-H., Hsu, J.L., Shen, J., Ye, Y., Wu, X. and Hou, M.-F. 2021. Human Ribonuclease 1 Serves as a Secretory Ligand of Ephrin A4 Receptor and Induces Breast Tumor Initiation. *Nature Communications*. **12**(1).

Love, M.I., Huber, W. and Anders, S. 2014. Moderated Estimation of Fold Change and Dispersion for RNA-seq Data with DESeq2. *Genome Biology*. **15**(12), p.550.

Mallick, S., Benson, R., Hakim, A. and Rath, G.K. 2016. Management of Glioblastoma after recurrence: a Changing Paradigm. *Journal of the Egyptian National Cancer Institute*. **28**(4), pp.199–210.

Mohammed, S., M, D. and T, A. 2022. Survival and Quality of Life Analysis in Glioblastoma Multiforme with Adjuvant chemoradiotherapy: a Retrospective Study. *Reports of Practical Oncology and Radiotherapy*. **27**(6), pp.1026–1036.

Nahm, J., Sinha, M., Highberg Schumann, E., Vu, M., Hsu, S.H. and Zhu, J.-J. 2023. Overall Survival in Patients with Recurrent Glioblastomas with Combination Chemotherapy and Tumor Treating Fields (TTF). *Journal of Clinical Oncology*. **41**(16_suppl), pp.e14057–e14057.

Qazi, M.A., Vora, P., Venugopal, C., Sidhu, S.S., Moffat, J., Swanton, C. and Singh, S.K. 2017. Intratumoral heterogeneity: Pathways to Treatment Resistance and Relapse in Human Glioblastoma. *Annals of Oncology*. **28**(7), pp.1448–1456.

Rabah, N., Ait Mohand, F.-E. and Kravchenko-Balasha, N. 2023. Understanding Glioblastoma Signaling, Heterogeneity, Invasiveness, and Drug Delivery Barriers. *International Journal of Molecular Sciences*. **24**(18), pp.14256–14256.

Rosati, D., Palmieri, M., Brunelli, G., Morrione, A., Iannelli, F., Frullanti, E. and Giordano, A. 2024. Differential Gene Expression Analysis Pipelines and Bioinformatic Tools for the Identification of Specific biomarkers: a Review. *Computational and Structural Biotechnology Journal*. **23**.

Spiteri, I., Caravagna, G., Cresswell, G.D., Vatsiou, A., Nichol, D., Acar, A., Ermini, L., Chkhaidze, K., Werner, B., Mair, R., Brognaro, E., Verhaak, R., Sanguinetti, G., Piccirillo, S., Watts, C. and Sottoriva, A. 2019. Evolutionary Dynamics of Residual Disease in Human Glioblastoma. *Annals of Oncology*. **30**(3), pp.456–463.

Stupp, R., Mason, W.P., van den Bent, M.J., Weller, M., Fisher, B., Taphoorn, M.J.B., Belanger, K., Brandes, A.A., Marosi, C., Bogdahn, U., Curschmann, J., Janzer, R.C., Ludwin, S.K., Gorlia, T., Allgeier, A., Lacombe, D., Cairncross, J.G., Eisenhauer, E. and Mirimanoff, R.O. 2005. Radiotherapy plus Concomitant and Adjuvant Temozolomide for Glioblastoma. *New England Journal of Medicine*. **352**(10), pp.987–996.

Stylli, S.S. 2020. Novel Treatment Strategies for Glioblastoma. *Cancers*. **12**(10), p.2883.

Sui, X., Wang, Y. and Liu, H. 2021. hsa_circ_0101119 Facilitates the Progression of Cervical Cancer via an Interaction with EIF4A3 to Inhibit TCEAL6 Expression. *Molecular Medicine Reports*. **24**(3).

Tanner, G., Barrow, R., Ajaib, S., Al-Jabri, M., Ahmed, N., Pollock, S., Finetti, M., Rippaus, N., Bruns, A.F., Syed, K., Poulter, J.A., Matthews, L., Hughes, T., Wilson, E., Johnson, C., Varn, F.S., Brüning-Richardson, A., Hogg, C., Droop, A. and Gusnanto, A. 2024. IDHwt Glioblastomas Can Be Stratified by Their Transcriptional Response to Standard treatment, with Implications for Targeted Therapy. *Genome Biology*. **25**(1).

Thomas, M.P.H., Shoaib Ajaib, Tanner, G., Bulpitt, A.J. and Stead, L.F. 2024. GBMPurity: a Machine Learning Tool for Estimating Glioblastoma Tumour Purity from Bulk RNA-seq Data. *BioRxiv* *(Cold Spring Harbor Laboratory)*.

Urbanowicz, R.J., Olson, R.S., Schmitt, P., Meeker, M. and Moore, J.H. 2018. Benchmarking relief-based Feature Selection Methods for Bioinformatics Data Mining. *Journal of Biomedical Informatics*. **85**, pp.168–188.

Wang, Q., Hu, B., Hu, X., Kim, H., Squatrito, M., Scarpace, L., deCarvalho, A.C., Lyu, S., Li, P., Li, Y., Barthel, F., Cho, H.J., Lin, Y.-H., Satani, N., Martinez-Ledesma, E., Zheng, S., Chang, E., Sauvé, C.-E.G., Olar, A. and Lan, Z.D. 2017. Tumor Evolution of Glioma-Intrinsic Gene Expression Subtypes Associates with Immunological Changes in the Microenvironment. *Cancer Cell*. **32**(1), pp.42-56.e6.

Watson, D.S. 2021. Interpretable Machine Learning for Genomics. *Human Genetics*. **141**.

Xiong, A., Zhang, J., Chen, Y., Zhang, Y. and Yang, F. 2022. Integrated single-cell Transcriptomic Analyses Reveal That GPNMB-high Macrophages Promote PN-MES Transition and Impede T Cell Activation in GBM. *EBioMedicine*. **83**, pp.104239–104239.

Zhang, X., Xu, S., Hu, C., Fang, K., Zhou, J., Guo, Z., Zhu, G. and Li, L. 2020. LncRNA ST8SIA6-AS1 Promotes Hepatocellular Carcinoma Progression by Regulating MAGEA3 and DCAF4L2 Expression. *Biochemical and Biophysical Research Communications*. **533**(4), pp.1039–1047.
