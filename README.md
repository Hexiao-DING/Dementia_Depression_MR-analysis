<!--
README redesigned for GitHub presentation.
Design principle: make the project immediately understandable, data-transparent, reproducible, and visually clean.
-->

<div align="center">

# 🧬 Dementia–Depression Mendelian Randomization Analysis Pipeline

### Causal inference pipeline for dementia–depression comorbidity, peripheral molecular traits, CSF metabolites, and SES-adjusted genetic evidence

<p>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-F2C94C?style=for-the-badge" alt="License: MIT"></a>
  <a href="https://www.r-project.org/"><img src="https://img.shields.io/badge/R-%E2%89%A5%204.0-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R >= 4.0"></a>
  <img src="https://img.shields.io/badge/Platform-Windows%20%7C%20Linux-6C757D?style=for-the-badge" alt="Platform: Windows | Linux">
  <img src="https://img.shields.io/badge/Status-Production%20Ready-2EA44F?style=for-the-badge" alt="Status: Production Ready">
</p>

<p>
  <img src="https://img.shields.io/badge/Methods-UVMR%20%7C%20MVMR%20%7C%20Mediation-5B7FA6?style=flat-square" alt="Methods">
  <img src="https://img.shields.io/badge/QC-MR--PRESSO%20%7C%20MRlap%20%7C%20FDR-C9956C?style=flat-square" alt="Quality control">
  <img src="https://img.shields.io/badge/Data-pQTL%20%7C%20mQTL%20%7C%20GWAS%20Summary%20Statistics-7A9E7E?style=flat-square" alt="Data layers">
  <img src="https://img.shields.io/badge/Ancestry-Mostly%20European-9B89AC?style=flat-square" alt="Most data are European ancestry">
</p>

<p>
  <a href="#-why-this-pipeline">Why this pipeline</a> ·
  <a href="#-quick-start">Quick start</a> ·
  <a href="#-data-atlas-most-data-are-european-ancestry">Data atlas</a> ·
  <a href="#-methodological-core">Methods</a> ·
  <a href="#-quality-control-and-result-filtering">Quality control</a> ·
  <a href="#-citation-and-references">Citation</a>
</p>

<img width="1846" height="1041" alt="Dementia–Depression Mendelian Randomization analysis framework" src="https://github.com/user-attachments/assets/4137db74-5fbc-40d0-8f71-a9d91a9aed03" />

</div>

---

> [!IMPORTANT]
> This repository contains the **analysis scripts and reproducible pipeline structure**, but it does **not** redistribute GWAS, pQTL, mQTL, LD reference, or MRlap reference files. Users must download the required public summary statistics and reference panels from the listed sources and configure local paths before running the full analysis.

---

## ✨ Project Snapshot

<table>
<tr>
<td width="25%" align="center"><b>🧠 Disease focus</b><br>Dementia, cognitive performance, depression-related traits, and dementia–depression comorbidity</td>
<td width="25%" align="center"><b>🧬 Molecular layer</b><br>Plasma proteome, inflammatory proteins, circulating metabolic biomarkers, and CSF metabolites</td>
<td width="25%" align="center"><b>🔬 Causal design</b><br>UVMR, MVMR, mediation MR, reverse MR, and sensitivity validation</td>
<td width="25%" align="center"><b>🛡️ Bias control</b><br>SES adjustment, MR-PRESSO, MRlap, heterogeneity testing, F-statistics, and FDR correction</td>
</tr>
</table>

| Item | Summary |
|:--|:--|
| **Primary scientific aim** | Identify genetically supported molecular pathways linking peripheral biology, central metabolic mediators, dementia phenotypes, and depression-related phenotypes. |
| **Exposure layer** | Circulating human plasma proteome, circulating metabolic biomarkers, and circulating inflammatory proteins. |
| **Mediator layer** | Cerebrospinal fluid metabolomics / metabQTLs. |
| **Outcome layer** | Dementia, Alzheimer's disease, cognitive performance, vascular dementia, Lewy body dementia, frontotemporal dementia, undefined dementia, depressive disorder, major depressive disorder, and mixed anxiety/depression. |
| **Covariate layer** | Socioeconomic status covariates: education, income, and occupation. |
| **Main outputs** | `uvmr_comprehensive_results.csv`, `mvmr_comprehensive_results.csv`, and `mediation_comprehensive_results.csv`. |
| **Recommended use** | Hypothesis generation, causal candidate prioritisation, molecular pathway screening, and robust MR-based triangulation. |

---

## 📚 Table of Contents

- [Why this pipeline](#-why-this-pipeline)
- [Key Features](#-key-features)
- [Quick Start](#-quick-start)
- [Data Atlas — Most data are European ancestry](#-data-atlas-most-data-are-european-ancestry)
- [Required File Schema](#-required-file-schema)
- [Instrument Selection](#-instrument-selection)
- [Methodological Core](#-methodological-core)
- [Installation](#-installation)
- [Project Structure](#-project-structure)
- [Running the Analysis](#-running-the-analysis)
- [Output Files and Key Columns](#-output-files-and-key-columns)
- [Quality Control and Result Filtering](#-quality-control-and-result-filtering)
- [Advanced Features](#-advanced-features)
- [Troubleshooting](#-troubleshooting)
- [Citation and References](#-citation-and-references)
- [Acknowledgements](#-acknowledgements)

---

## 🎯 Why this Pipeline

Dementia and depression are clinically intertwined, biologically complex, and strongly influenced by systemic, central, and social determinants. This pipeline is designed for a layered causal-inference question:

> **Which peripheral molecular traits and CSF metabolic mediators show genetically supported causal relevance to dementia and depression phenotypes after rigorous instrument selection, sensitivity analysis, sample-overlap correction, and socioeconomic adjustment?**

The pipeline combines three complementary MR modules:

| Module | Purpose | Main interpretation |
|:--|:--|:--|
| **Univariable MR** | Tests total causal associations between each exposure and each outcome. | Whether a molecular trait is genetically associated with a disease phenotype through an MR framework. |
| **Multivariable MR** | Estimates independent effects after adjusting for mediators and mandatory SES covariates. | Whether an exposure remains associated with the outcome after accounting for correlated traits and SES. |
| **Mediation MR** | Decomposes exposure → mediator → outcome pathways. | Whether a CSF metabolite or intermediate trait plausibly carries part of the causal signal. |

---

## ✨ Key Features

| Feature | Implementation |
|:--|:--|
| **Comprehensive MR methods** | UVMR with IVW, Weighted Median, MR-Egger, Weighted Mode, and Simple Mode; MVMR with IVW, Median, Egger, and Lasso. |
| **Strict instrument selection** | Genome-wide significant SNPs, LD clumping, minimum SNP count, F-statistic screening, and conditional F-statistics for MVMR. |
| **Sensitivity analysis** | MR-PRESSO outlier detection, Cochran's Q heterogeneity testing, reverse MR, direction-concordance checks, and method-concordance validation. |
| **Sample-overlap correction** | MRlap implementation with local LD scores and HapMap3 reference files. |
| **Mandatory SES adjustment** | Education, income, and occupation are treated as required socioeconomic covariates in MVMR. |
| **Multiple-testing correction** | Benjamini-Hochberg FDR correction and q-value reporting. |
| **Scalable output design** | Three comprehensive result tables for UVMR, MVMR, and mediation analysis. |
| **Reproducible local setup** | Local PLINK, 1000 Genomes EUR reference panel, and explicit directory-based file configuration. |

---

## 🚀 Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/Hexiao-DING/Dementia_Depression_MR-analysis.git
cd Dementia_Depression_MR-analysis
```

### 2. Install required R packages

```r
source("Install_All_Packages.R")
```

### 3. Place required data files

| Data type | Local directory |
|:--|:--|
| Exposure pQTLs / GWAS | `Standardized Circulating human plasma proteome_Data/` |
| Mediator mQTLs | `Cerebrospinal fluid metabolomics_Data/` |
| Outcome GWAS | `Outcomes/` |
| SES covariates | `Covariates_SES/` |
| LD reference panel | `LD_Reference/` |
| MRlap reference files | `MRlap_Reference/` |

### 4. Configure paths

Update data, PLINK, LD-reference, and MRlap-reference paths in `Main analysis.R`, especially the path block around lines 61–72.

### 5. Run the full analysis

```r
source("Main analysis.R")
```

### 6. Generate filtered summaries

```r
source("Results_Filter_Helper.R")
generate_summary_report()
```

---

## 🗂️ Data Atlas — Data are European Ancestry

### A. Exposure Layer — pQTLs and GWAS Summary Statistics

| Exposure group | Required files | Preprocessing / strategy | Source | Citation |
|:--|:--|:--|:--|:--|
| **Standardized circulating human plasma proteome data** | `GSCT005806_GRCh37.tsv.gz`; `GCST90240120_GRCh37.tsv.gz` → `GCST90243401_GRCh37.tsv.gz` | Standardize pQTLs into a homogeneous format using `MungeSumstats version 1.14.1`. | [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/29875488) | Sun, B. B. et al. *Nature*. 2018;558:73–79. |
| **Standardized circulating metabolic biomarkers data** | `GCST90301941.tsv` → `GCST90302173.tsv` | Standardize GWAS summary statistics into a homogeneous format using `MungeSumstats version 1.14.1`. | [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/38448586) | Karjalainen, M. K. et al. *Nature*. 2024;628:130–138. |
| **Circulating inflammatory proteins data** | `GCST90274758.tsv.gz` → `GCST90274848.tsv.gz` | Used as inflammatory-protein pQTL exposure files. | [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/37563310) | Zhao, J. H. et al. *Nature Immunology*. 2023;24:1540–1551. |

> [!TIP]
> The file name `GSCT005806_GRCh37.tsv.gz` is kept exactly as provided in the original README. When downloading from GWAS Catalog, verify the corresponding accession and local file name before running the scripts.

### B. Mediator Layer — CSF Metabolomics / metabQTLs

| Mediator group | Required files | Source | Citation |
|:--|:--|:--|:--|
| **Cerebrospinal fluid metabolomics data** | `GCST90025999_buildGRCh37.tsv.gz` → `GCST90026336_buildGRCh37.tsv.gz` | [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33437055) | Panyard, D. J. et al. *Communications Biology*. 2021;4:63. |

### C. Outcome Layer — Dementia and Cognitive Traits

| Trait | Local file | Source | Citation |
|:--|:--|:--|:--|
| **Dementia** | `Finn-b-F5_DEMENTIA.tsv.gz` | [FinnGen F5_DEMENTIA](https://r10.risteys.finngen.fi/endpoints/F5_DEMENTIA) | Kurki, M. I. et al. *Nature*. 2023. |
| **Alzheimer's disease** | `GCST90012877.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33589840) | Schwartzentruber, J. et al. *Nature Genetics*. 2021. |
| **Cognitive performance** | `GCST006572.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/30038396) | Lee, J. J. et al. *Nature Genetics*. 2018. |
| **Vascular dementia** | `Finn-b-F5_VASCDEM.tsv.gz` | [FinnGen F5_VASCDEM](https://r9.risteys.finngen.fi/endpoints/F5_VASCDEM) | Kurki, M. I. et al. *Nature*. 2023. |
| **Lewy body dementia** | `GCST90001390.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33589841) | Chia, R. et al. *Nature Genetics*. 2021. |
| **Frontotemporal dementia** | `GCST90558311.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/40280976) | Pottier, C. et al. *Nature Communications*. 2025. |
| **Undefined dementia** | `Finn_b_F5_DEMNAS.tsv.gz` | [FinnGen F5_DEMNAS](https://r12.risteys.finngen.fi/endpoints/F5_DEMNAS) | Kurki, M. I. et al. *Nature*. 2023. |

### D. Outcome Layer — Depression-Related Traits

| Trait | Local file | Source | Citation |
|:--|:--|:--|:--|
| **Depressive disorder** | `GCST90476922.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/39024449) | Verma, A. et al. *Science*. 2024. |
| **Major depressive disorder** | `GCST90468123.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/39789286) | Loya, N. et al. *Nature Genetics*. 2025. |
| **Mixed anxiety and depression** | `GCST90225526.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/37259642) | Brasher, M. S. et al. *Genes, Brain and Behavior*. 2023. |

### E. Covariate Layer — Socioeconomic Status

| SES covariate | Local file | Source | Citation | Required for MVMR |
|:--|:--|:--|:--|:--:|
| **Education** | `Education_GCST003676.txt.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/27225129) | Okbay, A. et al. *Nature*. 2016;533:539–542. | ✅ |
| **Income** | `Income_GCST90566700.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/40360725) | Xia, C. et al. *Molecular Psychiatry*. 2025. | ✅ |
| **Occupation** | `Occupation_GCST90566702.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/40360725) | Xia, C. et al. *Molecular Psychiatry*. 2025. | ✅ |

### F. Local Reference Data

| Resource | Required files / directories | Purpose | Source |
|:--|:--|:--|:--|
| **PLINK 1.9** | Executable available in system `PATH`, or full path specified in scripts. | Local LD clumping. | [PLINK website](https://www.cog-genomics.org/plink/) |
| **1000 Genomes EUR reference panel** | `1000G.EUR.QC.bed`; `1000G.EUR.QC.bim`; `1000G.EUR.QC.fam` | EUR-ancestry LD reference for clumping. | [PLINK Resources](https://www.cog-genomics.org/plink/1.9/resources); [1000 Genomes Project](https://www.internationalgenome.org/data/) |
| **MRlap LD scores** | `1000G_Phase3_ldscores/` | LD score files required by MRlap. | [MRlap GitHub](https://github.com/n-mounier/MRlap); [LDSC Resources](https://alkesgroup.broadinstitute.org/LDSCORE/) |
| **MRlap HapMap3 files** | `1000G_Phase3_hm3/` | HapMap3 variant information required by MRlap. | [MRlap GitHub](https://github.com/n-mounier/MRlap) |

---

## 🧾 Required File Schema

The pipeline accepts flexible column aliases, but each summary-statistics file should contain the following core fields.

| Field type | Accepted column names |
|:--|:--|
| **Variant ID** | `SNP`, `rsid`, `rs_id`, `MarkerName`, `variant` |
| **Effect allele** | `effect_allele`, `EA`, `A1` |
| **Other allele** | `other_allele`, `OA`, `A2` |
| **Effect size** | `beta`, `BETA`, `b`, `Beta` |
| **Standard error** | `se`, `SE`, `standard_error` |
| **P-value** | `pval`, `P`, `p_value`, `Pval` |

Recommended optional fields:

| Field type | Accepted column names |
|:--|:--|
| **Chromosome** | `chr`, `CHR`, `chromosome` |
| **Position** | `pos`, `BP`, `base_pair_location`, `POS` |
| **Sample size** | `n`, `N`, `samplesize` |
| **Effect allele frequency** | `eaf`, `EAF`, `freq` |

Supported file formats: `.tsv`, `.tsv.gz`, `.txt`, `.txt.gz`.

---

## 🎯 Instrument Selection

| Parameter | Default |
|:--|:--|
| **Genome-wide significance threshold** | `P < 5 × 10⁻⁸` |
| **LD clumping threshold** | `r² < 0.001` |
| **LD clumping window** | `10,000 kb` |
| **Minimum SNP count** | `≥ 3 SNPs` per analysis |
| **Weak instrument rule** | Prefer `F-statistic > 10` |

For MVMR, the combined instrument set follows the project rule:

> Genetic instruments are selected as the union of SNPs reaching genome-wide significance in **either** the exposure GWAS **or** the mediator GWAS.

Implementation: `select_combined_ivs_mvmr()`.

---

## 🔬 Methodological Core

### Reference standard

This pipeline follows the methodological logic of:

> Ye, C. J., Liu, D., Chen, M. L., *et al*. **Mendelian randomization evidence for the causal effect of mental well-being on healthy aging.** *Nature Human Behaviour*. 2024;8:1798–1809. DOI: `10.1038/s41562-024-01905-9`.
>
> GitHub: [https://github.com/yechaojie/mental_aging](https://github.com/yechaojie/mental_aging)

### Implementation checklist

| Domain | Component | Status | Notes |
|:--|:--|:--:|:--|
| **Causal estimation** | UVMR | ✅ | IVW, Weighted Median, MR-Egger, Weighted Mode, Simple Mode. |
| **Causal estimation** | MVMR | ✅ | IVW, Median, Egger, Lasso. |
| **Mediation** | Exposure → mediator → outcome decomposition | ✅ | Direct effect, indirect effect, and mediation proportion with 95% CI. |
| **Instrument strength** | F-statistics | ✅ | Used for UVMR instrument strength screening. |
| **Instrument strength** | Conditional F-statistics | ✅ | Used for MVMR instrument strength assessment. |
| **Horizontal pleiotropy** | MR-PRESSO | ✅ | Outlier detection and correction. |
| **Heterogeneity** | Cochran's Q | ✅ | IVW heterogeneity testing. |
| **Reverse causality** | Reverse MR | ✅ | Bidirectionality screening. |
| **Pathway logic** | Direction concordance | ✅ | Ensures exposure–mediator–outcome direction consistency. |
| **Sensitivity consistency** | Method concordance | ✅ | Validates IVW signals against sensitivity methods. |
| **Sample overlap** | MRlap | ✅ | Corrects IVW bias under sample overlap, weak instruments, and winner's curse. |
| **Confounding control** | SES adjustment | ✅ | Education, income, and occupation required in MVMR. |
| **Multiple testing** | FDR correction | ✅ | Benjamini-Hochberg and q-value reporting. |

**Total compliance:** ✅ **15/15 core requirements**.

---

## 💻 Installation

### Prerequisites

| Requirement | Recommendation |
|:--|:--|
| **R** | `≥ 4.0.0` |
| **RAM** | `8+ GB`; `16+ GB` recommended for MRlap. |
| **Internet connection** | Required for package installation and MRlap reference data. |
| **Disk space** | Sufficient space for GWAS files, reference panels, intermediate files, and results. |
| **PLINK** | PLINK 1.9 configured locally. |

### Automated installation

```r
source("Install_All_Packages.R")
```

### Manual installation

```r
install.packages(c(
  "data.table", "readr", "tibble",
  "TwoSampleMR", "ieugwasr",
  "MVMR", "MRPRESSO"
))

if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("n-mounier/MRlap")
```

### Verify installation

```r
source("Install_All_Packages.R")
source("Test_Covariate_Loading.R")
source("Demo_Test_Analysis.R")
```

---

## 🧭 Project Structure

```text
.
├── Main analysis.R                       # Main analysis script
├── Install_All_Packages.R                # Package installation
├── MungeSumstats_MR_data preprocessing.R # Data preprocessing
├── Results_Filter_Helper.R               # Result filtering utilities
├── Demo_Test_Analysis.R                  # Demo/test analysis
├── Test_Covariate_Loading.R              # Covariate testing
├── START_HERE_v2.5.txt                   # Getting started guide
│
├── Standardized Circulating human plasma proteome_Data/
│   └── [Exposure GWAS / pQTL files]
│
├── Cerebrospinal fluid metabolomics_Data/
│   └── [Mediator mQTL files]
│
├── Outcomes/
│   └── [Outcome GWAS files]
│
├── Covariates_SES/
│   ├── Education_GWAS.tsv
│   ├── Income_GWAS.tsv
│   └── Occupation_GWAS.tsv
│
├── LD_Reference/
│   ├── 1000G.EUR.QC.bed
│   ├── 1000G.EUR.QC.bim
│   └── 1000G.EUR.QC.fam
│
├── MRlap_Reference/
│   ├── 1000G_Phase3_ldscores/
│   └── 1000G_Phase3_hm3/
│
└── results_trial/
    ├── uvmr_comprehensive_results.csv
    ├── mvmr_comprehensive_results.csv
    └── mediation_comprehensive_results.csv
```

---

## ▶️ Running the Analysis

### 1. Data preprocessing

```r
source("MungeSumstats_MR_data preprocessing.R")
```

Purpose: harmonize GWAS/pQTL/mQTL summary statistics into a consistent MR-ready format.

### 2. Main causal analysis

```r
source("Main analysis.R")
```

Purpose: run the complete UVMR, MVMR, mediation, reverse MR, sensitivity, and correction modules.

### 3. Result filtering and summary report

```r
source("Results_Filter_Helper.R")
generate_summary_report()
```

Optional custom filtering:

```r
library(data.table)

uvmr <- fread("results_trial/uvmr_comprehensive_results.csv")

robust <- uvmr[
  grepl("IVW", method) &
    q_value < 0.05 &
    concordance_validated == TRUE &
    F_statistic > 10
]
```

---

## 📊 Output Files and Key Columns

All primary results are saved in `results_trial/`.

### Main result tables

| File | Approximate columns | Contents | Main use |
|:--|:--:|:--|:--|
| `uvmr_comprehensive_results.csv` | ~35 | All exposure–outcome pairs, five MR methods, MRlap correction, FDR q-values, and method-concordance validation. | Identify total causal associations. |
| `mvmr_comprehensive_results.csv` | ~20 | Covariate-adjusted independent effects, four MVMR methods, SES adjustment information, and FDR q-values. | Identify independent associations after covariate adjustment. |
| `mediation_comprehensive_results.csv` | ~32 | Exposure–mediator–outcome paths, reverse MR results, direction concordance, mediation proportions, and 95% CI. | Prioritise mechanistic pathways. |

### UVMR columns

| Column | Meaning |
|:--|:--|
| `b`, `se`, `pval` | Standard MR effect estimate, standard error, and P-value. |
| `q_value` | FDR-corrected P-value. |
| `concordance_validated` | Whether the IVW result is supported by sensitivity methods. |
| `F_statistic` | Instrument strength metric. |
| `mrlap_corrected_effect` | MRlap-corrected effect estimate. |
| `mrlap_corrected_pval` | MRlap-corrected P-value. |
| `mrlap_diff_pval` | Significance of difference between standard and corrected estimates. |

### MVMR columns

| Column | Meaning |
|:--|:--|
| `beta_mvmr`, `se_mvmr`, `pval_mvmr` | MVMR effect estimate, standard error, and P-value. |
| `q_value` | FDR-corrected P-value. |
| `adjusted_for` | Must include SES covariates. |
| `n_covariates` | Must be `3`: education, income, and occupation. |
| `method` | MVMR method, including IVW, Median, Egger, or Lasso. |
| `conditional_F_statistic` | MVMR instrument strength metric. |

### Mediation columns

| Column | Meaning |
|:--|:--|
| `beta_exp_med` | Exposure → mediator effect. |
| `beta_med_out_direct` | Mediator → outcome direct effect. |
| `beta_exp_out_total` | Exposure → outcome total effect. |
| `bidirectional` | Must be `No_Unidirectional` for a valid unidirectional pathway. |
| `direction_concordant` | Must be `TRUE` for direction-consistent mediation. |
| `mediation_proportion` | Estimated mediated proportion with 95% CI. |

---

## 🛡️ Quality Control and Result Filtering

### UVMR — robust causal associations

```r
robust_uvmr <- uvmr[
  grepl("IVW", method) &
    q_value < 0.05 &
    concordance_validated == TRUE &
    F_statistic > 10 &
    (is.na(Q_pval) | Q_pval > 0.05)
]
```

Interpretation: IVW estimates are treated as robust only when they are FDR-significant, instrument strength is acceptable, heterogeneity is not evident, and at least one sensitivity method supports the same direction and significance pattern.

### MVMR — independent effects after SES adjustment

```r
independent_mvmr <- mvmr[
  method == "MVMR-IVW" &
    q_value < 0.05 &
    n_covariates == 3 &
    grepl("Education", adjusted_for) &
    grepl("Income", adjusted_for) &
    grepl("Occupation", adjusted_for)
]
```

Interpretation: an MVMR result should not be interpreted as an SES-adjusted independent effect unless all three SES covariates are present.

### Mediation — valid pathways

```r
valid_mediation <- mediation[
  pval_exp_med < 0.05 &
    pval_med_out_direct < 0.05 &
    bidirectional == "No_Unidirectional" &
    direction_concordant == TRUE
]
```

A mediation pathway is considered valid only when all four criteria are satisfied:

- Exposure → mediator association is significant.
- Mediator → outcome direct association is significant.
- Reverse MR does not indicate bidirectionality.
- Direction concordance is preserved across the pathway.

### QC checklist

- [ ] F-statistics > 10 for all instruments.
- [ ] Cochran's Q P-value > 0.05 when requiring no heterogeneity.
- [ ] MR-PRESSO outlier test passed.
- [ ] Reverse MR shows no bidirectionality.
- [ ] Method concordance validated.
- [ ] MVMR adjusted for education, income, and occupation.
- [ ] FDR q-value < 0.05 for significant results.
- [ ] MRlap correction applied when sample overlap is detected or suspected.

---

## 💡 Advanced Features

### MRlap sample-overlap correction

MRlap provides bias correction for inverse-variance-weighted Mendelian randomization when exposure and outcome GWAS may have sample overlap. In this pipeline, MRlap is used to support correction for sample overlap, weak instruments, and winner's curse.

| Required MRlap asset | Local path | Purpose |
|:--|:--|:--|
| LD score files | `MRlap_Reference/1000G_Phase3_ldscores/` | LD score input for MRlap. |
| HapMap3 files | `MRlap_Reference/1000G_Phase3_hm3/` | HapMap3 variant filtering and harmonisation. |

Configuration example in `Main analysis.R`:

```r
mrlap_ldscores_path <- "path/to/1000G_Phase3_ldscores"
mrlap_hm3_path <- "path/to/1000G_Phase3_hm3"
```

MRlap auto-detection features in this pipeline:

- Extracts sample sizes from GWAS files.
- Estimates sample overlap using LDSC-related information.
- Corrects for winner's curse and weak instruments.
- Reduces the need for manual overlap configuration.

### Local PLINK and 1000 Genomes EUR LD clumping

The pipeline uses local PLINK for LD clumping with a European reference panel.

```r
ld_ref_path <- "path/to/1000G.EUR.QC"  # without .bed/.bim/.fam extension
```

Required files:

| File | Description |
|:--|:--|
| `1000G.EUR.QC.bed` | Binary genotype file. |
| `1000G.EUR.QC.bim` | SNP information file. |
| `1000G.EUR.QC.fam` | Sample information file. |

---

## 🔧 Troubleshooting

| Issue | Likely cause | Recommended solution |
|:--|:--|:--|
| **PLINK not found** | PLINK is not installed or not in `PATH`. | Install PLINK 1.9 and add it to `PATH`, or specify the full executable path in the scripts. |
| **1000G.EUR.QC files missing** | LD reference panel is absent or misconfigured. | Download `.bed`, `.bim`, and `.fam` files and update the LD reference path. |
| **MRlap reference data missing** | LD score or HapMap3 reference folder is absent. | Download `1000G_Phase3_ldscores/` and `1000G_Phase3_hm3/` and update paths in `Main analysis.R`. |
| **MRlap fails to detect sample sizes** | GWAS files lack `N` or case-control sample size columns. | Ensure files contain `N`, `N_cases`, and/or `N_controls` where required. |
| **MVMR fails with SES covariates** | SES files are missing or misnamed. | Verify education, income, and occupation files are present in `Covariates_SES/`. |
| **No instruments found** | Threshold too strict or column mapping failed. | Check `P < 5e-8`, LD clumping parameters, and P-value column names. |
| **MR-PRESSO errors** | Too few instruments. | Ensure at least four SNPs are available for MR-PRESSO. |
| **Memory issues with MRlap** | Large GWAS files or insufficient RAM. | Increase RAM allocation, use a machine with more memory, or run smaller batches. |
| **F-statistics < 10** | Weak instruments. | Exclude weak instruments or reassess instrument-selection thresholds. |

Getting help:

1. Read `START_HERE_v2.5.txt`.
2. Review comments in the analysis scripts.
3. Run `Demo_Test_Analysis.R` before the full analysis.
4. Open an issue on [GitHub](https://github.com/Hexiao-DING/Dementia_Depression_MR-analysis/issues).

---

## 📚 Citation and References

### If you use this pipeline

Please cite the methodological references and the original GWAS/QTL resources used in your analysis.

### Core MR methods

| Method / resource | Citation |
|:--|:--|
| **TwoSampleMR / MR-Base** | Hemani, G., Zheng, J., Elsworth, B., *et al*. The MR-Base platform supports systematic causal inference across the human phenome. *eLife*. 2018;7:e34408. DOI: `10.7554/eLife.34408`. |
| **MVMR** | Sanderson, E., Glymour, M. M., Holmes, M. V., *et al*. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. *International Journal of Epidemiology*. 2019;48:713–727. DOI: `10.1093/ije/dyy262`. |
| **MR-PRESSO** | Verbanck, M., Chen, C. Y., Neale, B., & Do, R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases. *Nature Genetics*. 2018;50:693–698. DOI: `10.1038/s41588-018-0099-7`. |
| **MRlap** | Mounier, N. & Kutalik, Z. Bias correction for inverse variance weighting Mendelian randomization. *Genetic Epidemiology*. 2023;47:314–331. DOI: `10.1002/gepi.22522`. GitHub: [https://github.com/n-mounier/MRlap](https://github.com/n-mounier/MRlap). |
| **Reference standard** | Ye, C. J., Liu, D., Chen, M. L., *et al*. Mendelian randomization evidence for the causal effect of mental well-being on healthy aging. *Nature Human Behaviour*. 2024;8:1798–1809. DOI: `10.1038/s41562-024-01905-9`. GitHub: [https://github.com/yechaojie/mental_aging](https://github.com/yechaojie/mental_aging). |

### Data standardisation and reference resources

| Resource | Citation / link |
|:--|:--|
| **MungeSumstats** | Murphy, A. E., Schilder, B. M., Skene, N. G., *et al*. MungeSumstats: a Bioconductor package for the standardization and quality control of many GWAS summary statistics. *Bioinformatics*. 2021;37:4593–4596. DOI: `10.1093/bioinformatics/btab665`. |
| **PLINK** | Purcell, S., Neale, B., Todd-Brown, K., *et al*. PLINK: a tool set for whole-genome association and population-based linkage analyses. *American Journal of Human Genetics*. 2007;81:559–575. |
| **1000 Genomes Project** | The 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature*. 2015;526:68–74. |
| **Benjamini-Hochberg FDR** | Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*. 1995;57:289–300. |

### GWAS / QTL data sources

| Data layer | Citation |
|:--|:--|
| **Plasma proteome pQTLs** | Sun, B. B., Maranville, J. C., Peters, J. E., Stacey, D., Staley, J. R., Blackshaw, J., *et al*. Genomic atlas of the human plasma proteome. *Nature*. 2018;558:73–79. |
| **Circulating metabolic biomarkers** | Karjalainen, M. K., Karthikeyan, S., Oliver-Williams, C., Sliz, E., Allara, E., Fung, W. T., *et al*. Genome-wide characterization of circulating metabolic biomarkers. *Nature*. 2024;628:130–138. |
| **Inflammatory proteins pQTLs** | Zhao, J. H., Stacey, D., Eriksson, N., Macdonald-Dunlop, E., Hedman, Å. K., Kalnapenkis, A., *et al*. Genetics of circulating inflammatory proteins identifies drivers of immune-mediated disease risk and therapeutic targets. *Nature Immunology*. 2023;24:1540–1551. |
| **CSF metabolomics / metabQTLs** | Panyard, D. J., Kim, K. M., Darst, B. F., Deming, Y. K., Zhong, X., Wu, Y., *et al*. Cerebrospinal fluid metabolomics identifies 19 brain-related phenotype associations. *Communications Biology*. 2021;4:63. |
| **FinnGen dementia endpoints** | Kurki, M. I., Karjalainen, J., Palta, P., Sipilä, T. P., Kristiansson, K., Donner, K. M., *et al*. FinnGen provides genetic insights from a well-phenotyped isolated population. *Nature*. 2023. |
| **Alzheimer's disease** | Schwartzentruber, J., Cooper, S., Liu, J. Z., Barrio-Hernandez, I., Bello, E., Kumasaka, N., *et al*. Genome-wide meta-analysis, fine-mapping and integrative prioritization implicate new Alzheimer's disease risk genes. *Nature Genetics*. 2021. |
| **Cognitive performance** | Lee, J. J., Wedow, R., Okbay, A., Kong, E., Maghzian, O., Zacher, M., *et al*. Gene discovery and polygenic prediction from a genome-wide association study of educational attainment in 1.1 million individuals. *Nature Genetics*. 2018. |
| **Lewy body dementia** | Chia, R., Sabir, M. S., Bandres-Ciga, S., Saez-Atienzar, S., Reynolds, R. H., Gustavsson, E., *et al*. Genome sequencing analysis identifies new loci associated with Lewy body dementia and provides insights into its genetic architecture. *Nature Genetics*. 2021. |
| **Frontotemporal dementia** | Pottier, C. et al. Frontotemporal dementia GWAS resource associated with `GCST90558311`. *Nature Communications*. 2025. |
| **Depressive disorder** | Verma, A. et al. Depression-related GWAS resource associated with `GCST90476922`. *Science*. 2024. |
| **Major depressive disorder** | Loya, N. et al. Major depressive disorder GWAS resource associated with `GCST90468123`. *Nature Genetics*. 2025. |
| **Mixed anxiety and depression** | Brasher, M. S. et al. Mixed anxiety/depression GWAS resource associated with `GCST90225526`. *Genes, Brain and Behavior*. 2023. |
| **Educational attainment** | Okbay, A., Beauchamp, J. P., Fontana, M. A., Lee, J. J., Pers, T. H., Rietveld, C. A., *et al*. Genome-wide association study identifies 74 loci associated with educational attainment. *Nature*. 2016;533:539–542. |
| **Income and occupation SES GWAS** | Xia, C. et al. Deciphering the influence of socioeconomic status on brain structure: insights from Mendelian randomization. *Molecular Psychiatry*. 2025. |

---

## 🙏 Acknowledgements

| Role | Name | Affiliation |
|:--|:--|:--|
| **Contributor** | Hexiao Ding | The Hong Kong Polytechnic University |
| **Contributor** | Jing Lan | The Hong Kong Polytechnic University |
| **Contributor** | Na Li | Sichuan University |
| **Supervisor** | Dr. Jung Sun Yoo | The Hong Kong Polytechnic University |

**Institutions**

- The Hong Kong Polytechnic University
- Sichuan University

**Version:** 2.7  
**Last Updated:** 2025-12-10

### Methodological acknowledgements

This pipeline implements and integrates methods from:

- **TwoSampleMR / MR-Base** for two-sample Mendelian randomization.
- **MVMR** for multivariable Mendelian randomization.
- **MR-PRESSO** for horizontal pleiotropy and outlier detection.
- **MRlap** for sample-overlap, weak-instrument, and winner's-curse correction.
- **MungeSumstats** for GWAS summary-statistics standardisation.

---

<div align="center">

**Designed for reproducible molecular epidemiology, neuropsychiatric genetics, and dementia–depression causal pathway prioritisation.**

</div>
