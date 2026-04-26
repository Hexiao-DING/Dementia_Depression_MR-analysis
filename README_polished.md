<!--
README polished for GitHub presentation.
Design goal: keep all original scientific, technical, data-source, path, output, QC, citation, and acknowledgement information, while adding stronger visual hierarchy, navigation, reproducibility notes, and supplementary methodological references.
-->

<div align="center">

# 🧬 Dementia–Depression Mendelian Randomization Analysis Pipeline

### A reproducible causal-inference workflow for dementia–depression comorbidity, peripheral molecular traits, CSF metabolites, and SES-adjusted genetic evidence

<p>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://www.r-project.org/"><img src="https://img.shields.io/badge/R-%3E%3D%204.0-276DC3?logo=r&logoColor=white" alt="R >= 4.0"></a>
  <img src="https://img.shields.io/badge/Platform-Windows%20%7C%20Linux-lightgrey" alt="Platform: Windows | Linux">
  <img src="https://img.shields.io/badge/Status-Production%20Ready-success" alt="Status: Production Ready">
  <img src="https://img.shields.io/badge/Ancestry-Mostly%20European-blueviolet" alt="Most data are European ancestry">
  <img src="https://img.shields.io/badge/Workflow-UVMR%20%7C%20MVMR%20%7C%20Mediation-0A7E8C" alt="Workflow: UVMR | MVMR | Mediation">
</p>

<p>
  <a href="#-overview">Overview</a> ·
  <a href="#-quick-start">Quick Start</a> ·
  <a href="#-data-requirements-most-data-are-eu-populations">Data</a> ·
  <a href="#-methodological-compliance">Methods</a> ·
  <a href="#-quality-control">Quality Control</a> ·
  <a href="#-citation-and-references">Citation</a>
</p>

<img width="1846" height="1041" alt="Dementia–Depression MR analysis framework" src="https://github.com/user-attachments/assets/4137db74-5fbc-40d0-8f71-a9d91a9aed03" />

</div>

---

## 🌟 At a Glance

| Dimension | What this pipeline provides |
|:--|:--|
| **Scientific question** | Genetic causal evidence for dementia–depression comorbidity and molecular mediating pathways |
| **Exposure layer** | Circulating human plasma proteome, circulating metabolic biomarkers, circulating inflammatory proteins |
| **Mediator layer** | Cerebrospinal fluid metabolomics / metabQTLs |
| **Outcome layer** | Dementia subtypes, cognitive performance, depression-related traits |
| **Covariate layer** | Socioeconomic status: education, income, occupation |
| **Core methods** | UVMR, MVMR, mediation MR, reverse MR, MRlap correction, MR-PRESSO, heterogeneity testing, FDR correction |
| **Main output** | `uvmr_comprehensive_results.csv`, `mvmr_comprehensive_results.csv`, `mediation_comprehensive_results.csv` |
| **Primary use case** | Prioritising peripheral-central molecular pathways and causal candidates in dementia–depression comorbidity |

> **Data-access note**  
> GWAS, pQTL, mQTL, LD reference, and MRlap reference files are **not included** in this repository. Users must download the required files from the listed public resources and configure local paths before running the full pipeline.

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Scientific Rationale](#-scientific-rationale)
- [Features](#-features)
- [Data Requirements — Most data are EU populations](#-data-requirements-most-data-are-eu-populations)
- [File Format Specification](#file-format-specification)
- [IV Selection Parameters](#-iv-selection-parameters)
- [Output Files](#-output-files)
- [Methodological Compliance](#-methodological-compliance)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Usage](#-usage)
- [Output Structure](#-output-structure)
- [Quality Control](#-quality-control)
- [Advanced Features](#-advanced-features)
- [Troubleshooting](#-troubleshooting)
- [Citation and References](#-citation-and-references)
- [Acknowledgements](#-acknowledgements)

---

## 🎯 Overview

This pipeline provides a **comprehensive, methodologically rigorous** framework for Mendelian Randomization (MR) analysis in dementia–depression comorbidity research. It implements univariable MR (UVMR), multivariable MR (MVMR), and mediation analysis following best practices established in high-impact publications.

The workflow is designed for large-scale summary-statistics analyses linking **peripheral molecular exposures**, **central nervous system metabolic mediators**, **neuropsychiatric outcomes**, and **socioeconomic covariates**. It supports a layered causal architecture in which circulating proteins, metabolites, and inflammatory factors are evaluated for their potential causal relevance to dementia and depression phenotypes, while mediation and SES-adjusted analyses are used to strengthen biological interpretability.

---

## 🧠 Scientific Rationale

Dementia and depression frequently co-occur across the life course and may share inflammatory, metabolic, vascular, and neurobiological mechanisms. Observational associations are vulnerable to confounding, reverse causality, disease-stage effects, and selection bias. Mendelian randomization uses genetic variants as instrumental variables to estimate causal effects from GWAS summary statistics, thereby providing a complementary line of evidence for prioritising mechanistic candidates and potential intervention targets.

This repository focuses on three connected questions:

1. **Peripheral molecular drivers** — whether circulating proteins, metabolites, or inflammatory proteins show genetic evidence of causal relevance to dementia- or depression-related outcomes.
2. **Central mediating biology** — whether CSF metabolomic traits statistically mediate exposure–outcome pathways.
3. **SES-independent evidence** — whether putative effects remain after multivariable adjustment for education, income, and occupation.

---

## 🧭 Analysis Framework

```mermaid
flowchart TD
    A[GWAS / pQTL / mQTL summary statistics] --> B[Data standardisation with MungeSumstats]
    B --> C[Instrument selection]
    C --> C1[Genome-wide significant SNPs\nP < 5 × 10^-8]
    C --> C2[LD clumping\nr² < 0.001, 10,000 kb]
    C --> C3[F-statistic > 10]

    C --> D[UVMR]
    C --> E[MVMR]

    D --> D1[IVW]
    D --> D2[Weighted median]
    D --> D3[MR-Egger]
    D --> D4[Weighted mode]
    D --> D5[Simple mode]
    D --> D6[MRlap correction]

    E --> E1[MVMR-IVW]
    E --> E2[MVMR median]
    E --> E3[MVMR Egger]
    E --> E4[MVMR Lasso]
    E --> E5[SES-adjusted effects]
    E --> E6[Conditional F-statistics]

    D --> F[Sensitivity analyses]
    E --> F
    F --> F1[MR-PRESSO]
    F --> F2[Cochran's Q]
    F --> F3[Reverse MR]
    F --> F4[Method concordance]

    F --> G[Mediation analysis]
    G --> G1[Exposure → Mediator]
    G --> G2[Mediator → Outcome]
    G --> G3[Direct / indirect effects]
    G --> G4[Mediation proportion]
    G --> G5[Direction concordance]

    G --> H[Multiple testing correction]
    H --> H1[FDR / Benjamini-Hochberg]
    H --> H2[q-values]
    H --> I[Comprehensive result tables]
```

<details>
<summary><strong>Text-based framework fallback</strong></summary>

```text
┌───────────────────────────────────────────────────────────┐
│                   MR Analysis Pipeline                    │
└───────────────────────────────────────────────────────────┘
        ┌───────────────────────────────────────┐
        │  Step 1: Instrument Selection         │
        │  • Genome-wide significant SNPs       │
        │  • LD clumping                        │
        │  • F-statistic > 10                   │
        └───────────────┬───────────────────────┘
                        │
        ┌───────────────┴───────────────┐
        │                               │
        ▼                               ▼
┌──────────────────┐          ┌────────────────────┐
│  UVMR Analysis   │          │  MVMR Analysis     │
│  • IVW           │          │  • IVW             │
│  • Weighted      │          │  • Median          │
│    Median        │          │  • Egger           │
│  • MR-Egger      │          │  • Lasso           │
│  • Weighted      │          │  • SES-adjusted    │
│    Mode          │          │  • Conditional F   │
│  • Simple Mode   │          └──────────┬─────────┘
│  • MRlap         │                     │
│    correction    │                     │
└────────┬─────────┘                     │
         │                               │
         └───────────────┬───────────────┘
                         │
                         ▼
        ┌───────────────────────────────────────┐
        │  Step 2: Sensitivity Analyses         │
        │  • MR-PRESSO outlier detection        │
        │  • Cochran's Q (heterogeneity)        │
        │  • Reverse MR (bidirectionality)      │
        │  • Method concordance validation      │
        └───────────────┬───────────────────────┘
                        │
                        ▼
        ┌───────────────────────────────────────┐
        │  Step 3: Mediation Analysis           │
        │  • Exposure → Mediator → Outcome      │
        │  • Direct vs indirect effects         │
        │  • Mediation proportion               │
        │  • Direction concordance check        │
        └───────────────┬───────────────────────┘
                        │
                        ▼
        ┌───────────────────────────────────────┐
        │  Step 4: Multiple Testing Correction  │
        │  • FDR (Benjamini-Hochberg)           │
        │  • q-values                           │
        └───────────────────────────────────────┘
```

</details>

---

## 🔬 Key Objectives

* 🔬 **Causal Inference**: Establish causal relationships between exposures and outcomes using genetic instruments.
* 📊 **Multivariable Analysis**: Adjust for confounders and mediators using MVMR.
* 🔄 **Mediation Analysis**: Identify intermediate pathways in causal chains.
* ✅ **Rigorous Validation**: Comprehensive sensitivity analyses and quality control.
* 📈 **Sample Overlap Correction**: Automatic correction using MRlap.
* 🎯 **SES Adjustment**: Mandatory socioeconomic status covariate adjustment.

---

## ✨ Features

| Feature | Included | Notes |
|:--|:--:|:--|
| 🔬 **Comprehensive MR Methods** | ✅ | UVMR with 5 methods and MVMR with 4 methods |
| 🛡️ **Rigorous QC** | ✅ | MR-PRESSO outlier detection, heterogeneity testing, reverse MR |
| 📊 **Sample Overlap Correction** | ✅ | Automatic MRlap implementation |
| 🎯 **SES Adjustment** | ✅ | Mandatory socioeconomic status covariate adjustment in MVMR |
| 🔄 **Mediation Analysis** | ✅ | Complete pathway decomposition |
| ✅ **Method Concordance** | ✅ | Validation across multiple MR methods |
| 📈 **FDR Correction** | ✅ | Benjamini-Hochberg multiple testing correction |
| 🔍 **Combined IV Selection** | ✅ | Standard-compliant instrument selection for MVMR |
| 📚 **Well-Documented** | ✅ | Comprehensive documentation and examples |

---

## 📋 Data Requirements [Most data are EU populations]

### 🔹 Required genome-wide association studies (GWAS) Summary Statistics, protein quantitative trait loci (pQTLs), and metabolite quantitative trait loci (mQTLs)

> ⚠️ **IMPORTANT:** You need to provide these files. They are **NOT included** in the repository.

> 🧬 **Ancestry note:** Most data sources are European-population datasets. The pipeline therefore uses an EUR LD reference panel by default. Downstream interpretation should avoid over-generalising results to non-European populations unless appropriate ancestry-matched data and LD references are used.

### Data Layer Map

| Layer | Data type | Biological meaning | Main folder |
|:--|:--|:--|:--|
| **Exposure** | pQTLs / GWAS | Peripheral proteins, metabolites, inflammatory proteins | `Standardized Circulating human plasma proteome_Data/` |
| **Mediator** | metabQTLs | CSF metabolomic traits | `Cerebrospinal fluid metabolomics_Data/` |
| **Outcome** | GWAS | Dementia, cognitive performance, depression traits | `Outcomes/` |
| **Covariate** | GWAS | Socioeconomic status indicators | `Covariates_SES/` |
| **LD reference** | 1000G EUR | Local LD clumping | `LD_Reference/` |
| **MRlap reference** | LD scores / HapMap3 | Sample-overlap correction | `MRlap_Reference/` |

---

### 1️⃣ Exposure pQTLs/GWAS — Multiple Exposure Files

#### (1) Standardized Circulating Human Plasma Proteome Data (pQTLs)

* **Files:** `GSCT005806_GRCh37.tsv.gz`, `GCST90240120_GRCh37.tsv.gz` → `GCST90243401_GRCh37.tsv.gz`
* **Strategy:** Standardize pQTLs into homogeneous format using `MungeSumstats version 1.14.1`.
* **Source:** [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/29875488)
* **Citation:**

    > Sun, B. B., Maranville, J. C., Peters, J. E., Stacey, D., Staley, J. R., Blackshaw, J., ... & Butterworth, A. S. (2018). Genomic atlas of the human plasma proteome. *Nature*, 558(7708), 73-79.

#### (2) Standardized Circulating Metabolic Biomarkers Data (GWAS)

* **Files:** `GCST90301941.tsv` → `GCST90302173.tsv`
* **Strategy:** Standardize GWAS summary statistics into homogeneous format using `MungeSumstats version 1.14.1`.
* **Source:** [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/38448586)
* **Citation:**

    > Karjalainen, M. K., Karthikeyan, S., Oliver-Williams, C., Sliz, E., Allara, E., Fung, W. T., ... & Kettunen, J. (2024). Genome-wide characterization of circulating metabolic biomarkers. *Nature*, 628(8006), 130-138.

#### (3) Circulating Inflammatory Proteins Data (pQTLs)

* **Files:** `GCST90274758.tsv.gz` → `GCST90274848.tsv.gz`
* **Source:** [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/37563310)
* **Citation:**

    > Zhao, J. H., Stacey, D., Eriksson, N., Macdonald-Dunlop, E., Hedman, Å. K., Kalnapenkis, A., ... & Peters, J. E. (2023). Genetics of circulating inflammatory proteins identifies drivers of immune-mediated disease risk and therapeutic targets. *Nature Immunology*, 24(9), 1540-1551.

---

### 2️⃣ Mediator mQTLs — Mediator Files

#### (1) Cerebrospinal Fluid Metabolomics Data (metabQTLs)

* **Files:** `GCST90025999_buildGRCh37.tsv.gz` → `GCST90026336_buildGRCh37.tsv.gz`
* **Source:** [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33437055)
* **Citation:**

    > Panyard, D. J., Kim, K. M., Darst, B. F., Deming, Y. K., Zhong, X., Wu, Y., ... & Lu, Q. (2021). Cerebrospinal fluid metabolomics identifies 19 brain-related phenotype associations. *Communications Biology*, 4(1), 63.

---

### 3️⃣ Outcome GWAS — Outcome Files

#### (1) Dementia Types *(7 traits in total, GWAS)*

| Trait | File | Source & Citation |
| :--- | :--- | :--- |
| **Dementia** | `Finn-b-F5_DEMENTIA.tsv.gz` | [FinnGen](https://r10.risteys.finngen.fi/endpoints/F5_DEMENTIA) — Kurki et al. (2023) *Nature* |
| **Alzheimer's** | `GCST90012877.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33589840) — Schwartzentruber et al. (2021) *Nature Genetics* |
| **Cognitive Perf.** | `GCST006572.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/30038396) — Lee et al. (2018) *Nature Genetics* |
| **Vascular Dem.** | `Finn-b-F5_VASCDEM.tsv.gz` | [FinnGen](https://r9.risteys.finngen.fi/endpoints/F5_VASCDEM) — Kurki et al. (2023) *Nature* |
| **Lewy Body Dem.** | `GCST90001390.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33589841) — Chia et al. (2021) *Nature Genetics* |
| **Frontotemporal** | `GCST90558311.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/40280976) — Pottier et al. (2025) *Nature Communications* |
| **Undefined Dem.** | `Finn_b_F5_DEMNAS.tsv.gz` | [FinnGen](https://r12.risteys.finngen.fi/endpoints/F5_DEMNAS) — Kurki et al. (2023) *Nature* |

#### (2) Depression Types *(3 traits in total, GWAS)*

| Trait | File | Source & Citation |
| :--- | :--- | :--- |
| **Depressive Dis.** | `GCST90476922.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/39024449) — Verma et al. (2024) *Science* |
| **Major DD** | `GCST90468123.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/39789286) — Loya et al. (2025) *Nature Genetics* |
| **Mixed Anxiety/Dep** | `GCST90225526.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/37259642) — Brasher et al. (2023) *Genes, Brain and Behavior* |

---

### 4️⃣ Covariate GWAS — Social Economic Status *(SES, 3 traits in total, GWAS)*

#### (1) Education

* **File:** `Education_GCST003676.txt.gz`
* **Source:** [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/27225129)
* **Citation:**

    > Okbay, A. et al. (2016). Genome-wide association study identifies 74 loci associated with educational attainment. *Nature*, 533(7604), 539-542.

#### (2) Income

* **File:** `Income_GCST90566700.tsv.gz`
* **Source:** [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/40360725)
* **Citation:**

    > Xia, C. et al. (2025). Deciphering the influence of socioeconomic status on brain structure: insights from Mendelian randomization. *Molecular Psychiatry*, 1-14.

#### (3) Occupation

* **File:** `Occupation_GCST90566702.tsv.gz`
* **Source:** [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/40360725)
* **Citation:**

    > Xia, C. et al. (2025). Deciphering the influence of socioeconomic status on brain structure: insights from Mendelian randomization. *Molecular Psychiatry*, 1-14.

---

### File Format Specification

**Required columns** with flexible naming:

| Concept | Accepted column names |
|:--|:--|
| SNP identifier | `SNP`, `rsid`, `rs_id`, `MarkerName`, `variant` |
| Effect allele | `effect_allele`, `EA`, `A1` |
| Other allele | `other_allele`, `OA`, `A2` |
| Effect estimate | `beta`, `BETA`, `b`, `Beta` |
| Standard error | `se`, `SE`, `standard_error` |
| P-value | `pval`, `P`, `p_value`, `Pval` |

**Optional but recommended**:

| Concept | Accepted column names |
|:--|:--|
| Chromosome | `chr`, `CHR`, `chromosome` |
| Position | `pos`, `BP`, `base_pair_location`, `POS` |
| Sample size | `n`, `N`, `samplesize` |
| Effect allele frequency | `eaf`, `EAF`, `freq` |

**Supported formats:** `.tsv`, `.tsv.gz`, `.txt`, `.txt.gz`

---

## 🎯 IV Selection Parameters

Following genome-wide standards:

| Parameter | Value | Purpose |
|:--|:--|:--|
| **P-value threshold** | `5×10⁻⁸` | Genome-wide significant instruments |
| **LD clumping** | `r² < 0.001`, window = `10,000 kb` | Reduce correlated instruments |
| **Minimum SNPs** | 3 per analysis | Ensure minimal instrument support |
| **Instrument strength** | F-statistic > 10 | Reduce weak-instrument bias |

---

## 📊 Output Files

### Main Results (3 CSV files)

1. **`uvmr_comprehensive_results.csv`** (~35 columns)

    * All exposure-outcome pairs
    * 5 MR methods per pair
    * MRlap sample overlap correction
    * FDR q-values
    * Method concordance validation

2. **`mvmr_comprehensive_results.csv`** (~20 columns)

    * Covariate-adjusted independent effects
    * 4 MVMR methods
    * SES adjustment information
    * FDR q-values

3. **`mediation_comprehensive_results.csv`** (~32 columns)

    * All exposure-mediator-outcome paths
    * Reverse MR results
    * Direction concordance
    * Mediation proportions with 95% CI

---

## 🔬 Methodological Compliance

### Reference Standard

This pipeline implements methods following the rigorous standards established in:

> **Ye CJ, Liu D, Chen ML, *et al*. Mendelian randomization evidence for the causal effect of mental well-being on healthy aging.**  
> *Nature Human Behaviour*. 2024.  
> DOI: [10.1038/s41562-024-01905-9](https://doi.org/10.1038/s41562-024-01905-9)  
> **GitHub**: [https://github.com/yechaojie/mental_aging](https://github.com/yechaojie/mental_aging)

### Implementation Checklist

| Component | Status | Description |
|:----------|:------:|:------------|
| **UVMR Methods** | ✅ | IVW, Weighted Median, MR-Egger, Weighted Mode, Simple Mode |
| **MVMR Methods** | ✅ | IVW, Median, Egger, Lasso |
| **MR-PRESSO** | ✅ | Outlier detection and correction |
| **MRlap** | ✅ | Sample overlap correction with automatic detection |
| **SES Covariate Adjustment** | ✅ | Education, Income, Occupation (mandatory for MVMR) |
| **Combined IV Selection** | ✅ | SNPs significant in either exposure or mediator GWAS |
| **FDR Correction** | ✅ | Benjamini-Hochberg multiple testing correction |
| **Method Concordance** | ✅ | Validation across sensitivity methods |
| **Reverse MR** | ✅ | Bidirectionality testing |
| **Direction Concordance** | ✅ | Pathway consistency validation |
| **F-statistics** | ✅ | Instrument strength assessment |
| **Cochran's Q** | ✅ | Heterogeneity testing |
| **Conditional F-statistics** | ✅ | MVMR instrument strength |
| **Mediation Analysis** | ✅ | Complete pathway decomposition |
| **Mediation Proportion** | ✅ | Indirect effect quantification with 95% CI |

**Total Compliance:** ✅ **15/15 Core Requirements** (100%)

---

## 💻 Installation

### Prerequisites

| Requirement | Minimum | Recommended |
|:--|:--|:--|
| **R** | ≥ 4.0.0 | Current stable R version |
| **RAM** | 8+ GB | 16+ GB for MRlap |
| **Internet Connection** | Required | Required for package installation and MRlap reference data |
| **Disk Space** | Depends on GWAS size | Sufficient space for GWAS files and results |

### PLINK and LD Reference Panel (EUR Population)

This pipeline uses **local PLINK** for LD clumping with European (EUR) population reference data.

> ⚠️ **IMPORTANT:** You must configure PLINK and the 1000 Genomes EUR reference panel locally before running the analysis.

#### Required Software

* **PLINK 1.9**: [Download from PLINK website](https://www.cog-genomics.org/plink/)
  * Add PLINK to your system PATH, or specify the full path in the analysis scripts.

#### Required Reference Files (1000G EUR)

Download and configure the following 1000 Genomes European population QC'd reference files:

| File | Description |
|:-----|:------------|
| `1000G.EUR.QC.bed` | Binary genotype file |
| `1000G.EUR.QC.bim` | SNP information file |
| `1000G.EUR.QC.fam` | Sample information file |

**Download Source:** These files can be obtained from:

* [PLINK Resources](https://www.cog-genomics.org/plink/1.9/resources)
* [1000 Genomes Project](https://www.internationalgenome.org/data/)

**Configuration:** Update the LD reference panel path in `Main analysis.R`:

```r
# Example configuration
ld_ref_path <- "path/to/1000G.EUR.QC"  # Without file extension
```

### Automated Installation (Recommended)

```bash
# Clone the repository
git clone https://github.com/Hexiao-DING/Dementia_Depression_MR-analysis.git

cd Dementia_Depression_MR-analysis

# Run automated installation script in R
Rscript -e 'source("Install_All_Packages.R")'
```

Alternative inside an R session:

```r
source("Install_All_Packages.R")
```

### Manual Installation

```r
# Install CRAN packages
install.packages(c(
  "data.table", "readr", "tibble", 
  "TwoSampleMR", "ieugwasr", 
  "MVMR", "MRPRESSO"
))

# Install MRlap from GitHub
if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

devtools::install_github("n-mounier/MRlap")
```

### Verify Installation

```r
# Test package loading
source("Install_All_Packages.R")

# Test covariate file loading
source("Test_Covariate_Loading.R")

# Run demo analysis before full analysis
source("Demo_Test_Analysis.R")
```

---

## 🚀 Quick Start

### Step 1: Prepare Data

1. **Place GWAS files** in designated directories:

    * Exposures: `Standardized Circulating human plasma proteome_Data/`
    * Mediators: `Cerebrospinal fluid metabolomics_Data/`
    * Outcomes: `Outcomes/`

2. **Download SES Covariates**. These are **REQUIRED for MVMR**:

    * Education
    * Income
    * Occupation
    * *Place in:* `Covariates_SES/`

3. **Configure paths** in `Main analysis.R`, especially lines 61-72.

### Step 2: Run Analysis

```r
# Full analysis; runtime may be hours depending on data size
source("Main analysis.R")
```

### Step 3: Filter Results

```r
# Automated filtering and summary report
source("Results_Filter_Helper.R")

generate_summary_report()
```

---

## 📚 Usage

### Project Structure

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
│   └── [Exposure GWAS files]
│
├── Cerebrospinal fluid metabolomics_Data/
│   └── [Mediator GWAS files]
│
├── Outcomes/
│   └── [Outcome GWAS files]
│
├── Covariates_SES/
│   ├── Education_GWAS.tsv
│   ├── Income_GWAS.tsv
│   └── Occupation_GWAS.tsv
│
├── LD_Reference/                         # LD reference panel
│   ├── 1000G.EUR.QC.bed
│   ├── 1000G.EUR.QC.bim
│   └── 1000G.EUR.QC.fam
│
├── MRlap_Reference/                      # MRlap reference data
│   ├── 1000G_Phase3_ldscores/
│   └── 1000G_Phase3_hm3/
│
└── results_trial/
    ├── uvmr_comprehensive_results.csv
    ├── mvmr_comprehensive_results.csv
    └── mediation_comprehensive_results.csv
```

### Running the Analysis

#### 1. Data Preprocessing

```r
source("MungeSumstats_MR_data preprocessing.R")

# Harmonizes GWAS data for MR analysis
```

#### 2. Main Analysis

```r
source("Main analysis.R")

# Runs complete UVMR, MVMR, and mediation analyses
```

#### 3. Result Filtering

```r
source("Results_Filter_Helper.R")

# Generate summary report
generate_summary_report()

# Or custom filtering
library(data.table)

uvmr <- fread("results_trial/uvmr_comprehensive_results.csv")

# Filter for robust causal associations
robust <- uvmr[
  grepl("IVW", method) & 
  q_value < 0.05 &                      # FDR significant
  concordance_validated == TRUE &       # Validated by sensitivity
  F_statistic > 10                      # Strong instruments
]
```

---

## 📊 Output Structure

### Result Files

All results are saved in `results_trial/`:

| File | Contents | Columns | Description |
|:-----|:---------|:--------|:------------|
| **uvmr_comprehensive_results.csv** | UVMR results | ~35 | All exposure-outcome pairs with 5 MR methods |
| **mvmr_comprehensive_results.csv** | MVMR results | ~20 | Covariate-adjusted independent effects |
| **mediation_comprehensive_results.csv** | Mediation results | ~32 | All exposure-mediator-outcome pathways |

### Key Output Columns

#### UVMR Results

| Column | Description |
|:-------|:------------|
| `b`, `se`, `pval` | Standard MR effect estimates |
| `q_value` | FDR-corrected p-value |
| `concordance_validated` | Validation by sensitivity methods |
| `F_statistic` | Instrument strength |
| `mrlap_corrected_effect` | MRlap-corrected effect |
| `mrlap_corrected_pval` | MRlap-corrected p-value |
| `mrlap_diff_pval` | Significance of correction |

#### MVMR Results

| Column | Description |
|:-------|:------------|
| `beta_mvmr`, `se_mvmr`, `pval_mvmr` | MVMR effect estimates |
| `q_value` | FDR-corrected p-value |
| `adjusted_for` | **CRITICAL**: Must include SES covariates |
| `n_covariates` | **CRITICAL**: Must be 3 (Education, Income, Occupation) |
| `method` | MVMR method (IVW, Median, Egger, Lasso) |
| `conditional_F_statistic` | MVMR instrument strength |

#### Mediation Results

| Column | Description |
|:-------|:------------|
| `beta_exp_med` | Exposure → Mediator effect |
| `beta_med_out_direct` | Mediator → Outcome direct effect |
| `beta_exp_out_total` | Exposure → Outcome total effect |
| `bidirectional` | **CRITICAL**: Must be `No_Unidirectional` |
| `direction_concordant` | **CRITICAL**: Must be `TRUE` |
| `mediation_proportion` | Mediation proportion from 0 to 1 with 95% CI |

---

## 🎯 Quality Control

### Result Filtering Criteria

Following Ye *et al*. 2024 standards:

#### UVMR — Causal Associations

```r
# IVW estimates considered causal only if:
# - Same direction AND significance as ≥1 sensitivity method

robust_uvmr <- uvmr[
  grepl("IVW", method) &
  q_value < 0.05 &                      # FDR significant
  concordance_validated == TRUE &       # Validated by sensitivity
  F_statistic > 10 &                    # Strong instruments
  (is.na(Q_pval) | Q_pval > 0.05)       # No heterogeneity
]
```

#### MVMR — Independent Effects

```r
# Must be adjusted for SES covariates

independent_mvmr <- mvmr[
  method == "MVMR-IVW" &
  q_value < 0.05 &
  n_covariates == 3 &                   # CRITICAL: must be 3
  grepl("Education", adjusted_for) &    # Must include Education
  grepl("Income", adjusted_for) &       # Must include Income
  grepl("Occupation", adjusted_for)     # Must include Occupation
]
```

#### Mediation — Valid Pathways

**All 4 criteria must be met**:

```r
valid_mediation <- mediation[
  pval_exp_med < 0.05 &                         # Criterion 1: Exp→Med significant
  pval_med_out_direct < 0.05 &                  # Criterion 2: Med→Out significant
  bidirectional == "No_Unidirectional" &        # Criterion 3: No reverse causality
  direction_concordant == TRUE                  # Criterion 4: Direction consistent
]
```

### Quality Control Checklist

- [ ] F-statistics > 10 for all instruments
- [ ] Cochran's Q p-value > 0.05, indicating no strong heterogeneity signal
- [ ] MR-PRESSO outlier test passed
- [ ] Reverse MR shows no bidirectionality
- [ ] Method concordance validated
- [ ] MVMR adjusted for all 3 SES covariates
- [ ] FDR q-value < 0.05 for significant results
- [ ] MRlap correction applied if sample overlap is detected

---

## 💡 Advanced Features

### MRlap Sample Overlap Correction

MRlap provides bias correction for inverse variance weighting Mendelian randomization when there is sample overlap between exposure and outcome GWAS.

#### Original Repository

* **GitHub:** [https://github.com/n-mounier/MRlap](https://github.com/n-mounier/MRlap)
* **Citation:** Mounier N, Kutalik Z. Bias correction for inverse variance weighting Mendelian randomization. *Genetic Epidemiology*. 2023;47(4):314-331. DOI: [10.1002/gepi.22522](https://doi.org/10.1002/gepi.22522)

#### Required Reference Data

> ⚠️ **IMPORTANT:** MRlap requires local LD scores and HapMap3 variant information files. You must download and configure these before running MRlap correction.

| Directory | Contents | Description |
|:----------|:---------|:------------|
| `1000G_Phase3_ldscores/` | LD score files | Pre-computed LD scores from 1000 Genomes Phase 3 |
| `1000G_Phase3_hm3/` | HapMap3 SNP list | HapMap3 variant information for filtering |

**Download Source:** Reference data can be obtained from:

* [MRlap GitHub Repository](https://github.com/n-mounier/MRlap) — See the README for download links
* [LDSC Resources](https://alkesgroup.broadinstitute.org/LDSCORE/) — Original LD score files

**Configuration:** Update the MRlap reference paths in `Main analysis.R`:

```r
# Example configuration
mrlap_ldscores_path <- "path/to/1000G_Phase3_ldscores"
mrlap_hm3_path <- "path/to/1000G_Phase3_hm3"
```

#### Auto-Detection Features

MRlap automatically:

* ✅ Extracts sample sizes from GWAS files
* ✅ Estimates sample overlap using LDSC
* ✅ Corrects for winner's curse and weak instruments
* ✅ No manual overlap configuration needed

### MVMR Combined IV Selection

Following standard practice:

> "Genetic instruments were the combination of SNPs, which had genome-wide significance in **either** the GWAS of each exposure **or** the GWAS of each mediator."

**Implementation:** `select_combined_ivs_mvmr()` function

### Comprehensive Validation

* **FDR Correction:** Controls Type I error in multiple testing.
* **Method Concordance:** Validates IVW with sensitivity analyses.
* **Reverse MR:** Excludes bidirectional relationships.
* **Direction Checking:** Ensures pathway consistency.
* **Heterogeneity Testing:** Cochran's Q for IVW.
* **Outlier Detection:** MR-PRESSO for robust estimates.

---

## 🔧 Troubleshooting

### Common Issues and Solutions

| Issue | Solution |
|:------|:---------|
| **PLINK not found** | Ensure PLINK is installed and added to system PATH, or specify full path in scripts |
| **1000G.EUR.QC files missing** | Download EUR reference panel files `.bed`, `.bim`, `.fam` and update path in scripts |
| **MRlap reference data missing** | Download `1000G_Phase3_ldscores` and `1000G_Phase3_hm3` from MRlap GitHub |
| **MRlap fails to detect sample sizes** | Ensure GWAS files contain `N` or `N_cases`/`N_controls` columns |
| **MVMR fails with SES covariates** | Verify all 3 SES files are in `Covariates_SES/` directory |
| **No instruments found** | Check p-value threshold, default `5e-8`, and LD clumping parameters |
| **MR-PRESSO errors** | Ensure sufficient instruments, minimum 4 SNPs |
| **Memory issues with MRlap** | Increase RAM allocation or use smaller datasets |
| **F-statistics < 10** | Remove weak instruments or use more lenient thresholds |

### Getting Help

1. Check `START_HERE_v2.5.txt` for detailed instructions.
2. Review script comments and documentation.
3. Run `Demo_Test_Analysis.R` to verify installation.
4. Open an issue on [GitHub](https://github.com/Hexiao-DING/Dementia_Depression_MR-analysis/issues).

---

## 📚 Citation and References

### If You Use This Pipeline

Please cite the methodological references:

#### TwoSampleMR

> Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic causal inference across the human phenome. *eLife*. 2018;7:e34408. DOI: [10.7554/eLife.34408](https://doi.org/10.7554/eLife.34408)

#### MVMR

> Sanderson E, Glymour MM, Holmes MV, Kang H, Morrison J, Munafò MR, Palmer T, Schooling CM, Wallace C, Zhao Q, Davey Smith G. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. *International Journal of Epidemiology*. 2019;48(3):713-727. DOI: [10.1093/ije/dyy262](https://doi.org/10.1093/ije/dyy262)

#### MR-PRESSO

> Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases. *Nature Genetics*. 2018;50(5):693-698. DOI: [10.1038/s41588-018-0099-7](https://doi.org/10.1038/s41588-018-0099-7)

#### MRlap

> Mounier N, Kutalik Z. Bias correction for inverse variance weighting Mendelian randomization. *Genetic Epidemiology*. 2023;47(4):314-331. DOI: [10.1002/gepi.22522](https://doi.org/10.1002/gepi.22522)
>
> **GitHub:** [https://github.com/n-mounier/MRlap](https://github.com/n-mounier/MRlap)

#### Reference Standard

> Ye CJ, Liu D, Chen ML, Kong LJ, Dou C, Wang YY, et al. Mendelian randomization evidence for the causal effect of mental well-being on healthy aging. *Nature Human Behaviour*. 2024;8:1798-1809. DOI: [10.1038/s41562-024-01905-9](https://doi.org/10.1038/s41562-024-01905-9)
>
> GitHub: [https://github.com/yechaojie/mental_aging](https://github.com/yechaojie/mental_aging)

### Additional Methodological and Data-Processing References

These references are recommended for users who want to report or extend the pipeline in a manuscript:

> Burgess S, Butterworth A, Thompson SG. Mendelian randomization analysis with multiple genetic variants using summarized data. *Genetic Epidemiology*. 2013;37(7):658-665. DOI: [10.1002/gepi.21758](https://doi.org/10.1002/gepi.21758)

> Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. *International Journal of Epidemiology*. 2015;44(2):512-525. DOI: [10.1093/ije/dyv080](https://doi.org/10.1093/ije/dyv080)

> Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator. *Genetic Epidemiology*. 2016;40(4):304-314. DOI: [10.1002/gepi.21965](https://doi.org/10.1002/gepi.21965)

> Hartwig FP, Davey Smith G, Bowden J. Robust inference in summary data Mendelian randomization via the zero modal pleiotropy assumption. *International Journal of Epidemiology*. 2017;46(6):1985-1998. DOI: [10.1093/ije/dyx102](https://doi.org/10.1093/ije/dyx102)

> Murphy AE, Schilder BM, Skene NG. MungeSumstats: a Bioconductor package for the standardization and quality control of many GWAS summary statistics. *Bioinformatics*. 2021;37(23):4593-4596. DOI: [10.1093/bioinformatics/btab665](https://doi.org/10.1093/bioinformatics/btab665)

> The 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature*. 2015;526:68-74. DOI: [10.1038/nature15393](https://doi.org/10.1038/nature15393)

> Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*. 1995;57(1):289-300. DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

---

## 🙏 Acknowledgements

| Role | Name | Affiliation |
|:-----|:-----|:------------|
| **Contributors** | Hexiao Ding | The Hong Kong Polytechnic University |
| **Contributors** | Jing Lan | The Hong Kong Polytechnic University |
| **Contributors** | Na Li | Sichuan University |
| **Supervisor** | Dr. Jung Sun Yoo | The Hong Kong Polytechnic University |

**Institutions:**

* The Hong Kong Polytechnic University
* Sichuan University

**Version:** 2.7  
**Last Updated:** 2025-12-10

### Methodological Acknowledgments

This pipeline implements methods from:

* **TwoSampleMR:** MR-Base platform
* **MVMR:** Multivariable MR package
* **MR-PRESSO:** Outlier detection and correction
* **MRlap:** Sample overlap correction ([GitHub](https://github.com/n-mounier/MRlap))

---

<div align="center">

### 🧬 Built for reproducible causal inference in neuropsychiatric comorbidity research

**Dementia · Depression · Molecular Epidemiology · Mendelian Randomization · Mediation · SES Adjustment**

</div>
