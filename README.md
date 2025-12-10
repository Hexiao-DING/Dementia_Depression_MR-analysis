<div align="center">

# 🧬 Dementia–Depression Mendelian Randomization Analysis Pipeline

**A comprehensive, methodologically rigorous pipeline for causal inference in dementia–depression comorbidity**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D%204.0-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux-lightgrey)]()
[![Status](https://img.shields.io/badge/Status-Production%20Ready-success)]()

[Documentation](#-overview) • [Quick Start](#-quick-start) • [Methods](#-methodological-compliance) • [Citation](#-citation-and-references)

</div>

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Features](#-features)
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

### Key Objectives

- 🔬 **Causal Inference**: Establish causal relationships between exposures and outcomes using genetic instruments
- 📊 **Multivariable Analysis**: Adjust for confounders and mediators using MVMR
- 🔄 **Mediation Analysis**: Identify intermediate pathways in causal chains
- ✅ **Rigorous Validation**: Comprehensive sensitivity analyses and quality control
- 📈 **Sample Overlap Correction**: Automatic correction using MRlap
- 🎯 **SES Adjustment**: Mandatory socioeconomic status covariate adjustment

### Analysis Framework

```
┌───────────────────────────────────────────────────────────┐
│                    MR Analysis Pipeline                   │
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

---

## ✨ Features

- 🔬 **Comprehensive MR Methods**: UVMR (5 methods) and MVMR (4 methods)
- 🛡️ **Rigorous QC**: MR-PRESSO outlier detection, heterogeneity testing, reverse MR
- 📊 **Sample Overlap Correction**: Automatic MRlap implementation
- 🎯 **SES Adjustment**: Mandatory socioeconomic status covariate adjustment in MVMR
- 🔄 **Mediation Analysis**: Complete pathway decomposition
- ✅ **Method Concordance**: Validation across multiple MR methods
- 📈 **FDR Correction**: Multiple testing correction with Benjamini-Hochberg
- 🔍 **Combined IV Selection**: Standard-compliant instrument selection for MVMR
- 📚 **Well-Documented**: Comprehensive documentation and examples

- ## 📋 Data Requirements [All data are EU populations]

### 🔹 Required genome-wide association studies (GWAS) Summary Statistics, protein quantitative trait loci (pQTLs), and metabolite quantitative trait loci (mQTLs)
> ⚠️ You need to provide these files (they are **NOT included** in the repository).

---

### 1️⃣ Exposure pQTLs/GWAS — Multiple Exposure Files

#### (1) **Standardized Circulating Human Plasma Proteome Data (pQTLs); EUR Population N=3301**
- `GSCT005806_GRCh37.tsv.gz`
- `GCST90240120_GRCh37.tsv.gz` → `GCST90243401_GRCh37.tsv.gz`
- Data pre-processing strategies: Standardize pQTLs into homogeneous format using `MungeSumstats version 1.14.1`.
- Source: `https://www.ebi.ac.uk/gwas/publications/29875488`
- **Citation**: Sun, B. B., Maranville, J. C., Peters, J. E., Stacey, D., Staley, J. R., Blackshaw, J., ... & Butterworth, A. S. (2018). Genomic atlas of the human plasma proteome. Nature, 558(7708), 73-79.
#### (2) **Standardized Circulating Metabolic Biomarkers Data (GWAS); EUR Population N= **
- `GCST90301941.tsv` → `GCST90302173.tsv`
- Data pre-processing strategies: Standardize GWAS summary statistics into homogeneous format using `MungeSumstats version 1.14.1`.
- Source: `https://www.ebi.ac.uk/gwas/publications/38448586`
- **Citation**: Karjalainen, M. K., Karthikeyan, S., Oliver-Williams, C., Sliz, E., Allara, E., Fung, W. T., ... & Kettunen, J. (2024). Genome-wide characterization of circulating metabolic biomarkers. Nature, 628(8006), 130-138.
#### (3) **Circulating Inflammatory Proteins Data (pQTLs); EUR Population N= **
- `GCST90274758.tsv.gz` → `GCST90274848.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/37563310`
- **Citation**: Zhao, J. H., Stacey, D., Eriksson, N., Macdonald-Dunlop, E., Hedman, Å. K., Kalnapenkis, A., ... & Peters, J. E. (2023). Genetics of circulating inflammatory proteins identifies drivers of immune-mediated disease risk and therapeutic targets. Nature immunology, 24(9), 1540-1551.

---

### 2️⃣ Mediator mQTLs — Mediator Files
## ✨ Features

#### (1) **Cerebrospinal Fluid Metabolomics Data (metabQTLs); EUR Population N= **
- `GCST90025999_buildGRCh37.tsv.gz` → `GCST90026336_buildGRCh37.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/33437055`
- **Citation**: Panyard, D. J., Kim, K. M., Darst, B. F., Deming, Y. K., Zhong, X., Wu, Y., ... & Lu, Q. (2021). Cerebrospinal fluid metabolomics identifies 19 brain-related phenotype associations. Communications biology, 4(1), 63.

---

### 3️⃣ Outcome GWAS — Outcome Files

#### (1) Dementia Types *(7 traits in total, GWAS)*
##### Dementia → `Finn-b-F5_DEMENTIA.tsv.gz`; EUR Population N= 
- Source: `https://r10.risteys.finngen.fi/endpoints/F5_DEMENTIA`
- **Citation**: Kurki, M. I., Karjalainen, J., Palta, P., Sipilä, T. P., Kristiansson, K., Donner, K. M., ... & Waring, J. (2023). FinnGen provides genetic insights from a well-phenotyped isolated population. Nature, 613(7944), 508-518.
##### Alzheimer's Disease → `GCST90012877.tsv.gz`; EUR Population N= 
- Source: `https://www.ebi.ac.uk/gwas/publications/33589840`
- **Citation**: Schwartzentruber, J., Cooper, S., Liu, J. Z., Barrio-Hernandez, I., Bello, E., Kumasaka, N., ... & Bassett, A. (2021). Genome-wide meta-analysis, fine-mapping and integrative prioritization implicate new Alzheimer’s disease risk genes. Nature genetics, 53(3), 392-402.
##### Cognitive Performance → `GCST006572.tsv.gz`; EUR Population N= 
- Source: `https://www.ebi.ac.uk/gwas/publications/30038396`
- **Citation**: Lee, J. J., Wedow, R., Okbay, A., Kong, E., Maghzian, O., Zacher, M., ... & Cesarini, D. (2018). Gene discovery and polygenic prediction from a genome-wide association study of educational attainment in 1.1 million individuals. Nature genetics, 50(8), 1112-1121.
##### Vascular Dementia → `Finn-b-F5_VASCDEM.tsv.gz`; EUR Population N= 
- Source: `https://r9.risteys.finngen.fi/endpoints/F5_VASCDEM`
- **Citation**: Kurki, M. I., Karjalainen, J., Palta, P., Sipilä, T. P., Kristiansson, K., Donner, K. M., ... & Waring, J. (2023). FinnGen provides genetic insights from a well-phenotyped isolated population. Nature, 613(7944), 508-518.
##### Lewy Body Dementia → `GCST90001390.tsv.gz`; EUR Population N= 
- Source: `https://www.ebi.ac.uk/gwas/publications/33589841`
- **Citation**: Chia, R., Sabir, M. S., Bandres-Ciga, S., Saez-Atienzar, S., Reynolds, R. H., Gustavsson, E., ... & Dickson, D. W. (2021). Genome sequencing analysis identifies new loci associated with Lewy body dementia and provides insights into its genetic architecture. Nature genetics, 53(3), 294-303.
##### Frontotemporal Dementia → `GCST90558311.tsv.gz`; EUR Population N= 
- Source: `https://www.ebi.ac.uk/gwas/publications/40280976`
- **Citation**: Pottier, C., Küçükali, F., Baker, M., Batzler, A., Jenkins, G. D., van Blitterswijk, M., ... & Rademakers, R. (2025). Deciphering distinct genetic risk factors for FTLD-TDP pathological subtypes via whole-genome sequencing. Nature communications, 16(1), 3914.
##### Undefined Dementia → `Finn_b_F5_DEMNAS.tsv.gz`; EUR Population N= 
- Source: `https://r12.risteys.finngen.fi/endpoints/F5_DEMNAS`
- **Citation**: Kurki, M. I., Karjalainen, J., Palta, P., Sipilä, T. P., Kristiansson, K., Donner, K. M., ... & Waring, J. (2023). FinnGen provides genetic insights from a well-phenotyped isolated population. Nature, 613(7944), 508-518.

#### (2) Depression Types *(3 traits in total, GWAS)*
##### Depressive Disorders → `GCST90476922.tsv.gz`; EUR Population N= 
- Source: `https://www.ebi.ac.uk/gwas/publications/39024449`
- **Citation**: Verma, A., Huffman, J. E., Rodriguez, A., Conery, M., Liu, M., Ho, Y. L., ... & Liao, K. P. (2024). Diversity and scale: Genetic architecture of 2068 traits in the VA Million Veteran Program. Science, 385(6706), eadj1182.

##### Major Depressive Disorders → `GCST90468123.tsv.gz`; EUR Population N= 
- Source: `https://www.ebi.ac.uk/gwas/publications/39789286`
- **Citation**: Loya, H., Kalantzis, G., Cooper, F., & Palamara, P. F. (2025). A scalable variational inference approach for increased mixed-model association power. Nature Genetics, 57(2), 461-468.

##### Mixed Anxiety and Depressive Disorder → `GCST90225526.tsv.gz`; EUR Population N= 
- Source: `https://www.ebi.ac.uk/gwas/publications/37259642`
- **Citation**: Brasher, M. S., Mize, T. J., Thomas, A. L., Hoeffer, C. A., Ehringer, M. A., & Evans, L. M. (2023). Testing associations between human anxiety and genes previously implicated by mouse anxiety models. Genes, Brain and Behavior, 22(6), e12851.

---

### 4️⃣ Covariate GWAS — Social Economic Status *(SES, 3 traits in total, GWAS)*

#### (1) Education → `Education_GCST003676.txt.gz`; EUR Population N= 
- Source: `https://www.ebi.ac.uk/gwas/publications/27225129`
- **Citation**: Okbay, A., Beauchamp, J. P., Fontana, M. A., Lee, J. J., Pers, T. H., Rietveld, C. A., ... & Rustichini, A. (2016). Genome-wide association study identifies 74 loci associated with educational attainment. Nature, 533(7604), 539-542.

#### (2) Income → `Income_GCST90566700.tsv.gz`; EUR Population N=   
- Source: `https://www.ebi.ac.uk/gwas/publications/40360725`
- **Citation**: Xia, C., Lu, Y., Zhou, Z., Marchi, M., Kweon, H., Ning, Y., ... & Hill, W. D. (2025). Deciphering the influence of socioeconomic status on brain structure: insights from Mendelian randomization. Molecular Psychiatry, 1-14.

#### (3) Occupation → `Occupation_GCST90566702.tsv.gz`; EUR Population N= 
- Source:  `https://www.ebi.ac.uk/gwas/publications/40360725`
- **Citation**: Xia, C., Lu, Y., Zhou, Z., Marchi, M., Kweon, H., Ning, Y., ... & Hill, W. D. (2025). Deciphering the influence of socioeconomic status on brain structure: insights from Mendelian randomization. Molecular Psychiatry, 1-14.

---

### File Format

**Required columns** (flexible naming):
- SNP identifier: `SNP`, `rsid`, `rs_id`, `MarkerName`, `variant`
- Effect allele: `effect_allele`, `EA`, `A1`
- Other allele: `other_allele`, `OA`, `A2`
- Beta: `beta`, `BETA`, `b`, `Beta`
- Standard error: `se`, `SE`, `standard_error`
- P-value: `pval`, `P`, `p_value`, `Pval`

**Optional but recommended**:
- Chromosome: `chr`, `CHR`, `chromosome`
- Position: `pos`, `BP`, `base_pair_location`, `POS`
- Sample size: `n`, `N`, `samplesize`
- Effect allele frequency: `eaf`, `EAF`, `freq`

**Supported formats**: `.tsv`, `.tsv.gz`, `.txt`, `.txt.gz`

---

## 🎯 IV Selection Parameters

Following genome-wide standards:

- **P-value threshold**: 5×10⁻⁸ (genome-wide significance)
- **LD clumping**: r² < 0.001, window = 10,000 kb
- **Minimum SNPs**: 3 per analysis

---

## 📊 Output Files

### Main Results (3 CSV files)

1. **uvmr_comprehensive_results.csv** (~35 columns)
   - All exposure-outcome pairs
   - 5 MR methods per pair
   - MRlap sample overlap correction
   - FDR q-values
   - Method concordance validation

2. **mvmr_comprehensive_results.csv** (~20 columns)
   - Covariate-adjusted independent effects
   - 4 MVMR methods
   - SES adjustment information
   - FDR q-values

3. **mediation_comprehensive_results.csv** (~32 columns)
   - All exposure-mediator-outcome paths
   - Reverse MR results
   - Direction concordance
   - Mediation proportions with 95% CI

---

## 🔬 Methodological Compliance

### Reference Standard

This pipeline implements methods following the rigorous standards established in:

> **Ye CJ, Liu D, Chen ML, *et al*. Mendelian randomization evidence for the causal effect of mental well-being on healthy aging.**  
> *Nature Human Behaviour*. 2024.  
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

**Total Compliance**: ✅ **15/15 Core Requirements** (100%)

---

## 💻 Installation

### Prerequisites

- **R** ≥ 4.0.0
- **RAM**: 8+ GB (16+ GB recommended for MRlap)
- **Internet Connection**: Required for package installation and MRlap reference data
- **Disk Space**: Sufficient space for GWAS files and results

### Automated Installation (Recommended)

```r
# Clone the repository
git clone https://github.com/Hexiao-DING/Dementia_Depression_MR-analysis.git
cd Dementia_Depression_MR-analysis

# Run automated installation script
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
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("n-mounier/MRlap")
```

### Verify Installation

```r
# Test package loading
source("Install_All_Packages.R")

# Test covariate file loading
source("Test_Covariate_Loading.R")

# Run demo analysis (recommended before full analysis)
source("Demo_Test_Analysis.R")
```

---

## 🚀 Quick Start

### Step 1: Prepare Data

1. **Place GWAS files** in designated directories:
   - Exposures: `Standardized Circulating human plasma proteome_Data/`
   - Mediators: `Cerebrospinal fluid metabolomics_Data/`
   - Outcomes: `Outcomes/`

2. **Download SES Covariates** (REQUIRED for MVMR):
   - Education
   - Income
   - Occupation
   - Place in: `Covariates_SES/`

3. **Configure paths** in `Main analysis.R` (lines 61-72)

### Step 2: Run Analysis

```r
# Full analysis (may take hours depending on data size)
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

```
.
├── Main analysis.R                          # Main analysis script
├── Install_All_Packages.R                   # Package installation
├── MungeSumstats_MR_data preprocessing.R    # Data preprocessing
├── Results_Filter_Helper.R                  # Result filtering utilities
├── Demo_Test_Analysis.R                     # Demo/test analysis
├── Test_Covariate_Loading.R                 # Covariate testing
├── START_HERE_v2.5.txt                      # Getting started guide
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
  q_value < 0.05 &                    # FDR significant
  concordance_validated == TRUE &     # Validated by sensitivity
  F_statistic > 10                    # Strong instruments
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
| `beta_med_out_direct` | Mediator → Outcome (direct) effect |
| `beta_exp_out_total` | Exposure → Outcome (total) effect |
| `bidirectional` | **CRITICAL**: Must be "No_Unidirectional" |
| `direction_concordant` | **CRITICAL**: Must be TRUE |
| `mediation_proportion` | Mediation proportion (0-1) with 95% CI |

---

## 🎯 Quality Control

### Result Filtering Criteria

Following Ye *et al*. 2024 standards:

#### UVMR - Causal Associations

```r
# IVW estimates considered causal only if:
# - Same direction AND significance as ≥1 sensitivity method
robust_uvmr <- uvmr[
  grepl("IVW", method) &
  q_value < 0.05 &                    # FDR significant
  concordance_validated == TRUE &     # Validated by sensitivity
  F_statistic > 10 &                  # Strong instruments
  (is.na(Q_pval) | Q_pval > 0.05)    # No heterogeneity
]
```

#### MVMR - Independent Effects

```r
# Must be adjusted for SES covariates
independent_mvmr <- mvmr[
  method == "MVMR-IVW" &
  q_value < 0.05 &
  n_covariates == 3 &                 # CRITICAL: must be 3
  grepl("Education", adjusted_for) &  # Must include Education
  grepl("Income", adjusted_for) &     # Must include Income
  grepl("Occupation", adjusted_for)  # Must include Occupation
]
```

#### Mediation - Valid Pathways

**All 4 criteria must be met**:

```r
valid_mediation <- mediation[
  pval_exp_med < 0.05 &                        # Criterion 1: Exp→Med significant
  pval_med_out_direct < 0.05 &                 # Criterion 2: Med→Out significant
  bidirectional == "No_Unidirectional" &       # Criterion 3: No reverse causality
  direction_concordant == TRUE                 # Criterion 4: Direction consistent
]
```

### Quality Control Checklist

- [ ] F-statistics > 10 for all instruments
- [ ] Cochran's Q p-value > 0.05 (no heterogeneity)
- [ ] MR-PRESSO outlier test passed
- [ ] Reverse MR shows no bidirectionality
- [ ] Method concordance validated
- [ ] MVMR adjusted for all 3 SES covariates
- [ ] FDR q-value < 0.05 for significant results
- [ ] MRlap correction applied (if sample overlap detected)

---

## 💡 Advanced Features

### MRlap Auto-Detection

MRlap automatically:

- ✅ Extracts sample sizes from GWAS files
- ✅ Estimates sample overlap using LDSC
- ✅ Corrects for winner's curse and weak instruments
- ✅ No manual configuration needed!

### MVMR Combined IV Selection

Following standard practice:

> "Genetic instruments were the combination of SNPs, which had genome-wide significance in **either** the GWAS of each exposure **or** the GWAS of each mediator"

**Implementation**: `select_combined_ivs_mvmr()` function

### Comprehensive Validation

- **FDR Correction**: Controls Type I error in multiple testing
- **Method Concordance**: Validates IVW with sensitivity analyses
- **Reverse MR**: Excludes bidirectional relationships
- **Direction Checking**: Ensures pathway consistency
- **Heterogeneity Testing**: Cochran's Q for IVW
- **Outlier Detection**: MR-PRESSO for robust estimates

---

## 🔧 Troubleshooting

### Common Issues and Solutions

| Issue | Solution |
|:------|:---------|
| **MRlap fails to detect sample sizes** | Ensure GWAS files contain `N` or `N_cases`/`N_controls` columns |
| **MVMR fails with SES covariates** | Verify all 3 SES files are in `Covariates_SES/` directory |
| **No instruments found** | Check p-value threshold (default: 5e-8) and LD clumping parameters |
| **MR-PRESSO errors** | Ensure sufficient instruments (minimum 4 SNPs) |
| **Memory issues with MRlap** | Increase RAM allocation or use smaller datasets |
| **F-statistics < 10** | Remove weak instruments or use more lenient thresholds |

### Getting Help

1. Check `START_HERE_v2.5.txt` for detailed instructions
2. Review script comments and documentation
3. Run `Demo_Test_Analysis.R` to verify installation
4. Open an issue on [GitHub](https://github.com/Hexiao-DING/Dementia_Depression_MR-analysis/issues)

---

## 📚 Citation and References

### If You Use This Pipeline

Please cite the methodological references:

#### TwoSampleMR

```
Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports 
systematic causal inference across the human phenome. eLife. 2018;7:e34408.
DOI: 10.7554/eLife.34408
```

#### MVMR

```
Sanderson E, Glymour MM, Holmes MV, et al. An examination of multivariable 
Mendelian randomization in the single-sample and two-sample summary data settings. 
Int J Epidemiol. 2019;48(3):713-727.
DOI: 10.1093/ije/dyy262
```

#### MR-PRESSO

```
Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal 
pleiotropy in causal relationships inferred from Mendelian randomization between 
complex traits and diseases. Nat Genet. 2018;50(5):693-698.
DOI: 10.1038/s41588-018-0099-7
```

#### MRlap

```
Mounier N, Kutalik Z. Bias correction for inverse variance weighting 
Mendelian randomization. Genet Epidemiol. 2023;47(4):314-331.
DOI: 10.1002/gepi.22522
```

#### Reference Standard

```
Ye CJ, Liu D, Chen ML, et al. Mendelian randomization evidence for the causal 
effect of mental well-being on healthy aging. Nat Hum Behav. 2024.
DOI: [To be updated]
GitHub: https://github.com/yechaojie/mental_aging
```

---

## 🙏 Acknowledgements

| Role | Name | Affiliation |
|:-----|:-----|:------------|
| **Contributors** | Hexiao Ding | The Hong Kong Polytechnic University |
| **Contributors** | Jing Lan | The Hong Kong Polytechnic University |
| **Contributors** | Na Li | Sichuan University |
| **Supervisor** | Dr. Jung Sun Yoo | The Hong Kong Polytechnic University |

**Institutions**:
- The Hong Kong Polytechnic University
- Sichuan University

**Version**: 2.6  
**Last Updated**: 2025-11-23

### Methodological Acknowledgments

This pipeline implements methods from:
- **TwoSampleMR**: MR-Base platform
- **MVMR**: Multivariable MR package
- **MR-PRESSO**: Outlier detection and correction
- **MRlap**: Sample overlap correction



