# Dementia & Depression Traits and Molecular Epidemiology 
# R code for MR-analysis
This Mendelian Randomization pipeline provides a comprehensive and methodologically rigorous framework for causal inference in complex traits. It implements a complete analytical workflow integrating univariable MR (UVMR), multivariable MR (MVMR) with socioeconomic status covariate adjustment, and formal mediation analysis. The pipeline incorporates robust sensitivity analyses for pleiotropy, heterogeneity, and sample overlap using MR-Egger, MR-PRESSO, and MRlap, while ensuring methodological validity through FDR correction, reverse MR testing, and concordance validation. Designed for high-dimensional molecular epidemiology studies, it enables reliable discovery of causal protein-disease and metabolite-disease relationships with comprehensive confounder control.

[![R Version](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)](https://www.r-project.org/)
[![Compliance](https://img.shields.io/badge/Compliance-100%25-brightgreen)](https://github.com/yechaojie/mental_aging)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)]()

---

## 📋 Overview

This pipeline provides a complete implementation of three complementary Mendelian Randomization (MR) analyses:

1. **UVMR** (Univariable MR) - Individual exposure-outcome associations
2. **MVMR** (Multivariable MR) - Independent effects with confounder adjustment
3. **Mediation Analysis** - Two-step MR for pathway analysis

### Key Features

✅ **All 23 methodological requirements** implemented  
✅ **4 required R packages**: TwoSampleMR, MVMR, MRPRESSO, MRlap  
✅ **Complete sensitivity analyses**: Heterogeneity, pleiotropy, outliers, sample overlap  
✅ **FDR correction** (Benjamini-Hochberg)  
✅ **SES covariate adjustment** (Education, Income, Occupation)  
✅ **Reverse MR** for bidirectionality testing  
✅ **Comprehensive English documentation**  
✅ **Server deployment ready** with demo testing

---

## 🚀 Quick Start

### Local Testing

```r
# 1. Install all required packages
source("Install_All_Packages.R")

# 2. Quick test
source("Demo_Test_Analysis.R")

# 3. Run full analysis
source("Main analysis.R")
```

---

## 📊 Analyses Implemented

### 1. UVMR (Univariable MR)

**Methods** (5):
- Inverse Variance Weighted (IVW) - main method
- MR-Egger - pleiotropy detection
- Weighted Median - robust to outliers
- Weighted Mode
- Simple Mode

**Sensitivity Analyses**:
- Heterogeneity test (Cochran's Q)
- Pleiotropy test (MR-Egger intercept)
- Outlier detection (MR-PRESSO)
- Sample overlap correction (MRlap)
- Instrument strength (F-statistic)

**Validation**:
- FDR correction (Benjamini-Hochberg)
- Method concordance checking

### 2. MVMR (Multivariable MR)

**Methods** (4):
- MVMR-IVW - main method
- MVMR-Median - sensitivity
- MVMR-Egger - pleiotropy detection
- MVMR-Lasso - variable selection

**Key Features**:
- **SES covariate adjustment** (Education, Income, Occupation)
- Combined IV selection (exposure + covariate SNPs)
- Conditional F-statistics
- Returns covariate-adjusted independent effects

### 3. Mediation Analysis

**Two-Step Framework**:
- Step 1: Exposure → Mediator (UVMR)
- Reverse MR: Mediator → Exposure (bidirectionality test)
- Step 2: [Mediator + Exposure] → Outcome (MVMR)
- Total: Exposure → Outcome

**Outputs**:
- Mediation proportion with 95% CI (Delta method)
- Direction concordance checking
- Validation of all 4 mediation criteria

---

## 📦 Required R Packages

| Package | Version | Purpose | Installation |
|---------|---------|---------|--------------|
| TwoSampleMR | ≥0.5.7 | UVMR framework | `install.packages("TwoSampleMR")` |
| MVMR | ≥0.4 | Multivariable MR | `install.packages("MVMR")` or GitHub |
| MRPRESSO | ≥1.0 | Outlier detection | `install.packages("MRPRESSO")` or GitHub |
| MRlap | ≥0.0.3 | Sample overlap | `devtools::install_github("n-mounier/MRlap")` |
| ieugwasr | Latest | LD clumping | `install.packages("ieugwasr")` |
| data.table | Latest | Data manipulation | `install.packages("data.table")` |

**One-Command Installation**:
```r
source("Install_All_Packages.R")
```

---

## 🗂️ Project Structure

```
MR_pipeline_demo/
│
├── 📜 Core Analysis Scripts
│   ├── Main analysis.R              ← Main analysis (UVMR + MVMR + Mediation)
│   ├── Install_All_Packages.R         ← Automated package installation
│   ├── Demo_Test_Analysis.R           ← Quick test (5-15 min)
│   ├── Results_Filter_Helper.R        ← Result filtering & summarization
│   └── Test_Covariate_Loading.R       ← Covariate validation
│
├── ⚙️ Configuration
│   └── sample_overlap_config.csv      ← MRlap config (auto-detect by default, empty file)
│
├── 📚 Documentation
│   ├── README.md                      ← This file
│
├── 📁 Data Directories (NOT in GitHub - too large)
│   ├── Covariates_SES                 ← Social-economic statu indicators (Education, Income, Occupation; with descriptions)
│   ├── Exposure Data Directories      ← GWAS data, human proteins (A list and descriptions)
│   ├── Mediator Data Directories      ← GWAS data, CSFs (A list and descriptions)
│   ├── Outcome Data Directories       ← GWAS data, Dementia types and Depression types (A list and descriptions)
│   └── local_clump                    ← LD reference panel
│
└── 📊 Output (NOT in GitHub - generated by analysis)
    └── results_trial/
        ├── uvmr_comprehensive_results.csv
        ├── mvmr_comprehensive_results.csv
        ├── mediation_comprehensive_results.csv
        └── demo_test/                 ← Demo output
```

---

## 📋 Data Requirements [All data are EU populations]

### 🔹 Required GWAS Summary Statistics
> ⚠️ You need to provide these files (they are **NOT included** in the repository).

---

### 1️⃣ Exposure GWAS — Multiple Exposure Files

#### (1) **Standardized Circulating Human Plasma Proteome Data**
- `GSCT005806_GRCh37.tsv.gz`
- `GCST90240120_GRCh37.tsv.gz` → `GCST90243401_GRCh37.tsv.gz`
- Data pre-processing strategies: 

- Source: `https://www.ebi.ac.uk/gwas/publications/29875488`
- **Citation**: Sun, B. B., Maranville, J. C., Peters, J. E., Stacey, D., Staley, J. R., Blackshaw, J., ... & Butterworth, A. S. (2018). Genomic atlas of the human plasma proteome. Nature, 558(7708), 73-79.
#### (2) **Standardized Circulating Metabolic Biomarkers Data**
- `GCST90301941.tsv` → `GCST90302173.tsv`
- Data pre-processing strategies: 

- Source: `https://www.ebi.ac.uk/gwas/publications/29875488](https://www.ebi.ac.uk/gwas/publications/38448586`
- **Citation**: Karjalainen, M. K., Karthikeyan, S., Oliver-Williams, C., Sliz, E., Allara, E., Fung, W. T., ... & Kettunen, J. (2024). Genome-wide characterization of circulating metabolic biomarkers. Nature, 628(8006), 130-138.
#### (3) **Circulating Inflammatory Proteins Data**
- `GCST90274758.tsv.gz` → `GCST90274848.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/29875488](https://www.ebi.ac.uk/gwas/publications/37563310`
- **Citation**: Zhao, J. H., Stacey, D., Eriksson, N., Macdonald-Dunlop, E., Hedman, Å. K., Kalnapenkis, A., ... & Peters, J. E. (2023). Genetics of circulating inflammatory proteins identifies drivers of immune-mediated disease risk and therapeutic targets. Nature immunology, 24(9), 1540-1551.

---

### 2️⃣ Mediator GWAS — Mediator Files

#### (1) **Cerebrospinal Fluid Metabolomics Data**
- `GCST90025999_buildGRCh37.tsv.gz` → `GCST90026336_buildGRCh37.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/29875488](https://www.ebi.ac.uk/gwas/publications/37563310](https://www.ebi.ac.uk/gwas/publications/33437055`
- **Citation**: Panyard, D. J., Kim, K. M., Darst, B. F., Deming, Y. K., Zhong, X., Wu, Y., ... & Lu, Q. (2021). Cerebrospinal fluid metabolomics identifies 19 brain-related phenotype associations. Communications biology, 4(1), 63.

---

### 3️⃣ Outcome GWAS — Outcome Files

#### (1) Dementia Types *(7 traits in total)*
##### Dementia → `Finn-b-F5_DEMENTIA.tsv.gz`
- Source: `https://r10.risteys.finngen.fi/endpoints/F5_DEMENTIA`
- **Citation**: Kurki, M. I., Karjalainen, J., Palta, P., Sipilä, T. P., Kristiansson, K., Donner, K. M., ... & Waring, J. (2023). FinnGen provides genetic insights from a well-phenotyped isolated population. Nature, 613(7944), 508-518.
##### Alzheimer's Disease → `GCST90012877.tsv.gz`
- Source: `https://r10.risteys.finngen.fi/endpoints/F5_DEMENTIA](https://www.ebi.ac.uk/gwas/publications/33589840`
- **Citation**: Schwartzentruber, J., Cooper, S., Liu, J. Z., Barrio-Hernandez, I., Bello, E., Kumasaka, N., ... & Bassett, A. (2021). Genome-wide meta-analysis, fine-mapping and integrative prioritization implicate new Alzheimer’s disease risk genes. Nature genetics, 53(3), 392-402.
##### Cognitive Performance → `GCST006572.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/30038396`
- **Citation**: Lee, J. J., Wedow, R., Okbay, A., Kong, E., Maghzian, O., Zacher, M., ... & Cesarini, D. (2018). Gene discovery and polygenic prediction from a genome-wide association study of educational attainment in 1.1 million individuals. Nature genetics, 50(8), 1112-1121.
##### Vascular Dementia → `Finn-b-F5_VASCDEM.tsv.gz`
- Source: `https://r9.risteys.finngen.fi/endpoints/F5_VASCDEM`
- **Citation**: Kurki, M. I., Karjalainen, J., Palta, P., Sipilä, T. P., Kristiansson, K., Donner, K. M., ... & Waring, J. (2023). FinnGen provides genetic insights from a well-phenotyped isolated population. Nature, 613(7944), 508-518.
##### Lewy Body Dementia → `GCST90001390.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/33589841`
- **Citation**: Chia, R., Sabir, M. S., Bandres-Ciga, S., Saez-Atienzar, S., Reynolds, R. H., Gustavsson, E., ... & Dickson, D. W. (2021). Genome sequencing analysis identifies new loci associated with Lewy body dementia and provides insights into its genetic architecture. Nature genetics, 53(3), 294-303.
##### Frontotemporal Dementia → `GCST90558311.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/40280976`
- **Citation**: Pottier, C., Küçükali, F., Baker, M., Batzler, A., Jenkins, G. D., van Blitterswijk, M., ... & Rademakers, R. (2025). Deciphering distinct genetic risk factors for FTLD-TDP pathological subtypes via whole-genome sequencing. Nature communications, 16(1), 3914.
##### Undefined Dementia → `Finn_b_F5_DEMNAS.tsv.gz`
- Source: `https://r12.risteys.finngen.fi/endpoints/F5_DEMNAS`
- **Citation**: Kurki, M. I., Karjalainen, J., Palta, P., Sipilä, T. P., Kristiansson, K., Donner, K. M., ... & Waring, J. (2023). FinnGen provides genetic insights from a well-phenotyped isolated population. Nature, 613(7944), 508-518.

#### (2) Depression Types *(3 traits in total)*
##### Depressive Disorders → `GCST90476922.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/39024449`
- **Citation**: Verma, A., Huffman, J. E., Rodriguez, A., Conery, M., Liu, M., Ho, Y. L., ... & Liao, K. P. (2024). Diversity and scale: Genetic architecture of 2068 traits in the VA Million Veteran Program. Science, 385(6706), eadj1182.

##### Major Depressive Disorders → `GCST90468123.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/39789286`
- **Citation**: Loya, H., Kalantzis, G., Cooper, F., & Palamara, P. F. (2025). A scalable variational inference approach for increased mixed-model association power. Nature Genetics, 57(2), 461-468.

##### Mixed Anxiety and Depressive Disorder → `GCST90225526.tsv.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/37259642`
- **Citation**: Brasher, M. S., Mize, T. J., Thomas, A. L., Hoeffer, C. A., Ehringer, M. A., & Evans, L. M. (2023). Testing associations between human anxiety and genes previously implicated by mouse anxiety models. Genes, Brain and Behavior, 22(6), e12851.

---

### 4️⃣ Covariate GWAS — Social Economic Status (SES)

#### (1) Education → `Education_GCST003676.txt.gz`
- Source: `https://www.ebi.ac.uk/gwas/publications/27225129`
- **Citation**: Okbay, A., Beauchamp, J. P., Fontana, M. A., Lee, J. J., Pers, T. H., Rietveld, C. A., ... & Rustichini, A. (2016). Genome-wide association study identifies 74 loci associated with educational attainment. Nature, 533(7604), 539-542.

#### (2) Income → `Income_GCST90566700.tsv.gz`  
- Source: `https://www.ebi.ac.uk/gwas/publications/40360725`
- **Citation**: Xia, C., Lu, Y., Zhou, Z., Marchi, M., Kweon, H., Ning, Y., ... & Hill, W. D. (2025). Deciphering the influence of socioeconomic status on brain structure: insights from Mendelian randomization. Molecular Psychiatry, 1-14.

#### (3) Occupation → `Occupation_GCST90566702.tsv.gz`
- Source:  `https://www.ebi.ac.uk/gwas/publications/40360725`
- **Citation**: Xia, C., Lu, Y., Zhou, Z., Marchi, M., Kweon, H., Ning, Y., ... & Hill, W. D. (2025). Deciphering the influence of socioeconomic status on brain structure: insights from Mendelian randomization. Molecular Psychiatry, 1-14.

---

### GWAS File Format

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

Based on:
> Ye CJ, Liu D, Chen ML, et al. **Mendelian randomization evidence for the causal effect of mental well-being on healthy aging**. *Nature Human Behaviour*. 2024.  
> GitHub: https://github.com/yechaojie/mental_aging

### Implementation Checklist

| Component | Status |
|-----------|--------|
| UVMR with 5 methods | ✅ |
| MVMR with 4 methods | ✅ |
| MR-PRESSO outlier detection | ✅ |
| MRlap sample overlap correction | ✅ |
| SES covariate adjustment | ✅ |
| Combined IV selection for MVMR | ✅ |
| FDR correction (Benjamini-Hochberg) | ✅ |
| Method concordance validation | ✅ |
| Reverse MR (bidirectionality) | ✅ |
| Direction concordance checking | ✅ |
| F-statistics & Cochran's Q | ✅ |
| Conditional F-statistics (MVMR) | ✅ |
| **Total: 23/23 requirements** | ✅ **100%** |


## 🎓 Citation

If you use this pipeline in your research, please cite:

### Methods
```
Hemani G, et al. The MR-Base platform supports systematic causal inference 
across the human phenome. eLife. 2018. (TwoSampleMR)

Sanderson E, et al. An examination of multivariable Mendelian randomization 
in the single-sample and two-sample summary data settings. 
Int J Epidemiol. 2019. (MVMR)

Verbanck M, et al. Detection of widespread horizontal pleiotropy in causal 
relationships inferred from Mendelian randomization between complex traits 
and diseases. Nat Genet. 2018. (MR-PRESSO)

Mounier N, Kutalik Z. Bias correction for inverse variance weighting 
Mendelian randomization. Genet Epidemiol. 2023. (MRlap)
```

---

## 🔧 Installation

### Prerequisites

- R ≥ 4.0.0
- 8+ GB RAM (16+ GB recommended for MRlap)
- Internet connection (for package installation and MRlap reference data)

### Install All Packages

```r
# Automated installation (recommended)
source("Install_All_Packages.R")

# Or manual installation
install.packages(c("data.table", "readr", "tibble", "TwoSampleMR", 
                   "ieugwasr", "MVMR", "MRPRESSO"))
devtools::install_github("n-mounier/MRlap")
```

---

## 📚 Usage

### Step 1: Prepare Data

1. **Place GWAS files** in designated directories:
   - Exposures: `Standardized Circulating human plasma proteome_Data/` etc.
   - Mediators: `Cerebrospinal fluid metabolomics_Data/`
   - Outcomes: `Outcomes/`

2. **Download SES covariates** (REQUIRED for MVMR):
   - Education
   - Income
   - Occupation
   - Place in: `Covariates_SES/`

3. **Configure paths** in `Main Analysis.R` (lines 61-72)

### Step 2: Test Installation

```r
# Verify packages
source("Install_All_Packages.R")

# Test covariate files
source("Test_Covariate_Loading.R")

# Quick demo test (recommended before full analysis)
source("Demo_Test_Analysis.R")
```

### Step 3: Run Analysis

```r
# Full analysis (may take hours depending on data size)
source("Main Analysis.R")
```

### Step 4: Filter Results

```r
# Automated filtering
source("Results_Filter_Helper.R")
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

All results saved in `results_trial/`:

| File | Contents | Columns |
|------|----------|---------|
| `uvmr_comprehensive_results.csv` | UVMR results | ~35 |
| `mvmr_comprehensive_results.csv` | MVMR results | ~20 |
| `mediation_comprehensive_results.csv` | Mediation results | ~32 |

### Key Columns

**UVMR**:
- Standard MR: `b`, `se`, `pval`
- FDR: `q_value`
- Validation: `concordance_validated`, `F_statistic`
- MRlap: `mrlap_corrected_effect`, `mrlap_corrected_pval`, `mrlap_diff_pval`

**MVMR**:
- Effect: `beta_mvmr`, `se_mvmr`, `pval_mvmr`, `q_value`
- **Adjustment**: `adjusted_for` (must include SES covariates), `n_covariates` (must be 3)
- Methods: `method` (MVMR-IVW, Median, Egger, Lasso)

**Mediation**:
- Pathways: `beta_exp_med`, `beta_med_out_direct`, `beta_exp_out_total`
- **Validation**: `bidirectional` (must be "No_Unidirectional"), `direction_concordant` (must be TRUE)
- Proportion: `mediation_proportion` (0-1)

---

## 🎯 Quality Control

### Result Filtering Criteria

Following Ye et al. 2024 standards:

**UVMR - Causal Associations**:
```r
# IVW estimates considered causal only if:
# - Same direction AND significance as ≥1 sensitivity method
robust_uvmr <- uvmr[
  grepl("IVW", method) &
  q_value < 0.05 &
  concordance_validated == TRUE &
  F_statistic > 10 &
  (is.na(Q_pval) | Q_pval > 0.05)  # No heterogeneity
]
```

**MVMR - Independent Effects**:
```r
# Must be adjusted for SES covariates
independent_mvmr <- mvmr[
  method == "MVMR-IVW" &
  q_value < 0.05 &
  n_covariates == 3  # CRITICAL: must be 3
]
```

**Mediation - Valid Pathways** (all 4 criteria):
```r
valid_mediation <- mediation[
  pval_exp_med < 0.05 &                        # Criterion 1
  pval_med_out_direct < 0.05 &                 # Criterion 2
  bidirectional == "No_Unidirectional" &       # Criterion 3
  direction_concordant == TRUE                 # Criterion 4
]
```

---

## 💡 Advanced Features

### MRlap Auto-Detection

MRlap automatically:
- Extracts sample sizes from GWAS files
- Estimates sample overlap using LDSC
- Corrects for winner's curse and weak instruments
- No manual configuration needed!

### MVMR Combined IV Selection

Following standard practice:
> "Genetic instruments were the combination of SNPs, which had genome-wide significance in **either** the GWAS of each exposure **or** the GWAS of each mediator"

Implementation: `select_combined_ivs_mvmr()` function

### Comprehensive Validation

- **FDR correction**: Controls Type I error in multiple testing
- **Method concordance**: Validates IVW with sensitivity analyses
- **Reverse MR**: Excludes bidirectional relationships
- **Direction checking**: Ensures pathway consistency

---


## 🌟 Acknowledgments

**Implements methods from**:
- [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)
- [MVMR](https://github.com/WSpiller/MVMR)
- [MR-PRESSO](https://github.com/rondolab/MR-PRESSO)
- [MRlap](https://github.com/n-mounier/MRlap)

---

- **Version**: 2.5 FINAL  
- **Last Updated**: 2025-10-29
- **Contributor**: Hexiao Ding (PolyU); Jing Lan (PolyU); Na Li (SCU)
- **Superviosr**: Dr. Jung Sun Yoo (PolyU)
- **Institute**: The Hong Kong Polytechnic University, Sichuan University



