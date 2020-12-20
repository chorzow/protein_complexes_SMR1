This repository contains results of the semester project "Detecting novel molecular events in proteomics data for genetic diagnostics" in Bioinformatics Institute.

# Aim
The aim of this project was to improve diagnostics of Mendelian disorders via proteomics. 

# Objectives
1. To review existing literature on the subject.
2. To acquire expertise in different methods of detection of expression outliers detection
3. To modify these methods for detection of aberrant expression of protein complexes
4. To acquire expertise in gene set enrichment analysis (GSEA)
5. To perform benchmarking of the different combinations of methods and sorting parameters used in GSEA
6. To find optimal combination of the aforementioned methods and parameters

# Methods
We analyzed data of 147 patients, for 61 of which one of the genes containing disease-causing mutation was determined earlier. We used CORUM database of protein complexes avaliable at http://mips.helmholtz-muenchen.de/corum/#download. For detection of aberrant expression, two methods were used. First of them, LIMMA, implies normalization on the control sample. Second, PROTRIDER, normalizes the data using autoencoder. To determine the enriched protein complexes and their function, we used gene set enrichment analysis (GSEA) implemented in R package `clusterProfiler`. It is a computational method that determines whether an a priori defined set of genes shows statistically significant, concordant differences between two biological states (e.g. phenotypes) - see more detailed information at https://www.gsea-msigdb.org/gsea/index.jsp

# System requirements
The commands and examples mentioned in this README have been tested on x86_64 Ubuntu 18.04 LTS with Intel(R) Core(TM) i7-3630QM CPU, 8 Gb system memory using R version 3.6.3 with the following packages:
- `data.table` v1.13.0
- `tidyr` v1.1.2
- `dplyr` v1.0.2
- `clusterProfiler` v3.14.3
- `ggplot2` v3.3.2
- `viridis` v0.5.1

The aforementioned packages are also required to launch the code on your computer.

# Data requirements
1. CORUM complexes database. Can be downloaded from [here](http://mips.helmholtz-muenchen.de/corum/#download).
2. Aberrant protein expression data to analyze. This should be the data obtained with either LIMMA or PROTRIDER  method. Columns in such files should include the following columns:
* proteome_ID - IDs of patients
For each dysregulated protein, the following columns should be specified:
* geneID - gene IDs
* UNIPROT_ID - UNIPROT IDs 
* PROTEIN_ZSCORE - Z-scores 
* PROTEIN_FC - Fold change
* PROTEIN_LOG2FC - log2(fold change)
* PROTEIN_PADJ - adjusted p-value
* PROTEIN_outlier - TRUE if the protein is considered as outlier; FALSE if not
* validated - TRUE if the protein has the validation in the previously obtained data.

# Usage
The main script `SMR1` is divided into several sections. The following overview will cover them consequently.

## 0. Attaching packages
This section installs the packages necessary for the analysis. The `viridis` package was chosen because of aesthetic preferences. 

## 1. Function definitions
Three main functions were used in the analysis. `data_prep` imports the CORUM database and aberrant expression results, filters them and returns the data frame prepared for GSEA, filtered CORUM database and name of the sorting variable chosen for ordering the results in GSEA. `gene_set_enrichment` performs GSEA on the data frame returned by `data_prep` and returns the results. `Enrichment` is the wrapper function for the two functions described above. 

### data_prep
This function takes four arguments as input:
* `path`: absolute path to folder with files including last slash. Do not specify the files themselves. **Example:** `'/home/chorzow/BI/SMR1/'`
* `expression_results`: filename of aberrant expression results. Should be provided in quotes. **Example:** `'expression_results.tsv'`
* `cor`: file containing the database of protein complexes. Should be provided in quotes. **Default:** `'CORUMcoreComplexes.txt'`
* `sorting_variable`: parameter by which sorting in downstream analysis (GSEA) will be performed. Should be one of `'fc'`, `'log2fc'`, `'zscore'`, `'padj'`.

The function returns a list containing several data frames, stored at the following indexes: 1 - data frame prepared for the GSEA; 2 - filtered CORUM database; 3 - chosen sorting variable that will be used for GSEA. You can extract the data frames by assigning the elements of this list to various variables (See Example).

**Example:** 
```r
dp <- data_prep(path = '/home/chorzow/BI/SMR1/, 
                     expression_results = 'expression_results_LIMMA.tsv', sorting_variable = 'zscore')
expression_data <- dp[[1]]
corum <- dp[[2]]
sort_var <- dp[[3]]
```
### gene_set_enrichment
This function performs the Gene Set Enrichment Analysis (GSEA). It takes four arguments as input:
* gsea_data: data file with prepared results (obtained from the data_prep function)
* corum: file containing the filtered database of protein complexes (obtained from the data_prep function)
* sorting_variable: parameter by which sorting in GSEA will be performed (obtained from the data_prep function)
* detection_method: method of aberrant expression detection. Should be specified in quotes and chosen according to the method used for obtaining aberrant expression results data. For now, `'PROTRIDER'` and `'LIMMA'` are supported.

The function returns a data frame containing the final results of GSEA, ordered by patient ID.

**Example:**
```r
res <- gene_set_enrichment(gsea_data = dp[[1]],
                           corum = dp[[2]],
                           sorting_variable = dp[[3]],
                           detection_method = 'LIMMA')
```

### Enrichment
Since it could be inconvenient to call `data_prep` and `gene_set_enrichment` every time, we implemented a wrapper function for the previous two called `Enrichment`. This function takes four arguments as input:
* `path`: absolute path to folder with files including last slash. Do not specify the files themselves. **Example:** `'/home/chorzow/BI/SMR1/'`
* `expression_results`: filename of aberrant expression results. Should be provided in quotes. **Example:**`'expression_results_LIMMA.tsv'`
* `sorting_variable`: parameter by which sorting in GSEA will be performed. Should be one of `'fc'`, `'log2fc'`, `'zscore'`, `'padj'`. **Default:** `'zscore'`.
* `detection_method`: method of aberrant expression detection. Should be specified in quotes and chosen according to the method used for obtaining aberrant expression results data. For now, `'PROTRIDER'` and `'LIMMA'` are supported.

**Example:**
```r
protrider_fc <- Enrichment(path = '/home/chorzow/BI/SMR1/', 
             expression_results = 'Protrider_results_autoencoder_normalised.tsv',
             sorting_variable = 'fc', detection_method = 'PROTRIDER')
```

## 2. Example analysis
In this section, we performed benchmarking of different aberrant detection methods with different sorting variables in GSEA. After that, we visualized number of protein complexes detected in each patient using different combinations of sorting variables and aberrant expression detection methods:

![Number of detected complexes](https://github.com/chorzow/protein_complexes_SMR1/blob/main/images/Compl_num.png)

We noticed that the number of detected complexes was too big if we used adjusted p-value as a sorting variable. That is why we did not include this variable in the further analysis. We then visualized the sensitivity of each combination by plotting -log10 of adjusted p-value over different sorting variables:

![Sensitivity benchmarking](https://github.com/chorzow/protein_complexes_SMR1/blob/main/images/Sensitivity.png)

Finally, we created some fancy cnetplots to visualize dysregulated protein complexes in some patients and extent of the dysregulation. For ethical reasons, we will not specify patient IDs neither in code nor on plots. 

![cnetplot example](https://github.com/chorzow/protein_complexes_SMR1/blob/main/images/patient_1.png)

# References

Zhou, B., Yan, Y., Wang, Y. et al. Quantitative proteomic analysis of prostate tissue specimens identifies deregulated protein complexes in primary prostate cancer. Clin Proteom 16, 15 (2019). https://doi.org/10.1186/s12014-019-9236-2
Jiang, L., Wang M., Lin S. et. al. A Quantitative Proteome Map of the Human Body. bioRxiv 797373; doi: https://doi.org/10.1101/797373
Kremer, L., Bader, D., Mertes, C. et al. Genetic diagnosis of Mendelian disorders via RNA sequencing. Nat Commun 8, 15824 (2017). https://doi.org/10.1038/ncomms15824
Vicente A. Yépez, Christian Mertes, Michaela F. Mueller et al. Detection of aberrant events in RNA sequencing data, 03 January 2020, PROTOCOL (Version 1) available at Protocol Exchange [+https://doi.org/10.21203/rs.2.19080/v1+]
