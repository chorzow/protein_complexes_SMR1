#--------0. ATTACHING PACKAGES--------
pckgs <- c('data.table', 'tidyr', 'dplyr', 'clusterProfiler', 'ggplot2', 'viridis')
for(i in pckgs){
  if(!require(i, character.only = T)){
    install.packages(i, dependencies = T)
    library(i)
  }
}


#--------1. FUNCTIONS DEFINITION--------
  data_prep <- function(path, expression_results, cor = 'CORUMcoreComplexes.txt', sorting_variable){
    # path: absolute path to folder with files including last slash. Do not specify the files themselves. Example: '/home/chorzow/SMR1/'
    # expression_results: filename of aberrant expression results. Should be provided in quotes. Example: 'expression_results.tsv'
    # cor: file containing the database of protein complexes. Should be provided in quotes. Default: 'CORUMcoreComplexes.txt'
    # sorting_variable: parameter by which sorting in downstream analysis (GSEA) will be performed. Should be one of 'fc', 'log2fc', 'zscore', 'padj'.
    if (sorting_variable == 'log2fc'){
      sorting_variable <- 'PROTEIN_LOG2FC'
    } else if (sorting_variable == 'zscore'){
      sorting_variable <- 'PROTEIN_ZSCORE'
    } else if (sorting_variable == 'padj'){
      sorting_variable <- 'PROTEIN_PADJ'
    } else if (sorting_variable == 'fc'){
      sorting_variable <- 'PROTEIN_FC'
    } else {
      stop("Wrong sorting variable!")  # associating chosen sorting variable with corresponding column names in aberrant expression file. Change names if necessary.
    }
    #--------1.1 READING CORUM DATABASE--------
    corum <- fread(paste0(path, 'CORUMcoreComplexes.txt')) %>% 
      filter(Organism == 'Human') %>%  # filter out unnecessary organisms
      rename(geneID = 'subunits(Gene name)') %>%  # rename one of the variables for convenience
      distinct()  # deduplicate
    
    corum <- separate_rows(corum, geneID) %>% 
      distinct() %>% 
      filter(!geneID %in% c('', ' ', '1', 'A')) %>%  # remove broken geneIDs
      as.data.table()  # convert to data.table
    corum[, N_subunits := .N, by = .(ComplexName)]  # create a variable with the number of subunits in each complex
    corum <- corum[N_subunits > 1]  # filter out all complexes with one subunit
    corum$geneID <- toupper(corum$geneID) 
    
    #--------2. READING AND FILTERING PROTEOMICS DATA--------
    prot_raw <- fread(paste0(path, expression_results))  # read raw aberrant expression data
    filt1 <- c(which(names(prot_raw) == 'proteome_ID'),
           which(names(prot_raw) == 'geneID'),
           which(names(prot_raw) == 'UNIPROT_ID'),
           which(names(prot_raw) == sorting_variable),
           which(names(prot_raw) == 'PROTEIN_outlier'),
           which(names(prot_raw) == 'validated'))
    prot <- select(prot_raw, filt1)  # select only necessary columns
    prot <- prot[!is.na(prot$UNIPROT_ID)  & !is.na(prot[,which(colnames(prot_raw) == sorting_variable)]) ]  # remove possible NAs
    prot <- prot[ !duplicated(prot) ]  # deduplicate
    prot <- prot[order(get(noquote(sorting_variable)))]  # sort by chosen sorting_variable
    prot <- prot[ !duplicated(prot[, c("proteome_ID", "geneID")]) ]  # deduplicate
    
    prot_res <- merge(corum, prot, by = "geneID", allow.cartesian=TRUE)  # map CORUM proteins to aberrant expression data by gene ID
    prot_res <- prot_res[ !duplicated(prot_res) ]  # deduplicate
    prot_res[, N_quantified_subunits:= .N, by = .(ComplexName, proteome_ID)]  # create a variable with the number of complex subunits detected in the data
    prot_res <- prot_res[ N_quantified_subunits > 1 ]  # filter out observations with only one dysregulated subunit
    prot_res <- prot_res[ N_quantified_subunits >= N_subunits / 2 ]  # filter out observations where number of dysregulated subunits is less than a half of total number of subunits in complex 
    
    #--------3. PREPARING DATA FOR GSEA--------
    filt2 <- c(which(names(prot_res) == 'proteome_ID'),
           which(names(prot_res) == 'geneID'),
           which(names(prot_res) == sorting_variable))
    gsea_p <- select(prot_res, filt2)  # select only necessary columns
    gsea_p <- gsea_p[order(get(noquote(sorting_variable)), decreasing = T)]  # sort by chosen sorting_variable
    gsea_p <- gsea_p[ !duplicated(gsea_p[,  c("proteome_ID", "geneID")]) ]  # deduplicate by patient ID and gene ID
    
    return (list(gsea_p, corum, sorting_variable))
  }
  
  gene_set_enrichment <- function(gsea_data, corum, sorting_variable, detection_method){  # performs GSEA itself 
    # gsea_data: data file with prepared results (obtained from the data_prep function)
    # corum: file containing the database of protein complexes (obtained from the data_prep function)
    # sorting_variable: parameter by which sorting in GSEA will be performed (obtained from the data_prep function)
    # detection_method: method of aberrant expression detection. Should be specified in quotes and chosen according to the method used for obtaining aberrant expression results data. For now, 'PROTRIDER' and 'LIMMA' are supported.
    
    res_GSEA <- data.frame()  # make empty data frame to store the results
    
    for (sample in unique( gsea_data$proteome_ID) ){  # for each unique patient ID in gsea_data
      samp <- gsea_data[proteome_ID == sample]  # add patient ID to the list
      geneList <- samp[[sorting_variable]]  
      names(geneList) <- as.character(samp$geneID)
      geneList <- sort(geneList, decreasing = TRUE)  # prepare gene list for GSEA
      ewp <- GSEA(geneList, TERM2GENE = corum[, c(2,18)], pvalueCutoff = 0.1, verbose=FALSE)  # GSEA from the clusterProfiler package
      
      ewp <- ewp@result  # extract the results
      ewp$proteome_ID <- rep(sample, nrow(ewp))
      ewp$Analysis <- rep(toupper(detection_method), nrow(ewp))
      ewp$Method <- rep("GSEA", nrow(ewp))
      ewp$Dataset <- rep("CORUM", nrow(ewp))
      res_GSEA <- rbind(res_GSEA, ewp)  # add variables to the resulting data frame and store it in res_GSEA
    }
  
    res_GSEA <- as.data.table( res_GSEA ) %>%  #convert to data.table
      filter(p.adjust < 0.05) %>%  # filter statistically significant enrichments
      mutate(ID = NULL)
    res_g <- res_GSEA[, c("proteome_ID", "Analysis", "Method", "Dataset", "Description", "NES", "p.adjust", "core_enrichment" )]  # choose necessary columns
    setnames(res_g, "core_enrichment", "genes")  # rename one variable for convenience
    
    res_g <- res_g[order(proteome_ID)]  # order by patient ID
    return (res_g)
  }
  
  Enrichment <- function(path, expression_results, sorting_variable = 'zscore', detection_method){  # wrapper function for data_prep and gene_set_enrichment
    # path: absolute path to folder with files including last slash. Do not specify the files themselves. Example: '/home/chorzow/project_data/'
    # expression_results: filename of aberrant expression results. Should be provided in quotes. Example: 'expression_results.tsv'
    # sorting_variable: parameter by which sorting in GSEA will be performed. Should be one of 'fc', 'log2fc', 'zscore', 'padj'. Default: 'zscore'.
    # detection_method: method of aberrant expression detection. Should be specified in quotes and chosen according to the method used for obtaining aberrant expression results data. For now, 'PROTRIDER' and 'LIMMA' are supported.
    prot_data <- data_prep(path = path, 
                           expression_results = expression_results,
                           sorting_variable = sorting_variable)
    return (gene_set_enrichment(prot_data[[1]], prot_data[[2]], prot_data[[3]], detection_method))
  }
  
#--------2. EXAMPLE ANALYSIS--------
  #0. Obtain the results for all combinations of sorting variables and aberrant detection methods
  
  #NOTE: in this example, all analysis was done for all protein complexes that are present in the aberrant expression.
  #  If you want to make the analysis for the patients for which the disease-causing gene is known, please uncomment all commented lines before the "Benchmark" section.
    
  path_to_files <- '/home/chorzow/BI/SMR1/'
  # protann_known <- fread(paste0(path_to_files, 'Proteomics_annotation.tsv')) %>% filter(!is.na(KNOWN_MUTATION))

  protrider_fc_known <- Enrichment(path = path_to_files, 
             expression_results = '3b_Protrider_results_autoencoder_normalised.tsv',
             sorting_variable = 'fc', detection_method = 'PROTRIDER') %>% 
    # filter(proteome_ID %in% protann_known$proteome_ID) %>% 
    mutate(sort_var = 'fc')
  
  protrider_log2fc_known <- Enrichment(path = path_to_files, 
                                   expression_results = '3b_Protrider_results_autoencoder_normalised.tsv',
                                   sorting_variable = 'log2fc', detection_method = 'PROTRIDER') %>% 
    # filter(proteome_ID %in% protann_known$proteome_ID) %>% 
    mutate(sort_var = 'log2fc')
  
  protrider_zscore_known <- Enrichment(path = '/media/chorzow/Data/current_projects/bioinf_institute/SMR1/', 
                                   expression_results = '3b_Protrider_results_autoencoder_normalised.tsv',
                                   sorting_variable = 'zscore', detection_method = 'PROTRIDER') %>% 
    # filter(proteome_ID %in% protann_known$proteome_ID) %>% 
    mutate(sort_var = 'zscore')
  
  protrider_padj_known <- Enrichment(path = '/media/chorzow/Data/current_projects/bioinf_institute/SMR1/',
                                     expression_results = '3b_Protrider_results_autoencoder_normalised.tsv',
                                     sorting_variable = 'padj', detection_method = 'PROTRIDER') %>% 
    # filter(proteome_ID %in% protann_known$proteome_ID) %>% 
    mutate(sort_var = 'padj')
  
  limma_fc_known <- Enrichment(path = '/media/chorzow/Data/current_projects/bioinf_institute/SMR1/', 
                                   expression_results = '3a_LIMMA_results_annotated.tsv',
                                   sorting_variable = 'fc', detection_method = 'LIMMA') %>% 
    # filter(proteome_ID %in% protann_known$proteome_ID) %>% 
    mutate(sort_var = 'fc')
  
  limma_log2fc_known <- Enrichment(path = '/media/chorzow/Data/current_projects/bioinf_institute/SMR1/', 
                               expression_results = '3a_LIMMA_results_annotated.tsv',
                               sorting_variable = 'log2fc', detection_method = 'LIMMA') %>% 
    # filter(proteome_ID %in% protann_known$proteome_ID) %>% 
    mutate(sort_var = 'log2fc')
  
  limma_zscore_known <- Enrichment(path = '/media/chorzow/Data/current_projects/bioinf_institute/SMR1/', 
                               expression_results = '3a_LIMMA_results_annotated.tsv',
                               sorting_variable = 'zscore', detection_method = 'LIMMA') %>% 
    # filter(proteome_ID %in% protann_known$proteome_ID) %>% 
    mutate(sort_var = 'zscore')
  
  limma_padj_known <- Enrichment(path = '/media/chorzow/Data/current_projects/bioinf_institute/SMR1/', 
                                   expression_results = '3a_LIMMA_results_annotated.tsv',
                                   sorting_variable = 'padj', detection_method = 'LIMMA') %>% 
    # filter(proteome_ID %in% protann_known$proteome_ID) %>% 
    mutate(sort_var = 'padj')
  
  #--------BENCHMARK--------
  
  # NOTE: Sensitivity benchmarking did not include sorting by p(adjusted), because the number of detected complexes obtained using this variable is too big and therefore useless for the biological analysis.
  # however, to show the abundance of complexes detected by p(adjusted), we showed this variable on the first plot.
  
  # 1. Number of detected complexes by each sorting variable
  complex_number <- rbind(protrider_fc_known, protrider_log2fc_known, protrider_padj_known, protrider_zscore_known, limma_fc_known, limma_log2fc_known, limma_padj_known, limma_zscore_known) %>% 
    select(1, 2, 9) %>% rename(Detection = Analysis) %>% 
    mutate_if(is.character, as.factor) %>% 
    group_by(proteome_ID, sort_var) %>% 
    mutate(N_complex = n()) %>%
    unique()  # resulting data frame with numbers of detected complexes for each sorting variable
  
  # 2. Sensitivity of different sorting algorithms: -log10(padj) by sorting variables
  
  sensitivity <- rbind(protrider_fc_known, protrider_log2fc_known, protrider_zscore_known, limma_fc_known, limma_log2fc_known, limma_zscore_known) %>% 
    select(1,2,7,9) %>% 
    mutate_if(is.character, as.factor) %>% 
    mutate(p.adjust = -log10(p.adjust)) %>% 
    rename(log10padj = p.adjust, Detection = Analysis)  # resulting data frame with -log10(p-values) for all complexes
 
  # 3. Visualizing the results
  
  ggplot(complex_number, aes(x = toupper(sort_var), y = N_complex, fill = Detection)) + 
    geom_boxplot(position = position_dodge(0.9)) + geom_jitter(size = 0.3, position = position_jitterdodge(dodge.width = 0.9)) + scale_fill_viridis(discrete = T, alpha = 0.6) + 
    theme_bw() + xlab('Sorting variable') + ylab('Number of complexes') + labs(fill = 'Detection\nmethod') + 
    theme(plot.title = element_text(hjust = 0.5, size = 16)) + ggtitle('Number of detected protein complexes') 
  
  ggplot(sensitivity, aes(toupper(sort_var), log10padj)) + 
    geom_boxplot(aes(fill = Detection)) + labs(fill = 'Detection\nmethod') + 
    scale_fill_viridis(discrete = T, alpha = 0.6) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16)) + 
    xlab('Sorting variable') + ylab('-log10(padj)') + ggtitle('Sensitivity benchmarking')
  
  
 #--------cnetplots--------  
PlotCORUMComplexes <- function(path, expression_results, sorting_variable, ID){  # perform part of the analysis and draw a cnetplot for specific patient
  # path: absolute path to folder with files including last slash. Do not specify the files themselves. Example: '/home/chorzow/project_data/'
  # expression_results: filename of aberrant expression results. Should be provided in quotes. Example: 'expression_results.tsv'
  # sorting_variable: parameter by which sorting in GSEA will be performed. Should be one of 'fc', 'log2fc', 'zscore', 'padj'.
  # ID: patient ID for which cnetplot will be drawn. Should be specified in quotes. Example: 'exPatient1'
  
  dp <- data_prep(path = path, 
                     expression_results = expression_results, sorting_variable = sorting_variable)
  expression_data <- dp[[1]]
  corum <- dp[[2]]
  samp <- expression_data[proteome_ID == ID]
  geneList <- samp[[3]]
  names(geneList) <- as.character(samp$geneID)
  geneList <- sort(geneList, decreasing = T)
  
  ewp <- GSEA(geneList, TERM2GENE = corum[, c(2, 18)], pvalueCutoff = 0.1, verbose = F)
  
  return (cnetplot(ewp, foldChange = geneList, layout = 'kk'))
}
  
# NOTE: for ethical reasons, real patient ID was substitited with 'exPatient1'

p3 <- PlotCORUMComplexes(path = path_to_files,
                         expression_results = '3b_Protrider_results_autoencoder_normalised.tsv',
                         sorting_variable = 'zscore', ID = 'P96993')
p3 + scale_color_gradient(low='red', high='blue') + labs(color = 'Z-score') + theme(plot.title = element_text(hjust = 0.5, size = 18), 
                                                                                    plot.subtitle = element_text(hjust = 0.5, size = 16)) + 
  ggtitle('Dysregulated protein complexes', subtitle = 'exPatient1')

