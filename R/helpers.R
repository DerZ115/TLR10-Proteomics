libraries <- c(
  "tidyverse", "magrittr", "readxl", "stringr", "MSstats",
  "SummarizedExperiment", "MsCoreUtils", "DEqMS", "qvalue", "ggrepel",
  "janitor", "ComplexHeatmap", "envalysis", "RColorBrewer",
  "patchwork", "circlize", "vsn", "ggforce"
)

lapply(libraries, library, character.only = TRUE)


#' Read DIA Report
#'
#' This function reads a DIA (Data Independent Acquisition) report from a specified file path,
#' processes the data, and returns a `SummarizedExperiment` object containing the MS1 and MS2
#' quantification data, annotation data, and sample data.
#'
#' @param path A character string specifying the file path to the DIA report.
#' @return A `SummarizedExperiment` object containing:
#' \describe{
#'   \item{MS1}{A matrix of MS1 quantification data.}
#'   \item{MS2}{A matrix of MS2 quantification data.}
#'   \item{rowData}{A data frame of annotations including protein names, protein groups, and peptide counts.}
#'   \item{colData}{A data frame of sample data including cell type, condition, and replicate information.}
#' }
#' @import dplyr
#' @import tidyr
#' @import SummarizedExperiment
#' @importFrom readr read_delim
#' @importFrom stringr str_split_i str_remove str_replace
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
#' # Example usage:
#' # result <- read_DIA_report("path/to/DIA_report.txt")
read_DIA_report <- function(path) {
  data <- read.delim(path, header = TRUE, sep = "\t", 
                     dec = ",", na.strings = c("NaN", "NA", "")) %>%
    mutate(PG.ProteinNames = str_split_i(PG.ProteinNames, ";", 1),
           PG.Genes = str_split_i(PG.Genes, ";", 1)) %>%
    select(c(PG.ProteinGroups,
             PG.Genes,
             PG.ProteinNames,
             PG.NrOfStrippedSequencesIdentified..Experiment.wide.,
             matches("MS[12]Quantity"))) %>%
    group_by(PG.Genes) %>%
    summarise(PG.ProteinNames = dplyr::first(PG.ProteinNames),
              PG.ProteinGroups = paste(PG.ProteinGroups, sep = ";", collapse = ""),
              across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA, sum(., na.rm = TRUE)))) %>%
    ungroup() %>%
    mutate(PG.Genes = coalesce(PG.Genes, PG.ProteinNames)) %>%
    set_colnames(c("Gene.Name", "Protein.Name", "Protein.Groups", "Peptides",
                   colnames(.)[5:ncol(.)] %>% 
                     str_remove("X\\.\\d+\\.\\.") %>% 
                     str_replace("\\.raw\\.PG\\.", "_") %>%
                     str_remove("Quantity") %>%
                     str_replace("_WO_", "_Dark_"))) %>%
  column_to_rownames("Gene.Name")
  
  quant_data_MS1 <- data %>% 
    select(matches("^(?:ADIPO_|OSTEO_)?ASC_(?:TERT_)?(?:WT|TLR10_LOV|TLR4_LOV)_(?:SN|WCL)_(?:Dark|Light|LPS)_\\d+(?:_\\d{8,})?_MS1$")) %>%
    set_colnames(str_remove(colnames(.), "_MS1$") %>%
                   str_remove("_\\d{8,}")) %>%
    as.matrix()
  
  quant_data_MS2 <- data %>% 
    select(matches("^(?:ADIPO_|OSTEO_)?ASC_(?:TERT_)?(?:WT|TLR10_LOV|TLR4_LOV)_(?:SN|WCL)_(?:Dark|Light|LPS)_\\d+(?:_\\d{8,})?_MS2$")) %>%
    set_colnames(str_remove(colnames(.), "_MS2$") %>%
                   str_remove("_\\d{8,}")) %>%
    as.matrix()
  
  annotation_data <- data %>%
    select(c(Protein.Name, Protein.Groups, Peptides))
  
  sample_data <- colnames(quant_data_MS1) %>% 
    data.frame(Sample.Name = .) %>% 
    separate_wider_regex(Sample.Name, c(Cell.Type = ".*", 
                                        "_[:alpha:]*_", 
                                        Condition = "[:alpha:]*", 
                                        "_", 
                                        Replicate = "\\d+"),
                         cols_remove = FALSE) %>%
    column_to_rownames("Sample.Name") %>%
    mutate(Cell.Type = factor(Cell.Type, unique(Cell.Type)),
           Condition = factor(Condition, unique(Condition)))
  
  SummarizedExperiment(list(MS1 = quant_data_MS1, MS2 = quant_data_MS2), rowData = annotation_data, colData = sample_data)
}


#' Plot Missing Values
#'
#' This function creates bar plots showing the count of missing values in MS1 and MS2 data frames.
#'
#' @param MS1.df A data frame containing MS1 quantification data.
#' @param MS2.df A data frame containing MS2 quantification data.
#' @return A combined plot of missing values for MS1 and MS2 data frames.
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @export
#'
#' @examples
#' # Example usage:
#' # p <- plot_missing(MS1.df, MS2.df)
#' # print(p)
plot_missing <- function(MS1.df, MS2.df) {
  missing_values_MS1 <- MS1.df %>%
    mutate(na_vals = rowSums(is.na(.)))
  missing_values_MS2 <- MS2.df %>%
    mutate(na_vals = rowSums(is.na(.)))
  
  xmax <- max(c(missing_values_MS1$na_vals, missing_values_MS2$na_vals))
  
  p1 <- ggplot(missing_values_MS1 %>% dplyr::count(na_vals), aes(x = na_vals, y = n)) + 
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust=-0.3) +
    scale_x_continuous(limits = c(-0.5, xmax + 0.5), breaks = seq(0, xmax, by = 1),
                       expand = expansion(add=0)) +
    scale_y_continuous(expand = expansion(mult=c(0, 0.15))) +
    labs(title = "MS1", x = "# Missing", y = "Count") +
    theme_publish()
  
  p2 <- ggplot(missing_values_MS2 %>% dplyr::count(na_vals), aes(x = na_vals, y = n)) + 
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust=-0.3) +
    scale_x_continuous(limits = c(-0.5, xmax + 0.5), breaks = seq(0, xmax, by = 1),
                       expand = expansion(add=0)) +
    scale_y_continuous(expand = expansion(mult=c(0, 0.15))) +
    labs(title = "MS2", x = "# Missing", y = "Count") +
    theme_publish()
  
  p1 / p2 + plot_layout(axis_titles = "collect")
}


#' Filter Genes with Too Many Missing Values
#'
#' This function filters out genes/proteins that have too many missing values in both MS1 and MS2 data frames.
#'
#' @param MS1.df A data frame containing MS1 quantification data.
#' @param MS2.df A data frame containing MS2 quantification data.
#' @param full.groups An integer specifying the minimum number of sample groups required to keep a gene.
#' @return A vector of gene names that pass the filter.
#' @import dplyr
#' @import tidyr
#' @export
#'
#' @examples
#' # Example usage:
#' # filtered_genes <- filter_too_many_missing(MS1.df, MS2.df, full.groups = 1)
filter_too_many_missing <- function(MS1.df, MS2.df, full.groups = 1) {
  keep1 <- MS1.df %>% pivot_longer(!Gene.Name, names_to = "Sample.Name", ) %>%
    mutate(value = !is.na(value), Sample.Name = str_remove(Sample.Name, "_\\d+$")) %>%
    group_by(Gene.Name, Sample.Name) %>%
    summarise(Group.Count = sum(value), .groups = "drop") %>%
    ungroup() %>%
    pivot_wider(names_from = Sample.Name, values_from = Group.Count) %>%
    filter(rowSums(across(!Gene.Name, ~ . == 3)) >= full.groups) %>%
    pull(Gene.Name)
  
  keep2 <- MS2.df %>% pivot_longer(!Gene.Name, names_to = "Sample.Name", ) %>%
    mutate(value = !is.na(value), Sample.Name = str_remove(Sample.Name, "_\\d+$")) %>%
    group_by(Gene.Name, Sample.Name) %>%
    summarise(Group.Count = sum(value), .groups = "drop") %>%
    ungroup() %>%
    pivot_wider(names_from = Sample.Name, values_from = Group.Count) %>%
    filter(rowSums(across(!Gene.Name, ~ . == 3)) >= full.groups) %>%
    pull(Gene.Name)
  
  intersect(keep1, keep2)
}


#' Plot Intensity Boxplot
#'
#' This function creates boxplots of log2 intensities for MS1 and MS2 data frames.
#'
#' @param MS1.df A data frame containing MS1 quantification data.
#' @param MS2.df A data frame containing MS2 quantification data.
#' @return A combined plot of intensity boxplots for MS1 and MS2 data frames.
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @export
#'
#' @examples
#' # Example usage:
#' # p <- plot_intensity_boxplot(MS1.df, MS2.df)
#' # print(p)
plot_intensity_boxplot <- function(MS1.df, MS2.df) {
  MS1.plot <- MS1.df %>%
    pivot_longer(!Gene.Name, names_to = "Sample", 
                 values_to = "Log2Intensity", values_drop_na = TRUE)
  
  MS2.plot <- MS2.df %>%
    pivot_longer(!Gene.Name, names_to = "Sample", 
                 values_to = "Log2Intensity", values_drop_na = TRUE)
  
  p1 <- ggplot(MS1.plot, aes(x = Sample, y = Log2Intensity)) +
    geom_boxplot() +
    labs(x = "", title = "MS1") +
    theme(axis.text.x = element_blank())
  
  p2 <- ggplot(MS2.plot, aes(x = Sample, y = Log2Intensity)) +
    geom_boxplot() +
    labs(x = "", title = "MS2") +
    theme(axis.text.x = element_blank())
  
  p1 / p2 + plot_layout(axes = "collect")
}
  

#' Get MNAR Genes
#'
#' This function identifies genes that are missing not at random (MNAR) in the data frame. 
#' It does this by checking if a gene/protein is not found in all replicates of a sample.
#'
#' @param df A data frame containing quantification data.
#' @return A logical vector indicating whether each gene is MNAR.
#' @import dplyr
#' @import tidyr
#' @export
#'
#' @examples
#' # Example usage:
#' # MNAR <- get_MNAR(df)
get_MNAR <- function(df) {
  MNAR_genes <- df %>%
    pivot_longer(!Gene.Name, names_to = "Sample.Name") %>%
    mutate(Sample.Name = str_remove(Sample.Name, "_\\d+$")) %>%
    group_by(Gene.Name, Sample.Name) %>%
    summarise(mnar = all(is.na(value)), .groups = "drop") %>%
    filter(mnar) %>%
    pull(Gene.Name) %>%
    unique()
  
  MNAR <- rownames(data.filtered) %in% MNAR_genes
  
  names(MNAR) <- rownames(data.filtered)
  
  MNAR
}


#' Plot Missing Values Heatmap
#'
#' This function creates a heatmap showing the missing values in MS1 and MS2 data frames.
#'
#' @param MS1.df A data frame containing MS1 quantification data.
#' @param MS2.df A data frame containing MS2 quantification data.
#' @param mnar A logical vector indicating whether each gene is MNAR.
#' @param samples A data frame containing sample information.
#' @param colors_rows A list of colors for row annotations.
#' @param colors_columns A list of colors for column annotations.
#' @return A heatmap of missing values.
#' @import dplyr
#' @import ComplexHeatmap
#' @import circlize
#' @export
#'
#' @examples
#' # Example usage:
#' # plot_missing_heatmap(MS1.df, MS2.df, mnar, samples)
plot_missing_heatmap <- function(MS1.df, MS2.df, mnar, samples, colors_rows = NULL, colors_columns = NULL) {
  missing_values_MS1 <- MS1.df %>%
    set_colnames(colnames(.) %>% str_remove("_MS1$")) %>%
    mutate(across(!Gene.Name, ~ is.na(.))) %>%
    mutate(na_vals = rowSums(across(!Gene.Name)), across(!c(Gene.Name, na_vals), as.numeric))
    
    
  missing_values_MS2 <- MS2.df %>%
    set_colnames(colnames(.) %>% str_remove("_MS2$")) %>%
    mutate(across(!Gene.Name, ~ is.na(.))) %>%
    mutate(na_vals = rowSums(across(!Gene.Name)), across(!c(Gene.Name, na_vals), ~2 * as.numeric(.x)))
  
  missing_values_hm <- bind_rows(missing_values_MS1,
                                 missing_values_MS2) %>%
    group_by(Gene.Name) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("Gene.Name") %>%
    filter(na_vals > 0) %>%
    select(!na_vals) %>%
    mutate(across(everything(), ~ case_match(.x,
                                             0 ~ "Neither",
                                             1 ~ "MS1",
                                             2 ~ "MS2",
                                             3 ~ "Both"))) %>%
    as.matrix()
  
  MNAR.hm <- mnar[names(mnar) %in% rownames(missing_values_hm)]
  
  dist_lookup <- matrix(c(0, 1, 1, 2, 
                          1, 0, 2, 1,  
                          1, 2, 0, 1,
                          2, 1, 1, 0),
                        nrow = 4,
                        dimnames = list(c("Neither", "MS1", "MS2", "Both"),
                                        c("Neither", "MS1", "MS2", "Both")))
  
  Heatmap(missing_values_hm, col = c(Neither="black", MS1="mediumblue", MS2="red3", Both="gray"),  
          name = "Missing Value", cluster_rows = TRUE, show_row_dend = FALSE, 
          clustering_distance_rows = function(x, y) {
            sum(dist_lookup[matrix(c(x, y), ncol = 2)])
          },
          cluster_row_slices = FALSE, show_column_names = FALSE, 
          show_row_names = FALSE, split = MNAR.hm, row_title = NULL, 
          column_title = NULL, cluster_columns = FALSE, 
          left_annotation = rowAnnotation(MNAR = MNAR.hm,
                                          col = colors_rows,
                                          show_annotation_name = FALSE),
          top_annotation = HeatmapAnnotation(Celltype = samples$Cell.Type,
                                             Condition = samples$Condition,
                                             col = colors_columns))
}


#' Plot Missing Values Density
#'
#' This function creates density plots showing the distribution of missing values in MS1 and MS2 data frames.
#'
#' @param MS1.df A data frame containing MS1 quantification data.
#' @param MS2.df A data frame containing MS2 quantification data.
#' @return A combined plot of missing values density for MS1 and MS2 data frames.
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @export
#'
#' @examples
#' # Example usage:
#' # p <- plot_missing_density(MS1.df, MS2.df)
#' # print(p)
plot_missing_density <- function(MS1.df, MS2.df) {
  missing_values_density_MS1 <- MS1.df %>%
    mutate(na_vals = rowSums(is.na(.)) > 0) %>%
    pivot_longer(!c(Gene.Name, na_vals), names_to = "Sample.Name", values_to = "Intensity") 
  
  missing_values_density_MS2 <- MS2.df %>%
    mutate(na_vals = rowSums(is.na(.)) > 0) %>%
    pivot_longer(!c(Gene.Name, na_vals), names_to = "Sample.Name", values_to = "Intensity") 
  
  missing_values_density <- bind_rows(list(missing_values_density_MS1,
                                           missing_values_density_MS2), 
                                      .id = "MS") %>%
    na.omit()
  
  p1 <- ggplot(missing_values_density, aes(x = Intensity, linetype = na_vals, color = MS)) +
    geom_density(linewidth = 1) +
    labs(x = "Log2Intensity", y = "Density") + 
    scale_linetype_manual(values = c("FALSE" = 1, "TRUE" = 2), labels = c("No", "Yes")) +
    theme_publish() + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(linewidth = 0.5))
  
  p2 <- ggplot(missing_values_density, aes(x = Intensity, linetype = na_vals, color = MS)) +
    stat_ecdf(linewidth = 1) +
    labs(x = "Log2Intensity", y = "Cumulative Fraction", color = "MS", linetype = "Has missing values") + 
    scale_linetype_manual(values = c("FALSE" = 1, "TRUE" = 2), labels = c("No", "Yes")) +
    theme_publish() +
    theme(panel.grid.major = element_line(linewidth = 0.5))
  
  p1 / p2 + plot_layout(axes = "collect")
  
}


#' Plot Imputation Density
#'
#' This function creates density plots comparing the intensity distributions before and after imputation.
#'
#' @param before A SummarizedExperiment object containing data before imputation.
#' @param after A SummarizedExperiment object containing data after imputation.
#' @return A combined plot of intensity density before and after imputation.
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' # Example usage:
#' # p <- plot_imputation_density(before, after)
#' # print(p)
plot_imputation_density <- function(before, after) {
  samples <- colData(after) %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample.Name")
  
  before.plot <- assays(before) %>%
    lapply(as.data.frame) %>%
    lapply(rownames_to_column, var = "Gene.Name") %>% 
    bind_rows(.id = "MS") %>%
    mutate(MS = as.factor(MS), na_vals = rowSums(is.na(across(where(is.numeric)))) > 0) %>%
    pivot_longer(!c(Gene.Name, MS, na_vals), names_to = "Sample.Name", values_to = "Intensity") %>%
    na.omit() %>%
    left_join(samples, by = "Sample.Name")
  
  after.plot <- assays(after) %>%
    lapply(as.data.frame) %>%
    lapply(rownames_to_column, var = "Gene.Name") %>%
    bind_rows(.id = "MS") %>%
    mutate(MS = as.factor(MS), na_vals = rowSums(is.na(across(where(is.numeric)))) > 0) %>%
    pivot_longer(!c(Gene.Name, MS, na_vals), names_to = "Sample.Name", values_to = "Intensity") %>%
    left_join(samples, by = "Sample.Name")
  
  
  p1 <- ggplot(before.plot, aes(x = Intensity,
                                color = MS)) +
    geom_density(aes(group = interaction(Cell.Type, Condition, sep = ":")),
                 linewidth = 1) +
    labs(title = "Before", x = "Log2Intensity", y = "Density", color = "MS") + 
    theme_publish() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(linewidth = 0.5))
  
  p2 <- ggplot(after.plot, aes(x = Intensity, 
                               color = MS)) +
    geom_density(aes(group = interaction(Cell.Type, Condition, sep = ":")),
                 linewidth = 1) +
    labs(title = "After", x = "Log2Intensity", y = "Density", color = "MS") + 
    theme_publish() +
    theme(panel.grid.major = element_line(linewidth = 0.5))
  
  p1 / p2 + plot_layout(axis_titles = "collect")
}


#' Plot PCA
#'
#' This function performs PCA on the data and creates PCA plots.
#'
#' @param data A SummarizedExperiment object containing quantification data.
#' @param n An integer specifying the number of top variable genes to use for PCA.
#' @param maxPC An integer specifying the maximum number of principal components to plot.
#' @param center A logical value indicating whether to center the data.
#' @param scale A logical value indicating whether to scale the data.
#' @param plot_all A logical value indicating whether to plot all principal components.
#' @return PCA plots.
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' # Example usage:
#' # plot_pca(data, n = 500, maxPC = 2)
plot_pca <- function(data, n = 500, maxPC = 2, center = TRUE, scale = TRUE, plot_all = F) {
  samples <- colData(data) %>% 
    as.data.frame() %>%
    rownames_to_column("Sample.Name")
  
  
  data.pca <- assays(data) %>%
    lapply(as.data.frame) %>%
    lapply(rownames_to_column, var = "Gene.Name") %>%
    bind_rows(.id = "MS") %>%
    unite(rowname, c("Gene.Name", "MS")) %>%
    column_to_rownames() %>%
    as.matrix() %>%
    t() %>%
    prcomp(center = center, scale. = scale) 
  
  data.pca.plot <- data.pca$x %>%
    as.data.frame() %>%
    rownames_to_column("Sample.Name") %>%
    left_join(samples, by = "Sample.Name")
  
  data.pca.var <- as.data.frame(t(summary(data.pca)$importance)) %>%
    set_colnames(c("sd", "var_prop", "var_cumprop")) %>%
    mutate(across(c(var_prop, var_cumprop), ~ .x * 100), Component = seq.int(1, nrow(.), 1))
  
  p1 <- ggplot(data.pca.var, aes(x = Component)) +
    geom_bar(aes(y = var_prop), stat = "identity") +
    geom_line(aes(y = var_cumprop)) +
    geom_point(aes(y = var_cumprop)) +
    scale_x_continuous(limits = c(0.5, nrow(data.pca.var) + 0.5), 
                       breaks = seq(1, nrow(data.pca.var), by = 1),
                       expand = expansion(add=0)) +
    scale_y_continuous(expand = expansion(mult=c(0, 0.05))) +
    ylab("Variance (%)") +
    theme_publish() +
    theme(panel.grid.major.y = element_line(linewidth = 0.5))
  
  print(p1)
  
  
  p2 <- ggplot(data.pca.plot, aes(x = PC1, y = PC2, color = interaction(Cell.Type, Condition, sep = ":"),
                             shape = Replicate)) +
    geom_point(size = 3) +
    labs(color = "Celltype:Condition") +
    theme_publish() +
    theme(panel.grid.major = element_line(linewidth = 0.5),
          legend.position = "right")
  
  print(p2)

  if (plot_all) {
    
    p3 <- ggplot(data.pca.plot, aes(color = interaction(Cell.Type, Condition, sep = ":"))) +
      geom_autopoint(aes(shape = Replicate), size = 2) +
      geom_autodensity(position = "identity", fill = NA) +
      facet_matrix(vars(paste0("PC", seq.int(1, maxPC, 1))), layer.diag = 2) +
      labs(color = "Celltype:Condition") +
      theme_bw() +
      theme(panel.grid.major = element_line(linewidth = 0.5))
    
    print(p3)
  }
    
    
    
  
}


#' Fit DEqMS Model
#'
#' This function fits a limma & DEqMS model to the quantification data.
#'
#' @param data A SummarizedExperiment object containing quantification data.
#' @param contrast_list A list of contrasts to test.
#' @return A DEqMS model fit object.
#' @import dplyr
#' @import limma
#' @import DEqMS
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' # Example usage:
#' # fit <- fit_DEqMS_model(data, contrast_list)
fit_DEqMS_model <- function(data, contrast_list) {
  samples <- colData(data) %>%
    as.data.frame() %>%
    rownames_to_column("Sample.Name") %>%
    bind_rows(., ., .id = "MS") %>%
    mutate(class = factor(paste(Cell.Type, Condition, sep = ".")),
           Sample.Name = paste(Sample.Name, MS, sep = "_MS")) %>%
    column_to_rownames("Sample.Name")
  
  
  data.merged <- merge(as.data.frame(assay(data, "MS1")),
                       as.data.frame(assay(data, "MS2")),
                       by = "row.names", suffixes = c("_MS1", "_MS2")) %>%
    column_to_rownames("Row.names") %>%
    as.matrix()
  

  DE_design <- model.matrix(~ 0 + samples$class + samples$Replicate + samples$MS)
  colnames(DE_design) <- c(levels(samples$class), "Replicate2", "Replicate3", "MS2")
  
  fit1 <- lmFit(data.merged, design = DE_design)
  cont_matrix <- makeContrasts(contrasts = contrast_list,
                               levels = DE_design)
  fit2 <- contrasts.fit(fit1, contrasts = cont_matrix)
  fit3 <- eBayes(fit2, trend = TRUE)
  fit3$count <- rowData(data)[rownames(fit3$coefficients),]$Peptides
  fit4 <- spectraCounteBayes(fit3)
  
  fit4
}


#' Plot Heatmap
#'
#' This function creates a heatmap of significant proteins.
#'
#' @param data.sig A SummarizedExperiment object containing significant data.
#' @param color_scale A color scale for the heatmap.
#' @param colors_columns A list of colors for column annotations.
#' @param qvalues A vector of q-values for the data.
#' @param title A character string specifying the title of the heatmap.
#' @param max_rows An integer specifying the maximum number of rows to display.
#' @param cluster_rows A logical value indicating whether to cluster rows.
#' @param cluster_cols A logical value indicating whether to cluster columns.
#' @param scale_by_row A logical value indicating whether to scale data by row.
#' @return A heatmap of significant data.
#' @import dplyr
#' @import ComplexHeatmap
#' @import circlize
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' # Example usage:
#' # plot_heatmap(data.sig, color_scale, colors_columns)
plot_heatmap <- function(data.sig, color_scale, colors_columns, qvalues = NA, 
                         title = NA, max_rows = 20, cluster_rows = TRUE, 
                         cluster_cols = TRUE, scale_by_row = TRUE) {
  hmap.annotation.row <- rowData(data.sig) %>%
    as.data.frame()
  
  hmap.annotation.col <- colData(data.sig) %>% 
    as.data.frame() %>%
    dplyr::select(Cell.Type, Condition)
  
  data.hm <- Reduce("+", assays(data.sig)) / 2
  
  if (scale_by_row) {
    data.hm %<>%
      t() %>%
      scale() %>%
      t()
  }
  
  if (nrow(data.hm) <= max_rows) {
    hm <- Heatmap(data.hm, col = color_scale,  
                  name = title, cluster_rows = cluster_rows, 
                  show_row_dend = TRUE, show_column_names = FALSE, show_row_names = TRUE, 
                  row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                  top_annotation = HeatmapAnnotation(Celltype = hmap.annotation.col$Cell.Type,
                                                     Condition = hmap.annotation.col$Condition,
                                                     col = colors_columns))
    
    draw(hm)
  } else {
    
    hm1 <- Heatmap(data.hm, col = color_scale,  
                   name = title, cluster_rows = cluster_rows, 
                   show_row_dend = FALSE, show_column_names = FALSE, show_row_names = FALSE, 
                   row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                   top_annotation = HeatmapAnnotation(Celltype = hmap.annotation.col$Cell.Type,
                                                      Condition = hmap.annotation.col$Condition,
                                                      col = colors_columns))
    
    data.top <- data.hm %>%
      as.data.frame() %>%
      mutate(qvalue = qvalues[rownames(data.hm)]) %>%
      arrange(qvalue) %>%
      head(max_rows) %>%
      select(!qvalue) %>%
      as.matrix()
    
    hmap.annotation.row %<>%
      filter(rownames(.) %in% rownames(data.top))
    
    hm2 <- Heatmap(data.top, col = color_scale,  
                   name = title, cluster_rows = cluster_rows, 
                   show_row_dend = TRUE, show_column_names = FALSE, show_row_names = TRUE, 
                   row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                   top_annotation = HeatmapAnnotation(Celltype = hmap.annotation.col$Cell.Type,
                                                      Condition = hmap.annotation.col$Condition,
                                                      col = colors_columns))
    
    draw(hm1)
    draw(hm2)
  }
}


#' Plot Volcano
#'
#' This function creates a volcano plot of differential expression results.
#'
#' @param res A data frame containing differential expression results.
#' @param title A character string specifying the title of the plot.
#' @param lfc_limit A numeric value specifying the limit for log fold change.
#' @return A volcano plot.
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @export
#'
#' @examples
#' # Example usage:
#' # plot_volcano(res, title = "Volcano Plot")
plot_volcano <- function(res, title = "", lfc_limit = NA) {
  res.plot <- res %>%
    filter(!is.na(sca.qvalue)) %>%
    arrange(sca.qvalue) %>%
    rowwise() %>%
    mutate(threshold = sca.qvalue < 0.05 & abs(logFC) > 0.58,
           out_of_bounds = ifelse(is.na(lfc_limit), 0, (abs(logFC) > lfc_limit) * sign(logFC)),
           logFC_capped = ifelse(out_of_bounds != 0, lfc_limit * sign(logFC), logFC)) %>%
    ungroup()
  
  res.sig <- res %>%
    filter(sca.qvalue < 0.05 & abs(logFC) >= 0.58)
  

  res.plot %<>% 
    mutate(genelabels = ifelse(threshold, gene, ""))

  
  maxFC <- max(abs(res.plot$logFC))
  if (!is.na(lfc_limit) && lfc_limit < maxFC) {
    xlim <- lfc_limit
  } else {
    xlim <- maxFC * 1.04
  }
  
  ggplot(res.plot, aes(x = logFC_capped, y = sca.qvalue)) +
    geom_point(data = subset(res.plot, out_of_bounds == 0), aes(colour=threshold), alpha = 0.5) +
    geom_point(data = subset(res.plot, out_of_bounds == -1), aes(colour=threshold), shape = "\u25c4", size=2) +
    geom_point(data = subset(res.plot, out_of_bounds == 1), aes(colour=threshold), shape = "\u25BA", size=2) +
    geom_hline(yintercept = 0.05, linetype = 2) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = 2) +
    geom_text_repel(aes(label = genelabels)) +
    scale_x_continuous(breaks = scales::pretty_breaks(), limits = c(-xlim, xlim), expand = expansion(0.01)) +
    scale_y_continuous(trans = c("log10", "reverse"), breaks = scales::log_breaks(), labels = scales::scientific) +
    ggtitle(title) +
    xlab("log2 Fold Change") +
    ylab("q-value") +
    theme_publish() +
    theme(legend.position = "none")
}