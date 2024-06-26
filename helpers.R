
libraries <- c("tidyverse", "magrittr", "readxl", "stringr","MSstats", 
               "SummarizedExperiment", "MsCoreUtils", "DEqMS", "ggrepel",
               "janitor", "ComplexHeatmap", "envalysis", "RColorBrewer", 
               "patchwork")

lapply(libraries, library, character.only = TRUE)

read_DIA_report <- function(path) {
  data <- read.delim(path, header = TRUE, sep = "\t", dec = ",", na.strings = c("NaN", "NA", "")) %>%
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
    select(matches("^(?:ADIPO_|OSTEO_)?ASC_(?:WT|TLR10_LOV)_(?:SN|WCL)_(?:Dark|Light)_\\d+_MS1$")) %>%
    set_colnames(str_remove(colnames(.), "_MS1$")) %>%
    as.matrix()
  
  quant_data_MS2 <- data %>% 
    select(matches("^(?:ADIPO_|OSTEO_)?ASC_(?:WT|TLR10_LOV)_(?:SN|WCL)_(?:Dark|Light)_\\d+_MS2$")) %>%
    set_colnames(str_remove(colnames(.), "_MS2$")) %>%
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
           Condition = factor(Condition, c("Dark", "Light")))
  
  SummarizedExperiment(list(MS1 = quant_data_MS1, MS2 = quant_data_MS2), rowData = annotation_data, colData = sample_data)
}


plot_missing <- function(MS1.df, MS2.df) {
  missing_values_MS1 <- MS1.df %>%
    mutate(na_vals = rowSums(is.na(.)))
  missing_values_MS2 <- MS2.df %>%
    mutate(na_vals = rowSums(is.na(.)))
  
  p1 <- ggplot(missing_values_MS1 %>% dplyr::count(na_vals), aes(x = na_vals, y = n)) + 
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust=-0.3) +
    scale_x_continuous(breaks = seq(0, max(c(missing_values_MS1$na_vals, 
                                             missing_values_MS2$na_vals)), 
                                    by = 1)) +
    scale_y_continuous(expand = expansion(mult=c(0, 0.15))) +
    labs(title = "MS1", x = "# Missing", y = "Count")
  
  p2 <- ggplot(missing_values_MS2 %>% dplyr::count(na_vals), aes(x = na_vals, y = n)) + 
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust=-0.3) +
    scale_x_continuous(breaks = seq(0, max(c(missing_values_MS1$na_vals, 
                                             missing_values_MS2$na_vals)), 
                                    by = 1)) +
    scale_y_continuous(expand = expansion(mult=c(0, 0.15))) +
    labs(title = "MS2", x = "# Missing", y = "Count")
  
  p1 / p2 + plot_layout(axes = "collect")
}


filter_too_many_missing <- function(MS1.df, MS2.df) {
  keep1 <- MS1.df %>% pivot_longer(!Gene.Name, names_to = "Sample.Name", ) %>%
    mutate(value = !is.na(value), Sample.Name = str_remove(Sample.Name, "_\\d+$")) %>%
    group_by(Gene.Name, Sample.Name) %>%
    summarise(Group.Count = sum(value), .groups = "drop") %>%
    ungroup() %>%
    pivot_wider(names_from = Sample.Name, values_from = Group.Count) %>%
    filter(if_any(everything(), ~ . == 3)) %>%
    pull(Gene.Name)
  
  keep2 <- MS2.df %>% pivot_longer(!Gene.Name, names_to = "Sample.Name", ) %>%
    mutate(value = !is.na(value), Sample.Name = str_remove(Sample.Name, "_\\d+$")) %>%
    group_by(Gene.Name, Sample.Name) %>%
    summarise(Group.Count = sum(value), .groups = "drop") %>%
    ungroup() %>%
    pivot_wider(names_from = Sample.Name, values_from = Group.Count) %>%
    filter(if_any(everything(), ~ . == 3)) %>%
    pull(Gene.Name)
  
  intersect(keep1, keep2)
}
  

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


plot_missing_heatmap <- function(MS1.df, MS2.df, mnar, colors_rows = NULL, colors_columns = NULL) {
  missing_values_MS1 <- MS1.df %>%
    column_to_rownames("Gene.Name") %>%
    mutate(across(everything(), ~ is.na(.))) %>%
    mutate(na_vals = rowSums(.))
    
    
  missing_values_MS2 <- MS2.df %>%
    column_to_rownames("Gene.Name") %>%
    mutate(across(everything(), ~ is.na(.))) %>%
    mutate(na_vals = rowSums(.))
  
  missing_values_hm <- merge(missing_values_MS1, missing_values_MS2, 
                             by = "row.names", suffixes = c("_MS1", "_MS2")) %>%
    filter() %>%
    mutate(across(!c(Gene.Name, na_vals), ~ is.na(.))) %>%
    select(!na_vals) %>%
    column_to_rownames("Gene.Name") %>%
    mutate(across(everything(), ~ as.integer(.x))) %>%
    as.matrix()
  
  MNAR.hm <- mnar[names(mnar) %in% rownames(missing_values_hm)]
  
  Heatmap(missing_values_hm, col = c("black", "gray"),  name = "Missing Value", 
          cluster_rows = TRUE, show_row_dend = FALSE, cluster_row_slices = FALSE,
          left_annotation = rowAnnotation(MNAR = MNAR.hm,
                                          col = colors_rows,
                                          show_annotation_name = FALSE),
          top_annotation = HeatmapAnnotation(Celltype = colData(data.filtered)$Cell.Type,
                                             Condition = colData(data.filtered)$Condition,
                                             col = colors_columns),
          split = MNAR.hm, show_column_names = FALSE, show_row_names = FALSE, 
          row_title = NULL, column_title = NULL, cluster_columns = FALSE, 
          heatmap_legend_param = list(labels = c("Yes", "No")))
}


plot_missing_density <- function(df) {
  missing_values_density <- df %>%
    mutate(na_vals = rowSums(is.na(.)) > 0) %>%
    pivot_longer(!c(Gene.Name, na_vals), names_to = "Sample.Name", values_to = "Intensity") %>%
    na.omit()
  
  p1 <- ggplot(missing_values_density, aes(x = Intensity, color = na_vals)) +
    geom_density() +
    labs(x = "Log2Intensity", y = "Density", color = "Has missing values") + 
    scale_color_discrete(labels = c("No", "Yes")) +
    theme_publish() + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(linewidth = 0.5))
  
  p2 <- ggplot(missing_values_density, aes(x = Intensity, color = na_vals)) +
    stat_ecdf() +
    labs(x = "Log2Intensity", y = "Cumulative Fraction", color = "Has missing values") + 
    scale_color_discrete(labels = c("No", "Yes")) +
    theme_publish() +
    theme(panel.grid.major = element_line(linewidth = 0.5))
  
  p1 / p2 + plot_layout(axes = "collect")
  
}


plot_imputation_density <- function(df_before, df_after, columndata) {
  samples <- columndata %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample.Name")
  
  missing_values_density <- df_before %>%
    mutate(na_vals = rowSums(is.na(.)) > 0) %>%
    pivot_longer(!c(Gene.Name, na_vals), names_to = "Sample.Name", values_to = "Intensity") %>%
    na.omit()
  
  imputed_plot <- df_after %>%
    pivot_longer(!Gene.Name, names_to = "Sample.Name", values_to = "Intensity") %>%
    left_join(samples, by = "Sample.Name", copy = TRUE)
  
  p1 <- ggplot(left_join(missing_values_density, samples, by = "Sample.Name"), aes(x = Intensity, color = interaction(Cell.Type, Condition, sep = ":"))) +
    geom_density() +
    labs(title = "Before", x = "Log2Intensity", y = "Density", color = "Celltype:Condition") + 
    theme_publish() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(linewidth = 0.5))
  
  p2 <- ggplot(imputed_plot, aes(x = Intensity, color = interaction(Cell.Type, Condition, sep = ":"))) +
    geom_density() +
    labs(title = "After", x = "Log2Intensity", y = "Density", color = "Celltype:Condition") + 
    theme_publish() +
    theme(panel.grid.major = element_line(linewidth = 0.5))
  
  p1 / p2
}
