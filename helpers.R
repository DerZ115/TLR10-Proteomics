
libraries <- c("tidyverse", "magrittr", "readxl", "stringr","MSstats", 
               "SummarizedExperiment", "MsCoreUtils", "DEqMS", "qvalue", "ggrepel",
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
    labs(title = "MS1", x = "# Missing", y = "Count") +
    theme_publish()
  
  p2 <- ggplot(missing_values_MS2 %>% dplyr::count(na_vals), aes(x = na_vals, y = n)) + 
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust=-0.3) +
    scale_x_continuous(breaks = seq(0, max(c(missing_values_MS1$na_vals, 
                                             missing_values_MS2$na_vals)), 
                                    by = 1)) +
    scale_y_continuous(expand = expansion(mult=c(0, 0.15))) +
    labs(title = "MS2", x = "# Missing", y = "Count") +
    theme_publish()
  
  p1 / p2 + plot_layout(axis_titles = "collect")
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


plot_missing_density <- function(MS1.df, MS2.df) {
  missing_values_density_MS1 <- MS1.df %>%
    mutate(na_vals = rowSums(is.na(.)) > 0) %>%
    pivot_longer(!c(Gene.Name, na_vals), names_to = "Sample.Name", values_to = "Intensity") 
  
  missing_values_density_MS2 <- MS2.df %>%
    mutate(na_vals = rowSums(is.na(.)) > 0) %>%
    pivot_longer(!c(Gene.Name, na_vals), names_to = "Sample.Name", values_to = "Intensity") 
  
  missing_values_density <- bind_rows(list(MS1=missing_values_density_MS1,
                                           MS2=missing_values_density_MS2), 
                                      .id = "Source") %>%
    na.omit()
  
  p1 <- ggplot(missing_values_density, aes(x = Intensity, color = interaction(Source, na_vals, sep = ":"))) +
    geom_density(linewidth = 1) +
    labs(x = "Log2Intensity", y = "Density", color = "Has missing values") + 
    scale_color_discrete(labels = c("No", "Yes")) +
    theme_publish() + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(linewidth = 0.5))
  
  p2 <- ggplot(missing_values_density, aes(x = Intensity, color = interaction(Source, na_vals, sep = ":"))) +
    stat_ecdf(linewidth = 1) +
    labs(x = "Log2Intensity", y = "Cumulative Fraction", color = "Has missing values") + 
    scale_color_discrete(labels = c("MS1:No", "MS2:No", "MS1:Yes", "MS2:Yes")) +
    theme_publish() +
    theme(panel.grid.major = element_line(linewidth = 0.5))
  
  p1 / p2 + plot_layout(axes = "collect")
  
}


plot_imputation_density <- function(before, after) {
  samples <- colData(after) %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample.Name")
  
  before.plot <- assays(before) %>%
    lapply(as.data.frame) %>%
    lapply(rownames_to_column, var = "Gene.Name") %>% 
    bind_rows(.id = "MS") %>%
    mutate(na_vals = rowSums(is.na(across(where(is.numeric)))) > 0) %>%
    pivot_longer(!c(Gene.Name, MS, na_vals), names_to = "Sample.Name", values_to = "Intensity") %>%
    na.omit() %>%
    left_join(samples, by = "Sample.Name")
  
  after.plot <- assays(after) %>%
    lapply(as.data.frame) %>%
    lapply(rownames_to_column, var = "Gene.Name") %>%
    bind_rows(.id = "MS") %>%
    mutate(na_vals = rowSums(is.na(across(where(is.numeric)))) > 0) %>%
    pivot_longer(!c(Gene.Name, MS, na_vals), names_to = "Sample.Name", values_to = "Intensity") %>%
    left_join(samples, by = "Sample.Name")
  
  
  p1 <- ggplot(before.plot, aes(x = Intensity, 
                                color = interaction(Cell.Type, Condition, sep = ":"), 
                                linetype = MS)) +
    geom_density(linewidth = 1) +
    scale_linetype_manual(values = c(MS1 = 3, MS2 = 1)) +
    labs(title = "Before", x = "Log2Intensity", y = "Density", color = "Celltype:Condition") + 
    theme_publish() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(linewidth = 0.5))
  
  p2 <- ggplot(after.plot, aes(x = Intensity, 
                               color = interaction(Cell.Type, Condition, sep = ":"),
                               linetype = MS)) +
    geom_density(linewidth = 1) +
    scale_linetype_manual(values = c(MS1 = 3, MS2 = 1)) +
    labs(title = "After", x = "Log2Intensity", y = "Density", color = "Celltype:Condition") + 
    theme_publish() +
    theme(panel.grid.major = element_line(linewidth = 0.5))
  
  p1 / p2 + plot_layout(axis_titles = "collect")
}


plot_pca <- function(data) {
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
    prcomp() %$% x %>%
    as.data.frame() %>%
    rownames_to_column("Sample.Name") %>%
    left_join(samples, by = "Sample.Name")
  
  ggplot(data.pca, aes(x = PC1, y = PC2, color = interaction(Cell.Type, Condition, sep = ":"))) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = Sample.Name), size=3) +
    labs(color = "Celltype:Condition") +
    theme_publish() +
    theme(panel.grid.major = element_line(linewidth = 0.5))
}


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