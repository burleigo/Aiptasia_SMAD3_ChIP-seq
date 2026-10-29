# SMAD3 ChIP-seq and RNA-seq integration — Pearson correlation between SMAD3 binding and target gene expression
# For each putative SMAD3 target gene, Pearson's r is calculated between the VST SMAD3 binding at its representative peak and its VST expression across the 18 paired samples.
# The peak-gene pair and AIPGENE gene_id for each target come from Peak_to_gene_correlation.R (ALL_contrasts_target_genes.csv)

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggnewscale)
library(ggh4x)

data_dir      <- "/Users/oliviaburleigh/Documents/SEQ_DATA"
targets_file  <- file.path(data_dir, "GO_Enrichment/topGO_PeakToGene/topGO_ds_promoterOverlap/Intersection_gene_lists/ALL_contrasts_target_genes.csv")
rna_vst_file  <- file.path(data_dir, "DESEQ2/output_files/normalized_counts_vst_interaction.csv")
chip_vst_file <- file.path(data_dir, "DiffBind/Diffbind_ds/CSV_files/chip_vst_interaction.csv")
beta_dir      <- file.path(data_dir, "BETA/output_c0.05")
plot_dir        <- file.path(data_dir, "ChIP-RNA_Correlation/Integration_plots")

output_dir <- file.path(data_dir, "ChIP-RNA_Correlation/CSV_files")
beta_output_dir <- file.path(data_dir, "BETA/beta_validation")

# Contrasts of interest for this integration of datasets
Between_symState_contrasts <- c("CS_vs_CA", "HS6_vs_HA6", "HS18_vs_HA18")

cond_levels <- c("CA", "HA6", "HA18", "CS", "HS6", "HS18")

###############################################################################
# VARIANCE-STABILIZED COUNTS (non-blind)
# RNA:  DESeq2 interaction model, rows = AIPGENE gene_id (saved in DESeq2_RNAseq.R)
# ChIP: DiffBind/DESeq2 interaction model, row i = peak_id i (saved in DiffBind script as chip_vst_interaction.csv)

rna_vst  <- as.matrix(read.csv(rna_vst_file,  row.names = 1, check.names = FALSE))
chip_vst <- as.matrix(read.csv(chip_vst_file, row.names = 1, check.names = FALSE))

###############################################################################
# SAMPLE PAIRING
# Match each ChIP library to its RNA library by condition + replicate

parse_sample <- function(ids) {
  clean <- toupper(gsub("[^A-Za-z0-9]", "", sub("_RNA$", "", ids)))
  cond  <- str_extract(clean, "^(HA18|HS18|HA6|HS6|CA|CS)")
  tibble(id = ids, condition = cond, rep = str_extract(str_remove(clean, cond), "[0-9]+$"))
}

pair_map <- inner_join(parse_sample(colnames(chip_vst)) %>% dplyr::rename(chip_id = id),
                       parse_sample(colnames(rna_vst))  %>% dplyr::rename(rna_id = id),
                       by = c("condition", "rep")) %>%
  mutate(condition = factor(condition, levels = cond_levels)) %>%
  arrange(condition, rep)

nrow(pair_map) # 18 paired samples (6 conditions x 3 replicates)

###############################################################################
# PEARSON CORRELATION PER TARGET GENE (one row per gene per contrast)

targets <- read_csv(targets_file, show_col_types = FALSE) %>%
  filter(contrast %in% Between_symState_contrasts)

nrow(targets) # 666 gene-contrast rows (625 unique genes)

calc_cor <- function(gene_id, peak_id, method) {
  cor(chip_vst[peak_id, pair_map$chip_id], rna_vst[gene_id, pair_map$rna_id], method = method)
}

targets <- targets %>%
  mutate(pearson_r = unname(mapply(calc_cor, gene_id, peak_id, "pearson"))) %>%
  arrange(match(contrast, Between_symState_contrasts), desc(abs(pearson_r)))

write_csv(targets, file.path(output_dir, "Between_symState_targets_Pearson.csv"))


###############################################################################
# BETA TARGET PREDICTIONS
# One file per contrast and direction ({contrast}_uptarget.txt, {contrast}_downtarget.txt); refseqID = AIPGENE gene_id
# If a gene is predicted in both directions, the prediction with the lower rank product is kept

beta <- expand_grid(contrast = Between_symState_contrasts, beta_direction = c("Up", "Down")) %>%
  mutate(path = file.path(beta_dir, paste0(contrast, "_", tolower(beta_direction), "target.txt"))) %>%
  pmap_dfr(function(contrast, beta_direction, path) {
    read_tsv(path, col_types = cols(.default = "c")) %>%
      transmute(gene_id            = str_trim(refseqID),
                contrast, beta_direction,
                beta_rank_product  = as.numeric(rank_product),
                beta_reg_potential = as.numeric(reg_potential))
  }) %>%
  group_by(gene_id, contrast) %>%
  slice_min(beta_rank_product, n = 1, with_ties = FALSE) %>%
  ungroup()

###############################################################################
# FILE S3: PUTATIVE SMAD3 TARGET GENES + BETA PREDICTIONS

Pearson_plus_BETA <- targets %>%
  left_join(beta, by = c("gene_id", "contrast")) %>%
  transmute(`Gene ID`                   = gene_id,
            `LOC ID`                    = LOC_id,
            `Gene description`          = gene_desc,
            `Contrast`                  = str_replace_all(contrast, "_", " "),
            `SMAD3 regulatory mode`     = smad3_mode,
            `ChIP log2FC`               = round(chip_log2FC, 2),
            `RNA log2FC`                = round(rna_log2FC, 2),
            `RNA padj`                  = signif(rna_padj, 3),
            `Pearson r`                 = round(pearson_r, 3),
            `Peak ID`                   = peak_id,
            `Peak distance to TSS (bp)` = dist_to_tss,
            `BETA target`               = if_else(!is.na(beta_direction), "Yes", "No"),
            `BETA regulatory potential` = beta_reg_potential,
            `BETA rank product`         = signif(beta_rank_product, 3))

write_csv(Pearson_plus_BETA, file.path(beta_output_dir, "File_S3_Putative_SMAD3_targeted_genes.csv"), na = "")

###############################################################################
# SUMMARY: number of target genes also predicted as targets by BETA

BETA_summary <- Pearson_plus_BETA %>%
  group_by(Contrast, `SMAD3 regulatory mode`) %>%
  summarise(n_genes     = n(),
            n_in_beta   = sum(`BETA target` == "Yes"),
            pct_in_beta = round(100 * n_in_beta / n_genes, 1),
            .groups = "drop")

print(BETA_summary)
  # Contrast     `SMAD3 regulatory mode` n_genes n_in_beta pct_in_beta
  # 1 CS vs CA     Activation                   95        85        89.5
  # 2 CS vs CA     Repression                   81        69        85.2
  # 3 HS18 vs HA18 Activation                  144       131        91  
  # 4 HS18 vs HA18 Repression                  302       253        83.8
  # 5 HS6 vs HA6   Activation                   16        10        62.5
  # 6 HS6 vs HA6   Repression                   28        19        67.9

write_csv(BETA_summary, file.path(beta_output_dir, "BETA_summary_by_contrast_and_mode.csv"))


###############################################################################
# SCATTER PLOTS: SMAD3 BINDING VS EXPRESSION FOR SELECTED TARGET GENES 
# Genes of interest are selected from File S3 by Gene ID and contrast
# peak, gene ID, regulatory mode and Pearson r all come this file

cond_palette <- c(CA = "#B2BEB5", HA6 = "#71797E", HA18 = "#36454F",
                  CS = "#E49B0F", HS6 = "#FF7F50", HS18 = "#C04000")

# One color per functional category
category_palette <- c("Mucociliary surface"                = "#76D7C4",
                      "Carbon/nitrogen exchange"           = "#7FB3D2",
                      "Retinoid nuclear receptors"         = "#CCFF99",
                      "Aquaporins"                         = "#5499E9",
                      "Immunity"                           = "pink1",
                      "Sterol/lipid trafficking"           = "#BB8FCE",
                      "TGF-\U03B2/SMAD pathway regulation" = "coral",
                      "GPCR signaling"                     = "#F7DC6F",
                      "RNA modification"                   = "#F5B041",
                      "Chromatin regulation"               = "#AAB7B8")

# Selected target genes: File S3 Gene ID + contrast, short gene name for the panel label, and functional category
selected_genes <- tribble(
  ~gene_id,        ~Contrast,      ~gene_name, ~category,
  # CS vs CA, activation (Fig 6C)
  "AIPGENE28493",  "CS vs CA",     "HNF4G",    "Retinoid nuclear receptors",
  "AIPGENE19169",  "CS vs CA",     "GBP1",     "Immunity",
  "AIPGENE26763",  "CS vs CA",     "GS",       "Carbon/nitrogen exchange",
  "AIPGENE3689",   "CS vs CA",     "RXRA",     "Retinoid nuclear receptors",
  "AIPGENE8831",   "CS vs CA",     "Mucin",    "Mucociliary surface",
  # HS18 vs HA18, activation (Fig 6D)
  "AIPGENE27506",  "HS18 vs HA18", "PPP3CA",   "Immunity",
  "AIPGENE17680",  "HS18 vs HA18", "STRA6",    "Retinoid nuclear receptors",
  "AIPGENE17420",  "HS18 vs HA18", "AMT1",     "Carbon/nitrogen exchange",
  "AIPGENE15748",  "HS18 vs HA18", "AQP5",     "Aquaporins",
  "AIPGENE7917",   "HS18 vs HA18", "AQP3",     "Aquaporins",
  # CS vs CA, repression (Fig 7C)
  "AIPGENE3167",   "CS vs CA",     "NPC2",     "Sterol/lipid trafficking",
  "AIPGENE11655",  "CS vs CA",     "HTR1A",    "GPCR signaling",
  "AIPGENE18931",  "CS vs CA",     "CD109",    "Immunity",
  "AIPGENE18788",  "CS vs CA",     "GLTP",     "Sterol/lipid trafficking",
  "AIPGENE2250",   "CS vs CA",     "ZCCHC4",   "RNA modification",
  # HS18 vs HA18, repression (Fig 7D)
  "AIPGENE10787",  "HS18 vs HA18", "ADGRD1",   "GPCR signaling",
  "AIPGENE18720",  "HS18 vs HA18", "HRH2",     "GPCR signaling",
  "AIPGENE15844",  "HS18 vs HA18", "TACR1",    "GPCR signaling",
  "AIPGENE10187",  "HS18 vs HA18", "C2",       "Immunity",
  "AIPGENE19496",  "HS18 vs HA18", "GSK3",     "TGF-\U03B2/SMAD pathway regulation",
  "AIPGENE10266",  "HS18 vs HA18", "KAT6A",    "Chromatin regulation")

selected_genes <- selected_genes %>%
  inner_join(Pearson_plus_BETA, by = c("gene_id" = "Gene ID", "Contrast")) %>%
  mutate(panel_label = paste0(gene_name, "\n", gene_id),
         r_label     = sprintf("r = %.2f", `Pearson r`))

nrow(selected_genes) # 21

# Binding and expression values for each selected gene across the 18 paired samples
plot_data <- selected_genes %>%
  rowwise() %>%
  reframe(panel_label, Contrast, `SMAD3 regulatory mode`,
          condition = pair_map$condition,
          chip      = chip_vst[`Peak ID`, pair_map$chip_id],
          rna       = rna_vst[gene_id, pair_map$rna_id])

# One figure per contrast and regulatory mode; panels ordered by |r|
make_scatter <- function(comp, mode, plot_title) {
  genes <- selected_genes %>%
    filter(Contrast == comp, `SMAD3 regulatory mode` == mode) %>%
    arrange(desc(abs(`Pearson r`))) %>%
    mutate(panel_label = factor(panel_label, levels = panel_label))
  
  df <- plot_data %>%
    filter(Contrast == comp, `SMAD3 regulatory mode` == mode) %>%
    mutate(panel_label = factor(panel_label, levels = levels(genes$panel_label)))
  
  # r goes in the empty corner: top-left for activation, top-right for repression
  r_x     <- if (mode == "Activation") -Inf else Inf
  r_hjust <- if (mode == "Activation") -0.1 else 1.1
  
  ggplot(df, aes(x = chip, y = rna)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey50", linewidth = 0.8, linetype = "dashed") +
    geom_point(aes(fill = condition), shape = 21, size = 4, color = "white", stroke = 0.4) +
    geom_text(data = genes, aes(label = r_label), x = r_x, y = Inf,
              hjust = r_hjust, vjust = 1.5, size = 5, inherit.aes = FALSE) +
    scale_fill_manual(values = cond_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +   # empty band at the top of each panel for the r label
    facet_wrap2(~ panel_label, scales = "free", ncol = 2,
                strip = strip_themed(background_x = elem_list_rect(fill = unname(category_palette[genes$category])))) +
    labs(title = plot_title, x = "ChIP binding (VST)", y = "Expression (VST)") +
    theme_bw(base_size = 18) +
    theme(panel.grid      = element_blank(),
          strip.text      = element_text(face = "bold", size = 14, lineheight = 0.9),
          axis.text       = element_text(size = 16),
          axis.title      = element_text(size = 21),
          plot.title      = element_text(face = "bold", size = 16),
          legend.position = "none")
}

# 4 panels 
scatter_figs <- tribble(
  ~comp,          ~mode,        ~plot_title,                          ~file_name,
  "CS vs CA",     "Activation", "CS vs CA (upregulated genes)",       "CS_vs_CA_activation_scatter.pdf",
  "HS18 vs HA18", "Activation", "HS18 vs HA18 (upregulated genes)",   "HS18_vs_HA18_activation_scatter.pdf",
  "CS vs CA",     "Repression", "CS vs CA (downregulated genes)",     "CS_vs_CA_repression_scatter.pdf",
  "HS18 vs HA18", "Repression", "HS18 vs HA18 (downregulated genes)", "HS18_vs_HA18_repression_scatter.pdf")

# Make graphs
for (i in seq_len(nrow(scatter_figs))) {
  fig     <- scatter_figs[i, ]
  p       <- make_scatter(fig$comp, fig$mode, fig$plot_title)
  n_genes <- sum(selected_genes$Contrast == fig$comp & selected_genes$`SMAD3 regulatory mode` == fig$mode)
  ggsave(file.path(plot_dir, fig$file_name), p, width = 6.5, height = 2.8 * ceiling(n_genes / 2) + 1.5)
}

###############################################################################
# Generate legend separately (conditions and functional category) 

legend_plot <- ggplot() +
  geom_point(data = tibble(condition = factor(names(cond_palette), levels = names(cond_palette))),
             aes(x = 1, y = 1, fill = condition), shape = 21, size = 5, color = "white") +
  scale_fill_manual(values = cond_palette, name = "Condition") +
  new_scale_fill() +
  geom_point(data = tibble(category = factor(names(category_palette), levels = names(category_palette))),
             aes(x = 1, y = 1, fill = category), shape = 22, size = 5) +
  scale_fill_manual(values = category_palette, name = "Functional category",
                    labels = function(x) str_wrap(x, 25)) +
  theme_void() +
  theme(legend.title      = element_text(face = "bold", size = 14),
        legend.text       = element_text(size = 12, lineheight = 0.8),
        legend.key.height = unit(1.75, "lines"))

legend_grob <- ggplotGrob(legend_plot)
legend_grob <- legend_grob$grobs[[which(legend_grob$layout$name %in% c("guide-box-right", "guide-box"))[1]]]

quartz(type = "pdf", file = file.path(plot_dir, "scatter_legend.pdf"), width = 3.5, height = 7)
grid::grid.draw(legend_grob)
dev.off()





