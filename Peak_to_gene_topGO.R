# SMAD3 ChIP-seq and RNA-seq integration — GO enrichment of putative SMAD3 target genes
# Target genes (significant DiffBind peak in promoter-proximal region AND is a DEG) are split by direction and tested with topGO (weight01, Fisher) against all expressed genes with GO annotations
# ChIP_up_RNA_up   = putative activation targets
# ChIP_up_RNA_down = putative repression targets

library(readr)
library(dplyr)
library(stringr)
library(topGO)
library(GO.db)
library(ggplot2)
library(ggprism)
library(GenomicRanges)

data_dir    <- "/Users/oliviaburleigh/Documents/SEQ_DATA"
rnaseq_dir  <- file.path(data_dir, "DESEQ2/output_files/full_results_for_topGO")
gene2go_tsv <- file.path(data_dir, "GO_Enrichment/Aiptasia_gene2go.tsv")
conv_file   <- file.path(data_dir, "Aiptasia_Genome/LOC_to_AIPGENE_conversion.csv")

output_root <- file.path(data_dir, "GO_Enrichment/topGO_PeakToGene/topGO_ds_promoterOverlap")
tsv_dir     <- file.path(output_root, "TSV_files")
plot_dir    <- file.path(output_root, "GO_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

node_size  <- 5
top_nodes  <- 20

# only testing direction subsets with at least 50 target genes otherwise not enough statistical power for GO enrichment
min_genes  <- 50   

ontologies <- c("BP", "MF", "CC")

# Contrasts of interest for this integration of datasets
Between_symState_contrasts <- c("CS_vs_CA", "HS6_vs_HA6", "HS18_vs_HA18")

###############################################################################
# GO MAPPING (LOC_id -> GO terms)

gene2go <- read_tsv(gene2go_tsv, col_names = c("LOC_id", "GO_ID"), col_types = "cc") %>%
  filter(!is.na(LOC_id), !is.na(GO_ID)) %>%
  group_by(LOC_id) %>%
  summarise(GO_ID = paste(unique(GO_ID), collapse = ";"), .groups = "drop")

gene2go <- setNames(strsplit(gene2go$GO_ID, ";"), gene2go$LOC_id)

###############################################################################
# topGO BACKGROUND (AIPGENE -> LOC_id)
# All genes tested in DESeq2 that have GO annotations

conv <- read_csv(conv_file, show_col_types = FALSE) %>% distinct(AIPGENE, LOC_ID) %>% dplyr::rename(LOC_id = LOC_ID)
names(conv)

# All of the RNA-seq *_ALL_results.csv files have the same list of genes from the dds_int model (see RNA_DESeq2.R)
rna_background <- read_csv(file.path(rnaseq_dir, "CS_vs_CA_ALL_results.csv"), show_col_types = FALSE) %>%
  inner_join(conv, by = c("gene_id" = "AIPGENE")) %>%
  pull(LOC_id) %>%
  unique()

go_background <- intersect(rna_background, names(gene2go))
length(go_background)
  # RNA background for topGO: 12,893 expressed genes with GO annotations

###############################################################################
# MAKE topGO FUNCTION
# Will run topGO for each ontology (BP, MF, CC), contrast, and direction (up- or down-regulated)

run_topgo <- function(target_genes, ontology) {
  gene_list <- factor(as.integer(go_background %in% target_genes), levels = c(0, 1))
  names(gene_list) <- go_background
  
  go_data <- new("topGOdata", ontology = ontology, allGenes = gene_list,
                 annot = annFUN.gene2GO, gene2GO = gene2go, nodeSize = node_size)
  go_result <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
  
  GenTable(go_data, weightFisher = go_result, orderBy = "weightFisher",
           topNodes = min(top_nodes, length(score(go_result)))) %>%
    # GenTable reports very small p-values as "< 1e-30"; remove "<" so these terms are kept
    mutate(weightFisher = as.numeric(sub("<", "", weightFisher)))
}

###############################################################################
# GO PLOT FUNCTION
# Uses the full GO term name from GO.db (GenTable truncates long names)

manual_relabel <- c("GO:0043491" = "PI3K/protein kinase B signal transduction")

plot_topgo <- function(go_tbl, plot_title) {
  go_sig <- go_tbl %>%
    filter(weightFisher < 0.05) %>%
    mutate(full_term  = coalesce(AnnotationDbi::Term(GO.ID), Term),
           full_term  = ifelse(GO.ID %in% names(manual_relabel), manual_relabel[GO.ID], full_term),
           GeneRatio  = Significant / Annotated,
           neglog10_p = -log10(weightFisher),
           label      = str_wrap(full_term, width = 35)) %>%
    arrange(weightFisher)
  
  if (nrow(go_sig) == 0) return(NULL)
  
  ggplot(go_sig, aes(x = GeneRatio, y = factor(label, levels = rev(unique(label))), fill = neglog10_p)) +
    geom_col() +
    scale_fill_gradient(low = "grey80", high = "red", name = "-log10(p)") +
    labs(x = "Gene Ratio (Significant / Annotated)", y = NULL, title = str_wrap(plot_title, width = 40)) +
    theme_prism() +
    theme(axis.text.y  = element_text(size = 14, lineheight = 0.9),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 15),
          plot.title   = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text  = element_text(size = 15))
}

###############################################################################
# RUN topGO AND SAVE PLOTS/FILES
# Only CS_vs_CA and HS18_vs_HA18 have direction subsets with >= 50 genes (HS6_vs_HA6: 16 up, 28 down)

for (comp in Between_symState_contrasts) {
  outdir <- file.path(tsv_dir, comp)
  
  intersection_tbl <- read_tsv(file.path(outdir, paste0(comp, "_ChIP_RNA_intersection.tsv")), show_col_types = FALSE) %>%
    mutate(direction = paste0("ChIP_", if_else(chip_log2FC > 0, "up", "down"), "_RNA_", tolower(rna_direction)))
  
  for (direction in unique(intersection_tbl$direction)) {
    target_genes <- unique(intersection_tbl$LOC_id[intersection_tbl$direction == direction])
    if (length(target_genes) < min_genes) next
    message(comp, " | ", direction, " | n = ", length(target_genes))
    
    for (ont in ontologies) {
      go_tbl <- run_topgo(target_genes, ont) %>%
        mutate(contrast = comp, direction = direction, ontology = ont,
               n_background = length(go_background), n_target_genes = length(target_genes))
      
      file_prefix <- paste0(comp, "_", direction, "_topGO_", ont)
      write_tsv(go_tbl, file.path(outdir, paste0(file_prefix, ".tsv")))
      
      p <- plot_topgo(go_tbl, paste(comp, direction, ont, sep = " | "))
      if (!is.null(p)) {
        ggsave(file.path(plot_dir, paste0(file_prefix, ".pdf")),
               p, width = 7, height = max(4, 0.45 * sum(go_tbl$weightFisher < 0.05, na.rm = TRUE) + 2))
      }
    }
  }
}
