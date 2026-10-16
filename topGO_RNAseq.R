# Gene ontology analysis using topGO - RNA-seq DEGs 
# weight01 / Fisher, up- and down-regulated genes tested separately
# Background = all expressed and GO-annotated genes in each contrast

library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(topGO)
library(ggplot2)
library(ggprism)

base_dir    <- "/Users/oliviaburleigh/Documents/SEQ_DATA"
deseq2_dir  <- file.path(base_dir, "DESEQ2/output_files/full_results_for_topGO")
output_root <- file.path(base_dir, "GO_Enrichment/topGO_RNAseq")

padj_cutoff <- 0.05
lfc_cutoff  <- 1
node_size   <- 5
top_nodes   <- 30
ontologies  <- c("BP", "MF", "CC")

contrasts <- c("CS_vs_CA", "HS6_vs_HA6", "HS18_vs_HA18",
               "HA6_vs_CA", "HA18_vs_CA", "HS6_vs_CS", "HS18_vs_CS",
               "Interaction_6h", "Interaction_18h")


###############################################################################
# GO-term mapping (AIPGENE ID -> LOC GENE ID -> GO terms)

map_df <- read_csv(file.path(base_dir, "Aiptasia_Genome/Aiptasia_AIPGENE_LOC_mapping_GCA_GCF.csv"),
                   show_col_types = FALSE) %>%
  transmute(AIPGENE = as.character(AIPGENE), LOC_id = as.character(LOC_ID)) %>%
  filter(!is.na(AIPGENE), AIPGENE != "", !is.na(LOC_id), LOC_id != "") %>%
  distinct(AIPGENE, .keep_all = TRUE)

gene2go_raw <- read_tsv(file.path(base_dir, "GO_Enrichment/Aiptasia_gene2go.tsv"),
                        col_names = c("LOC_id", "GO_terms"),
                        col_types = cols(.default = col_character())) %>%
  filter(!is.na(LOC_id), LOC_id != "", !is.na(GO_terms), GO_terms != "")

geneID2GO <- setNames(strsplit(gene2go_raw$GO_terms, ";"), gene2go_raw$LOC_id)


###############################################################################
# Functions

run_topgo <- function(geneList, ontology) {
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = node_size)
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  GenTable(GOdata, weightFisher = result, orderBy = "weightFisher",
           topNodes = min(top_nodes, length(score(result)))) %>%
    # GenTable reports very small p-values as "< 1e-30"; strip "<" so they aren't lost
    mutate(weightFisher = as.numeric(sub("<", "", weightFisher)))
}

# Bar plot of significant terms (distinct GO.IDs + make.unique stop terms with
# identical wrapped labels from stacking into one bar)
plot_topgo <- function(tab, title) {
  go_sig <- tab %>%
    filter(weightFisher < 0.05) %>%
    distinct(GO.ID, .keep_all = TRUE) %>%
    mutate(GeneRatio  = Significant / Annotated,
           neglog10_p = -log10(weightFisher)) %>%
    arrange(weightFisher) %>%
    slice_head(n = top_nodes) %>%
    mutate(label = make.unique(str_wrap(Term, width = 50)),
           label = factor(label, levels = rev(label)))

  if (nrow(go_sig) == 0) return(NULL)

  ggplot(go_sig, aes(x = GeneRatio, y = label, fill = neglog10_p)) +
    geom_col() +
    scale_fill_gradient(low = "grey80", high = "blue", name = "-log10(p)") +
    labs(x = "GeneRatio (Significant / Annotated)", y = NULL, title = title) +
    theme_prism() +
    theme(axis.text.y  = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"),
          plot.title   = element_text(size = 14, face = "bold"))
}


###############################################################################
# Run topGO for each contrast 

all_results <- list()

for (comp in contrasts) {
  outdir <- file.path(output_root, comp)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  deg <- read_csv(file.path(deseq2_dir, paste0(comp, "_ALL_results.csv")), show_col_types = FALSE) %>%
    left_join(map_df, by = c("gene_id" = "AIPGENE")) %>%
    filter(!is.na(LOC_id)) %>%
    arrange(padj) %>%
    distinct(LOC_id, .keep_all = TRUE)

  background <- intersect(deg$LOC_id, names(geneID2GO))
  sig <- deg %>% filter(!is.na(padj), padj < padj_cutoff, !is.na(log2FoldChange))
  foreground <- list(up   = intersect(sig$LOC_id[sig$log2FoldChange >  lfc_cutoff], background),
                     down = intersect(sig$LOC_id[sig$log2FoldChange < -lfc_cutoff], background))

  for (direction in names(foreground)) {
    geneList <- factor(as.integer(background %in% foreground[[direction]]))
    names(geneList) <- background

    for (ont in ontologies) {
      tab <- run_topgo(geneList, ont) %>%
        mutate(ontology = ont, comparison = comp, direction = direction,
               nUniverse = length(background), nSelected = length(foreground[[direction]]))

      prefix <- paste0("topGO_", comp, "_", ont, "_", direction)
      write_tsv(tab, file.path(outdir, paste0(prefix, "_weight01_fisher.tsv")))

      p <- plot_topgo(tab, paste(comp, ont, direction, sep = " | "))
      if (!is.null(p)) ggsave(file.path(outdir, paste0(prefix, "_colorplot_GeneRatio.pdf")), p, width = 9, height = 9)

      all_results[[prefix]] <- tab
    }
  }
}

# All contrasts & ontologies & DEG directions in one table 
write_tsv(bind_rows(all_results), file.path(output_root, "topGO_ALL_contrasts_ALL_ontologies_split.tsv"))
