# RNA-Seq — Differential gene expression analysis using DESeq2

library(tximport)
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(vegan)

base_dir    <- "/Users/oliviaburleigh/Documents/SEQ_DATA/DESEQ2"
output_dir  <- file.path(base_dir, "output_files")
deg_root    <- file.path(output_dir, "pairwise_deg_lists")
topgo_root  <- file.path(output_dir, "full_results_for_topGO")
summary_dir <- file.path(output_dir, "summaries")
pca_dir     <- file.path(output_dir, "pca_plots")
volcano_dir <- file.path(output_dir, "volcano_plots")

for (d in c(output_dir, deg_root, topgo_root, summary_dir, pca_dir, volcano_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

padj_cutoff <- 0.05
lfc_cutoff  <- 1

###############################################################################
# Load data

# Configuration file 
samples <- read.csv(file.path(base_dir, "config_RNA_DESeq2.csv"), stringsAsFactors = FALSE)
samples$Condition <- factor(samples$Condition)
samples$Tissue    <- factor(samples$Tissue)
samples$Treatment <- factor(samples$Treatment, levels = c("Control", "Heat_6", "Heat_18"))

# Abundance files from Kallisto pseudoalignment
abundance_files <- file.path(base_dir, "host_abundance_tsv", paste0("host_", samples$SampleID, "_abundance.tsv"))
names(abundance_files) <- samples$SampleID
txi <- tximport(abundance_files, type = "kallisto", txOut = TRUE)

# Annotation file
annots <- read.csv(file.path(base_dir, "aip_tabulated_annots.csv"), stringsAsFactors = FALSE) %>%
  dplyr::select(Query, Hit.description)

###############################################################################
# Functions for saving results

save_all <- function(res, comparison) {
  df <- as.data.frame(res) %>% rownames_to_column("gene_id")
  write.csv(df, file.path(topgo_root, paste0(comparison, "_ALL_results.csv")), row.names = FALSE)

  out <- file.path(deg_root, comparison)
  dir.create(out, showWarnings = FALSE)
  df   <- left_join(df, annots, by = c("gene_id" = "Query"))
  sig  <- df %>% filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)
  up   <- sig %>% filter(log2FoldChange > 0)
  down <- sig %>% filter(log2FoldChange < 0)
  write.csv(sig,  file.path(out, paste0(comparison, "_SIG_all.csv")),  row.names = FALSE)
  write.csv(up,   file.path(out, paste0(comparison, "_SIG_up.csv")),   row.names = FALSE)
  write.csv(down, file.path(out, paste0(comparison, "_SIG_down.csv")), row.names = FALSE)

  data.frame(comparison = comparison, n_sig = nrow(sig), n_up = nrow(up), n_down = nrow(down))
}

###############################################################################
# PCA (blind VST) - Fig. S5A

dds_null <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)
dds_null <- dds_null[rowSums(counts(dds_null)) >= 10, ]
vsd_null <- vst(dds_null, blind = TRUE)

pcaData    <- plotPCA(vsd_null, intgroup = c("Tissue", "Treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$ColorGroup <- factor(paste(pcaData$Tissue, pcaData$Treatment, sep = "_"),
                             levels = c("Apo_Control", "Apo_Heat_6", "Apo_Heat_18",
                                        "Sym_Control", "Sym_Heat_6", "Sym_Heat_18"))
pcaData$HeatStatus <- factor(ifelse(pcaData$Treatment == "Control", "Control", "Heat"),
                             levels = c("Control", "Heat"))
pcaData$name <- sub("_RNA$", "", pcaData$name)

color_palette <- c("Apo_Control" = "#B2BEB5", "Apo_Heat_6" = "#71797E", "Apo_Heat_18" = "#36454F",
                   "Sym_Control" = "#E49B0F", "Sym_Heat_6" = "#FF7F50", "Sym_Heat_18" = "#C04000")
pretty_labels <- c("Apo_Control" = "CA", "Apo_Heat_6" = "HA6", "Apo_Heat_18" = "HA18",
                   "Sym_Control" = "CS", "Sym_Heat_6" = "HS6", "Sym_Heat_18" = "HS18")

pca <- ggplot(pcaData, aes(PC1, PC2, color = ColorGroup, shape = HeatStatus)) +
  geom_point(size = 5) +
  scale_color_manual(values = color_palette, labels = pretty_labels) +
  scale_shape_manual(values = c("Control" = 16, "Heat" = 17)) +
  geom_text_repel(aes(label = name), size = 4, max.overlaps = Inf,
                  box.padding = 0.4, point.padding = 0.5, force = 2, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        plot.title   = element_text(size = 20),
        axis.title   = element_text(size = 20),
        axis.text    = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 18)) +
  labs(title = "PCA: RNA-seq variance-stabilized expression", color = "Condition", shape = "Treatment")

print(pca)
ggsave(file.path(pca_dir, "PCA_RNA_blind_square_byShape.pdf"), pca, width = 7, height = 7)

###############################################################################
# PERMANOVA + BETADISPER (on blind VST distances)

d <- dist(t(assay(vsd_null)))

set.seed(42)
adonis2(d ~ Tissue * Treatment, data = samples, permutations = 999, by = "terms")
    #                 Df SumOfSqs      R2       F Pr(>F)
    #Tissue            1   136883 0.45488 20.5435  0.001 ***
    #Treatment         2    61267 0.20360  4.5975  0.004 **
    #Tissue:Treatment  2    22815 0.07582  1.7120  0.124
    #Residual         12    79957 0.26571
    #Total            17   300922 1.00000

permutest(betadisper(d, samples$Tissue),    permutations = 999) #p = 0.187
permutest(betadisper(d, samples$Treatment), permutations = 999) #p = 0.212

###############################################################################
# DESeq2 Interaction model (~Tissue * Treatment)

dds_int <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Tissue * Treatment)
dds_int <- dds_int[rowSums(counts(dds_int)) >= 10, ]
dds_int <- DESeq(dds_int)

resultsNames(dds_int)
    #"Intercept"  "Tissue_Sym_vs_Apo"  "Treatment_Heat_18_vs_Control"
    #"Treatment_Heat_6_vs_Control"  "TissueSym.TreatmentHeat_18"  "TissueSym.TreatmentHeat_6"

# Model-fit VST, used for ChIP/RNA correlation plots
vsd_int <- vst(dds_int, blind = FALSE)
write.csv(as.data.frame(assay(vsd_int)), file.path(output_dir, "normalized_counts_vst_interaction.csv"))

# Between symbiotic states
res_CS_vs_CA     <- results(dds_int, name = "Tissue_Sym_vs_Apo", alpha = padj_cutoff)
res_HS6_vs_HA6   <- results(dds_int, contrast = list(c("Tissue_Sym_vs_Apo", "TissueSym.TreatmentHeat_6")),  alpha = padj_cutoff)
res_HS18_vs_HA18 <- results(dds_int, contrast = list(c("Tissue_Sym_vs_Apo", "TissueSym.TreatmentHeat_18")), alpha = padj_cutoff)

# Heat contrasts within aposymbiotic animals
res_HA6_vs_CA   <- results(dds_int, name = "Treatment_Heat_6_vs_Control",  alpha = padj_cutoff)
res_HA18_vs_CA  <- results(dds_int, name = "Treatment_Heat_18_vs_Control", alpha = padj_cutoff)
res_HA6_vs_HA18 <- results(dds_int, contrast = list("Treatment_Heat_6_vs_Control", "Treatment_Heat_18_vs_Control"),
                           alpha = padj_cutoff)

# Heat contrasts within symbiotic animals
res_HS6_vs_CS   <- results(dds_int, contrast = list(c("Treatment_Heat_6_vs_Control",  "TissueSym.TreatmentHeat_6")),  alpha = padj_cutoff)
res_HS18_vs_CS  <- results(dds_int, contrast = list(c("Treatment_Heat_18_vs_Control", "TissueSym.TreatmentHeat_18")), alpha = padj_cutoff)
res_HS6_vs_HS18 <- results(dds_int, contrast = list(c("Treatment_Heat_6_vs_Control",  "TissueSym.TreatmentHeat_6"),
                                                    c("Treatment_Heat_18_vs_Control", "TissueSym.TreatmentHeat_18")),
                           alpha = padj_cutoff)

# Interaction terms
res_Interaction_6h  <- results(dds_int, name = "TissueSym.TreatmentHeat_6",  alpha = padj_cutoff)
res_Interaction_18h <- results(dds_int, name = "TissueSym.TreatmentHeat_18", alpha = padj_cutoff)

# Save results + DEG counts (Table S2A)
int_summary <- rbind(
  save_all(res_CS_vs_CA,        "CS_vs_CA"),
  save_all(res_HS6_vs_HA6,      "HS6_vs_HA6"),
  save_all(res_HS18_vs_HA18,    "HS18_vs_HA18"),
  save_all(res_HA6_vs_CA,       "HA6_vs_CA"),
  save_all(res_HA18_vs_CA,      "HA18_vs_CA"),
  save_all(res_HA6_vs_HA18,     "HA6_vs_HA18"),
  save_all(res_HS6_vs_CS,       "HS6_vs_CS"),
  save_all(res_HS18_vs_CS,      "HS18_vs_CS"),
  save_all(res_HS6_vs_HS18,     "HS6_vs_HS18"),
  save_all(res_Interaction_6h,  "Interaction_6h"),
  save_all(res_Interaction_18h, "Interaction_18h"))

print(int_summary)
write.csv(int_summary, file.path(summary_dir, "DEG_int_summary_counts.csv"), row.names = FALSE)

###############################################################################
# BETA input file format (gene_id, log2FC, padj) for symbiotic state contrasts

for (comp in c("CS_vs_CA", "HS6_vs_HA6", "HS18_vs_HA18")) {
  df <- read.csv(file.path(topgo_root, paste0(comp, "_ALL_results.csv")))
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), c("gene_id", "log2FoldChange", "padj")]
  write.table(df, paste0("RNA_", comp, "_beta.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

###############################################################################
# Volcano plots - Figs. S5B–H

make_volcano <- function(res, title) {
  df <- as.data.frame(res) %>%
    filter(!is.na(padj), padj > 0) %>%
    mutate(neglog10 = -log10(padj),
           sig = case_when(padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
                           padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
                           TRUE ~ "NS"))

  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = title, x = "log2 Fold Change", y = "-log10 adjusted p-value", color = "DE status") +
    theme_bw(base_size = 16)
}

volcano_list <- list(
  "CS vs CA"     = res_CS_vs_CA,
  "HS6 vs HA6"   = res_HS6_vs_HA6,
  "HS18 vs HA18" = res_HS18_vs_HA18,
  "HA6 vs CA"    = res_HA6_vs_CA,
  "HA18 vs CA"   = res_HA18_vs_CA,
  "HS6 vs CS"    = res_HS6_vs_CS,
  "HS18 vs CS"   = res_HS18_vs_CS)

# Shared axis limits across panels
ymax <- max(sapply(volcano_list, function(r) max(-log10(r$padj[r$padj > 0]), na.rm = TRUE)))
xmax <- max(sapply(volcano_list, function(r) max(abs(r$log2FoldChange), na.rm = TRUE)))

v <- lapply(names(volcano_list), function(nm) make_volcano(volcano_list[[nm]], nm))

panel <- (v[[1]] | v[[2]] | v[[3]]) /
         (v[[4]] | v[[5]] | plot_spacer()) /
         (v[[6]] | v[[7]] | plot_spacer()) +
  plot_layout(guides = "collect") &
  coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax)) &
  theme(legend.position = "bottom")

ggsave(file.path(volcano_dir, "volcano_panel_FigS5.pdf"), panel, width = 16, height = 12)


###############################################################################
# Other DESeq2 models (not used in manuscript)

# Additive (~Tissue + Treatment)
dds_add <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Tissue + Treatment)
dds_add <- DESeq(dds_add[rowSums(counts(dds_add)) >= 10, ])
summary(results(dds_add, name = "Tissue_Sym_vs_Apo",             alpha = padj_cutoff)) #9012 up, 8563 down
summary(results(dds_add, name = "Treatment_Heat_6_vs_Control",   alpha = padj_cutoff)) #6086 up, 4578 down
summary(results(dds_add, name = "Treatment_Heat_18_vs_Control",  alpha = padj_cutoff)) #3659 up, 3052 down
summary(results(dds_add, contrast = c("Treatment", "Heat_6", "Heat_18"), alpha = padj_cutoff)) #1599 up, 1836 down

# LRT: interaction (full ~Tissue * Treatment vs reduced ~Tissue + Treatment)
dds_lrt <- DESeq(dds_int, test = "LRT", reduced = ~ Tissue + Treatment)
sum(results(dds_lrt)$padj < 0.05, na.rm = TRUE) #1523
