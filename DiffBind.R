# Differential binding scripts using DiffBind

########################## diffbind_counts_ds.R ##########################
#!/usr/bin/env Rscript

library(BiocParallel)
library(DiffBind)

# Using original BAM files for DiffBind

# Load Samplesheet
df <- read.csv("/nfs4/IB/Weis_Lab/burleigo/chip_project/chipseq/diffbind/DiffBind_samplesheet_hpc_ds.csv")

# Create DBA Object
dbObj <- dba(sampleSheet = df)

# DBA Count
# Summits = 200 creates uniform 401bp peaks centered on the binding site
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
register(MulticoreParam(workers = ncores, tasks = 1))  

dbObj_counts <- dba.count(dbObj, 
    minOverlap = 2,
    summits = 200, 
    bParallel = TRUE)

# Normalize by library size 
dbObj_counts <- dba.normalize(dbObj_counts, 
    method = DBA_DESEQ2,
    library = DBA_LIBSIZE_BACKGROUND, 
    background = TRUE)

# Save RDS
saveRDS(dbObj_counts, "dbObj_counts_norm_ds.rds")
print("Done with dba.count")
##############################################################################


################################  Diffbind.sh ################################

library(DiffBind)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(vegan)
library(ggVennDiagram)

setwd("/Users/oliviaburleigh/Documents/SEQ_DATA/DiffBind/Diffbind_ds")
outdir   <- "CSV_files"
plot_dir <- "Volcano_plots"
pca_dir  <- "PCA_plots"
for (d in c(outdir, plot_dir, pca_dir)) dir.create(d, showWarnings = FALSE)

###############################################################################
# Load DiffBind object (counted + normalized on HPC)

dbObj_counts <- readRDS("dbObj_counts_norm_ds.rds")
dba.show(dbObj_counts)

peaks_gr <- dba.peakset(dbObj_counts, bRetrieve = TRUE)
length(peaks_gr) #17,569 peaks in universe

###############################################################################
# Export all peaks used in DiffBind analysis
# BETA: BED with score = mean normalized count across samples
# MEME/HOMER: peaks ranked by mean affinity (all + top 1000)

peaks_df_all <- as.data.frame(peaks_gr)
count_cols   <- setdiff(colnames(peaks_df_all), c("seqnames", "start", "end", "width", "strand"))

peaks_df_all <- peaks_df_all %>%
  mutate(peak_id       = paste0("peak_", row_number()),
         mean_affinity = rowMeans(across(all_of(count_cols))),
         start0        = start - 1)  

# BETA input
write.table(data.frame(peaks_df_all$seqnames, peaks_df_all$start0, peaks_df_all$end,
                       peaks_df_all$peak_id, round(peaks_df_all$mean_affinity, 3), "."),
            "SMAD3_consensus_peaks_LibSize_norm_counts.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Motif analysis input
peaks_meme <- peaks_df_all %>%
  arrange(desc(mean_affinity)) %>%
  dplyr::select(seqnames, start0, end, peak_id, mean_affinity)

write.table(peaks_meme, file.path(outdir, "MEME_peaks_ranked_all.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(head(peaks_meme, 1000), file.path(outdir, "MEME_peaks_top1000.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


###############################################################################
# Differential binding with additive model (~Tissue + Treatment) - used in Fig. 4A

dbObj_add <- dba.contrast(dbObj_counts, design = "~Tissue + Treatment")
dbObj_add <- dba.analyze(dbObj_add,
                         bBlacklist = FALSE, bGreylist = FALSE,
                         bParallel = FALSE, method = DBA_DESEQ2)

dba.show(dbObj_add, bContrast = TRUE, th = 0.05)
#   Factor   Group Samples  Group2 Samples2 DB.DESeq2
#1  Tissue     Apo       9     Sym        9      9003
#2  Treatment Control      6 Heat_18        6         0
#3  Treatment Control      6  Heat_6        6       304
#4  Treatment  Heat_6      6 Heat_18        6        17

###############################################################################
# Differential binding with interaction model (~Tissue * Treatment) 

dbObj_analysis <- dba.contrast(dbObj_counts, design = "~Tissue * Treatment")
dbObj_analysis <- dba.contrast(dbObj_analysis, contrast = "TissueSym.TreatmentHeat_18")
dbObj_analysis <- dba.contrast(dbObj_analysis, contrast = "TissueSym.TreatmentHeat_6")
dbObj_analysis <- dba.analyze(dbObj_analysis,
                              bBlacklist = FALSE, bGreylist = FALSE,
                              bParallel = FALSE, method = DBA_DESEQ2)

dba.show(dbObj_analysis, bContrast = TRUE, th = 0.05)
    #   Factor                      Group Samples  Group2 Samples2 DB.DESeq2
    #1  Tissue                        Apo       9     Sym        9      1448
    #2  Treatment                 Control       6 Heat_18        6        51
    #3  Treatment                 Control       6  Heat_6        6       105
    #4  Treatment                  Heat_6       6 Heat_18        6       138
    #5  Name TissueSym.TreatmentHeat_18                                   41
    #6  Name  TissueSym.TreatmentHeat_6                                   97

dds      <- dbObj_analysis$DESeq2$DEdata
peaks_df <- as.data.frame(dba.peakset(dbObj_analysis, bRetrieve = TRUE))

resultsNames(dds)
    #"Intercept"  "Tissue_Sym_vs_Apo"  "Treatment_Heat_18_vs_Control"
    #"Treatment_Heat_6_vs_Control"  "TissueSym.TreatmentHeat_18"  "TissueSym.TreatmentHeat_6"

write.csv(peaks_df, file.path(outdir, "DiffBind_peaks_coordinates_ds.csv"), row.names = FALSE)

###############################################################################
# Define contrasts to extract from DESeqDataSet (DDS) object

alpha <- 0.05

# Between symbiotic states
res_CS_vs_CA     <- results(dds, name = "Tissue_Sym_vs_Apo", alpha = alpha)
res_HS6_vs_HA6   <- results(dds, contrast = list(c("Tissue_Sym_vs_Apo", "TissueSym.TreatmentHeat_6")), alpha = alpha)
res_HS18_vs_HA18 <- results(dds, contrast = list(c("Tissue_Sym_vs_Apo", "TissueSym.TreatmentHeat_18")), alpha = alpha)

# Heat contrasts within aposymbiotic animals
res_HA6_vs_CA   <- results(dds, name = "Treatment_Heat_6_vs_Control", alpha = alpha)
res_HA18_vs_CA  <- results(dds, name = "Treatment_Heat_18_vs_Control", alpha = alpha)
res_HA6_vs_HA18 <- results(dds, contrast = c("Treatment", "Heat_6", "Heat_18"), alpha = alpha)

# Heat contrasts within symbiotic animals
res_HS6_vs_CS   <- results(dds, contrast = list(c("Treatment_Heat_6_vs_Control",  "TissueSym.TreatmentHeat_6")), alpha = alpha)
res_HS18_vs_CS  <- results(dds, contrast = list(c("Treatment_Heat_18_vs_Control", "TissueSym.TreatmentHeat_18")), alpha = alpha)
res_HS6_vs_HS18 <- results(dds, contrast = list(c("Treatment_Heat_6_vs_Control",  "TissueSym.TreatmentHeat_6"),
                                                c("Treatment_Heat_18_vs_Control", "TissueSym.TreatmentHeat_18")), alpha = alpha)

# Interaction terms
res_interaction_6  <- results(dds, name = "TissueSym.TreatmentHeat_6",  alpha = alpha)
res_interaction_18 <- results(dds, name = "TissueSym.TreatmentHeat_18", alpha = alpha)


###############################################################################
# Export DiffBind results from interaction model (full table + significant peaks for each contrast)

write_report <- function(res, out_prefix, interaction = FALSE) {
  df <- cbind(chr = peaks_df$seqnames, peaks_df[, c("start", "end")], as.data.frame(res))
  labels <- if (interaction) c("Stronger_in_Sym", "Stronger_in_Apo") else c("Up", "Down")
  df$direction <- ifelse(df$log2FoldChange > 0, labels[1], labels[2])

  sig <- df[!is.na(df$padj) & df$padj < alpha, ]
  sig <- sig[order(sig$padj), ]

  write.csv(df,  file.path(outdir, paste0(out_prefix, "_FULL.csv")), row.names = FALSE)
  write.csv(sig, file.path(outdir, paste0(out_prefix, "_SIG_padj", alpha, ".csv")), row.names = FALSE)
  cat(out_prefix, ":", nrow(sig), "significant\n")
}

write_report(res_CS_vs_CA,       "CS_vs_CA")        #1448
write_report(res_HS6_vs_HA6,     "HS6_vs_HA6")      #280
write_report(res_HS18_vs_HA18,   "HS18_vs_HA18")    #5636
write_report(res_HA6_vs_CA,      "HA6_vs_CA")       #183
write_report(res_HA18_vs_CA,     "HA18_vs_CA")      #102
write_report(res_HA6_vs_HA18,    "HA6_vs_HA18")     #474
write_report(res_HS6_vs_CS,      "HS6_vs_CS")       #63
write_report(res_HS18_vs_CS,     "HS18_vs_CS")      #0
write_report(res_HS6_vs_HS18,    "HS6_vs_HS18")     #0
write_report(res_interaction_6,  "Interaction_6h",  interaction = TRUE) #97
write_report(res_interaction_18, "Interaction_18h", interaction = TRUE) #41

###############################################################################
# PCA (blind VST)

vsd_blind <- vst(dds, blind = TRUE, fitType = "local")
pcaData   <- plotPCA(vsd_blind, intgroup = c("Tissue", "Treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$ColorGroup <- factor(paste(pcaData$Tissue, pcaData$Treatment, sep = "_"),
                             levels = c("Apo_Control", "Apo_Heat_6", "Apo_Heat_18",
                                        "Sym_Control", "Sym_Heat_6", "Sym_Heat_18"))

pcaData$HeatStatus <- factor(ifelse(pcaData$Treatment == "Control", "Control", "Heat"),
                             levels = c("Control", "Heat"))

color_palette <- c("Apo_Control" = "#B2BEB5", "Apo_Heat_6" = "#71797E", "Apo_Heat_18" = "#36454F",
                   "Sym_Control" = "#E49B0F", "Sym_Heat_6" = "#FF7F50", "Sym_Heat_18" = "#C04000")

graph_labels <- c("Apo_Control" = "CA", "Apo_Heat_6" = "HA6", "Apo_Heat_18" = "HA18",
                   "Sym_Control" = "CS", "Sym_Heat_6" = "HS6", "Sym_Heat_18" = "HS18")

pca <- ggplot(pcaData, aes(PC1, PC2, color = ColorGroup, shape = HeatStatus)) +
  geom_point(size = 5) +
  scale_color_manual(values = color_palette, labels = graph_labels) +
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
  labs(title = "PCA: Smad3 ChIP-seq binding affinity", color = "Condition", shape = "Treatment")

print(pca)
ggsave(file.path(pca_dir, "PCA_ds_square_byShape.pdf"), pca, width = 7, height = 7)

###############################################################################
# PERMANOVA + BETADISPER (on blind VST distances)

d_ds <- dist(t(assay(vsd_blind)))

set.seed(42)
adonis2(d_ds ~ Tissue * Treatment, data = colData(dds), permutations = 999, by = "terms")
    #                 Df SumOfSqs      R2      F Pr(>F)
    #Tissue            1  1608992 0.22056 5.3194  0.007 **
    #Treatment         2  1356961 0.18601 2.2431  0.038 *
    #Tissue:Treatment  2   699375 0.09587 1.1561  0.286
    #Residual         12  3629709 0.49756
    #Total            17  7295037 1.00000

permutest(betadisper(d_ds, colData(dds)$Tissue),    permutations = 999) #p = 0.477
permutest(betadisper(d_ds, colData(dds)$Treatment), permutations = 999) #p = 0.691

###############################################################################
# Volcano plots
# Figure 4C-E: between symbiotic states
# Figure S4: heat contrasts within each symbiotic state

sig_levels <- c("Up", "Down", "NS")
sig_colors <- c(Up = "red", Down = "blue", NS = "grey70")

make_volcano <- function(res, title, show_total = TRUE) {
  df <- as.data.frame(res) %>%
    filter(!is.na(padj)) %>%
    mutate(neglog10 = -log10(padj),
           sig = case_when(padj < alpha & log2FoldChange > 0 ~ "Up",
                           padj < alpha & log2FoldChange < 0 ~ "Down",
                           TRUE ~ "NS"),
           sig = factor(sig, levels = sig_levels))

  n_up   <- sum(df$sig == "Up")
  n_down <- sum(df$sig == "Down")

  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(size = 2, alpha = 0.5) +
    scale_color_manual(name = "Differential\nbinding", values = sig_colors,
                       limits = sig_levels, drop = FALSE) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    annotate("text", x = Inf,  y = Inf, label = n_up,   color = "red",  hjust = 1.5,  vjust = 1.5, size = 5) +
    annotate("text", x = -Inf, y = Inf, label = n_down, color = "blue", hjust = -0.5, vjust = 1.5, size = 5) +
    labs(title = if (show_total) paste0(title, " (", n_up + n_down, ")") else title, x = "Log2 Fold Change", y = "-log10(FDR)") +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank())
}

# Figure 4C-E
main_panel <- (make_volcano(res_CS_vs_CA,     "CS vs CA",     show_total = FALSE) |
               make_volcano(res_HS6_vs_HA6,   "HS6 vs HA6",   show_total = FALSE) |
               make_volcano(res_HS18_vs_HA18, "HS18 vs HA18", show_total = FALSE)) +
  plot_layout(guides = "collect", axis_titles = "collect")

ggsave(file.path(plot_dir, "volcano_main_3across.pdf"), main_panel, width = 13, height = 5)

# Figure S4: apo (A–C) top row, sym (D–F) bottom row
supp_panel <- make_volcano(res_HA6_vs_CA,   "HA6 vs CA") +
              make_volcano(res_HA18_vs_CA,  "HA18 vs CA") +
              make_volcano(res_HA6_vs_HA18, "HA6 vs HA18") +
              make_volcano(res_HS6_vs_CS,   "HS6 vs CS") +
              make_volcano(res_HS18_vs_CS,  "HS18 vs CS") +
              make_volcano(res_HS6_vs_HS18, "HS6 vs HS18") +
  plot_layout(ncol = 3, guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(file.path(plot_dir, "volcano_supplement_FigS4.pdf"), supp_panel, width = 15, height = 10)


###############################################################################
# Additive model bubble plot (Fig. 4A) 

bubble_df <- data.frame(
  contrast = c("Sym vs Apo", "Control vs Heat 6h", "Control vs Heat 18h", "Heat 6h vs Heat 18h"),
  effect   = c("Symbiotic state", "Treatment", "Treatment", "Treatment"),
  n_peaks  = add_counts$DB.DESeq2[c(1, 3, 2, 4)])
bubble_df$contrast <- factor(bubble_df$contrast, levels = rev(bubble_df$contrast))

effect_colors <- c("Symbiotic state" = "#C04000", "Treatment" = "#36454F")

bubble <- ggplot(bubble_df, aes(x = 1, y = contrast, color = effect)) +
  geom_point(data = subset(bubble_df, n_peaks > 0), aes(size = n_peaks)) +
  geom_point(data = subset(bubble_df, n_peaks == 0), shape = 1, size = 3, stroke = 1) +
  geom_text(aes(label = format(n_peaks, big.mark = ",")), nudge_x = 0.5, hjust = 0, size = 5, show.legend = FALSE) +
  scale_color_manual(values = effect_colors, name = "Effect type") +
  scale_size_area(max_size = 20, guide = "none") +
  scale_x_continuous(limits = c(0.5, 2)) +
  labs(x = "Additive model\n(~ Symbiotic state + Treatment)", y = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())

print(bubble)
ggsave(file.path(plot_dir, "Fig4A_additive_bubbleplot.pdf"), bubble, width = 7, height = 4.5)

###############################################################################
# Venn of symbiotic state contrasts (interaction model) (Fig. 4B)

get_sig <- function(res) rownames(res)[!is.na(res$padj) & res$padj < alpha]

venn <- ggVennDiagram(list(
  "CS vs CA"     = get_sig(res_CS_vs_CA),
  "HS6 vs HA6"   = get_sig(res_HS6_vs_HA6),
  "HS18 vs HA18" = get_sig(res_HS18_vs_HA18)),
  label_alpha = 0, label = "count", label_size = 5, set_size = 5, edge_size = 0.5) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme(plot.margin = margin(10, 20, 10, 20))

# White count labels on dark regions
venn$layers[[4]]$data$label_col <- ifelse(as.numeric(venn$layers[[4]]$data$count) > 3000, "white", "black")
venn$layers[[4]]$mapping$colour <- aes(colour = label_col)$colour
venn$layers[[4]]$aes_params$colour <- NULL
venn <- venn + scale_color_identity()

print(venn)
ggsave(file.path(plot_dir, "Fig4B_Venn_tissue_effects.pdf"), venn, width = 5, height = 5)
         
