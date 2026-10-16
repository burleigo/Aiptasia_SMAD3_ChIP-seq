# CHIPSEEKER — genomic annotation of MSPC consensus peaks (c = 2)
# Promoter region = TSS ± 3 kb 
 
library(txdbmaker)
library(ChIPseeker)
library(ggplot2)
library(ggrepel)
 
consensus_bed <- "/Users/oliviaburleigh/Documents/SEQ_DATA/MSPC/MSPC_output_ds/Consensus_bed"
output_dir    <- "/Users/oliviaburleigh/Documents/SEQ_DATA/ChIPseeker/Plots"
 
cond_order <- c("CA", "HA6", "HA18", "CS", "HS6", "HS18")
 
consensus_files <- file.path(consensus_bed, paste0(cond_order, "_consensus_peaks.bed"))
names(consensus_files) <- cond_order
all_consensus_peak_file <- file.path(consensus_bed, "All_consensus_peaks.bed")
 
txdb <- makeTxDbFromGFF("/Users/oliviaburleigh/Documents/SEQ_DATA/ChIPseeker/GCA_001417965.1_Aiptasia_genome_1.1_genomic.gtf",
                        format = "gtf")
 
 
###############################################################################
# Annotations per condition (not used in manuscript)
 
peakAnnoList_3kb <- lapply(consensus_files, annotatePeak, TxDb = txdb,
                           tssRegion = c(-3000, 3000), verbose = FALSE)
 
p_bar <- plotAnnoBar(peakAnnoList_3kb)
ggsave(file.path(output_dir, "Fig_AnnotBar_3kb.pdf"), p_bar, width = 10, height = 5)
 
p_tss <- plotDistToTSS(peakAnnoList_3kb)
ggsave(file.path(output_dir, "Fig_distToTSS_3kb.pdf"), p_tss, width = 10, height = 5)
 
 
###############################################################################
# Pie chart, all conditions combined (Fig. 5A)
 
peak_annots_all_3kb <- annotatePeak(all_consensus_peak_file, TxDb = txdb,
                                    tssRegion = c(-3000, 3000), verbose = FALSE)
 
# Remove 5' UTR slice since it's negligible (0.03%)

anno_df <- droplevels(subset(peak_annots_all_3kb@annoStat, Feature != "5' UTR"))
anno_df$Feature <- factor(anno_df$Feature, levels = anno_df$Feature)
 
# Label positions and text
csum <- rev(cumsum(rev(anno_df$Frequency)))
anno_df$pos   <- anno_df$Frequency / 2 + c(csum[-1], 0)
anno_df$label <- ifelse(anno_df$Frequency < 2,
                        sprintf("%.1f%%", anno_df$Frequency),
                        sprintf("%.0f%%", anno_df$Frequency))
 
anno_colors <- c(
  "Promoter (<=1kb)"   = "gold3",
  "Promoter (1-2kb)"   = "gold2",
  "Promoter (2-3kb)"   = "#F5D76E",
  "3' UTR"             = "#238B45",
  "1st Exon"           = "#08519C",
  "Other Exon"         = "#6BAED6",
  "1st Intron"         = "firebrick",
  "Other Intron"       = "red",
  "Downstream (<=300)" = "black",
  "Distal Intergenic"  = "#969696")

# slices >= 10% labeled inside; smaller ones labeled outside
thr <- 10   
 
p_pie <- ggplot(anno_df, aes(x = "", y = Frequency, fill = Feature)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = anno_colors) +
  coord_polar(theta = "y", clip = "off") +
  geom_text(data = subset(anno_df, Frequency >= thr),
            aes(x = 1.21, y = pos, label = label), size = 6.5) +
  geom_text_repel(data = subset(anno_df, Frequency < thr),
                  aes(x = 1.5, y = pos, label = label), size = 5.5, nudge_x = 0.35,
                  segment.size = 0.4, min.segment.length = 0) +
  guides(fill = guide_legend(title = NULL, nrow = 5, byrow = TRUE)) +
  theme_void() +
  theme(legend.position    = "bottom",
        legend.text        = element_text(size = 12),
        legend.key.size    = unit(0.7, "cm"),
        legend.box.spacing = unit(-25, "pt"),
        legend.margin      = margin(t = 1),
        plot.margin        = margin(0, 0, 10, 0))
 
print(p_pie)
ggsave(file.path(output_dir, "Fig_annoPie_allcond_3kb.pdf"), p_pie, width = 6, height = 6)
