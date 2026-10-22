# SMAD3 ChIP-seq and RNA-seq integration — peak-to-gene assignment
# Differentially bound peaks (DiffBind, p-adj < 0.05) are assigned to every gene whose promoter window (TSS ± 3 kb) contains the peak midpoint, then crossed with DEGs (DESeq2, p-adj < 0.05, absVal LFC > 1) from the same contrast.

library(readr)
library(dplyr)
library(tidyr)
library(GenomicRanges)

data_dir      <- "/Users/oliviaburleigh/Documents/SEQ_DATA"
diffbind_dir  <- file.path(data_dir, "DiffBind/Diffbind_ds/CSV_files")
rnaseq_dir    <- file.path(data_dir, "DESEQ2/output_files/full_results_for_topGO")
feature_table <- file.path(data_dir, "ChIPpeakAnno/GCF_001417965.1_Aiptasia_genome_1.1_feature_table.txt")
assembly_rpt  <- file.path(data_dir, "ChIPpeakAnno/GCF_001417965.1_Aiptasia_genome_1.1_assembly_report.txt")
conv_file     <- file.path(data_dir, "Aiptasia_Genome/LOC_to_AIPGENE_conversion.csv")

output_root      <- file.path(data_dir, "GO_Enrichment/topGO_PeakToGene/topGO_ds_promoterOverlap")
tsv_dir          <- file.path(output_root, "TSV_files")
intersection_dir <- file.path(output_root, "Intersection_gene_lists")
for (d in c(tsv_dir, intersection_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Promoter-proximal region defined as ± 3 kb of TSS for a given gene
promoter_bp  <- 3000   

chip_alpha   <- 0.05
rna_padj_cut <- 0.05
rna_lfc_cut  <- 1

contrasts <- c("CS_vs_CA", "HS6_vs_HA6", "HS18_vs_HA18",
               "HA6_vs_CA", "HA18_vs_CA", "HS6_vs_CS", "HS18_vs_CS",
               "Interaction_6h", "Interaction_18h")

# Contrasts of interest for this integration of datasets
Between_symState_contrasts <- c("CS_vs_CA", "HS6_vs_HA6", "HS18_vs_HA18")

###############################################################################
# GENOME ANNOTATION: Identify genomic locations of TSS's and promoter-proximal windows for each gene
# GenBank (DiffBind peaks) -> RefSeq (feature table) scaffold names

asm_lines <- read_lines(assembly_rpt)
hdr_idx   <- which(grepl("^#?\\s*Sequence-Name\\b", asm_lines))[1]
gb_to_rs  <- read_delim(assembly_rpt, delim = "\t", skip = hdr_idx - 1,
                        col_types = cols(.default = col_character()), trim_ws = TRUE) %>%
  transmute(gb = `GenBank-Accn`, rs = `RefSeq-Accn`) %>%
  filter(!is.na(gb), !is.na(rs), gb != "", rs != "")

ft <- read_tsv(feature_table, show_col_types = FALSE)
names(ft)[1] <- sub("^#\\s*", "", names(ft)[1])

gene_df <- ft %>%
  filter(feature == "gene") %>%
  transmute(seqnames  = genomic_accession,
            strand    = if_else(strand %in% c("+", "-"), strand, "*"),
            tss       = if_else(strand == "+", as.integer(start), as.integer(end)),
            LOC_id    = coalesce(symbol, paste0("GeneID:", GeneID), name),
            gene_desc = name) %>%
  drop_na(seqnames, tss, LOC_id)

tss_gr      <- GRanges(gene_df$seqnames, IRanges(gene_df$tss, gene_df$tss), strand = gene_df$strand)
promoter_gr <- promoters(tss_gr, upstream = promoter_bp, downstream = promoter_bp)

###############################################################################
# RNA-SEQ RESULTS (AIPGENE -> LOC_id)

conv <- read_csv(conv_file, show_col_types = FALSE) %>% distinct(AIPGENE, LOC_ID) %>% rename(LOC_id = LOC_ID)

load_rna <- function(comp) {
  read_csv(file.path(rnaseq_dir, paste0(comp, "_ALL_results.csv")), show_col_types = FALSE) %>%
    inner_join(conv, by = c("gene_id" = "AIPGENE")) %>%
    # For instances where there were multiple AIPGENE IDs assigned to one LOC_id, I kept the one with the highest baseMean value
    group_by(LOC_id) %>%
    slice_max(baseMean, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(rna_sig       = !is.na(padj) & padj < rna_padj_cut & abs(log2FoldChange) > rna_lfc_cut,
           rna_direction = case_when(rna_sig & log2FoldChange > 0 ~ "Up",
                                     rna_sig & log2FoldChange < 0 ~ "Down",
                                     TRUE ~ "NS"))
}

###############################################################################
# PEAK-TO-GENE MAPPING (peak midpoint within any promoter window)

map_peaks_to_genes <- function(peak_df) {
  seqs <- if_else(peak_df$chr %in% gb_to_rs$gb, gb_to_rs$rs[match(peak_df$chr, gb_to_rs$gb)], peak_df$chr)
  mid  <- floor((peak_df$start + peak_df$end) / 2)
  peak_gr <- GRanges(seqs, IRanges(mid, mid))

  # Drop scaffolds absent from the annotation (avoids seqlevel warnings)
  keep    <- seqs %in% as.character(seqnames(promoter_gr))
  hits    <- findOverlaps(peak_gr[keep], promoter_gr, ignore.strand = TRUE)
  q       <- which(keep)[queryHits(hits)]
  s       <- subjectHits(hits)

  tibble(peak_id     = peak_df$peak_id[q],
         chip_padj   = peak_df$padj[q],
         chip_log2FC = peak_df$log2FoldChange[q],
         LOC_id      = gene_df$LOC_id[s],
         gene_desc   = gene_df$gene_desc[s],
         dist_to_tss = distance(peak_gr[q], tss_gr[s], ignore.strand = TRUE))
}

###############################################################################
# PER-CONTRAST INTERSECTION (significant DiffBind peak in promoter-proximal region AND is a DEG)

summary_rows <- list()
all_genes    <- list()

for (comp in contrasts) {
  outdir <- file.path(tsv_dir, comp)
  dir.create(outdir, showWarnings = FALSE)

  peak_df <- read.csv(file.path(diffbind_dir, paste0(comp, "_FULL.csv"))) %>%
    mutate(peak_id = row_number())
  rna <- load_rna(comp)

  p2g_sig  <- map_peaks_to_genes(peak_df) %>% filter(!is.na(chip_padj), chip_padj < chip_alpha)
  rna_degs <- rna$LOC_id[rna$rna_sig]

  intersection_tbl <- p2g_sig %>%
    filter(LOC_id %in% rna_degs) %>%
    left_join(rna %>% dplyr::select(LOC_id, rna_log2FC = log2FoldChange, rna_padj = padj, rna_direction),
              by = "LOC_id") %>%
    arrange(chip_padj)

  write_tsv(intersection_tbl, file.path(outdir, paste0(comp, "_ChIP_RNA_intersection.tsv")))
  all_genes[[comp]] <- intersection_tbl %>% mutate(contrast = comp)

  summary_rows[[comp]] <- tibble(
    contrast                = comp,
    n_sig_peaks             = sum(!is.na(peak_df$padj) & peak_df$padj < chip_alpha),
    n_sig_peaks_in_promoter = n_distinct(p2g_sig$peak_id),
    n_genes_near_sig        = n_distinct(p2g_sig$LOC_id),
    n_rna_universe          = nrow(rna),
    n_rna_degs              = length(rna_degs),
    n_chip_rna_genes        = n_distinct(intersection_tbl$LOC_id))
}

master_summary <- bind_rows(summary_rows)
print(master_summary)
  # contrast        n_sig_peaks n_sig_peaks_in_promoter n_genes_near_sig
  # CS_vs_CA               1448                     653              864
  # HS6_vs_HA6              282                     113              135
  # HS18_vs_HA18           5636                    2642             3260

write_tsv(master_summary, file.path(output_root, "master_contrast_summary.tsv"))

###############################################################################
# ACTIVATION / REPRESSION TARGET GENE COUNTS
# One row per gene; if a gene has several peaks, the lowest-padj peak sets the ChIP peak direction

all_genes <- bind_rows(all_genes) %>%
  distinct(contrast, LOC_id, gene_desc, chip_log2FC, rna_log2FC, rna_padj, rna_direction, dist_to_tss, .keep_all = TRUE)
write_csv(all_genes, file.path(intersection_dir, "ALL_contrasts_intersection_genes.csv"))
nrow(all_genes) 
  # 752 gene-peak-contrast rows

direction_summary <- all_genes %>%
  group_by(contrast, LOC_id) %>%
  slice_min(chip_padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(chip_dir = if_else(chip_log2FC > 0, "up", "down")) %>%
  count(contrast, chip_dir, rna_direction) %>%
  mutate(category = paste0("ChIP_", chip_dir, "_RNA_", tolower(rna_direction))) %>%
  dplyr::select(contrast, category, n) %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0)

print(direction_summary)
write_csv(direction_summary, file.path(output_root, "ChIP_RNA_direction_summary_uniqueGenes.csv"))
    # contrast      ChIP_up_RNA_up (activation)  ChIP_up_RNA_down (repression)
    # CS_vs_CA        95                           81
    # HS6_vs_HA6      16                           28
    # HS18_vs_HA18   144                          302

###############################################################################
# CONTRASTS OF INTEREST (between symbiotic states at each temperature) 
# gene overlap across contrasts

Between_symState_deg_overlap <- all_genes %>%
  filter(contrast %in% Between_symState_contrasts) %>%
  distinct(LOC_id, gene_desc, contrast, rna_direction) %>%
  arrange(LOC_id, contrast) %>%
  group_by(LOC_id, gene_desc) %>%
  summarise(n_contrasts   = n_distinct(contrast),
            contrasts_sig = paste(contrast, collapse = ";"),
            directions    = paste(rna_direction, collapse = ";"),
            .groups = "drop") %>%
  arrange(desc(n_contrasts), gene_desc)

write_csv(Between_symState_deg_overlap, file.path(intersection_dir, "Between_symState_DEG_peak_overlap.csv"))

count(Between_symState_deg_overlap, n_contrasts)
# 625 genes total now
# 587 in 1 contrast, 35 in 2 contrasts, 3 in all 3 contrasts

# See how many of these target genes are uncharacterized
mean(grepl("uncharacterized", Between_symState_deg_overlap$gene_desc, ignore.case = TRUE)) #49.6% (310 / 625)

