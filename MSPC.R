# Multiple Sample Peak Calling (MSPC)

library(rmspc)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ComplexUpset)   # version 1.3.3

narrowPeak_dir <- "/Users/oliviaburleigh/Documents/SEQ_DATA/MSPC/DS_Smad3_Peaks"
subset_dir     <- "/Users/oliviaburleigh/Documents/SEQ_DATA/MSPC/DS_Smad3_Peaks/subset_narrowPeak"
output_dir     <- "/Users/oliviaburleigh/Documents/SEQ_DATA/MSPC/MSPC_output_ds"
output_dir_c3  <- file.path(output_dir, "output_c3")

consensus_dir    <- file.path(output_dir, "Consensus_bed")
consensus_dir_c3 <- file.path(output_dir_c3, "Consensus_bed")
confirmed_dir    <- file.path(output_dir, "Confirmed_bed")

for (d in c(consensus_dir, consensus_dir_c3, confirmed_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

condition_order <- c("CA", "HA6", "HA18", "CS", "HS6", "HS18")

set_colors <- c(CA = "#71797E", HA6 = "#71797E", HA18 = "#71797E",
                CS = "#CC5500", HS6 = "#CC5500", HS18 = "#CC5500")


##############################################
# Run MSPC 
# c = 2: peak in 2/3 replicates (DiffBind input + Fig. 3)
# c = 3: peak in 3/3 replicates (Fig. S3)

dir.create(subset_dir, showWarnings = FALSE, recursive = TRUE)

# NarrowPeak files from MACS3 output
all_np_files <- list.files(narrowPeak_dir,
                           pattern = "ds_macs3_PE_.*_peaks\\.narrowPeak$",
                           full.names = TRUE)

for (condition in condition_order) {
  cat("Processing:", condition, "\n")

  cond_files <- all_np_files[grepl(
    paste0("ds_macs3_PE_", condition, "[0-9]+_peaks\\.narrowPeak$"),
    basename(all_np_files))]

  # Keep only the 5 columns MSPC uses
  subset_files <- c()
  for (file in cond_files) {
    df  <- read.table(file)
    out <- file.path(subset_dir,
                     paste0(gsub("ds_macs3_PE_|_peaks\\.narrowPeak", "", basename(file)),
                            "_subset.narrowPeak"))
    write.table(df[, c(1:4, 8)], out,
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    subset_files <- c(subset_files, out)
  }

  mspc(input               = subset_files,
       replicateType       = "Biological",
       stringencyThreshold = 1e-5,
       weakThreshold       = 1e-2,
       c                   = 2,
       outputPath          = file.path(output_dir, condition),
       GRanges             = TRUE)

  mspc(input               = subset_files,
       replicateType       = "Biological",
       stringencyThreshold = 1e-5,
       weakThreshold       = 1e-2,
       c                   = 3,
       outputPath          = file.path(output_dir_c3, condition),
       GRanges             = TRUE)
}


##############################################################################
# Copy MSPC output files so theyre in their own folders 
    # (easier to transfer and use for downstream applications)
# Consensus beds (c=2, c=3) for overlap analysis 
# Confirmed beds (c=2) for DiffBind on HPC

for (condition in condition_order) {
  file.copy(file.path(output_dir, condition, "ConsensusPeaks.bed"),
            file.path(consensus_dir, paste0(condition, "_consensus_peaks.bed")),
            overwrite = TRUE)

  file.copy(file.path(output_dir_c3, condition, "ConsensusPeaks.bed"),
            file.path(consensus_dir_c3, paste0(condition, "_consensus_peaks_c3.bed")),
            overwrite = TRUE)

  for (rep in 1:3) {
    sample <- paste0(condition, rep)
    file.copy(file.path(output_dir, condition, paste0(sample, "_subset"), "Confirmed.bed"),
              file.path(confirmed_dir, paste0(sample, "_Confirmed.bed")),
              overwrite = TRUE)
  }
}

##############################################################################
# Build overlap matrix of consensus peaks 

build_overlap <- function(bed_dir, suffix) {
  peaks <- list()
  for (cond in condition_order) {
    df <- read.table(file.path(bed_dir, paste0(cond, suffix)))
    peaks[[cond]] <- setNames(df[, 1:3], c("seqnames", "start", "end"))
    cat(cond, ":", nrow(peaks[[cond]]), "peaks\n")
  }

  loci <- do.call(rbind, peaks) %>%
    arrange(seqnames, start) %>%
    group_by(seqnames) %>%
    mutate(grp = cumsum(start > cummax(lag(end, default = end[1])) + 1)) %>%
    group_by(seqnames, grp) %>%
    summarise(start = min(start), end = max(end), .groups = "drop") %>%
    dplyr::select(seqnames, start, end) %>%
    mutate(locus_id = row_number())

  mat <- matrix(0L, nrow = nrow(loci), ncol = length(condition_order),
                dimnames = list(loci$locus_id, condition_order))

  for (cond in condition_order) {
    cp <- peaks[[cond]]
    for (i in seq_len(nrow(loci))) {
      if (any(cp$seqnames == loci$seqnames[i] &
              cp$end      >= loci$start[i] &
              cp$start    <= loci$end[i])) {
        mat[i, cond] <- 1L
      }
    }
  }

  cbind(loci, as.data.frame(mat))
}

#############################################################################
# Consensus peaks overlap summary (c=2)

overlap_df <- build_overlap(consensus_dir, "_consensus_peaks.bed")
    #CA : 4789 peaks
    #HA6 : 4212 peaks
    #HA18 : 3842 peaks
    #CS : 5161 peaks
    #HS6 : 1626 peaks
    #HS18 : 1394 peaks

nrow(overlap_df) #17,598 merged loci across all conditions

# Loci per condition
colSums(overlap_df[, condition_order])
    #CA  HA6 HA18   CS  HS6 HS18
    #4777 4178 3833 5152 1616 1388

# Number of conditions sharing each locus
# 15,147 / 17,598 (86%) unique to one condition; 13 peaks in all 6
table(rowSums(overlap_df[, condition_order]))
    #1     2     3     4     5     6
    #15147  1807   458   134    39    13

write.csv(overlap_df,
          file.path(consensus_dir, "peak_overlap_matrix_with_coords.csv"),
          row.names = FALSE)

##################################################################
# Define condition intersections of interest for ComplexUpset plot

focal_intersections <- list(
  c("CA"), c("CS"), c("HA6"), c("HA18"), c("HS6"), c("HS18"),
  c("CA", "HA6", "HA18"),                     
  c("CS", "HS6", "HS18"),                      
  c("HA6", "HA18"),                           
  c("HS6", "HS18"),                           
  c("CA", "CS"),                               
  c("HA6", "HA18", "HS6", "HS18"),             
  c("CA", "HA6", "HA18", "CS", "HS6", "HS18")  
)

# Size of each intersection 
for (inter in focal_intersections) {
  present <- rowSums(overlap_df[, inter, drop = FALSE]) == length(inter)
  absent  <- rowSums(overlap_df[, setdiff(condition_order, inter), drop = FALSE]) == 0
  cat(paste(inter, collapse = "+"), ":", sum(present & absent), "\n")
}

##########################################################
# ComplexUpset plot function (used for Fig. 3 and Fig. S3)

set_bottom_to_top <- c("HS18", "HS6", "CS", "HA18", "HA6", "CA")

set_queries <- lapply(set_bottom_to_top, function(s) upset_query(set = s, fill = set_colors[[s]]))

plot_upset <- function(df) {
  df[set_bottom_to_top] <- lapply(df[set_bottom_to_top], as.logical)
  
  p <- upset(
    df,
    intersect          = set_bottom_to_top,
    intersections      = focal_intersections,
    sort_intersections = "descending",
    sort_sets          = FALSE,
    queries            = set_queries,
    name               = "",
    width_ratio        = 0.4,
    height_ratio       = 0.6,
    matrix    = intersection_matrix(geom = geom_point(size = 6), segment = geom_segment(linewidth = 1.8)),
    set_sizes = upset_set_size() + ylab("Peaks per condition"),
    base_annotations = list(
      "Number of peak loci" = (
        ggplot(mapping = aes(x = intersection)) +
          geom_bar(fill = "grey25", width = 0.75) +
          geom_text_repel(
            stat = "count",
            aes(label = after_stat(count), y = after_stat(count)),
            direction    = "y",
            nudge_y      = 100,
            segment.size = 0.3,
            box.padding  = 0.4,
            size         = 7
          ) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
          ylab("Number of peak loci")
      )
    )) &
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 24),
          axis.text  = element_text(size = 18))
  
  p + theme(axis.text.y = element_text(size = 20))
}

##############################################
# Fig. 3 — c=2 focal UpSet

p <- plot_upset(overlap_df)
print(p)

ggsave(file.path(consensus_dir, "Upset_consensus_peak_c=2_focal_complexupset.pdf"),
       p, width = 13, height = 13)

##############################################
# Consensus peaks overlap summary (c=3)

overlap_df_c3 <- build_overlap(consensus_dir_c3, "_consensus_peaks_c3.bed")
#CA : 527 peaks
#HA6 : 748 peaks
#HA18 : 950 peaks
#CS : 901 peaks
#HS6 : 259 peaks
#HS18 : 205 peaks

nrow(overlap_df_c3) #3309 merged loci

colSums(overlap_df_c3[, condition_order])
#CA  HA6 HA18   CS  HS6 HS18
#527  747  950  901  259  205

# 3095 / 3309 (94%) unique to one condition
table(rowSums(overlap_df_c3[, condition_order]))
#1    2    3    4    5
#3095  167   32   11    4

# Loci unique to each condition
colSums(overlap_df_c3[rowSums(overlap_df_c3[, condition_order]) == 1, condition_order])
#CA  HA6 HA18   CS  HS6 HS18
#459  642  892  767  181  154

# Size of each focal intersection (c=3)
for (inter in focal_intersections) {
  present <- rowSums(overlap_df_c3[, inter, drop = FALSE]) == length(inter)
  absent  <- rowSums(overlap_df_c3[, setdiff(condition_order, inter), drop = FALSE]) == 0
  cat(paste(inter, collapse = "+"), ":", sum(present & absent), "\n")
}

##############################################
# Fig. S3 - c=3 UpSet plot

p_c3 <- plot_upset(overlap_df_c3)
print(p_c3)

ggsave(file.path(consensus_dir_c3, "Upset_consensus_peak_c=3_focal_complexupset.pdf"),
       p_c3, width = 13, height = 13)
