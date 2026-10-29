# Aiptasia gene ID conversion — AIPGENE (Baumgarten et al. 2015 gene models; GenBank GCA_001417965.1) to LOC (NCBI RefSeq GCF_001417965.1)
# RNA-seq was quantified against AIPGENE gene models, while ChIP-seq peaks and GO terms use RefSeq LOC IDs. 
# Both annotations are on the same genome assembly, so genes are matched by genomic coordinates.
# Each AIPGENE gene is assigned to the RefSeq gene it overlaps by the most base pairs on the same strand.
# Output is used by DESeq2 topGO, peak-to-gene, and Pearson correlation scripts.

library(readr)
library(dplyr)
library(GenomicRanges)

data_dir        <- "/Users/oliviaburleigh/Documents/SEQ_DATA/Aiptasia_Genome"
gca_feature_tbl <- file.path(data_dir, "GCA_001417965.1_Aiptasia_genome_1.1_feature_table.txt")
gcf_feature_tbl <- file.path(data_dir, "GCF_001417965.1_Aiptasia_genome_1.1_feature_table.txt")
assembly_rpt    <- file.path(data_dir, "GCF_001417965.1_Aiptasia_genome_1.1_assembly_report.txt")

###############################################################################
# GENBANK (GCA) GENES: AIPGENE locus tags on LJWW scaffolds

gca <- read_tsv(gca_feature_tbl, show_col_types = FALSE)
names(gca)[1] <- sub("^#\\s*", "", names(gca)[1])

gca_genes <- gca %>%
  filter(feature == "gene") %>%
  transmute(chr_gca  = genomic_accession,
            start    = as.integer(start),
            end      = as.integer(end),
            strand   = strand,
            AIPGENE  = sub("^AC249_", "", locus_tag),   # AC249_AIPGENE13054 -> AIPGENE13054 (matches RNA-seq IDs)
            name_gca = name)

nrow(gca_genes) # 26,789

###############################################################################
# REFSEQ (GCF) GENES: LOC IDs on NW scaffolds

gcf <- read_tsv(gcf_feature_tbl, show_col_types = FALSE)
names(gcf)[1] <- sub("^#\\s*", "", names(gcf)[1])

gcf_genes <- gcf %>%
  filter(feature == "gene") %>%
  transmute(chr_gcf  = genomic_accession,
            start    = as.integer(start),
            end      = as.integer(end),
            strand   = strand,
            LOC_ID   = if_else(grepl("^LOC", symbol), symbol, paste0("LOC", GeneID)),
            name_gcf = name)

nrow(gcf_genes) # 24,862

###############################################################################
# SCAFFOLD NAMES: GenBank (LJWW) -> RefSeq (NW) from the NCBI assembly report

asm_lines <- read_lines(assembly_rpt)
hdr_idx   <- which(grepl("^#?\\s*Sequence-Name\\b", asm_lines))[1]

gb_to_rs <- read_delim(assembly_rpt, delim = "\t", skip = hdr_idx - 1,
                       col_types = cols(.default = col_character()), trim_ws = TRUE) %>%
  transmute(chr_gca = `GenBank-Accn`, chr_gcf = `RefSeq-Accn`) %>%
  filter(chr_gca != "", chr_gcf != "") %>%
  distinct()

gca_genes <- gca_genes %>%
  inner_join(gb_to_rs, by = "chr_gca")

###############################################################################
# COORDINATE OVERLAP (same strand; keep the RefSeq (LOC) gene with the largest overlap)

gca_gr <- GRanges(gca_genes$chr_gcf, IRanges(gca_genes$start, gca_genes$end), strand = gca_genes$strand)
gcf_gr <- GRanges(gcf_genes$chr_gcf, IRanges(gcf_genes$start, gcf_genes$end), strand = gcf_genes$strand)

hits <- findOverlaps(gca_gr, gcf_gr, ignore.strand = FALSE)

overlaps <- tibble(gca_idx    = queryHits(hits),
                   gcf_idx    = subjectHits(hits),
                   overlap_bp = width(pintersect(gca_gr[queryHits(hits)], gcf_gr[subjectHits(hits)]))) %>%
  group_by(gca_idx) %>%
  slice_max(overlap_bp, n = 1, with_ties = FALSE) %>%
  ungroup()

gene_mapping <- tibble(AIPGENE    = gca_genes$AIPGENE[overlaps$gca_idx],
                       LOC_ID     = gcf_genes$LOC_ID[overlaps$gcf_idx],
                       name_gcf   = gcf_genes$name_gcf[overlaps$gcf_idx],
                       name_gca   = gca_genes$name_gca[overlaps$gca_idx],
                       chr_gcf    = gca_genes$chr_gcf[overlaps$gca_idx],
                       strand     = gca_genes$strand[overlaps$gca_idx],
                       overlap_bp = overlaps$overlap_bp)

# AIPGENE genes matched to a RefSeq (LOC) gene
nrow(gene_mapping) # 21,257

# Unique RefSeq (LOC) genes (several AIPGENEs can map to one LOC)
n_distinct(gene_mapping$LOC_ID) # 19,027

###############################################################################
# SAVE FILES

# Full table (with gene descriptions and overlap size)
write_csv(gene_mapping, file.path(data_dir, "Aiptasia_AIPGENE_LOC_mapping_GCA_GCF.csv"))

# Simplified version of conversion table used downstream
write_csv(gene_mapping %>% distinct(AIPGENE, LOC_ID, name_gcf),
          file.path(data_dir, "LOC_to_AIPGENE_conversion.csv"))
