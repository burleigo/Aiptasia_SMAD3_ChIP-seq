#KALLISTO SCRIPTS
#Using Aiptasia + S.linuchae transcriptome for reference. 
#Aiptasia transcriptome file obtained from Aranda Lab: aip.genome_models.no_isoforms.no_duplication.mRNA.fa

#Obtained S.linuchae transcript info from NCBI. Contains other RNA types as well, not just mRNA 
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/221/635/GCA_905221635.1_Slin_CCMP2456/Slin_CCMP2456_rna_from_genomic.fna

$ cat aip.genome_models.no_isoforms.no_duplication.mRNA.fa Slin_CCMP2456_rna_from_genomic.fna > combined_RNA_Aip_Slin.fna

#Newest Kallisto version wasn't working so had to download older version v0.46.2
$ wget https://github.com/pachterlab/kallisto/releases/download/v0.46.2/kallisto_linux-v0.46.2.tar.gz
$ tar -xvzf kallisto_linux-v0.46.2.tar.gz

Location of Kallisto program: 
/nfs4/IB/Weis_Lab/burleigo/aiptasia_genome/kallisto/kallisto/kallisto


############### 1. MAKE COMBINED HOST+SYM KALLISTO INDEX (index_kallisto.sh) ############ 
#!/bin/bash
#SBATCH --job-name=kallisto_index
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00

/nfs4/IB/Weis_Lab/burleigo/aiptasia_genome/kallisto/kallisto/kallisto index -i Aip_Slin_kallisto_index.idx combined_RNA_Aip_Slin.fna


############### 2. KALLISTO ALIGNMENT OF PAIRED READS (kallisto_alignment_paired.sh) ################
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-18
#SBATCH --output=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/kallisto_map/logs/%x_%a.out
#SBATCH --error=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/kallisto_map/logs/%x_%a.err

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config_RNA.txt"
sample="$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' ${config})"

#paired-end samples used 
input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/polyG_trim/trimmed_polyG_output"

base_output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/kallisto_map/paired_aligned_output"
output_dir="${base_output_dir}/${sample}"
mkdir -p "${output_dir}"

#Using combined host+sym index
index="/nfs4/IB/Weis_Lab/burleigo/aiptasia_genome/Aip_Slin_kallisto_index.idx"

kallisto="/nfs4/IB/Weis_Lab/burleigo/aiptasia_genome/kallisto/kallisto/kallisto"

"${kallisto}" quant -i "${index}" \
    -o "${output_dir}" \
    -b 100 -t 8 \
    "${input_dir}/${sample}_R1_val_1_val_1.148bp_5prime_val_1.fq.gz" \
    "${input_dir}/${sample}_R2_val_2_val_2.148bp_5prime_val_2.fq.gz"


############### 3. PSEUDOALIGNMENT SUMMARY (summarize_kallisto_runinfo.sh) ################
#!/bin/bash

base="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/kallisto_map"
out="${base}/kallisto_runinfo_summary.tsv"

echo -e "sample\tn_processed\tn_pseudoaligned\tp_pseudoaligned\tp_unique" > "${out}"

for json in "${base}"/paired_aligned_output/paired_*_run_info.json; do
  [[ -f "$json" ]] || continue

  sample="$(basename "$json" _run_info.json | sed 's/^paired_//')"

  n_processed=$(sed -n 's/.*"n_processed":[[:space:]]*\([0-9]\+\).*/\1/p' "$json")
  n_pseudo=$(sed -n 's/.*"n_pseudoaligned":[[:space:]]*\([0-9]\+\).*/\1/p' "$json")
  p_pseudo=$(sed -n 's/.*"p_pseudoaligned":[[:space:]]*\([0-9.]\+\).*/\1/p' "$json")
  p_unique=$(sed -n 's/.*"p_unique":[[:space:]]*\([0-9.]\+\).*/\1/p' "$json")

  echo -e "${sample}\t${n_processed}\t${n_pseudo}\t${p_pseudo}\t${p_unique}" >> "${out}"
done

#### Output: kallisto_runinfo_summary.tsv
#sample	n_processed	n_pseudoaligned	p_pseudoaligned	p_unique
#CA1_RNA	27569921	17900031	64.9	62.3
#CA2_RNA	28525502	18106487	63.5	60.8
#CA3_RNA	24084017	13427181	55.8	53.6
#CS1_RNA	32185809	22302063	69.3	66.9
#CS2_RNA	33791193	23950682	70.9	68.8
#CS3_RNA	40524684	26989391	66.6	64.4
#HA181_RNA	23903337	11977084	50.1	47.5
#HA182_RNA	45965593	29663095	64.5	61.8
#HA183_RNA	25012240	16277581	65.1	62.2
#HA61_RNA	20696274	9883792	47.8	45.0
#HA62_RNA	27125689	12868706	47.4	45.2
#HA63_RNA	31992634	13382668	41.8	39.0
#HS181_RNA	33462011	23624385	70.6	68.0
#HS182_RNA	55019302	36494411	66.3	63.8
#HS183_RNA	29970797	21480069	71.7	69.2
#HS61_RNA	32203775	19536542	60.7	58.6
#HS62_RNA	25426394	16474324	64.8	62.3
#HS63_RNA	32309241	22566906	69.8	67.2


############### 4. SEPARATE AND REMOVE SYMBIONT FROM HOST READS (sep_hostsym_paired_reads.sh) ################
#!/bin/bash
#SBATCH --job-name=split_host_sym_paired_only
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-18
#SBATCH --output=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/kallisto_map/logs/%x_%a.out
#SBATCH --error=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/kallisto_map/logs/%x_%a.err

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config_RNA.txt"
sample="$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' ${config})"

base="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/kallisto_map"
in="${base}/paired_aligned_output/abundance_tsv/paired_${sample}_abundance.tsv"
out_base="${base}/paired_only_split"
host_out="${out_base}/host/host_${sample}_abundance.tsv"
sym_out="${out_base}/sym/sym_${sample}_abundance.tsv"
qc_out="${out_base}/qc/${sample}_host_sym_fraction.tsv"

mkdir -p "$(dirname "${host_out}")" "$(dirname "${sym_out}")" "$(dirname "${qc_out}")"

# Host targets: contains "AIPGENE" in header
awk 'BEGIN{FS=OFS="\t"} NR==1 || $1 ~ /^AIPGENE/ {print}' "${in}" > "${host_out}"
# Sym targets: contains "SLIN2456" in header
awk 'BEGIN{FS=OFS="\t"} NR==1 || $1 ~ /SLIN2456/ {print}' "${in}" > "${sym_out}"
# QC: fraction of est_counts assigned host vs sym (est_counts is column 4)

awk 'BEGIN{FS=OFS="\t"}
NR==1{next}
{
  c=$4
  total+=c
  if($1 ~ /^AIPGENE/) host+=c
  else if($1 ~ /SLIN2456/) sym+=c
}
END{
  hf=(total?host/total:"NA"); sf=(total?sym/total:"NA")
  print "Sample","total_est_counts","host_est_counts","sym_est_counts","host_frac","sym_frac"
  print "'"${sample}"'", total, host, sym, hf, sf
}' "${in}" > "${qc_out}"



############### 5. READ COUNT SUMMARY FOR HOST AND SYMBIONT (tsv_read_counts.sh) ################
#!/bin/bash

base="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/kallisto_map"
paired_dir="${base}/paired_aligned_output/abundance_tsv"
host_dir="${base}/paired_only_split/host"

out="${base}/kallisto_estcounts_paired_vs_host.tsv"

echo -e "sample\tpaired_total_est_counts\thost_total_est_counts\thost_fraction" > "${out}"

for f in "${paired_dir}"/paired_*_abundance.tsv; do
  [[ -f "$f" ]] || continue

  sample="$(basename "$f" _abundance.tsv | sed 's/^paired_//')"

  host_file="${host_dir}/host_${sample}_abundance.tsv"

  paired_sum=$(awk -F'\t' 'NR>1 {sum+=$4} END {print sum+0}' "$f")
  host_sum=$(awk -F'\t' 'NR>1 {sum+=$4} END {print sum+0}' "$host_file")

  frac=$(awk -v h="$host_sum" -v p="$paired_sum" 'BEGIN{if(p>0) print h/p; else print "NA"}')

  echo -e "${sample}\t${paired_sum}\t${host_sum}\t${frac}" >> "${out}"
done


#### Output: estcounts_paired_vs_host.tsv ####
#sample	paired_total_est_counts	host_total_est_counts	host_fraction
#CA1_RNA	1.79e+07	1.7899e+07	0.999944
#CA2_RNA	1.81065e+07	1.81049e+07	0.999912
#CA3_RNA	1.34272e+07	1.34267e+07	0.999963
#CS1_RNA	2.23021e+07	2.16922e+07	0.972653
#CS2_RNA	2.39507e+07	2.18699e+07	0.913122
#CS3_RNA	2.69894e+07	2.5203e+07	0.933811
#HA181_RNA	1.19771e+07	1.19751e+07	0.999833
#HA182_RNA	2.96631e+07	2.96428e+07	0.999316
#HA183_RNA	1.62776e+07	1.62771e+07	0.999969
#HA61_RNA	9.88379e+06	9.88135e+06	0.999753
#HA62_RNA	1.28687e+07	1.28681e+07	0.999953
#HA63_RNA	1.33827e+07	1.33814e+07	0.999903
#HS181_RNA	2.36244e+07	2.34637e+07	0.993198
#HS182_RNA	3.64944e+07	3.50993e+07	0.961772
#HS183_RNA	2.14801e+07	2.03703e+07	0.948334
#HS61_RNA	1.95365e+07	1.89622e+07	0.970604
#HS62_RNA	1.64743e+07	1.5962e+07	0.968903
#HS63_RNA	2.25669e+07	2.20262e+07	0.97604


