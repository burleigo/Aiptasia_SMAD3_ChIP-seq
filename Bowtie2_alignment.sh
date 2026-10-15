# Bowtie2 alignment scripts for ChIP-seq samples.

################ 1. Build Bowtie2 index for combined reference (bowtie2_build.sh) ############## 
#!/bin/bash
#SBATCH --job-name=bowtie2_index
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00

/local/cluster/bin/bowtie2-build --threads 8 Aip_Slin_combined_genome.fna Aip_Slin_index
################################################################################################


################ 2. Bowtie2 Alignment to combined genome (Bowtie2_alignment.sh) ######################
#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH --array=1-72
#SBATCH --output=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/logs/%x_%a.out
#SBATCH --error=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/logs/%x_%a.err

#align trimmed reads in fastq format, output is a bam file (binary sequence alignment map)
config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample="$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' ${config})"
 
input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/polyG_trim/trimmed_polyG_output"
output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/paired_aligned_output"
genome_index="/nfs4/IB/Weis_Lab/burleigo/aiptasia_genome/combined_genome/Aip_Slin_index"
 
bowtie2="/local/cluster/bin/bowtie2"
samtools="/local/cqls/opt/x86_64/bin/samtools"
 
mkdir -p "$output_dir"
 
"${bowtie2}" --threads 16 \
    --no-unal \
    --score-min L,16,1 \
    --local \
    -L 16 \
    -x "${genome_index}" \
    -1 "${input_dir}/${sample}_R1_val_1_val_1.148bp_5prime_val_1.fq.gz" \
    -2 "${input_dir}/${sample}_R2_val_2_val_2.148bp_5prime_val_2.fq.gz" \
    | "${samtools}" view -b -q 30 - \
    | "${samtools}" sort -o "${output_dir}/${sample}_bt2_HostSym.sorted.bam"
 
"${samtools}" index "${output_dir}/${sample}_bt2_HostSym.sorted.bam"
#samtools -q 30 is filtering to keep only reads with a MAPQ quality score over 30, which should be almost all uniquely mapped reads
################################################################################################



############# 3. Remove symbiont-mapped reads to get host-aligned reads (separate_hostsym_reads.sh) ##########################
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
#SBATCH --array=19-72
#SBATCH --output=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/logs/%x_%a.out
#SBATCH --error=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/logs/%x_%a.err

# Host:LJWW, Sym:CAJNDV

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample="$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' ${config})"

paired_input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/paired_aligned_output"
output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/sorted_reads"

samtools="/local/cqls/opt/x86_64/bin/samtools"

# Remove SLIN-matching reads. Output is just host reads. Using these moving forward.
samtools view -h ${paired_input_dir}/${sample}'_bt2_HostSym.sorted.bam' \
    | awk '!/CAJNDV/' | samtools view -b -o ${output_dir}/${sample}'_bt2_paired_HOST.sorted.bam'
 
# Remove Exaiptasia-matching reads. Output is just symbiont reads
samtools view -h ${paired_input_dir}/${sample}'_bt2_HostSym.sorted.bam' \
    | awk '!/LJWW/' | samtools view -b -o ${output_dir}/${sample}'_bt2_paired_SYM.sorted.bam'
##########################################################################################################################



################## 4. Remove PCR duplicate reads from BAM files (bam_clean_dedup.sh) #######################
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --array=19-72
#SBATCH --output=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/logs2/%x_%a.out
#SBATCH --error=/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/logs2/%x_%a.err

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample="$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' ${config})"
 
input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/sorted_reads"
output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/clean_reads"
stats_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/stats"
tmp_dir="/fs1/home/zoo/burleigo/tmp_picard"
 
samtools="/local/cqls/opt/x86_64/bin/samtools"
picard="/local/cqls/opt/x86_64/bin/picard"
 
input_bam="${input_dir}/${sample}_bt2_paired_HOST.sorted.bam"
dedup_bam="${output_dir}/${sample}_bt2_paired_HOST.dedup.bam"
 
mkdir -p "$output_dir" "$stats_dir" "$tmp_dir"
 
"${picard}" MarkDuplicates \
    I="${input_bam}" \
    O="${dedup_bam}" \
    M="${stats_dir}/${sample}_dup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    TMP_DIR="${tmp_dir}"
 
"${samtools}" index "${dedup_bam}"
"${samtools}" flagstat "${dedup_bam}" > "${stats_dir}/${sample}_paired_HOST_dedup.flagstat.txt"
"${samtools}" idxstats "${dedup_bam}" > "${stats_dir}/${sample}_paired_HOST_dedup.idxstats.txt"
############################################################################################################



################ 5. Read count summary for BAM files (bam_generate_read_counts.sh) ###################
#!/bin/bash
#SBATCH --job-name=read_counts
#SBATCH --output=logs/bam_read_counts_%A_%a.out
#SBATCH --error=logs/bam_read_counts_%A_%a.err
#SBATCH --array=19-72

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample=$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' "${config}")
 
samtools="/local/cqls/opt/x86_64/bin/samtools"
aligned="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads"
 
tmp_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/Nreads_stats/tmp_bam_summaries"
summary_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/Nreads_stats/summary"
summary_file="${summary_dir}/bam_read_counts_summary.tsv"
header="sample\tbt2_reads\tbt2_host_reads\tbt2_host_dedup_reads"
 
mkdir -p "$tmp_dir"
 
if [ "$SLURM_ARRAY_TASK_ID" -eq 19 ]; then
  sbatch --dependency=afterok:${SLURM_ARRAY_JOB_ID} \
         --job-name=combine_bam_counts \
         --output=logs/combine_bam_counts_%j.out \
         --wrap="mkdir -p ${summary_dir}; printf '${header}\n' > ${summary_file}; cat ${tmp_dir}/*.tsv >> ${summary_file}"
fi
 
bt2="${aligned}/paired_aligned_output/${sample}_bt2_HostSym.sorted.bam"
bt2_host="${aligned}/sorted_reads/${sample}_bt2_paired_HOST.sorted.bam"
bt2_host_dedup="${aligned}/clean_reads/${sample}_bt2_paired_HOST.dedup.bam"
 
count_reads_bam() {
  "${samtools}" view -c -F 0x900 "$1"
}
 
echo -e "${sample}\t$(count_reads_bam "$bt2")\t$(count_reads_bam "$bt2_host")\t$(count_reads_bam \
"$bt2_host_dedup")" > "${tmp_dir}/${sample}.tsv"
####################################################################################################

