# MACS3 Peak calling scripts

############# 1. Calculate effective Aiptasia genome size (effective_gen_size.sh) ################
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G

grep -v "^>" GCA_001417965.1_Aiptasia_genome_1.1_genomic.fna | grep -o '[A-Z]' | wc -l
###############################################################################################



############### 2. Call SMAD3 peaks against input with MACS3 using subsampled BAM files(subsample_macs3_p0.05.sh) #############
#!/bin/bash
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --array=55-72
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample=$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' "${config}")
basename=$(echo "${sample}" | sed 's/^Smad3_//' | sed 's/_ChIP$//')
input_sample="Input_${basename}_ChIP"

macs3="/local/cqls/opt/x86_64/bin/macs3"

input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/downsampled_bams"
output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/chipseq/macs3_peakcall/ds_p0.05_PE_peakcall_output"
mkdir -p "${output_dir}"

CHIP_BAM="${input_dir}/${sample}_bt2_paired_HOST.dedup.ds.bam"
INPUT_BAM="${input_dir}/${input_sample}_bt2_paired_HOST.dedup.ds.bam"

${macs3} callpeak \
  -t ${CHIP_BAM} \
  -c ${INPUT_BAM} \
  -f BAMPE \
  -g 2.06e8 \
  -n ds_macs3_PE_${basename} \
  --outdir ${output_dir} \
  -B \
  -p 0.05 \
  --call-summits \
  --bw 300 \
  --keep-dup all \
  --SPMR

echo "Finished MACS3 for ${basename}"
###############################################################################################


########################## 3. Count peaks per sample (count_peaks.sh) ############################
#!/bin/bash

input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/chipseq/macs3_peakcall/ds_p0.05_PE_peakcall_output"
output_file="peak_counts_ds.txt"

echo -e "Sample\tPeaks" > "$output_file"

for file in ${input_dir}/*_peaks.xls; do
    sample=$(basename "$file" _peaks.xls | sed 's/^ds_macs3_PE_//')
    count=$(grep -c '^LJWW' "$file")
    echo -e "${sample}\t${count}" >> "$output_file"
done
###############################################################################################



################### 4. Calculate fraction of reads in peaks (FRiP) (FRiP.sh) ###################
#!/bin/bash
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.err
#SBATCH --array=55-72
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample=$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' "${config}")
basename=$(echo "${sample}" | sed 's/^Smad3_//' | sed 's/_ChIP$//')

samtools="/local/cqls/opt/x86_64/bin/samtools"

bam_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/downsampled_bams"
peak_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/chipseq/macs3_peakcall/ds_p0.05_PE_peakcall_output"
results_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/chipseq/diagnostics/ds_frip_results"
mkdir -p "${results_dir}"

bam="${bam_dir}/${sample}_bt2_paired_HOST.dedup.ds.bam"
peaks="${results_dir}/${basename}_peaks.bed"

cut -f1-3 "${peak_dir}/ds_macs3_PE_${basename}_peaks.narrowPeak" | sort -u > "${peaks}"

peak_count=$(wc -l < "${peaks}")
total_reads=$("${samtools}" view -c -F 2308 "${bam}")
reads_in_peaks=$("${samtools}" view -c -F 2308 -L "${peaks}" "${bam}")
frip_pct=$(echo "scale=2; 100 * ${reads_in_peaks} / ${total_reads}" | bc)

echo -e "${sample}\t${total_reads}\t${reads_in_peaks}\t${peak_count}\t${frip_pct}" \
    > "${results_dir}/${basename}_frip.txt"

cat ${results_dir}/*_frip.txt > all_samples_frip.txt
###############################################################################################
