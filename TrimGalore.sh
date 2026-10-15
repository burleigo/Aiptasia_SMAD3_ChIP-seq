
# TrimGalore
######################################## raw_fastqc.sh ############################################ 
#!/bin/bash

echo "Starting FastQC at $(date)"
#module load fastqc
/local/cqls/opt/x86_64/bin/fastqc -t 4 /nfs4/IB/Weis_Lab/burleigo/raw_seq_data/BMK250108-CN204-ZX01-0101/BMK_DATA_20250220102315_1/*.gz -o ./QCOutput/fastqc/

echo "FastQC finished at $(date)"
####################################################################################################



################## 1. Adapter Trim (submit_trimgalore.sh/script_trimgalore.sh) #####################
############### a) submit_trimgalore.sh ###############
#!/bin/bash
#SBATCH --job-name=trim_galore_job  
#SBATCH --output=logs/%x_%j.out    
#SBATCH --error=logs/%x_%j.err     
#SBATCH --time=04:00:00           
#SBATCH --mem=16G                   
#SBATCH --cpus-per-task=4          

# Input and output directories
input_dir="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/BMK250108-CN204-ZX01-0101/BMK_DATA_20250220102315_1"

output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/adapter_trim/adapter_trimmed_fastq"

# Loop through all R1 files in the input directory
for R1 in "$input_dir"/*_R1.fq.gz; do
    R2="${R1/_R1.fq.gz/_R2.fq.gz}"

    sbatch --export=R1="$R1",R2="$R2",OUTPUT_DIR="$output_dir" script_trimgalore.sh
done

############### b) script_trimgalore.sh ###############
#!/bin/bash

# Trimming adapters from all samples using default adapter detection instead of providing them
# Will be calling this script in submit_trimgalore.sh file to access every file and run trimgalore

R1=${R1}
R2=${R2}
OUTPUT_DIR=${OUTPUT_DIR}

# Extract sample name
sample=$(basename "$R1" _R1.fq.gz)

#Run TrimGalore
/local/cluster/bin/trim_galore --paired "$R1" "$R2" \
    --output_dir "$OUTPUT_DIR" \
    --nextseq 20 \
    --retain_unpaired

#Output: *_R1/2_val_2.fq.gz
#Trimming mode: paired-end
#Trim Galore version: 0.6.7
#Cutadapt version: 4.9
#Number of cores used for trimming: 1
#Quality encoding type selected: ASCII+33
####################################################################################################



##################### 2. Poly A/T Trim (config_trimgalore_polyAT.sh) ########################## 
#!/bin/bash
#SBATCH --output=/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/logs/%x_%a.out
#SBATCH --error=/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/logs/%x_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-72

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample="$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' ${config})"

input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/adapter_trim/adapter_trimmed_fastq"

output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/polyAT_trimmed_fastq"

/local/cluster/bin/trim_galore \
    --paired "${input_dir}/${sample}_R1_val_1.fq.gz" "${input_dir}/${sample}_R2_val_2.fq.gz" \
    --fastqc \
    --output_dir "${output_dir}" \
    --nextseq 20 \
    -a " AAAAAAAAAA -a TTTTTTTTTT -n 3" \
    -a2 " AAAAAAAAAA -a TTTTTTTTTT -n 3" \
    --retain_unpaired
####################################################################################################
	
	
	
################# 3. Hard clip to 148 bp (clip_config_trimgalore.sh) ######################## 	
#!/bin/bash
#SBATCH --output=/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/clip_logs/%x_%a.out
#SBATCH --error=/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/clip_logs/%x_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-72

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample="$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' ${config})"

input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/polyAT_trimmed_fastq"

output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/clipped_fastq"

/local/cluster/bin/trim_galore \
    --paired "${input_dir}/${sample}_R1_val_1_val_1.fq.gz" "${input_dir}/${sample}_R2_val_2_val_2.fq.gz" \
    --fastqc \
    --output_dir "${output_dir}" \
    --hardtrim5 148 \
    --retain_unpaired
####################################################################################################


  
###################### 4. PolyG Trim (config_trimgalore_polyG.sh) ########################
#!/bin/bash
#SBATCH --output=/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/polyG_trim/logs/%x_%a.out
#SBATCH --error=/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/polyG_trim/logs/%x_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-72

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample="$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' ${config})"

input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/clipped_fastq"
output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/polyG_trim/trimmed_polyG_output"

/local/cluster/bin/trim_galore \
    --paired "${input_dir}/${sample}_R1_val_1_val_1.148bp_5prime.fq.gz" "${input_dir}/${sample}_R2_val_2_val_2.148bp_5prime.fq.gz" \
    --fastqc \
    --output_dir "${output_dir}" \
    --nextseq 20 \
    -a " GGGGGGGGG -a GGGGGGGGG -n 3" \
    -a2 " GGGGGGGGGGG -a GGGGGGGGGG -n 3" \
    --retain_unpaired
####################################################################################################	

  
	
############# 5. Read count summary at each trimming step (generate_read_counts.sh) ################ 
#!/bin/bash
#SBATCH --job-name=read_counts
#SBATCH --output=logs/read_counts_%A_%a.out
#SBATCH --error=logs/read_counts_%A_%a.err
#SBATCH --array=1-72

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample=$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' "${config}")

tmp_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/Nreads_stats/tmp_fastq_summaries"
summary_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/Nreads_stats/summary"
summary_file="${summary_dir}/read_counts_summary.tsv"
header="sample\traw\tadapter_trim_paired\tadapter_trim_unpaired\tpolyAT_trim_paired\tpolyAT_trim_unpaired\tclipped_paired\tclipped_unpaired\tpolyG_trim_paired\tpolyG_trim_unpaired"

mkdir -p "$tmp_dir"

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
  sbatch --dependency=afterok:${SLURM_ARRAY_JOB_ID} \
         --job-name=combine_read_counts \
         --output=logs/combine_read_counts_%j.out \
         --wrap="mkdir -p ${summary_dir}; printf '${header}\n' > ${summary_file}; cat ${tmp_dir}/*.tsv >> ${summary_file}"
fi

raw="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/BMK250108-CN204-ZX01-0101/BMK_DATA_20250220102315_1/${sample}_R1.fq.gz"
adapter_trimmed_paired="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/adapter_trim/adapter_trimmed_fastq/${sample}_R1_val_1.fq.gz"
adapter_trimmed_unpaired="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/adapter_trim/adapter_trimmed_fastq/unpaired/${sample}_R2_unpaired_2.fq.gz"
polyAT_trimmed_paired="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/polyAT_trimmed_fastq/${sample}_R1_val_1_val_1.fq.gz"
polyAT_trimmed_unpaired="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/polyAT_trimmed_fastq/unpaired_fastq/${sample}_val_unpaired.fq.gz"
clipped_paired="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/clipped_fastq/${sample}_R1_val_1_val_1.148bp_5prime.fq.gz"
clipped_unpaired="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/try2_polyAT_trim/unpaired_clipped_fastq/${sample}_val_unpaired.148bp_5prime.fq.gz"
polyG_trimmed_paired="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/polyG_trim/trimmed_polyG_output/${sample}_R1_val_1_val_1.148bp_5prime_val_1.fq.gz"
polyG_trimmed_unpaired="/nfs4/IB/Weis_Lab/burleigo/chip_project/trimmedreads/trim_galore/all_samples/polyG_trim/unpaired_trimmed_polyG_output/${sample}_val_unpaired.148bp_5prime_trimmed.fq.gz"

count_reads() {
  if [ -f "$1" ]; then
    zcat "$1" | awk 'NR % 4 == 1' | wc -l
  else
    echo "0"
  fi
}

echo -e "${sample}\t$(count_reads "$raw")\t$(count_reads "$adapter_trimmed_paired")\t$(count_reads "$adapter_trimmed_unpaired")\t$(count_reads "$polyAT_trimmed_paired")\t$(count_reads "$polyAT_trimmed_unpaired") \
\t$(count_reads "$clipped_paired")\t$(count_reads "$clipped_unpaired")\t$(count_reads "$polyG_trimmed_paired")\t$(count_reads "$polyG_trimmed_unpaired")" > "${tmp_dir}/${sample}.tsv"	
####################################################################################################

