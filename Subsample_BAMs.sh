#Subsample BAMs to target depth. Any samples below the threshold will not be subsampled.

############################# subsample_bams.sh #############################
#!/bin/bash
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.err
#SBATCH --array=28
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Subsample ChIP-seq BAM files before MACS3 peak calling
# SMAD3 ChIP target depth: 25 million reads 
# Input target depth: 40 million reads  
# Mock target: 20 million reads  

config="/nfs4/IB/Weis_Lab/burleigo/raw_seq_data/config.txt"
sample=$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} '$1==ArrayTaskID {print $2}' "${config}")

samtools="/local/cqls/opt/x86_64/bin/samtools"

input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/clean_reads"
output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/align_and_index/alignedreads/downsampled_bams_median"
mkdir -p "${output_dir}"

SEED=42

# Make downsample_bam function
# Required arguments: $1=input_bam  $2=output_bam  $3=target_reads

downsample_bam() {
    local in_bam="$1"
    local out_bam="$2"
    local target="$3"

    # Count total mapped reads
    local total
    total=$("${samtools}" view -c -F 2308 "${in_bam}")
    echo "Total reads: ${total}"

    if (( total <= target )); then
        echo "  Below target (${target}). Copying as-is."
        cp "${in_bam}" "${out_bam}"
        cp "${in_bam}.bai" "${out_bam}.bai" 2>/dev/null || true
    else
        # Calculate fraction and format for samtools -s SEED.FRAC
        local subsample_arg
		
        subsample_arg=$(awk -v t="${target}" -v n="${total}" -v s="${SEED}" \
            'BEGIN { printf "%d.%07d", s, (t/n)*10000000 }')
			
        echo "Downsampling to ${target} reads (samtools -s ${subsample_arg})"
        "${samtools}" view -b -s "${subsample_arg}" "${in_bam}" \
            | "${samtools}" sort -o "${out_bam}"
    fi

    # Index the output
    "${samtools}" index "${out_bam}"

    # Verify
    local out_count
    out_count=$("${samtools}" view -c -F 2308 "${out_bam}")
    echo "Output reads: ${out_count}"
}

SMAD3_TARGET=25000000
INPUT_TARGET=40000000
MOCK_TARGET=20000000

# Run downsample_bam function for each sample type
if [[ "${sample}" == Smad3_* ]]; then
    in_bam="${input_dir}/${sample}_bt2_paired_HOST.dedup.bam"
    out_bam="${output_dir}/${sample}_bt2_paired_HOST.dedup.ds.bam"
    downsample_bam "${in_bam}" "${out_bam}" "${SMAD3_TARGET}"

elif [[ "${sample}" == Input_* ]]; then
    in_bam="${input_dir}/${sample}_bt2_paired_HOST.dedup.bam"
    out_bam="${output_dir}/${sample}_bt2_paired_HOST.dedup.ds.bam"
    downsample_bam "${in_bam}" "${out_bam}" "${INPUT_TARGET}"

elif [[ "${sample}" == Mock_* ]]; then
    in_bam="${input_dir}/${sample}_bt2_paired_HOST.dedup.bam"
    out_bam="${output_dir}/${sample}_bt2_paired_HOST.dedup.ds.bam"
    downsample_bam "${in_bam}" "${out_bam}" "${MOCK_TARGET}"

else
    echo "Sample ${sample} is not a Smad3/Input/Mock sample, skipping."
fi

