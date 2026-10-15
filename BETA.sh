# BETA script

# Adjusted genome fasta and refgene file format to match BETA program requirements 

################## beta_format_genome.sh ########################### 
#!/bin/bash

in=/nfs4/IB/Weis_Lab/burleigo/aiptasia_genome/Aiptasia_genome.fa
out=/nfs4/IB/Weis_Lab/burleigo/aiptasia_genome/Aiptasia_genome_beta_format.fa

awk '
/^>/{
  h=substr($0,2)
  gsub(/\./,"_",h)
  if (h !~ /^chr/) h="chr"h
  print ">"h
  next
}
{print}
' "$in" > "$out"
########################################################################
    # $ head Aiptasia_genome_beta_format.fa
    # >chrLJWW01000001_1 Exaiptasia diaphana isolate CC7 scaffold1, whole genome shotgun sequence
    # GATGGGGTTATAGTCAATGGTCTGATAATGTAGTGTATTTGATAATTGGcgatgggcttcctgtgtgtaCATGTCAACAG


# File format for Aiptasia_refgene_final.txt required for BETA to run properly
      #name	chrom	strand	txStart	txEnd	name2
      #AIPGENE13054	chrLJWW01000001_1	+	564	2385	AIPGENE13054


################## BETA.sh ########################
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00

BETA="/home/zoo/burleigo/beta_wrapper.py"
PEAK_FILE="/nfs4/IB/Weis_Lab/burleigo/chip_project/chipseq/BETA/input_files/SMAD3_consensus_peaks_LibSize_norm_counts.bed"
GENOME_FASTA="/nfs4/IB/Weis_Lab/burleigo/aiptasia_genome/Aiptasia_genome_beta_format.fa"
REFGENE="/nfs4/IB/Weis_Lab/burleigo/chip_project/chipseq/BETA/Aiptasia_refgene_final.txt"
input_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/chipseq/BETA/input_files"
output_dir="/nfs4/IB/Weis_Lab/burleigo/chip_project/chipseq/BETA/output_0.05"

export PYTHONPATH="$HOME/beta_patched:$PYTHONPATH"

for expr_file in $input_dir/RNA_*_beta.txt; do
  contrast=$(basename "$expr_file" _beta.txt)
  contrast=${contrast#RNA_}                     # strip leading "RNA_"

  echo "Running BETA plus for: $contrast"

  $BETA plus \
    -p "$PEAK_FILE" \
    -e "$expr_file" \
    -k O \
    --info 1,2,3 \
    -r "$REFGENE" \
    --gs "$GENOME_FASTA" \
    --gname2 \
    -d 3000 \
    --pn 17569 \
    --df 0.05 \
    -n "$contrast" \
    -c 0.05 \
    -o ${output_dir}/"output_${contrast}"

  echo "Done: $contrast"
done
