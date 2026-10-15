

#### effective_gen_size.sh ####
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G

grep -v "^>" GCA_001417965.1_Aiptasia_genome_1.1_genomic.fna | grep -o '[A-Z]' | wc -l
