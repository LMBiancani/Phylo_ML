#!/bin/bash
#SBATCH --job-name="IQconcat"
#SBATCH --time=172:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node
#SBATCH --mem=200G
#SBATCH --mail-user="biancani@uri.edu" #CHANGE TO user email address
#SBATCH --mail-type=ALL

# Update Path:
OUTPUT=/data/schwartzlab/Biancani/Phylo_ML/output
iqtree_exe=/data/schwartzlab/Biancani/Software/iqtree-2.1.2-Linux/bin/iqtree2
amas=/data/schwartzlab/Biancani/Software/AMAS/amas/AMAS.py
aligned_loci_path=/data/schwartzlab/Biancani/PlacentalPolytomy/output/01_SISRS_loci_filtered

date
module purge
module load Python/3.7.4-GCCcore-8.3.0

mkdir -p $OUTPUT/all_loci
cd $OUTPUT/all_loci
pwd

mkdir alignmentGroups
cd alignmentGroups
pwd
ls ${aligned_loci_path} | rev | cut -f1 -d/ | rev | split -l 5500 - aligned_loci_list_
ls aligned_loci_list_* > array_list.txt
cd ..
pwd

mkdir iqtree_full_concattree
cd iqtree_full_concattree
pwd

files=$(cat $(sed -n p ../alignmentGroups/array_list.txt | awk '{print "../alignmentGroups/"$0}') | awk '{print "../alignments/"$0}' | paste -sd" ")
echo "Running AMAS"
python3 ${amas} concat -f fasta -d dna --out-format fasta --part-format raxml -i ${files} -t concatenated.fasta -p partitions.txt
	${iqtree_exe} -nt 20 -s concatenated.fasta -spp partitions.txt -pre inference -m MFP -bb 1000 -alrt 1000
date

