#!/bin/bash
#SBATCH --job-name="prep"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-user="biancani@uri.edu" #CHANGE TO user email address
#SBATCH --mail-type=ALL

OUTPUT=/data/schwartzlab/Biancani/Phylo_ML/output
aligned_loci_path=/data/schwartzlab/Biancani/PlacentalPolytomy/output/01_SISRS_loci_filtered

date
mkdir -p $OUTPUT
cd $OUTPUT
pwd

mkdir alignmentGroups
cd alignmentGroups
ls ${aligned_loci_path} | rev | cut -f1 -d/ | rev | split -l 250 - aligned_loci_list_
ls aligned_loci_list_* > array_list.txt
cd ..

mkdir alignments3
mkdir iqtree_genetrees2
mkdir iqtree_genetrees3
mkdir phylomad_assessment
mkdir rate_assessment
mkdir iqtree_concattree
mkdir astral_tree

date

