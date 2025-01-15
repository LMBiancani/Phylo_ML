#!/bin/bash
#SBATCH --job-name="IQLiu"
#SBATCH --time=96:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=2   # number of nodes
#SBATCH --ntasks-per-node=24   # processor core(s) per node
#SBATCH --exclusive
#SBATCH --mem=250G


cd $SLURM_SUBMIT_DIR

date

#Path to IQTREE executable. Modify with path to your own executable.
iqtree_exe="/home/aknyshov/alex_data/andromeda_tools/iqtree-2.1.2-Linux/bin/iqtree2"

#Concatenate input fasta files and prepare partitions ahead of IQTree run
python3 /home/aknyshov/alex_data/andromeda_tools/AMAS/amas/AMAS.py concat -f fasta -d dna --out-format fasta --part-format raxml -i *fas -c 20 -t concatenatedTrain.fasta -p partitionsTrain.txt

#Run IQtree. Flags: -nt: use 20 CPU cores -spp: specifies partition file but allows partitions to have different evolutionary speeds -pre: specifies prefix for output files -m: determine best fit model immediately followed by tree reconstruction -bb: sets 1000 bootstrap replicates  -alrt: sets 1000 replicates to perform SH-like approximate likelihood test (SH-aLRT)
${iqtree_exe} -nt 20 -s concatenatedTrain.fasta -spp partitionsTrain.txt -pre inferenceEmpirical -m MFP -bb 1000 -alrt 1000


date
