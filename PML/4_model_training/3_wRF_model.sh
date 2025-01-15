#!/bin/bash
#SBATCH --job-name="wrfmodel"
#SBATCH --time=72:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=2   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node
#SBATCH --mem-per-cpu=12G


cd $SLURM_SUBMIT_DIR

date


module load scikit-learn/1.1.2-foss-2022a
module load matplotlib/3.5.2-foss-2022a
module load SHAP/0.42.1-foss-2022a
module load treeinterpreter/0.2.3-foss-2022a

cd wRF_model
pwd
python3 ../train_random_forest.py -i wRFtrain_tab.tsv -f auto -d 15000 --msl 3 --mss 2 -e 6250 -t 20 
 

date
