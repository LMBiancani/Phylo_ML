#!/bin/bash
#SBATCH --job-name="FastSP"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=2   # processor core(s) per node

module purge
module load all/Java/17.0.2 

fastsp="/home/aknyshov/alex_data/andromeda_tools/FastSP/FastSP.jar"


for i in ../simulations/*/*/1
do	
	cd ${i}
	pwd
	rm fastsp_output.csv
	>fastsp_output.csv
	cat alignmentGroups/array_list.txt | while read fileline
	do
		#echo ${fileline}
		cat alignmentGroups/${fileline} | while read line
		do
		#	echo ${line}
		#	echo ${line}","$(java -jar ${fastsp} -r alignments2/${line} -e alignments3/${line} | grep "SP-Score" | cut -f2 -d" ")  >> fastsp_output.csv
			echo ${line}","$(java -jar ${fastsp} -r alignments2/${line} -e alignments3/${line} | grep "SP-Score" | sed 's/SP-Score //g')  >> fastsp_output.csv
		done
	done
cd ../../../../3_feature_assessment/
done
