#!/bin/sh

# Run by cron

source ~/.bash_profile
	# allows running through cron

module load python/intel/2.7.12
	# module load needed only on cluster
python /home/yp19/mut_effect_estimation/plavskin_code/pipeline_v2/MLE_runner.py -s /home/yp19/setup_file.csv
	# make sure to modify directories as needed


