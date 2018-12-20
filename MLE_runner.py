#!/usr/bin/python

# Estimates mutational effects and proportion petites for each strain
	# based on differences in average growth rates between mean test
	# and reference strain GRs in a group within a field

import sys, os, getopt
import csv
import subprocess
import math
#import numpy
#import re 
from shutil import copyfile
from cluster_wrangler import cluster_functions
from mlestimator import mle_functions
from mlestimator import mle_run_functions
from mlestimator import mle_sim_functions
#import warnings as war
#war.filterwarnings("ignore", message="numpy.dtype size changed")

usage = '''

Usage:

    MLE_runner_slurm_GRdiff.py -s setup_file

Defaulting to ~/setup_data_MLE.csv

'''

# test

# code below necessary to avoid errors in reading csv files with large fields
	# from user1251007 on stackoverflow
maxInt = sys.maxsize
decrement = True

while decrement:
    # decrease the maxInt value by factor 10 
    # as long as the OverflowError occurs.

    decrement = False
    try:
        csv.field_size_limit(maxInt)
    except OverflowError:
        maxInt = int(maxInt/10)
        decrement = True

def Check_Input():
	# gets input passed to program, which should contain filename of
		# setup file
	if len(sys.argv) <1 :
		print >> sys.stderr, 'Setup file needed.'
		print >> sys.stderr, usage
	else:
		opts, args = getopt.gnu_getopt(sys.argv[1:], "s:h", ["--help"])

	setup_file = '~/setup_data.csv'
	try:
		for o, a in opts:
			if o == "-h" or o == "--help":
				print usage
				exit()
			if o == '-s':
				setup_file = a
	except:
		print >> sys.stderr, 'Setup file needed.'
		print >> sys.stderr, usage
		exit()

	return (setup_file)

##############################################################################
# main

setup_file = Check_Input()
	# setup_file is given as input when running the code, and contains
		# info about directories containing code and data

initial_parameter_list = cluster_functions.parse_setup_file(setup_file)

# path of file that determines whether MLE_runner already running
currently_running_checkfile = os.path.join(initial_parameter_list['composite_data_folder'],\
	'MLE_running.txt')

# Loop through all directories in composite_data_dir
for experiment_folder_name in os.walk(initial_parameter_list['composite_data_folder']).next()[1]:

	experiment_path = \
		os.path.join(initial_parameter_list['composite_data_folder'], \
			experiment_folder_name)
	complete_checkfile = \
		os.path.join(experiment_path,'processing_complete.txt')

	# Only run MLE_runner if it's not already running and if this
		# folder hasn't been processed yet
	if not os.path.isfile(currently_running_checkfile) and \
		not os.path.isfile(complete_checkfile):

		# If setup_file doesn't exist in current directory, copy it to
			# experiment_path
		local_setup_file = os.path.join(experiment_path,'setup_file.csv')
		if not os.path.isfile(local_setup_file):
			copyfile(setup_file, local_setup_file)

		# Read setup_file from experiment_path
		parameter_list = cluster_functions.parse_setup_file(local_setup_file)
		cluster_parameters = cluster_functions.ClusterParameters(parameter_list)
		# get general info necessary to run the rest of code
		mle_parameters = mle_functions.MLEParameters(parameter_list)
		sim_parameters = mle_sim_functions.SimParameters(parameter_list)

		# create MLE_running.txt so no new instances of MLE_runner.py run
		open(currently_running_checkfile,'w+').close()

		# Get names of necessary folders, and, if this hasn't been done
			# before, create those folders
		cluster_folders = cluster_functions.FolderManager(cluster_parameters, \
			experiment_folder_name)
		mle_folders = mle_functions.FolderManager(cluster_parameters, \
			cluster_folders, experiment_folder_name)

		##########
		# Run data analysis

#		try:
#			phenotype_file_list = subprocess.check_output('ls -lrt '
#				+ os.path.join(('"' + experiment_path + '"'),'*_phenotype_file.csv'),shell=True)
#		except subprocess.CalledProcessError:
#			phenotype_file_list = ''

		############### ??? TO DO ??? ###############
#		rep_list = ['1','2']
		rep_list = ['1']
		# somehow create starting_parameter_val_file
		# as output of data analysis, create phenotype_file
		############### ??? TO DO ??? ###############

		

		############### ??? TO DO ??? ###############
		# Update file types being taken in

		# Check for data analysis complete_file before running anything in this loop
		############### ??? TO DO ??? ###############

	#	current_input_data_prefix = os.path.join(cluster_folders.experiment_path, \
	#		(current_growth_condition + '_'))
	
		##########
		# Run MLE - fixed and unfixed
		# Run 2 reps - one 'fast' with standard grid_power and L,
			# and a second one to ensure same results when
			# grid_power and L are increased
		rep_completeness_tracker = cluster_functions.CompletenessTracker(rep_list)
		for rep_index, rep in enumerate(rep_list):
			rep_float = float(rep)

			# set L and gridpower values for FFT
#			current_L = parameter_list["starting_L"]*pow(1.5,rep_float-1)
				# returns starting_L for rep1 and starting_L*1.5 for rep2
#			current_gridpower = parameter_list["starting_gridpower"]+(rep_float-1)

			# set memory and time for current run
			cluster_parameters.set_current_time(rep_float*cluster_parameters.starting_time)
			cluster_parameters.set_current_mem(rep_float*cluster_parameters.starting_mem)

			# set subfolder to be in current rep
			mle_folders.set_current_output_subfolder('rep_' + rep)

#				current_petite_file = os.path.join(cluster_folders.get_path('experiment_path'), \
#					current_growth_condition + '_petite_GR_data.csv')

#				# For 'sim', mat_phenotype_file needs to include the job iteration
#				if rep == 'sim':
#					########### FIX PATH AND ALSO SLURM PART HERE!!!!!!
#					current_phenotype_file = os.path.join(cluster_folders.get_path('experiment_path'), \
#						(current_growth_condition + '_phenotype_file_' + '${SLURM_ARRAY_TASK_ID}.csv'))
#				else:
#					current_phenotype_file = os.path.join(cluster_folders.get_path('experiment_path'), \
#						(current_growth_condition + '_phenotype_file.csv'))

			# use additional_code_run_keys and values to specify where input
				# data comes from (and any other extra information that
				# doesn't come from setup_file)
			additional_code_run_keys = []
			additional_code_run_values = []			
#			additional_code_run_keys = ['L', 'gridpower']
#			additional_code_run_values = [current_L, current_gridpower]
			#output_id_string_start = '_'.join([current_growth_condition, rep])
			output_id_string_start = rep

			# run MLE and LL profile creation, as well as CI
				# identification, for the current model
			model_loop_completeness = mle_run_functions.loop_over_models(mle_parameters, \
				cluster_parameters, cluster_folders, mle_folders, \
				experiment_path, additional_code_run_keys, \
				additional_code_run_values, output_id_string_start, \
				sim_parameters)
			# run model comparisons
			if model_loop_completeness:
				model_comparison_completeness = \
					mle_run_functions.loop_over_model_comparisons(mle_folders, \
						sim_parameters, cluster_folders, cluster_parameters, \
						output_id_string_start, additional_code_run_keys, \
						additional_code_run_values)
			# update current rep completeness
			current_rep_completeness = \
				model_loop_completeness and model_comparison_completeness
			rep_completeness_tracker.switch_key_completeness(rep, \
				current_rep_completeness)

		rep_completeness = rep_completeness_tracker.get_completeness()
		if rep_completeness:
			open(complete_checkfile,'a').close()




		# remove MLE_running.txt so MLE_runner.py can run again
		os.remove(currently_running_checkfile)




