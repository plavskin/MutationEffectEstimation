#!/usr/bin/python

# Estimates mutational effects and proportion petites for each strain
	# based on differences in average growth rates between mean test
	# and reference strain GRs in a group within a field

import sys, os, getopt
import csv
import subprocess
import math
import numpy
from scipy.stats import chi2
import re
from shutil import copyfile
import Cluster_Functions
import MLE_Functions

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

initial_parameter_list = Cluster_Functions.parse_setup_file(setup_file)

# path of file that determines whether MLE_runner already running
currently_running_checkfile = os.path.join(initial_parameter_list['home_folder'],\
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
		parameter_list = Cluster_Functions.parse_setup_file(local_setup_file)
		cluster_parameters = Cluster_Functions.ClusterParameters(parameter_list)
		# get general info necessary to run the rest of code
		mle_parameters = MLE_Functions.MLEParameters(parameter_list)

		# create MLE_running.txt so no new instances of MLE_runner.py run
#		open(currently_running_checkfile,'w+').close()

		# Get names of necessary folders, and, if this hasn't been done
			# before, create those folders
		cluster_folders = Cluster_Functions.FolderManager(cluster_parameters, \
			experiment_folder_name)
		mle_folders = MLE_Functions.FolderManager(cluster_parameters, \
			cluster_folders, experiment_folder_name)

		##########
		# Run data analysis

		try:
			phenotype_file_list = subprocess.check_output('ls -lrt '
				+ os.path.join(('"' + experiment_path + '"'),'*_phenotype_file.csv'),shell=True)
		except subprocess.CalledProcessError:
			phenotype_file_list = ''
		growth_condition_list = re.findall(os.path.join(experiment_path,'(.+?)_phenotype_file.csv'),
			phenotype_file_list,re.MULTILINE)

		############### ??? TO DO ??? ###############
#		rep_list = ['1','2']
		rep_list = ['1']
		# somehow create starting_parameter_val_file
		# as output of data analysis, create phenotype_file
		############### ??? TO DO ??? ###############

		# Run all steps for one phenotype (growth condition) at a time,
			# but loop through modes separately for each step
		for current_growth_condition in growth_condition_list:

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

			for rep_index, rep in enumerate(rep_list):
				rep_float = float(rep)

				# set L and gridpower values for FFT
				current_L = parameter_list["starting_L"]*pow(1.5,rep_float-1)
					# returns starting_L for rep1 and starting_L*1.5 for rep2
				current_gridpower = parameter_list["starting_gridpower"]+(rep_float-1)

				# set memory and time for current run
				cluster_parameters.set_current_time(rep_float*cluster_parameters.starting_time)
				cluster_parameters.set_current_mem(rep_float*cluster_parameters.starting_mem)

				# set subfolder to be in current rep
				mle_folders.set_current_output_subfolder('rep_' + rep)

				current_petite_file = os.path.join(cluster_folders.experiment_path, \
					current_growth_condition + '_petite_GR_data.csv')

				# For 'sim', mat_phenotype_file needs to include the job iteration
				if rep == 'sim':
					########### FIX PATH AND ALSO SLURM PART HERE!!!!!!
					current_phenotype_file = os.path.join(cluster_folders.experiment_path, \
						(current_growth_condition + '_phenotype_file_' + '${SLURM_ARRAY_TASK_ID}.csv'))
				else:
					current_phenotype_file = os.path.join(cluster_folders.experiment_path, \
						(current_growth_condition + '_phenotype_file.csv'))

				print(mle_parameters.mode_list)
				for current_mode in mle_parameters.mode_list:

					print('##################################################################################')
					print(current_mode)

					# Use parameter lists appropriate to current_mode, but
						# adding 'unfixed' parameter regime to the list
					output_id_string = (current_growth_condition + '_' + \
						current_mode + '_' + rep)
					mle_parameters.set_mode(current_mode,output_id_string)

					##### RUN MLE #####
					output_id_string = current_mode
					mle_parameters.set_mode(current_mode,output_id_string)
					MLE_summary_file_path = os.path.join(experiment_path, \
						'_'.join([output_id_string,'MLE_file.csv']))
					# use additional_code_run_keys and values to specify where input
						# data comes from (and any other extra information that
						# doesn't come from setup_file)
					additional_code_run_keys = ['phenotype_file','petite_file','L', \
						'gridpower']
					additional_code_run_values = [current_phenotype_file, \
						current_petite_file,current_L,current_gridpower]
					# run MLE for current set of parameters
					MLE_Functions.run_MLE(mle_parameters, cluster_parameters, cluster_folders, mle_folders, \
						additional_code_run_keys, additional_code_run_values)
					# if all parameters for this mode are complete, update mode completeness
					# this also updates completeness across modes
					current_mode_mle_complete_status = \
						mle_parameters.check_completeness_within_mode()

#					if sum(MLE_completeness_test_list) == len(parameters_to_loop_over):
#					#if sum(MLE_completeness_test_list) >= 24:
#						# run MLE on every non-fixed parameter and in unfixed mode
#						
#						fixed_parameter_list = parameters_to_loop_over[1:]
#						#fixed_parameter_list = parameters_to_loop_over[1:24]
#							# doesn't include 'unfixed'
#
#						current_ML_combination_completefile = MLE_combiner(
#							completefile_folder,current_growth_condition,
#							current_mode,rep,current_folder,MLE_summary_file,
#							point_numbers_to_loop_over,fixed_parameter_list,
#							current_MLE_output_subfolder,LL_profile_folder,lambda_SNM,current_parameter_list)

		# remove MLE_running.txt so MLE_runner.py can run again
#		os.remove(currently_running_checkfile)




