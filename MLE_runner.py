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
import pandas
import Cluster_Functions

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

class Parameter(object):
	# processes a single parameter passed to it
    def __init__(self, name, value, explanation, type):
        self.name = name
        self.explanation = explanation
        self.type = type
        self._parse_value(value, type)
    def _parse_value(self, value, type):
    	# determines what type the parameter belongs to
        if   type == "str":
            self.value = value
        elif type == "list":
            self.value = value.split(";")
        elif type == "float":
            self.value = float(value)
        elif type == "float_list":
            self.value = [float(x) for x in value.split(";")]
        else:
            raise ValueError("invalid data type: " + type)
    def concat_to_string(self):
    	# returns a string with parameter info
        concatenated_string = ""
        concatenated_string += self.name
        if self.type == "str":
            value_str = self.value
        elif self.type == "list":
            value_str = ','.join(self.value)
        elif self.type == "float":
            value_str = str(self.value)
        elif self.type == "float_list":
            value_str = ','.join(str(element) for element in self.value)
        concatenated_string += '\n'+value_str
        concatenated_string += '\n'+self.type
        concatenated_string += '\n'+self.explanation
        return(concatenated_string)

class Parameters(object):
	# creates a list of arbitrary parameters
    def __init__(self, parameters):
        self.parameters = {}
        for p in parameters:
            self.parameters[p.name] = p
    def __getitem__(self, name):
        if name in self.parameters:
            return self.parameters[name].value
        else:
            return 'this parameter does not exist'
    def print_values(self):
        for name in self.parameters:
            print(self.parameters[name].concat_to_string())

def parse_setup_file(filename):
	# parses setup file and creates a list of the parameters therein
    f = open(filename)
    lines = f.readlines()
    f.close()
    lines.pop(0)
    parameters = []
    for l in lines:
        spl = l.rstrip().split(",")
        if len(spl) != 4:
            raise ValueError("needs 4 values: " + l)
        parameters.append(Parameter(spl[0], spl[1], spl[2], spl[3].lstrip().rstrip()))
    return Parameters(parameters)

##############################################################################
# main

setup_file = Check_Input()
	# setup_file is given as input when running the code, and contains
		# info about directories containing code and data

initial_parameter_list = parse_setup_file(setup_file)

# path of file that determines whether MLE_runner already running
currently_running_checkfile = os.path.join(initial_parameter_list['home_folder'], \
	initial_parameter_list['username'],'MLE_running.txt')

# Loop through all directories in composite_data_dir
for experiment_folder_name in os.walk(initial_parameter_list['composite_data_folder']).next()[1]:

	complete_checkfile = \
		os.path.join(initial_parameter_list['composite_data_folder'], \
			experiment_folder_name,'processing_complete.txt')

	# Only run MLE_runner if it's not already running and if this
		# folder hasn't been processed yet
	if not os.path.isfile(currently_running_checkfile) and \
		not os.path.isfile(complete_checkfile):

		# If setup_file doesn't exist in current directory, copy it to
			# current_folder
		local_setup_file = os.path.join(current_folder,'setup_file.csv')
		if not os.path.isfile(local_setup_file):
			copyfile(setup_file, local_setup_file)

		# Read setup_file from current_folder
		parameter_list = parse_setup_file(local_setup_file)
		cluster_parameters = Cluster_Functions.ClusterParameters(parameter_list)
		# get general info necessary to run the rest of code

		# create MLE_running.txt so no new instances of MLE_runner.py run
		open(currently_running_checkfile,'a').close()

		# Get names of necessary folders, and, if this hasn't been done
			# before, create those folders
		cluster_folders = Cluster_Functions.FolderManager(cluster_parameters, \
			experiment_folder_name)

		##########
		# Run data analysis

		try:
			phenotype_file_list = subprocess.check_output('ls -lrt '
				+ current_folder + '/' + '*_phenotype_file.csv',shell=True)
		except subprocess.CalledProcessError:
			phenotype_file_list = ''
		growth_condition_list = re.findall(current_folder + '/(.+?)_phenotype_file.csv',
			phenotype_file_list,re.MULTILINE)
#		growth_condition_list = ['SC']
		############### ??? TO DO ??? ###############
		# somehow get growth_condition_list, possibly as output of data analysis?
#		rep_list = ['1','2']
		rep_list = ['1']
		# somehow create starting_parameter_val_file
		# as output of data analysis, create phenotype_file
		############### ??? TO DO ??? ###############

		# Run all steps for one phenotype (growth condition) at a time,
			# but loop through modes separately for each step
		for current_growth_condition in growth_condition_list:

			############### ??? TO DO ??? ###############
			# Check for data analysis complete_file before running anything in this loop
			############### ??? TO DO ??? ###############

			current_phenotype_file = current_folder + '/' \
				+ current_growth_condition + '_phenotype_file.csv'
			current_petite_file = current_folder + '/' \
				+ current_growth_condition + '_petite_GR_data.csv'
				# contains 4 columns, with data for each individual
					# measured strain in each row.
					# 1: observed phenotype;
					# 2: s.e. of observed phenotype;
					# 3. known SNP number in strain;
					# 4. proportion of genome sequenced to satisfactory
						# depth to make mutation calls
		
			##########
			# Run MLE - fixed and unfixed
			# Run 2 reps - one 'fast' with standard grid_power and L,
				# and a second one to ensure same results when
				# grid_power and L are increased

			for rep_index, rep in enumerate(rep_list):
				rep_float = float(rep)
				current_L = starting_L*pow(1.5,rep_float-1)
					# returns starting_L for rep1 and starting_L*1.5 for rep2
				current_gridpower = starting_gridpower+(rep_float-1)
				current_MLE_time = MLE_time*rep_float
				current_MLE_memory = MLE_memory*rep_float

				current_MLE_output_subfolder = MLE_output_folder + '/rep_' + rep
				if not os.path.isdir(current_MLE_output_subfolder):
					os.makedirs(current_MLE_output_subfolder)

				for current_mode_index, current_mode in enumerate(mode_list):

					print('##################################################################################')
					print(current_mode)

					# Use parameter lists appropriate to current_mode, but
						# adding 'unfixed' parameter regime to the list
					current_parameter_list = parameters_by_mode[current_mode_index]
					current_parameter_max_list = numpy.array(
						parameter_max_by_mode[current_mode_index])
					current_parameter_min_list = numpy.array(
							parameter_min_by_mode[current_mode_index])
					current_start_values = numpy.array(
							parameter_start_by_mode[current_mode_index])
					current_parameter_profile_ub_list = numpy.array(
						parameter_profile_ub_by_mode[current_mode_index])
					current_parameter_profile_lb_list = numpy.array(
							parameter_profile_lb_by_mode[current_mode_index])

					# find the total number of parameters, including fixed ones
					total_parameter_number = len(current_parameter_list)

					current_profile_points_list = numpy.array(
						profile_points_by_mode[current_mode_index])
					if len(current_profile_points_list)==1:
						current_profile_points_list = \
							numpy.repeat(current_profile_points_list,total_parameter_number)

					# Determine which parameters are always fixed
						# (the ones where min and max values are the same)
					pernamently_fixed_parameter_bool = \
						current_parameter_max_list==current_parameter_min_list

					parameters_to_loop_over_bool = numpy.invert(pernamently_fixed_parameter_bool)* \
						numpy.invert((current_profile_points_list < 1))

					non_permafixed_parameters = [item for (item,bool_val) in \
						zip(current_parameter_list,parameters_to_loop_over_bool) \
						if bool_val]
					parameters_to_loop_over = ['unfixed'] + non_permafixed_parameters
					point_numbers_to_loop_over = numpy.append([1], \
						current_profile_points_list[ \
							numpy.invert(pernamently_fixed_parameter_bool)])

	#				starting_parameter_val_file = (current_folder + '/'
	#					+ current_growth_condition + '_' + current_mode
	#					+ '_MLE_starting_vals.csv')

					MLE_summary_file = (current_folder + '/' + current_mode + '_'
						+ current_growth_condition + '_' + rep + '_MLE_file.csv')

					MLE_completefile_list = []

					for current_fixed_parameter in parameters_to_loop_over:

						current_fixed_parameter_bool_list = numpy.empty_like(pernamently_fixed_parameter_bool)
						numpy.copyto(current_fixed_parameter_bool_list,pernamently_fixed_parameter_bool)


						if current_fixed_parameter == 'unfixed':
							current_profile_point_num = 1
						else:
							# find index of current_fixed_parameter in
								# parameter list, and add it to bool
								# list of fixed parameters
							current_fixed_parameter_index = current_parameter_list.index(current_fixed_parameter)
							current_fixed_parameter_bool_list[current_fixed_parameter_index] = True
							current_profile_point_num = current_profile_points_list[current_fixed_parameter_index]

						#current_fixed_param_min = \
						#	current_parameter_min_list[current_fixed_parameter_index]
						#current_fixed_param_max = \
						#	current_parameter_max_list[current_fixed_parameter_index]

						print(current_fixed_parameter)

						current_MLE_completefile = ML_Estimator(
							completefile_folder,trackfile_folder,
							current_profile_point_num,current_mode,mode_list,username,
							current_MLE_time,current_MLE_memory,slurm_folder,
							current_fixed_parameter,
							current_fixed_parameter_bool_list,
							current_MLE_output_subfolder,
							current_growth_condition,current_L,current_gridpower,
							current_phenotype_file,current_petite_file,rep,strain_number,
							parallel_processors_max,output_dir_name,max_memory,
							max_time,current_start_values,
							current_parameter_max_list,
							current_parameter_min_list,
							pernamently_fixed_parameter_bool,
							current_parameter_list,total_parameter_number,
							ms_positions,current_parameter_profile_ub_list,
							current_parameter_profile_lb_list)

						MLE_completefile_list.append(current_MLE_completefile)

						print('done')


					# Check that MLE is complete for all parameters in this folder
					MLE_completeness_test_list = []
					for MLE_completefile_to_test in MLE_completefile_list:
						if os.path.isfile(MLE_completefile_to_test):
							MLE_completeness_test_list.append(1)
					print('completeness done')

					if sum(MLE_completeness_test_list) == len(parameters_to_loop_over):
					#if sum(MLE_completeness_test_list) >= 24:
						# run MLE on every non-fixed parameter and in unfixed mode
						
						fixed_parameter_list = parameters_to_loop_over[1:]
						#fixed_parameter_list = parameters_to_loop_over[1:24]
							# doesn't include 'unfixed'

						current_ML_combination_completefile = MLE_combiner(
							completefile_folder,current_growth_condition,
							current_mode,rep,current_folder,MLE_summary_file,
							point_numbers_to_loop_over,fixed_parameter_list,
							current_MLE_output_subfolder,LL_profile_folder,lambda_SNM,current_parameter_list)







#						# plot pdf for multiple mutations
#						use_single_mutation = False
#						if not current_mode == 'known_muts':							
#							current_pdf_completefile = pdf_Plotter(completefile_folder,trackfile_folder,
#								pdf_points,current_mode,username,
#								slurm_folder,current_folder,current_growth_condition,
#								current_L,current_gridpower,MLE_summary_file,
#								current_phenotype_file,rep,lambda_SNM,epilogue_folder,output_dir_name,max_memory,max_time,use_single_mutation)
#						# plot pdf for single mutational effect
#						use_single_mutation = True
#						single_mut_L = 'NaN'
#						single_mut_gridpower = current_gridpower+8
#						single_mut_pdf_points = pdf_points*5
#						current_pdf_completefile = pdf_Plotter(completefile_folder,trackfile_folder,
#							single_mut_pdf_points,current_mode,username,
#							slurm_folder,current_folder,current_growth_condition,
#							single_mut_L,single_mut_gridpower,MLE_summary_file,
#							current_phenotype_file,rep,lambda_SNM,epilogue_folder,output_dir_name,max_memory,max_time,use_single_mutation)
#
#						if current_mode == 'mixedgauss':
#							# also plot 'combined' single pdf whose
#								# poisson mixture is equivalent to
#								# mixedgauss distribution
#							use_single_mutation = True
#							single_mut_L = 'NaN'
#							single_mut_gridpower = current_gridpower+8
#							single_mut_pdf_points = pdf_points*5
#							current_pdf_completefile = pdf_Plotter(completefile_folder,trackfile_folder,
#								single_mut_pdf_points,'combo_poissonmut_mixedgauss',username,
#								slurm_folder,current_folder,current_growth_condition,
#								single_mut_L,single_mut_gridpower,MLE_summary_file,
#								current_phenotype_file,rep,lambda_SNM,epilogue_folder,output_dir_name,max_memory,max_time,use_single_mutation)
#						elif current_mode == 'mixed':
#							# also plot 'combined' single pdf whose
#								# poisson mixture is equivalent to
#								# mixedgauss distribution
#							use_single_mutation = True
#							single_mut_L = 'NaN'
#							single_mut_gridpower = current_gridpower+8
#							single_mut_pdf_points = pdf_points*5
#							current_pdf_completefile = pdf_Plotter(completefile_folder,trackfile_folder,
#								single_mut_pdf_points,'combo_poissonmut_mixed',username,
#								slurm_folder,current_folder,current_growth_condition,
#								single_mut_L,single_mut_gridpower,MLE_summary_file,
#								current_phenotype_file,rep,lambda_SNM,epilogue_folder,output_dir_name,max_memory,max_time,use_single_mutation)
#
#			#########
#			# Compare MLE results between the two reps?
#
#			#########
#			# Run LR: unfixed (and fixed?) AFTER MLE COMPLETE FOR ALL PARAMETERS
#
#			# Run Simulation to create simulated phenotype files
#			rep = 'sim'
#			# doesn't make sense to run simulation on knownmuts mode?
#			mode_list_random_muts = mode_list[:]
#			mode_list_random_muts.remove('knownmuts')
#
#			parameters_by_mode_random_muts = parameters_by_mode[:]
#			del parameters_by_mode_random_muts[mode_list.index('knownmuts')]
#
#			for current_mode_index, current_mode in enumerate(mode_list_random_muts):
#
#				# Determine which parameters are always fixed
#					# (the ones where min and max values are the same)
#				current_parameter_max_list = numpy.array(parameter_max_by_mode[current_mode_index])
#				current_parameter_min_list = numpy.array(parameter_min_by_mode[current_mode_index])
#				pernamently_fixed_parameter_bool = \
#					current_parameter_max_list==current_parameter_min_list
#
#				fixed_parameter_list = parameters_by_mode_random_muts[current_mode_index]
#
#				MLE_summary_file = (current_folder + '/' + current_mode + '_'
#							+ current_growth_condition + '_1_MLE_file.csv')
#
#				if os.path.isfile(MLE_summary_file):
#
#					simulation_setup_completefile = Simulation_Setup(sim_repeats,current_mode,
#						username,output_dir_name,slurm_folder,user_email,strain_number,
#						MLE_summary_file,current_growth_condition,completefile_folder,trackfile_folder,
#						sim_output_folder,current_phenotype_file,epilogue_folder,max_memory,
#						max_time,lambda_SNM)
#
#		#			for current_mode in mode_list:
#
#		#				Simulation_Setup(sim_repeats,current_mode,
#		#					mode_list,username,output_dir_name,slurm_folder,user_email,
#		#					strain_number,MLE_output_file,current_growth_condition,
#		#					completefile_folder,trackfile_folder,sim_output_folder)
#					if os.path.isfile(simulation_setup_completefile):
#						current_fixed_parameter = 'unfixed'
#						current_fixed_param_min = 0;
#						current_fixed_param_max = 0;
#
#						current_MLE_output_subfolder = MLE_output_folder + '/rep_' + rep
#						if not os.path.isdir(current_MLE_output_subfolder):
#							os.makedirs(current_MLE_output_subfolder)
#
#						starting_parameter_val_file = MLE_summary_file
#							# start parameter search at MLE values being simulated
#
#						current_MLE_phenotype_file = sim_output_folder + '/phenotype_table_'\
#							 + current_mode + '_' + current_growth_condition + '_'
#							# ML_Estimator will add id# and .csv
#
#						sim_MLE_completefile = ML_Estimator(
#										completefile_folder,trackfile_folder,
#										sim_repeats,current_mode,mode_list,username,
#										MLE_time,MLE_memory,slurm_folder,
#										current_fixed_parameter,current_fixed_param_min,
#										current_fixed_param_max,current_MLE_output_subfolder,
#										current_growth_condition,starting_L,
#										starting_gridpower,starting_parameter_val_file,
#										current_MLE_phenotype_file,rep,strain_number,
#										parallel_processors_max,epilogue_folder,output_dir_name,max_memory,max_time,pernamently_fixed_parameter_bool)
#
#						if os.path.isfile(sim_MLE_completefile):
#
#							sim_process_completefile = Simulation_Processor(
#								completefile_folder,current_mode,rep,
#								current_growth_condition,MLE_summary_file,sim_repeats,
#								fixed_parameter_list,MLE_combined_sim_folder,
#								current_MLE_output_subfolder,current_fixed_parameter)

						############### ??? TO DO ??? ###############
						# Currently, Simulation_MLE only runs an 'unfixed' MLE on simulated data - correct if necessary
						# Simulation_MLE is NOT ACTUALLY FINISHED
						############### ??? TO DO ??? ###############
		#				Simulation_MLE()




				

			# Compile LR results and calculate p-value

############### ??? TO DO ??? ###############
# Remove temporary folder on scratch?
# Remove epilogue folder
# Back up current_folder on archive
# !!! write complete_checkfile !!!
############### ??? TO DO ??? ###############

		# remove MLE_running.txt so MLE_runner.py can run again
		os.remove(currently_running_checkfile)




