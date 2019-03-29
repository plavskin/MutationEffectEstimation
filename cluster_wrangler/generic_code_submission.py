#!/usr/bin/python

'''
Generic function to read in setup_file and run code specified therein
'''

import csv
import cluster_functions
import os
import sys
import getopt
from shutil import copyfile

usage = '''

Usage:

	generic_code_submission -s setup_file

'''

class CodeSubParameters(cluster_functions.InputParameterHolder):
	'''
	Holds parameters related to code being submitted
	'''
	def __init__(self, parameter_list):
		required_key_list = ['code_path', 'code_to_run', 'job_number', \
			'output_file_prefix', 'output_file_extension', 'module']
		super(CodeSubParameters, self).__init__(parameter_list, \
			required_key_list)
	def set_code(self, code_name):
		current_option_dict = {'current_code': code_name}		
		code_list = self.input_val_dict['code_to_run']
		code_idx = code_list.index(code_name)
		number_of_codes = len(code_list)
		self.current_option_dict = self._select_sublist(self.input_val_dict, \
			current_option_dict, code_idx, number_of_codes)
	def get_option(self, key):
		return(self.current_option_dict[key])
	def get_option_dict(self):
		return(self.current_option_dict)


class GenericCodeSubmitter(cluster_functions.CodeSubmitter):
	'''
	Submits info to cluster_wrangler.cluster_functions.job_flow_handler
	to run code that performs Maximum Likelihood Estimation
	'''
	def __init__(self, cluster_parameters, cluster_folders, \
		code_sub_parameters, job_name, output_path):
		self.code_sub_parameters = code_sub_parameters
		completefile = \
			os.path.join(cluster_folders.get_path('completefile_path'), \
				'generic_completefile.txt')
		job_numbers = [x + 1 for x in \
			list(range(code_sub_parameters.get_option('job_number')))]
		module = code_sub_parameters.get_option('module')
		parallel_processors = cluster_parameters.get_input_option('parallel_processors')
		code_name = code_sub_parameters.get_option('code_to_run')
		output_file_prefix = code_sub_parameters.get_option('output_file_prefix')
		output_extension = code_sub_parameters.get_option('output_file_extension')
		additional_beginning_lines_in_job_sub = []
		additional_end_lines_in_job_sub = []
		initial_sub_time = cluster_parameters.get_input_option('current_time')
		initial_sub_mem = cluster_parameters.get_input_option('current_mem')
		self.within_batch_counter_call = \
			cluster_parameters.get_batch_counter_call()
		self.output_file = os.path.join(output_path, output_file_prefix + '_' + \
			self.within_batch_counter_call + '.' + output_extension)
		code_path = code_sub_parameters.get_option('code_path')
		# run __init__ from parent class, which in turn runs
			# _create_code_run_input_lists
		super(GenericCodeSubmitter, self).__init__(cluster_parameters, \
			cluster_folders, completefile, job_name, \
			job_numbers, module, parallel_processors, \
			output_extension, code_name, \
			additional_beginning_lines_in_job_sub, \
			additional_end_lines_in_job_sub, initial_sub_time, \
			initial_sub_mem, output_path, output_file_prefix, code_path)
	def _create_code_run_input_lists(self):
		'''
		Creates list of keys and their values to be submitted to
		external code being run
		'''
		input_dict = self.code_sub_parameters.get_option_dict()
		self.key_list = list(input_dict.keys()) + ['external_counter', 'output_file']
		self.value_list = list(input_dict.values()) + \
			[self.within_batch_counter_call, self.output_file]

def set_csv_fieldsize_limit():
	'''
	Necessary to avoid errors in reading csv files with large fields
	'''
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

def check_input():
	# gets input passed to program, which should contain filename of
		# setup file
	if len(sys.argv) <1 :
		print >> sys.stderr, 'Setup file needed.'
		print >> sys.stderr, usage
	else:
		opts, args = getopt.gnu_getopt(sys.argv[1:], "s:h", ["--help"])
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

def main():
	set_csv_fieldsize_limit()
	setup_file = check_input()
	initial_parameter_list = cluster_functions.parse_setup_file(setup_file)
	# path of file that determines whether MLE_runner already running
	# check if any current processes for this user include MLE_runner
	code_run_counts = int(subprocess.check_output(
			('ps aux | egrep ^' + initial_parameter_list['username'] + \
				' | grep MLE_runner.py | grep -v "grep" | wc -l'), shell=True))
	code_currently_running = code_run_counts > 1
	# Loop through all directories in composite_data_dir
	for experiment_folder_name in \
		os.walk(initial_parameter_list['composite_data_path']).next()[1]:
		experiment_path = \
			os.path.join(initial_parameter_list['composite_data_path'], \
				experiment_folder_name)
		complete_checkfile = \
			os.path.join(experiment_path,'processing_complete.txt')
		# Only run MLE_runner if it's not already running and if this
			# folder hasn't been processed yet
		if not code_currently_running and not os.path.isfile(complete_checkfile):
			# create MLE_running.txt so no new instances of MLE_runner.py run
			# If setup_file doesn't exist in current directory, copy it to
				# experiment_path
			local_setup_file = os.path.join(experiment_path,'setup_file.csv')
			if not os.path.isfile(local_setup_file):
				copyfile(setup_file, local_setup_file)
			# Read setup_file from experiment_path
			parameter_list = cluster_functions.parse_setup_file(local_setup_file)
			cluster_parameters = cluster_functions.ClusterParameters(parameter_list)
			code_sub_parameters = CodeSubParameters(parameter_list)
			# set up necessary folders
			cluster_folders = \
				cluster_functions.FolderManager(cluster_parameters, \
					experiment_folder_name)
			code_to_run_list = code_sub_parameters.get_input_option('code_to_run')
			code_completeness_tracker = \
				cluster_functions.CompletenessTracker(code_to_run_list)
			for current_code in code_to_run_list:
				current_output_folder = current_code + '_results'
				current_output_path = \
					os.path.join(experiment_path, current_output_folder)
				if not os.path.isdir(current_output_path):
					os.makedirs(current_output_path)
				code_sub_parameters.set_code(current_code)
				code_submitter = GenericCodeSubmitter(cluster_parameters, \
					cluster_folders, code_sub_parameters, experiment_folder_name, \
					current_output_path)
				code_submitter.run_job_submission()
				code_completefile = code_submitter.get_completefile_path()
				code_completeness_tracker.update_key_status(current_code, \
					code_completefile)
			experiment_completeness = code_completeness_tracker.get_completeness()
			if experiment_completeness:
				open(complete_checkfile,'w+').close()
						
main()