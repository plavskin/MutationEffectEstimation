#!/usr/bin/python

# Contains objects needed for running Maximum Likelihood Estimation

import os
import numpy
import Cluster_Functions
import copy

class FolderManager(object):
	def __init__(self, cluster_parameters, cluster_folders, \
		experiment_folder_name):
		self.experiment_folder_name = experiment_folder_name
		self.experiment_path = \
			os.path.join(cluster_parameters.composite_data_path,experiment_folder_name)
		self.sim_output_path = os.path.join(cluster_parameters.temp_storage_path, \
			experiment_folder_name,'simulated_phenotypes')
		self.MLE_output_path = os.path.join(cluster_parameters.temp_storage_path, \
			experiment_folder_name,'MLE_output')
		self.completefile_folder = cluster_folders.completefile_path
		self.LL_profile_path = os.path.join(self.experiment_path,'/LL_profiles')
	#	self.epilogue_path = os.path.join(cluster_parameters.home_path,cluster_parameters.username,'mut_effect_epilogue_files',experiment_folder_name)
		self.MLE_combined_sim_path = os.path.join(self.experiment_path,'/MLE_sim_outputs')
		self._set_up_folders()
	def _set_up_folders(self):
		setup_complete_file = \
			os.path.join(self.completefile_folder, \
			'folder_setup_complete.txt')
		if not os.path.isfile(setup_complete_file):
			new_directory_list = (self.sim_output_path,self.MLE_output_path, \
				self.LL_profile_path,self.MLE_combined_sim_path)
			for current_new_directory in new_directory_list:
				if not os.path.isdir(current_new_directory):
					os.makedirs(current_new_directory)
			open(setup_complete_file,'a').close()
	def set_current_output_subfolder(self,current_subfolder):
		# set (and if necessary, create) a subfolder to write temp output to
		self.current_output_subfolder = os.path.join(self.MLE_output_path, \
			current_subfolder)
		if not os.path.isdir(self.current_output_subfolder):
			os.makedirs(self.current_output_subfolder)


class MLEParameters(object):
	############### ??? TO DO ??? ###############
	# Check that each mode or parameter is in the list once?
		# (currently only pays attention to the first time mode or parameter listed)
	############### ??? TO DO ??? ###############
	def __init__(self,parameter_list):
		self.sim_repeats_by_mode = parameter_list["simulation_repeats_by_mode"]
		self.profile_points_by_mode = parameter_list["profile_points_by_mode"]
		self.mode_list = parameter_list["mode_list"]
		self.parameters_by_mode = parameter_list["parameters_by_mode"]
		self.min_parameter_vals_by_mode = \
			parameter_list["min_parameter_vals_by_mode"]
		self.max_parameter_vals_by_mode = \
			parameter_list["max_parameter_vals_by_mode"]
		self.starting_parameter_vals_by_mode = \
		 	parameter_list["starting_parameter_vals_by_mode"]
		self.profile_lower_limits_by_mode = \
			parameter_list["profile_lower_limits_by_mode"]
		self.profile_upper_limits_by_mode = \
			parameter_list["profile_upper_limits_by_mode"]
		self.ms_positions = parameter_list["multistart_positions"]
		self.multistart_grid_parameters = parameter_list["multistart_grid_parameters"]
		self.logspace_profile_list_by_mode = parameter_list["logspace_profile_list"]
		self.parallel_processors = parameter_list["parallel_processors"]
		self.mode_completeness_tracker = CompletenessTracker(self.mode_list)
		self.all_modes_complete = False
	def _retrieve_current_values(self,list_by_mode,mode_idx,current_mode,param_num):
		# retrieves the appropriate list from a list of parameter lists
			# by mode
		output_list = list_by_mode[mode_idx]
		# ensure length of outputs is same as length of parameter_list
		output_list_length = len(output_list)
		if output_list_length == 1:
			output_list_trimmed = output_list*param_num
		elif output_list_length > param_num:
			print('Warning! More elements in list than parameters in mode ' + \
				current_mode)
			print(output_list)
			print('Replacing all values with NAs')
			output_list_trimmed = [numpy.NaN]*param_num
		elif output_list_length < param_num:
			print('Warning! Fewer elements in list than parameters in mode ' + \
				current_mode)
			print(output_list)
			print('Replacing all values with NAs')
			output_list_trimmed = [numpy.NaN]*param_num
		else:
			output_list_trimmed = output_list
		return(output_list_trimmed)
	def _id_parameters_to_loop_over(self):
		# identify which parameters need to be looped through in MLE
			# i.e. fitted parameters that MLE needs to be performed on
		# include 'unfixed' parameter, in which case no parameter is fixed
		parameters_to_loop_over_bool = \
			numpy.invert(self.current_permafixed_parameter_bool)* \
			numpy.invert((self.current_profile_point_list < 1))
		non_permafixed_parameters = [item for (item,bool_val) in \
			zip(self.current_parameter_list,parameters_to_loop_over_bool) \
			if bool_val]
		self.current_parameters_to_loop_over = ['unfixed'] + non_permafixed_parameters
		# include 1 profile point for 'unfixed' setting
			# (i.e. only perform 'unfixed' MLE once, since no need to
			# get likelihood profile)
		self.point_numbers_to_loop_over = numpy.append([1], \
			self.current_profile_point_list[ \
			parameters_to_loop_over_bool])
	def set_mode(self,mode_name,output_identifier):
		# for all MLE_parameter attributes, retrieves the parameter or
			# list of parameters corresponding to the current mode
		self.current_mode = mode_name
		mode_idx = self.mode_list.index(mode_name)
		self.current_sim_repeats = self.sim_repeats_by_mode[mode_idx]
		self.current_ms_positions = self.ms_positions[mode_idx]
		self.current_ms_grid_parameters = self.multistart_grid_parameters[mode_idx]
		self.current_logspace_profile_list = self.logspace_profile_list_by_mode[mode_idx]
		self.current_parameter_list = self.parameters_by_mode[mode_idx]
		# find the total number of parameters, including fixed ones
		self.total_param_num = len(self.current_parameter_list)
		# create lists, of length total_param_num, of settings for each
			# parameter in current_parameter_list
		# if input in file was incorrect length relative to number of
			# parameters, corresponding list is just NaN
		self.current_min_parameter_val_list = \
			numpy.array(self._retrieve_current_values(self.min_parameter_vals_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_max_parameter_val_list = \
			numpy.array(self._retrieve_current_values(self.max_parameter_vals_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_start_parameter_val_list = \
			numpy.array(self._retrieve_current_values(self.starting_parameter_vals_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_profile_point_list = \
			numpy.array(self._retrieve_current_values(self.profile_points_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_profile_lower_limit_list = \
			numpy.array(self._retrieve_current_values(self.profile_lower_limits_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_profile_upper_limit_list = \
			numpy.array(self._retrieve_current_values(self.profile_upper_limits_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		# identify list of parameters that are permanently fixed
		self.current_permafixed_parameter_bool = \
			self.current_max_parameter_val_list == self.current_min_parameter_val_list
		# identify parameters MLE needs to be performed on
		self._id_parameters_to_loop_over()
		# set up completefile tracker for these parameters
		self.parameter_completeness_tracker = \
			CompletenessTracker(self.current_parameters_to_loop_over)
		self.current_mode_complete = False
		# output_identifier is a string that will be included in filenames
		self.output_identifier = output_identifier
	def get_fitted_parameter_list(self):
		# get list of parameters to loop over for current mode
		return(self.current_parameters_to_loop_over)
	def set_parameter(self,parameter_name):
		# set current parameter, number of likelihood profile points
			# for it, and create a temporary list of fixed parameters
			# that includes it
		self.current_fixed_parameter = parameter_name
		if self.current_fixed_parameter == 'unfixed':
			self.current_profile_point_num = 1
			# list of fixed parameters is unchanged from default
			self.current_tempfixed_parameter_bool = \
				copy.copy(self.current_permafixed_parameter_bool)
			# other properties need to be NaN
		else:
			# find index of current_fixed_parameter in parameter list
			current_fixed_parameter_idx = self.current_parameter_list.index(parameter_name)
			# temporarily fix current parameter
			self.current_tempfixed_parameter_bool = \
				copy.copy(self.current_permafixed_parameter_bool)
			self.current_tempfixed_parameter_bool[current_fixed_parameter_idx] = \
				True
			self.current_profile_point_num = \
				self.current_profile_point_list[current_fixed_parameter_idx]
		self.output_id_parameter = self.output_identifier + '_' + self.current_fixed_parameter
		# identify how many dimensions are being used in multistart
		self.current_ms_grid_dimensions = sum([x == self.current_fixed_parameter \
			for x in self.current_ms_grid_parameters])
		self.current_parallel_processors = min(self.parallel_processors,
			self.current_ms_positions**self.current_ms_grid_dimensions)
	def update_parameter_completeness(self, completefile):
		# checks whether jobs for current parameter are all complete
		self.parameter_completeness_tracker.update_key_status( \
			self.current_fixed_parameter,completefile)
	def check_completeness_within_mode(self):
		# checks whether all parameters within mode are complete
		# change mode completeness status accordingly
		self.parameter_completeness_tracker.check_completeness()
		self.current_mode_complete = \
			self.parameter_completeness_tracker.get_completeness()
		self.mode_completeness_tracker.switch_key_completeness( \
			self.current_mode, self.current_mode_complete)
		return(self.current_mode_complete)
	def check_completeness_across_modes(self):
		self.mode_completeness_tracker.check_completeness()
		self.all_modes_complete = self.mode_completeness_tracker.get_completeness()
		return(self.all_modes_complete)

class SubmissionStringProcessor(object):
	# creates a code submission string appropriate to the programming
		# environment (module) being used
	def __init__(self, module, key_list, value_list, code_name):
		self._submission_string_generator(key_list,value_list,module,code_name)
	def get_code_run_input(self):
		return(self.code_sub_input_processor)
	def _submission_string_generator(self, key_list, value_list, module, \
		code_name):
		if module == 'matlab':
			code_sub_input_processor = \
			Cluster_Functions.MatlabInputProcessor(code_name)
		elif module == 'r':
			code_sub_input_processor = \
			Cluster_Functions.RInputProcessor(code_name)
		else:
			raise ValueError('%s is not an available module at the moment' \
				% (module))
		converted_key_string = \
			code_sub_input_processor.convert_mixed_list(key_list)
		converted_value_string = \
			code_sub_input_processor.convert_mixed_list(value_list)
		code_input_arg_list = [converted_key_string, converted_value_string]
		code_sub_input_processor.set_code_run_argument_string(code_input_arg_list)
		self.code_sub_input_processor = code_sub_input_processor

class MLEstimation(object):
	def __init__(self, mle_parameters, cluster_parameters, cluster_folders, \
		mle_folders, additional_code_run_keys, additional_code_run_values):
		self.mle_parameters = copy.deepcopy(mle_parameters)
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
		self.cluster_folders = copy.deepcopy(cluster_folders)
		self.mle_folders = copy.deepcopy(mle_folders)
		self.completefile = os.path.join(cluster_folders.completefile_path, \
			'_'.join(['MLE',mle_parameters.output_id_parameter,'completefile.txt']))
		self.job_name = '-'.join([mle_folders.experiment_folder_name,'MLE', \
			mle_parameters.output_id_parameter])
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self.output_filename = '_'.join(['data',mle_parameters.output_id_parameter])
		self.output_extension = 'csv'
		self.output_path = os.path.join(mle_folders.current_output_subfolder,'csv_output')
		self.module = 'matlab'
		self.code_name = '_'.join(['MLE',mle_parameters.current_mode])
		self.additional_beginning_lines_in_sbatch = []
		self.additional_end_lines_in_sbatch = []
			# don't include lines specific to matlab parallelization here
		self._create_code_run_input()
	def get_completefile_path(self):
		return(self.completefile)
	def _create_code_run_input(self):
		key_list = ['external_counter','combined_fixed_parameter_array', \
			'combined_min_array','combined_max_array','combined_length_array', \
			'combined_position_array','combined_start_values_array', \
			'parameter_list','csv_output_prename','output_folder', \
			'parallel_processors','ms_positions','combined_profile_ub_array', \
			'combined_profile_lb_array','ms_grid_parameter_array', \
			'combined_logspace_parameters']
		value_list = ['${' + self.cluster_parameters.within_batch_counter + '}', \
			self.mle_parameters.current_tempfixed_parameter_bool, \
			self.mle_parameters.current_min_parameter_val_list, \
			self.mle_parameters.current_max_parameter_val_list, \
			self.mle_parameters.current_profile_point_list, \
			['${' + self.cluster_parameters.within_batch_counter + '}'], \
				# if combined_position_array has length=1, MLE programs
					# interpret it as an array of the correct length
					# with the same value repeated
			self.mle_parameters.current_start_parameter_val_list, \
			self.mle_parameters.current_parameter_list, \
			self.output_filename, self.output_path, \
			self.mle_parameters.current_parallel_processors, \
			self.mle_parameters.current_ms_positions, \
			self.mle_parameters.current_profile_upper_limit_list, \
			self.mle_parameters.current_profile_lower_limit_list, \
			self.mle_parameters.current_ms_grid_parameters, \
			self.mle_parameters.current_logspace_profile_list]
		# take values from self.additional_code_run_keys and
			# self.additional_code_run_values where
			# self.additional_code_run_keys isn't already in key_list,
			# and add remaining values to key_list and value_list
		for current_key, current_val in \
			zip(self.additional_code_run_keys,self.additional_code_run_values):
			if not current_key in key_list:
				key_list.append(current_key)
				value_list.append(current_val)
		# process key_list and value_list into a submission string
		submission_string_processor = \
			SubmissionStringProcessor(self.module, key_list, value_list, \
				self.code_name)
		self.code_run_input = submission_string_processor.get_code_run_input()
	def run_job_submission(self):
		# handles submission of the job
		job_name = self.job_name
		job_numbers = [x + 1 for x in \
			list(range(self.mle_parameters.current_profile_point_num))]
		initial_time = self.cluster_parameters.current_time
		initial_mem = self.cluster_parameters.current_mem * \
			self.mle_parameters.current_parallel_processors
		cluster_parameters = self.cluster_parameters
		output_folder = self.output_path
		output_extension = self.output_extension
		output_filename = self.output_filename
		cluster_job_submission_folder = self.cluster_folders.cluster_job_submission_path
		experiment_folder = self.mle_folders.experiment_path
		module = self.module
		code_run_input = self.code_run_input
		additional_beginning_lines_in_sbatch = self.additional_beginning_lines_in_sbatch
		additional_end_lines_in_sbatch = self.additional_end_lines_in_sbatch
		parallel_processors = self.mle_parameters.current_parallel_processors
		completefile_path = self.completefile
		# set up and run batch jobs
		Cluster_Functions.job_flow_handler(job_name, job_numbers, initial_time, \
			initial_mem, cluster_parameters, output_folder, output_extension, \
			output_filename, cluster_job_submission_folder, experiment_folder, \
			module, code_run_input, additional_beginning_lines_in_sbatch, \
			additional_end_lines_in_sbatch, parallel_processors, completefile_path)

class CompletenessTracker(object):
	# Keeps a running tally of keys (which can be parameters, modes, etc) and checks their completeness
	def __init__(self,key_list):
		self._create_completeness_dict(key_list)
		self.completeness_status = False
	def _create_completeness_dict(self,key_list):
		# creates dictionary with parameter names as keys and False as values
		completeness_list = [False]*len(key_list)
		self.completeness_dict = dict(zip(key_list,completeness_list))
	def update_key_status(self, key, completefile_path):
		# checks whether key has associated completefile and changes its status accordingly
		if os.path.isfile(completefile_path):
			self.completeness_dict[key] = True
	def switch_key_completeness(self, key, value_bool):
		# switches a key value to value_bool
		# useful for cases when completefiles not kept track of
		self.completeness_dict[key] = value_bool
	def check_completeness(self):
		# checks whethere all parameters completed
		if all(self.completeness_dict.values()):
			self.completeness_status = True
		else:
			self.completeness_status = False
	def get_completeness(self):
		# returns completeness status
		return(self.completeness_status)

########################################################################

def loop_over_modes(mle_parameters, cluster_parameters, cluster_folders, \
	mle_folders, experiment_path):
	# Handles all of MLE across modes, including confidence
		# interval identification
	# May be too rigid, so should potentially be moved to main code
	for current_mode in mle_parameters.mode_list:
		##### RUN MLE #####
		output_id_string = current_mode
		mle_parameters.set_mode(current_mode,output_id_string)
		MLE_summary_file_path = os.path.join(experiment_path, \
			'_'.join([output_id_string,'MLE_file.csv']))
		# use additional_code_run_keys and values to specify where input
			# data comes from (and any other extra information that
			# doesn't come from setup_file)
		additional_code_run_keys = []
		additional_code_run_values = []
		# run MLE for current set of parameters
		run_MLE(mle_parameters, cluster_parameters, cluster_folders, mle_folders, \
			additional_code_run_keys, additional_code_run_values)
		# if all parameters for this mode are complete, update mode completeness
		# this also updates completeness across modes
		current_mode_mle_complete_status = \
			mle_parameters.check_completeness_within_mode()
	#	if current_mode_complete_status:
			##### RUN ASYMPTOTIC CI IDENTIFICATION #####

			# if asymptotic CI identification is complete:
			#	- identify position at which sims need to happen
			#		- run sims
			#			- get CIs from sims

def run_MLE(mle_parameters, cluster_parameters, cluster_folders, mle_folders, \
	additional_code_run_keys, additional_code_run_values):
	# Handles all of MLE for a particular mode, including confidence
		# interval identification
	# Loops through parameters for particular mode, finding ML
		# parameters at every fixed parameter value
	# Compiles results together to find asymptotic CI values (based on
		# chi-sq distribution of 2*log(LR))
	# Runs through simulations to find simulation-based CI values
	# mle_parameters must have the mode already set
	parameters_to_loop_over = mle_parameters.get_fitted_parameter_list()
	for current_fixed_parameter in parameters_to_loop_over:
		# set current parameter
		mle_parameters.set_parameter(current_fixed_parameter)
		# create MLEstimation object
		ml_estimator = MLEstimation(mle_parameters, cluster_parameters, \
			cluster_folders, mle_folders, additional_code_run_keys, \
			additional_code_run_values)
		# submit and track current set of jobs
		ml_estimator.run_job_submission()
		# track completeness within current mode
		mle_completefile = ml_estimator.get_completefile_path()
		mle_parameters.update_parameter_completeness(mle_completefile)









