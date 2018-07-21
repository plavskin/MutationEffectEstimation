#!/usr/bin/python

# Contains objects needed for running Maximum Likelihood Estimation

import os
import numpy
import Cluster_Functions
import copy
import pandas
from scipy.stats import chi2
import csv

class FolderManager(object):
	def __init__(self, cluster_parameters, cluster_folders, \
		experiment_folder_name):
		self.path_dict = {}
		self.path_dict['experiment_folder_name'] = experiment_folder_name
		self.path_dict['experiment_path'] = \
			os.path.join(cluster_parameters.composite_data_path, \
				experiment_folder_name)
		self.path_dict['sim_output_path'] = \
			os.path.join(cluster_parameters.temp_storage_path, \
				experiment_folder_name, 'simulated_phenotypes')
		self.path_dict['MLE_output_path'] = \
			os.path.join(cluster_parameters.temp_storage_path, \
				experiment_folder_name,'MLE_output')
		self.path_dict['completefile_folder'] = \
			cluster_folders.get_path('completefile_path')
		self.path_dict['LL_profile_path'] = \
			os.path.join(self.path_dict['experiment_path'],'LL_profiles')
		self.path_dict['CI_bound_path'] = \
			os.path.join(self.path_dict['experiment_path'],'CI_bounds')
	#	self.epilogue_path = os.path.join(cluster_parameters.home_path,cluster_parameters.username,'mut_effect_epilogue_files',experiment_folder_name)
		self.path_dict['MLE_combined_sim_path'] = \
			os.path.join(cluster_parameters.temp_storage_path, \
				'MLE_sim_outputs')
		self._set_up_folders()
	def _set_up_folders(self):
		setup_complete_file = \
			os.path.join(self.path_dict['completefile_folder'], \
			'mle_folder_setup_complete.txt')
		if not os.path.isfile(setup_complete_file):
			for current_folder_key, current_path in self.path_dict.iteritems():
				if not os.path.isdir(current_path):
					os.makedirs(current_path)
			open(setup_complete_file,'a').close()
	def get_path(self, folder_name):
		return(self.path_dict[folder_name])
	def set_current_output_subfolder(self,current_subfolder):
		# set (and if necessary, create) a subfolder to write temp output to
		if current_subfolder:
			self.path_dict['current_output_subfolder'] = \
				os.path.join(self.path_dict['MLE_output_path'], \
					current_subfolder)
			if not os.path.isdir(self.path_dict['current_output_subfolder']):
				os.makedirs(self.path_dict['current_output_subfolder'])
		else:
			self.path_dict['current_output_subfolder'] = self.path_dict['MLE_output_path']


class MLEParameters(object):
	############### ??? TO DO ??? ###############
	# Check that each mode or parameter is in the list once?
		# (currently only pays attention to the first time mode or parameter listed)
	############### ??? TO DO ??? ###############
	def __init__(self,parameter_list):
		self.runtime_percentile = parameter_list["runtime_percentile"]
		self.CI_pval_by_mode = parameter_list['CI_pval_by_mode']
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
		self.current_CI_pval = self.CI_pval_by_mode[mode_idx]
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
	def get_fitted_parameter_list(self, include_unfixed):
		# get list of parameters to loop over for current mode
		parameters_to_return = copy.copy(self.current_parameters_to_loop_over)
		if not include_unfixed:
			parameters_to_return.remove('unfixed')
		return(parameters_to_return)
	def get_complete_parameter_list(self):
		# get list of parameters in current mode
		return(self.current_parameter_list)
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
			self.current_fixed_parameter_idx = self.current_parameter_list.index(parameter_name)
			# temporarily fix current parameter
			self.current_tempfixed_parameter_bool = \
				copy.copy(self.current_permafixed_parameter_bool)
			self.current_tempfixed_parameter_bool[self.current_fixed_parameter_idx] = \
				True
			self.current_profile_point_num = \
				self.current_profile_point_list[self.current_fixed_parameter_idx]
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
		self.current_mode_complete = \
			self.parameter_completeness_tracker.get_completeness()
		self.mode_completeness_tracker.switch_key_completeness( \
			self.current_mode, self.current_mode_complete)
		return(self.current_mode_complete)
	def check_completeness_across_modes(self):
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
		self.completefile = \
			os.path.join(cluster_folders.get_path('completefile_path'), \
				'_'.join(['MLE',mle_parameters.output_id_parameter, \
					'completefile.txt']))
		self.job_name = '-'.join([mle_folders.get_path('experiment_folder_name'),'MLE', \
			mle_parameters.output_id_parameter])
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self.within_batch_counter_call = \
			'${' + self.cluster_parameters.within_batch_counter + '}'
		self.output_path = mle_folders.get_path('current_output_subfolder')
		self.output_filename = _generate_filename(self.output_path, \
			self.within_batch_counter_call, mle_parameters.output_identifier, \
			mle_parameters.current_fixed_parameter, 'data')
		self.output_file_label = _generate_file_label('data', \
			mle_parameters.output_identifier, \
			mle_parameters.current_fixed_parameter)
		self.module = 'matlab'
		self.output_extension = 'csv'
		self.code_name = '_'.join(['MLE',mle_parameters.current_mode])
		self.additional_beginning_lines_in_job_sub = []
		self.additional_end_lines_in_job_sub = []
			# don't include lines specific to matlab parallelization here
		self._create_code_run_input()
	def get_output_path(self):
		return(self.output_path)
	def get_completefile_path(self):
		return(self.completefile)
	def _create_code_run_input(self):
		key_list = ['external_counter','combined_fixed_parameter_array', \
			'combined_min_array','combined_max_array','combined_length_array', \
			'combined_position_array','combined_start_values_array', \
			'parameter_list','output_file', \
			'parallel_processors','ms_positions','combined_profile_ub_array', \
			'combined_profile_lb_array','ms_grid_parameter_array', \
			'combined_logspace_parameters']
		value_list = [self.within_batch_counter_call, \
			self.mle_parameters.current_tempfixed_parameter_bool, \
			self.mle_parameters.current_min_parameter_val_list, \
			self.mle_parameters.current_max_parameter_val_list, \
			self.mle_parameters.current_profile_point_list, \
			[self.within_batch_counter_call], \
				# if combined_position_array has length=1, MLE programs
					# interpret it as an array of the correct length
					# with the same value repeated
			self.mle_parameters.current_start_parameter_val_list, \
			self.mle_parameters.current_parameter_list, \
			self.output_filename, \
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
		output_file_label = self.output_file_label
		cluster_job_submission_folder = \
			self.cluster_folders.get_path('cluster_job_submission_path')
		experiment_folder = self.mle_folders.get_path('experiment_path')
		module = self.module
		code_run_input = self.code_run_input
		additional_beginning_lines_in_job_sub = self.additional_beginning_lines_in_job_sub
		additional_end_lines_in_job_sub = self.additional_end_lines_in_job_sub
		parallel_processors = self.mle_parameters.current_parallel_processors
		completefile_path = self.completefile
		# set up and run batch jobs
		Cluster_Functions.job_flow_handler(job_name, job_numbers, initial_time, \
			initial_mem, cluster_parameters, output_folder, output_extension, \
			output_file_label, cluster_job_submission_folder, experiment_folder, \
			module, code_run_input, additional_beginning_lines_in_job_sub, \
			additional_end_lines_in_job_sub, parallel_processors, completefile_path)

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
	def get_key_completeness(self, key):
		return(self.completeness_dict[key])
	def _check_completeness(self):
		# checks whethere all parameters completed
		if all(self.completeness_dict.values()):
			self.completeness_status = True
		else:
			self.completeness_status = False
	def get_completeness(self):
		# checks and returns completeness status
		self._check_completeness()
		return(self.completeness_status)

class LLWarning(object):
	# stores warnings for LLProfile objects
	def __init__(self):
		self.non_monotonic = False
	def set_non_monotonic(self):
		self.non_monotonic = True
	def get_warning_line(self):
		warning_list = []
		if self.non_monotonic:
			warning_list.append('non-monotonic LL profile before or after ML')
		warning_string = ';'.join(warning_list)
#		if warning_list:
#			warning_string = 'warning! ' + warning_string
		return(warning_string)

class LLProfile(object):
	# Gets, holds, and updates log likelihood profile
	def __init__(self, mle_parameters, datafile_path, LL_profile_folder, additional_param_df):
#		unfixed_ll_param_df, true_max_param_df
		self.warnings = LLWarning()
		self.profile_points = mle_parameters.current_profile_point_num
		self.output_identifier = mle_parameters.output_identifier
		self.current_fixed_parameter = mle_parameters.current_fixed_parameter
		self.datafile_path = datafile_path
		self.fixed_param = mle_parameters.current_fixed_parameter
#		self.fixed_param_idx = mle_parameters.current_fixed_parameter_idx
#		self.parameter_list = mle_parameters.current_parameter_list
			# the percentile of runtimes returned by LLprofile
#		self.unfixed_ll_param_df = unfixed_ll_param_df
#		self.true_max_param_df = true_max_param_df
		self.additional_param_df = copy.copy(additional_param_df)
		self.max_LL = None
		self.LL_file = os.path.join(LL_profile_folder, \
			('_'.join(['LL_file', self.output_identifier, \
				self.current_fixed_parameter]) + '.csv'))
		self.LL_df = pandas.DataFrame()
#		self.non_parameter_columns = ['LL','runtime_in_secs']
#		self.non_parameter_column_num = len(self.non_parameter_columns)
#		self.LL_matrix = numpy.array([])
#		self.LL_df = pandas.DataFrame(\
#			columns = (self.non_parameter_columns + self.parameter_list))
	def _set_ML_params(self, ml_param_df):
		self.ML_params = copy.copy(ml_param_df)
		self.max_LL = ml_param_df.iloc[0]['LL']
		self.fixed_param_MLE_val = ml_param_df.iloc[0][self.fixed_param]
#	def _check_and_update_ML(self, ml_param_df):
		# checks whether ml_param_df contains a higher LL value
			# that current self.max_LL; if so, update self.ML_params and
			# self.max_LL, and add a warning to this LL profile
#		current_LL = ml_param_df['LL'][0]
#		if current_LL > self.max_LL:
#			self._set_ML_params(ml_param_df)
#			self.warnings.set_non_unfixed_ML(fixed_param)
	def _add_vals(self, ll_param_df):
		# adds values to self.LL_df from current_param_datafile
		# if current LL is max, updates max_LL, ML_params, and
			# fixed_param_MLE_val
		if not self.LL_df.empty:
			# concatenate data frames, dropping any duplicate rows and
				# resetting indices
			self.LL_df = self.LL_df.append(ll_param_df).drop_duplicates().reset_index(drop=True)
		else:
			self.LL_df = ll_param_df
#		if len(self.LL_matrix) == 0:
#			self.LL_matrix = ll_param_df
#		else:
#			self.LL_matrix = vstack(self.LL_matrix,ll_param_df)
#	def _populate_LL_matrix(self):
	def _id_max_LL(self):
		# identifies and sets the parameter values corresponding to the
			# max likelihood
		# id row of LL_df corresponding to max LL
#		ml_params = pandas.DataFrame(self.LL_df.loc[self.LL_df['LL'].idxmax()]).transpose()
		ml_param_df = self.LL_df.iloc[[self.LL_df['LL'].idxmax()]]
		### ml_param_df is a series, needs to be a df
			# idxmax OK here because indices are reset during appending
		# set this row to be the ML parameters
		self._set_ML_params(ml_param_df)
	def _populate_LL_df(self):
		# fill in data in LL_df
		for current_pp in range(1,(self.profile_points+1)):
			current_datafile = _generate_filename(self.datafile_path, \
				str(current_pp), self.output_identifier, \
				self.current_fixed_parameter, 'data')
			if os.path.isfile(current_datafile):
				current_data_df = _get_MLE_params(current_datafile)
				# add current_data_np to LL array and update max likelihood value
				self._add_vals(current_data_df)
	def _sort_by_profiled_param(self):
		# returns order of rows in self.LL_df sorted by the
			# parameter being profiled
		self.LL_df_sorted = self.LL_df.sort_values(self.fixed_param)
#		parameter_indices_sorted = self.LL_matrix[:, \
#			self.fixed_param_idx + self.non_parameter_column_num].argsort()
#		return(numpy.array(parameter_indices_sorted))
#	def get_LL_profile(self):
#		# gets a numpy array in which the first column is the parameter
#			# being profiled and the second column are the
#			# correspinding LL values, sorted by the profiled parameter
#		sorted_indices = self._sort_by_profiled_param()
#		LL_profile = self.LL_matrix[sorted_indices[:,numpy.newaxis], \
#			[self.fixed_param_idx + self.non_parameter_column_num,0]]
#		return(LL_profile)
	def _check_monotonicity(self):
		### PROBABLY MOVE THIS TO CI SECTION TO AVOID REDUNDANCY
		# check whether LL profile is monotonic before and after ML parameter value
		# In simple and accurately estimated LL landscape, LL
			# profile expected to increase monotonically up until
			# max LL, then decrease monotonically after; if this
			# isn't the case, throw a warning
#		LL_profile = self.get_LL_profile()
#		y_vals = LL_profile[:,1]
#		x_vals = LL_profile[:,0]
		x_vals = self.LL_df_sorted[self.fixed_param]
		y_vals = self.LL_df_sorted['LL']
		y_diffs = numpy.diff(y_vals)
		increasing_LL_section = y_diffs[x_vals[:-1] < self.fixed_param_MLE_val]
		decreasing_LL_section = y_diffs[x_vals[:-1] >= self.fixed_param_MLE_val]
			# remember indexing is different for y_diffs and y_vals,
				# it's correct here
		monotonicity_state = numpy.all(increasing_LL_section >= 0) and \
			numpy.all(decreasing_LL_section <= 0)
		if not monotonicity_state:
			self.warnings.set_non_monotonic()
	def _write_LL_df(self):
#		LL_array_sorted = \
#			self.LL_matrix[self.LL_matrix[:,\
#				self.fixed_param_idx + self.non_parameter_column_num].argsort()]
#		LL_df = pandas.DataFrame(
#			columns = (self.non_parameter_columns + self.parameter_list),
#			data = LL_array_sorted)
#		if not os.path.isfile(self.LL_file):
		# write LL_df to file regardless of whether file already
			# exists, since LL_df may have been updated
		self.LL_df_sorted.to_csv(path_or_buf=self.LL_file,index=False)
#	def _read_LL_matrix(self):
	def _read_LL_df(self):
		# reads self.LL_df from a pre-recorded file
		self.LL_df = pandas.read_csv(filepath_or_buffer=self.LL_file)
		self._sort_by_profiled_param()
#		self.LL_matrix = LL_df.as_matrix()
	def _set_LL_df(self):
		# gets LL_df, either from a pre-recorded file, or from MLE outputs
		# if true_max_param_df is non-empty, add it to LL_df and update
			# max LL-related parameters
		if os.path.isfile(self.LL_file):
#			self._read_LL_matrix()
			self._read_LL_df()
		else:
			self._populate_LL_df()
		if not self.additional_param_df.empty:
			self._add_vals(self.additional_param_df)
		self._sort_by_profiled_param()
	def run_LL_profile_compilation(self):
		# compiles and writes LL_profile;
			# identifies and sets max LL parameters;
			# checks monotonicity
		self._set_LL_df()
		self._write_LL_df()
		self._id_max_LL()
		self._check_monotonicity()
	def get_LL_df(self):
		# returns sorted LL_df
		return(self.LL_df_sorted)
	def get_fixed_param_MLE_val(self):
		# returns MLE value of current fixed parameter
		return(self.fixed_param_MLE_val)
	def get_max_LL(self):
		return(self.max_LL)
	def get_ML_params(self):
		return(self.ML_params)
	def get_warnings(self):
		warning_line = self.warnings.get_warning_line()
		return(warning_line)
	def get_runtime(self, runtime_display_percentile):
		# finds the runtime, in hours, corresponding to
			# self.runtime_display_percentile
		percentile_as_decimal = runtime_display_percentile/100
		time_quantile_in_seconds = \
			self.LL_df_sorted.runtime_in_secs.quantile(percentile_as_decimal)
		time_quantile_in_hours = time_quantile_in_seconds/3600
		return(time_quantile_in_hours)

class CIWarning(object):
	# stores warnings for OneSidedCIBound objects
	def __init__(self, CI_side, CI_type):
		self.points_to_create_CI = 0
		self.all_points_within_CI_bound = False
		self.CI_side = CI_side
		self.CI_type = CI_type
	def set_all_points_within_CI_bound(self):
		self.all_points_within_CI_bound = True
	def set_points_to_create_CI(self,points_to_create_CI):
		self.points_to_create_CI = points_to_create_CI
	def set_no_output_file(self):
		self.no_output_file = True
	def get_warning_line(self):
		warning_list = []
		if self.all_points_within_CI_bound:
			warning_list.append('All points used to calculate ' + self.CI_side + \
				 ' CI for ' + self.CI_type + \
				 ' CI are between expected CI bound and MLE.')
		if self.no_output_file:
			warning_list.append('Although the curve-fitting job was completed, there was no output file created for ' + \
				self.CI_side + 'CI for ' + self.CI_type)
		if set_points_to_create_CI == 0:
			warning_list.append('It seems '  + self.CI_side + ' CI for ' + \
				self.CI_type + ' CI was not set.')
		elif set_points_to_create_CI == 1:
			warning_list.append('No profile points besides MLE on ' + \
				self.CI_side + ' side of CI for ' + self.CI_type + ' CI.')
		elif set_points_to_create_CI == 2:
			warning_list.append('Only one profile point besides MLE on ' + \
				self.CI_side + \
				' side of CI for ' + self.CI_type + \
				'CI; CI boundary calculated as linear interpolation of p-val between to two points')
		warning_string = ';'.join(warning_list)
		return(warning_string)

class OneSidedCIBound(object):
	# Stores data for CI bound on one side of MLE
	def __init__(self, p_val, profile_side, LL_df, df, fixed_param_MLE_val,
		fixed_param, max_LL, CI_type, mle_folders, cluster_parameters, cluster_folders):
		self.cdf_bound = 1-p_val
		self.points_to_fit_curve = 3
		self.profile_side = profile_side
		self.df = df
		self.CI_type = CI_type
			# CI_type can be either 'asymptotic' or 'sim'
		self.fixed_param = fixed_param
		self.fixed_param_MLE_val = fixed_param_MLE_val
		self.max_LL = max_LL
		self.cluster_parameters = cluster_parameters
		self.cluster_folders = cluster_folders
		self.mle_folders = mle_folders
		self.module = 'matlab'
		self.CI_bound_set = False
		self.output_prename = '-'.join(['CI_bound', self.profile_side, \
			self.CI_type])
		self.CI_bound_name = _generate_file_label(self.output_prename, \
			self.output_identifier, self.fixed_param)
		self.CI_bound_output_file = \
			_generate_filename(self.mle_folders.get_path('CI_bound_path'), '1', \
				self.output_identifier, self.fixed_param, self.output_prename)
		self.CI_bound_fit_file = \
			_generate_filename(self.mle_folders.get_path('CI_bound_path'), '1', \
				self.output_identifier, self.fixed_param, \
				(self.output_prename + '_fit_file'))
		self.completefile = os.path.join(cluster_folders.get_path('completefile_path'), \
			'_'.join([self.CI_bound_name,'completefile.txt']))
		self.additional_beginning_lines_in_job_sub = []
		self.additional_end_lines_in_job_sub = []
		self.additional_code_run_keys = []
		self.additional_code_run_values = []
		self._select_profile_side(LL_df)
		self.warning = CIWarning(profile_side)
		# need to run _asymptotic_p_val_calc and
			# _find_CI_proximal_LL_points even if curve fitting to LL
			# points has already occurred, since these functions throw
			# important warnings
		self._asymptotic_p_val_calc()
		self._find_CI_proximal_LL_points()
	def _select_profile_side(self, LL_df):
		# selects data from correct side of likelihood profile
		if self.profile_side == 'lower':
			LL_df_one_side = LL_df[(LL_df[self.fixed_param] <= self.fixed_param_MLE_val)]
			self.default_CI_bound = float("-inf")
		elif self.profile_side == 'upper':
			LL_df_one_side = LL_df[(LL_df[self.fixed_param] >= self.fixed_param_MLE_val)]
			self.default_CI_bound = float("inf")
		else:
			print(('Error! Underfined profile side ' + self.profile_side))
		self.one_sided_LL_df = LL_df_one_side.sort_values(self.fixed_param)
	def _asymptotic_p_val_calc(self):
		# calculates p-values for every point in LL_profile, following asymptotic assumption
		x_vals = self.one_sided_LL_df[self.fixed_param]
		y_vals = self.one_sided_LL_df['LL']
		D_vals = 2*(self.max_LL-y_vals)
		cdf_vals = 0.5 + 0.5*chi2.cdf(D_vals,self.df)
			# rather than calculating p-vals for the whole distribution, we
				# calculate a 'reflected' cdf for the lower side of the
				# distribution so that the same fitting algorithm could be
				# used to identify the x-value closest to the desired
				# p-value cutoff
		self.one_sided_LL_df['cdf_vals'] = cdf_vals
	def _id_proximal_points(self, target_y, y_vals, num_points):
		# returns indices of up to num_points points closest to target_y
		dist_to_target = numpy.absolute(target_y-y_vals)
		sorted_indices = numpy.argsort(dist_to_target)
		closest_indices = sorted_indices[0:num_points]
		return(closest_indices)
	def _set_CI_bound(self, CI_bound):
		# sets self.CI_bound and changed self.CI_bound_set to true
		self.CI_bound = CI_bound
		self.CI_bound_set = True
	def _find_CI_proximal_LL_points(self):
		# identify points that are most proximal to
			# conf int cutoff
		# if there are 3 such points, there's no problem
			# if there are 2 such points, run linear CI estimation and throw warning
			# if there is 1 or 0 such points, set CI bound to -/+Inf, throw warning
			# if all points within CI bound, throw warning
		number_profile_points = self.one_sided_LL_df.shape[0]
		if number_profile_points <= 1:
			self._set_CI_bound(self.default_CI_bound)
		elif number_profile_points == 2:
			self.code_name = 'Linear_Bound_Finder'
		elif number_profile_points > 2:
			self.code_name = 'Quadratic_Bound_Finder'
		self.warning.set_points_to_create_CI(number_profile_points)
		if numpy.max(profile_points) < self.cdf_bound:
			self.warning.set_all_points_within_CI_bound()
		# get indices of closest points to CI bound
		CI_bound_proximal_indices = \
			self._id_proximal_points(self.cdf_bound, self.one_sided_LL_df['LL'], \
				self.points_to_fit_curve)
		# create a new df with only points closest to CI bound
		self.CI_bound_proximal_points = self.one_sided_LL_df[CI_bound_proximal_indices]
	def _create_code_run_input(self):
		key_list = ['cdf_bound', \
			'mle_param_val', \
			'parameter_values',\
			'cdf_vals', \
			'output_file', \
			'fit_file']
		value_list = [self.cdf_bound, \
			self.fixed_param_MLE_val, \
			self.CI_bound_proximal_points[self.fixed_param], \
			self.CI_bound_proximal_points['cdf_vals'], \
			self.CI_bound_output_file, \
			self.CI_bound_fit_file]
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
	def _run_CI_finder_submission(self):
		# handles submission of the job
		self.job_name = \
			'-'.join([self.mle_folders.get_path('experiment_folder_name'), \
				self.CI_bound_name])
		job_numbers = [1]
		initial_time = 30
		initial_mem = 1024
		cluster_parameters = self.cluster_parameters
		output_folder = self.mle_folders.get_path('CI_bound_path')
		output_extension = '.csv'
		output_file_label = self.CI_bound_name
		cluster_job_submission_folder = \
			self.cluster_folders.get_path('cluster_job_submission_path')
		experiment_folder = self.mle_folders.get_path('experiment_path')
		module = self.module
		code_run_input = self.code_run_input
		additional_beginning_lines_in_job_sub = self.additional_beginning_lines_in_job_sub
		additional_end_lines_in_job_sub = self.additional_end_lines_in_job_sub
		parallel_processors = 1
		completefile_path = self.completefile
		# set up and run batch jobs
		Cluster_Functions.job_flow_handler(job_name, job_numbers, initial_time, \
			initial_mem, cluster_parameters, output_folder, output_extension, \
			output_file_label, cluster_job_submission_folder, experiment_folder, \
			module, code_run_input, additional_beginning_lines_in_job_sub, \
			additional_end_lines_in_job_sub, parallel_processors, completefile_path)
	def find_CI_bound(self):
		# check if CI bound has been IDed
		# if it has, read it in, assign it to self
		# if not, submit job
		# look for completefile
		if not self.CI_bound_set:
			if os.path.isfile(self.completefile):
				# if completefile exists, look for CI bound file
				# if CI bound file doesn't exist, set CI bound to
					# negative infinity, add error message
				# otherwsise read CI bound file
				if os.path.isfile(self.CI_bound_output_file):
					with open(self.CI_bound_output_file, 'rU') as \
						CI_bound_contents:
						self._set_CI_bound(list(csv.reader(CI_bound_contents))[0])
				else:
					self._set_CI_bound(self.default_CI_bound)
					self.warning.set_no_output_file()
			else:
				self._run_CI_finder_submission()
	def get_CI_bound(self):
		# returns CI_bound if it exists, otherwise returns None
		if self.CI_bound_set:
			return(self.CI_bound)
		else:
			return(None)
	def get_CI_bound_warning(self):
		# returns the warning string for current CI bound
		warning_string = self.warning.get_warning_line()
		return(warning_string)

class TwoSidedCI(object):
	# compiles two-sided CI
	def __init__(self, p_val, LL_df, deg_freedom, fixed_param_MLE_val, \
		fixed_param, max_LL, CI_type, mle_folders, cluster_parameters, \
		cluster_folders):
		self.CI_sides = ['lower','upper']
		self.CI_object_dictionary = dict()
		self.CI_dictionary = dict()
		self.CI_completeness_tracker = CompletenessTracker(self.CI_sides)
		self.CI_complete = False
		self.CI_warning_list = []
		for current_CI_side in self.CI_sides:
			self.CI_object_dictionary[current_CI_side] = OneSidedCIBound(p_val, \
				current_CI_side, LL_df, deg_freedom, fixed_param_MLE_val, \
				fixed_param, max_LL, CI_type, mle_folders, cluster_parameters, \
				cluster_folders)
	def find_CI(self):
		# identifies confidence interval
		for current_CI_side in self.CI_sides:
			self.CI_object_dictionary[current_CI_side].find_CI_bound()
			current_CI_bound = self.CI_object_dictionary[current_CI_side].get_CI_bound()
			if current_CI_bound:
				self.CI_dictionary[current_CI_side] = current_CI_bound
				self.CI_completeness_tracker.switch_key_completeness(current_CI_side, \
					True)
				self.CI_complete = \
					self.CI_completeness_tracker.get_completeness()
				current_CI_warning = \
					self.CI_object_dictionary[current_CI_side].get_CI_bound_warning()
				self.CI_warning_list.append(current_CI_warning)
	def get_CI(self):
		# returns confidence interval dictionary if it's complete,
			# otherwise, returns None
		if self.CI_complete:
			return(self.CI_dictionary)
		else:
			return(None)
	def get_CI_warning(self):
		# returns a combined warning line for both CI bounds
		combined_warning_line = ';'.join(self.CI_warning_list)
		return(combined_warning_line)

class SingleParamResultSummary(object):
	# stores MLE estimates and CIs for a single parameter
	def __init__(self, fixed_param_MLE, max_LL, warning_string, \
		CI_dict):
		self.content_dict = {'param_MLE': fixed_param_MLE, \
			'max_LL': max_LL}
		# Add contents of CI_dict to self.content_dict
		# CI_dict contains a dictionary of CIs; each key is the type of
			# CI calculation (e.g. 'asymptotic' or 'simulation-based');
			# each value is a dictionary of CI bounds, with the key
			# being the bound name and the value being the bound val
		for current_CI_type, current_CI_bound_dict in CI_dict.iteritems():
			self._set_CI(current_CI_bound_dict, current_CI_type)
		# add warnings
		self.content_dict['warnings'] = warning_string
	def _set_CI(self, CI_bound_dict, CI_type):
		# adds keys for each element of CI_bound_dict, combining key
			# name with CI_type
		for current_key, current_value in CI_bound_dict.iteritems():
			new_key = '_'.join([CI_type, 'CI_bound', current_key])
			self.content_dict[new_key] = current_value
	def get_contents(self):
		return(self.content_dict)

class CombinedResultWarning(object):
	# stores warnings for CombinedResultSummary objects
	def __init__(self):
		self.unfixed_file_missing = False
		self.unfixed_not_ML = False
	def set_unfixed_file_missing(self):
		self.unfixed_file_missing = True
	def set_non_unfixed_ML(self,fixed_param):
		self.unfixed_not_ML = True
		self.parameter_LL_containing_ML = fixed_param
	def get_non_unfixed_ML_status(self):
		return(self.unfixed_not_ML)
	def get_warning_line(self):
		warning_list = []
		if self.unfixed_file_missing:
			warning_list.append('no unfixed condition file found')
		if self.unfixed_not_ML:
			warning_list.append( \
				'unfixed condition did not find true ML! ML found in fixed ' + \
				self.parameter_LL_containing_ML + ' search')
		warning_string = ';'.join(warning_list)
		return(warning_string)


class CombinedResultSummary(object):
	# stores, and writes, summary of MLE results
	# Only to be run after initial MLE across profile points complete!
	# Loops through parameters within a mode, checks max_LL in each
	# Once all LL_profiles are complete, if max_LL in at least one
		# LL_profile is not unfixed_LL, throws warning, recalculates
		# LL_profiles with new max_LL
	# Writes LL_profiles, gets lower and upper asymptotic CI
	# If asymptotic CI identification is complete:
	#	- identify position at which sims need to happen
	#		- run sims
	#			- get CIs from sims
	# keep track of asymptotic and sim CIs using completeness tracker across parameters within mode
	# once CIs complete (either asymptotic only or asymptotic and sim, depending on settings),
	#	create summary file
	# keep track of summary files for each mode; once those are complete, stop running current folder
	def __init__(self, mle_folders, mle_parameters, cluster_parameters, \
		cluster_folders):
		self.mle_datafile_path = mle_folders.get_path('current_output_subfolder')
		self.mle_parameters = copy.deepcopy(mle_parameters)
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
		self.cluster_folders = copy.deepcopy(cluster_folders)
		self.LL_profile_folder = mle_folders.get_path('LL_profile_path')
		self.runtime_percentile = mle_parameters.runtime_percentile
		self.pval = mle_parameters.current_CI_pval
		self._create_combined_output_file(mle_folders)
		self.max_LL = None
		self.warnings = CombinedResultWarning()
		self.completeness_tracker = CompletenessTracker(['initialization', \
			'asymptotic_CIs', 'simulation-based_CIs'])
		self.runtime_quant_list = numpy.array([])
		self.true_max_param_df = pandas.DataFrame()
		self.combined_results_df = pandas.DataFrame()
		# set MLE results from 'unfixed' parameter (i.e. fitting all
			# unfixed params together)
		self.unfixed_mle_file = _generate_filename(self.mle_datafile_path, \
			'1', mle_parameters.output_identifier, 'unfixed', 'data')
		self._set_unfixed_param_data()
		self._check_and_update_ML(self.unfixed_ll_param_df,'unfixed')
	def _create_combined_output_file(self, mle_folders):
		# creates the name of the combined output file for the results
		experiment_path = mle_folders.get_path('experiment_path')
		self.combined_results_output_file = os.path.join(experiment_path, \
			('_'.join(['MLE_output', self.mle_parameters.output_identifier]) + \
				'.csv'))
	def _set_unfixed_param_data(self):
		# gets data from self.unfixed_mle_file and uses it to update
			# self.max_LL or, if file is not there, to create a warning
		if os.path.isfile(self.unfixed_mle_file):
			self.unfixed_ll_param_df = _get_MLE_params(self.unfixed_mle_file)
			self._check_and_update_ML(self.unfixed_ll_param_df,'unfixed')
		else:
			self.unfixed_ll_param_df = pandas.DataFrame()
			self.warnings.set_unfixed_file_missing()
	def _set_combined_ML_params(self, ml_param_df):
		self.true_max_param_df = ml_param_df
		self.max_LL = ml_param_df.iloc[0]['LL']
	def _check_and_update_ML(self, ml_param_df, fixed_param):
		# checks whether ml_param_np_array contains a higher LL value
			# that current self.max_LL; if so, update self.ML_params and
			# self.max_LL, and add a warning to these combined results
		current_LL = ml_param_df.iloc[0]['LL']
		if current_LL > self.max_LL:
			if self.max_LL:
				# only set warning about surpassing unfixed file max_LL
					# if self.max_LL is not None
				self.warnings.set_non_unfixed_ML(fixed_param)
			self._set_combined_ML_params(ml_param_df)
	def _update_LL_profile(self, non_profile_max_params, fixed_parameter):
		# updates LL_profile for current_fixed_parameter with known max
			# parameters (non_profile_max_params)
		# returns any warnings, runtime_quantile, and ML_params from
			# LL_profile
		# first set current parameter
		self.mle_parameters.set_parameter(fixed_parameter)
		# create an LLProfile for current parameter
		ll_profile = LLProfile(self.mle_parameters, \
			self.mle_datafile_path, self.LL_profile_folder, \
			non_profile_max_params)
		ll_profile.run_LL_profile_compilation()
		warning_line = ll_profile.get_warnings()
		ML_params = ll_profile.get_ML_params()
		runtime_quantile = ll_profile.get_runtime(self.runtime_percentile)
		return({'warnings':warning_line, 'runtime_quantile':runtime_quantile, \
			'ML_params':ML_params})
	def _create_initial_LL_profiles(self):
		# creates LL_profiles (if necessary) and gets their max parameters
		parameters_to_loop_over = self.mle_parameters.get_fitted_parameter_list(False)
		for current_fixed_parameter in parameters_to_loop_over:
			ll_profile_outputs = \
				self._update_LL_profile(self.unfixed_ll_param_df,
					current_fixed_parameter)
			# check if LL of current_ll_profile is higher than current max_LL, and update true_max_param_df accordingly
			current_ML_params = ll_profile_outputs['ML_params']
			self._check_and_update_ML(current_ML_params, current_fixed_parameter)
			# get runtime
			runtime_quantile = ll_profile_outputs['runtime_quantile']
			self.runtime_quant_list = numpy.append(self.runtime_quant_list,runtime_quantile)
	def _correct_LL_profiles(self):
		# updates LL_profiles
		parameters_to_loop_over = mle_parameters.get_fitted_parameter_list(False)
		ml_params_to_include = pandas.concat(self.unfixed_ll_param_df, \
			self.true_max_param_df)
		for current_fixed_parameter in parameters_to_loop_over:
			self._update_LL_profile(self.unfixed_ll_param_df,
				ml_params_to_include)
	def _get_runtime_CI(self):
		# creates 'confidence intervals' across parameters for the
			# self.runtime_percentile-th estimates within a parameter's
			# estimated profile points
		runtime_CI_bounds = {}
		runtime_CI_bounds['lower'] = \
			numpy.percentile(self.runtime_quant_list,self.pval/2*100)
		runtime_CI_bounds['upper'] = \
			numpy.percentile(self.runtime_quant_list,(1-self.pval/2)*100)
		runtime_CIs = {'asymptotic':{'lower':numpy.nan,'upper':numpy.nan}, \
			'simulation-based':runtime_CI_bounds}
		return(runtime_CIs)
	def _add_line_to_combined_df(self,fixed_param_name, fixed_param_MLE, \
		warning_string, CI_dict):
		# creates self.combined_summary or adds a line to it
		current_results_line = SingleParamResultSummary(fixed_param_MLE, \
			self.max_LL, warning_string, CI_dict)
		current_results_dict = current_results_line.get_contents()
		if self.combined_results_df.empty:
			# create dataframe using currently passed line
			self.combined_results_df = \
				pandas.DataFrame(current_results_dict, index = [fixed_param_name])
		else:
			# check that all keys in current_results_dict are columns
				# in self.combined_results_df; if not, create these
				# columns before updating data
			dict_keys = current_results_dict.keys()
			df_columns = list(self.combined_results_df.columns.values)
			new_columns = set(dict_keys).symmetric_difference(df_columns)
			self.combined_results_df = \
				self.combined_results_df.append(pandas.DataFrame(columns=new_columns), \
					sort = False)
			# insert current_results_dict into correct row in df
			current_results_series = pandas.Series(current_results_dict)
			self.combined_results_df.loc[fixed_param_name] = \
				current_results_series
	def _write_combined_df(self):
		# writes self.combined_results_df to file
		self.combined_results_df.to_csv(path_or_buf=self.combined_results_output_file, \
			index=True)
	def _read_combined_df(self):
		# reads self.combined_results_df from a pre-recorded file
		self.combined_results_df = \
			pandas.read_csv(filepath_or_buffer=self.combined_results_output_file, \
				index_col = 0)
		### READ in MLE, add LL, convert to true_max_param_df
	def _create_combined_df(self):
		# creates combined_results_df and adds line for runtime and general warnings
		# include a line containing the mean self.runtime_quantile-th
			# time within the mode, as well as the self.pval-based CI
			# on this time, and the general warning line that applies
			# the current mode
		time_string = '_'.join(['avg', (str(self.runtime_percentile) + 'th'), \
			'percentile', 'runtime', 'across', 'profile', 'points', 'in', 'hrs'])
		# get runtime mean and CIs
		mean_runtime_quant = numpy.mean(self.runtime_quant_list)
		runtime_CIs = self._get_runtime_CI()
		# get warning line for current mode
		mode_warning_line = self.warnings.get_warning_line()
		# make a line to write to combined_results_df
		self._add_line_to_combined_df(time_string, mean_runtime_quant, \
			mode_warning_line, runtime_CIs)
		# append a DataFrame where row names are parameter names
		mle_df = self.true_max_param_df.transpose()
		# rename column to 'param_MLE' and remove 'LL' row
		mle_df.columns=['param_MLE']
		mle_df.drop(index=['LL','runtime_in_secs'],inplace = True)
		self.combined_results_df = \
			self.combined_results_df.append(mle_df, sort = False)
	def initialize_combined_results(self):
		# check whether initial mle has been completed
		mle_complete = self.mle_parameters.check_completeness_within_mode()
		if mle_complete:
			# check whether initialization has been run
			self.completeness_tracker.update_key_status('initialization', \
				self.combined_results_output_file)
			init_complete = \
				self.completeness_tracker.get_key_completeness('initialization')
			if init_complete:
				self._read_combined_df()
			else:
				# create initial LL profiles and check whether one of
					# them contains a LL higher than the 'unfixed' mode
				self._create_initial_LL_profiles()
				# create self,combined_results_df, using parameter
					# names as index names; calculate runtime
					# parameters and write line to
					# self.combined_results_df, including general
					# warnings for the current mode
				self._create_combined_df()
				# update LL profiles if an ML point was identified
					# that's not unfixed_ll_parameters
				non_unfixed_ML_identified = \
					self.warnings.get_non_unfixed_ML_status()
				if non_unfixed_ML_identified:
					self._correct_LL_profiles()
				self._write_combined_df()
	def generate_asymptotic_CIs(self):
		# generate and record asymptotic CIs for fitted parameters
		self._read_combined_df()
		parameters_to_loop_over = mle_parameters.get_fitted_parameter_list(False)
		for current_fixed_parameter in parameters_to_loop_over:
			# first set current parameter
			self.mle_parameters.set_parameter(fixed_parameter)
			# create an LLProfile for current parameter
			ll_profile = LLProfile(self.mle_parameters, \
				self.mle_datafile_path, self.LL_profile_folder, \
				non_profile_max_params)
			ll_profile.run_LL_profile_compilation()





########################################################################
def _generate_filename(output_file_path, profile_point_as_str, \
	output_id, parameter_name, output_prename):
	# creates a filename to which output file of MLE will be written,
		# to be read by LLProfile
	output_file_label = _generate_file_label(output_prename, output_id, \
		parameter_name)
	output_file = '_'.join([output_file_label, profile_point_as_str]) + '.csv'
	output_filename = os.path.join(output_file_path,output_file)
	return(output_filename)

def _generate_file_label(output_prename, output_id, parameter_name):
	# creates a file label that can be used by Cluster_Functions to
		# track the completeness of the file
	# the file label doesn't include the path, profile point, or file
		# extension
	output_file_label = '_'.join([output_prename, output_id, parameter_name])
	return(output_file_label)

def _get_MLE_params(current_param_datafile):
	# get MLE param values and run info from output (csv) file
	current_param_df = pandas.read_csv(current_param_datafile)
	return(current_param_df)
#	with open(current_param_datafile, 'rU') as \
#		current_param_datafile_contents:
#		current_param_data = list(csv.reader(current_param_datafile_contents))
#	current_param_np = numpy.array([float(i) for i in current_param_data[0]])
#	return(current_param_np)

def loop_over_modes(mle_parameters, cluster_parameters, cluster_folders, \
	mle_folders, experiment_path, additional_code_run_keys, \
	additional_code_run_values, output_id_string_start):
	# Handles all of MLE across modes, including confidence
		# interval identification
	for current_mode in mle_parameters.mode_list:
		##### RUN MLE #####
		output_id_string = '_'.join([output_id_string_start, current_mode])
		mle_parameters.set_mode(current_mode, output_id_string)
		MLE_summary_file_path = os.path.join(experiment_path, \
			'_'.join([output_id_string,'MLE_file.csv']))
		# run MLE for current set of parameters
		run_MLE(mle_parameters, cluster_parameters, cluster_folders, mle_folders, \
			additional_code_run_keys, additional_code_run_values)
		# generate LL_profiles and MLE_output file, and identify CIs
		current_combined_results = CombinedResultSummary(mle_folders, \
			mle_parameters, cluster_parameters, cluster_folders)
		current_combined_results.initialize_combined_results()

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
	parameters_to_loop_over = mle_parameters.get_fitted_parameter_list(True)
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
		# if all parameters for this mode are complete, update mode completeness
		# this also updates completeness across modes
		current_mode_mle_complete_status = \
			mle_parameters.check_completeness_within_mode()
#		if current_mode_complete_status:
#			datafile_path = mle_estimator.get_output_path()
#			LL_profile_folder = mle_folders.LL_profile_path
#			current_ll_profile = LLProfile(mle_parameters, datafile_path, LL_profile_folder, CI_p_val)
#			current_ll_profile.run_LL_profile_compilation()

def id_CI(mle_folders, mle_parameters, cluster_parameters, cluster_folders, \
	fixed_param, CI_completeness_tracker, additional_param_df, datafile_path, \
	LL_profile_folder, p_val):
	# get LL_profile
	# when LL_profile complete, get lower and upper asymptotic CI
	# if asymptotic CI identification is complete:
	#	- identify position at which sims need to happen
	#		- run sims
	#			- get CIs from sims
	# keep track of asymptotic and sim CIs using completeness tracker across parameters within mode
	# once CIs complete (either asymptotic only or asymptotic and sim, depending on settings),
	#	create summary file
	# keep track of summary files for each mode; once those are complete, stop running current folder
	##
	# intialize list of warnings for current LL profile and CIs
	CI_and_profile_warning_list = []
	# create LL profile
	current_ll_profile = LLProfile(mle_parameters, datafile_path, LL_profile_folder, additional_param_df)
	current_ll_profile.run_LL_profile_compilation()
	# get data frame with LL profile, warnings, runtime, and MLE vals
	current_LL_df = current_ll_profile.get_LL_df()
	current_fixed_param_MLE_val = current_ll_profile.get_fixed_param_MLE_val()
	current_max_LL = current_ll_profile.get_max_LL()
	LL_profile_warnings = current_ll_profile.get_warnings()
	CI_and_profile_warning_list.append(LL_profile_warnings)
#	runtime_quantile = current_ll_profile.get_runtime(runtime_percentile)

	###### BEFORE RUNNING ANY CI COMPUTATION, NEED TO MAKE SURE THAT MAX_LL VAL IDENTIFIED IN ALL LL PROFILES IS THE SAME!
		###### IF NOT, FORCE RE-RUN OF LL PROFILE COMPUTATION WITH HIGHEST MAX_LL VALUE

	# get asymptotic CI
	deg_freedom = 1
		# 1 df for chi-square test for LL profile comparisons
	CI_type = 'asymptotic'
	asymptotic_CI = TwoSidedCI(p_val, current_LL_df, deg_freedom, current_fixed_param_MLE_val, \
		current_fixed_param, current_max_LL, CI_type, mle_folders, cluster_parameters, \
		cluster_folders)
	asymptotic_CI.find_CI()
	asymptotic_CI_dict = asymptotic_CI.get_CI()
	if asymptotic_CI_dict:
		asymptotic_CI_warnings = asymptotic_CI.get_CI_warning()
		CI_and_profile_warning_list.append(asymptotic_CI_warnings)
		CI_completeness_tracker.switch_key_completeness(fixed_param, True)
		###### Really, completeness tracker should only be switched once sim CI completed, or if sims not required


##### WHEN COMBINING DATA INTO OUTPUT TABLE, REMEMBER TO INCLUDE CALCULATION OF TRUE MAX LL
			







