#!/usr/bin/python

# Contains objects needed for running Maximum Likelihood Estimation

import os
import numpy
from cluster_wrangler import cluster_functions
import copy
import pandas
import csv
from mle_filenaming_functions import generate_file_label, generate_filename
import mle_CI_functions

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
		self.scaling_arrays_by_mode = \
			parameter_list["scaling_arrays_by_mode"]
		self.ms_positions = parameter_list["multistart_positions"]
		self.multistart_grid_parameters = parameter_list["multistart_grid_parameters"]
		self.logspace_profile_list_by_mode = parameter_list["logspace_profile_list"]
		self.parallel_processors = parameter_list["parallel_processors"]
		self.mode_completeness_tracker = cluster_functions.CompletenessTracker(self.mode_list)
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
		self.current_scaling_val_list = \
			numpy.array(self._retrieve_current_values(self.scaling_arrays_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		# identify list of parameters that are permanently fixed
		self.current_permafixed_parameter_bool = \
			self.current_max_parameter_val_list == self.current_min_parameter_val_list
		# identify parameters MLE needs to be performed on
		self._id_parameters_to_loop_over()
		# set up completefile tracker for these parameters
		self.parameter_completeness_tracker = \
			cluster_functions.CompletenessTracker(self.current_parameters_to_loop_over)
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
		self.current_ms_grid_dimensions = sum([x != self.current_fixed_parameter \
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
		job_submission_manager = cluster_parameters.get_job_sub_manager()
		self.within_batch_counter = \
			job_submission_manager.get_within_batch_counter()
		self.within_batch_counter_call = \
			'${' + self.within_batch_counter + '}'
		self.output_path = mle_folders.get_path('current_output_subfolder')
		self.output_filename = generate_filename(self.output_path, \
			self.within_batch_counter_call, mle_parameters.output_identifier, \
			mle_parameters.current_fixed_parameter, 'data')
		self.output_file_label = generate_file_label('data', \
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
			'combined_logspace_parameters','datafile_path', \
			'output_id_parameter', 'combined_scaling_array']
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
			self.mle_parameters.current_logspace_profile_list,
			self.output_path, self.mle_parameters.output_id_parameter, \
			self.mle_parameters.current_scaling_val_list]
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
			cluster_functions.SubmissionStringProcessor(self.module, key_list, value_list, \
				self.code_name)
		self.code_run_input = submission_string_processor.get_code_run_input()
	def run_job_submission(self):
		# handles submission of the job
		job_name = self.job_name
		job_numbers = [x + 1 for x in \
			list(range(self.mle_parameters.current_profile_point_num))]
		initial_time = self.cluster_parameters.current_time
		initial_mem = self.cluster_parameters.current_mem
			# initial_mem scaled for multiple parallel clusters at a later stage
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
		cluster_functions.job_flow_handler(job_name, job_numbers, initial_time, \
			initial_mem, cluster_parameters, output_folder, output_extension, \
			output_file_label, cluster_job_submission_folder, experiment_folder, \
			module, code_run_input, additional_beginning_lines_in_job_sub, \
			additional_end_lines_in_job_sub, parallel_processors, completefile_path)

class LLProfile(object):
	# Gets, holds, and updates log likelihood profile
	def __init__(self, mle_parameters, datafile_path, LL_profile_folder, additional_param_df):
		self.profile_points = mle_parameters.current_profile_point_num
		self.output_identifier = mle_parameters.output_identifier
		self.current_fixed_parameter = mle_parameters.current_fixed_parameter
		self.datafile_path = datafile_path
		self.fixed_param = mle_parameters.current_fixed_parameter
		self.warning_line = ''
		self.additional_param_df = copy.copy(additional_param_df)
		self.max_LL = None
		self.ML_params = None
		self.fixed_param_MLE_val = None
		self.LL_file = os.path.join(LL_profile_folder, \
			('_'.join(['LL_file', self.output_identifier, \
				self.current_fixed_parameter]) + '.csv'))
		self.LL_df = pandas.DataFrame()
	def _set_ML_params(self, ml_param_df):
		self.ML_params = copy.copy(ml_param_df)
		self.max_LL = ml_param_df.iloc[0]['LL']
		self.fixed_param_MLE_val = ml_param_df.iloc[0][self.fixed_param]
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
	def _id_max_LL(self):
		# identifies and sets the parameter values corresponding to the
			# max likelihood
		# id row of LL_df corresponding to max LL
		if not self.LL_df.empty:
			ml_param_df = self.LL_df.iloc[[self.LL_df['LL'].idxmax()]]
				# idxmax OK here because indices are reset during appending
			# set this row to be the ML parameters
			self._set_ML_params(ml_param_df)
	def _populate_LL_df(self):
		# fill in data in LL_df
		for current_pp in range(1,(self.profile_points+1)):
			current_datafile = generate_filename(self.datafile_path, \
				str(current_pp), self.output_identifier, \
				self.current_fixed_parameter, 'data')
			if os.path.isfile(current_datafile):
				current_data_df = _get_MLE_params(current_datafile)
				# add current_data_np to LL array and update max likelihood value
				self._add_vals(current_data_df)
	def _sort_by_profiled_param(self):
		# returns order of rows in self.LL_df sorted by the
			# parameter being profiled
		if not self.LL_df.empty:
			self.LL_df_sorted = self.LL_df.sort_values(self.fixed_param)
		else:
			self.LL_df_sorted = self.LL_df
	def _write_LL_df(self):
		# write LL_df to file regardless of whether file already
			# exists, since LL_df may have been updated
		self.LL_df_sorted.to_csv(path_or_buf=self.LL_file,index=False)
	def _read_LL_df(self):
		# reads self.LL_df from a pre-recorded file
		try:
			self.LL_df = pandas.read_csv(filepath_or_buffer=self.LL_file)
		except pandas.io.common.EmptyDataError:
			# if file is empty, just create an empty df for self.LL_df
			self.LL_df = pandas.DataFrame()
		self._sort_by_profiled_param()
	def _set_LL_df(self):
		# gets LL_df, either from a pre-recorded file, or from MLE outputs
		# if true_max_param_df is non-empty, add it to LL_df and update
			# max LL-related parameters
		if os.path.isfile(self.LL_file):
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
	def run_CI(self, pval, mle_folders, cluster_parameters, cluster_folders):
		# when LL_profile complete, get lower and upper asymptotic CI
		# if asymptotic CI identification is complete:
		#	- identify position at which sims need to happen
		#		- run sims
		#			- get CIs from sims
		# get asymptotic CI
		deg_freedom = 1
			# 1 df for chi-square test for LL profile comparisons
		CI_type = 'asymptotic'
		self.asymptotic_CI = mle_CI_functions.TwoSidedCI(pval, self.LL_df_sorted, deg_freedom, \
			self.fixed_param_MLE_val, self.fixed_param, self.max_LL, CI_type, \
			mle_folders, cluster_parameters, cluster_folders, self.output_identifier)
		self.asymptotic_CI.find_CI()
		self.asymptotic_CI_dict = self.asymptotic_CI.get_CI()
		if self.asymptotic_CI_dict:
			self.warning_line = self.asymptotic_CI.get_CI_warning()
	def get_LL_df(self):
		# returns sorted LL_df
		return(self.LL_df_sorted)
	def get_fixed_param_MLE_val(self):
		# returns MLE value of current fixed parameter
		return(self.fixed_param_MLE_val)
	def get_asymptotic_CI(self):
		return(self.asymptotic_CI_dict)
	def get_max_LL(self):
		return(self.max_LL)
	def get_ML_params(self):
		return(self.ML_params)
	def get_warnings(self):
		return(self.warning_line)
	def get_runtime(self, runtime_display_percentile):
		# finds the runtime, in hours, corresponding to
			# self.runtime_display_percentile
		percentile_as_decimal = runtime_display_percentile/100
		try:
			time_quantile_in_seconds = \
				self.LL_df_sorted.runtime_in_secs.quantile(percentile_as_decimal)
			time_quantile_in_hours = time_quantile_in_seconds/3600
		except AttributeError:
			time_quantile_in_hours = numpy.NaN
		return(time_quantile_in_hours)

class SingleParamResultSummary(object):
	# stores MLE estimates and CIs for a single parameter
	def __init__(self):
		self.content_dict = {}
	def _get_keys_by_searchstring(self, searchstring):
		# returns keys of self.content_dict containing searchstring
		searchstring_keys = [key for key in self.content_dict.keys() if \
			searchstring in key]
		return(searchstring_keys)
	def read_from_df_line(self, df, index_name):
		# reads content from line index_name in dataframe df
		if index_name in df.index.values:
			self.content_dict = df.loc[index_name].to_dict()
		else:
			key_list = list(df.columns.values)
			self.content_dict = {current_key:numpy.NaN for current_key in key_list}
	def set_max_LL(self, max_LL):
		self.content_dict['max_LL'] = max_LL
	def set_param_MLE(self, param_MLE):
		self.content_dict['param_MLE'] = param_MLE
	def set_warnings(self, warnings):
		self.content_dict['warnings'] = warnings
	def set_CI(self, CI_type, CI_bound_dict):
		# adds keys for each element of CI_bound_dict, combining key
			# name with CI_type
		# key in CI_bound_dict is the bound name and the value is the
			# bound value
		for current_key, current_value in CI_bound_dict.iteritems():
			new_key = '_'.join([CI_type, 'CI_bound', current_key])
			self.content_dict[new_key] = current_value
	def check_column_filledness_by_keyword(self, searchstring):
		# check whether any keys containing keyword exist, and whether
			# all of them contain non-NaN values
		searchstring_keys = self._get_keys_by_searchstring(searchstring)
		if not searchstring_keys:
			columns_filled = False
		else:
			columns_filled = \
				not numpy.any([numpy.isnan(self.content_dict[key]) for \
					key in searchstring_keys])
		return(columns_filled)
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
		self.mle_folders = copy.deepcopy(mle_folders)
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
		self.cluster_folders = copy.deepcopy(cluster_folders)
		self.LL_profile_folder = mle_folders.get_path('LL_profile_path')
		self.runtime_percentile = mle_parameters.runtime_percentile
		self.pval = mle_parameters.current_CI_pval
		self._create_combined_output_file()
		self.max_LL = None
		self.warnings = CombinedResultWarning()
		self.completeness_tracker = cluster_functions.CompletenessTracker(['initialization', \
			'asymptotic_CIs', 'simulation-based_CIs'])
		self.runtime_quant_list = numpy.array([])
		self.true_max_param_df = pandas.DataFrame()
		self.combined_results_df = pandas.DataFrame()
		# set completefiles for asymptotic and sim-based CIs
		self.asymptotic_CI_completefile = \
			os.path.join(cluster_folders.get_path('completefile_path'), \
				'_'.join(['asymoptotic_CI',mle_parameters.output_identifier, \
					'completefile.txt']))
		self.sim_CI_completefile = \
			os.path.join(cluster_folders.get_path('completefile_path'), \
				'_'.join(['simulation-based_CI', \
					mle_parameters.output_identifier, 'completefile.txt']))
		# set MLE results from 'unfixed' parameter (i.e. fitting all
			# unfixed params together)
		self.unfixed_mle_file = generate_filename(self.mle_datafile_path, \
			'1', mle_parameters.output_identifier, 'unfixed', 'data')
	def _create_combined_output_file(self):
		# creates the name of the combined output file for the results
		experiment_path = self.mle_folders.get_path('experiment_path')
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
			if not current_ML_params.empty:
				self._check_and_update_ML(current_ML_params, current_fixed_parameter)
			# get runtime
			runtime_quantile = ll_profile_outputs['runtime_quantile']
			self.runtime_quant_list = numpy.append(self.runtime_quant_list,runtime_quantile)
	def _correct_LL_profiles(self):
		# updates LL_profiles
		parameters_to_loop_over = self.mle_parameters.get_fitted_parameter_list(False)
		ml_params_to_include = pandas.concat([self.unfixed_ll_param_df, \
			self.true_max_param_df])
		for current_fixed_parameter in parameters_to_loop_over:
			self._update_LL_profile(ml_params_to_include, current_fixed_parameter)
	def _get_runtime_CI(self):
		# creates 'confidence intervals' across parameters for the
			# self.runtime_percentile-th estimates within a parameter's
			# estimated profile points
		runtime_CI_bounds = {}
		if not self.runtime_quant_list.size == 0:
			runtime_CI_bounds['lower'] = \
				numpy.percentile(self.runtime_quant_list,self.pval/2*100)
			runtime_CI_bounds['upper'] = \
				numpy.percentile(self.runtime_quant_list,(1-self.pval/2)*100)
		else:
			runtime_CI_bounds['lower'] = numpy.NaN
			runtime_CI_bounds['upper'] = numpy.NaN
		return(runtime_CI_bounds)
	def _add_line_to_combined_df(self,fixed_param_name, current_results_dict):
		# creates self.combined_summary or adds a line to it
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
					sort=False)
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
		runtime_CI = self._get_runtime_CI()
		# get warning line for current mode
		mode_warning_line = self.warnings.get_warning_line()
		# combine data into dictionary that can be added to df
		current_results_line = SingleParamResultSummary()
		current_results_line.set_max_LL(self.max_LL)
		current_results_line.set_param_MLE(mean_runtime_quant)
		current_results_line.set_warnings,(mode_warning_line)
		current_results_line.set_CI('asymptotic', runtime_CI)
		current_results_dict = current_results_line.get_contents()
		# make a line to write to combined_results_df
		self._add_line_to_combined_df(time_string, current_results_dict)
		# append a DataFrame where row names are parameter names
		if not self.true_max_param_df.empty:
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
			self._set_unfixed_param_data()
			# check whether initialization has been run
			self.completeness_tracker.update_key_status('initialization', \
				self.combined_results_output_file)
			init_complete = \
				self.completeness_tracker.get_key_completeness('initialization')
			if init_complete:
				self._read_combined_df()
			else:
#				self._check_and_update_ML(self.unfixed_ll_param_df,'unfixed')
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
		init_complete = \
			self.completeness_tracker.get_key_completeness('initialization')
		if init_complete:
			self._read_combined_df()
			# check whether asymptotic CI has been completed
			self.completeness_tracker.update_key_status('asymptotic_CI', \
				self.asymptotic_CI_completefile)
			asymptotic_CIs_complete = \
				self.completeness_tracker.get_key_completeness('asymptotic_CIs')
			if not asymptotic_CIs_complete:
				parameters_to_loop_over = \
					self.mle_parameters.get_fitted_parameter_list(False)
				asymptotic_CI_completeness_tracker = \
					cluster_functions.CompletenessTracker(parameters_to_loop_over)
				for current_fixed_parameter in parameters_to_loop_over:
					# check whether this asymptotic CI has been identified
					# read line for current_fixed_parameter from
						# combined_results_df
					current_results_line = SingleParamResultSummary()
					current_results_line.read_from_df_line(self.combined_results_df, \
						current_fixed_parameter)
					# identify columns containing the word 'asymptotic_CI'
						# and check that such columns exist and that none
						# of them contain NaN)
					current_asymptotic_CI_complete = \
						current_results_line.check_column_filledness_by_keyword('asymptotic_CI')
					if current_asymptotic_CI_complete:
						asymptotic_CI_completeness_tracker.switch_key_completeness(current_fixed_parameter,True)
					else:
						# first set current parameter
						self.mle_parameters.set_parameter(current_fixed_parameter)
						# create an LLProfile for current parameter
						ll_profile = LLProfile(self.mle_parameters, \
							self.mle_datafile_path, self.LL_profile_folder, \
							pandas.DataFrame())
						ll_profile.run_LL_profile_compilation()
						ll_profile.run_CI(self.pval, self.mle_folders, \
							self.cluster_parameters, self.cluster_folders)
						asymptotic_CI_dict = ll_profile.get_asymptotic_CI()
							# if jobs to calculate both CI bounds have not yet been
								# completed, returns None
						if asymptotic_CI_dict:
							# set the new confidence interval
							current_results_line.set_CI('asymptotic', asymptotic_CI_dict)
							# replace previous warning entry--all warnings are
								# recalculated individually when a CI bound is
								# initialized, so having a CI_bound completely
								# identified previously will not cause a problem
								# in skipping warnings
							current_warning = ll_profile.get_warnings()
							current_results_line.set_warnings(current_warning)
							# replace line in combined_results_df with updated line
								# that has asymptotic CI and warnings
							current_results_dict = current_results_line.get_contents()
							self._add_line_to_combined_df(current_fixed_parameter, current_results_dict)
							asymptotic_CI_completeness_tracker.switch_key_completeness(current_fixed_parameter,True)
				self._write_combined_df()
				asymptotic_CIs_just_completed = \
					asymptotic_CI_completeness_tracker.get_completeness()
				if asymptotic_CIs_just_completed:
					open(self.asymptotic_CI_completefile,'a').close()


			# keep track of asymptotic and sim CIs using completeness tracker across parameters within mode
		# once CIs complete (either asymptotic only or asymptotic and sim, depending on settings),
		#	create summary file
		# keep track of summary files for each mode; once those are complete, stop running current folder
		##





########################################################################

def _get_MLE_params(current_param_datafile):
	# get MLE param values and run info from output (csv) file
	current_param_df = pandas.read_csv(current_param_datafile)
	return(current_param_df)

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
		current_combined_results.generate_asymptotic_CIs()

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








