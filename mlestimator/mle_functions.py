#!/usr/bin/python

# Contains objects needed for running Maximum Likelihood Estimation

import os
import numpy
from cluster_wrangler import cluster_functions
import copy
import pandas
import csv
from mlestimator.mle_filenaming_functions import generate_file_label, generate_filename
import mle_CI_functions
import warnings as war
war.filterwarnings("ignore", message="numpy.dtype size changed")

class FolderManager(object):
	def __init__(self, cluster_parameters, cluster_folders, \
		experiment_folder_name):
		self.path_dict = {}
		self.experiment_folder_name = experiment_folder_name
		self.path_dict['experiment_path'] = \
			os.path.join(cluster_parameters.composite_data_path, \
				experiment_folder_name)
	#	self.path_dict['sim_output_path'] = \
	#		os.path.join(cluster_parameters.temp_storage_path, \
	#			experiment_folder_name, 'simulated_phenotypes')
		self.path_dict['MLE_output_path'] = \
			os.path.join(cluster_parameters.temp_storage_path, \
				experiment_folder_name,'MLE_output')
		self.path_dict['completefile_folder'] = \
			cluster_folders.get_path('completefile_path')
		self.path_dict['LL_list_path'] = \
			os.path.join(self.path_dict['experiment_path'],'LL_profiles')
		self.path_dict['CI_bound_path'] = \
			os.path.join(self.path_dict['experiment_path'],'CI_bounds')
		self.setup_complete_file = \
			os.path.join(self.path_dict['completefile_folder'], \
			'mle_folder_setup_complete.txt')
	#	self.epilogue_path = os.path.join(cluster_parameters.home_path,cluster_parameters.username,'mut_effect_epilogue_files',experiment_folder_name)
	#	self.path_dict['MLE_combined_sim_path'] = \
	#		os.path.join(cluster_parameters.temp_storage_path, \
	#			'MLE_sim_outputs')
		self.path_dict['sim_output_path'] = \
			os.path.join(cluster_parameters.temp_storage_path, \
				experiment_folder_name, 'simulated_phenotypes')
		self.path_dict['sim_profile_fixed_pt_folder'] = \
			os.path.join(self.path_dict['LL_list_path'], \
				'sim_profile_fixed_points')
		self.path_dict['sim_profile_fixed_pt_folder'] = \
			os.path.join(self.path_dict['LL_list_path'], \
				'sim_profile_fixed_points')
		self.path_dict['sim_output_list_folder'] = \
			os.path.join(self.path_dict['LL_list_path'], \
				'sim_output_list_folder')
		self.path_dict['key_organizer_folder'] = \
			os.path.join(self.path_dict['experiment_path'],'key_organizers')
		# create folders in path_dict
		self._set_up_folders()
		# set up organizer files for sim
		self.hypothesis_key_organizer_file = \
			os.path.join(self.path_dict['key_organizer_folder'], \
				'hypothesis_key_organizer.csv')
		self.sim_key_organizer_file = \
			os.path.join(self.path_dict['key_organizer_folder'], \
				'sim_key_organizer.csv')
	def _set_up_folders(self):
		if not os.path.isfile(self.setup_complete_file):
			for current_folder_key, current_path in self.path_dict.iteritems():
				if not os.path.isdir(current_path):
					os.makedirs(current_path)
			open(self.setup_complete_file,'a').close()
	def get_experiment_folder_name(self):
		return(self.experiment_folder_name)
	def get_hypothesis_key_organizer_file(self):
		return(self.hypothesis_key_organizer_file)
	def get_sim_key_organizer_file(self):
		return(self.sim_key_organizer_file)
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
	# mle_sim_functions.SimParameters inherits from this list
	############### ??? TO DO ??? ###############
	# Check that each mode or parameter is in the list once?
		# (currently only pays attention to the first time mode or parameter listed)
	############### ??? TO DO ??? ###############
	def __init__(self,parameter_list):
		self.input_datafile_keys = parameter_list["input_datafile_keys"]
		self.input_datafile_values = parameter_list["input_datafile_values"]
		self.runtime_percentile = parameter_list["runtime_percentile"]
		self.CI_pval_by_mode = parameter_list['CI_pval_by_mode']
#		self.sim_repeats_by_mode = parameter_list["simulation_repeats_by_mode"]
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
		self.x_tolerance_by_mode = \
			parameter_list["x_tolerance_by_mode"]
		self.fun_tolerance_by_mode = \
			parameter_list["fun_tolerance_by_mode"]
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
		self.current_x_tolerance = self.x_tolerance_by_mode[mode_idx]
		self.current_fun_tolerance = self.fun_tolerance_by_mode[mode_idx]
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

class MLEstimation(cluster_functions.CodeSubmitter):
	'''
	Submits info to cluster_wrangler.cluster_functions.job_flow_handler
	to run matlab code that performs Maximum Likelihood Estimation
	'''
	def __init__(self, mle_parameters, cluster_parameters, cluster_folders, \
		mle_folders, additional_code_run_keys, additional_code_run_values, \
		input_data_folder):
		self.mle_parameters = copy.deepcopy(mle_parameters)
		experiment_folder = mle_folders.get_path('experiment_path')
		completefile = \
			os.path.join(cluster_folders.get_path('completefile_path'), \
				'_'.join(['MLE',mle_parameters.output_id_parameter, \
					'completefile.txt']))
		job_name = '-'.join([mle_folders.get_experiment_folder_name,'MLE', \
			mle_parameters.output_id_parameter])
		job_numbers = [x + 1 for x in \
			list(range(mle_parameters.current_profile_point_num))]
		module = 'matlab'
		parallel_processors = mle_parameters.current_parallel_processors
		output_extension = 'csv'
		code_name = '_'.join(['MLE',mle_parameters.current_mode])
		additional_beginning_lines_in_job_sub = []
		additional_end_lines_in_job_sub = []
		initial_sub_time = cluster_parameters.current_time
		initial_sub_mem = cluster_parameters.current_mem
		self.within_batch_counter_call = \
			cluster_parameters.get_batch_counter_call()
		output_path = mle_folders.get_path('current_output_subfolder')
		output_file_label = generate_file_label('data', \
			mle_parameters.output_identifier, \
			mle_parameters.current_fixed_parameter)
		self.output_filename = generate_filename(output_path, \
			self.within_batch_counter_call, mle_parameters.output_identifier, \
			mle_parameters.current_fixed_parameter, 'data')
		# set up input_datafile_keys and input_datafile_paths
			# attributes, which will be used by
			# _create_code_run_input_lists
		self.input_datafile_keys = mle_parameters.input_datafile_keys
		self.input_datafile_paths = \
			[os.path.join(input_data_folder, \
				current_input_datafile) for current_input_datafile in \
				mle_parameters.input_datafile_values]
		# run __init__ from parent class, which in turn runs
			# _create_code_run_input_lists
		super(cluster_wrangler.cluster_functions.CodeSubmitter, \
			self).__init__(cluster_parameters, \
			cluster_folders, completefile, job_name, \
			job_numbers, module, parallel_processors, \
			experiment_folder, output_extension, code_name, \
			additional_beginning_lines_in_job_sub, \
			additional_end_lines_in_job_sub, initial_sub_time, \
			initial_sub_mem, output_path, output_file_label)
	def _create_code_run_input_lists():
		'''
		Creates list of keys and their values to be submitted to
		external code being run
		'''
		if (len(self.input_datafile_keys) == len(self.input_datafile_paths)) \
			and (len(self.additional_code_run_keys) == \
			len(self.additional_code_run_values)):
			key_list = ['external_counter','combined_fixed_parameter_array', \
				'combined_min_array','combined_max_array','combined_length_array', \
				'combined_position_array','combined_start_values_array', \
				'parameter_list','output_file', \
				'parallel_processors','ms_positions','combined_profile_ub_array', \
				'combined_profile_lb_array','ms_grid_parameter_array', \
				'combined_logspace_parameters', \
				'output_id_parameter', 'combined_scaling_array', 'tolx_val', \
				'tolfun_val'] + self.input_datafile_keys + \
				self.additional_code_run_keys
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
				self.mle_parameters.output_id_parameter, \
				self.mle_parameters.current_scaling_val_list, \
				self.mle_parameters.current_x_tolerance, \
				self.mle_parameters.current_fun_tolerance] + \
				self.input_datafile_paths + self.additional_code_run_values
		else:
			raise RuntimeError('input_datafile_paths or ' + \
				'additional_code_run_values is not the same length as its ' + \
				'corresponding list of keys in MLEstimation class!')

class BoundAbuttingPointRemover(object):
	"""
	Identifies and removes points from LL_df in which at least one
	parameter value abutts a parameter from a list of boundary values
	"""
	def __init__(self, LL_df_prefilter, x_tolerance, bound_array, \
		scaling_array, logspace_parameters, parameter_names, fixed_param):
		self.LL_df_prefilter = LL_df_prefilter
		self.x_tolerance = x_tolerance
		self.bound_array = bound_array
		self.scaling_array = scaling_array
		self.logspace_parameters = logspace_parameters
		# allow fixed_param to abut a boundary
		self.parameter_names = [x for x in parameter_names if x != fixed_param]
		self.fixed_param = fixed_param
		self.parameter_val_df = LL_df_prefilter[parameter_names]
		self.num_rows = LL_df_prefilter.shape[0]
		self.scaling_matrix = np.tile(scaling_array, (self.num_rows, 1))
		# get abs val of difference between rescaled
			# self.parameter_val_df and rescaled bound_array
		self._get_scaled_array_diff()
		self._remove_abutting_points()
	def _rescale_df(self, df):
		# rescales df by converting necessary columns to logspace and
		# then multipying by scaling_array
		logspace_df = copy.copy(df)
		if any(self.logspace_parameters):
			logspace_df[self.logspace_parameters] = \
				np.log(logspace_df[self.logspace_parameters])
		scaled_df = logspace_df * self.scaling_matrix
		return(scaled_df)
	def _get_scaled_array_diff(self):
		# rescales columns in parameter_val_df and finds abs val of the
		# difference between each row and a rescaled comparison_array
		bound_matrix = np.tile(self.bound_array, (self.num_rows, 1))
		bound_df = pd.DataFrame(bound_matrix, \
			index = self.parameter_val_df.index.values, \
			columns = self.parameter_names)
		bound_df_rescaled = self._rescale_df(bound_df)
		parameter_val_df_rescaled = self._rescale_df(self.parameter_val_df)
		self.LL_df_diff = abs(parameter_val_df_rescaled - bound_df_rescaled)
	def _remove_abutting_points(self):
		# identify indices to remove from df, remove them, and save
		# removed vals of fixed_param as self.removed_param_vals
		comparison_df = self.LL_df_diff < self.x_tolerance
		indices_to_remove_bool = comparison_df.any(axis = 'columns')
		indices_to_remove = list(compress(comparison_df.index.values, \
			indices_to_remove_bool))
		self.removed_param_vals = \
			np.array(self.LL_df_prefilter[self.fixed_param].loc[indices_to_remove])
		self.LL_df = self.LL_df_prefilter.drop(indices_to_remove)
	def get_LL_df(self):
		return(self.LL_df)
	def get_removed_param_vals(self):
		return(self.removed_param_vals)

class LLHolder(object):
	'''
	Populates, holds, and saves a dataframe containing outputs of MLE
	'''
	def __init__(self, profile_points, output_identifier, fixed_param, datafile_path, LL_list_folder):
		self.profile_points = profile_points
		self.output_identifier = output_identifier
		self.datafile_path = datafile_path
		self.fixed_param = fixed_param
		self.warning_line = ''
		self.LL_file = os.path.join(LL_list_folder, \
			('_'.join(['LL_file', self.output_identifier, \
				self.fixed_param]) + '.csv'))
		self.LL_df = pandas.DataFrame()
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
	def _populate_LL_df(self):
		''' Fills in data in LL_df '''
		for current_pp in range(1,(self.profile_points+1)):
			current_datafile = generate_filename(self.datafile_path, \
				str(current_pp), self.output_identifier, \
				self.fixed_param, 'data')
			if os.path.isfile(current_datafile):
				current_data_df = _get_MLE_params(current_datafile)
				# add current_data_np to LL array and update max likelihood value
				self._add_vals(current_data_df)
	def _sort_by_profiled_param(self):
		'''
		Returns order of rows in self.LL_df sorted by the parameter
		being profiled
		'''
		if not self.LL_df.empty:
			self.LL_df_sorted = self.LL_df.sort_values(self.fixed_param)
		else:
			self.LL_df_sorted = self.LL_df
	def _write_LL_df(self):
		'''
		Writes LL_df to file regardless of whether file already exists
		'''
		self.LL_df_sorted.to_csv(path_or_buf=self.LL_file,index=False)
	def _read_LL_df(self):
		''' Reads self.LL_df from a pre-recorded file '''
		try:
			self.LL_df = pandas.read_csv(filepath_or_buffer=self.LL_file)
		except pandas.io.common.EmptyDataError:
			# if file is empty, just create an empty df for self.LL_df
			self.LL_df = pandas.DataFrame()
		self._sort_by_profiled_param()
	def _set_LL_df(self):
		'''
		Gets LL_df, either from a pre-recorded file, or from MLE outputs
		'''
		if os.path.isfile(self.LL_file):
			self._read_LL_df()
		else:
			self._populate_LL_df()
		self._sort_by_profiled_param()
	def remove_bound_abutting_points(self, x_tolerance, \
		parameter_max_vals, parameter_min_vals, scaling_array, \
		logspace_parameters, parameter_names):
		# identifies points that abutt parameter_max_vals or
		# parameter_min_vals at a parameter point, and removes them
		# from LL_df
		bound_abutting_point_remover_min = \
			BoundAbuttingPointRemover(self.LL_df_sorted, x_tolerance, \
				parameter_min_vals, scaling_array, logspace_parameters, \
				parameter_names, self.fixed_param)
		LL_df_minfilter = bound_abutting_point_remover_min.get_LL_df()
		min_removed_param_vals = \
			bound_abutting_point_remover_min.get_removed_param_vals()
		bound_abutting_point_remover_max = \
			BoundAbuttingPointRemover(LL_df_minfilter, x_tolerance, \
				parameter_max_vals, scaling_array, logspace_parameters, \
				parameter_names, fixed_param)
		LL_df_postfilter = bound_abutting_point_remover_max.get_LL_df()
		max_removed_param_vals = \
			bound_abutting_point_remover_max.get_removed_param_vals()
		self.removed_param_vals = {'lower' : min_removed_param_vals, \
			'upper' : max_removed_param_vals}
		self._check_bound_abutting_point_warning(fixed_param)
		return(LL_df_postfilter)
	def _check_bound_abutting_point_warning(self, fixed_param):
		for key, val in self.removed_param_vals.iteritems():
			if not val.size == 0:
				current_list_as_str = ';'.join([str(i) for i in val])
				current_warning_string = \
					'datapoints for the following values of ' + \
					fixed_param + 'were removed for abutting the ' + \
					key + ' parameter bounds: ' + current_list_as_str
				self.warning_line = self.warning_line + current_warning_string
	def run_LL_list_compilation(self):
		''' Compiles and writes LL_df '''
		self._set_LL_df()
		self._write_LL_df()
	def get_LL_df(self):
		''' Returns sorted LL_df '''
		return(self.LL_df_sorted)
	def get_LL_file(self):
		''' Returns the filepath containing the LL list '''
		return(self.LL_file)


class LLProfile(LLHolder):
	# Gets, holds, and updates log likelihood profile
	def __init__(self, mle_parameters, datafile_path, LL_profile_folder, additional_param_df):
		super(LLHolder, \
			self).__init__(mle_parameters.current_profile_point_num, \
			mle_parameters.output_identifier, \
			mle_parameters.current_fixed_parameter, datafile_path, \
			LL_profile_folder)
		self.additional_param_df = copy.copy(additional_param_df)
		self.max_LL = None
		self.ML_params = None
		self.fixed_param_MLE_val = None
		self.LL_file = os.path.join(LL_profile_folder, \
			('_'.join(['LL_file', self.output_identifier, \
				self.fixed_param]) + '.csv'))
		self.LL_df = pandas.DataFrame()
	def _set_ML_params(self, ml_param_df):
		self.ML_params = copy.copy(ml_param_df)
		self.max_LL = ml_param_df.iloc[0]['LL']
		self.fixed_param_MLE_val = ml_param_df.iloc[0][self.fixed_param]
	def _id_max_LL(self):
		# identifies and sets the parameter values corresponding to the
			# max likelihood
		# id row of LL_df corresponding to max LL
		if not self.LL_df.empty:
			ml_param_df = self.LL_df.iloc[[self.LL_df['LL'].idxmax()]]
				# idxmax OK here because indices are reset during appending
			# set this row to be the ML parameters
			self._set_ML_params(ml_param_df)
	def _set_LL_df(self):
		'''
		Gets LL_df, either from a pre-recorded file, or from MLE outputs
		If true_max_param_df is non-empty, add it to LL_df and update
		max LL-related parameters
		'''
		if os.path.isfile(self.LL_file):
			self._read_LL_df()
		else:
			self._populate_LL_df()
		if not self.additional_param_df.empty:
			self._add_vals(self.additional_param_df)
		self._sort_by_profiled_param()
	def run_LL_list_compilation(self):
		'''
		Compiles and writes LL_df;
		Identifies and sets max LL parameters
		'''
		self._set_LL_df()
		self._write_LL_df()
		self._id_max_LL()
	def run_CI(self, pval, mle_folders, cluster_parameters, cluster_folders, mle_parameters, CI_type):
		# when LL_profile complete, get lower and upper asymptotic CI
		# if asymptotic CI identification is complete:
		#	- identify position at which sims need to happen
		#		- run sims
		#			- get CIs from sims
		# get asymptotic CI
		deg_freedom = 1
			# 1 df for chi-square test for LL profile comparisons
		LL_df_sorted_cleaned = \
			self.remove_bound_abutting_points(self.mle_parameters.current_x_tolerance, \
				self.mle_parameters.current_max_parameter_val_list, \
				self.mle_parameters.current_min_parameter_val_list, \
				self.mle_parameters.current_scaling_val_list, \
				self.mle_parameters.current_logspace_profile_list, \
				self.mle_parameters.current_parameter_list,
				self.mle_parameters.current_fixed_parameter)
		self.CI = mle_CI_functions.TwoSidedCI(pval, \
			LL_df_sorted_cleaned, deg_freedom, self.fixed_param_MLE_val, \
			self.fixed_param, CI_type, mle_folders, \
			cluster_parameters, cluster_folders, self.output_identifier, \
			mle_parameters)
		self.CI.find_CI()
		self.CI_dict = self.CI.get_CI()
		if self.CI_dict:
			self.warning_line = self.CI.get_CI_warning()
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


########################################################################

def _get_MLE_params(current_param_datafile):
	# get MLE param values and run info from output (csv) file
	current_param_df = pandas.read_csv(current_param_datafile)
	return(current_param_df)

def run_MLE(mle_parameters, cluster_parameters, cluster_folders, mle_folders, \
	additional_code_run_keys, additional_code_run_values, \
	include_unfixed_parameter, input_data_folder):
	# Handles all of MLE for a particular mode, including confidence
		# interval identification
	# Loops through parameters for particular mode, finding ML
		# parameters at every fixed parameter value
	# Compiles results together to find asymptotic CI values (based on
		# chi-sq distribution of 2*log(LR))
	# Runs through simulations to find simulation-based CI values
	# mle_parameters must have the mode already set
	parameters_to_loop_over = \
		mle_parameters.get_fitted_parameter_list(include_unfixed_parameter)
	for current_fixed_parameter in parameters_to_loop_over:
		# set current parameter
		mle_parameters.set_parameter(current_fixed_parameter)
		# create MLEstimation object
		ml_estimator = MLEstimation(mle_parameters, cluster_parameters, \
			cluster_folders, mle_folders, additional_code_run_keys, \
			additional_code_run_values, input_data_folder)
		# submit and track current set of jobs
		ml_estimator.run_job_submission()
		# track completeness within current mode
		mle_completefile = ml_estimator.get_completefile_path()
		mle_parameters.update_parameter_completeness(mle_completefile)
		# if all parameters for this mode are complete, update mode completeness
		# this also updates completeness across modes
		mle_parameters.check_completeness_within_mode()








