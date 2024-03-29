#!/usr/bin/python

# Contains objects needed for running Maximum Likelihood Estimation

import os
import numpy as np
from cluster_wrangler import cluster_functions
import copy
import pandas as pd
import csv
from mlestimator.mle_filenaming_functions import generate_file_label, generate_filename
import mle_CI_functions
import warnings as war
war.filterwarnings("ignore", message="numpy.dtype size changed")

class FolderManager(object):
	def __init__(self, cluster_parameters, cluster_folders, \
		experiment_folder_name):
		self.path_dict = {}
		temp_storage_folder = \
			cluster_parameters.get_input_option('temp_storage_folder')
		self.experiment_folder_name = experiment_folder_name
		self.path_dict['experiment_path'] = \
			os.path.join(\
				cluster_parameters.get_input_option('composite_data_path'), \
				experiment_folder_name)
		self.path_dict['MLE_output_path'] = \
			os.path.join(\
				temp_storage_folder, \
				experiment_folder_name,'MLE_output')
		self.path_dict['completefile_folder'] = \
			cluster_folders.get_path('completefile_path')
		self.path_dict['LL_list_path'] = \
			os.path.join(self.path_dict['experiment_path'],'LL_profiles')
		self.path_dict['CI_bound_path'] = \
			os.path.join(\
				temp_storage_folder, \
				experiment_folder_name, 'CI_bounds')
		self.setup_complete_file = \
			os.path.join(self.path_dict['completefile_folder'], \
			'mle_folder_setup_complete.txt')
	#	self.epilogue_path = os.path.join(cluster_parameters.home_path,cluster_parameters.username,'mut_effect_epilogue_files',experiment_folder_name)
	#	self.path_dict['MLE_combined_sim_path'] = \
	#		os.path.join(cluster_parameters.temp_storage_folder, \
	#			'MLE_sim_outputs')
		self.path_dict['sim_output_path'] = \
			os.path.join(\
				temp_storage_folder, \
				experiment_folder_name, 'simulated_phenotypes')
		self.path_dict['sim_profile_fixed_pt_folder'] = \
			os.path.join(self.path_dict['LL_list_path'], \
				'sim_profile_fixed_points')
#		self.path_dict['sim_output_list_folder'] = \
#			os.path.join(self.path_dict['LL_list_path'], \
#				'sim_output_list_folder')
		self.path_dict['sim_output_list_folder'] = \
			os.path.join(temp_storage_folder, \
				experiment_folder_name, 'sim_output_list_folder')
		self.path_dict['key_organizer_home_folder'] = \
			os.path.join(self.path_dict['experiment_path'],'key_organizers')
		self.path_dict['key_organizer_folder'] = \
			os.path.join(temp_storage_folder, \
				experiment_folder_name, 'key_organizers')
		self.path_dict['mle_finder_folder'] = \
			os.path.join(\
				cluster_parameters.get_input_option('pipeline_path'), \
				'mle_finder')
		self.path_dict['CI_finder_folder'] = \
			os.path.join(\
				cluster_parameters.get_input_option('pipeline_path'), \
				'CI_finder')
		# create folders in path_dict
		self._set_up_folders()
		# set up organizer files for sim
		self.hypothesis_key_organizer_file = \
			os.path.join(self.path_dict['key_organizer_folder'], \
				'hypothesis_key_organizer.csv')
		self.sim_key_organizer_file = \
			os.path.join(self.path_dict['key_organizer_folder'], \
				'sim_key_organizer.csv')
		self.path_dict['current_output_subfolder'] = \
			self.path_dict['MLE_output_path']
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
		self.path_dict['current_output_subfolder'] = \
			os.path.join(self.path_dict['MLE_output_path'], \
				current_subfolder)
		if not os.path.isdir(self.path_dict['current_output_subfolder']):
			os.makedirs(self.path_dict['current_output_subfolder'])

class MLEParameters(cluster_functions.InputParameterHolder):
	def __init__(self, parameter_input):
		self.model_list = parameter_input["model"]
		self.model_completeness_tracker = \
			cluster_functions.CompletenessTracker(self.model_list)
		self.all_models_complete = False
		required_mle_key_list = ['input_datafile_keys', \
			'input_datafile_values', 'model', 'parameter_list', \
			'top_level_parameters', 'permafixed_parameters', \
			'starting_parameter_vals', 'min_parameter_vals', \
			'max_parameter_vals', 'multistart_positions', \
			'multistart_grid_parameters', 'logspace_profile_parameters', \
			'scaling_array', 'x_tolerance', 'fun_tolerance', \
			'LL_calculator', 'pre_MLE_function', 'post_MLE_function', \
			'gradient_specification', 'bound_abut_safe_parameters']
		required_CI_key_list = ['profile_point_num_list', \
			'profile_lower_limits', 'profile_upper_limits', 'CI_pval', \
			'runtime_percentile']
		complete_required_key_list = \
			required_mle_key_list + required_CI_key_list
		super(MLEParameters, self).__init__(parameter_input, \
			complete_required_key_list)
	def _id_parameters_to_loop_over(self):
		# identify which parameters need to be looped through in MLE
			# i.e. fitted parameters that MLE needs to be performed on
		# include 'unfixed' parameter, in which case no parameter is fixed
		parameters_to_loop_over_bool = \
			np.invert(self.current_option_dict['permafixed_parameter_bool'])* \
			np.invert((self.current_option_dict['profile_point_num_list'] < 1))
		non_permafixed_parameters = [item for (item, bool_val) in \
			zip(self.current_option_dict['parameter_list'], \
				parameters_to_loop_over_bool) if bool_val]
		self.current_option_dict['parameters_to_loop_over'] = ['unfixed'] + non_permafixed_parameters
		# include 1 profile point for 'unfixed' setting
			# (i.e. only perform 'unfixed' MLE once, since no need to
			# get likelihood profile)
		self.point_numbers_to_loop_over = np.append([1], \
			self.current_option_dict['profile_point_num_list']\
				[parameters_to_loop_over_bool])
	def set_model(self, model_name, output_identifier):
		# for all MLE_parameter attributes, retrieves the parameter or
			# list of parameters corresponding to the current model
		self.current_model_complete = False
		# output_identifier is a string that will be included in filenames
		current_option_dict = {'output_identifier': output_identifier, \
			'current_model': model_name}		
		model_idx = self.model_list.index(model_name)
		number_of_models = len(self.input_val_dict['model'])
		self.current_option_dict = self._select_sublist(self.input_val_dict, \
			current_option_dict, model_idx, number_of_models)
		# check that each parameter for current model is in parameter_list once
		if len(self.current_option_dict['parameter_list']) > \
			len(set(self.current_option_dict['parameter_list'])):
			raise AttributeError('Parameter list for model ' + model_name + \
				'contains non-unique parameters: ' + \
				str(self.current_option_dict['parameter_list']))
		# identify list of parameters that are permanently fixed
		permafixed_parameter_bool = \
			np.array([x in self.current_option_dict['permafixed_parameters'] \
				for x in self.current_option_dict['parameter_list']], \
				dtype=bool)
		self.current_option_dict['permafixed_parameter_bool'] = \
			permafixed_parameter_bool
		# for permafixed parameters, set lower and upper profile pt limits to starting_val
			# profile_point_num needs to be reset too, but that can happen when specifying
		self.current_option_dict['profile_lower_limits']\
			[permafixed_parameter_bool] = \
			self.current_option_dict['starting_parameter_vals']\
				[permafixed_parameter_bool]
		self.current_option_dict['profile_upper_limits']\
			[permafixed_parameter_bool] = \
			self.current_option_dict['starting_parameter_vals']\
				[permafixed_parameter_bool]
		# identify parameters MLE needs to be performed on
		self._id_parameters_to_loop_over()
		# set up completefile tracker for these parameters
		self.parameter_completeness_tracker = \
			cluster_functions.CompletenessTracker(\
				self.current_option_dict['parameters_to_loop_over'])
	def get_current_option_dict(self):
		return(self.current_option_dict)
	def get_option(self, key):
		return(self.current_option_dict[key])
	def get_fitted_parameter_list(self, include_unfixed):
		# get list of parameters to loop over for current model
		parameters_to_return = \
			copy.copy(self.current_option_dict['parameters_to_loop_over'])
		if not include_unfixed:
			parameters_to_return.remove('unfixed')
		return(parameters_to_return)
	def get_complete_parameter_list(self):
		# get list of parameters in current model
		return(self.current_option_dict['parameter_list'])
	def set_parameter(self,parameter_name):
		# set current parameter, number of likelihood profile points
			# for it, and create a temporary list of fixed parameters
			# that includes it
		self.current_option_dict['fixed_parameter'] = parameter_name
		permafixed_parameter_bool = \
			copy.copy(self.current_option_dict['permafixed_parameter_bool'])
		if parameter_name == 'unfixed':
			current_profile_point_num = 1
			# list of fixed parameters is unchanged from default
			self.current_option_dict['tempfixed_parameter_bool'] = \
				permafixed_parameter_bool
		else:
			# find index of current_fixed_parameter in parameter list
			current_fixed_parameter_idx = \
				self.current_option_dict['parameter_list'].index(parameter_name)
			# temporarily fix current parameter
			self.current_option_dict['tempfixed_parameter_bool'] = \
				permafixed_parameter_bool
			self.current_option_dict['tempfixed_parameter_bool']\
				[current_fixed_parameter_idx] = True
			current_profile_point_num = \
				self.current_option_dict['profile_point_num_list']\
					[current_fixed_parameter_idx]
		# set 'profile_point_num' option to number determined above
		self.current_option_dict['profile_point_num'] = current_profile_point_num
		# for permafixed parameters, also reset their profile point to
			# the number determined above
		self.current_option_dict['profile_point_num_list']\
			[permafixed_parameter_bool] = \
				current_profile_point_num
		self.current_option_dict['output_id_parameter'] = \
			self.current_option_dict['output_identifier'] + '_' + parameter_name
		# identify how many dimensions are being used in multistart
		self.current_option_dict['ms_grid_dimensions'] = \
			sum([x != parameter_name for x in \
				self.current_option_dict['multistart_grid_parameters']])
		self.current_option_dict['parallel_processors'] = \
			min(self.input_val_dict['parallel_processors'],
				self.current_option_dict['multistart_positions']**\
					self.current_option_dict['ms_grid_dimensions'])
	def update_parameter_completeness(self, completefile):
		# checks whether jobs for current parameter are all complete
		self.parameter_completeness_tracker.update_key_status( \
			self.current_option_dict['fixed_parameter'], completefile)
	def check_completeness_within_model(self):
		# checks whether all parameters within model are complete
		# change model completeness status accordingly
		self.current_model_complete = \
			self.parameter_completeness_tracker.get_completeness()
		self.model_completeness_tracker.switch_key_completeness( \
			self.current_option_dict['current_model'], self.current_model_complete)
		return(self.current_model_complete)
	def check_completeness_across_models(self):
		self.all_models_complete = self.model_completeness_tracker.get_completeness()
		return(self.all_models_complete)

class MLEstimation(cluster_functions.CodeSubmitter):
	'''
	Submits info to cluster_wrangler.cluster_functions.job_flow_handler
	to run code that performs Maximum Likelihood Estimation
	'''
	def __init__(self, mle_parameters, cluster_parameters, cluster_folders, \
		mle_folders, additional_code_run_keys, additional_code_run_values, \
		input_data_folder):
		self.mle_parameters = copy.deepcopy(mle_parameters)
		completefile = \
			os.path.join(cluster_folders.get_path('completefile_path'), \
				'_'.join(['MLE',mle_parameters.get_option('output_id_parameter'), \
					'completefile.txt']))
		job_name = '-'.join([mle_folders.get_experiment_folder_name(), 'MLE', \
			mle_parameters.get_option('output_id_parameter')])
		job_numbers = [x + 1 for x in \
			list(range(mle_parameters.get_option('profile_point_num')))]
		module = mle_parameters.get_option('module')
		parallel_processors = mle_parameters.get_option('parallel_processors')
		output_extension = 'csv'
		code_name = 'MLE_finder'
		additional_beginning_lines_in_job_sub = []
		additional_end_lines_in_job_sub = []
		initial_sub_time = cluster_parameters.get_input_option('current_time')
		initial_sub_mem = cluster_parameters.get_input_option('current_mem')
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self.within_batch_counter_call = \
			cluster_parameters.get_batch_counter_call()
		output_path = mle_folders.get_path('current_output_subfolder')
		output_file_label = generate_file_label('data', \
			mle_parameters.get_option('output_identifier'), \
			mle_parameters.get_option('fixed_parameter'))
		self.output_filename = generate_filename(output_path, \
			self.within_batch_counter_call, \
			mle_parameters.get_option('output_identifier'), \
			mle_parameters.get_option('fixed_parameter'), 'data')
		self.checkpoint_filename = generate_filename(output_path, \
				self.within_batch_counter_call, \
				mle_parameters.get_option('output_identifier'), \
				mle_parameters.get_option('fixed_parameter'), 'checkpoint')
		# set up input_datafile_keys and input_datafile_paths
			# attributes, which will be used by
			# _create_code_run_input_lists
		self.input_datafile_keys = mle_parameters.get_option('input_datafile_keys')
		self.input_datafile_paths = \
			[os.path.join(input_data_folder, \
				current_input_datafile) for current_input_datafile in \
				mle_parameters.get_option('input_datafile_values')]
		# run __init__ from parent class, which in turn runs
			# _create_code_run_input_lists
		code_path = mle_folders.get_path('mle_finder_folder')
		super(MLEstimation, self).__init__(cluster_parameters, \
			cluster_folders, completefile, job_name, \
			job_numbers, module, parallel_processors, \
			output_extension, code_name, \
			additional_beginning_lines_in_job_sub, \
			additional_end_lines_in_job_sub, initial_sub_time, \
			initial_sub_mem, output_path, output_file_label, code_path)
	def _create_code_run_input_lists(self):
		'''
		Creates list of keys and their values to be submitted to
		external code being run
		'''
		if (len(self.input_datafile_keys) == len(self.input_datafile_paths)) \
			and (len(self.additional_code_run_keys) == \
			len(self.additional_code_run_values)):
			mle_param_dict = self.mle_parameters.get_current_option_dict() 
			self.key_list = ['external_counter', 'combined_position_array', \
				'output_file', 'pause_at_end', 'output_path', \
				'output_file_label', 'checkpoint_file', 'write_checkpoint'] + \
				list(mle_param_dict.keys()) + \
				self.input_datafile_keys + self.additional_code_run_keys
			self.value_list = [self.within_batch_counter_call, \
				[self.within_batch_counter_call], \
					# if combined_position_array has length=1, MLE programs
						# interpret it as an array of the correct length
						# with the same value repeated
				self.output_filename, \
				self.cluster_parameters.get_input_option('pause_at_end'), \
				self.output_path, self.output_file_label, \
				self.checkpoint_filename, True] + \
				list(mle_param_dict.values()) + \
				self.input_datafile_paths + self.additional_code_run_values
		else:
			raise RuntimeError('input_datafile_paths or ' + \
				'additional_code_run_values is not the same length as its ' + \
				'corresponding list of keys in MLEstimation class!')

class BoundAbuttingPointRemover(object):
	"""
	Identifies and removes points from LL_df in which either at least
	one or all the parameter values abutt a parameter from a list of
	boundary values
	"""
	def __init__(self, LL_df_prefilter, x_tolerance, bound_array, \
		scaling_array, logspace_parameters, parameter_names, \
		row_removal_criteria = 'any'):
		self.LL_df_prefilter = LL_df_prefilter
		self.x_tolerance = x_tolerance
		# allow fixed_param to abut a boundary, so remove it from self.parameter_names and other related lists
		self.parameter_names = parameter_names
		self.logspace_parameters = logspace_parameters
		self.bound_array = bound_array
		self.scaling_array = scaling_array
		self.num_rows = LL_df_prefilter.shape[0]
		self.scaling_matrix = np.tile(self.scaling_array, (self.num_rows, 1))
		self.row_removal_criteria = row_removal_criteria
		# get abs val of difference between rescaled
			# self.parameter_val_df and rescaled bound_array
		if len(LL_df_prefilter.index) > 0:
			self._set_up_parameter_val_df(LL_df_prefilter)
			self._get_scaled_array_diff()
			self._remove_abutting_points()
		else:
			self.LL_df = LL_df_prefilter
	def _set_up_parameter_val_df(self, LL_df_prefilter):
		'''
		Sets up a dataframe from LL_df_prefilter containing columns
		corresponding to all the parameter names in self.parameter_names
		'''
		LL_df_cols = LL_df_prefilter.columns
		new_columns = \
			list(set(self.parameter_names).difference(set(LL_df_cols)))
		old_columns_to_keep = \
			list(set(self.parameter_names).intersection(set(LL_df_cols)))
		self.parameter_val_df = copy.copy(LL_df_prefilter[old_columns_to_keep])
		for current_col in new_columns:
			self.parameter_val_df[current_col] = np.nan
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
		parameter_val_df_rescaled = \
			self._rescale_df(self.parameter_val_df[self.parameter_names])
		self.LL_df_diff = abs(parameter_val_df_rescaled - bound_df_rescaled)
		# since self.LL_df_diff will throw NaN when Inf or -Inf values
			# are being compared, also use an equality df for comparing
			# parameter_val_df and bound_df (which will have True when
			# values are both Inf or -Inf)
		self.LL_df_equality = parameter_val_df_rescaled == bound_df_rescaled
	def _remove_abutting_points(self):
		# identify indices to remove from df, remove them, and save
		# removed vals of fixed_param as self.removed_param_vals
		comparison_df = (self.LL_df_diff < self.x_tolerance) | self.LL_df_equality
		if self.row_removal_criteria == 'any':
			indices_to_remove_bool = comparison_df.any(axis = 'columns')
		elif self.row_removal_criteria == 'all':
			indices_to_remove_bool = comparison_df.all(axis = 'columns')
		else:
			raise ValueError('Row removal criteria is \'' + \
				self.row_removal_criteria + \
				'\' but must be \'any\' or \'all\'')
		self.indices_to_remove = list(np.compress(indices_to_remove_bool, \
			comparison_df.index.values))
		self.LL_df = self.LL_df_prefilter.drop(self.indices_to_remove)
	def get_LL_df(self):
		return(self.LL_df)
	def get_removed_indices(self):
		return(self.indices_to_remove)

class BoundAbuttingPointRemoverParamAware(BoundAbuttingPointRemover):
	"""
	Identifies and removes points from LL_df in which at least one
	parameter value abutts a parameter from a list of boundary values,
	excluding the 'safe' parameters
	"""
	def __init__(self, LL_df_prefilter, x_tolerance, bound_array, \
		scaling_array, logspace_parameters, parameter_names_original, \
		bound_abut_safe_parameters, fixed_param, row_removal_criteria = 'any'):
		self.fixed_param = fixed_param
		self.bound_abut_safe_parameters = bound_abut_safe_parameters
		# allow bound_abut_safe_parameters to abut a boundary, so remove them from
		# self.parameter_names and other related lists
		parameter_names_filtered = \
			self._remove_safe_parameters(parameter_names_original, \
				parameter_names_original)
		logspace_parameters_filtered = \
			self._remove_safe_parameters(logspace_parameters, logspace_parameters)
		bound_array_filtered = \
			self._remove_safe_parameters(parameter_names_original, bound_array)
		scaling_array_filtered = \
			self._remove_safe_parameters(parameter_names_original, scaling_array)
		super(BoundAbuttingPointRemoverParamAware, \
			self).__init__(LL_df_prefilter, x_tolerance, bound_array_filtered, \
			scaling_array_filtered, logspace_parameters_filtered, \
			parameter_names_filtered, row_removal_criteria)
	def _remove_safe_parameters(self, array_to_find_safe_parameters_in, array_to_filter):
		'''
		Removes indices corresponding to self.bound_abut_safe_parameters in
		array_to_find_safe_parameters_in from array_to_filter
		'''
		filtered_array = []
		for counter, value in enumerate(array_to_filter):
			current_param = array_to_find_safe_parameters_in[counter]
			if current_param not in self.bound_abut_safe_parameters:
				filtered_array.append(value)
		return(filtered_array)
	def get_removed_param_vals(self):
		if len(self.LL_df_prefilter.index) > 0:
			if self.fixed_param == 'unfixed':
				parameters_to_return = self.parameter_names
			else:
				parameters_to_return = self.fixed_param
			removed_param_vals = \
				np.array(self.LL_df_prefilter[parameters_to_return].loc[self.indices_to_remove])
		else:
			removed_param_vals = np.array([])
		return(removed_param_vals)

class LLHolder(object):
	'''
	Populates, holds, and saves a dataframe containing outputs of MLE
	'''
	def __init__(self, mle_parameters, datafile_path, LL_list_folder):
		self.mle_parameters = mle_parameters
		self.profile_points = mle_parameters.get_option('profile_point_num')
		self.output_identifier = mle_parameters.get_option('output_identifier')
		self.fixed_param = mle_parameters.get_option('fixed_parameter')
		self.bound_abut_safe_parameters = \
			copy.copy(mle_parameters.get_option('bound_abut_safe_parameters'))
		if not self.fixed_param == 'unfixed':
			self.bound_abut_safe_parameters.append(self.fixed_param)
		self.datafile_path = datafile_path
		self.warning_line = ''
		self.LL_file = os.path.join(LL_list_folder, \
			('_'.join(['LL_file', self.output_identifier, \
				self.fixed_param]) + '.csv'))
		self.LL_file_unfiltered = os.path.join(LL_list_folder, \
			('_'.join(['pre_cleanup_LL_file', self.output_identifier, \
				self.fixed_param]) + '.csv'))
		self.LL_df = pd.DataFrame()
		self.x_tolerance = mle_parameters.get_option('x_tolerance')
		self.cdf_tolerance = 10^-10
		self.parameter_max_vals = mle_parameters.get_option('max_parameter_vals')
		self.parameter_min_vals = mle_parameters.get_option('min_parameter_vals')
		self.scaling_array = mle_parameters.get_option('scaling_array')
		self.logspace_parameters = mle_parameters.get_option('logspace_profile_parameters')
		self.parameter_names = mle_parameters.get_option('parameter_list')
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
				# check that file not empty (can result from errors in
					# MLE code)
				if not os.stat(current_datafile).st_size == 0:
					try:
						current_data_df = pd.read_csv(current_datafile)
					except:
						war.warn('issue with file '+str(current_datafile))
					# add current_data_np to LL array and update max likelihood value
					self._add_vals(current_data_df)
	def _sort_by_profiled_param(self):
		'''
		Returns order of rows in self.LL_df sorted by the parameter
		being profiled
		'''
		if not self.LL_df.empty and not (self.fixed_param == 'unfixed'):
			self.LL_df_sorted = self.LL_df.sort_values(self.fixed_param)
		else:
			self.LL_df_sorted = self.LL_df
	def _write_LL_df(self):
		'''
		Writes LL_df to file regardless of whether file already exists
		'''
		self.LL_df_cleaned.to_csv(path_or_buf = self.LL_file , index = False)
		# only write unfiltered file if it hasn't been written before,
			# otherwise pre-filtration values get deleted
		if not os.path.isfile(self.LL_file_unfiltered):
			self.LL_df_sorted.to_csv(path_or_buf = self.LL_file_unfiltered, \
				index = False)
	def _read_LL_df(self):
		''' Reads self.LL_df from a pre-recorded file '''
		try:
			self.LL_df = pd.read_csv(filepath_or_buffer=self.LL_file)
		except pd.io.common.EmptyDataError:
			# if file is empty, just create an empty df for self.LL_df
			self.LL_df = pd.DataFrame()
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
		self._remove_bound_abutting_points()
		self._remove_duplicates()
	def _remove_bound_abutting_points(self):
		'''
		Identifies points that abutt parameter_max_vals or
		parameter_min_vals at a parameter point, and removes them from
		LL_df_sorted
		'''
		bound_abutting_point_remover_min = \
			BoundAbuttingPointRemoverParamAware(self.LL_df_sorted, \
				self.x_tolerance, \
				self.parameter_min_vals, self.scaling_array, \
				self.logspace_parameters, self.parameter_names, \
				self.bound_abut_safe_parameters, self.fixed_param)
		LL_df_minfilter = bound_abutting_point_remover_min.get_LL_df()
		min_removed_param_vals = \
			bound_abutting_point_remover_min.get_removed_param_vals()
#		print('##########################')
#		print(self.parameter_max_vals)
#		print(self.parameter_min_vals)
#		print(self.parameter_names)
#		print(min_removed_param_vals)
		bound_abutting_point_remover_max = \
			BoundAbuttingPointRemoverParamAware(LL_df_minfilter, \
				self.x_tolerance, \
				self.parameter_max_vals, self.scaling_array, \
				self.logspace_parameters, self.parameter_names, \
				self.bound_abut_safe_parameters, self.fixed_param)
		LL_df_cleaned = bound_abutting_point_remover_max.get_LL_df()
		max_removed_param_vals = \
			bound_abutting_point_remover_max.get_removed_param_vals()
		self.removed_param_vals = {'lower' : min_removed_param_vals, \
			'upper' : max_removed_param_vals}
		self._check_bound_abutting_point_warning(self.fixed_param)
		self.LL_df_cleaned = LL_df_cleaned
	def _remove_duplicates(self):
		'''
		Identifies points with LL within provided fun_tolerance (or cdf
		within cdf_tolerance set in this class) and all parameter values
		within parameter tolerance of each other and removes them
		'''
		if not self.LL_df_cleaned.empty:
			# identify rows with an LL or cdf value within fun_tolerance of
		# another LL or cdf value
			if 'LL' in self.LL_df_cleaned.columns:
				current_tolerance = self.mle_parameters.get_option('fun_tolerance')
				diff_list = self.LL_df_cleaned.LL.diff()
			elif 'cdf_vals' in self.LL_df_cleaned.columns:
				current_tolerance = self.cdf_tolerance
				diff_list = self.LL_df_cleaned.cdf_vals.diff()
			diff_indices_within_tolerance = \
				[i for i, x in enumerate(diff_list) if abs(x) < current_tolerance]
			diff_rownames_within_tolerance = \
				self.LL_df_cleaned.iloc[diff_indices_within_tolerance].index.values
			rows_to_remove_unfiltered = []
			# loop through identified rows and add any rows in which every
			# parameter value is within parameter tolerance of an identified
			# row to rows_to_remove_unfiltered
			for current_idx in diff_rownames_within_tolerance:
				current_row_param_vals = \
					self.LL_df_cleaned[self.parameter_names].loc[\
						current_idx].values
				current_row_abutting_point_remover = \
					BoundAbuttingPointRemover(self.LL_df_cleaned, \
						self.x_tolerance, current_row_param_vals, \
						self.scaling_array, self.logspace_parameters, \
						self.parameter_names, row_removal_criteria = 'all')
				current_indices_to_remove = \
					current_row_abutting_point_remover.get_removed_indices()
				rows_to_remove_unfiltered.extend(current_indices_to_remove)
			rows_to_remove = \
				set(rows_to_remove_unfiltered) - \
				set(diff_rownames_within_tolerance)
			self.LL_df_cleaned = self.LL_df_cleaned.drop(rows_to_remove)
	def _check_bound_abutting_point_warning(self, fixed_param):
		for key, val in self.removed_param_vals.iteritems():
			if not val.size == 0:
				current_list_as_str = ';'.join([str(i) for i in val])
				current_warning_string = \
					'datapoints for the following values of ' + \
					fixed_param + 'were removed for abutting the ' + \
					key + ' parameter bounds: ' + current_list_as_str
#				print(current_warning_string)
				self.warning_line = self.warning_line + current_warning_string
	def run_LL_list_compilation(self):
		''' Compiles and writes LL_df '''
		self._set_LL_df()
		self._write_LL_df()
	def get_LL_df(self):
		''' Returns sorted LL_df '''
		return(self.LL_df_cleaned)
	def get_LL_file(self):
		''' Returns the filepath containing the LL list '''
		return(self.LL_file)

class LLProfile(LLHolder):
	# Gets, holds, and updates log likelihood profile
	def __init__(self, mle_parameters, datafile_path, LL_profile_folder, additional_param_df):
		super(LLProfile, \
			self).__init__(mle_parameters, datafile_path, LL_profile_folder)
		self.additional_param_df = copy.copy(additional_param_df)
		self.max_LL = None
		self.ML_params = None
		self.fixed_param_MLE_val = None
		self.LL_file = os.path.join(LL_profile_folder, \
			('_'.join(['LL_file', self.output_identifier, \
				self.fixed_param]) + '.csv'))
		self.LL_df = pd.DataFrame()
	def _set_ML_params(self, ml_param_df):
		self.ML_params = copy.copy(ml_param_df)
		if 'LL' in self.LL_df.columns:
			self.max_LL = ml_param_df.iloc[0]['LL']
		self.fixed_param_MLE_val = ml_param_df.iloc[0][self.fixed_param]
	def _id_max_LL(self):
		# identifies and sets the parameter values corresponding to the
			# max likelihood
		# id row of LL_df corresponding to max LL
		if not self.LL_df.empty:
			if 'LL' in self.LL_df.columns:
				ml_param_df = self.LL_df.loc[[self.LL_df['LL'].idxmax()]]
#				print(ml_param_df)
					# idxmax OK here because indices are reset during appending
				# set this row to be the ML parameters
			elif 'cdf_vals' in self.LL_df.columns:
				ml_param_df = \
					self.LL_df.iloc[[int((abs(self.LL_df['cdf_vals'] - 0)).idxmin())]]
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
		self._remove_bound_abutting_points()
		self._remove_duplicates()
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
		self.CI = mle_CI_functions.TwoSidedCI(pval, \
			self.LL_df_cleaned, deg_freedom, self.fixed_param_MLE_val, \
			self.fixed_param, CI_type, mle_folders, \
			cluster_parameters, cluster_folders, self.output_identifier)
		self.CI.find_CI()
		self.CI_dict = self.CI.get_CI()
		if self.CI_dict:
			self.warning_line = self.CI.get_CI_warning()
	def get_fixed_param_MLE_val(self):
		# returns MLE value of current fixed parameter
		return(self.fixed_param_MLE_val)
	def get_CI(self):
		return(self.CI_dict)
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
				self.LL_df_cleaned.runtime_in_secs.quantile(percentile_as_decimal)
			time_quantile_in_hours = time_quantile_in_seconds/3600
		except AttributeError:
			time_quantile_in_hours = np.NaN
		return(time_quantile_in_hours)


########################################################################

def run_MLE(mle_parameters, cluster_parameters, cluster_folders, mle_folders, \
	additional_code_run_keys, additional_code_run_values, \
	include_unfixed_parameter, input_data_folder):
	# Handles all of MLE for a particular model, including confidence
		# interval identification
	# Loops through parameters for particular model, finding ML
		# parameters at every fixed parameter value
	# Compiles results together to find asymptotic CI values (based on
		# chi-sq distribution of 2*log(LR))
	# Runs through simulations to find simulation-based CI values
	# mle_parameters must have the model already set
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
		# track completeness within current model
		mle_completefile = ml_estimator.get_completefile_path()
		mle_parameters.update_parameter_completeness(mle_completefile)
		# if all parameters for this model are complete, update model completeness
		# this also updates completeness across models
		mle_parameters.check_completeness_within_model()








