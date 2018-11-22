#!/usr/bin/python

# Contains functions necessary for running simulations and combining their results

import pandas as pd
import numpy as np
import os
import warnings as war
import mle_functions
from scipy.stats import norm
from mlestimator import mle_filenaming_functions
import cluster_wrangler.cluster_functions
import copy
from shutil import copyfile
from mlestimator.mle_filenaming_functions import generate_filename
war.filterwarnings("ignore", message="numpy.dtype size changed")

#####################################################################
# Need to:
# 1. Simulate data (using correct sim_mode, correct sim_params; don't forget to random-seed or seed with sim #)
# 2. MLE data (using independent mle_mode, correct starting_params, correct fixed_params)
	# (potentially with multiple independent mle_modes, starting_params, fixed_params for same data - need to be able to differentiate in output)
# 3. Compile a list of likelihoods from a given simulation * mle combo
# 4. Potentially compare likelihoods from different simulation * mle combos and create lrt lists
# 5. Potentially compare either likelihood list or lrt list to benchmark likelihood/LRT value, and calculate proportion of list BELOW benchmark
# 6. 	a. For sim CI calculations, the results of #5 (using lrt comparison to LR given 'true' parameters) becomes likelihood profile point
#		b. For model mode comparisons, the results of #5
#			(coming from sims with SIMPLE model and estimation with simple model vs complex model!?!?) are compared to LR given real data,
#			with the comparison p-val becoming the p-val of extra parameters

# FOR NOW, MAYBE DON'T WRITE IN FRAMEWORK FOR MODEL COMPARISON!?
#####################################################################


#####################################################################
class SimParameters(mle_functions.MLEParameters):
	'''
	This class inherits from mle_functions.MLEParameters, but allows
	modification of the mode and parameter variable assignments
	initially made for running hypothesis testing MLE:
		1.	Selecting only one of the modes orignally specified to be
			kept by the object, to simplify completeness tracking
		2.	Changing the starting_param_vals for that mode, if necessary
		3.	If a parameter in a simulation mle mode is meant to be
			fixed, it is forced to be treated as a fixed parameter by
			setting the min and max values, as well as the lower and
			upper profile limits for that parameter,to its starting
			value, and the mle is then only run for that parameter
			-	The number of profile points is reset to the number of
				simulation_repeats
	It also reads in and processes parameters specific to simulations
	'''
	def __init__(self, parameter_list):
		super(SimParameters, self).__init__(parameter_list)
		self.unmod_input_datafile_names = self.input_datafile_values
		# delete self.input_datafile_values so that they can't be passed
			# to downstream code without setting the sim first
		del self.input_datafile_values
		self.sim_CI_parameters_by_mode = \
			parameter_list["sim_CI_parameters_by_mode"]
		self.sim_repeats_by_mode = \
			parameter_list["simulation_repeats_by_mode"]
		self.sim_time = parameter_list["sim_time"]
		self.sim_mem = parameter_list["sim_mem"]
		# reset profile_points_by_mode to number of sim repeats by mode
		self.profile_points_by_mode = self.sim_repeats_by_mode
		# set default parameter_to_select for looping over parameters
		self.parameter_to_select = 'unfixed'
	def set_sim(self, sim_key, within_batch_counter_call):
		'''
		Sets names of input datafiles to retrieve sim data from
		'''
		self.sim_key = sim_key
		self.input_datafile_values = []
		for current_input_datafile_name in self.unmod_input_datafile_names:
			new_input_datafile_label = _generate_sim_file_label(sim_key, \
				current_input_datafile_name)
			new_input_datafile_val = \
				_generate_sim_filename(within_batch_counter_call, \
					new_input_datafile_label)
			self.input_datafile_values.append(new_input_datafile_val)
	def set_mode(self, mode_name, output_identifier):
		mode_idx = self.mode_list.index(mode_name)
		self.current_sim_rep_num = self.sim_repeats_by_mode[mode_idx]
		self.current_sim_CI_parameters = \
			self.sim_CI_parameters_by_mode[mode_idx]
		super(SimParameters, self).set_mode(mode_name, \
				output_identifier)
	def _id_parameters_to_loop_over(self):
		'''
		Replaces same function in mle_functions.MLEParameters
		Since there are no profile calculations, and any fixed
		parameters in the MLE are fixed by setting their min and max
		values equal to their starting vals, run all MLEs as
		parameter_to_select (defaults to 'unfixed', but should be
		changed to the name of the fixed parameter for the sim mode if
		respecify_for_hypothesis_testing has been run with a specified
		fixed_param) and set point_numbers_to_loop_over to the number of
		sim reps
		'''
		self.current_parameters_to_loop_over = [self.parameter_to_select]
		self.point_numbers_to_loop_over = np.array(self.current_sim_rep_num)
	def _select_mode_to_keep(self, mode_name):
		'''
		Keeps only one of the original modes specified in the object,
		resets mode_completeness_tracker accordingly
		'''
		# only allows for one mode to be kept
		mode_idx = self.mode_list.index(mode_name)
		self.CI_pval_by_mode = [self.CI_pval_by_mode[mode_idx]]
		self.sim_repeats_by_mode = [self.sim_repeats_by_mode[mode_idx]]
		self.ms_positions = [self.ms_positions[mode_idx]]
		self.multistart_grid_parameters = \
			[self.multistart_grid_parameters[mode_idx]]
		self.logspace_profile_list_by_mode = \
			[self.logspace_profile_list_by_mode[mode_idx]]
		self.parameters_by_mode = [self.parameters_by_mode[mode_idx]]
		self.x_tolerance_by_mode = [self.x_tolerance_by_mode[mode_idx]]
		self.fun_tolerance_by_mode = [self.fun_tolerance_by_mode[mode_idx]]
		self.min_parameter_vals_by_mode = \
			[self.min_parameter_vals_by_mode[mode_idx]]
		self.max_parameter_vals_by_mode = \
			[self.max_parameter_vals_by_mode[mode_idx]]
		self.starting_parameter_vals_by_mode = \
			[self.starting_parameter_vals_by_mode[mode_idx]]
		self.profile_points_by_mode = [self.profile_points_by_mode[mode_idx]]
		self.profile_lower_limits_by_mode = \
			[self.profile_lower_limits_by_mode[mode_idx]]
		self.profile_upper_limits_by_mode = \
			[self.profile_upper_limits_by_mode[mode_idx]]
		self.scaling_arrays_by_mode = [self.scaling_arrays_by_mode[mode_idx]]
		self.sim_repeats_by_mode = \
			[self.sim_repeats_by_mode[mode_idx]]
		self.sim_CI_parameters_by_mode = \
			[self.sim_CI_parameters_by_mode[mode_idx]]
		# reset mode_completeness_tracker
		self.mode_completeness_tracker = cluster_wrangler.cluster_functions.CompletenessTracker(self.mode_list)
		self.all_modes_complete = False
	def _change_starting_param_vals(self, mode_name, new_starting_param_vals):
		mode_idx = self.mode_list.index(mode_name)
		self.starting_parameter_vals_by_mode[mode_idx] = \
			new_starting_param_vals
	def _set_permafixed_parameter(self, mode_name, fixed_param):
		'''
		Sets fixed_param to 'fixed' (i.e. not fitted by the MLE) in
		mode_name by setting the min_parameter_val and max_parameter_val
		associated with that parameter equal to its starting_param_val
		Sets 'parameter_to_select' to fixed_param; this will be the only
		parameter looped over during mle with these SimParameters
		'''
		mode_idx = self.mode_list.index(mode_name)
		current_parameter_list = self.parameters_by_mode[mode_idx]
		fixed_param_idx = current_parameter_list.index(fixed_param)
		# run _retrieve_current_values on min, max, and starting
			# parameter vals, as well as profile lower and upper
			# limits, for this mode in case any of these lists
			# were inputted as single elements (i.e. all elements of list
			# are the same)
		self.min_parameter_vals_by_mode[mode_idx] = \
			np.array(self._retrieve_current_values(self.min_parameter_vals_by_mode,\
				mode_idx,mode_name,len(current_parameter_list)))
		self.max_parameter_vals_by_mode[mode_idx] = \
			np.array(self._retrieve_current_values(self.max_parameter_vals_by_mode,\
				mode_idx,mode_name,len(current_parameter_list)))
		self.starting_parameter_vals_by_mode[mode_idx] = \
			np.array(self._retrieve_current_values(self.starting_parameter_vals_by_mode,\
				mode_idx,mode_name,len(current_parameter_list)))
		self.profile_lower_limits_by_mode[mode_idx] = \
			np.array(self._retrieve_current_values(self.profile_lower_limits_by_mode,\
				mode_idx,mode_name,len(current_parameter_list)))
		self.profile_upper_limits_by_mode[mode_idx] = \
			np.array(self._retrieve_current_values(self.profile_upper_limits_by_mode,\
				mode_idx,mode_name,len(current_parameter_list)))
		# change min and max vals for fixed_param
		new_min_max_val = \
			self.starting_parameter_vals_by_mode[mode_idx][fixed_param_idx]
		self.min_parameter_vals_by_mode[mode_idx][fixed_param_idx] = \
			new_min_max_val
		self.max_parameter_vals_by_mode[mode_idx][fixed_param_idx] = \
			new_min_max_val
		self.profile_lower_limits_by_mode[mode_idx][fixed_param_idx] = \
			new_min_max_val
		self.profile_upper_limits_by_mode[mode_idx][fixed_param_idx] = \
			new_min_max_val
		self.parameter_to_select = fixed_param
	def change_profile_pt_number(self, new_profile_pt_num):
		'''
		Allows profile_points_by_mode and current_profile_points_list to
		be altered to new_profile_pt_num, which is an int
		'''
		old_profile_pts_by_mode = self.profile_points_by_mode
		new_profile_pts_by_mode = []
		for current_pts_list in old_profile_pts_by_mode:
			new_pts_list = [new_profile_pt_num] * len(current_pts_list)
			new_profile_pts_by_mode.append(new_pts_list)
		self.profile_points_by_mode = new_profile_pts_by_mode
		if hasattr(self, 'current_mode'):
			mode_idx = self.mode_list.index(self.current_mode)
			self.current_profile_point_list = \
				np.array(self._retrieve_current_values(self.profile_points_by_mode,\
					mode_idx, self.current_mode, self.total_param_num))
	def respecify_for_hypothesis_testing(self, mode_name, fixed_param, starting_vals):
		self._select_mode_to_keep(mode_name)
		self._change_starting_param_vals(mode_name, starting_vals)
		if fixed_param != 'unfixed':
			self._set_permafixed_parameter(mode_name, fixed_param)
		
class SimKeyHolder(object):
	'''
	Stores individual simulation keys as a pandas Series
	'''
	def __init__(self, mode, starting_param_vals):
		self.params = \
			pd.Series([mode, starting_param_vals],
				index = ['mode', 'starting_param_vals'])
	def get_params(self):
		return(self.params)

class HypothesisKeyHolder(object):
	'''
	Stores individual hypothesis keys as a pandas Series
	'''
	def __init__(self, mode, fixed_param, starting_param_vals, sim_key):
		### ??? ###
		# Insert a test that fixed_param is one string element?
		### ??? ###
		self.params = \
			pd.Series([mode, fixed_param, starting_param_vals, sim_key],
				index = ['mode', 'fixed_param', 'starting_param_vals', 'sim_key'])
	def get_params(self):
		return(self.params)

class KeyOrganizer(object):
	'''
	Organizes keys and checks if key with current properties already
	exists; index in self.key_df serves as key
	Keys may not contain lists, with the exception of
	'starting_param_vals' column, which should contain a numpy array
	'''
	def __init__(self, key_file, sim_parameters):
		self.key_file = key_file
		self.param_val_colname = 'starting_param_vals'
		self.sim_parameters = copy.deepcopy(sim_parameters)
		self._read_key_file()
	def _read_key_file(self):
		# if self.key_file exists, reads it in;
		# otherwise, creates one
		if os.path.isfile(self.key_file):
			try:
				self.key_df = \
					pd.read_csv(filepath_or_buffer = self.key_file, \
						index_col = 'idx', \
						converters={self.param_val_colname: _csv_list_converter})
			except pd.io.common.EmptyDataError:
				# if file is empty, just create an empty df for self.LL_df
				self._create_key_df()
		else:
			self._create_key_df()
	def _create_key_df(self):
		''' Creates new key_df '''
		self.key_df = pd.DataFrame()
	def _find_matching_key(self, key_holder):
		'''
		Takes in a KeyHolder object and identifies the row in
		self.key_df that matches it, if any do
		'''
		matching_keys = []
		# identify indices that are within x_tolerance of key_holder
			# for any columns that contain parameter values (first
			# filter by columns containing the correct mode to avoid
			# nonsensical parameter value comparisons)
		if len(self.key_df.index) > 0:
			current_mode = key_holder.get_params()['mode']
			key_df_current_mode = \
				self.key_df.loc[self.key_df['mode'] == current_mode]
		else:
			key_df_current_mode = self.key_df
		if len(key_df_current_mode.index) > 0:
			bound_array = list(key_holder.get_params()[self.param_val_colname])
			param_df = \
				pd.DataFrame(data = \
					np.vstack(key_df_current_mode[self.param_val_colname].values), \
					columns = self.sim_parameters.get_complete_parameter_list(), \
					index = key_df_current_mode.index.values.tolist())
			self.sim_parameters.set_mode(current_mode, '')
			bound_abutting_point_identifier = \
				mle_functions.BoundAbuttingPointRemover(param_df, \
					self.sim_parameters.current_x_tolerance, bound_array, \
					self.sim_parameters.current_scaling_val_list, \
					self.sim_parameters.current_logspace_profile_list,
					self.sim_parameters.get_complete_parameter_list(), 'all')
			bound_abutting_point_indices = bound_abutting_point_identifier.get_removed_indices()
		else:
			bound_abutting_point_indices = []
		if len(bound_abutting_point_indices) > 0:
			columns_to_keep = \
				[current_col for current_col in key_df_current_mode.columns if \
					current_col != self.param_val_colname]
			subset_of_key_df = \
				key_df_current_mode.loc[bound_abutting_point_indices][columns_to_keep]
			for i in subset_of_key_df.index.values.tolist():
				new_dict = key_holder.get_params().to_dict()
				old_dict = subset_of_key_df.loc[i].to_dict()
				match_bool_list = []
				match_bool = all([new_dict[j] == old_dict[j] for j in old_dict.keys()])
				if match_bool:
					matching_keys.append(i)		
		if matching_keys:
			key_val = matching_keys[0]
		else:
			key_val = None
		return(key_val)
	def get_key(self, key_holder):
		'''
		Takes in a KeyHolder object and determines whether it's already
		in self.key_df
		If not, adds it, and writes updated file to self.key_file
		Either way, returns the key for this set of KeyHolder
		(i.e. the index of the line in key_df holding these params)
		'''
		current_key = self._find_matching_key(key_holder)
		if current_key is None:
			print('---- ... before appending ... ----')
			print(self.key_df)
			self.key_df = \
				self.key_df.append(key_holder.get_params(), sort = False, \
					ignore_index = True)
			print('---- ... after appending ... ----')
			print(self.key_df)
			current_key = self._find_matching_key(key_holder)
			self.write_key_df()
		return(current_key)
	def set_key(self, key_holder, key_idx):
		''' Sets index key_idx in self.key_df to key_holder '''
		# check that all keys in current_results_dict are columns
			# in self.combined_results_df; if not, create these
		dict_keys = key_holder.get_params().to_dict()
		df_columns = list(self.key_df.columns.values)
		new_columns = set(dict_keys).symmetric_difference(df_columns)
		self.key_df = \
			self.key_df.append(pd.DataFrame(columns=new_columns), \
				sort=False)
		# insert key_holder into correct row in df
		self.key_df.loc[key_idx] = pd.Series(key_holder.get_params())
		self.write_key_df()
	def write_key_df(self):
		''' Writes self.key_df to self.key_file '''
		self.key_df.to_csv(path_or_buf = self.key_file, \
			index = True, index_label = 'idx')

class HypothesisTestingInfo(object):
	'''
	Holds information for performing MLE on a pair of null and
	alternative hypotheses in a pandas DF
	Hypotheses are 'H0' (null) and 'H1' (alternative)
	Need to first run set_hypothesis_parameters method with both H0 and
	H1, then run set_up_sim_key method
	'''
	def __init__(self):
		self.hypotheses = ['H1', 'H0']
		self.hypothesis_info = \
			pd.DataFrame(columns = \
				['mode', 'fixed_param', 'starting_param_vals', 'sim_key', \
					'hypothesis_key'],
				index = self.hypotheses)
	def set_hypothesis_parameters(self, Hnum, hypothesis_key_holder, \
		hypothesis_key_organizer):
		'''
		Sets the mode, fixed_param (None if no parameter is fixed) and
		starting_param_vals (None if default should be used) that will
		be used for hypothesis Hnum ('H1' or 'H0')
		'''
		if Hnum in self.hypotheses:
			self.hypothesis_info.loc[Hnum] = \
				pd.Series(hypothesis_key_holder.get_params())
			current_h_key = \
				hypothesis_key_organizer.get_key(hypothesis_key_holder)
			self.hypothesis_info.at[Hnum, 'hypothesis_key'] = current_h_key
		else:
			raise ValueError("invalid hypothesis (Hnum): " + Hnum + \
				'; Hnum may only be one of the following: ' + \
				', '.join(self.hypotheses))
	def set_up_sim_key(self):
		'''
		Check that hypothesis testing info hypotheses all have a single
		common sim_key (i.e. are performing MLE with different parameter
		settings on the same set of simulated data)
		'''
		unique_sim_key_vals = \
			self.hypothesis_info['sim_key'].unique()
		if len(unique_sim_key_vals) != 1:
			raise ValueError("All hypotheses being compared must be using " + \
				"the same simulation data. The following " + \
				"HypothesisTestingInfo object lists multiple unique " + \
				"sim_keys:\n" + self.hypothesis_info.to_string())
		else:
			self.sim_key = unique_sim_key_vals[0]
	def _get_hypothesis_info(self, Hnum, desired_var):
		if Hnum in self.hypotheses:
			return(self.hypothesis_info.loc[Hnum][desired_var])
		else:
			raise ValueError("invalid hypothesis (Hnum): " + Hnum + \
				'; Hnum may only be one of the following: ' + \
				', '.join(self.hypotheses))
	def get_hypothesis_key_holder(self, Hnum):
		'''
		Returns HypothesisKeyHolder object corresponding to values in
		self.hypothesis_info.loc[Hnum]
		'''
		current_mode = self.get_mode(Hnum)
		current_fixed_param = self.get_fixed_param(Hnum)
		current_starting_param_vals = self.get_starting_param_vals(Hnum)
		current_sim_key = self.sim_key
		current_hypothesis_key = HypothesisKeyHolder(current_mode, \
			current_fixed_param, current_starting_param_vals, current_sim_key)
		return(current_hypothesis_key)
	def get_sim_key(self):
		return(self.sim_key)
	def get_mode(self, Hnum):
		return(self._get_hypothesis_info(Hnum, 'mode'))
	def get_fixed_param(self, Hnum):
		return(self._get_hypothesis_info(Hnum, 'fixed_param'))
	def get_starting_param_vals(self, Hnum):
		return(self._get_hypothesis_info(Hnum, 'starting_param_vals'))
	def get_hypothesis_key(self, Hnum):
		return(self._get_hypothesis_info(Hnum, 'hypothesis_key'))
	def get_hypotheses(self):
		return(self.hypotheses)

class SimPreparer(object):
	'''
	Parent class for classes that run through hypotheses in a
	HypothesisTestingInfo class object and perform some actions (e.g.
	LLRCalculator, SimRunner)
	'''
	def __init__(self, output_id_prefix, sim_parameters, hypothesis_testing_info, \
		cluster_parameters, cluster_folders, sim_folders, 
		additional_code_run_keys, additional_code_run_values):
		self.sim_key = hypothesis_testing_info.get_sim_key()
		self.output_id_sim = '_'.join([output_id_prefix, str(self.sim_key)])
			# output_id_sim needs to include sim key but not sim
				# number; fixed parameter will be added to filenames
				# by SimParameters, and sim number will be added in MLE
		self.within_batch_counter_call = \
			cluster_parameters.get_batch_counter_call()
		self.sim_parameters = copy.deepcopy(sim_parameters)
		self.sim_folders = sim_folders
		self.cluster_parameters = cluster_parameters
		self.cluster_folders = cluster_folders
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self.sim_datafile_path = sim_folders.get_path('sim_output_path')
		self.hypothesis_testing_info = hypothesis_testing_info
		self.hypothesis_list = ['H1','H0']
		self._set_up_sim_parameters()
	def _set_up_sim_parameters(self):
		'''
		Creates a dictionary with a SimParameters object corresponding
		to each hypothesis in hypothesis_testing_info; respecifies
		SimParameters object for hypothesis testing, and sets its sim,
		mode, and fixed_param
		'''
		# It's important to keep sim_parameters object as attribute of
			# LRCalculator so that the objects can be modified to keep
			# track of modes and parameters
		self.sim_param_dict = dict()
		for current_Hnum in self.hypothesis_list:
			# respecify SimParameter object for each hypothesis
			current_sim_params = copy.deepcopy(self.sim_parameters)
			current_mode = self.hypothesis_testing_info.get_mode(current_Hnum)
			current_fixed_param = \
				self.hypothesis_testing_info.get_fixed_param(current_Hnum)
			current_starting_param_vals = \
				self.hypothesis_testing_info.get_starting_param_vals(current_Hnum)
			current_sim_params.respecify_for_hypothesis_testing(current_mode, \
				current_fixed_param, current_starting_param_vals)
			# set sim, mode, fixed_param for current_sim_params
			current_h_key = \
				self.hypothesis_testing_info.get_hypothesis_key(current_Hnum)
			current_output_id = '_'.join([self.output_id_sim, current_Hnum, \
				str(current_h_key)])
			current_sim_params.set_sim(self.sim_key, \
				self.within_batch_counter_call)
			current_sim_params.set_mode(current_mode, current_output_id)
			current_sim_params.set_parameter(current_fixed_param)
			self.sim_param_dict[current_Hnum] = current_sim_params

class LLRCalculator(SimPreparer):
	'''
	For self.data (which can be real or come from a sim) and fixed_param_val:
		1. Calculates H1: 'unfixed' LL (all unknown parameters freely
			estimated)
		2. Calculates H0: LL at fixed_param_val
		3. Calculates LL(H0)-LL(H1)
		4. Calculates Deviance as -2*LLR
	'''
	def __init__(self, output_id_prefix, \
		sim_parameters, hypothesis_testing_info, cluster_parameters, \
		cluster_folders, sim_folders, additional_code_run_keys, \
		additional_code_run_values):
		super(LLRCalculator, self).__init__(output_id_prefix, sim_parameters, \
			hypothesis_testing_info, cluster_parameters, cluster_folders, \
			sim_folders, additional_code_run_keys, additional_code_run_values)
		self.LL_list_folder = sim_folders.get_path('sim_output_list_folder')
		self.mle_datafile_path = sim_folders.get_path('current_output_subfolder')
		self._generate_LLR_filename()
		self.LL_list_dict = dict()
		# create CompletenessTracker object to track whether each LL
			# list is complete
		self.LL_list_completeness_tracker = \
			cluster_wrangler.cluster_functions.CompletenessTracker(self.hypothesis_list)
		self.llr_completeness = False
	def _generate_LLR_filename(self):
		''' Generates name for combined LLR file '''
		hypothesis_key_string = ''
		for current_Hnum in self.hypothesis_list:
			current_h_key = \
				self.hypothesis_testing_info.get_hypothesis_key(current_Hnum)
			hypothesis_key_string = '_'.join([hypothesis_key_string, \
				current_Hnum, str(current_h_key)])
		LLR_file_name = 'LLR_file' + '_' + str(self.output_id_sim) + \
			hypothesis_key_string + '.csv'
		self.LLR_file = os.path.join(self.LL_list_folder, LLR_file_name)
	def _run_MLE(self, current_sim_parameters):
		'''
		Runs MLE for getting LL values to be used in hypothesis testing
		'''
		# run MLE; completeness of current_sim_parameters will
			# automatically be updated
		include_unfixed_parameter = True
		input_data_folder = self.sim_folders.get_path('sim_output_path')
		mle_functions.run_MLE(current_sim_parameters, self.cluster_parameters, \
			self.cluster_folders, self.sim_folders, self.additional_code_run_keys, \
			self.additional_code_run_values, include_unfixed_parameter, \
			input_data_folder)
	def _compile_LL_list(self, Hnum, current_sim_parameters):
		'''
		Compiles and writes LL_list for each hypothesis
		'''
		ll_list = mle_functions.LLHolder(current_sim_parameters, \
			self.mle_datafile_path, self.LL_list_folder)
		ll_list.run_LL_list_compilation()
			# compiles ll_list and writes it to a file
		# get sorted list of LL values with any values in which a
			# parameter (except the fixed parameter) abutted a min or
			# max bound removed from the list
		current_ll_df = ll_list.get_LL_df()
		current_ll_file = ll_list.get_LL_file()
		self.LL_list_dict[Hnum] = current_ll_df
		self.LL_list_completeness_tracker.update_key_status(Hnum, \
			current_ll_file)
	def _find_LLs(self, Hnum):
		'''
		Respecifies sim_parameters for current hypothesis being tested;
		runs MLE and, once that is complete, compiles LL_list for
		hypothesis
		'''
		current_sim_parameters = self.sim_param_dict[Hnum]
		self._run_MLE(current_sim_parameters)
		mle_complete = current_sim_parameters.check_completeness_within_mode()
		if mle_complete:
			self._compile_LL_list(Hnum, current_sim_parameters)
	def _calculate_LLR(self):
		'''
		Matches MLE results for each hypothesis by simulation point and
		calculates log(likelihood ratio) (LLR) for each point
		'''
		abbreviated_LL_df_dict = {}
		for current_Hnum in self.hypothesis_list:
			current_full_df = self.LL_list_dict[current_Hnum]
			current_fixed_param = \
				self.hypothesis_testing_info.get_fixed_param(current_Hnum)
			if len(current_full_df.index) == 0:
				current_full_df = \
					pd.DataFrame(columns = ['LL', 'point_num', current_fixed_param])
			if current_fixed_param is 'unfixed':
				abbreviated_LL_df_dict[current_Hnum] = \
					current_full_df[['LL', 'point_num']]
			else:
				abbreviated_LL_df_dict[current_Hnum] = \
					current_full_df[['LL', 'point_num', current_fixed_param]]
		self.LLR_df = pd.merge(abbreviated_LL_df_dict['H1'], \
			abbreviated_LL_df_dict['H0'], how='outer', left_on = 'point_num', \
			right_on = 'point_num', suffixes = (['_H1','_H0']))
		self.LLR_df['LLR'] = self.LLR_df['LL_H0'] - self.LLR_df['LL_H1']
		self.LLR_df['deviance'] = -2 * self.LLR_df['LLR']
		self.LLR_df.to_csv(path_or_buf = self.LLR_file)
		self.llr_completeness = True
	def run_LLR(self):
		'''
		Runs through steps to submit MLE jobs on simulations, compile
		LL lists, and compute LLR
		'''
		for current_Hnum in self.hypothesis_list:
			self._find_LLs(current_Hnum)
		ll_list_completeness = \
			self.LL_list_completeness_tracker.get_completeness()
		if ll_list_completeness:
			self._calculate_LLR()
	def get_LLR_filepath(self):
		return(self.LLR_file)
	def get_LLR_completeness(self):
		return(self.llr_completeness)
	def get_deviances(self):
		'''
		If llr calculation has been completed, returns np array of
		values in deviance column of self.LLR_df (unless self.LLR_df is
		empty, in which an np array with a single NaN value is returned)
		'''
		if self.llr_completeness:
			if self.LLR_df.empty:
				return(np.array([np.nan]))
			else:
				return(self.LLR_df['deviance'].values)
		else:
			raise RuntimeError('Requestion LLR values before LLR ' + \
				'calculation has been completed')
	def get_LL(self, Hnum):
		return(self.LL_list_dict[Hnum])

class Simulator(cluster_wrangler.cluster_functions.CodeSubmitter):
	'''
	Submits info to cluster_wrangler.cluster_functions.job_flow_handler
	to run matlab code that performs simulations based on original input
	data
	'''
	def __init__(self, sim_parameters, cluster_parameters, cluster_folders, \
		sim_folders, additional_code_run_keys, additional_code_run_values):
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self.sim_parameters = copy.deepcopy(sim_parameters)
		experiment_folder = sim_folders.get_path('experiment_path')
		completefile = \
			os.path.join(cluster_folders.get_path('completefile_path'), \
				'_'.join(['sim',str(sim_parameters.sim_key), \
					'completefile.txt']))
		job_name = '-'.join([sim_folders.get_experiment_folder_name(), \
			'sim', str(sim_parameters.sim_key)])
		job_numbers = [x + 1 for x in \
			list(range(sim_parameters.current_profile_point_num))]
		module = 'matlab'
		parallel_processors = 1
		output_extension = 'csv'
		code_name = '_'.join(['sim',sim_parameters.current_mode])
		additional_beginning_lines_in_job_sub = []
		additional_end_lines_in_job_sub = []
		initial_sub_time = sim_parameters.sim_time
		initial_sub_mem = sim_parameters.sim_mem
		# need to set an output_file_label for the
			# purposes of tracking job completion
		# use last file in sim_parameters.input_datafile_values
		sim_output_path = sim_folders.get_path('sim_output_path')
		output_file_label = sim_parameters.input_datafile_values[-1]
		self.within_batch_counter_call = \
			cluster_parameters.get_batch_counter_call()
		# set up input_datafile_keys and input_datafile_paths
			# attributes, which will be used by
			# _create_code_run_input_lists
		self.input_datafile_keys = sim_parameters.input_datafile_keys
		self.input_datafile_paths = \
			[os.path.join(sim_output_path, current_input_datafile) for \
				current_input_datafile in sim_parameters.input_datafile_values]
		# Need to add original data and also include that in code input
		self._create_original_input_lists(sim_parameters, sim_output_path)
		# run __init__ from parent class, which in turn runs
			# _create_code_run_input_lists
		super(Simulator, self).__init__(cluster_parameters, \
			cluster_folders, completefile, job_name, \
			job_numbers, module, parallel_processors, \
			experiment_folder, output_extension, code_name, \
			additional_beginning_lines_in_job_sub, \
			additional_end_lines_in_job_sub, initial_sub_time, \
			initial_sub_mem, sim_output_path, output_file_label)
	def _create_original_input_lists(self, sim_parameters, sim_output_path):
		''' Creates lists of original data to include in code input '''
		temp_sim_parameters = copy.deepcopy(sim_parameters)
		temp_sim_parameters.set_sim('original', str(1))
		self.original_input_datafile_keys = ['original_' + current_key for \
			current_key in temp_sim_parameters.input_datafile_keys]
		self.original_input_datafile_paths = \
			[os.path.join(sim_output_path, current_input_datafile) for \
				current_input_datafile in \
					temp_sim_parameters.input_datafile_values]
	def _create_code_run_input_lists(self):
		'''
		Creates list of keys and their values to be submitted to
		external code being run
		'''
		if (len(self.original_input_datafile_keys) == \
			len(self.original_input_datafile_paths)) and \
			(len(self.input_datafile_keys) == len(self.input_datafile_paths)) \
				and (len(self.additional_code_run_keys) == \
				len(self.additional_code_run_values)):
			self.key_list = ['external_counter','combined_start_values_array', \
				'parameter_list'] + self.original_input_datafile_keys + \
				self.input_datafile_keys + self.additional_code_run_keys
			self.value_list = [self.within_batch_counter_call, \
				self.sim_parameters.current_start_parameter_val_list, \
				self.sim_parameters.current_parameter_list] + \
				self.original_input_datafile_paths + self.input_datafile_paths \
				+ self.additional_code_run_values
		else:
			raise RuntimeError('original_input_datafile_paths, ' + \
				'input_datafile_paths, or additional_code_run_values is not' + \
				' the same length as its corresponding list of keys in ' + \
				'Simulator class!')
	def get_original_input_datafile_paths(self):
		return(self.original_input_datafile_paths)

class SimRunner(SimPreparer):
	'''
	Runs simulations with parameters corresponding to H1 hypothesis in
	hypothesis_testing_info
	'''
	def __init__(self, output_id_prefix, sim_parameters, hypothesis_testing_info, \
		cluster_parameters, cluster_folders, sim_folders, additional_code_run_keys, \
		additional_code_run_values):
		# create CompletenessTracker object to track whether each LL
			# list is complete
		super(SimRunner, self).__init__(output_id_prefix, sim_parameters, \
			hypothesis_testing_info, cluster_parameters, cluster_folders, \
			sim_folders, additional_code_run_keys, additional_code_run_values)
		self.sim_completeness = False
		# create Simulator object
		print('Profile points to run')
		print(self.sim_param_dict['H1'].current_profile_point_num)
		self.simulator = Simulator(self.sim_param_dict['H1'], \
			self.cluster_parameters, self.cluster_folders, \
			self.sim_folders, additional_code_run_keys, \
			additional_code_run_values)
	def _copy_original_files(self):
		''' Copies original datafiles to simulation folder '''
		experiment_path = self.sim_folders.get_path('experiment_path')
		input_files = [os.path.join(experiment_path, current_datafile_val) \
			for current_datafile_val in \
				self.sim_parameters.unmod_input_datafile_names]
		output_files = self.simulator.get_original_input_datafile_paths()
		if len(input_files) == len(output_files):
			for (current_input_file, current_output_file) in \
				zip(input_files, output_files):
				if not os.path.isfile(current_output_file):
					copyfile(current_input_file, current_output_file)
			self.sim_completeness = True
		else:
			raise RuntimeError('Cannot create copy of original data, ' + \
				'more output files than input files listed')
	def _submit_and_track_sim_job(self):
		''' Submit and track current set of jobs '''
		self.simulator.run_job_submission()
		# track completeness within current mode
		sim_completefile = self.simulator.get_completefile_path()
		self.sim_param_dict['H1'].update_parameter_completeness(sim_completefile)
		# if all parameters for this mode are complete, update mode completeness
		# this also updates completeness across modes
		self.sim_completeness = \
			self.sim_param_dict['H1'].check_completeness_within_mode()
	def submit_sim(self):
		'''
		If sim_key is 'original', copies original datafiles to
		destination of 'original' datafiles in simulation folder;
		otherwise, submits code to run sim with parameters corresponding
		to H1
		'''
		if self.sim_key == 'original':
			self._copy_original_files()
		else:
			self._submit_and_track_sim_job()
	def get_sim_completeness(self):
		'''
		Returns bool of whether or not all sim jobs have run
		'''
		return(self.sim_completeness)

class FixedPointCDFvalCalculator(object):
	'''
	Performs simulations and MLEs to calculate cumulative density
	function val given parameters in hypothesis_testing_info
	To calculate cdf val at a given point:
		1. 	a.	Calculate H1: the hypothesis with the higher number of
				degrees of freedom; for comparisons of hypothesis in one
				of which a parameter is set to a predetermined value,
				this corresponds to the 'unfixed' LL (all unknown
				parameters freely estimated)
			b.	Calculate H0: the hypothesis with the lower number of
				degrees of freedom; LL at fixed point in parameter space
			c. 	'original' log likelihood ratio 'LLR' is LL(H0)-LL(H1)
				(actually to make the cdf vals easier to understand,
				it's easier to work with the deviance, -2*LLR)
			d. 	MLE parameter values for all parameters estimated in (b)
				will be used as sim starting values!
		2. 	Run sim_number simulations using starting parameter values
			from (1d)
		3. Repeat and record (1) but with sim data rather than real data
	'''
	def __init__(self, mode_dict, fixed_param_dict, \
		fixed_param_val_dict, sim_parameters, sim_folders, \
		additional_code_run_keys, additional_code_run_values, \
		output_id_prefix, output_file, cluster_folders, cluster_parameters):
		self.cluster_folders = cluster_folders
		self.cluster_parameters = cluster_parameters
		self.sim_folders = copy.deepcopy(sim_folders)
		self.hypothesis_key_organizer = \
			KeyOrganizer(sim_folders.get_hypothesis_key_organizer_file(), \
				sim_parameters)
		self.sim_key_organizer = \
			KeyOrganizer(sim_folders.get_sim_key_organizer_file(), \
				sim_parameters)
		self.hypotheses = ['H1','H0']
		self._set_up_sim_params_by_data_type(sim_parameters)
		original_data_completeness_tracker = cluster_wrangler.cluster_functions.CompletenessTracker(['sim','LLR'])
		simulated_data_completeness_tracker = cluster_wrangler.cluster_functions.CompletenessTracker(['sim','LLR'])
		self.completeness_tracker_dict = \
			{'original': original_data_completeness_tracker, \
			'simulated': simulated_data_completeness_tracker}
		self.llr_calculator_dict = dict()
		self.hypothesis_testing_info_dict = dict()
		self.mode_dict = mode_dict
		self.fixed_param_dict = fixed_param_dict
		self.fixed_param_val_dict = fixed_param_val_dict
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self.output_id_prefix = output_id_prefix
		self.output_file = output_file
		self.cdf_val_calc_complete = False
	def _set_up_sim_params_by_data_type(self, sim_parameters):
		self.sim_param_dict = {}
		self.sim_param_dict['simulated'] = copy.deepcopy(sim_parameters)
		self.sim_param_dict['original'] = copy.deepcopy(sim_parameters)
		self.sim_param_dict['original'].change_profile_pt_number(1)
	def _get_starting_params(self, Hnum, current_sim_parameters):
		current_sim_parameters.set_mode(self.mode_dict[Hnum], self.output_id_prefix)
		current_sim_parameters.set_parameter(self.fixed_param_dict[Hnum])
		parameter_names = current_sim_parameters.get_complete_parameter_list()
		starting_vals = \
			copy.copy(current_sim_parameters.current_start_parameter_val_list)
		current_fixed_param = current_sim_parameters.current_fixed_parameter
		if current_fixed_param != 'unfixed':
			fixed_param_idx = parameter_names.index(self.fixed_param_dict[Hnum])
			starting_vals[fixed_param_idx] = self.fixed_param_val_dict[Hnum]
		return(np.array(starting_vals))
	def _set_up_hypothesis_testing_info(self, data_type, H1_key_holder, \
		H0_key_holder):
		'''
		Sets up HypothesisTestingInfo object with H1_key_holder and
		H0_key_holder in self.hypothesis_testing_info_dict[data_type]
		'''
		current_hypothesis_testing_info = HypothesisTestingInfo()
		current_hypothesis_testing_info.set_hypothesis_parameters('H1', \
			H1_key_holder, self.hypothesis_key_organizer)
		current_hypothesis_testing_info.set_hypothesis_parameters('H0', \
			H0_key_holder, self.hypothesis_key_organizer)
		current_hypothesis_testing_info.set_up_sim_key()
		self.hypothesis_testing_info_dict[data_type] = \
			current_hypothesis_testing_info
	def _generate_hypothesis_key_holder(self, Hnum, sim_key, \
		starting_param_vals):
		''' Generates HypothesisKeyHolder for Hnum with sim_key '''
		current_key_holder = HypothesisKeyHolder(self.mode_dict[Hnum], \
			self.fixed_param_dict[Hnum], starting_param_vals, sim_key)
		return(current_key_holder)
	def _set_up_original_data(self):
		''' Sets up hypothesis_testing_info for original data '''
		sim_key = 'original'
		H1_starting_param_vals = \
			self._get_starting_params('H1', self.sim_param_dict['original'])
		H1_key_holder = self._generate_hypothesis_key_holder('H1', sim_key, \
			H1_starting_param_vals)
		H0_starting_param_vals = self._get_starting_params('H0', \
			self.sim_param_dict['simulated'])
		H0_key_holder = self._generate_hypothesis_key_holder('H0', sim_key, \
			H0_starting_param_vals)
		# create sim_key_holder and get sim_key from
			# sim_key_organizer
		sim_key_holder = SimKeyHolder(self.mode_dict['H1'], H1_starting_param_vals)
		self.sim_key_organizer.set_key(sim_key_holder, sim_key)
		print('----- setting key ----')
		print(sim_key_holder.get_params())
		print(self.sim_key_organizer.key_df)
		# set up hypothesis testing info for original data
		self._set_up_hypothesis_testing_info('original', H1_key_holder, \
			H0_key_holder)
	def _get_original_MLE_param_vals(self, Hnum):
		'''
		Returns MLE parameters from MLE on original data for
		current_Hnum
		'''
		# set mode and fixed parameter in sim_parameters, get list
			# of parameter names
		self.sim_param_dict['original'].set_mode(self.mode_dict[Hnum], \
			self.output_id_prefix)
		self.sim_param_dict['original'].set_parameter(self.fixed_param_dict[Hnum])
		parameter_names = self.sim_param_dict['original'].get_complete_parameter_list()
		current_original_LL_df = \
			self.llr_calculator_dict['original'].get_LL(Hnum)
		current_original_MLE_param_vals = \
			current_original_LL_df.loc[0][parameter_names].tolist()
		return(np.array(current_original_MLE_param_vals))
	def _set_up_sim_data(self):
		''' Sets up hypothesis_testing_info for sim data '''
		# get MLE parameter values from LL_df, and use them as new
			# starting vals
		# only H1 data will be used for running simulations, but
			# starting_param_vals for both will be used as MLE starting
			# vals
		H1_starting_param_vals = self._get_original_MLE_param_vals('H1')
		# create sim_key_holder and get sim_key from sim_key_organizer
		sim_key_holder = SimKeyHolder(self.mode_dict['H1'], \
			H1_starting_param_vals)
		print('----- before getting key ... ----')
		print(self.sim_key_organizer.key_df)
		sim_key = self.sim_key_organizer.get_key(sim_key_holder)
		print('----- getting key ----')
		print(sim_key_holder.get_params())
		print(self.sim_key_organizer.key_df)
		print(sim_key)
		# create hypothesis_key_holder
		H1_key_holder = self._generate_hypothesis_key_holder('H1', sim_key, \
			H1_starting_param_vals)
		# now repeat for H0, using sim_key determined above
		H0_starting_param_vals = self._get_original_MLE_param_vals('H0')
		# create hypothesis_key_holder
		H0_key_holder = self._generate_hypothesis_key_holder('H0', sim_key, \
			H0_starting_param_vals)
		# pass hypothesis_key_holder to self.hypothesis_testing_info_dict['sim']
		self._set_up_hypothesis_testing_info('simulated', H1_key_holder, \
			H0_key_holder)
	def _run_sim(self, current_hypothesis_testing_info, \
		current_completeness_tracker, current_sim_parameters):
		'''
		Runs simulations based on H1 parameters in
		current_hypothesis_testing_info
		'''
		sim_runner = SimRunner(self.output_id_prefix, current_sim_parameters, \
			current_hypothesis_testing_info, self.cluster_parameters, \
			self.cluster_folders, self.sim_folders, \
			self.additional_code_run_keys, self.additional_code_run_values)
		sim_runner.submit_sim()
		sim_completeness = sim_runner.get_sim_completeness()
		current_completeness_tracker.switch_key_completeness('sim', \
			sim_completeness)
	def _run_LLR_calc(self, current_hypothesis_testing_info, \
		current_completeness_tracker, current_sim_parameters, \
		data_type):
		'''
		Estimates log likelihood ratios for given
		current_hypothesis_testing_info based on sim outputs
		'''
		self.llr_calculator_dict[data_type] = \
			LLRCalculator(self.output_id_prefix, current_sim_parameters, \
				current_hypothesis_testing_info, self.cluster_parameters, \
				self.cluster_folders, self.sim_folders, \
				self.additional_code_run_keys, self.additional_code_run_values)
		self.llr_calculator_dict[data_type].run_LLR()
		llr_completeness = \
			self.llr_calculator_dict[data_type].get_LLR_completeness()
		current_completeness_tracker.switch_key_completeness('LLR', \
			llr_completeness)
	def _run_sim_and_LLR(self, data_type):
		'''
		Based on current_hypothesis_testing_info, runs simulations (if
		necessary) and LLR calculation
		'''
		current_hypothesis_testing_info = \
			self.hypothesis_testing_info_dict[data_type]
		current_completeness_tracker = \
			self.completeness_tracker_dict[data_type]
		print('#################################')
		print('data type is ' + data_type)
		print('original pt num: ' + str(self.sim_param_dict['original'].current_profile_point_num))
		print('simulated pt num: ' + str(self.sim_param_dict['simulated'].current_profile_point_num))
		print('current pt num: ' + str(self.sim_param_dict[data_type].current_profile_point_num))
		current_sim_parameters = self.sim_param_dict[data_type]
		data_type_complete = current_completeness_tracker.get_completeness()
		if not data_type_complete:
			# run sim if it has not been completed
			if not current_completeness_tracker.get_key_completeness('sim'):
				self._run_sim(current_hypothesis_testing_info, \
					current_completeness_tracker, current_sim_parameters)
			# run LLR if sim has been completed
				# (need to re-check completion status in case sims got
					# completed in above run)
			if current_completeness_tracker.get_key_completeness('sim'):
				self._run_LLR_calc(current_hypothesis_testing_info, \
					current_completeness_tracker, current_sim_parameters, \
					data_type)
	def _calculate_cdf_val(self):
		'''
		Calculates proportion of sim data deviances above deviance of
		original data
		'''
		# sim_deviance_array and original_deviance_array are np arrays
		if len(self.original_deviance_array) > 1:
			raise RuntimeError('More than one deviance value returned for ' + \
				'original data')
		sim_deviance_array = \
			self.llr_calculator_dict['simulated'].get_deviances()
		if (len(sim_deviance_array) == 0) or \
			(len(self.original_deviance_array) == 0) or \
			np.all(np.isnan(self.original_deviance_array)) or \
			np.all(np.isnan(sim_deviance_array)):
			self.cdf_val = np.nan
		else:
			original_deviance = self.original_deviance_array[0]
			num_sim_deviances_above_original_deviance = \
				sum(np.greater(original_deviance, sim_deviance_array))
			num_sim_deviances = sum(~np.isnan(sim_deviance_array))
			self.cdf_val = \
				num_sim_deviances_above_original_deviance/num_sim_deviances
		self.cdf_val_calc_complete = True
		_write_fixed_pt_output(self.fixed_param_dict['H0'], \
			self.fixed_param_val_dict['H0'], self.cdf_val, self.output_file)
	def run_fixed_pt_cdf_val_estimation(self):
		'''
		Determine the cumulative density function val of the hypothesis
		comparison being performed
		If MLE on original data for one of the hypotheses fails, returns
		NaN
		'''
		self._set_up_original_data()
		self._run_sim_and_LLR('original')
		original_completeness = \
			self.completeness_tracker_dict['original'].get_completeness()
		if original_completeness:
			# check that deviance is non-NaN for original data; if it is
				# NaN, there's no point running sims
			self.original_deviance_array = \
				self.llr_calculator_dict['original'].get_deviances()
			if (len(self.original_deviance_array) == 0) or \
				np.all(np.isnan(self.original_deviance_array)):
				self.cdf_val = np.nan
				self.completeness_tracker_dict['simulated'].switch_key_completeness('sim', True)
				self.completeness_tracker_dict['simulated'].switch_key_completeness('LLR', True)
			else:
				self._set_up_sim_data()
				self._run_sim_and_LLR('simulated')
				simulated_completeness = \
					self.completeness_tracker_dict['simulated'].get_completeness()
				if simulated_completeness:
					self._calculate_cdf_val()
	def get_cdf_val(self):
		return(self.cdf_val)
	def get_cdf_val_completeness(self):
		return(self.cdf_val_calc_complete)

#####################################################################
# Pipeline for sim-LRT-based p-val

# Inputs:
# optimal_param_vector is fixed_param_mle
# fixed_parameter_indices include parameter for which LRT is being calculated (just as in asymptotic search)
# original_LL - likelihood value at MLE
# sim_number
# cutoff_pval
# parameter lower and upper bounds
# idealized_cutoff - x-val of asymptotic CI bound


class OneSidedSimProfiler(object):
	'''
	Calculates three sim-based cumulative density values on one side of
	MLE of a fixed parameter with a given mode
	To determine which point to calculate cumulative density vals at:
		1. 	Calculate cdf val at asymptotic CI value
		2. 	Taking cdf val from (1) as a cdf val on a normal curve,
			determine at which point on fixed parameter axis
			1/2*(target p-val) would be, and calculate the cdf
			val at that point
		3. 	a. 	If point (2) has a higher empirical cdf val than the
				target cdf val, interpolate between points (2) and (1)
				to id the best candidate for the target cdf val
			b. 	If point (2) has a lower empirical cdf val than the
				target cdf val, determine at which point on the fixed
				parameter axis 1/4*(target p-val) would be, taking the
				cdf val from (2) as a cdf val on a normal curve, and
				calculate the cdf val at that point
	'''
	def __init__(self, fixed_param_mle, asymptotic_CI_val, mode, fixed_param, \
		sim_parameters, sim_folders, additional_code_run_keys, \
		additional_code_run_values, output_id_prefix, cdf_bound, profile_pt_list, \
		cluster_folders, cluster_parameters, profile_path):
		self.mode_dict = {'H1': mode, 'H0': mode}
		self.fixed_param_dict = {'H1': 'unfixed', 'H0': fixed_param}
		self.sim_parameters = sim_parameters
		self.sim_folders = sim_folders
		self.cluster_folders = cluster_folders
		self.cluster_parameters = cluster_parameters
		self.profile_path = profile_path
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self.output_id_prefix = output_id_prefix
		self.asymptotic_CI_val = asymptotic_CI_val
		self.fixed_param_mle = fixed_param_mle
		self.cdf_bound = cdf_bound
		self.profile_pt_list = profile_pt_list
		self.profile_pt_num = len(self.profile_pt_list)
		self.fixed_point_df = \
			pd.DataFrame(\
				{'fixed_param_val_dict': [{'H1': self.fixed_param_mle, \
					'H0': np.nan}] * self.profile_pt_num, \
				'completeness': [False] * self.profile_pt_num, \
				'cdf_val': [np.nan] * self.profile_pt_num}, \
				index = self.profile_pt_list)
		self.side_completeness = False
#	def _convert_profile_pt_num(self, profile_pt):
#		'''
#		Converts profile_pt to an index that will be unique across both
#		sides of the sim-based profile
#		'''
#		if self.fixed_param_mle < self.asymptotic_CI_val:
#			# upper side of profi;le
#			profile_pt_converted = profile_pt + self.profile_points_per_side
#		elif self.fixed_param_mle > self.asymptotic_CI_val:
#			# lower side of profile, re-order point indices
#			profile_pt_converted = self.profile_points_per_side - profile_pt + 1
#		else:
#			raise ValueError('fixed_param_mle is neither greater than nor ' + \
#				'less than asymptotic_CI_val for sim-based CI estimation ' + \
#				'on ' + self.output_id_prefix + '; cannot detect profile side')
#		return(profile_pt_converted)
	def _fixed_param_val_finder(self, current_fixed_param_val, current_cdf_val, \
		target_cdf_val):
		'''
		Finds value of fixed_param that corresponds to target_cdf_val,
		assuming density function of fixed parameter probabilities is
		normal
		'''
		current_z_val = norm.ppf(current_cdf_val)
		z_val_scaler = \
			(current_fixed_param_val - self.fixed_param_mle) / current_z_val
		target_z_val = norm.ppf(target_cdf_val)
		target_fixed_param_val = \
			target_z_val * z_val_scaler + self.fixed_param_mle
		return(target_fixed_param_val)
	def _select_first_point(self):
		'''
		Assigns self.asymptotic_CI_val to fixed parameter val of 1st
		profile point
		'''
		current_fixed_param_val_dict = \
			copy.copy(self.fixed_point_df.loc[self.profile_pt_list[0]]['fixed_param_val_dict'])
		current_fixed_param_val_dict['H0'] = self.asymptotic_CI_val
		self.fixed_point_df.at[self.profile_pt_list[0], \
			'fixed_param_val_dict'] = current_fixed_param_val_dict
	def _select_second_point(self):
		'''
		Assumes cdf on first point is on normal distribution with center
		at mle val of fixed parameter, and assign second point to fixed
		parameter value where p-val is expected to be 1/2*(target p-val)
		'''
		first_fixed_param_val = \
			self.fixed_point_df.loc[self.profile_pt_list[0]]['fixed_param_val_dict']['H0']
		first_cdf_val = self.fixed_point_df.loc[self.profile_pt_list[0]]['cdf_val']
		pval_scaler = 0.5
		target_cdf_val = 1 - pval_scaler * (1 - self.cdf_bound)
		second_pt_fixed_param_val = \
			self._fixed_param_val_finder(first_fixed_param_val, first_cdf_val, \
				target_cdf_val)
		current_fixed_param_val_dict = \
			copy.copy(self.fixed_point_df.loc[self.profile_pt_list[1]]['fixed_param_val_dict'])
		current_fixed_param_val_dict['H0'] = second_pt_fixed_param_val
		self.fixed_point_df.at[self.profile_pt_list[1], 'fixed_param_val_dict'] = \
			current_fixed_param_val_dict
	def _select_third_point(self):
		'''
		If second fixed param val has higher cdf than cdf_bound, uses
		linear interpolation to get close to cdf_bound based on the two
		known values (or, in the unlikely case that cdf at both first
		and second point is above cdf_bound, extrapolates to cdf_bound;
		this will most likely occur due to sampling error in sim or poor
		estimation of asymptotic CI due to insufficient spacing of LL
		profile points)
		Otherwise, proceed as in second point selection, but re-scaling
		target p-val by 1/4 instead of 1/2
		'''
		second_fixed_param_val = \
			self.fixed_point_df.loc[self.profile_pt_list[1]]['fixed_param_val_dict']['H0']
		second_cdf_val = self.fixed_point_df.loc[self.profile_pt_list[1]]['cdf_val']
		if second_cdf_val > self.cdf_bound:
			# interpolate (or if necessary extrapolate) between first
				# and second fixed param val on x-axis to point
				# corresponding to cdf_bound y-val
			first_fixed_param_val = \
				self.fixed_point_df.loc[self.profile_pt_list[0]]['fixed_param_val_dict']['H0']
			first_cdf_val = self.fixed_point_df.loc[self.profile_pt_list[0]]['cdf_val']
			line_fit = \
				np.polyfit([first_fixed_param_val, second_fixed_param_val],
					[first_cdf_val, second_cdf_val], 1)
			interpolator = poly1d(line_fit)
			third_pt_fixed_param_val = interpolator(self.cdf_bound)
		else:
			pval_scaler = 0.25
			target_cdf_val = 1 - pval_scaler * (1 - self.cdf_bound)
			third_pt_fixed_param_val = \
				self._fixed_param_val_finder(second_fixed_param_val, \
					second_cdf_val, target_cdf_val)
		current_fixed_param_val_dict = \
			copy.copy(self.fixed_point_df.loc[self.profile_pt_list[2]]['fixed_param_val_dict'])
		current_fixed_param_val_dict['H0'] = third_pt_fixed_param_val
		self.fixed_point_df.at[self.profile_pt_list[2], 'fixed_param_val_dict'] = \
			current_fixed_param_val_dict
	def _run_fixed_pt_calc(self, profile_pt):
#		output_id_fixed_pt = '_'.join([self.output_id_prefix, 'fixed_pt_cdf'])
		output_file = generate_filename(self.profile_path, str(profile_pt), \
			self.output_id_prefix, self.fixed_param_dict['H0'], 'data')
		current_fixed_param_val_dict = \
			self.fixed_point_df.loc[profile_pt]['fixed_param_val_dict']
		fixed_point_pval_calc = FixedPointCDFvalCalculator(self.mode_dict, \
			self.fixed_param_dict, current_fixed_param_val_dict, \
			self.sim_parameters, self.sim_folders, \
			self.additional_code_run_keys, self.additional_code_run_values, \
			self.output_id_prefix, output_file, self.cluster_folders, \
			self.cluster_parameters)
		fixed_point_pval_calc.run_fixed_pt_cdf_val_estimation()
		current_fixed_pt_complete = \
			fixed_point_pval_calc.get_cdf_val_completeness()
		self.fixed_point_df.at[profile_pt, 'completeness'] = \
			current_fixed_pt_complete
		if current_fixed_pt_complete:
			current_cdf_val = fixed_point_pval_calc.get_cdf_val()
			self.fixed_point_df.at[profile_pt, 'cdf_val'] = \
				current_cdf_val
		return(current_fixed_pt_complete)
	def get_sim_profile_side_completeness(self):
		return(self.side_completeness)
	def run_sim_profiler(self):
		'''
		Determines three fixed pts at which to run sims, and runs sims
		and mle on them
		'''
		self._select_first_point()
		first_fixed_pt_complete = self._run_fixed_pt_calc(self.profile_pt_list[0])
		if first_fixed_pt_complete:
			self._select_second_point()
			second_fixed_pt_complete = self._run_fixed_pt_calc(self.profile_pt_list[1])
			if second_fixed_pt_complete:
				self._select_third_point()
				third_fixed_pt_complete = self._run_fixed_pt_calc(self.profile_pt_list[2])
				if third_fixed_pt_complete:
					self.side_completeness = True

class TwoSidedProfiler(object):
	'''
	Runs OneSidedSimProfiler on both upper and lower sides of mle value
	of fixed parameter
	'''
	def __init__(self, fixed_param_mle, asymptotic_CI_dict, mode, fixed_param, \
		sim_parameters, sim_folders, additional_code_run_keys, \
		additional_code_run_values, output_id_prefix, cdf_bound, cluster_folders, \
		cluster_parameters):
		self.mode = mode
		self.fixed_param = fixed_param
		self.sim_parameters = sim_parameters
		self.sim_folders = sim_folders
		self.cluster_folders = cluster_folders
		self.cluster_parameters = cluster_parameters
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self.output_id_prefix = output_id_prefix
		self.asymptotic_CI_dict = asymptotic_CI_dict
		self.fixed_param_mle = fixed_param_mle
		self.cdf_bound = cdf_bound
		self.profile_points_per_side = 3
		self.profile_pt_list = \
			np.array(range(0, (self.profile_points_per_side * 2 + 1)))
		self.profile_path = sim_folders.get_path('sim_profile_fixed_pt_folder')
		self.CI_sides = ['lower', 'upper']
		self.completeness_tracker = \
			cluster_wrangler.cluster_functions.CompletenessTracker(self.CI_sides)
		self._create_profile_pt_list_dict()
		self.completefile = \
			os.path.join(cluster_folders.get_path('completefile_path'), \
				'_'.join(['param_profile_', sim_parameters.output_identifier, \
					fixed_param, 'completefile.txt']))
	def _create_profile_pt_list_dict(self):
		lower_profile_pts = \
			self.profile_pt_list[range(1, (self.profile_points_per_side + 1))]
		upper_profile_pts = \
			self.profile_pt_list[range((self.profile_points_per_side + 1), \
				(2 * self.profile_points_per_side + 1))]
		self.profile_pt_list_dict = {'lower': lower_profile_pts, \
			'upper': upper_profile_pts}
	def get_completefile(self):
		return(self.completefile)
	def run_profiler(self):
		# write output file for 0th pt
		profile_pt = 0
		output_file = generate_filename(self.profile_path, str(profile_pt), \
			self.output_id_prefix, self.fixed_param, 'data')
		_write_fixed_pt_output(self.fixed_param, self.fixed_param_mle, profile_pt, \
			output_file)
		# run through sides
		for current_profile_side in self.CI_sides:
			current_profile_pt_list = \
				self.profile_pt_list_dict[current_profile_side]
			asymptotic_CI_val = self.asymptotic_CI_dict[current_profile_side]
			current_side_sim_profiler = \
				OneSidedSimProfiler(self.fixed_param_mle, asymptotic_CI_val, \
					self.mode, self.fixed_param, self.sim_parameters, \
					self.sim_folders, self.additional_code_run_keys, \
					self.additional_code_run_values, self.output_id_prefix, \
					self.cdf_bound, current_profile_pt_list, \
					self.cluster_folders, self.cluster_parameters, \
					self.profile_path)
			current_side_sim_profiler.run_sim_profiler()
			current_side_completeness = \
				current_side_sim_profiler.get_sim_profile_side_completeness()
			self.completeness_tracker.switch_key_completeness(current_profile_side, \
				current_side_completeness)
		if self.completeness_tracker.get_completeness():
			open(self.completefile,'a').close()



			
# For comparing models, don't run FixedPointIdentifier, just run FixedPointCDFvalCalculator on the two models
#####################################################################

def _csv_list_converter(instr):
	''' Converts values read in from csv cell to np array '''
	# solution from https://stackoverflow.com/questions/50934301/numpy-array-changes-to-string-when-writing-to-file
	return np.fromstring(instr[1:-1],sep=' ')

def _generate_sim_file_label(sim_key, input_datafile_name):
	sim_file_label = \
		mle_filenaming_functions.generate_file_label('sim', str(sim_key), \
			input_datafile_name)
	return(sim_file_label)

def _generate_sim_filename(within_batch_counter_call, sim_file_label):
	sim_file_label_split_by_dot = sim_file_label.split('.')
	sim_file_extension = sim_file_label_split_by_dot[-1]
	sim_file_label_no_extension = '.'.join(sim_file_label_split_by_dot[0:-1])
	new_sim_file_label = '_'.join([sim_file_label_no_extension, \
		within_batch_counter_call])
	sim_filename = '.'.join([new_sim_file_label, sim_file_extension])
	return(sim_filename)

def _write_fixed_pt_output(fixed_param, fixed_param_val, current_cdf_val, \
	current_output_file):
		'''
		Writes mlestimation.LLProfile-readable output containing
		current value of fixed_param and current cdf val
		'''
		output_df = pd.DataFrame({fixed_param: [fixed_param_val], \
			'cdf_vals': [current_cdf_val]})
		output_df.to_csv(path_or_buf = current_output_file, \
			index = False)

def generate_sim_based_profile_pts(mode, sim_parameters, sim_folders, \
	additional_code_run_keys, additional_code_run_values, output_id_prefix, \
	combined_results, cluster_folders, cluster_parameters):
		cdf_bound = 1 - sim_parameters.current_CI_pval
		for fixed_param in sim_parameters.current_sim_CI_parameters:
			asymptotic_CI_complete = \
				combined_results.check_step_completeness('asymptotic_CIs')
			if asymptotic_CI_complete:
				fixed_param_mle = \
					combined_results.get_param_mle_val(fixed_param)
				asymptotic_CI_dict = \
					combined_results.get_CI_dict(fixed_param, 'asymptotic')
				current_param_profiler = TwoSidedProfiler(fixed_param_mle, \
					asymptotic_CI_dict, mode, fixed_param, sim_parameters, \
					sim_folders, additional_code_run_keys, \
					additional_code_run_values, output_id_prefix, cdf_bound, \
					cluster_folders, cluster_parameters)
				current_param_profiler.run_profiler()
				profile_completefile = \
					current_param_profiler.get_completefile()
				sim_parameters.set_parameter(fixed_param)
				sim_parameters.update_parameter_completeness(profile_completefile)


