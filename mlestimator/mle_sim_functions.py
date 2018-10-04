#!/usr/bin/python

# Contains functions necessary for running simulations and combining their results

import pandas as pd
import numpy as np
import warnings as war
import mle_functions
from cluster_wrangler.cluster_functions import CompletenessTracker
import copy
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

# have a list of simulation key files, by mode and fixed sim parameters; each will have a unique key number
# before starting new set of sims, see if key for such a sim already exists (e.g. look through ones with correct sim mode);
	# if it does, use the key number corresponding to that sim for trackfiles, completefiles, etc (i.e. in sim output_id)

#####################################################################
class SimFolders(object):
	pass()
	# needs to have a get_experiment_path method, as MLEFolders does
	# needs to have all the same folder types as MLEFolders does
	# also needs:
	#	current_sim_folder

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
		super(mle_functions.MLEParameters, self).__init__(parameter_list)
		self.sim_CI_parameters_by_mode = \
			parameter_list["sim_CI_parameters_by_mode"]
		self.simulation_repeats_by_mode = \
			parameter_list["simulation_repeats_by_mode"]
		# reset profile_points_by_mode to number of sim repeats by mode
		self.profile_points_by_mode = self.simulation_repeats_by_mode
		# set default parameter_to_select for looping over parameters
		self.parameter_to_select = 'unfixed'
	def set_mode(self, mode_name, output_identifier):
		mode_idx = self.mode_list.index(mode_name)
		self.current_sim_rep_num = self.simulation_repeats_by_mode[mode_idx]
		self.current_sim_CI_parameters = \
			self.sim_CI_parameters_by_mode[mode_idx]
		super(mle_functions.MLEParameters, self).set_mode(mode_name, \
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
		self.point_numbers_to_loop_over = numpy.array(self.current_sim_rep_num)
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
		self.simulation_repeats_by_mode = \
			[self.simulation_repeats_by_mode[mode_idx]]
		self.sim_CI_parameters_by_mode = \
			[self.sim_CI_parameters_by_mode[mode_idx]]
		# reset mode_completeness_tracker
		self.mode_completeness_tracker = CompletenessTracker(self.mode_list)
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
	def respecify_for_hypothesis_testing(self, mode_name, fixed_param, starting_vals):
		self._select_mode_to_keep(mode_name)
		self._change_starting_param_vals(mode_name, starting_vals)
		if fixed_param:
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
	'''
	def __init__(self, key_file):
		self.key_file = key_file
		self._read_key_file()
	def _read_key_file(self):
		# if self.key_file exists, reads it in;
		# otherwise, creates one
		if os.path.isfile(self.LL_file):
			try:
				self.key_df = pd.read_csv(filepath_or_buffer=self.key_file)
			except pd.io.common.EmptyDataError:
				# if file is empty, just create an empty df for self.LL_df
				self._create_key_df()
		else:
			self._create_key_df()
	def _create_key_df(self):
		''' Creates new key_df '''
		self.key_df = pandas.DataFrame()
	def _find_matching_key(self, key_holder):
		'''
		Takes in a KeyHolder object and identifies the row in
		self.key_df that matches it, if any do
		'''
		current_keys = self.key_df.index.values.tolist()
		row_match_boolean_list = [key_holder.get_params().equals(self.key_df.loc[i]) for i in current_keys]
		matching_keys = self.key_df.index[row_match_boolean_list].tolist()
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
			self.key_df = \
				self.key_df.append(key_holder.get_params(), sort = False, \
					ignore_index = True)
			current_key = self._find_matching_key(key_holder)
			self.write_key_df()
		return(current_key)
	def set_key(self, key_holder, key_idx):
		''' Sets index key_idx in self.key_df to key_holder '''
		self.key_df.loc[key_idx] = key_holder
	def write_key_df(self):
		''' Writes self.key_df to self.key_file '''
		self.key_df.to_csv(path_or_buf = self.key_file, index = True)

class HypothesisTestingInfo(object):
	'''
	Holds information for performing MLE on a pair of null and
	alternative hypotheses in a pandas DF
	Hypotheses are 'H0' (null) and 'H1' (alternative)
	Need to first run set_hypothesis_parameters method with both H0 and
	H1, then run set_up_sim_key method
	'''
	def __init__(self):
		self.hypotheses = ['H0', 'H1']
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
		be used for hypothesis Hnum ('H0' or 'H1')
		'''
		if Hnum in self.hypotheses:
			self.hypothesis_info.loc[Hnum] = 
				hypothesis_key_holder.get_params()
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
		if len(unique_sim_key_vals) not 1:
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

# ID 3 points at which to find LRs--maybe, 2*(idealized_cutoff-optimal_param_vector(current_param))?
#	- one option is to ID these points stepwise: i.e. first try asymptotic bound, then go from there
#		- example flow:
#			a. measure LR p-val at asymptotic value
#			b. use this value and x-position of MLE to approximate a normal(? or t?) distribution,
#				and find second x-val to try based on where desired p-val would be on that distribution
#			c. Linear interpolate between the p-vals and x-vals measured so far to ID a third x-position to try
#			d. Quad-fit the 3 values
#		The problem with this approach is if the initial values randomly happen to go in wrong directions
#		(i.e. p-val increases, not decreases, with x-distance from optimum), linear interpolation will take
#		next value further from target
#		However, this could either be explicitly dealt with, or it could be assumed that if this happens, all
#		values are close enough together and to target for exact x-values to not matter much
#		Also need to deal with what to do when p-val at both points is identical...

class LLRCalculator(object):
	'''
	For self.data (which can be real or come from a sim) and fixed_param_val:
		1. Calculates H1: 'unfixed' LL (all unknown parameters freely
			estimated)
		2. Calculates H0: LL at fixed_param_val
		3. Calculates LL(H0)-LL(H1)
	'''
	def __init__(self, output_id_sim, sim_key, \
		sim_parameters, hypothesis_testing_info, cluster_parameters, \
		cluster_folders, sim_folders, additional_code_run_keys, \
		additional_code_run_values):
		self.output_id_sim = output_id_sim
			# output_id_sim needs to include sim key but not sim
				# number; fixed parameter will be added to filenames
				# by SimParameters, and sim number will be added in MLE
		self.sim_key = sim_key
		self.sim_parameters = copy.deepcopy(sim_parameters)
		self.sim_folders = sim_folders
		self.cluster_parameters = cluster_parameters
		self.cluster_folders = cluster_folders
		self.mle_datafile_path = sim_folders.get_path('current_output_subfolder')
		self.LL_list_folder = sim_folders.get_path('LL_list_path')
		### ??? ###
		# Each hypothesis will need a unique title, likely incorporating fixed param, fixed param val, and mode setting
		# LLR file should include, in order, titles for these parameters?
		### ??? ###
		self.hypothesis_testing_info = hypothesis_testing_info
		#self.hypothesis_list = self.hypothesis_testing_info.get_hypotheses()
		# Hard-code hypotheses to be H0 and H1
		self.hypothesis_list = ['H0','H1']
		self._generate_LLR_filename()
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self._set_up_sim_parameters()
		self.LL_list_dict = dict()
		# create CompletenessTracker object to track whether each LL
			# list is complete
		self.LL_list_completeness_tracker = CompletenessTracker(hypothesis_list)
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
	def _set_up_sim_parameters(self):
		'''
		Uses info in self.hypothesis_testing_info to respecify
		SimParameters object for hypothesis testing
		'''
		# It's important to keep sim_parameters object as attribute of
			# LRCalculator so that the objects can be modified to keep
			# track of modes and parameters
		self.sim_param_dict = dict()
		for current_Hnum in self.hypothesis_list:
			# respecify SimParameter object for each hypothesis
			current_sim_param = copy.deepcopy(self.sim_parameters)
			current_mode = self.hypothesis_testing_info.get_mode(current_Hnum)
			current_fixed_param = \
				self.hypothesis_testing_info.get_fixed_param(current_Hnum)
			current_starting_param_vals = \
				self.hypothesis_testing_info.get_starting_param_vals(current_Hnum)
			current_sim_params.respecify_for_hypothesis_testing(current_mode, \
				current_fixed_param, current_starting_param_vals)
			self.sim_param_dict[current_Hnum] = current_sim_params
	def _run_MLE(self, current_sim_parameters):
		'''
		Runs MLE for getting LL values to be used in hypothesis testing
		'''
		# run MLE; completeness of current_sim_parameters will
			# automatically be updated
		include_unfixed_parameter = False
		mle_functions.run_MLE(current_sim_parameters, self.cluster_parameters, \
			self.cluster_folders, self.sim_folders, \
			self.additional_code_run_keys, self.additional_code_run_values, \
			include_unfixed_parameter)
	def _compile_LL_list(self, Hnum, current_sim_parameters):
		'''
		Compiles and writes LL_list for each hypothesis
		'''
		ll_list = LLHolder(current_sim_parameters, \
			self.mle_datafile_path, self.LL_list_folder)
		ll_list.run_LL_list_compilation()
			# compiles ll_list and writes it to a file
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
		current_mode = self.hypothesis_testing_info.get_mode(Hnum)
		current_h_key = self.hypothesis_testing_info.get_hypothesis_key(Hnum)
		current_output_id = '_'.join([self.output_id_sim, Hnum, str(current_h_key)])
		current_sim_parameters = self.sim_param_dict[Hnum]
		current_sim_parameters.set_mode(current_mode, current_output_id)
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
			if current_fixed_param = 'unfixed':
				abbreviated_LL_df_dict[current_Hnum] = \
					current_full_df[['LL', 'point_num']]
			else:
				abbreviated_LL_df_dict[current_Hnum] = \
					current_full_df[['LL', 'point_num', current_fixed_param]]
		self.LLR_df = pd.merge(abbreviated_LL_df_dict['H0'], \
			abbreviated_LL_df_dict['H1'], how='outer', left_on = 'point_num', \
			right_on = 'point_num', suffixes = (['_H0','_H1']))
		self.LLR_df['LLR'] = self.LLR_df['LL_H0'] - self.LLR_df['LL_H1']
		self.LLR_df.to_csv(path_or_buf = self.LLR_file)
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
	def get_LLRs(self):
		return(self.LLR_df['LLR'])
	def get_LL(self, Hnum):
		return(self.LL_list_dict[Hnum])




######## ??? #######
# Before running LLRCalculator on original data, a copy of it needs to
# be created in the same place where simulation folders are stored;
# The 'unfixed' output file for this data should be placed in this folder
######## ??? #######

### sim_key == 'original' for original data

class FixedPointPvalCalculator(object):
	'''
	Performs simulations and MLEs to calculate p-val given parameters in
	hypothesis_testing_info
	To calculate p-val at a given point:
		1. 	a.	Calculate H1: the hypothesis with the higher number of
				degrees of freedom; for comparisons of hypothesis in one
				of which a parameter is set to a predetermined value,
				this corresponds to the 'unfixed' LL (all unknown
				parameters freely estimated)
			b.	Calculate H0: the hypothesis with the lower number of
				degrees of freedom; LL at fixed point in parameter space
			c. 	'True' log likelihood ratio 'LLR' is LL(H0)-LL(H1)
				(actually to make the p-vals easier to understand, it's
				easier to work with the negative of LLR values)
			d. 	MLE parameter values for all parameters estimated in (b)
				will be used as sim starting values!
		2. 	Run sim_number simulations using starting parameter values
			from (1d)
		3. Repeat and record (1) but with sim data rather than real data
	'''
	def __init__(self, H0_mode, H1_mode, H0_fixed_param, \
		H1_fixed_param, H0_fixed_param_val, H1_fixed_param_val, \
		sim_parameters, hypothesis_key_organizer, sim_folders,
		sim_key_organizer):
		self.sim_folders = copy.deepcopy(sim_folders)
		self.hypothesis_key_organizer = hypothesis_key_organizer
		self.sim_key_organizer = sim_key_organizer
		self.hypotheses = ['H0','H1']
		self.sim_parameters = copy.deepcopy(sim_parameters)
		self._set_up_original_data_hypothesis(H0_mode, H1_mode, \
			H0_fixed_param, H1_fixed_param, H0_fixed_param_vals, \
			H1_fixed_param_vals)
		self.original_data_completeness_tracker = CompletenessTracker(['sim','LLR'])
		self.sim_data_completeness_tracker = CompletenessTracker(['sim','LLR'])
		self.llr_calculator_dict = dict()
		self.hypothesis_testing_info_dict = dict()
	def _get_starting_params(self, current_mode, current_fixed_param, \
		current_fixed_param_vals):
		self.sim_parameters.set_mode(current_mode)
		self.sim_parameters.set_parameter(current_fixed_param)
		parameter_names = self.sim_parameters.get_complete_parameter_list()
		starting_vals = \
			copy.copy(self.sim_parameters.current_start_parameter_val_list())
		if current_fixed_param not 'unfixed':
			fixed_param_idx = parameter_names.index(current_fixed_param)
			starting_vals[fixed_param_idx] = current_fixed_param_vals
		return(starting_vals)
	def _set_up_hypothesis_testing_info(self, data_type, H0_key_holder, \
		H1_key_holder):
		'''
		Sets up HypothesisTestingInfo object with H0_key_holder and
		H1_key_holder in self.hypothesis_testing_info_dict[data_type]
		'''
		current_hypothesis_testing_info = HypothesisTestingInfo()
		current_hypothesis_testing_info.set_hypothesis_parameters('H0', \
			H0_key_holder, self.hypothesis_key_organizer)
		current_hypothesis_testing_info.set_hypothesis_parameters('H1', \
			H1_key_holder, self.hypothesis_key_organizer)
		current_hypothesis_testing_info..set_up_sim_key()
		self.hypothesis_testing_info_dict[data_type] = \
			current_hypothesis_testing_info
	def _set_up_original_data(self, H0_mode, H1_mode, \
			H0_fixed_param, H1_fixed_param, H0_fixed_param_vals, \
			H1_fixed_param_vals):
		''' Sets up hypothesis_testing_info for original data '''
		sim_key = 'original'
		H0_starting_param_vals = self._get_starting_params(H0_mode, H0_fixed_param, \
			H0_fixed_param_vals)
		H0_key_holder = self.HypothesisKeyHolder(H0_mode, H0_fixed_param, \
			H0_starting_param_vals, sim_key)
		H1_starting_param_vals = self._get_starting_params(H1_mode, H1_fixed_param, \
			H1_fixed_param_vals)
		H1_key_holder = self.HypothesisKeyHolder(H1_mode, H1_fixed_param, \
			H1_starting_param_vals, sim_key)
		# create sim_key_holder and get sim_key from
			# sim_key_organizer
		sim_key_holder = SimKeyHolder(H0_mode, H0_starting_param_vals)
		self.sim_key_organizer.set_key(sim_key_holder, sim_key)
		# set up hypothesis testing info for original data
		self._set_up_hypothesis_testing_info('original', H0_key_holder, \
			H1_key_holder)
	def _get_original_MLE_param_vals(self, current_Hnum, current_mode, \
		current_fixed_param):
		'''
		Returns MLE parameters from MLE on original data for
		current_Hnum
		'''
		# set mode and fixed parameter in sim_parameters, get list
			# of parameter names
		self.sim_parameters.set_mode(current_mode)
		self.sim_parameters.set_parameter(current_fixed_param)
		parameter_names = self.sim_parameters.get_complete_parameter_list()
		current_original_LL_df = \
			self.llr_calculator_dict['original'].get_LL(current_Hnum)
		current_original_MLE_param_vals = \
			current_original_LL_df.loc[0][parameter_names].tolist()
		return(current_original_MLE_param_vals)
	def _set_up_sim_data(self):
		self.sim_hypothesis_testing_info = HypothesisTestingInfo()
		# first, get H0 results from original data MLE, since simulation
			# will be based on MLE parameter values for those
		H0_mode = self.original_hypothesis_testing_info.get_mode('H0')
		H0_fixed_param = \
			self.original_hypothesis_testing_info.get_fixed_param('H0')
		# get MLE parameter values from LL_df, and use them as new
			# starting vals
		H0_starting_param_vals = \
			self._get_original_MLE_param_vals('H0', H0_mode, H0_fixed_param)
		# create sim_key_holder and get sim_key from
			# sim_key_organizer
		sim_key_holder = SimKeyHolder(H0_mode, H0_starting_param_vals)
		sim_key = self.sim_key_organizer.get_key(sim_key_holder)
		# create hypothesis_key_holder
		H0_key_holder = self.HypothesisKeyHolder(H0_mode, H0_fixed_param, \
			H0_starting_param_vals, sim_key)
		# now repeat for H1, using sim_key determined above
		H1_mode = self.original_hypothesis_testing_info.get_mode('H1')
		H1_fixed_param = \
			self.original_hypothesis_testing_info.get_fixed_param('H1')
		H1_starting_param_vals = \
			self._get_original_MLE_param_vals('H1', H1_mode, H1_fixed_param)
		# create hypothesis_key_holder
		H1_key_holder = self.HypothesisKeyHolder(H1_mode, H1_fixed_param, \
			H1_starting_param_vals, sim_key)
		# pass hypothesis_key_holder to self.hypothesis_testing_info_dict['sim']
		self._set_up_hypothesis_testing_info('sim', H0_key_holder, \
			H1_key_holder)
	def _run_sim(self, current_hypothesis_testing_info, \
		current_completeness_tracker):
		'''
		Runs simulations based on parameters in
		current_hypothesis_testing_info
		'''
	def _run_LLR_calc(self, current_hypothesis_testing_info, \
		current_completeness_tracker):
		'''
		Estimates log likelihood ratios for given
		current_hypothesis_testing_info based on sim outputs
		'''
	def _run_sim_and_LLR(self, current_hypothesis_testing_info, \
		current_completeness_tracker):
		'''
		Based on current_hypothesis_testing_info, runs simulations (if
		necessary) and LLR calculation
		'''
		self.original_llr_calculator = 

		### Need to allow for option that original data MLE (or LLR) ran but did not successfully complete
			# LLR value may be NaN or missing
		### Need to change the way data inputs are handled in mle_functions--make separate dict of data input files (read in from setup_file) that gets added to submission keys and values during mle batch submission
		return(LLR_list)
	def run_fixed_pt_pval_estimation(self):
		'''
		Determine the p-val of the hypothesis comparison being performed
		'''

		LLR_list_original = \
			self._run_sim_and_LLR(self.original_hypothesis_testing_info, \
			self.original_data_completeness_tracker)
		new_starting_vals = # read H0 MLE output file, extract values corresponding to parameter names of H0
			# does this work for comparing different modes????? think about this.
		self.sim_hypothesis_testing_info = copy.deepcopy(self.original_hypothesis_testing_info)




class FixedPointIdentifier(object):
	'''
	To determine which point to calculate p-vals at:
		1. 	Calculate p-val at asymptotic CI value
		2. 	Taking p-val from (1) as a p-val on a normal curve,
			determine at which point on fixed parameter axis
			1/2*(target p-val) would be, and calculate the p-val at
			that point
		3. 	a. 	If point (2) has a lower empirical p-val than the target
				p-val, interpolate between points (2) and (1) to id the
				best candidate for the target p-val
			b. 	If point (2) has a higher empirical p-val than the
				target p-val, determine at which point on the fixed
				parameter axis 1/4*(target p-val) would be, taking the
				p-val from (2) as a p-val on a normal curve, and
				calculate the p-val at that point
	'''
			
# For comparing models, don't run FixedPointIdentifier, just run FixedPointPvalCalculator on the two models
#####################################################################
