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
	Stores and compares individual simulation keys as a pandas Series
	'''
	def __init__(self, mode, fixed_params, fixed_param_vals):
		'''
		Sorts fixed_params alphabetically, and rearranges
		fixed_param_vals according to new fixed_param order, to
		enable easy comparison to sim keys that may have been entered
		in a different order
		'''
		sorted_indices = sorted(range(len(fixed_params)), \
			key=lambda k: fixed_params[k])
		fixed_params_sorted = [fixed_params[i] for i in sorted_indices]
		fixed_param_vals_sorted = np.array([fixed_param_vals[i] \
			for i in sorted_indices])
		self.params = \
			pd.Series([mode, fixed_params_sorted, fixed_param_vals_sorted],
				index = ['mode', 'fixed_params', 'fixed_param_vals'])
	def get_params(self):
		return(self.params)

class SimKeyOrganizer(object):
	'''
	Organizes simulation keys and checks if sim with current
	properties already exists; index in self.sim_key_df serves as
	simulation key
	'''
	def __init__(self, sim_key_file):
		self.sim_key_file = sim_key_file
		self._read_sim_key_file()
	def _read_sim_key_file(self):
		# if self.sim_key_file exists, reads it in;
		# otherwise, creates one
		if os.path.isfile(self.LL_file):
			try:
				self.sim_key_df = pd.read_csv(filepath_or_buffer=self.LL_file)
			except pd.io.common.EmptyDataError:
				# if file is empty, just create an empty df for self.LL_df
				self._create_sim_key_df()
		else:
			self._create_sim_key_df()
	def _create_sim_key_df(self):
		''' Creates new sim_key_df '''
		self.sim_key_df = pandas.DataFrame(columns = \
			['sim_mode', 'fixed_params', 'fixed_param_vals'])
	def _find_matching_key(self, sim_key_holder):
		'''
		Takes in a SimKeyHolder object and identifies the row in
		self.sim_key_df that matches it, if any do
		'''
		current_keys = self.sim_key_df.index.values.tolist()
		row_match_boolean_list = [sim_key_holder.get_params().equals(self.sim_key_df.loc[i]) for i in current_keys]
		matching_keys = self.sim_key_df.index[row_match_boolean_list].tolist()
		if matching_keys:
			key_val = matching_keys[0]
		else:
			key_val = None
		return(key_val)
	def get_sim_key(self, sim_key_holder):
		'''
		Takes in a SimKeyHolder object and determines whether it's already
		in self.sim_key_df
		If not, adds it, and writes updated file to self.sim_key_file
		Either way, returns the sim_key for this set of SimKeyHolder
		(i.e. the index of the line in sim_key_df holding these params)
		'''
		current_key = self._find_matching_key(sim_key_holder)
		if not current_key:
			self.sim_key_df = \
				self.sim_key_df.append(sim_key_holder.get_params(), sort = False, \
					ignore_index = True)
			current_key = self._find_matching_key(sim_key_holder)
			self.write_key_df()
		return(current_key)
	def write_key_df(self):
		''' Writes self.sim_key_df to self.sim_key_file '''
		self.sim_key_df.to_csv(path_or_buf = self.sim_key_file, index = True)

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

class HypothesisTestingInfo(object):
	'''
	Holds information for performing MLE on a pair of null and alternative
	hypotheses in a pandas DF
	Hypotheses are 'H0' (null) and 'H1' (alternative)
	'''
	def __init__(self):
		self.hypotheses = ['H0', 'H1']
		self.hypothesis_info = \
			pd.DataFrame(columns = \
				['mode', 'fixed_param', 'starting_param_vals'],
				index = self.hypotheses)
	def set_hypothesis_parameters(self, Hnum, mode, fixed_param, \
		starting_param_vals):
		'''
		Sets the mode, fixed_param (None if no parameter is fixed) and
		starting_param_vals (None if default should be used) that will
		be used for hypothesis Hnum ('H0' or 'H1')
		'''
		if Hnum in self.hypotheses:
			self.hypothesis_info.loc[Hnum] = \
				[mode, fixed_param, starting_param_vals]
		else:
			raise ValueError("invalid hypothesis (Hnum): " + Hnum + \
				'; Hnum may only be one of the following: ' + \
				', '.join(self.hypotheses))
	def _get_hypothesis_info(self, Hnum, desired_var):
		if Hnum in self.hypotheses:
			return(self.hypothesis_info.loc[Hnum][desired_var])
		else:
			raise ValueError("invalid hypothesis (Hnum): " + Hnum + \
				'; Hnum may only be one of the following: ' + \
				', '.join(self.hypotheses))
	def get_mode(self, Hnum):
		return(self._get_hypothesis_info(Hnum, 'mode'))
	def get_fixed_param(self, Hnum):
		return(self._get_hypothesis_info(Hnum, 'fixed_param'))
	def get_starting_param_vals(self, Hnum):
		return(self._get_hypothesis_info(Hnum, 'starting_param_vals'))
	def get_hypotheses(self):
		return(self.hypotheses)

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
		self.hypothesis_testing_info = hypothesis_testing_info
		self.additional_code_run_keys = additional_code_run_keys
		self.additional_code_run_values = additional_code_run_values
		self._set_up_sim_parameters()
		self.LL_list_dict = dict()
	def _set_up_sim_parameters(self):
		'''
		Uses info in self.hypothesis_testing_info to respecify
		SimParameters object for hypothesis testing
		'''
		# It's important to keep sim_parameters object as attribute of
			# LRCalculator so that the objects can be modified to keep
			# track of modes and parameters
		self.sim_param_dict = dict()
		self.completeness_tracker_dict = dict()
		for current_Hnum in self.hypothesis_testing_info.get_hypotheses():
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
			# create CompletenessTracker object for each hypothesis
			self.completeness_tracker_dict[current_Hnum] = \
				CompletenessTracker(['MLE', 'LL_list'])
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
	def _compile_LL_list(self, current_sim_parameters):
		'''
		Compiles and writes LL_list for each hypothesis
		'''
		ll_list = LLHolder(current_sim_parameters, \
			self.mle_datafile_path, self.LL_list_folder)
		ll_profile.run_LL_list_compilation()
		current_ll_df = ll_profile.get_LL_df()
		return(current_ll_df)
	def _find_LLs(self, Hnum):
		'''
		Respecifies sim_parameters for current hypothesis being tested;
		runs MLE and, once that is complete, compiles LL_list for
		hypothesis
		'''
		current_output_id = output_id_sim + '_' + Hnum
		current_mode = self.hypothesis_testing_info.get_mode(Hnum)
		current_sim_parameters = self.sim_param_dict[Hnum]
		current_sim_parameters.set_mode(current_mode, current_output_id)
		self._run_MLE(current_sim_parameters)
		mle_complete = current_sim_parameters.check_completeness_within_mode()
		if mle_complete:
			self.completeness_tracker_dict[Hnum].switch_key_completeness('MLE', \
				True)
			self.LL_list_dict[Hnum] = \
				self._compile_LL_list(current_sim_parameters)
			self.completeness_tracker_dict[Hnum].switch_key_completeness('LL_list', \
				True)
	def _calculate_LLR(self):
		# check for LL file
			# if not there, check for MLE completefile
				# if not there, run MLE
				# if there, compile LL

######## ??? #######
# Before running LLRCalculator on original data, a copy of it needs to
# be created in the same place where simulation folders are stored;
# The 'unfixed' output file for this data should be placed in this folder
######## ??? #######
class FixedPointPvalCalculator(object):
	'''
	To calculate p-val at a given point:
		1. 	a.	Calculate H1: 'unfixed' LL (all unknown parameters
				freely estimated)
			b.	Calculate H0: LL at fixed point in parameter space
			c. 	'True' log likelihood ratio 'LLR' is LL(H0)-LL(H1)
				(actually to make the p-vals easier to understand, it's
				easier to work with the negative of LLR values)
			d. 	MLE parameter values for all parameters estimated in (b)
				will be used as sim starting values!
		2. 	Run sim_number simulations using starting parameter values
			from (1d)
		3. Repeat and record (1) but with sim data rather than real data
	'''
	def __init__(self, fixed_param, mode):
		pass()

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
			
#####################################################################
