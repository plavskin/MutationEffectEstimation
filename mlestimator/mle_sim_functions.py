#!/usr/bin/python

# Contains functions necessary for running simulations and combining their results

import pandas as pd
import numpy as np
import warnings as war
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
		self.hypotheses = ['H0','H1']
		self.hypothesis_info = \
			pd.DataFrame(columns = \
				['mode', 'fixed_param', 'starting_param_vals'],
				index = self.hypotheses)
	def set_hypothesis_parameters(self, Hnum, mode, fixed_param, \
		starting_param_vals):
		'''
		Sets the mode, fixed_param ('unfixed' if none) and
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
	def get_hypothesis_info(self, Hnum, desired_var):
		if Hnum in self.hypotheses:
			return(self.hypothesis_info.loc[Hnum][desired_var])
		else:
			raise ValueError("invalid hypothesis (Hnum): " + Hnum + \
				'; Hnum may only be one of the following: ' + \
				', '.join(self.hypotheses))
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
	def __init__(self, output_id_sim, data_folder, sim_key, sim_number_list, \
		sim_parameters, hypothesis_testing_info):
		self.output_id_sim = output_id_sim
			# output_id_sim needs to include sim key but not sim number; fixed parameter will be added to filenames by mleparameters, and sim number will be added in MLE
		self.sim_number_list = sim_number_list
		self.sim_key = sim_key
		self.sim_parameters = copy.deepcopy(sim_parameters)
		self.hypothesis_testing_info = hypothesis_testing_info
		self.hypotheses = hypothesis_testing_info.get_hypotheses()
	def _run_MLE(self, )
	def _find_LLs(self, Hnum, current_fixed_param):
		# use mleparameters and set mode and fixed_param; change starting vals to whatever was used in simulation
		current_output_id = output_id_sim + '_' + Hnum
		self.mle_parameters.set_mode(mode, current_output_id)
		self.mle_parameters.set_parameter(current_fixed_param)
	def _calculate_LLR(self):

######## ??? #######
# Before running LLRCalculator on original data, a copy of it needs to
# be created in the same place where simulation folders are stored;
# The 'unfixed' output file for this data should be placed in this folder
######## ??? #######
class FixedPointPvalCalculator(object):
	'''
	To calculate p-val at a given point:
		1. 	a. Calculate H1: 'unfixed' LL (all unknown parameters
			freely estimated)
			b. Calculate H0: LL at fixed point in parameter space
			c. 'True' log likelihood ratio 'LLR' is LL(H0)-LL(H1)
				(actually to make the p-vals easier to understand, it's
				easier to work with the negative of LLR values)
			d. MLE parameter values for all parameters estimated in (b)
			will be used as sim starting values!
		2. Run sim_number simulations using starting parameter values
		from (1d)
		3. Repeat and record (1) but with sim data rather than real data
	'''
	def __init__(self, fixed_param, mode):



			
#####################################################################
