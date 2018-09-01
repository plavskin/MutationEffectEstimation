#!/usr/bin/python

# Contains functions necessary for running simulations and combining their results

import pandas as pd
import numpy as np

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

class SimParams(object):
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
	def _find_matching_key(self, sim_params):
		'''
		Takes in a SimParams object and identifies the row in
		self.sim_key_df that matches it, if any do
		'''
		current_keys = self.sim_key_df.index.values.tolist()
		row_match_boolean_list = [sim_params.equals(self.sim_key_df.loc[i]) for i in current_keys]
		matching_keys = self.sim_key_df.index[row_match_boolean_list].tolist()
		if matching_keys:
			key_val = matching_keys[0]
		else:
			key_val = None
		return(key_val)
	def get_sim_key(self, sim_params):
		'''
		Takes in a SimParams object and determines whether it's already
		in self.sim_key_df
		If not, adds it, and writes updated file to self.sim_key_file
		Either way, returns the sim_key for this set of SimParams
		(i.e. the index of the line in sim_key_df holding these params)
		'''
		current_key = self._find_matching_key(sim_params)
		if not current_key:
			self.sim_key_df = \
				self.sim_key_df.append(sim_params, sort = False, \
					ignore_index = True)
			current_key = self._find_matching_key(sim_params)
			self.write_key_df()
		return(current_key)
	def write_key_df(self):
		''' Writes self.sim_key_df to self.sim_key_file '''
		self.sim_key_df.to_csv(path_or_buf = self.sim_key_file, index = True)


			
