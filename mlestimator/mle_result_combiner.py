#!/usr/bin/python

'''
Contains function needed to initiate CI compilation and combine MLE
results into summary table
'''

import pandas as pd
import numpy as np
import warnings as war
import mle_functions
import mle_CI_functions
import mle_sim_functions
from cluster_wrangler import cluster_functions

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
			self.content_dict = {current_key:np.NaN for current_key in key_list}
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
				not np.any([np.isnan(self.content_dict[key]) for \
					key in searchstring_keys])
		return(columns_filled)
	def get_contents_by_searchstring(self, searchstring):
		''' Get dictionary elements whose keys contain searchstring '''
		searchstring_keys = self._get_keys_by_searchstring(searchstring)
		searchstring_dict = {key: self.content_dict[key] for key in searchstring_keys}
		return(searchsting_dict)
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
		self.LL_profile_folder = mle_folders.get_path('LL_list_path')
		self.runtime_percentile = mle_parameters.runtime_percentile
		self.pval = mle_parameters.current_CI_pval
		self._create_combined_output_file()
		self.max_LL = None
		self.warnings = CombinedResultWarning()
		self.completeness_tracker = cluster_functions.CompletenessTracker(['initialization', \
			'asymptotic_CIs', 'simulation-based_CIs'])
		self.runtime_quant_list = np.array([])
		self.true_max_param_df = pd.DataFrame()
		self.combined_results_df = pd.DataFrame()
		# set completefiles for asymptotic and sim-based CIs
		self.asymptotic_CI_completefile = \
			os.path.join(cluster_folders.get_path('completefile_path'), \
				'_'.join(['asymoptotic_CI', mle_parameters.output_identifier, \
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
			self.unfixed_ll_param_df = pd.DataFrame()
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
		ll_profile = mle_functions.LLProfile(self.mle_parameters, \
			self.mle_datafile_path, self.LL_profile_folder, \
			non_profile_max_params)
		ll_profile.run_LL_list_compilation()
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
			self.runtime_quant_list = np.append(self.runtime_quant_list,runtime_quantile)
	def _correct_LL_profiles(self):
		# updates LL_profiles
		parameters_to_loop_over = self.mle_parameters.get_fitted_parameter_list(False)
		ml_params_to_include = pd.concat([self.unfixed_ll_param_df, \
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
				np.percentile(self.runtime_quant_list,self.pval/2*100)
			runtime_CI_bounds['upper'] = \
				np.percentile(self.runtime_quant_list,(1-self.pval/2)*100)
		else:
			runtime_CI_bounds['lower'] = np.NaN
			runtime_CI_bounds['upper'] = np.NaN
		return(runtime_CI_bounds)
	def _add_line_to_combined_df(self,fixed_param_name, current_results_dict):
		# creates self.combined_summary or adds a line to it
		if self.combined_results_df.empty:
			# create dataframe using currently passed line
			self.combined_results_df = \
				pd.DataFrame(current_results_dict, index = [fixed_param_name])
		else:
			# check that all keys in current_results_dict are columns
				# in self.combined_results_df; if not, create these
				# columns before updating data
			dict_keys = current_results_dict.keys()
			df_columns = list(self.combined_results_df.columns.values)
			new_columns = set(dict_keys).symmetric_difference(df_columns)
			self.combined_results_df = \
				self.combined_results_df.append(pd.DataFrame(columns=new_columns), \
					sort=False)
			# insert current_results_dict into correct row in df
			current_results_series = pd.Series(current_results_dict)
			self.combined_results_df.loc[fixed_param_name] = \
				current_results_series
	def _write_combined_df(self):
		# writes self.combined_results_df to file
		self.combined_results_df.to_csv(path_or_buf=self.combined_results_output_file, \
			index=True)
	def _read_combined_df(self):
		# reads self.combined_results_df from a pre-recorded file
		self.combined_results_df = \
			pd.read_csv(filepath_or_buffer=self.combined_results_output_file, \
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
		mean_runtime_quant = np.mean(self.runtime_quant_list)
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
							pd.DataFrame())
						ll_profile.run_LL_list_compilation()
						ll_profile.run_CI(self.pval, self.mle_folders, \
							self.cluster_parameters, self.cluster_folders, \
							self.mle_parameters)
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
	def get_CI_dict(self, fixed_param, CI_type):
		'''
		Get a dictionary with 'lower' and 'upper' as keys for CI_type
		bounds
		'''
		current_results_line = SingleParamResultSummary()
		current_results_line.read_from_df_line(self.combined_results_df, \
			fixed_param)
		CI_bound_dict = \
			current_results_line.get_contents_by_searchstring(CI_type)
		CI_bound_names = ['lower','upper']
		CI_bound_dict_renamed = {}
		for current_bound_name in CI_bound_names:
			bound_name_dict = \
				current_results_line.get_contents_by_searchstring(current_bound_name)
			current_bound_key = list(set(CI_bound_dict.keys()) & set(bound_name_dict.keys()))
			current_bound_val_list = [val for key, val in CI_bound_dict if key in current_bound_key]
			if len(current_bound_val_list) > 1:
				raise ValueError('Too many values for ' + \
					current_bound_name + 'side of ' + CI_type + 'bound')
			else:
				current_bound_val = current_bound_val_list[0]
			CI_bound_dict_renamed[current_bound_name] = current_bound_val
		return(CI_bound_dict_renamed)
	def get_param_mle_val(self, fixed_param):
		return(self.content_dict.loc[fixed_param]['param_MLE'])


			# keep track of asymptotic and sim CIs using completeness tracker across parameters within mode
		# once CIs complete (either asymptotic only or asymptotic and sim, depending on settings),
		#	create summary file
		# keep track of summary files for each mode; once those are complete, stop running current folder
		##

