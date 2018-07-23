#!/usr/bin/python

# Contains functions necessary for identifying MLE Confidence Intervals

import os
import pandas as pd
import numpy as np
from mle_filenaming_functions import generate_file_label
import Cluster_Functions


class CIWarning(object):
	# stores warnings for OneSidedCIBound objects
	def __init__(self):
		self.points_to_create_CI = 0
		self.all_points_within_CI_bound = False
		self.non_monotonic = False
	def set_all_points_within_CI_bound(self):
		self.all_points_within_CI_bound = True
	def set_points_to_create_CI(self,points_to_create_CI):
		self.points_to_create_CI = points_to_create_CI
	def set_no_output_file(self):
		self.no_output_file = True
	def set_non_monotonic(self):
		self.non_monotonic = True
	def get_warning_line(self, CI_side, CI_type):
		warning_list = []
		if self.all_points_within_CI_bound:
			warning_list.append('All points used to calculate ' + CI_side + \
				 ' CI for ' + CI_type + \
				 ' CI are between expected CI bound and MLE.')
		if self.no_output_file:
			warning_list.append('Although the curve-fitting job was completed, there was no output file created for ' + \
				CI_side + 'CI for ' + CI_type)
		if self.non_monotonic:
			warning_list.append('non-monotonic LL profile on ' + CI_side + ' side of ML parameter value')
		if set_points_to_create_CI == 0:
			warning_list.append('It seems '  + CI_side + ' CI for ' + \
				CI_type + ' CI was not set.')
		elif set_points_to_create_CI == 1:
			warning_list.append('No profile points besides MLE on ' + \
				CI_side + ' side of CI for ' + CI_type + ' CI.')
		elif set_points_to_create_CI == 2:
			warning_list.append('Only one profile point besides MLE on ' + \
				CI_side + \
				' side of CI for ' + CI_type + \
				'CI; CI boundary calculated as linear interpolation of p-val between to two points')
		warning_string = ';'.join(warning_list)
		return(warning_string)

class OneSidedCIBound(object):
	# Stores data for CI bound on one side of MLE
	def __init__(self, p_val, LL_df, df, fixed_param_MLE_val,fixed_param, \
		max_LL, CI_type, mle_folders, cluster_parameters, cluster_folders):
		self.cdf_bound = 1-p_val
		self.points_to_fit_curve = 3
		self.df = df
		self.CI_type = CI_type
			# CI_type can be either 'asymptotic' or 'sim'
		self.default_CI_bound = None
		self.fixed_param = fixed_param
		self.fixed_param_MLE_val = fixed_param_MLE_val
		self.max_LL = max_LL
		self.cluster_parameters = cluster_parameters
		self.cluster_folders = cluster_folders
		self.mle_folders = mle_folders
		self.module = 'matlab'
		self.CI_bound_set = False
		self.additional_beginning_lines_in_job_sub = []
		self.additional_end_lines_in_job_sub = []
		self.additional_code_run_keys = []
		self.additional_code_run_values = []
		self._select_profile_side(LL_df)
		self.warning = CIWarning()
		self.one_sided_LL_df = pd.DataFrame()
		# need to run _asymptotic_p_val_calc and
			# _find_CI_proximal_LL_points even if curve fitting to LL
			# points has already occurred, since these functions throw
			# important warnings
		self._asymptotic_p_val_calc()
		self._find_CI_proximal_LL_points()
	def _set_output_filenames(self, profile_side):
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
	def _select_profile_side(self, LL_df):
		# selects data from correct side of likelihood profile
		print('error! Profile side not selected in parent class of OneSidedCIBound objects')
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
	def _check_monotonicity(self, x_vals, y_vals):
		# check whether LL profile is monotonic before and after ML parameter value
		# In simple and accurately estimated LL landscape, LL
			# profile expected to increase monotonically up until
			# max LL, then decrease monotonically after; if this
			# isn't the case, throw a warning
		print('error! Monotonicity not checked in parent class of OneSidedCIBound objects')
	def _id_proximal_points(self, target_y, y_vals, num_points):
		# returns indices of up to num_points points closest to target_y
		dist_to_target = np.absolute(target_y-y_vals)
		sorted_indices = np.argsort(dist_to_target)
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
		if np.max(profile_points) < self.cdf_bound:
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

class OneSidedCIBoundLower(OneSidedCIBound):
	def __init__(self, p_val, LL_df, df, fixed_param_MLE_val, fixed_param, \
		max_LL, CI_type, mle_folders, cluster_parameters, cluster_folders):
		super(OneSidedCIBound).__init__(p_val, LL_df, df, fixed_param_MLE_val, \
			fixed_param, max_LL, CI_type, mle_folders, cluster_parameters, \
			cluster_folders)
		self._set_output_filenames('upper')
#		self._set_output_filenames('lower')
	def _select_profile_side(self, LL_df):
		# selects data from correct side of likelihood profile
		self.default_CI_bound = float("-inf")
		LL_df_one_side = LL_df[(LL_df[self.fixed_param] <= self.fixed_param_MLE_val)]		
#		LL_df_one_side = LL_df[(LL_df[self.fixed_param] >= self.fixed_param_MLE_val)]
#		self.default_CI_bound = float("inf")
		self.one_sided_LL_df = LL_df_one_side.sort_values(self.fixed_param)
		self._check_monotonicity()
	def _check_monotonicity(self):
		# check whether LL profile is monotonic before or after ML parameter value
		# In simple and accurately estimated LL landscape, LL
			# profile expected to increase monotonically up until
			# max LL, then decrease monotonically after; if this
			# isn't the case, throw a warning
		y_diffs = np.diff(self.one_sided_LL_df['LL'])
		monotonicity_state = np.all(y_diffs >= 0)
#		monotonicity_state = np.all(y_diffs <= 0)
		if not monotonicity_state:
			self.warnings.set_non_monotonic()

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
		self.CI_object_dictionary = {}
		self.CI_object_dictionary['lower'] = OneSidedCIBoundLower(p_val, \
				LL_df, deg_freedom, fixed_param_MLE_val, fixed_param, max_LL, \
				CI_type, mle_folders, cluster_parameters, cluster_folders)
		self.CI_object_dictionary['upper'] = OneSidedCIBoundUpper(p_val, \
				LL_df, deg_freedom, fixed_param_MLE_val, fixed_param, max_LL, \
				CI_type, mle_folders, cluster_parameters, cluster_folders)
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
