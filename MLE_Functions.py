#!/usr/bin/python

# Contains objects needed for running Maximum Likelihood Estimation

import os
import numpy

class FolderManager(object):
	def __init__(self,cluster_parameters,experiment_folder_name):
		self.experiment_path = \
			os.path.join(cluster_parameters.composite_data_folder,experiment_folder_name)
		self.sim_output_path = os.path.join(cluster_parameters.temp_storage_path, \
			experiment_folder_name,'simulated_phenotypes')
		self.MLE_output_path = os.path.join(cluster_parameters.temp_storage_path, \
			experiment_folder_name,'MLE_output')
		self.LL_profile_path = os.path.join(self.experiment_path,'/LL_profiles')
	#	self.epilogue_path = os.path.join(cluster_parameters.home_path,cluster_parameters.username,'mut_effect_epilogue_files',experiment_folder_name)
		self.MLE_combined_sim_path = os.path.join(self.experiment_path,'/MLE_sim_outputs')
		self._set_up_folders()
	def _set_up_folders(self):
		setup_complete_file = os.path.join(self.completefile_folder,'folder_setup_complete.txt')
		if not os.path.isfile(setup_complete_file):
			new_directory_list = (self.sim_output_path,self.MLE_output_path, \
				self.LL_profile_path,self.MLE_combined_sim_path)
			for current_new_directory in new_directory_list:
				if not os.path.isdir(current_new_directory):
					os.makedirs(current_new_directory)
			open(setup_complete_file,'a').close()
	def set_current_output_subfolder(self,current_subfolder):
		# set (and if necessary, create) a subfolder to write temp output to
		self.current_output_subfolder = os.path.join(self.MLE_output_path, \
			current_subfolder)
		if not os.path.isdir(self.current_output_subfolder):
			os.makedirs(self.current_output_subfolder)


class MLEParameters(object):
	def __init__(self,parameter_list):
		self.sim_repeats_by_mode = parameter_list["sim_repeats"].split('|')
		self.profile_points_by_mode = parameter_list["profile_points_by_mode"].split('|')
		self.mode_list = parameter_list["mode_list"].split('|')
		self.parameters_by_mode = parameter_list["parameters_by_mode"].split('|')
		self.min_parameter_vals_by_mode = \
			parameter_list["min_parameter_vals_by_mode"].split('|')
		self.max_parameter_vals_by_mode = \
			parameter_list["max_parameter_vals_by_mode"].split('|')
		self.starting_parameter_vals_by_mode = \
			parameter_list["starting_parameter_vals_by_mode"].split('|')
		self.profile_lower_limits_by_mode = \
			parameter_list["profile_lower_limits_by_mode"].split('|')
		self.profile_upper_limits_by_mode = \
			parameter_list["profile_upper_limits_by_mode"].split('|')
	def _retrieve_current_values(list_by_mode,mode_idx,current_mode,param_num):
		# retrieves the appropriate list from a list of parameter lists
			# by mode
		if len(list_by_mode) > 1:
			output_list = list_by_mode[mode_idx].split(';')
		else:
			output_list = list_by_mode.split(';')
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
		# include 'unfixed' 'parameter', in which case no parameter is fixed
		parameters_to_loop_over_bool = \
			numpy.invert(self.current_permafixed_parameter_bool)* \
			numpy.invert((self.current_profile_point_list < 1))
		non_permafixed_parameters = [item for (item,bool_val) in \
			zip(self.current_parameter_list,parameters_to_loop_over_bool) \
			if bool_val]
		self.parameters_to_loop_over = ['unfixed'] + non_permafixed_parameters
		# include 1 profile point for 'unfixed' setting
			# (i.e. only perform 'unfixed' MLE once, since no need to
			# get likelihood profile)
		self.point_numbers_to_loop_over = numpy.append([1], \
			self.current_profile_point_list[ \
			numpy.invert(pernamently_fixed_parameter_bool)])
	def set_mode(self,mode_idx):
		# for all MLE_parameter attributes, retrieves the parameter or
			# list of parameters corresponding to the current mode
		self.current_mode = self._retrieve_current_values(self.mode_list,mode_idx)
		self.current_sim_repeats = self.sim_repeats_by_mode[mode_idx]
		self.current_parameter_list = self.parameters_by_mode[mode_idx].split(';')
		# find the total number of parameters, including fixed ones
		self.total_param_num = len(self.current_parameter_list)
		# create lists, of length total_param_num, of settings for each
			# parameter in current_parameter_list
		# if input in file was incorrect length relative to number of
			# parameters, corresponding list is just NaN
		self.current_min_parameter_val_list = \
			numpy.array(self._retrieve_current_values(min_parameter_vals_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_max_parameter_val_list = \
			numpy.array(self._retrieve_current_values(max_parameter_vals_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_start_parameter_val_list = \
			numpy.array(self._retrieve_current_values(starting_parameter_vals_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_profile_point_list = \
			numpy.array(self._retrieve_current_values(profile_points_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_profile_lower_limit_list = \
			numpy.array(self._retrieve_current_values(profile_lower_limits_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		self.current_profile_upper_limit_list = \
			numpy.array(self._retrieve_current_values(profile_upper_limits_by_mode,\
				mode_idx,self.current_mode,self.total_param_num))
		# identify list of parameters that are permanently fixed
		self.current_permafixed_parameter_bool = \
			self.current_max_parameter_val_list == self.current_min_parameter_val_list
		# identify parameters MLE needs to be performed on
		self._id_parameters_to_loop_over()
		# set up a list of completefiles and the current summary file
		self.current_completefile_list = []






