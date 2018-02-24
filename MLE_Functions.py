#!/usr/bin/python

# Contains objects needed for running Maximum Likelihood Estimation

import os

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