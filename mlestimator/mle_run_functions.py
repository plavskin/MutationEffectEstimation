#!/usr/bin/python

''' Contains function needed to run MLE '''

import os
from mlestimator.mle_functions import MLEParameters, run_MLE
from mlestimator.mle_result_combiner import CombinedResultSummary
from mlestimator.mle_sim_functions import generate_sim_based_profile_pts

def loop_over_modes(mle_parameters, cluster_parameters, cluster_folders, \
	mle_folders, experiment_path, additional_code_run_keys, \
	additional_code_run_values, output_id_string_start, sim_parameters, \
	sim_folders):
	# Handles all of MLE across modes, including confidence
		# interval identification
	for current_mode in mle_parameters.mode_list:
		##### RUN MLE #####
		output_id_string = '_'.join([output_id_string_start, current_mode])
		mle_parameters.set_mode(current_mode, output_id_string)
#		MLE_summary_file_path = os.path.join(experiment_path, \
#			'_'.join([output_id_string,'MLE_file.csv']))
		# run MLE for current set of parameters
		include_unfixed_param = True
		input_data_folder = mle_folders.get_path('experiment_path')
		run_MLE(mle_parameters, cluster_parameters, cluster_folders, \
			mle_folders, additional_code_run_keys, additional_code_run_values, \
			include_unfixed_param)
		# generate LL_profiles and MLE_output file, and identify
			# asymptotic CIs
		current_combined_results = CombinedResultSummary(mle_folders, \
			mle_parameters, cluster_parameters, cluster_folders, \
			sim_parameters, sim_folders)
		current_combined_results.generate_CIs('asymptotic')
		# run simulations, run MLE on those simulations, and generate
			# empirical sim-based cdf value profile points across each
			# parameter for which sim-based CIs are required
		output_id_string_sim = 'sim' + output_id_string
		sim_parameters.set_mode(current_mode, output_id_string_sim)
		generate_sim_based_profile_pts(current_mode, sim_parameters, \
			sim_folders, additional_code_run_keys, additional_code_run_values, \
			output_id_string_sim, current_combined_results, cluster_folders, \
			cluster_parameters)
		current_combined_results.generate_CIs('sim_based')


