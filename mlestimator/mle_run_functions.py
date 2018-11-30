#!/usr/bin/python

''' Contains function needed to run MLE '''

import os
import copy
from mlestimator.mle_functions import MLEParameters, run_MLE
from mlestimator.mle_result_combiner import CombinedResultSummary
from mlestimator.mle_sim_functions import generate_sim_based_profile_pts
from cluster_wrangler.cluster_functions import CompletenessTracker

def loop_over_modes(mle_parameters, cluster_parameters, cluster_folders, \
	mle_folders, experiment_path, additional_code_run_keys, \
	additional_code_run_values, output_id_string_start, sim_parameters):
	# Handles all of MLE across modes, including confidence
		# interval identification
	mode_completeness_tracker = CompletenessTracker(mle_parameters.mode_list)
	for current_mode in mle_parameters.mode_list:
		##### RUN MLE #####
		output_id_string = '_'.join([output_id_string_start, current_mode])
		mle_parameters.set_mode(current_mode, output_id_string)
		output_id_string_sim = 'sim_' + output_id_string
		sim_parameters.set_mode(current_mode, output_id_string_sim)
#		MLE_summary_file_path = os.path.join(experiment_path, \
#			'_'.join([output_id_string,'MLE_file.csv']))
		# run MLE for current set of parameters
		include_unfixed_param = True
		input_data_folder = mle_folders.get_path('experiment_path')
		run_MLE(mle_parameters, cluster_parameters, cluster_folders, \
			mle_folders, additional_code_run_keys, additional_code_run_values, \
			include_unfixed_param, input_data_folder)
		# generate LL_profiles and MLE_output file, and identify
			# asymptotic CIs
		current_combined_results = CombinedResultSummary(mle_folders, \
			mle_parameters, cluster_parameters, cluster_folders, \
			sim_parameters)
		current_combined_results.generate_CIs('asymptotic')
		# run simulations, run MLE on those simulations, and generate
			# empirical sim-based cdf value profile points across each
			# parameter for which sim-based CIs are required
		sim_folders = copy.deepcopy(mle_folders)
		sim_folders.set_current_output_subfolder('sim')
		sim_MLEs_completefile = generate_sim_based_profile_pts(current_mode, sim_parameters, \
			sim_folders, additional_code_run_keys, additional_code_run_values, \
			output_id_string_sim, current_combined_results, cluster_folders, \
			cluster_parameters)
		current_combined_results.initialize_sim_based_results(sim_MLEs_completefile)
		current_combined_results.generate_CIs('sim_based')
		# update completeness of current mode based on whether or not combined results have been completed
		current_mode_completeness = \
			current_combined_results.get_completeness()
		mode_completeness_tracker.switch_key_completeness(current_mode, \
			current_mode_completeness)
	mode_loop_completeness = mode_completeness_tracker.get_completeness()
	return(mode_loop_completeness)


