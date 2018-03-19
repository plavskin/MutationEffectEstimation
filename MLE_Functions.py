#!/usr/bin/python

# Contains objects needed for running Maximum Likelihood Estimation

import os
import numpy
import Cluster_Functions

class FolderManager(object):
	def __init__(self,cluster_parameters,experiment_folder_name):
		self.experiment_folder_name = experiment_folder_name
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
	############### ??? TO DO ??? ###############
	# Check that each mode or parameter is in the list once?
		# (currently only pays attention to the first time mode or parameter listed)
	############### ??? TO DO ??? ###############
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
		self.ms_positions = parameter_list["multistart_positions"].split('|')
		self.ms_grid_dimensions = parameter_list["multistart_grid_dimensions"].split('|')
		self.parallel_processors = int(parameter_list["parallel_processors"])
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
		self.current_parameters_to_loop_over = ['unfixed'] + non_permafixed_parameters
		# include 1 profile point for 'unfixed' setting
			# (i.e. only perform 'unfixed' MLE once, since no need to
			# get likelihood profile)
		self.point_numbers_to_loop_over = numpy.append([1], \
			self.current_profile_point_list[ \
			numpy.invert(pernamently_fixed_parameter_bool)])
	def set_mode(self,mode_name,output_identifier):
		# for all MLE_parameter attributes, retrieves the parameter or
			# list of parameters corresponding to the current mode
		self.current_mode = mode_name
		mode_idx = self.mode_list.index(mode_name)
		self.current_sim_repeats = self.sim_repeats_by_mode[mode_idx]
		self.current_ms_positions = int(self.ms_positions[mode_idx])
		self.mode_ms_grid_dimensions = int(self.ms_grid_dimensions[mode_idx])
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
		# output_identifier is a string that will be included in filenames
		self.output_identifier = output_identifier
	def set_parameter(self,parameter_name):
		# set current parameter, number of likelihood profile points
			# for it, and create a temporary list of fixed parameters
			# that includes it
		self.current_fixed_parameter = parameter_name
		if current_fixed_parameter == 'unfixed':
			self.current_profile_point_num = 1
			# list of fixed parameters is unchanged from default
			self.current_tempfixed_parameter_bool = \
				self.current_permafixed_parameter_bool
			# other properties need to be NaN
			self.current_ms_grid_dimensions = self.mode_ms_grid_dimensions
		else:
			# find index of current_fixed_parameter in parameter list
			current_fixed_parameter_idx = self.current_parameter_list.index(parameter_name)
			# temporarily fix current parameter
			self.current_tempfixed_parameter_bool = \
				self.current_permafixed_parameter_bool
			self.current_tempfixed_parameter_bool[current_fixed_parameter_idx] = \
				True
			self.current_profile_point_num = \
				self.current_profile_points_list[current_fixed_parameter_idx]
			self.current_ms_grid_dimensions = self.mode_ms_grid_dimensions-1
		self.output_id_parameter = self.output_identifier + '_' self.current_fixed_parameter
		self.current_parallel_processors = min(self.parallel_processors,
			self.current_ms_positions**self.current_ms_grid_dimensions)

class MLEstimation(object):
	def __init__(self,mle_parameters,cluster_folders,mle_folders,input_data_list):
		self.mle_parameters = mle_parameters
		self.cluster_folders = cluster_folders
		self.mle_folders = mle_folders
		self.completefile = os.path.join(cluster_folders.completefile_path, \
			('MLE_' + mle_parameters.output_identifier + '_completefile.txt'))
		self.job_name = mle_folders.experiment_folder_name + '-MLE-' + \
			mle_parameters.output_id_parameter
		self.input_data_list = input_data_list
		self.output_filename = 'data_'+mle_parameters.output_id_parameter
		self.output_extension = 'csv'
		self.output_path = os.path.join(mle_folders.MLE_output_path,'csv_output')


########################################################################

def ML_Estimator(completefile_folder,trackfile_folder,profile_points,
	current_mode,mode_list,username,MLE_time,MLE_memory,slurm_folder,
	current_fixed_parameter,current_fixed_parameter_bool_list,
	MLE_output_folder,growth_condition,current_L,current_gridpower,
	phenotype_file,petite_file,rep,strain_number,
	parallel_processors_max,output_dir_name,max_memory,max_time,
	start_values,parameter_max_list,parameter_min_list,
	pernamently_fixed_parameter_bool,current_parameter_list,
	total_parameter_number,ms_positions,parameter_profile_ub_list,
	parameter_profile_lb_list):
	# runs MLE


	if not os.path.isfile(current_MLE_completefile):
		# Update the trackfile; get list of jobs, if any, for which
			# simulation still needs to be run

		MLE_memory = MLE_memory*parallel_processors

		[new_jobs_to_submit,aborted_jobs_to_submit,aborted_newtimes_to_submit, \
			aborted_newmems_to_submit] = Trackfile_Processor(profile_points,
				current_mode,current_MLE_trackfile,current_MLE_completefile,
				username,current_job_name,slurm_folder,
				folder_with_csv_output_files,current_output_filename,
				current_output_extension,output_dir_name,
				max_memory,max_time,MLE_memory,MLE_time)

		if new_jobs_to_submit or aborted_jobs_to_submit:
			# write a sbatch file to submit a batch process
			current_sbatch_filename = slurm_folder + '/' + current_job_name + '.q'
			# when doing MLE in 'knownmuts' mode, LL for each strain
				# has to be calculated separately, so MLE code uses
				# parfor loop, making code run faster on multiple
				# processors

			if parallel_processors > 1:
				additional_beginning_lines_in_sbatch = ['if [ \"$SLURM_JOBTMP" == \"\" ]; then',\
				'    export SLURM_JOBTMP=/state/partition1/$USER/$$',\
				'    mkdir -p $SLURM_JOBTMP',\
				'fi',\
				'export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)']#,\
#					'export NTHREADS=$(cat $PBS_NODEFILE | wc -l)']
				additional_end_lines_in_sbatch = ['rm -rf $SLURM_JOBTMP/*']
			else:
				additional_beginning_lines_in_sbatch = []
				additional_end_lines_in_sbatch = []

			module_to_use = 'matlab'

			####### ??? TO DO ??? #######
			# need to consider/try parallelization of the strainwise stage
			# Fix parallelization code above to be up to date with SLURM
			####### ??? TO DO ??? #######

			parameter_list = '{\'' + '\',\''.join(current_parameter_list) + '\'}'
			combined_fixed_parameter_array = '[' + ','.join(map(str,(current_fixed_parameter_bool_list*1))) + ']'
				# current_fixed_parameter_bool_list should be a numpy
					# array of boolean values, so multiplying it by 1
					# produces integers of 1 or 0, which is input matlab
					# can take
			combined_min_array = '[' + ','.join(map(str,parameter_min_list)) + ']'
			combined_max_array = '[' + ','.join(map(str,parameter_max_list)) + ']'
			combined_profile_ub_array = '[' + ','.join(map(str,parameter_profile_ub_list)) + ']'
			combined_profile_lb_array = '[' + ','.join(map(str,parameter_profile_lb_list)) + ']'
			combined_length_array = '[' + ','.join([str(profile_points)]*total_parameter_number) + ']'
			combined_position_array = '[$SLURM_ARRAY_TASK_ID]'
				# if combined_position_array has length=1, MLE programs
					# interpret it as an array of the correct length
					# with the same value repeated
			combined_start_values_array = '[' + ','.join(map(str,start_values)) + ']'
			external_counter = '$SLURM_ARRAY_TASK_ID'
			assigned_parallel_processors = '$SLURM_CPUS_PER_TASK'

			csv_output_prename = growth_condition + '_data_' + current_mode + \
				'_' + current_fixed_parameter

			matlab_input_list = [combined_fixed_parameter_array,external_counter, \
				combined_min_array,combined_max_array,combined_length_array, \
				combined_position_array,combined_start_values_array, \
				parameter_list,('\''+csv_output_prename+'\''),('\''+MLE_output_folder+'\''), \
				('\''+mat_phenotype_file+'\''),assigned_parallel_processors,str(ms_positions),
				('\''+petite_file+'\''),combined_profile_ub_array,combined_profile_lb_array]

			code_run_string = ('\'MLE_' + current_mode + '(\'\"'
				+ '\"\",\"\"'.join(matlab_input_list) + '\"\");exit\"')

			if current_mode =='mixed':
				default_single_job_time = 3*MLE_time
				# 'mixed' mode requires additional time
			elif current_mode == 'mixedgauss':
				default_single_job_time = 2*MLE_time
			else:
				default_single_job_time = MLE_time
			default_single_job_mem = MLE_memory # in MB

			# submit aborted jobs one-by-one, new_jobs_to_submit as a batch job
			Job_List_Submitter(new_jobs_to_submit,aborted_jobs_to_submit,aborted_newtimes_to_submit,
				aborted_newmems_to_submit,current_sbatch_filename,
				parallel_processors,code_dir,user_email,module_to_use,
				code_run_string,current_job_name,
				additional_beginning_lines_in_sbatch,additional_end_lines_in_sbatch,
				default_single_job_time,default_single_job_mem)

	return current_MLE_completefile







