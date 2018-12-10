#!/usr/bin/python

"""
Contains objects needed for running jobs on cluster and keeping track
of their progress
"""

import os
import csv
import subprocess
import numpy as np
import copy
import re
import cluster_sub_functions
import warnings as war
war.filterwarnings("ignore", message="numpy.dtype size changed")

class Job(object):
	""" job object that stores properties of individual jobs in queue """
	def __init__(self, number, status, time, mem):
		self.number = number
		self.status = status
		self.time = time
		self.mem = mem
	def change_status(self,new_status):
		self.status = new_status
	def get_status(self):
		return(self.status)
	def get_time(self):
		return(self.time)
	def get_mem(self):
		return(self.mem)
	def change_mem(self,new_mem):
		self.mem = new_mem
	def change_time(self,new_time):
		self.time = new_time
	def extract_job_info(self):
		info_list = [self.number,self.status,self.time,self.mem]
		return(info_list)

class JobStatus(object):
	""" enum for job status """
	TO_PROCESS = 'To Process'
	PROCESSING = 'Processing'
	COMPLETED = 'Completed'
	ABORTED_FOREVER = 'Aborted Forever'
	ABORTED_TO_RESTART = 'Aborted to Restart'
	ERROR = 'Error'

class Parameter(object):
	""" processes a single parameter passed to it """
	def __init__(self, name, value, current_type, explanation):
		self.name = name
		self.explanation = explanation
		self.type = current_type
		self.value = value
		self._parse_value()
	def concat_to_string(self):
		""" returns a string with parameter info """
		concatenated_string = ""
		concatenated_string += self.name
		if self.type == "str":
			value_str = self.value
		elif self.type == "str_list":
			value_str = ','.join(self.value)
		elif self.type == "nested_str_list":
			value_str = '[['+'],['.join(','.join(current_sublist) for \
				current_sublist in self.value) + ']]'
		elif self.type == "float" or self.type == 'int' or self.type == 'bool':
			value_str = str(self.value)
		elif self.type == "float_list" or self.type == "int_list" or self.type == "bool_list":
			value_str = ','.join(str(element) for element in self.value)
		elif self.type == "nested_float_list" or self.type == "nested_int_list" or self.type == "nested_bool_list":
			value_str = '[['+'],['.join(','.join(str(element) for element in \
				current_sublist) for current_sublist in self.value) + ']]'
		concatenated_string += '\n'+value_str
		concatenated_string += '\n'+self.type
		concatenated_string += '\n'+self.explanation
		return(concatenated_string)
	def _parse_value(self):
		""" determines what type the parameter belongs to """
		if   self.type == "str":
			self.value = self.value
		elif self.type == "nested_str_list":
			self.value = [i.split(";") for i in self.value.split("|")]
		elif self.type == "str_list":
			self.value = self.value.split("|")
		elif self.type == "float":
			self.value = float(self.value)
		elif self.type == "nested_float_list":
			self.value = [[float(j) for j in i.split(";")] for i in self.value.split("|")]
		elif self.type == "float_list":
			self.value = [float(i) for i in self.value.split("|")]
		elif self.type == "int":
			self.value = int(self.value)
		elif self.type == "nested_int_list":
			self.value = [[int(j) for j in i.split(";")] for i in self.value.split("|")]
		elif self.type == "int_list":
			self.value = [int(i) for i in self.value.split("|")]
		elif self.type == "bool":
			self.value = bool(self.value)
		elif self.type == "nested_bool_list":
			self.value = [[bool(j) for j in i.split(";")] for i in self.value.split("|")]
		elif self.type == "bool_list":
			self.value = [bool(i) for i in self.value.split("|")]
		else:
			raise ValueError("invalid data type: " + self.type)

class Parameters(object):
	""" creates a list of arbitrary parameters """
	def __init__(self, parameters):
		self.parameters = {}
		for p in parameters:
			self.parameters[p.name] = p
	def __getitem__(self, name):
		if name in self.parameters:
			return self.parameters[name].value
		else:
			raise ValueError('this parameter does not exist: ' + name)
	def print_values(self):
		for name in self.parameters:
			print(self.parameters[name].concat_to_string())

class ClusterParameters(object):
	""" holds parameters general to the system """
	def __init__(self,parameter_list):
		self.code_path = parameter_list["code_folder"]
#		self.home_path = parameter_list["home_folder"]
		self.composite_data_path = parameter_list["composite_data_folder"]
		self.max_mem = parameter_list["max_mem"]
		self.max_time = parameter_list["max_time"]
		self.username = parameter_list["username"]
		self.user_email = parameter_list["user_email"]
		self.starting_mem = parameter_list["starting_mem"]
		self.starting_time = parameter_list["starting_time"]
		self.temp_storage_path = parameter_list["temp_storage_folder"]
		self.max_char_num = parameter_list["max_char_num"]
		self.max_jobs_per_batch = parameter_list["max_jobs_per_batch"]
		self.cluster_architecture = parameter_list["cluster_architecture"].lower()
		if self.cluster_architecture == 'macosx':
			self.pause_at_end = False
		else:
			self.pause_at_end = True
		self._set_cluster_architecture_properties()
		self.current_mem = self.starting_mem
		self.current_time = self.starting_time
	def _set_cluster_architecture_properties(self):
		""" sets properties related to the cluster architecture """
		if self.cluster_architecture == 'slurm':
			self.job_submission_manager = \
				cluster_sub_functions.SlurmManager(copy.deepcopy(self))
		elif self.cluster_architecture == 'macosx':
			self.job_submission_manager = \
				cluster_sub_functions.MacOSXManager(copy.deepcopy(self))
		# if cluster determines max_jobs_per_batch that is smaller than self.max_jobs_per_batch, update this value
		forced_max_jobs_per_batch = \
			self.job_submission_manager.get_max_jobs_per_batch()
		if forced_max_jobs_per_batch:
			self.max_jobs_per_batch = \
				int(round(min(forced_max_jobs_per_batch, self.max_jobs_per_batch)))
	def get_job_sub_manager(self):
		return(copy.deepcopy(self.job_submission_manager))
	def get_batch_counter_call(self):
		within_batch_counter = \
			self.job_submission_manager.get_within_batch_counter()
		within_batch_counter_call = \
			'${' + within_batch_counter + '}'
		return(within_batch_counter_call)
	def set_current_time(self, new_time):
		self.current_time = new_time
	def set_current_mem(self, new_mem):
		self.current_mem = new_mem

class FolderManager(object):
	""" Creates and holds paths used by cluster_wrangler """
	def __init__(self, cluster_parameters, experiment_folder_name):
		self.path_dict = {}
#		self.path_dict['experiment_path'] = \
#			os.path.join(cluster_parameters.composite_data_path,experiment_folder_name)
		self.path_dict['tempfolder_experiment_path'] = \
			os.path.join(cluster_parameters.temp_storage_path, experiment_folder_name)
#		self.path_dict['trackfile_path'] = \
#			os.path.join(self.path_dict['experiment_path'],'trackfiles')
#		self.path_dict['completefile_path'] = \
#			os.path.join(self.path_dict['experiment_path'],'completefiles')
		self.path_dict['trackfile_path'] = \
			os.path.join(self.path_dict['tempfolder_experiment_path'], \
				'trackfiles')
		self.path_dict['completefile_path'] = \
			os.path.join(self.path_dict['tempfolder_experiment_path'],\
				'completefiles')
		self.path_dict['cluster_job_submission_path'] = \
			os.path.join(self.path_dict['tempfolder_experiment_path'], \
				'cluster_job_submission_folder')
		self._set_up_folders()
	def _set_up_folders(self):
		setup_complete_file = os.path.join(self.path_dict['completefile_path'], \
			'cluster_folder_setup_complete.txt')
		if not os.path.isfile(setup_complete_file):
			for current_folder_key, current_path in self.path_dict.iteritems():
				if not os.path.isdir(current_path):
					os.makedirs(current_path)
			open(setup_complete_file,'a').close()
	def get_path(self, folder_name):
		return(self.path_dict[folder_name])

class JobParameters(object):
	""" Holds parameters of the job currently being run """
	def __init__(self, name, output_folder, output_extension, output_file_label, \
		cluster_job_submission_folder, experiment_folder, module, code_run_input, \
		additional_beginning_lines_in_job_sub, additional_end_lines_in_job_sub, \
		parallel_processors):
		self.name = name
		self.output_folder = output_folder
		self.output_extension = output_extension
		self.output_file_label = output_file_label
		self.cluster_job_submission_folder = cluster_job_submission_folder
		self.experiment_folder = experiment_folder
		self.module = module
		self.code_run_input = code_run_input
		self.additional_beginning_lines_in_job_sub = \
			self._string_to_list(additional_beginning_lines_in_job_sub)
		self.additional_end_lines_in_job_sub = \
			self._string_to_list(additional_end_lines_in_job_sub)
		self.parallel_processors = parallel_processors
	def _string_to_list(self, parameter):
		"""
		Checks whether parameter is entered as a string or list, and
		if a string, converts to a list holding that string; if
		neither, returns an error
		"""
		if isinstance(parameter, list):
			output_parameter = parameter
		elif isinstance(parameter, basestring):
			output_parameter = [parameter]
		else:
			print('Error: Non-string, non-list passed to JobParameters')
		return(output_parameter)

class CompletenessTracker(object):
	"""
	Keeps a running tally of keys (which can be parameters, modes, etc)
	and checks their completeness
	"""
	def __init__(self,key_list):
		if len(key_list) == 0:
			self.completeness_status = True
		else:
			self._create_completeness_dict(key_list)
			self.completeness_status = False
	def _create_completeness_dict(self,key_list):
		"""
		Creates dictionary with parameter names as keys and False
		as values
		"""
		key_list = [str(x) for x in key_list]
		completeness_list = [False]*len(key_list)
		self.completeness_dict = dict(zip(key_list,completeness_list))
	def update_key_status(self, key, completefile_path):
		"""
		Checks whether key has associated completefile and changes its
		status accordingly
		"""
		if os.path.isfile(completefile_path):
			self.completeness_dict[str(key)] = True
	def switch_key_completeness(self, key, value_bool):
		"""
		Switches a key value to value_bool
		Useful for cases when completefiles not kept track of
		"""
		self.completeness_dict[str(key)] = value_bool
	def get_key_completeness(self, key):
		return(self.completeness_dict[str(key)])
	def _check_completeness(self):
		""" Checks whethere all parameters completed """
		if all(self.completeness_dict.values()):
			self.completeness_status = True
		else:
			self.completeness_status = False
	def get_completeness(self):
		""" Checks and returns completeness status """
		if not self.completeness_status:
			self._check_completeness()
		return(self.completeness_status)

class JobListManager(object):
	"""
	Holds list of jobs corresponding to a single 'name'
	Updates current job status
	"""
	def __init__(self, jobs, job_parameters, cluster_parameters):
		self.jobs = {}
		for j in jobs:
			self.jobs[j.number] = copy.deepcopy(j)
		self.job_parameters = copy.deepcopy(job_parameters)
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
	def get_job_name(self):
		""" Returns the name of the stored jobs """
		job_name = self.job_parameters.name
		return(job_name)
#	def _get_status_list(self):
#		# gets current status of every job
#		# ??? Do I need this function? ???
#		status_list=[self.jobs[num].status for num in self.jobs]
#		return status_list
	def get_jobs_by_status(self, status):
		""" Gets list of jobs with a specific status """
		job_subset_list = []
		for num in self.jobs:
			if self.jobs[num].status == status:
				job_subset_list.append(num)
		return job_subset_list
	def get_job_times(self, number_list):
		""" Gets list of job times for all jobs in number_list """
		job_time_list = []
		for num in number_list:
			job_time_list.append(self.jobs[num].get_time())
		return(job_time_list)
	def get_job_mems(self, number_list):
		""" Gets list of job memories for all jobs in number_list """
		job_mem_list = []
		for num in number_list:
			job_mem_list.append(self.jobs[num].get_mem())
		return(job_mem_list)
	def batch_status_change(self, number_list, new_status):
		""" Changes the status of all jobs in number_list to new_status """
		for num in number_list:
			self.jobs[num].change_status(new_status)
	def extract_contents(self, required_status_list):
		"""
		Returns a list of job information lists (using
		Job.extract_job_info) for all jobs with a status that is in
		required_status_list
		"""
		joblist_contents = []
		for current_key,current_job in self.jobs.iteritems():
			current_status = current_job.get_status()
			if current_status in required_status_list:
				current_job_extract = current_job.extract_job_info()
				joblist_contents.append(current_job_extract)
		return(joblist_contents)
	def _batch_mem_change(self, number_list, new_mem):
		"""
		Changes the mem of all jobs in number_list to new_mem
		"""
		for num in number_list:
			self.jobs[num].change_mem(new_mem)
	def _batch_time_change(self, number_list, new_time):
		"""
		Changes the time of all jobs in number_list to new_time
		"""
		for num in number_list:
			self.jobs[num].change_time(new_time)
	def _just_completed_finder(self, jobs_just_finished):
		"""
		Identifies which jobs from a list have been successfully completed
		"""
		try:
			completed_files = subprocess.check_output('ls -lrt '
				+ '"' + self.job_parameters.output_folder + '"',shell=True)
		except subprocess.CalledProcessError:
			completed_files = ''
		completed_job_list = re.findall(' ' + self.job_parameters.output_file_label + '_(\d+?)\.'
			+ self.job_parameters.output_extension,completed_files,re.DOTALL)
		completed_job_int_list = [int(x) for x in completed_job_list]
		jobs_just_completed = list(set(completed_job_int_list) & set(jobs_just_finished))
		return(jobs_just_completed)
	def _extract_latest_errorfile(self, cluster_job_submission_folder, \
		job_sub_run_output_list, job_name,job_num):
		"""
		Identifies the latest error file associated with a job, if one
		exists, and extracts its contents
		"""
		missing_job_codes = re.findall(job_name + '.e(\d+?)-' + str(job_num) + \
			'$', job_sub_run_output_list, re.MULTILINE)
		if missing_job_codes:
			job_code = str(max(map(int,missing_job_codes)))
				# jobs are assigned consecutive code numbers on cluster,
					# so since we want to look at the error file associated
					# with the most recent run of current_missing_job,
					# look for max code
			latest_errorfile = cluster_job_submission_folder + '/' + job_name + '.e' + job_code \
				+ '-' + str(job_num)
			if os.path.isfile(latest_errorfile):
				latest_errorfile_contents = open(latest_errorfile).read()
			else:
				latest_errorfile_contents = ''
		else:
			latest_errorfile_contents = ''
		return(latest_errorfile_contents)
	def _aborted_job_processor(self, max_mem, max_time, time_multiplier,
		mem_multiplier, job_num):
		""" For aborted jobs, updates job memory or time, and changes job status """
		current_job = self.jobs[job_num]
		most_recent_time = current_job.time
		most_recent_mem = current_job.mem
		if (most_recent_time == max_time and time_multiplier > 1) or \
			(most_recent_mem == max_mem and mem_multiplier > 1):
			self.batch_status_change([job_num],JobStatus.ABORTED_FOREVER)
		else:
			self.batch_status_change([job_num],JobStatus.ABORTED_TO_RESTART)
			current_updated_mem = min(most_recent_mem*mem_multiplier,max_mem)
			current_updated_time = min(most_recent_time*time_multiplier,max_time)
			self._batch_mem_change([job_num],current_updated_mem)
			self._batch_time_change([job_num],current_updated_time)
	def _parse_job_submission_output(self):
		""" Gets list of files in cluster_job_submission output folder """
		try:
			job_sub_run_output_list = subprocess.check_output('ls -lrt \'' + \
				self.job_parameters.cluster_job_submission_folder + '\'',shell=True)
		except subprocess.CalledProcessError:
			job_sub_run_output_list = ''
		return(job_sub_run_output_list)
	def _missing_job_processor(self, missing_jobs, job_submission_manager):
		"""
		Looks through any jobs that are no longer processing but have
		not been successfully completed
		Checks for error files associated with these jobs, and updates
		their status to ABORTED_FOREVER, ABORTED_TO_RESTART, or
		TO_PROCESS
		If error file exists:
			1. check whether it contains info on exceeded time or
			mem limits; if so, resubmit with max time and memory
			2. Otherwise, check whether error file empty
				if not, add field to error list
				if so, restart field
		"""
		job_sub_run_output_list = \
			self._parse_job_submission_output()
		for current_missing_job in missing_jobs:
			latest_errorfile_contents = self._extract_latest_errorfile( \
				self.job_parameters.cluster_job_submission_folder,job_sub_run_output_list, \
				self.job_parameters.name,current_missing_job)
			# If errorfile is missing or empty; or if a cluster
				# error has occurred; but job is listed as
				# having started, return job to list of jobs
				# to process
			# If errorfile contains a time limit error or a
				# memory limit error, process as an aborted
				# job
			# If errorfile contains something else, report
				# as an error
			if latest_errorfile_contents:
				error_status_dict = \
					job_submission_manager.error_status_check(\
						latest_errorfile_contents)
				if error_status_dict['cluster_error_check']:
					self.batch_status_change([current_missing_job], \
						JobStatus.TO_PROCESS)
				elif error_status_dict['memory_limit_check']:
					# updated memory should be 1.5x times previous memory
						# allotment
					# If previous memory allotment was the max allowed memory,
						# abort job forever
					time_multiplier = 1
					mem_multiplier = 1.5
					self._aborted_job_processor(self.cluster_parameters.max_mem*self.job_parameters.parallel_processors, \
						self.cluster_parameters.max_time, time_multiplier, \
						mem_multiplier, current_missing_job)
				elif error_status_dict['time_limit_check']:
					# updated time should be 2x times previous time
						# allotment
					# If previous time allotment was the max allowed time,
						# abort job forever
					time_multiplier = 2
					mem_multiplier = 1
					self._aborted_job_processor(self.cluster_parameters.max_mem*self.job_parameters.parallel_processors, \
						self.cluster_parameters.max_time, time_multiplier, \
						mem_multiplier, current_missing_job)
				elif error_status_dict['unidentified_error_check']:
					self.batch_status_change([current_missing_job], \
						JobStatus.ERROR)
			else:
				# if 'missing' job has no error file, put it back in
					# queue of jobs to process
				self.batch_status_change([current_missing_job], \
					JobStatus.TO_PROCESS)
	def autoupdate_job_status(self, job_submission_manager):
		"""
		Updates the status of any jobs that were in the 'PROCESSING'
		status based on their current status on the cluster
		Checks for jobs that are still processing
		Among the rest, look for errors, completed jobs
		parses errors, moves jobs among 'PROCESSING',
		'ABORTED_FOREVER', 'ABORTED_TO_RESTART' statuses
		looks for lost jobs and moves those to 'TO_PROCESS' list
		Gets list of jobs that should be processing
		"""
		prev_processing_list = self.get_jobs_by_status(JobStatus.PROCESSING)
		current_processing_list = job_submission_manager.job_process_finder()
		# identify jobs that no longer have 'PROCESSING' status
		jobs_just_finished = list(set(prev_processing_list)-set(current_processing_list))
		# identify jobs that no longer have 'PROCESSING' status because
			# they've just been completed, and change their status
		jobs_just_completed = self._just_completed_finder(jobs_just_finished)
		self.batch_status_change(jobs_just_completed,JobStatus.COMPLETED)
		# identify jobs that no longer have 'PROCESSING' status because
			# they've been either dropped, aborted because of time/mem
			# limites, or because they've encountered an error, and 
			# change their status
		missing_jobs = list(set(jobs_just_finished)-set(jobs_just_completed))
		if missing_jobs:
			self._missing_job_processor(missing_jobs, job_submission_manager)
	def _get_jobs_to_submit(self):
		"""
		Gets list of job numbers that have to_process or
		aborted_to_restart status
		"""
		jobs_aborted_to_restart = \
			self.get_jobs_by_status(JobStatus.ABORTED_TO_RESTART)
		jobs_to_process = \
			self.get_jobs_by_status(JobStatus.TO_PROCESS)
		job_candidate_numbers = jobs_aborted_to_restart + jobs_to_process
		return(job_candidate_numbers)
	def _group_jobs_by_time_and_mem(self, times, mems, order):
		"""
		Groups jobs based on common time and mem requirement
			times, mems, order must all by np.arrays
		Reorders job_candidate_times and job_candidate_mems
		"""
		times_sorted = times[order]
		mems_sorted = mems[order]
		# identify positions where time and mem don't change from previous position
		time_unchanged = np.diff(times_sorted) == 0
		mem_unchanged = np.diff(mems_sorted) == 0
		properties_unchanged = time_unchanged*mem_unchanged
		split_indices = np.where(np.logical_not(properties_unchanged))[0]+1
		# Split data at split_indices
		split_order = np.split(order,split_indices)
		return(split_order)
	def _group_jobs_for_sub(self, job_num_list):
		"""
		Gets the order in which jobs should be considered for submission
			Then groups jobs based on common mem and time requirement
		Gets times  and memory amts for jobs that have to_process or
		aborted_to_restart status
		"""
		job_candidate_times = \
			np.array(self.get_job_times(job_num_list))
		job_candidate_mems = \
			np.array(self.get_job_mems(job_num_list))
		# sort first by job times, then by memory
			# (lexsort input looks backwards)
		new_job_order_indices = np.lexsort((job_candidate_mems,job_candidate_times))[::-1]
		# group sorted jobs by common memory and time
		grouped_job_order_indices = \
			self._group_jobs_by_time_and_mem(job_candidate_times,job_candidate_mems, \
				new_job_order_indices)
		job_num_list_grouped = [[job_num_list[j] for j in i] for i in \
			grouped_job_order_indices]
		return(job_num_list_grouped)
	def submit_jobs(self, job_submission_manager):
		""" Submits jobs and updates job statuses accordingly """
		job_sub_candidate_nums = self._get_jobs_to_submit()
		job_sub_candidate_nums_grouped = \
			self._group_jobs_for_sub(job_sub_candidate_nums)
		# while there is still space in the cluster queue and while the
			# max number of consecutive batch submissions has not been
			# surpassed, submit jobs in batches, starting with the
			# batches requiring the most time and memory
		batch_number_remaining = job_submission_manager.max_sub_batches_in_one_run
		# get amount of space available in queue
		space_in_queue = job_submission_manager.free_job_calculator()
		for current_batch in job_sub_candidate_nums_grouped:
			# initialize batch submission manager
			submission_manager = \
				copy.deepcopy(job_submission_manager.get_submission_manager())
			submission_manager.select_jobs_to_sub(current_batch, \
				space_in_queue, batch_number_remaining)
			# update batch_number_remaining and space_in_queue
			batch_number_remaining = submission_manager.get_batches_remaining()
			space_in_queue = submission_manager.get_space_left_in_queue()
			# get list of jobs whose corresponding strings (for the job
				# alone plus a comma) are longer than max_char_num, and
				# give these 'ERROR' status
			jobs_exceeding_max_char_num = submission_manager.get_error_jobs()
			self.batch_status_change(jobs_exceeding_max_char_num,JobStatus.ERROR)
			# get list of jobs which can be submitted for this batch
			current_jobs_to_submit = submission_manager.get_submitted_jobs()
			if len(current_jobs_to_submit) > 0:
				current_job_submission_strings = submission_manager.get_submission_string_list()
				# get time and memory requirements for these jobs
					# these should all be the same, but take the max anyway
				current_batch_time = max(self.get_job_times(current_jobs_to_submit))
				current_batch_mem = max(self.get_job_mems(current_jobs_to_submit))
				# submit jobs, one batch at a time
				for current_job_sub_string in current_job_submission_strings:
					job_submission_manager.create_and_submit_batch_job( \
						current_job_sub_string, current_batch_time, \
						current_batch_mem)
				# update status of submitted jobs
				self.batch_status_change(current_jobs_to_submit,JobStatus.PROCESSING)

class CompletefileManager(object):
	"""
	Handles checking whether current set of jobs is complete, and
	writing completefile
	"""
	def __init__(self, completefile_path):
		self.completefile_path = completefile_path
		self.incomplete_status_list = [JobStatus.TO_PROCESS, \
			JobStatus.PROCESSING,JobStatus.ABORTED_TO_RESTART]
			# If there are any jobs with one of these statuses, this
				# set of jobs is not yet complete
	def add_job_list_manager(self, job_list_manager):
		self.job_list_manager = job_list_manager
	def check_completeness(self):
		""" Checks whether jobs complete; returns completeness status """
		if os.path.isfile(self.completefile_path):
			self.completeness = True
		elif hasattr(self, 'job_list_manager'):
			self._check_job_status()
		else:
			self.completeness = False
		return(self.completeness)
	def _check_job_status(self):
		"""
		Checks to see whether any jobs still need to be run or are
		running in job_list_manager
		If jobs have all been completed, writes a completefile
		"""
		incomplete_job_list = []
		for current_job_status in self.incomplete_status_list:
			current_job_list = \
				self.job_list_manager.get_jobs_by_status(current_job_status)
			incomplete_job_list.extend(current_job_list)
		incomplete_job_number = len(incomplete_job_list)
		if incomplete_job_number == 0:
			self.completeness = True
			# write completeness file
			open(self.completefile_path,'a').close()
		else:
			self.completeness = False

class SubmissionStringProcessor(object):
	"""
	Creates a code submission string appropriate to the programming
	environment (module) being used
	"""
	def __init__(self, module, key_list, value_list, code_name):
		self._submission_string_generator(key_list,value_list,module,code_name)
	def get_code_run_input(self):
		return(self.code_sub_input_processor)
	def _submission_string_generator(self, key_list, value_list, module, \
		code_name):
		if module == 'matlab':
			code_sub_input_processor = \
				MatlabInputProcessor(code_name)
		elif module == 'r':
			code_sub_input_processor = \
				RInputProcessor(code_name)
		else:
			raise ValueError('%s is not an available module at the moment' \
				% (module))
		converted_key_string = \
			code_sub_input_processor.convert_mixed_list(key_list)
		converted_value_string = \
			code_sub_input_processor.convert_mixed_list(value_list)
		code_input_arg_list = [converted_key_string, converted_value_string]
		code_sub_input_processor.set_code_run_argument_string(code_input_arg_list)
		self.code_sub_input_processor = code_sub_input_processor

class TrackfileManager(object):
	"""
	Handles writing and reading of trackfiles, which store information
	on current job status of each job
	"""
	def __init__(self, job_parameters, cluster_parameters, trackfile_folder):
		self.trackfile_path = \
			os.path.join(trackfile_folder, \
				('trackfile_' + job_parameters.name + '.csv'))
		self.summaryfile_path = os.path.join(trackfile_folder, \
			('summary_trackfile_'+job_parameters.name+'.csv'))
		self.job_parameters = copy.deepcopy(job_parameters)
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
	def get_summaryfile_path(self):
		""" Returns summaryfile path """
		return(self.summaryfile_path)
	def get_trackfile_path(self):
		""" Returns trackfile path """
		return(self.trackfile_path)
	def read_trackfile(self):
		"""
		Reads a csv containing the status of each job in job_list, as
		well as the time and memory allocated to it, and returns a
		JobListManager object
		"""
		current_job_list = []
		with open(self.trackfile_path, 'rU') as trackfile_opened:
			trackfile_contents = list(csv.reader(trackfile_opened))
			for row_list in trackfile_contents[1:]:
				current_job = Job(int(row_list[0]), row_list[1], float(row_list[2]), \
					float(row_list[3]))
				current_job_list.append(current_job)
		return JobListManager(current_job_list,self.job_parameters,self.cluster_parameters)
	def write_output_files(self, job_list_manager):
		"""
		Writes trackfile and summaryfile for current object for each
		job in job_list_manager, which is a JobListManager object
		"""
		self._write_trackfile(job_list_manager)
		self._write_summaryfile(job_list_manager)
	def _get_complete_status_list(self):
		""" Gets entire list of possible statuses from JobStatus class """
		status_list = [getattr(JobStatus, attr) for attr in vars(JobStatus) \
			if not attr.startswith("__")]
		return(status_list)
	def _write_trackfile(self, job_list_manager):
		"""
		Writes a csv containing the status of each job in
		job_list_manager, as well as the time and memory allocated to it
		"""
		complete_status_list = self._get_complete_status_list()
		job_list_contents = job_list_manager.extract_contents(complete_status_list)
		with open(self.trackfile_path, 'wb') as trackfile_opened:
			trackfile_writer = csv.writer(trackfile_opened)
			trackfile_writer.writerow(['Job Number','Job Status','Time Allocated to Job','Memory Allocated to Job'])
			for current_job in job_list_contents:
				trackfile_writer.writerow(current_job)
	def _write_summaryfile(self, job_list_manager):
		"""
		Writes a csv file counting the number of jobs of each status in
		job_list_manager
		"""
		with open(self.summaryfile_path, 'wb') as summaryfile_opened:
			summaryfile_writer = csv.writer(summaryfile_opened)
			summaryfile_writer.writerow(['Status','# of jobs','job indices'])
			complete_status_list = self._get_complete_status_list()
			for status in complete_status_list:
				indices = job_list_manager.get_jobs_by_status(status)
				current_n = len(indices)
				indices_concat = ';'.join(str(x) for x in indices)
				summaryfile_writer.writerow([status,current_n,indices_concat])

class RInputProcessor(object):
	""" Generates proper submission string properties for submitting R jobs """
	def __init__(self, code_name):
		self.code_name = code_name
		self.set_code_run_argument_string([])

class MatlabInputProcessor(object):
	"""
	Generates proper submission string properties for submitting
	matlab jobs
	"""
	def __init__(self, code_name):
		self.code_name = code_name
		self.set_code_run_argument_string([])
	def set_code_run_argument_string(self, code_input_arguments):
		""" Gets a list of arguments and creates a submission string from them """
		self.code_run_arguments = '\"\",\"\"'.join(code_input_arguments)
	def set_full_code_run_string(self, cluster_architecture):
		if cluster_architecture == 'slurm':
			self.code_run_string = \
				'matlab -nodisplay -nosplash -nodesktop -r \'' + \
				self.code_name + '(\'\"' + self.code_run_arguments + \
				'\"\");exit\"'
		elif cluster_architecture == 'macosx':
			self.code_run_string = \
				'matlab -nodisplay -nosplash -nodesktop -r \'try ' + \
				self.code_name + '(\'\"' + self.code_run_arguments + \
				"\"\"); catch error_contents; fprintf('There was an error:\\n');" + \
				" fprintf(1,'The identifier was:\\n%s\\n',error_contents.identifier);" + \
				" fprintf(1,'The error message was:\\n%s\\n',error_contents.message); end; exit\""
	def get_code_run_string(self):
		return(self.code_run_string)
	def convert_mixed_list(self, current_list):
		"""
		Checks type of every element of current_list and converts it
		to a string
		Joins elements into single string
		"""
		converted_by_part_list = []
		for current_value in current_list:
			current_val_converted = self.convert_val(current_value)
			converted_by_part_list.append(current_val_converted)
		converted_list = '{' + ','.join(converted_by_part_list) + '}'
		return(converted_list)
	def convert_str_list(self,current_list):
		converted_strings = [self.convert_str(i) for i in current_list]
		converted_list = '{' + ','.join(converted_strings) + '}'
		return(converted_list)
	def convert_int_list(self,current_list):
		converted_ints = [self.convert_int(i) for i in current_list]
		converted_list = '[' + ','.join(converted_ints) + ']'
		return(converted_list)
	def convert_float_list(self,current_list):
		converted_floats = [self.convert_float(i) for i in current_list]
		converted_list = '[' + ','.join(converted_floats) + ']'
		return(converted_list)
	def convert_bool_list(self,current_list):
		converted_bools = [self.convert_bool(i) for i in current_list]
		converted_list = '[' + ','.join(converted_bools) + ']'
		return(converted_list)
	def convert_str(self, current_str):
		converted_str = '\'' + current_str + '\''
		return(converted_str)
	def convert_int(self, current_int):
		converted_int = str(current_int)
		return(converted_int)
	def convert_float(self, current_float):
		converted_float = str(current_float)
		return(converted_float)
	def convert_bool(self, current_bool):
		if current_bool:
			converted_bool = 'true'
		else:
			converted_bool = 'false'
		return(converted_bool)
	def convert_none(self):
		return('NaN')
	def convert_any_list(self, current_value):
		if all(isinstance(temp_val,bool) for temp_val in current_value):
			converted_value = self.convert_bool_list(current_value)
		elif all(isinstance(temp_val,int) for temp_val in current_value):
			converted_value = self.convert_int_list(current_value)
		elif all(isinstance(temp_val,float) for temp_val in current_value):
			converted_value = self.convert_float_list(current_value)
		elif all(isinstance(temp_val,str) for temp_val in current_value):
			converted_value = self.convert_str_list(current_value)
		else:
			converted_value = self.convert_mixed_list(current_value)
		return(converted_value)
	def convert_val(self,current_value):
		if current_value is None:
			converted_value = self.convert_none()
		elif isinstance(current_value, np.ndarray):
			current_value_listified = current_value.tolist()
			converted_value = self.convert_any_list(current_value_listified)
		elif isinstance(current_value, list):
			converted_value = self.convert_any_list(current_value)
		elif isinstance(current_value, bool):
			converted_value = self.convert_bool(current_value)
		elif isinstance(current_value, int):
			converted_value = self.convert_int(current_value)
		elif isinstance(current_value, float):
			converted_value = self.convert_float(current_value)
		elif isinstance(current_value, basestring):
			converted_value = self.convert_str(current_value)
		else:
			print('Error: Trying to convert list element of unrecognized type in submission string conversion:')
			print(current_value)
		return(converted_value)

class CodeSubmitter(object):
	'''
	Class that combines inputs necessary for submission of jobs using
	some external_function, and submits those jobs
	'''
	def __init__(self, cluster_parameters, cluster_folders, completefile, job_name, \
		job_numbers, module, parallel_processors, experiment_folder, output_extension, code_name, \
		additional_beginning_lines_in_job_sub, additional_end_lines_in_job_sub, \
		initial_sub_time, initial_sub_mem, output_path, output_file_label):
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
		self.cluster_folders = copy.deepcopy(cluster_folders)
		self.parallel_processors = parallel_processors
		# The following should likely be generated in the classes
			# inheriting from CodeSubmitter and then passed to
			# CodeSubmitter's __init__ method
		self.completefile = completefile
		self.job_name = job_name
		self.job_numbers = job_numbers
		self.module = module
		self.output_extension = output_extension
		self.code_name = code_name
		self.additional_beginning_lines_in_job_sub = additional_beginning_lines_in_job_sub
		self.additional_end_lines_in_job_sub = additional_end_lines_in_job_sub
		self.initial_sub_time = initial_sub_time
		self.initial_sub_mem = initial_sub_mem
		self.output_path = output_path
		self.output_file_label = output_file_label
		self.experiment_folder = experiment_folder
		self._create_code_run_input()
	def get_completefile_path(self):
		return(self.completefile)
	def _create_code_run_input_lists(self):
		'''
		Creates list of keys and their values to be submitted to
		external code being run
		'''
		self.key_list = []
		self.value_list = []
	def _create_code_run_input(self):
		self._create_code_run_input_lists()
		# process key_list and value_list into a submission string
		submission_string_processor = \
			SubmissionStringProcessor(self.module, \
				self.key_list, self.value_list, self.code_name)
		self.code_run_input = submission_string_processor.get_code_run_input()
	def run_job_submission(self):
		# handles submission of the job
		cluster_job_submission_folder = \
			self.cluster_folders.get_path('cluster_job_submission_path')
		trackfile_folder = \
			self.cluster_folders.get_path('trackfile_path')
		# set up and run batch jobs
		job_flow_handler(self.job_name, self.job_numbers, self.initial_sub_time, \
			self.initial_sub_mem, self.cluster_parameters, self.output_path, self.output_extension, \
			self.output_file_label, cluster_job_submission_folder, self.experiment_folder, \
			self.module, self.code_run_input, self.additional_beginning_lines_in_job_sub, \
			self.additional_end_lines_in_job_sub, self.parallel_processors, self.completefile,
			trackfile_folder)

#######################################################
def parse_setup_file(filename):
	"""
	Parses setup file and creates a list of the parameters therein
	"""
	f = open(filename)
	lines = f.readlines()
	f.close()
	lines.pop(0)
	parameters = []
	for l in lines:
		spl = l.rstrip().split(",")
		if len(spl) != 4:
			raise ValueError("needs 4 values: " + l)
		if not spl[1] == spl[2] == '':
			parameters.append(Parameter(spl[0], spl[1], spl[2], spl[3].lstrip().rstrip()))
	return Parameters(parameters)

def _create_job_list(job_name, job_numbers, initial_time, initial_mem, \
		job_parameters, cluster_parameters):
	"""
	Create a list of jobs sharing name, and number from 'numbers' list
	All new jobs will have 'TO_PROCESS' status
	"""
	current_job_list = []

	for n in job_numbers:
		if type(n) is not int:
			raise TypeError("numbers contains non-integers: " + l)
		current_job = Job(n, JobStatus.TO_PROCESS, initial_time, \
			initial_mem*job_parameters.parallel_processors)
		current_job_list.append(current_job)
	return JobListManager(current_job_list,job_parameters,cluster_parameters)

def job_flow_handler(job_name, job_numbers, initial_time, initial_mem, \
	cluster_parameters, output_folder, output_extension, output_file_label, \
	cluster_job_submission_folder, experiment_folder, module, code_run_input, \
	additional_beginning_lines_in_job_sub, additional_end_lines_in_job_sub, \
	parallel_processors, completefile_path, trackfile_folder):
	"""
	Handles entire flow of job, from creation of new trackfile to
	submission of jobs to updating trackfile
	"""
	completefile_manager = CompletefileManager(completefile_path)
	# check completeness and only run other code if job is not complete
	jobs_complete = completefile_manager.check_completeness()
	if not jobs_complete:
		########### ADD COMPLETENESS UPDATE BELOW
		job_parameters = JobParameters(job_name, output_folder, \
			output_extension, output_file_label, cluster_job_submission_folder, experiment_folder, \
			module, code_run_input, additional_beginning_lines_in_job_sub, \
			additional_end_lines_in_job_sub, parallel_processors)
		current_trackfile = TrackfileManager(job_parameters, \
			cluster_parameters, trackfile_folder)
		# check whether trackfile exists; if so, use it to get job data and
			# update statuses; otherwise, create new trackfile
		current_trackfile_path = current_trackfile.get_trackfile_path()
		# create job_submission_manager
		job_submission_manager = cluster_parameters.get_job_sub_manager()
		job_submission_manager.set_job_parameters(job_parameters)
		if os.path.isfile(current_trackfile_path):
			# retrieve job list, update it, submit jobs
			# job_list_manager is a JobListManager object
			job_list_manager = current_trackfile.read_trackfile()
			# update status of all jobs in job_list_manager
			job_list_manager.autoupdate_job_status(job_submission_manager)
		else:
			# create job list (JobListManager object)
			job_list_manager = _create_job_list(job_name, job_numbers, initial_time, \
				initial_mem, job_parameters, cluster_parameters)
		# submit whatever jobs possible
		job_list_manager.submit_jobs(job_submission_manager)
		# update trackfile and summaryfile
		current_trackfile.write_output_files(job_list_manager)
		# update completefile
		completefile_manager.add_job_list_manager(job_list_manager)
		completefile_manager.check_completeness()
