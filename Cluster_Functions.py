#!/usr/bin/python

# Contains objects needed for running jobs on cluster and keeping track
	# of their progress

import os
import csv
import subprocess
import numpy
import copy
import re

class Job(object):
	# job object that stores properties of individual jobs in queue
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
	# enum for job status
	TO_PROCESS = 1
	PROCESSING = 2
	COMPLETED = 3
	ABORTED_FOREVER = 4
	ABORTED_TO_RESTART = 5
	ERROR = 6

class ClusterParameters(object):
	# holds parameters general to the system
	def __init__(self,parameter_list):
		self.code_path = parameter_list["code_folder"]
		self.home_path = parameter_list["home_folder"]
		self.composite_data_path = parameter_list["composite_data_folder"]
		self.max_mem = parameter_list["max_mem"]
		self.max_time = parameter_list["max_time"]
		self.username = parameter_list["username"]
		self.user_email = parameter_list["user_email"]
		self.starting_mem = parameter_list["starting_mem"]
		self.starting_time = parameter_list["starting_time"]
		self.temp_storage_path = parameter_list["temp_storage_folder"]
		self.max_char_num = parameter_list["max_char_num"]
		self.cluster_architecture = parameter_list["cluster_architecture"].lower()
		self._set_cluster_architecture_properties()
		self.current_mem = self.starting_mem
		self.current_time = self.starting_time
	def _set_cluster_architecture_properties(self):
		# sets properties related to the cluster architecture
		if self.cluster_architecture == 'slurm':
			self.within_batch_counter = '$SLURM_ARRAY_TASK_ID'
		elif self.cluster_architecture == 'macosx':
			self.within_batch_counter = '$ARRAY_TASK_ID'
		else:
			print('Error! Did not enter recognized cluster architecture')
	def set_current_time(self,new_time):
		self.current_time = new_time
	def set_current_mem(self,new_mem):
		self.current_mem = new_mem

class FolderManager(object):
	def __init__(self,cluster_parameters,experiment_folder_name):
		self.experiment_path = \
			os.path.join(cluster_parameters.composite_data_path,experiment_folder_name)
		self.trackfile_path = os.path.join(self.experiment_path,'trackfiles')
		self.completefile_path = os.path.join(self.experiment_path,'completefiles')
		self.cluster_job_submission_path = os.path.join(cluster_parameters.temp_storage_path, \
			experiment_folder_name,'cluster_job_submission_folder')
		self._set_up_folders()
	def _set_up_folders(self):
		setup_complete_file = os.path.join(self.completefile_folder,'folder_setup_complete.txt')
		if not os.path.isfile(setup_complete_file):
			new_directory_list = (self.trackfile_path,self.completefile_path, \
				self.cluster_job_submission_path)
			for current_new_directory in new_directory_list:
				if not os.path.isdir(current_new_directory):
					os.makedirs(current_new_directory)
			open(setup_complete_file,'a').close()

class JobParameters(object):
	# holds parameters of the job currently being run
	def __init__(self, name, output_folder, output_extension, output_filename, \
		cluster_job_submission_folder,experiment_folder, module, code_run_string, \
		additional_beginning_lines_in_sbatch, additional_end_lines_in_sbatch, \
		parallel_processors):
		self.name = name
		self.output_folder = output_folder
		self.output_extension = output_extension
		self.output_filename = output_filename
		self.cluster_job_submission_folder = cluster_job_submission_folder
		self.experiment_folder = experiment_folder
		self.module = module
		self.module_path = self._get_latest_module_path(module)
		self.code_run_string = code_run_string
		self.additional_beginning_lines_in_sbatch = \
			self._string_to_list(additional_beginning_lines_in_sbatch)
		self.additional_end_lines_in_sbatch = \
			self._string_to_list(additional_end_lines_in_sbatch)
		self.parallel_processors = parallel_processors
	def _string_to_list(parameter):
		# checks whether parameter is entered as a string or list, and
			# if a string, converts to a list holding that string; if
			# neither, returns an error
		if isinstance(parameter, list):
			output_parameter = parameter
		elif isinstance(parameter, basestring):
			output_parameter = [parameter]
		else:
			print('Error! Non-string, non-list passed to JobParameters')
		return(output_parameter)
	def _get_latest_module_path(module):
		try:
			module_output = subprocess.check_output('module avail', shell=True)
		except subprocess.CalledProcessError:
			module_output = ''
		relevant_module_list = re.findall(' (' + module.lower() + os.sep + '.+?)\\n', module_output)
		# latest module is the last one in the list
		latest_module_path = re.findall('(' + module.lower() + os.sep + '\S+)', \
			relevant_module_list[-1])[0]
		return(latest_module_path)

class BatchSubmissionManager(object):
	# Holds parameters for current job submission session, and manages
		# which jobs from current batch are submitted
	# Algorithm here assumes space_in_queue is typically more limiting
		# than max_char_num and batches_remaining
	def __init__(self,job_list,max_char_num,space_in_queue,batches_remaining):
		self.space_in_queue = space_in_queue
		self.batches_remaining = batches_remaining
		self.job_list = job_list
		self.jobs_to_submit = []
		self.job_submission_string_list = ['']
		self.max_char_num = max_char_num
		self.remaining_chars = max_char_num
		self.jobs_to_error = []
	def get_error_jobs(self):
		# returns jobs whose corresponding strings (for the job alone
			# plus a comma) are longer than max_char_num
		return(self.jobs_to_error)
	def get_submission_string_list(self):
		return(self.job_submission_string_list)
	def get_submitted_jobs(self):
		return(self.jobs_to_submit)
	def get_batches_remaining(self):
		return(self.batches_remaining)
	def get_space_left_in_queue(self):
		return(self.space_in_queue)
	def _update_space_in_queue(self,job_number):
		self.space_in_queue = self.space_in_queue - job_number
	def _update_batches_remaining(self,batch_number):
		self.batches_remaining = batches_remaining - batch_number
	def _update_remaining_chars(self,char_number):
		self.remaining_chars = self.remaining_chars - char_number
	def _add_to_submission_string(self,new_string):
		self.job_submission_string_list[-1] = self.job_submission_string_list[-1]+new_string
	def _add_jobs_to_submission_list(self,job_sublist,job_string):
		self._add_to_submission_string(job_string)
		self.jobs_to_submit.extend(job_sublist)
		self._update_space_in_queue(len(job_sublist))
		self._update_remaining_chars(len(job_string))
	def _reset_batch(self):
		self.remaining_chars = max_char_num
		self.job_submission_string_list.append('')
	def _consecutive_parser(job_list, stepsize = 1):
		# Splits list of integers into list of lists of consecutive nums
		# Returns the list of consecutive integer lists, as well as a
			# list of lists of the indices of each member of each
			# consecutive integer lists in the original job_list
		np_job_list = numpy.array(job_list)
		sorted_indices = numpy.argsort(np_job_list)
		sorted_job_list = np_job_list[sorted_indices]
		# Find position of indices where the number is higher than
			# previous number + stepsize
		split_indices = numpy.where(numpy.diff(sorted_job_list) != stepsize)[0]+1
		# Split data at split_indices
		split_job_list = numpy.split(sorted_job_list,split_indices)
		return(split_job_list)
	def _job_list_string_converter(self,consecutive_job_list):
		# Converts a consecutive list of job numbers into a job
			# submission string for SLURM
		# test that consecutive_job_list is actually consecutive
		test_list = self._consecutive_parser(consecutive_job_list)
		if len(test_list) > 1:
			print('Error! BatchSubmissionManager._job_list_string_converter provided unparsed list')
		# Denote consecutive job sublists by dashes between the min
			# and max of the sublist; separate job sublists by commas
		current_job_num = len(consecutive_job_list)
		if current_job_num==1:
			current_job_string = (str(consecutive_job_list[0])+',')
		else:
			current_job_string = (str(min(consecutive_job_list))+'-'+
				str(max(consecutive_job_list))+',')
		return(current_job_string)
	def _job_list_parser(self,job_sublist):
		# Taking in a list of jobs, reorganize them so that consecutive
			# jobs can be called as intervals (e.g. '5-8' for
			# '5,6,7,8'), and find the optimal arrangement that
			# maximizes job number so that the character count for the
			# list of jobs doesn't exceed max_char_num
		# Assumes that max_char_num is sufficiently large
			# i.e. >(2*(max number of digits in a job number)+2)
		if self.space_in_queue < len(job_sublist):
			# reorder job_sublist and truncate
			job_sublist = numpy.sort(job_sublist)[0:self.space_in_queue]
		# create a string that can be used to submit those jobs
			# to SLURM, and count the number of characters in that string
		current_jobs_string = self._job_list_string_converter(job_sublist)
		if len(current_jobs_string) < self.remaining_chars:
			self._add_jobs_to_submission_list(job_sublist,current_job_string)
		else:
			# start new batch, try adding jobs to list again
			self._update_batches_remaining(1)
			if self.batches_remaining > 0:
				self._reset_batch()
				if len(current_jobs_string) < self.remaining_chars:
					self._add_jobs_to_submission_list(job_sublist,current_job_string)
				else:
					# try splitting up job_list one-by-one
					for current_job in job_sublist:
						current_jobs_string = self._job_list_string_converter([current_job])
						if len(current_jobs_string) < self.remaining_chars:
							self._add_jobs_to_submission_list(job_sublist,current_job_string)
						else:
							self.jobs_to_error.extend(current_job)
	def select_jobs_to_sub(self):
		# Based on the amount of space available in the queue, splits jobs
			# into ones that can be run now vs those that have to be run
			# later
		#################
		# SLURM-based cluster doesn't appear to slow down as a penalty for
			# submitting jobs one-by-one so don't reset
			# max_sbatches_in_one_run if submitting jobs one-by-one
		#################
		# Take the maximum number of jobs that can be submitted from
			# initial_job_list, and parse them to ensure that
			#	1. they're listed in the most efficient SLURM-readable
			#		format
			#	2. this list doesn't exceed the max number of chars
			#		allowed by SLURM
		################# maybe correct? #################	
		job_list_split = self._consecutive_parser(self.job_list)
		job_number_list = [len(i) for i in job_list_split]
		order_by_job_num = numpy.argsort(job_number_list)
		job_list_split_sorted = job_list_split[order_by_job_num]
		for current_job_list in job_list_split_sorted:
			if self.space_in_queue > 0 and self.batches_remaining > 0:
				self._job_list_parser(current_job_list)
		self._update_batches_remaining(1)

class JobListManager(object):
	# holds list of jobs corresponding to a single 'name'
	# updates current job status
	def __init__(self, jobs, job_parameters,cluster_parameters):
		self.jobs = {}
		for j in jobs:
			self.jobs[j.number] = copy.deepcopy(j)
		self.job_parameters = copy.deepcopy(job_parameters)
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
	def get_job_name(self):
		# returns the name of the jobs
		job_name = self.job_parameters.name
		return(job_name)
	def _get_status_list(self):
		# gets current status of every job
		# ??? Do I need this function? ???
		status_list=[self.jobs[num].status for num in self.jobs]
		return status_list
	def get_jobs_by_status(self,status):
		# gets list of jobs with a specific status
		job_subset_list = []
		for num in self.jobs:
			if self.jobs[num].status == status:
				job_subset_list.extend(num)
		return job_subset_list
	def get_job_times(self,number_list):
		# gets list of job times for all jobs in number_list
		job_time_list = []
		for num in number_list:
			job_time_list.extend(self.jobs[num].get_time())
		return(job_time_list)
	def get_job_mems(self,number_list):
		# gets list of job memories for all jobs in number_list
		job_mem_list = []
		for num in number_list:
			job_mem_list.extend(self.jobs[num].get_mem())
		return(job_mem_list)
	def batch_status_change(self, number_list, new_status):
		# changes the status of all jobs in number_list to new_status
		for num in number_list:
			self.jobs[num].change_status(new_status)
	def extract_contents(self, required_status_list):
		joblist_contents = []
		for current_job in self.jobs:
			current_status = curent_job.get_status()
			if current_status in required_status_list:
				current_job_extract = current_job.extract_job_info()
				joblist_contents.extend(current_job_extract)
		return(joblist_contents)
	def _batch_mem_change(self, number_list, new_mem):
		# changes the mem of all jobs in number_list to new_mem
		for num in number_list:
			self.jobs[num].change_mem(new_mem)
	def _batch_time_change(self, number_list, new_time):
		# changes the time of all jobs in number_list to new_time
		for num in number_list:
			self.jobs[num].change_time(new_time)
	def _job_process_finder(self):
		# Identify jobs that are still being run by SLURM
		# Return list of jobs currently being processed
		try:
			qstat_output = subprocess.check_output('squeue -u ' + self.cluster_parameters.username + ' -r -n '
				+ self.job_parameters.current_job_name + ' | egrep " PD | CG | R " ',shell=True)
		except subprocess.CalledProcessError:
			qstat_output = ''
		jobs_still_in_processing = re.findall('^\s+\d+_(\d+)\s',qstat_output,re.MULTILINE)
		return(jobs_still_in_processing)
	def _just_completed_finder(self,jobs_just_finished):
		# Identify which jobs from a list have been successfully completed
		try:
			completed_files = subprocess.check_output('ls -lrt '
				+ self.job_parameters.output_folder,shell=True)
		except subprocess.CalledProcessError:
			completed_files = ''
		completed_job_list = re.findall(' ' + self.job_parameters.output_filename + '_(\d+?)\.'
			+ self.job_parameters.output_extension,completed_files,re.DOTALL)
		jobs_just_completed = list(set(completed_job_list) & set(jobs_just_finished))
		return(jobs_just_completed)
	def _parse_sbatch_output(cluster_job_submission_folder):
		# gets list of files in cluster_job_submission output folder
		try:
			sbatch_run_output_list = subprocess.check_output('ls -lrt ' + cluster_job_submission_folder,shell=True)
		except subprocess.CalledProcessError:
			sbatch_run_output_list = ''
		return(sbatch_run_output_list)
	def _extract_latest_errorfile(cluster_job_submission_folder,sbatch_run_output_list,job_name,job_num):
		# identifies the latest error file associated with a job,
			# if one exists, and extracts its contents
		missing_job_codes = re.findall(job_name + '.e(\d+?)-' + job_num + '$',
			sbatch_run_output_list, re.MULTILINE)
		if missing_job_codes:
			job_code = str(max(map(int,missing_job_codes)))
				# jobs are assigned consecutive code numbers on cluster,
					# so since we want to look at the error file associated
					# with the most recent run of current_missing_job,
					# look for max code
			latest_errorfile = cluster_job_submission_folder + '/' + job_name + '.e' + job_code \
				+ '-' + job_num
			if os.path.isfile(latest_errorfile):
				latest_errorfile_contents = open(latest_errorfile).read()
			else:
				latest_errorfile_contents = ''
		else:
			latest_errorfile_contents = ''
		return(latest_errorfile_contents.lower())
	def _aborted_job_processor(self, max_mem, max_time, time_multiplier,
		mem_multiplier, job_num):
		# For aborted jobs, update job memory or time, and change job status
		current_job = self.jobs[job_num]
		most_recent_time = current_job.time
		most_recent_mem = current_job.mem
		if most_recent_time == max_time or most_recent_mem == max_mem:
			self.batch_status_change([job_num],JobStatus.ABORTED_FOREVER)
		else:
			self.batch_status_change([job_num],JobStatus.ABORTED_TO_RESTART)
			current_updated_mem = min(most_recent_mem*mem_multiplier,max_mem)
			current_updated_time = min(most_recent_time*time_multiplier,max_time)
			self._batch_mem_change([job_num],current_updated_mem)
			self._batch_time_change([job_num],current_updated_time)
	def _missing_job_processor(self,missing_jobs):
		# Looks through any jobs that are no longer processing but have
			# not been successfully completed
		# Checks for error files associated with these jobs, and updates
			# their status to ABORTED_FOREVER, ABORTED_TO_RESTART, or
			# TO_PROCESS
		# if error file exists:
			#	1. check whether it contains info on exceeded time or
			#		mem limits; if so, resubmit with max time and memory
			#	2. Otherwise, check whether error file empty
			#		if not, add field to error list
			#		if so, restart field
		sbatch_run_output_list = \
			self._parse_sbatch_output(self.job_parameters.cluster_job_submission_folder)
		for current_missing_job in missing_jobs:
			latest_errorfile_contents = self._extract_latest_errorfile( \
				self.job_parameters.cluster_job_submission_folder,sbatch_run_output_list, \
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
				time_limit_check = 'time limit' in latest_errorfile_contents
				memory_limit_check = 'memory limit' in latest_errorfile_contents
				cluster_error_check = 'bus error' in latest_errorfile_contents \
					or 'fatal error on startup' in latest_errorfile_contents \
					or 'reload' in latest_errorfile_contents \
					or 'MatlabException' in latest_errorfile_contents
				if cluster_error_check:
					self.batch_status_change([current_missing_job], \
						JobStatus.TO_PROCESS)
				elif memory_limit_check:
					# updated memory should be 1.5x times previous memory
						# allotment
					# If previous memory allotment was the max allowed memory,
						# abort job forever
					time_multiplier = 1
					mem_multiplier = 1.5
					self._aborted_job_processor(self.cluster_parameters.max_mem, \
						self.cluster_parameters.max_time, time_multiplier, \
						mem_multiplier, job_num)
				elif time_limit_check:
					# updated time should be 2x times previous time
						# allotment
					# If previous time allotment was the max allowed time,
						# abort job forever
					time_multiplier = 2
					mem_multiplier = 1
					self._aborted_job_processor(self.cluster_parameters.max_mem, \
						self.cluster_parameters.max_time, time_multiplier, \
						mem_multiplier, job_num)
				else:
					self.batch_status_change([current_missing_job], \
						JobStatus.ERROR)
			else:
				# if 'missing' job has no error file, put it back in
					# queue of jobs to process
				self.batch_status_change([current_missing_job], \
					JobStatus.TO_PROCESS)
	def autoupdate_job_status(self):
		# update the status of any jobs that were in the 'PROCESSING'
			# status based on their current status on the cluster
		# check for jobs that are still processing
		# among the rest, look for errors, completed
			# parse errors, move jobs among 'PROCESSING', 'ABORTED_FOREVER', 'ABORTED_TO_RESTART' statuses
			# look for lost jobs and move those to 'TO_PROCESS' list
		# get list of jobs that should be processing
		prev_processing_list = self.get_jobs_by_status(JobStatus.PROCESSING)
		current_processing_list = self._job_process_finder()
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
			self._missing_job_processor(missing_jobs)
	def _get_jobs_to_submit(self):
		# get list of job numbers that have to_process or
			# aborted_to_restart status
		jobs_aborted_to_restart = \
			self.get_jobs_by_status(JobStatus.ABORTED_TO_RESTART)
		jobs_to_process = \
			self.get_jobs_by_status(JobStatus.TO_PROCESS)
		job_candidate_numbers = jobs_aborted_to_restart + jobs_to_process
		return(job_candidate_numbers)
	def _group_jobs_by_time_and_mem(times,mems,order):
		# Group jobs based on common time and mem requirement
			# times, mems, order must all by numpy.arrays
		# reorder job_candidate_times and job_candidate_mems
		times_sorted = times[order]
		mems_sorted = mems[order]
		# identify positions where time and mem don't change from previous position
		time_unchanged = numpy.diff(times_sorted) == 0
		mem_unchanged = numpy.diff(mems_sorted) == 0
		properties_unchanged = time_unchanged*mem_unchanged
		split_indices = numpy.where(numpy.logical_not(properties_unchanged))[0]+1
		# Split data at split_indices
		split_order = numpy.split(order,split_indices)
		return(split_order)
	def _group_jobs_for_sub(self,job_num_list):
		# Get the order in which jobs should be considered for submission
			# Then group jobs based on common mem and time requirement
		# get times  and memory amts for jobs that have to_process or
			# aborted_to_restart status
		job_candidate_times = \
			numpy.array(self.get_job_times(job_num_list))
		job_candidate_mems = \
			numpy.array(self.get_job_mems(job_num_list))
		# sort first by job times, then by memory
			# (lexsort input looks backwards)
		new_job_order_indices = numpy.lexsort((job_candidate_mems,job_candidate_times))[::-1]
		# group sorted jobs by common memory and time
		grouped_job_order_indices = \
			self._group_jobs_by_time_and_mem(job_candidate_times,job_candidate_mems, \
				new_job_order_indices)
		job_num_list_grouped = [[job_num_list[j] for j in i] for i in \
			grouped_job_order_indices]
		return(job_num_list_grouped)
	def submit_jobs(self,cluster_submission_manager):
		# Submits jobs and updates job statuses accordingly
		job_sub_candidate_nums = self._get_jobs_to_submit()
		job_sub_candidate_nums_grouped = \
			self._group_jobs_for_sub(job_sub_candidate_nums)
		# while there is still space in the cluster queue and while the
			# max number of consecutive batch submissions has not been
			# surpassed, submit jobs in batches, starting with the
			# batches requiring the most time and memory
		batch_number_remaining = self.max_sub_batches_in_one_run()
		# get amount of space available in queue
		space_in_queue = cluster_submission_manager.free_job_calculator()
		for current_batch in job_sub_candidate_nums_grouped:
			# initialize batch submission manager
			submission_manager = BatchSubmissionManager(current_batch, \
				self.max_char_num, space_in_queue, batch_number_remaining)
			submission_manager.select_jobs_to_sub()
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
			current_job_submission_strings = submission_manager.get_submission_string_list()
			# get time and memory requirements for these jobs
				# these should all be the same, but take the max anyway
			current_batch_time = max(self.get_job_times(current_jobs_to_submit))
			current_batch_mem = max(self.get_job_mems(current_jobs_to_submit))
			# submit jobs
			cluster_submission_manager.create_and_submit_batch_job(job_list_string,job_time,job_mem)
			# update status of submitted jobs
			self.batch_status_change(current_jobs_to_submit,JobStatus.PROCESSING)

class CompletefileManager(object):
	# Handles checking if current set of jobs is complete, and writing completefile
	def __init__(self, completefile_path):
		self.completefile_path = completefile_path
		self.incomplete_status_list = [JobStatus.TO_PROCESS, \
			JobStatus.PROCESSING,JobStatus.ABORTED_TO_RESTART]
			# If there are any jobs with one of these statuses, this
				# set of jobs is not yet complete
	def add_job_list_manager(self,job_list_manager):
		self.job_list_manager = job_list_manager
	def check_completeness(self):
		# check whether jobs complete
		# return completeness status
		if os.path.isfile(self.completefile_path):
			self.completeness = True
		elif hasattr(self, 'job_list_manager'):
			self._check_job_status()
		else:
			self.completeness = False
		return(self.completeness)
	def _check_job_status(self):
		# checks to see whether any jobs still need to be run or are
			# running in job_list_manager
		# If jobs have all been completed, writes a completefile
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


class TrackfileManager(object):
	# Handles writing and reading of trackfile
	def __init__(self, job_parameters, cluster_parameters):
		self.trackfile_path = os.path.join(job_parameters.experiment_folder, \
			'trackfiles',('trackfile_'+job_parameters.name+'.csv'))
		self.summaryfile_path = os.path.join(job_parameters.experiment_folder, \
			'trackfiles',('summary_trackfile_'+job_parameters.name+'.csv'))
		self.job_parameters = copy.deepcopy(job_parameters)
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
	def get_summaryfile_path(self):
		# returns summaryfile path
		return(self.summaryfile_path)
	def get_trackfile_path(self):
		# returns trackfile path
		return(self.trackfile_path)
	def read_trackfile(self):
		# reads a csv containing the status of each job in job_list,
			# as well as the time and memory allocated to it, and 
			# returns a JobListManager object
		current_job_list = []
		with open(self.trackfile_path, 'rU') as trackfile_opened:
			trackfile_contents = list(csv.reader(trackfile_opened))
			for row in trackfile_contents[1:]:
				current_job = Job(*row)
				current_job_list.extend(current_job)
		return JobListManager(current_job_list,self.job_parameters,self.cluster_parameters)
	def write_output_files(self,job_list_manager):
		# writes trackfile and summaryfile for current object for each
			# job in job_list_manager, which is a JobListManager object
		self._write_trackfile(job_list_manager)
		self._write_summaryfile(job_list_manager)
	def _get_complete_status_list():
		# gets entire list of possible statuses from JobStatus class
		status_list = [getattr(JobStatus, attr) for attr in vars(JobStatus) \
			if not attr.startswith("__")]
		return(status_list)
	def _write_trackfile(self, job_list_manager):
		# writes a csv containing the status of each job in job_list_manager,
			# as well as the time and memory allocated to it
		complete_status_list = self._get_complete_status_list()
		job_list_contents = job_list_manager.extract_contents(complete_status_list)
		with open(self.trackfile_path, 'wb') as trackfile_opened:
			trackfile_writer = csv.writer(trackfile_opened)
			trackfile_writer.writerow(['Job Number','Job Status','Time Allocated to Job','Memory Allocated to Job'])
			for current_job in job_list_contents:
				trackfile_writer.writerow(current_job)
	def _write_summaryfile(self, job_list_manager):
		# writes a csv file counting the number of jobs of each status in job_list_manager
		with open(self.summaryfile_path, 'wb') as summaryfile_opened:
			summaryfile_writer = csv.writer(summaryfile_opened)
			summaryfile_writer.writerow(['Status','# of jobs','job indices'])
			complete_status_list = self._get_complete_status_list()
			for status in complete_status_list:
				indices = job_list_manager.get_jobs_by_status(status)
				current_n = len(indices)
				indices_concat = ';'.join(str(x) for x in indices)
				summaryfile_writer.writerow([status,current_n,indices_concat])

class SlurmManager(object):
	# Handles getting information from and passing information to the
		# slurm cluster system
	def __init__(self,job_parameters,cluster_parameters):
		self.job_parameters = copy.deepcopy(job_parameters)
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
		# at the request of hpc staff, don't use all available queue space
		self.max_job_proportion = 0.95
		self.sbatch_filename = os.path.join(self.job_parameters.cluster_job_submission_folder,\
			(self.job_parameters.name + '.q'))
		# add necessary lines for running multiple parallel matlab jobs
		if self.job_parameters.module == 'matlab' and \
			self.job_parameters.parallel_processors > 1:
			matlab_parallel_start_lines = ['if [ \"$SLURM_JOBTMP" == \"\" ]; then',\
				'    export SLURM_JOBTMP=/state/partition1/$USER/$$',\
				'    mkdir -p $SLURM_JOBTMP',\
				'fi',\
				'export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)']
			matlab_parallel_end_lines = ['rm -rf $SLURM_JOBTMP/*']
			self.job_parameters.additional_beginning_lines_in_sbatch.extend(matlab_parallel_start_lines)
			self.job_parameters.additional_end_lines_in_sbatch.extend(matlab_parallel_end_lines)
	def free_job_calculator(self):
		# gets the number of jobs that can still be submitted to
			# the routing queue
		# Get max number of jobs in a single array submission
		max_array_size_response_string = subprocess.check_output(
			'scontrol show config | grep MaxArraySize',shell=True)
		max_array_size_response_split = max_array_size_response_string.split(' = ')
		max_array_size = int(max_array_size_response_split[1].rstrip())
		# Get max number of jobs user can have in queue or running
		max_submit_response_string = subprocess.check_output(
			('sacctmgr list assoc format=user,maxsubmit where user='+self.cluster_parameters.username),
			shell=True)
		max_submit = int(re.findall((self.cluster_parameters.username+'\s+(\d+)'),max_submit_response_string)[0])
		# how many jobs are currently in default_queue for this user?
		jobs_in_queue = int(subprocess.check_output(
			('squeue -u '+self.cluster_parameters.username+' -r | egrep " PD | R | CG  " | wc -l'),
			shell=True))
		# find the max number of jobs you can run at one time
		max_allowed_jobs = round(min(max_array_size,max_submit)*self.max_job_proportion)
		# find max amount of jobs that can be added to queue without
			# making it overflow
		space_in_queue = max_allowed_jobs-jobs_in_queue
		return(space_in_queue)
	def _convert_time_format(time_in_mins):
		# convert time_in_mins from minutes to hh:mm:ss	
		hours = int(math.floor(float(time_in_mins)/60))
		if hours > 0:
			minutes = int(float(time_in_mins) % (float(hours)*60))
		else:
			minutes = int(float(time_in_mins))
		time_string = str(hours).zfill(2)+':'+ str(minutes).zfill(2)+':00'
		return(time_string)
	def _convert_mem_format(mem_in_mb):
		# convert single_job_mem from # of Mb to memory string
		# cluster only cares about units of GB
		mem_in_gb = math.ceil(float(mem_in_mb)/1024)
		mem_string = str(int(mem_in_gb))+'GB'
		return(mem_string)
	def _create_slurm_job(self,job_list_string,job_time,job_mem):
		# Writes files to submit to cluster queue
		# convert memory and time formats
		single_job_time_string = self._convert_time_format(job_time)
		single_job_mem_string = self._convert_mem_format(job_mem)
		# write submission file
		with open(self.sbatch_filename,'w') \
			as sbatch_job_file:
			sbatch_job_file.write('#!/bin/bash\n')
			sbatch_job_file.write('#SBATCH --job-name=' + \
				self.job_parameters.name + '\n')
				# name of current job
			sbatch_job_file.write('#SBATCH --output=' + \
				self.job_parameters.name + '.o%A-%a\n')
			sbatch_job_file.write('#SBATCH --error=' + \
				self.job_parameters.name + '.e%A-%a\n')
			sbatch_job_file.write('#SBATCH --time=' + \
				single_job_time_string + '\n')
				# amount of time allowed per job
			sbatch_job_file.write('#SBATCH --mem=' + \
				single_job_mem_string + '\n')
				# memory allocated to each job
			sbatch_job_file.write('#SBATCH --nodes=1\n')
			sbatch_job_file.write('#SBATCH --cpus-per-task=' + \
				str(self.job_parameters.parallel_processors) + '\n')
				# nodes and processors used per job
			sbatch_job_file.write('#SBATCH --array=' + job_list_string + '\n')
				# list of job IDs to be submitted
			sbatch_job_file.write('#SBATCH --mail-type=FAIL\n')
				# only sends email if job array fails
			sbatch_job_file.write('#SBATCH --mail-user=' + \
				self.cluster_parameters.user_email + '\n')
				# email that gets notification about aborted job
			# add any rows that need to be written for each particular file
			if self.job_parameters.additional_beginning_lines_in_sbatch:
				for additional_sbatch_beginning_row in \
					self.job_parameters.additional_beginning_lines_in_sbatch:
					sbatch_job_file.write(additional_sbatch_beginning_row + '\n')

			sbatch_job_file.write('cd ' + self.cluster_parameters.code_path + '\n')
				# cd into code directory
			sbatch_job_file.write('module purge\n')
			sbatch_job_file.write('module load ' + \
				self.job_parameters.module_path + '\n')
				# load module
			# write appropriate code-running line
			if self.job_parameters.module.lower() == 'matlab':
				sbatch_job_file.write('matlab -nodisplay -r '+
					code_run_string+'\n')
			elif self.job_parameters.module.lower() == 'r':
				sbatch_job_file.write('R CMD BATCH --vanilla ' + \
					code_run_string + '\n')
			# add any rows that need to be written at the end of each
				# particular file
			if self.job_parameters.additional_end_lines_in_sbatch:
				for additional_sbatch_end_row in \
					self.job_parameters.additional_end_lines_in_sbatch:
					sbatch_job_file.write(additional_sbatch_end_row + '\n')
			sbatch_job_file.write('\n\n')
				# need additional returns at end of shell scripts
	def _submit_slurm_job(self):
		# submits sbatch job
		# cd into sbatch directory
		os.chdir(self.job_parameters.cluster_job_submission_folder)			
		# Run .q file for sim
		subprocess.call('sbatch ' + \
			self.sbatch_filename,shell=True)
	def create_and_submit_batch_job(self,job_list_string,job_time,job_mem):
		# creates and submits batch job on slurm
		self._create_slurm_job(job_list_string,job_time,job_mem)
		self._submit_slurm_job()



#######################################################

def _create_job_list(job_name, job_numbers, initial_time, initial_mem, \
		job_parameters, cluster_parameters):
	# create a list of jobs sharing name, and number from 'numbers' list
	# all new jobs will have 'TO_PROCESS' status
	current_job_list = []
	for n in job_numbers:
		if type(n) is not int:
			raise TypeError("numbers contains non-integers: " + l)
		current_job = Job(n,JobStatus.TO_PROCESS,initial_time,initial_mem)
		current_job_list.extend(current_job)
	return JobListManager(current_job_list,job_parameters,cluster_parameters)

def job_flow_handler(job_name, job_numbers, initial_time, initial_mem, \
	cluster_parameters, output_folder, output_extension, output_filename, \
	cluster_job_submission_folder, experiment_folder, module, code_run_string, \
	additional_beginning_lines_in_sbatch, additional_end_lines_in_sbatch, \
	completefile_path):
	# Handles entire flow of job, from creation of new trackfile to
		# submission of jobs to updating trackfile
	completefile_manager = CompletefileManager(completefile_path)
	# check completeness and only run other code if job is not complete
	jobs_complete = completefile_manager.check_completeness()
	if not jobs_complete:
		########### ADD COMPLETENESS UPDATE BELOW
		job_parameters = JobParameters(job_name, output_folder, \
			output_extension, output_filename, cluster_job_submission_folder, experiment_folder, \
			module, code_run_string, additional_beginning_lines_in_sbatch, \
			additional_end_lines_in_sbatch, parallel_processors)
		current_trackfile = TrackfileManager(job_parameters, cluster_parameters)
		# check whether trackfile exists; if so, use it to get job data and
			# update statuses; otherwise, create new trackfile
		current_trackfile_path = current_trackfile.get_trackfile_path()
		if os.path.isfile(current_trackfile_path):
			# retrieve job list, update it, submit jobs
			# job_list_manager is a JobListManager object
			job_list_manager = current_trackfile.read_trackfile()
			# update status of all jobs in job_list_manager
			job_list_manager.autoupdate_job_status()
		else:
			# create job list (JobListManager object)
			job_list_manager = _create_job_list(job_name, job_numbers, initial_time, \
				initial_mem, job_parameters, cluster_parameters)
		# create cluster_submission_manager
		if cluster_parameters.cluster_architecture.lower() == 'slurm':
			cluster_submission_manager = SlurmManager(job_parameters,cluster_parameters)
		# submit whatever jobs possible
		job_list_manager.submit_jobs(cluster_submission_manager)
		# update trackfile and summaryfile
		current_trackfile.write_output_files(job_list_manager)
		# update completefile
		completefile_manager.add_job_list_manager(job_list_manager)
		completefile_manager.check_completeness()
