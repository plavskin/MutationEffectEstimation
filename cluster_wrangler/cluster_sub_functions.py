#!/usr/bin/python

# Contains submission managers for cluster

import copy
import os
import re
import subprocess
import math
import numpy as np
import time
import warnings as war
war.filterwarnings("ignore", message="numpy.dtype size changed")

class BatchSubmissionManager(object):
	"""
	Holds parameters for current job submission session, and manages
	which jobs from current batch are submitted
	Algorithm here assumes space_in_queue is typically more limiting
	than max_char_num and batches_remaining
	"""
	def __init__(self, max_char_num, max_jobs_per_batch):
		self.jobs_to_submit = []
		self.job_submission_string_list = ['']
		self.max_char_num = max_char_num
		self.remaining_chars = max_char_num
		self.max_jobs_per_batch = max_jobs_per_batch
		self.remaining_jobs_in_batch = max_jobs_per_batch
		self.jobs_to_error = []
	def get_error_jobs(self):
		"""
		Returns jobs whose corresponding strings (for the job alone
		plus a comma) are longer than max_char_num
		"""
		return(self.jobs_to_error)
	def get_submission_string_list(self):
		return(self.job_submission_string_list)
	def get_submitted_jobs(self):
		return(self.jobs_to_submit)
	def get_batches_remaining(self):
		return(self.batches_remaining)
	def get_space_left_in_queue(self):
		return(self.space_in_queue)
	def _update_space_in_queue(self, job_number):
		self.space_in_queue = self.space_in_queue - job_number
	def _update_batches_remaining(self, batch_number):
		self.batches_remaining = self.batches_remaining - batch_number
	def _update_remaining_chars(self, char_number):
		self.remaining_chars = self.remaining_chars - char_number
	def _update_remaining_jobs_in_batch(self, job_number):
		self.remaining_jobs_in_batch = self.remaining_jobs_in_batch - job_number
	def _add_to_submission_string(self, new_string):
		self.job_submission_string_list[-1] = self.job_submission_string_list[-1]+new_string
	def _add_jobs_to_submission_list(self, job_sublist, job_string):
		self._add_to_submission_string(job_string)
		self.jobs_to_submit.extend(job_sublist)
		self._update_space_in_queue(len(job_sublist))
		self._update_remaining_chars(len(job_string))
		self._update_remaining_jobs_in_batch(len(job_sublist))
	def _reset_batch(self):
		self._update_batches_remaining(1)
		self.remaining_chars = self.max_char_num
	#	if not self.job_submission_string_list[-1] == '':
	#		self.job_submission_string_list.append('')
		self.job_submission_string_list.append('')
		self.remaining_jobs_in_batch = self.max_jobs_per_batch
	def _job_list_string_converter(self, job_list):
		"""
		Converts a list of job numbers into a job submission string
		"""
		raise AttributeError(\
			'BatchSubmissionManager._job_list_string_converter is a ' + \
			'placeholder, please use a child class of BatchSubmissionManager' + \
			' that properly parses job_list intro a job submission string')
		return(job_string)
	def _update_job_sublist(self, job_sublist, current_job_sublist, current_job_string):
		"""
		Adds current_job_sublist to job_submission_list,
		removes current_job_sublist from job_list
		"""
		self._add_jobs_to_submission_list(current_job_sublist,current_job_string)
		job_sublist = [x for x in job_sublist if x not in current_job_sublist]
		return(job_sublist)
	def _job_list_parser(self, job_sublist):
		"""
		Taking in a list of jobs, reorganizes them so that consecutive
		jobs can be called as intervals (e.g. '5-8' for '5,6,7,8'), and
		finds the optimal arrangement that maximizes job number so that
		the character count for the list of jobs doesn't exceed
		max_char_num
		Assumes that max_char_num is sufficiently large
		i.e. >(2*(max number of digits in a job number)+2)
		"""
		if self.space_in_queue < len(job_sublist):
			# reorder job_sublist and truncate
			job_sublist = np.sort(job_sublist)[0:self.space_in_queue]
		while len(job_sublist) > 0 and self.batches_remaining > 0:
			if self.remaining_jobs_in_batch > 0:
				# move number of jobs that can be from job_sublist to current_job_sublist
				num_jobs_to_run = int(min(len(job_sublist),self.remaining_jobs_in_batch))
				current_job_sublist = job_sublist[0:num_jobs_to_run]
				# create a string that can be used to submit those jobs
					# to SLURM, and count the number of characters in that string
				current_job_string = self._job_list_string_converter(current_job_sublist)
				if len(current_job_string) <= self.remaining_chars:
					job_sublist = self._update_job_sublist(job_sublist, \
						current_job_sublist, current_job_string)
				else:
					# try splitting up job_list one-by-one
					for current_job in current_job_sublist:
						current_job_string = self._job_list_string_converter([current_job])
						# for single jobs, first check whether string
							# is longer than max chars; if it is, add
							# it to errors; otherwise, if length is
							# less than remaining chars, add it to
							# current batch; otherwise, start new batch
						if len(current_job_string) > self.max_char_num:
							self.jobs_to_error.append(current_job)
							job_sublist.remove(current_job)
						elif len(current_job_string) <= self.remaining_chars:
							job_sublist = self._update_job_sublist(job_sublist, \
								[current_job], current_job_string)
						else:
							self._reset_batch()
							job_sublist = self._update_job_sublist(job_sublist, \
								[current_job], current_job_string)
			else:
				self._reset_batch()
	def select_jobs_to_sub(self, job_list, space_in_queue, batches_remaining):
		"""
		Based on the amount of space available in the queue, splits
		jobs into ones that can be run now vs those that have to be run
		later
		Takes the maximum number of jobs that can be submitted from
		initial_job_list, and parse them to ensure that
			1. they're listed in the most efficient format
			2. this list doesn't exceed the max number of chars
			allowed
		"""
		self.space_in_queue = space_in_queue
		self.batches_remaining = batches_remaining
		self.job_list = job_list
		if self.space_in_queue > 0 and self.batches_remaining > 0:
			self._job_list_parser(self.job_list)
		self._update_batches_remaining(1)

class BatchSubmissionManagerSlurm(BatchSubmissionManager):
	"""
	Holds parameters for current job submission session, and manages
	which jobs from current batch are submitted, for the SLURM
	workload manager
	Algorithm here assumes space_in_queue is typically more limiting
	than max_char_num and batches_remaining
	"""
	def __init__(self, max_char_num, max_jobs_per_batch):
		super(BatchSubmissionManagerSlurm, \
			self).__init__(max_char_num, max_jobs_per_batch)	
	def _consecutive_parser(self, job_list, stepsize = 1):
		"""
		Splits list of integers into list of lists of consecutive nums
		Returns the list of consecutive integer lists, as well as a
		list of lists of the indices of each member of each consecutive
		integer lists in the original job_list
		"""
		np_job_list = np.array(job_list)
		sorted_indices = np.argsort(np_job_list)
		sorted_job_list = np_job_list[sorted_indices]
		# Find position of indices where the number is higher than
			# previous number + stepsize
		split_indices = np.where(np.diff(sorted_job_list) != stepsize)[0]+1
		# Split data at split_indices
		split_job_list = np.array(np.split(sorted_job_list,split_indices))
		return(split_job_list)
	def _job_list_string_converter(self, consecutive_job_list):
		"""
		Converts a consecutive list of job numbers into a job
		submission string for SLURM
		"""
		# test that consecutive_job_list is actually consecutive
		test_list = self._consecutive_parser(consecutive_job_list)
		if len(test_list) > 1:
			print('Error: BatchSubmissionManager._job_list_string_converter provided unparsed list')
		# Denote consecutive job sublists by dashes between the min
			# and max of the sublist; separate job sublists by commas
		current_job_num = len(consecutive_job_list)
		if current_job_num==1:
			current_job_string = (str(consecutive_job_list[0])+',')
		else:
			current_job_string = (str(min(consecutive_job_list))+'-'+
				str(max(consecutive_job_list))+',')
		return(current_job_string)
	def select_jobs_to_sub(self, job_list, space_in_queue, batches_remaining):
		"""
		Based on the amount of space available in the queue, splits
		jobs into ones that can be run now vs those that have to be run
		later
		Takes the maximum number of jobs that can be submitted from
		initial_job_list, and parse them to ensure that
			1. they're listed in the most efficient format
			2. this list doesn't exceed the max number of chars
			allowed
		"""
		self.space_in_queue = space_in_queue
		self.batches_remaining = batches_remaining
		self.job_list = job_list
		job_list_split = self._consecutive_parser(self.job_list)
		job_number_list = np.array([len(i) for i in job_list_split])
		order_by_job_num = np.argsort(job_number_list)
		job_list_split_sorted = job_list_split[order_by_job_num]
		for current_job_list in job_list_split_sorted:
			if self.space_in_queue > 0 and self.batches_remaining > 0:
				self._job_list_parser(current_job_list)
		self._update_batches_remaining(1)

class BatchSubmissionManagerMacOSX(BatchSubmissionManager):
	"""
	Holds parameters for current job submission session, and manages
	which jobs from current batch are submitted, for the SLURM
	workload manager
	Algorithm here assumes space_in_queue is typically more limiting
	than max_char_num and batches_remaining
	"""
	def __init__(self, max_char_num, max_jobs_per_batch):
		super(BatchSubmissionManagerMacOSX, \
			self).__init__(max_char_num, max_jobs_per_batch)
	def _job_list_string_converter(self, job_list):
		"""
		Converts a list of job numbers into a job submission string for
		MacOSX (a bash array)
		"""
		job_list_as_str = [str(x) for x in job_list]
		current_job_string = ' '.join(job_list_as_str) + ' '
		return(current_job_string)

class JobSubmissionManager(object):
	"""
	Base class that handles getting information from and passing
	information to the cluster system
	"""
	def __init__(self, cluster_parameters):
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
		self.submission_manager = \
			BatchSubmissionManager(cluster_parameters.max_char_num, \
				cluster_parameters.max_jobs_per_batch)
	def get_within_batch_counter(self):
		return(self.within_batch_counter)
	def set_job_parameters(self,job_parameters):
		self.job_parameters = copy.deepcopy(job_parameters)
	def get_max_jobs_per_batch(self):
		return(None)
	def get_submission_manager(self):
		return(self.submission_manager)
	def free_job_calculator(self):
		"""
		Gets the number of jobs that can still be submitted to
		the routing queue
		"""
		print('error: cannot run free_job_calculator on base JobSubmissionManager class.')
		space_in_queue = 0
		return(space_in_queue)
	def job_process_finder(self):
		"""
		Identifies jobs that are still being run
		Returns list of jobs currently being processed
		"""
		print('error: cannot run job_process_finder on base JobSubmissionManager class.')
		jobs_still_in_processing = []
		return(jobs_still_in_processing)
	def _create_submission_job(self, job_list_string, job_time, job_mem):
		""" Writes files to submit to cluster queue """
		# convert memory and time formats
		print('error: cannot create submission job in base JobSubmissionManager class.')
	def _submit_job(self):
		""" Submits job/job batch """
		# cd into sbatch directory
		# Run submission file
		print('error: cannot submit job from base JobSubmissionManager class.')
	def create_and_submit_batch_job(self, job_list_string, job_time, job_mem):
		""" Creates and submits batch/single job on computer/cluster """
		self._create_submission_job(job_list_string,job_time,job_mem)
		self._submit_job()
	def error_status_check(self, latest_errorfile_contents):
		""" Parses contents of job submission run error file """
		latest_errorfile_contents = latest_errorfile_contents.lower()
		error_status_dict = dict()
		error_status_dict['time_limit_check'] = \
			'due to time' in latest_errorfile_contents or \
			'job time limit' in latest_errorfile_contents
		error_status_dict['memory_limit_check'] = \
			'due to memory' in latest_errorfile_contents or \
			'job memory limit' in latest_errorfile_contents
		error_status_dict['cluster_error_check'] = \
			'bus error' in latest_errorfile_contents \
			or 'fatal error on startup' in latest_errorfile_contents \
			or 'reload' in latest_errorfile_contents \
			or 'matlabexception' in latest_errorfile_contents
		error_status_dict['unidentified_error_check'] = \
			len(latest_errorfile_contents) > self.empty_errorfile_size and \
			not any([error_status_dict['time_limit_check'], \
				error_status_dict['memory_limit_check'], error_status_dict['cluster_error_check']])
		return(error_status_dict)

class MacOSXManager(JobSubmissionManager):
	"""
	Handles getting information from and passing information to a
	MacOSX computer
	"""
	def __init__(self, cluster_parameters):
		super(MacOSXManager, self).__init__(cluster_parameters)
		self.submission_manager = \
			BatchSubmissionManagerMacOSX(cluster_parameters.max_char_num, \
				cluster_parameters.max_jobs_per_batch)
		# Don't use every processor on computer! (duh)
		self.free_processors = 1
		# On a personal computer, no slowing down as a penalty for
			# submitting jobs one-by-one so don't reset
			# max_sub_batches_in_one_run if submitting jobs one-by-one
		self.max_sub_batches_in_one_run = float('Inf')
		self.sh_filename_suffix = '.sh'
			# job_name needs to be part of this filename in order for job tracking to work
		# specify size of an empty errorfile on this cluster architecture
		self.empty_errorfile_size = 0
		self.within_batch_counter = 'ARRAY_TASK_ID'
	def set_job_parameters(self,job_parameters):
		self.job_parameters = copy.deepcopy(job_parameters)
		self.sh_filename_prefix = \
			os.path.join(self.job_parameters.cluster_job_submission_folder, \
			self.job_parameters.name)
	def free_job_calculator(self):
		"""
		Gets the number of jobs that can still be submitted to
		the routing queue
		"""
		# Get max number of processors user can use
		try:
			number_cpus = \
				int(subprocess.check_output('getconf _NPROCESSORS_ONLN',shell=True))
		except subprocess.CalledProcessError:
			number_cpus = \
				int(subprocess.check_output('getconf NPROCESSORS_ONLN',shell=True))
			# one of the above should work on MacOSX and linux machines
		# how many jobs are currently running on computer?
		# calculate this by assuming only jobs from module (e.g.
			# matlab) are relevant, i.e. count those
#		module_execution_path = \
#			subprocess.check_output('which ' + self.job_parameters.module, \
#				shell = True).rstrip()
		jobs_running = int(subprocess.check_output(
			('ps aux | grep -i ' + self.job_parameters.module + \
				' | grep -v "grep" | wc -l'),shell=True))
		# find the max number of jobs you can run at one time
		max_allowed_jobs = number_cpus - self.free_processors
		# find max amount of jobs that can be added to queue without
			# making it overflow
		space_in_queue = max_allowed_jobs - jobs_running
		return(space_in_queue)
	def job_process_finder(self):
		"""
		Identifies jobs that are still being run
		Returns list of jobs currently being processed
		"""
		try:
			jobs_running_list = subprocess.check_output(('ps aux | grep ' + \
				self.job_parameters.name + ' | grep -v "grep"'),shell=True)
		except subprocess.CalledProcessError:
			jobs_running_list = ''
		jobs_still_in_processing = \
			[int(k) for k in re.findall('\w+_(\d+)' + self.sh_filename_suffix, \
				jobs_running_list, re.MULTILINE)]
		return(jobs_still_in_processing)
	def _generate_filename_for_sub_job(self, prenumber_prefix):
		current_filename = \
			self.sh_filename_prefix + '.'+ str(prenumber_prefix) + \
			'1-\'${' + self.within_batch_counter + '}'
		return(current_filename)
	def _create_submission_job(self, job_number_string, *unused):
		""" Writes files to submit to cluster queue """
		# take first int in job_number_string as the required job number
		self.sh_filename = \
			self.sh_filename_prefix + self.sh_filename_suffix
		self.code_run_filename = \
			self.sh_filename_prefix + '_code_run_file' + self.sh_filename_suffix
		# generate filenames needed in scripts below
		output_filename = self._generate_filename_for_sub_job('o')
		error_filename = self._generate_filename_for_sub_job('e')
		screen_rc_filename = self._generate_filename_for_sub_job('rc')
		screen_out_filename = self._generate_filename_for_sub_job('screen-out')
		# get string that will submit actual code to run (e.g. MLE code)
		code_run_input = self.job_parameters.code_run_input
		code_run_input.set_full_code_run_string('macosx')
		code_run_string = code_run_input.get_code_run_string()
		# write code run file
		with open(self.code_run_filename,'w') \
			as code_run_file:
			code_run_file.write('#!/bin/bash\n')
			code_run_file.write(self.within_batch_counter + '=${1}\n')
			code_run_file.write('output_file=\'' + \
				output_filename + '\n')
			code_run_file.write('error_file=\'' + \
				error_filename + '\n')
			code_run_file.write('screen_out_file=\'' + \
				screen_out_filename + '\n')
			# cd into code directory
			code_run_file.write('cd \'' + self.job_parameters.code_path + '\'\n')
			# write appropriate code-running line
			code_run_file.write(code_run_string  + \
				' | tee "${output_file}"\n')
			code_run_file.write('sed -n -e "/[Ee][Rr][Rr][Oo][Rr]/,\$w ${error_file}" "${output_file}"\n')
			code_run_file.write('rm "${screen_out_file}"\n\n\n')
		# write submission file
		with open(self.sh_filename,'w') \
			as batch_job_file:
			batch_job_file.write('#!/bin/bash\n')
			batch_job_file.write('ARRAY_TASK_LIST=(' + job_number_string + ')\n')
			batch_job_file.write('for ' + self.within_batch_counter + \
				' in ${ARRAY_TASK_LIST[@]}\n')
			batch_job_file.write('do\n')
			batch_job_file.write('screen_rc_file=\'' + \
				screen_rc_filename + '\n')
			batch_job_file.write('code_run_file=\'' + \
				self.code_run_filename + '\'\n')
			batch_job_file.write('screen_out_file=\'' + \
				screen_out_filename + '\n')
			batch_job_file.write('cat << EOF >${screen_rc_file}\n')
			batch_job_file.write('logfile "${screen_out_file}"\n')
			batch_job_file.write('EOF\n')
			# add any rows that need to be written for each particular file
			if self.job_parameters.additional_beginning_lines_in_job_sub:
				for additional_sbatch_beginning_row in \
					self.job_parameters.additional_beginning_lines_in_job_sub:
					batch_job_file.write(additional_sbatch_beginning_row + '\n')
			# run code via detached screen, move any error message in output_file into error file
			batch_job_file.write('screen -d -m -c "${screen_rc_file}" -L sh ' + \
				'"${code_run_file}" ${' + self.within_batch_counter + '}\n')
			batch_job_file.write('rm "${screen_rc_file}"\n')
			batch_job_file.write('done\n')
			# add any rows that need to be written at the end of each
				# particular file
			if self.job_parameters.additional_end_lines_in_job_sub:
				for additional_sbatch_end_row in \
					self.job_parameters.additional_end_lines_in_job_sub:
					batch_job_file.write(additional_sbatch_end_row + '\n')
			batch_job_file.write('\n\n')
				# need additional returns at end of shell scripts
	def _submit_job(self):
		""" Submits sh job """
		# cd into sh directory
		os.chdir(self.job_parameters.cluster_job_submission_folder)
		subprocess.call('sh \'' + self.sh_filename + '\'', shell=True)
		# pause 6 seconds to allow submitted processes to start, so
			# that they will be detected when looking for any running
			# jobs
		time.sleep(6)

class SlurmManager(JobSubmissionManager):
	"""
	Handles getting information from and passing information to the
	slurm cluster system
	"""
	def __init__(self, cluster_parameters):
		super(SlurmManager, self).__init__(cluster_parameters)
		self.submission_manager = \
			BatchSubmissionManagerSlurm(cluster_parameters.max_char_num, \
				cluster_parameters.max_jobs_per_batch)
		# at the request of hpc staff, don't use all available queue space
		self.max_job_proportion = 0.95
		# specify size of an empty errorfile on this cluster architecture
		self.empty_errorfile_size = 0
		# SLURM-based cluster doesn't appear to slow down as a penalty for
			# submitting jobs one-by-one so don't reset
			# max_sub_batches_in_one_run if submitting jobs one-by-one
		self.max_sub_batches_in_one_run = float('Inf')
		self.within_batch_counter = 'SLURM_ARRAY_TASK_ID'
	def set_job_parameters(self,job_parameters):
		self.job_parameters = copy.deepcopy(job_parameters)
		self._get_latest_module_path()
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
			self.job_parameters.additional_beginning_lines_in_job_sub.extend(matlab_parallel_start_lines)
			self.job_parameters.additional_end_lines_in_job_sub.extend(matlab_parallel_end_lines)
	def get_max_jobs_per_batch(self):
		""" Gets max number of jobs in a single array submission """
		max_array_size_response_string = subprocess.check_output(
			'scontrol show config | grep MaxArraySize',shell=True)
		max_array_size_response_split = max_array_size_response_string.split(' = ')
		max_array_size = int(max_array_size_response_split[1].rstrip())
		# Get max number of jobs user can have in queue or running
		max_submit_response_string = subprocess.check_output(
			('sacctmgr list assoc format=user,maxsubmit where user='+self.cluster_parameters.username),
			shell=True)
		max_submit = int(re.findall((self.cluster_parameters.username+'\s+(\d+)'),max_submit_response_string)[0])
		# find the max number of jobs you can run at one time
		max_allowed_jobs = int(round(min(max_array_size,max_submit)*self.max_job_proportion))
		return(max_allowed_jobs)
	def free_job_calculator(self):
		"""
		Gets the number of jobs that can still be submitted to
		the routing queue
		"""
		max_allowed_jobs = self.get_max_jobs_per_batch()
		# how many jobs are currently in default_queue for this user?
		jobs_in_queue = int(subprocess.check_output(
			('squeue -u '+self.cluster_parameters.username+' -r | egrep " PD | R | CG  " | wc -l'),
			shell=True))
		# find max amount of jobs that can be added to queue without
			# making it overflow
		space_in_queue = max_allowed_jobs-jobs_in_queue
		return(space_in_queue)
	def job_process_finder(self):
		"""
		Identifies jobs that are still being run by SLURM
		Return list of jobs currently being processed
		"""
		try:
			jobs_running_list = subprocess.check_output('squeue -u ' + self.cluster_parameters.username + ' -r -n '
				+ self.job_parameters.name + ' | egrep " PD | CG | R " ',shell=True)
		except subprocess.CalledProcessError:
			jobs_running_list = ''
		jobs_still_in_processing = [int(k) for k in re.findall('^\s+\d+_(\d+)\s',jobs_running_list,re.MULTILINE)]
		return(jobs_still_in_processing)
	def _get_latest_module_path(self):
		"""
		Identifes the most up-to-date path for the module of the
		application you need to use
		"""
		try:
			module_output = subprocess.check_output('module avail', shell=True)
		except subprocess.CalledProcessError:
			module_output = ''
		relevant_module_list = re.findall(' (' + self.job_parameters.module + \
			os.sep + '.+?)\\n', module_output)
		# latest module is the last one in the sorted list
		relevant_module_list.sort()
		latest_module_path = re.findall('(' + self.job_parameters.module + \
			os.sep + '\S+)', relevant_module_list[-1])[0]
		self.module_path = latest_module_path
	def _convert_time_format(self, time_in_mins):
		""" Converts time_in_mins from minutes to hh:mm:ss """	
		hours = int(math.floor(float(time_in_mins)/60))
		if hours > 0:
			minutes = int(float(time_in_mins) % (float(hours)*60))
		else:
			minutes = int(float(time_in_mins))
		time_string = str(hours).zfill(2)+':'+ str(minutes).zfill(2)+':00'
		return(time_string)
	def _convert_mem_format(self, mem_in_mb):
		"""
		Converts single_job_mem from # of Mb to memory string
		Cluster only cares about units of GB
		"""
		mem_in_gb = math.ceil(float(mem_in_mb)/1024)
		mem_string = str(int(mem_in_gb))+'GB'
		return(mem_string)
	def _create_submission_job(self, job_list_string, job_time, job_mem):
		""" Writes files to submit to cluster queue """
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
			if self.job_parameters.additional_beginning_lines_in_job_sub:
				for additional_sbatch_beginning_row in \
					self.job_parameters.additional_beginning_lines_in_job_sub:
					sbatch_job_file.write(additional_sbatch_beginning_row + '\n')

			sbatch_job_file.write('cd ' + self.job_parameters.code_path + '\n')
				# cd into code directory
			sbatch_job_file.write('module purge\n')
			sbatch_job_file.write('module load ' + \
				self.module_path + '\n')
				# load module
			# write appropriate code-running line
			code_run_input = self.job_parameters.code_run_input
			code_run_input.set_full_code_run_string('slurm')
			code_run_string = code_run_input.get_code_run_string()
			sbatch_job_file.write(code_run_string + '\n')
			# add any rows that need to be written at the end of each
				# particular file
			if self.job_parameters.additional_end_lines_in_job_sub:
				for additional_sbatch_end_row in \
					self.job_parameters.additional_end_lines_in_job_sub:
					sbatch_job_file.write(additional_sbatch_end_row + '\n')
			sbatch_job_file.write('\n\n')
				# need additional returns at end of shell scripts
	def _submit_job(self):
		""" Submits sbatch job """
		# cd into sbatch directory
		os.chdir(self.job_parameters.cluster_job_submission_folder)			
		subprocess.call('sbatch ' + \
			self.sbatch_filename,shell=True)
