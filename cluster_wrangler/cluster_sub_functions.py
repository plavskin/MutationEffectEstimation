#!/usr/bin/python

# Contains submission managers for cluster

import copy
import os
import re
import subprocess
import math

class JobSubmissionManager(object):
	"""
	Base class that handles getting information from and passing
	information to the cluster system
	"""
	def __init__(self, cluster_parameters):
		self.cluster_parameters = copy.deepcopy(cluster_parameters)
	def get_within_batch_counter(self):
		return(self.within_batch_counter)
	def set_job_parameters(self,job_parameters):
		self.job_parameters = copy.deepcopy(job_parameters)
	def get_max_jobs_per_batch(self):
		return(None)
	def free_job_calculator(self):
		"""
		Gets the number of jobs that can still be submitted to
		the routing queue
		"""
		print('error! cannot run free_job_calculator on base JobSubmissionManager class.')
		space_in_queue = 0
		return(space_in_queue)
	def job_process_finder(self):
		"""
		Identifies jobs that are still being run
		Returns list of jobs currently being processed
		"""
		print('error! cannot run job_process_finder on base JobSubmissionManager class.')
		jobs_still_in_processing = []
		return(jobs_still_in_processing)
	def _create_submission_job(self, job_list_string, job_time, job_mem):
		""" Writes files to submit to cluster queue """
		# convert memory and time formats
		print('error! cannot create submission job in base JobSubmissionManager class.')
	def _submit_job(self):
		""" Submits job/job batch """
		# cd into sbatch directory
		# Run submission file
		print('error! cannot submit job from base JobSubmissionManager class.')
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
		### CURRENTLY CAN'T HANDLE PARALLEL MATLAB JOBS ###
		self.within_batch_counter = 'ARRAY_TASK_ID'
	def set_job_parameters(self,job_parameters):
		self.job_parameters = copy.deepcopy(job_parameters)
		self.sh_filename_prefix = os.path.join(self.job_parameters.cluster_job_submission_folder,\
			(self.job_parameters.name + '_'))
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
		jobs_running = int(subprocess.check_output(
			('ps aux | grep ' + self.job_parameters.module + \
				' | grep -v "grep" | wc -l'),shell=True))
		# find the max number of jobs you can run at one time
		max_allowed_jobs = number_cpus - self.free_processors
		# find max amount of jobs that can be added to queue without
			# making it overflow
		space_in_queue = max_allowed_jobs-jobs_running
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
	def _create_submission_job(self, job_number_string, *unused):
		""" Writes files to submit to cluster queue """
		# take first int in job_number_string as the required job number
		job_number = str(re.findall('\d+',job_number_string)[0])
		self.sh_filename = self.sh_filename_prefix + job_number + \
			self.sh_filename_suffix
		# write submission file
		with open(self.sh_filename,'w') \
			as batch_job_file:
			batch_job_file.write('#!/bin/bash\n')
			batch_job_file.write(\
				self.within_batch_counter + '=' + \
				job_number + '\n')
			batch_job_file.write('output_file=\'' + \
				os.path.join(self.job_parameters.cluster_job_submission_folder, \
				self.job_parameters.name) + '.o1-\'${' + \
				self.within_batch_counter + '}\n')
			batch_job_file.write('error_file=\'' + \
				os.path.join(self.job_parameters.cluster_job_submission_folder, \
				self.job_parameters.name) + \
				'.e1-\'${' + self.within_batch_counter + '}' + \
				'\n')
			# add any rows that need to be written for each particular file
			if self.job_parameters.additional_beginning_lines_in_job_sub:
				for additional_sbatch_beginning_row in \
					self.job_parameters.additional_beginning_lines_in_job_sub:
					batch_job_file.write(additional_sbatch_beginning_row + '\n')
			batch_job_file.write('cd \'' + self.cluster_parameters.code_path + '\'\n')
				# cd into code directory
			# write appropriate code-running line
			code_run_input = self.job_parameters.code_run_input
			code_run_input.set_full_code_run_string('macosx')
			code_run_string = code_run_input.get_code_run_string()
			batch_job_file.write(code_run_string  + '> ${output_file}\n')
			# move  any error message in output_file into error file
			batch_job_file.write('sed -n -e "/[Ee][Rr][Rr][Oo][Rr]/,\$w ${error_file}" "${output_file}"')
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
		# Run .q file for sim
		subprocess.call('sh \'' + self.sh_filename + '\'', shell=True)

class SlurmManager(JobSubmissionManager):
	"""
	Handles getting information from and passing information to the
	slurm cluster system
	"""
	def __init__(self, cluster_parameters):
		super(SlurmManager, self).__init__(cluster_parameters)
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
		# latest module is the last one in the list
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

			sbatch_job_file.write('cd ' + self.cluster_parameters.code_path + '\n')
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
		# Run .q file for sim
		subprocess.call('sbatch ' + \
			self.sbatch_filename,shell=True)
