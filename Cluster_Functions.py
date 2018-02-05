#!/usr/bin/python

# Contains objects needed for running jobs on cluster and keeping track
	# of their progress

from enum import Enum

class Job(object):
	# job object that stores properties of individual jobs in queue
	def __init__(self, number, status):
		self.number = number
		self.status = status
	def change_status(self,new_status):
		self.status = new_status

class JobStatus(object):
	# enum for job status
	TO_PROCESS = 1
	PROCESSING = 2
	COMPLETED = 3
	ABORTED_FOREVER = 4
	ABORTED_TO_RESTART = 5

class JobList(object):
	# holds list of jobs corresponding to a single 'name'
	# updates current job status
	def __init__(self, jobs, name):
		self.jobs = {}
		for j in jobs:
			self.jobs[j.number] = j
		self.name = name
	def _get_status(self, number):
		# gets status of individual job
		# ??? Do I need this function? ???
		if number in self.jobs:
			return self.jobs[number].status
		else:
			raise ValueError("Invalid job : " + number)
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
				job_subset_list.append(num)

	def batch_status_change(self, number_list, new_status):
		# changes the status of al jobs in number_list to new_status
		for num in number_list:
			self.jobs[num].change_status(new_status)


def create_job_list(name,numbers):
	# create a list of jobs sharing name, and number from 'numbers' list
	# all new jobs will have 'TO_PROCESS' status
	current_job_list = []
	for n in numbers:
		if type(n) is not int:
			raise TypeError("numbers contains non-integers: " + l)
		current_job = Job(n,JobStatus.TO_PROCESS)
		current_job_list.append(current_job)
	return JobList(current_job_list,name)


test_list.jobs[1].change_status(JobStatus.COMPLETED)





def Free_Job_Calculator(username):
	# gets the number of jobs that can still be submitted to
		# the routing queue

	# Get max number of jobs in a single array submission
	max_array_size_response_string = subprocess.check_output(
		'scontrol show config | grep MaxArraySize',shell=True)
	max_array_size_response_split = max_array_size_response_string.split(' = ')
	max_array_size = int(max_array_size_response_split[1].rstrip())

	# Get max number of jobs user can have in queue or running
	max_submit_response_string = subprocess.check_output(
		('sacctmgr list assoc format=user,maxsubmit where user='+username),
		shell=True)
	max_submit = int(re.findall((username+'\s+(\d+)'),max_submit_response_string)[0])
	#max_submit = 1

	# how many jobs are currently in default_queue for this user?
	jobs_in_queue = int(subprocess.check_output(
		('squeue -u '+username+' -r | egrep " PD | R | CG  " | wc -l'),
		shell=True))

	# find the max number of jobs you can run at one time
	max_allowed_jobs = min(max_array_size,max_submit)

#	# at the request of hpc staff, don't use all available queue space
	max_allowed_jobs = round(max_allowed_jobs*.9)

	# find max amount of jobs that can be added to queue without
		# making it overflow
	space_in_queue = max_allowed_jobs-jobs_in_queue

	return(space_in_queue)


def Aborted_Job_Processor(jobs_aborted_to_restart,jobs_aborted_newmems,
	jobs_aborted_newtimes,jobs_aborted_forever,jobs_aborted_to_restart_new,
	jobs_aborted_newtimes_new,jobs_aborted_newmems_new,max_memory,max_time,
	time_multiplier,memory_multiplier,current_missing_job,default_memory,
	default_time):
	
	# find indices of current_missing_job in jobs_aborted_to_restart
	indices_in_old_restart_list = [i for i, x in \
		enumerate(jobs_aborted_to_restart) if x == current_missing_job]
	if indices_in_old_restart_list:
		# use info from most recent job abort
		most_recent_restart_idx = max(indices_in_old_restart_list)
		most_recent_time = float(jobs_aborted_newtimes[most_recent_restart_idx])
		most_recent_mem = float(jobs_aborted_newmems[most_recent_restart_idx])
		del jobs_aborted_to_restart[most_recent_restart_idx]
		del jobs_aborted_newtimes[most_recent_restart_idx]
		del jobs_aborted_newmems[most_recent_restart_idx]
	else:
		most_recent_time = default_time
		most_recent_mem = default_memory

	current_updated_memory = min(most_recent_mem*memory_multiplier,max_memory)
	current_updated_time = min(most_recent_time*time_multiplier,max_time)

	if most_recent_time == max_time or most_recent_mem == max_memory:
		jobs_aborted_forever.extend([current_missing_job])
	else:
		jobs_aborted_to_restart_new.extend([current_missing_job])
		jobs_aborted_newtimes_new.extend([current_updated_time])
		jobs_aborted_newmems_new.extend([current_updated_memory])

	return(jobs_aborted_forever,jobs_aborted_to_restart_new,
		jobs_aborted_newtimes_new,jobs_aborted_newmems_new,
		jobs_aborted_to_restart,jobs_aborted_newtimes,jobs_aborted_newmems)

def Job_Process_Finder(username,current_job_name):
	# Identify jobs that are still being run by SLURM

	# count jobs currently being processed
	try:
		qstat_output = subprocess.check_output('squeue -u ' + username + ' -r -n '
			+ current_job_name + ' | egrep " PD | CG | R " ',shell=True)
	except subprocess.CalledProcessError:
		qstat_output = ''

	jobs_still_in_processing = re.findall('^\s+\d+_(\d+)\s',qstat_output,re.MULTILINE)

	return(jobs_still_in_processing)

def Trackfile_Processor(max_repeat,current_mode,current_trackfile,current_completefile,username,current_job_name,
	slurm_folder,current_output_folder,current_output_filename,current_output_extension,output_dir_name,
	max_memory,max_time,default_memory,default_time):
	# Creates and updates trackfiles, which keep track of jobs to be completed
	# writes separate trackfile for every mode, growth condition, rep, parameter

	# Check if trackfile exists, and if not, create it
	if not os.path.isfile(current_trackfile):
		# generate a string of numbers from 1 to max_repeat
		all_jobs_to_process = ';'.join(str(x) for x in range(1,max_repeat+1))

		# write trackfile:
			# rows 1-4: 'knownmuts','poissonmut','lambdafit','SSR_mixed';
			# cols 1-4: 'to process', 'processing', 'completed'
		with open(current_trackfile, 'wb') as csvfile:
			trackfile_writer = csv.writer(csvfile)
			trackfile_writer.writerow(['mode','to process','processing','completed','errors','aborted forever','aborted to restart','aborted restart times','aborted restart memory'])
			trackfile_writer.writerow([current_mode,all_jobs_to_process,'','','','','','',''])
			trackfile_writer.writerow(['# jobs','','','','','','','',''])

	# Read trackfile for current mode
	with open(current_trackfile, 'rU') as trackfile_opened:
		trackfile_contents = list(csv.reader(trackfile_opened))

	jobs_to_process_string = trackfile_contents[1][1]
	jobs_processing_string = trackfile_contents[1][2]
	jobs_completed_string = trackfile_contents[1][3]
	jobs_with_errors_string = trackfile_contents[1][4]
	jobs_aborted_forever_string = trackfile_contents[1][5]
	jobs_aborted_to_restart_string = trackfile_contents[1][6]
	jobs_aborted_newtimes_string = trackfile_contents[1][7]
	jobs_aborted_newmems_string = trackfile_contents[1][8]

	jobs_to_process = filter(None,jobs_to_process_string.split(';'))
	jobs_processing = filter(None,jobs_processing_string.split(';'))
	jobs_completed = filter(None,jobs_completed_string.split(';'))
	jobs_with_errors = filter(None,jobs_with_errors_string.split(';'))
	jobs_aborted_forever = filter(None,jobs_aborted_forever_string.split(';'))
	jobs_aborted_to_restart = filter(None,jobs_aborted_to_restart_string.split(';'))
	jobs_aborted_newtimes = filter(None,jobs_aborted_newtimes_string.split(';'))
	jobs_aborted_newmems = filter(None,jobs_aborted_newmems_string.split(';'))
		# 'filter' removes empty strings

	jobs_to_process_updated = jobs_to_process[:]
	jobs_processing_updated = []
	jobs_completed_updated = jobs_completed[:]
	jobs_with_errors_updated = jobs_with_errors[:]
	jobs_aborted_forever_updated = jobs_aborted_forever[:]
	jobs_aborted_to_restart_updated = jobs_aborted_to_restart[:]
	jobs_aborted_newtimes_updated = jobs_aborted_newtimes[:]
	jobs_aborted_newmems_updated = jobs_aborted_newmems[:]
	# make sure there are any jobs being processed right now

	jobs_aborted_to_restart_new = []
	jobs_aborted_newtimes_new = []
	jobs_aborted_newmems_new = []

	jobs_still_in_processing = Job_Process_Finder(username,current_job_name)
	jobs_processing_updated.extend(jobs_still_in_processing)

	# move jobs that are not still in processing into 'jobs_completed_updated'
	jobs_just_finished = list(set(jobs_processing)-set(jobs_still_in_processing))

	# In jobs that have just been finished, look for jobs that never
		# completed, and sort those into jobs that just need to be
		# restarted and error jobs

	if jobs_just_finished:
		
		try:
			completed_files = subprocess.check_output('ls -lrt '
				+ current_output_folder,shell=True)
		except subprocess.CalledProcessError:
			completed_files = ''

		completed_job_list = re.findall(' ' + current_output_filename + '_(\d+?)\.'
			+ current_output_extension,completed_files,re.DOTALL)

		true_jobs_just_finished = list(set(completed_job_list) & set(jobs_just_finished))
		jobs_completed_updated.extend(true_jobs_just_finished)

		missing_jobs = list(set(jobs_just_finished)-set(true_jobs_just_finished))
		if missing_jobs:

			# if error file exists:
				#	1. check whether it contains info on exceeded time or mem limits;
				#		if so, resubmit with max time and memory
				#	2. Otherwise, check whether error file empty
				#		if not, add field to error list
				#		if so, restart field

			try:
				sbatch_run_output_list = subprocess.check_output('ls -lrt ' + slurm_folder,shell=True)
			except subprocess.CalledProcessError:
				sbatch_run_output_list = ''
			
			for current_missing_job in missing_jobs:

				current_missing_job_codes = re.findall(
					current_job_name + '.e(\d+?)-' + current_missing_job + '$',
					sbatch_run_output_list,re.MULTILINE)
				if current_missing_job_codes:
					current_job_code = str(max(map(int,current_missing_job_codes)))
						# jobs are assigned consecutive code numbers on cluster,
							# so since we want to look at the error file associated
							# with the most recent run of current_missing_job,
							# look for max code
					current_errorfile = slurm_folder + '/' + current_job_name + '.e' + current_job_code \
						+ '-' + current_missing_job

					# If errorfile is missing or empty; or if a cluster
						# error has occurred; but job is listed as
						# having started, return job to list of jobs
						# to process
					# If errorfile contains a time limit error or a
						# memory limit error, process as an aborted
						# job
					# If errorfile contains something else, report
						# as an error
					if os.path.isfile(current_errorfile):
						current_errorfile_contents = open(current_errorfile).read()
					else:
						current_errorfile_contents = ''

					if len(current_errorfile_contents) == 0:
						jobs_to_process_updated.extend([current_missing_job])
					else:
						time_limit_check = 'time limit' in current_errorfile_contents.lower()
						memory_limit_check = 'memory limit' in current_errorfile_contents.lower()
						cluster_error_check = 'bus error' in current_errorfile_contents.lower() \
							or 'fatal error on startup' in current_errorfile_contents.lower() \
							or 'reload' in current_errorfile_contents.lower() \
							or 'MatlabException' in current_errorfile_contents.lower()

						if cluster_error_check:
							jobs_to_process_updated.extend([current_missing_job])
						elif memory_limit_check:
							# check whether job had been restarted
								# before; if so, updated memory should
								# be 1.5x times previous memory allotment
								# otherwise, it should be 1.5x the default
								# mem allotment
							# If previous mem allotment was the max
								# allowed memory, abort job forever
							time_multiplier = 1
							memory_multiplier = 1.5

							(jobs_aborted_forever_updated,jobs_aborted_to_restart_new,
								jobs_aborted_newtimes_new,jobs_aborted_newmems_new,
								jobs_aborted_to_restart_updated,
								jobs_aborted_newtimes_updated,
								jobs_aborted_newmems_updated) = \
								Aborted_Job_Processor(jobs_aborted_to_restart_updated,
									jobs_aborted_newmems_updated,jobs_aborted_newtimes_updated,
									jobs_aborted_forever_updated,jobs_aborted_to_restart_new,
									jobs_aborted_newtimes_new,jobs_aborted_newmems_new,
									max_memory,max_time,time_multiplier,
									memory_multiplier,current_missing_job,
									default_memory,default_time)

						elif time_limit_check:
							# check whether job had been restarted
								# before; if so, updated time should
								# be 2x times previous time allotment;
								# otherwise, it should be 2x the default
								# time allotment
							# If previous mem allotment was the max
								# allowed memory, abort job forever
							time_multiplier = 2
							memory_multiplier = 1

							(jobs_aborted_forever_updated,jobs_aborted_to_restart_new,
								jobs_aborted_newtimes_new,jobs_aborted_newmems_new,
								jobs_aborted_to_restart_updated,
								jobs_aborted_newtimes_updated,
								jobs_aborted_newmems_updated) = \
								Aborted_Job_Processor(jobs_aborted_to_restart_updated,
									jobs_aborted_newmems_updated,jobs_aborted_newtimes_updated,
									jobs_aborted_forever_updated,jobs_aborted_to_restart_new,
									jobs_aborted_newtimes_new,jobs_aborted_newmems_new,
									max_memory,max_time,time_multiplier,
									memory_multiplier,current_missing_job,
									default_memory,default_time)
						else:
							jobs_with_errors_updated.extend([current_missing_job])
				else:
					jobs_to_process_updated.extend([current_missing_job])

			jobs_aborted_to_restart_updated.extend(jobs_aborted_to_restart_new)
			jobs_aborted_newtimes_updated.extend(jobs_aborted_newtimes_new)
			jobs_aborted_newmems_updated.extend(jobs_aborted_newmems_new)

	aborted_jobs_to_submit = []
	aborted_newtimes_to_submit = []
	aborted_newmems_to_submit = []
	new_jobs_to_submit = []

	# If there's space in the queue and jobs_to_process isn't empty, submit jobs
	space_in_queue = Free_Job_Calculator(username)

	############### ??? TO DO ??? ###############
	# Need to group aborted jobs by combo of used memory and time, then loop through groups in Jobs_to_Run_Selector
	############### ??? TO DO ??? ###############

	# run aborted jobs through Jobs_to_Run_Selector one-by-one to avoid list truncation due to 
	for counter in range(0,len(jobs_aborted_to_restart_updated)):
		# check that aborted job is not currently processing
		if not jobs_aborted_to_restart_updated[counter] in \
			(jobs_still_in_processing + jobs_with_errors_updated + jobs_completed_updated):

			current_jobs_aborted_to_restart_updated = [jobs_aborted_to_restart_updated[counter]]
			current_jobs_aborted_newtimes_updated = [jobs_aborted_newtimes_updated[counter]]
			current_jobs_aborted_newmems_updated = [jobs_aborted_newmems_updated[counter]]

			[[current_aborted_jobs_to_submit,
					current_aborted_newtimes_to_submit,
					current_aborted_newmems_to_submit],\
				[current_jobs_aborted_to_restart_updated,
					current_jobs_aborted_newtimes_updated,
					current_jobs_aborted_newmems_updated]] = \
				Jobs_to_Run_Selector(username,
					[current_jobs_aborted_to_restart_updated,
						current_jobs_aborted_newtimes_updated,
						current_jobs_aborted_newmems_updated],
					space_in_queue)

				# should return empty output list if input lists are empty
			space_in_queue = space_in_queue - 1
			aborted_jobs_to_submit.extend(current_aborted_jobs_to_submit)
			aborted_newtimes_to_submit.extend(current_aborted_newtimes_to_submit)
			aborted_newmems_to_submit.extend(current_aborted_newmems_to_submit)
			jobs_aborted_to_restart_updated.extend(current_jobs_aborted_to_restart_updated)
			jobs_aborted_newtimes_updated.extend(current_jobs_aborted_newtimes_updated)
			jobs_aborted_newmems_updated.extend(current_jobs_aborted_newmems_updated)

			jobs_processing_updated.extend(current_aborted_jobs_to_submit)

	# select new jobs to submit
	jobs_processing_updated.extend(aborted_jobs_to_submit)
	new_space_in_queue = space_in_queue - len(aborted_jobs_to_submit)

	[[new_jobs_to_submit],[jobs_to_process_updated]] = \
		Jobs_to_Run_Selector(username,
			[jobs_to_process_updated],
			new_space_in_queue)
		# should return empty output list if input lists are empty

	jobs_processing_updated.extend(new_jobs_to_submit)

	if not jobs_to_process_updated and not jobs_processing_updated:
		# if no more jobs to process, write complete_file
		open(current_completefile,'a').close()

	# ensure jobs only counted once
	jobs_to_process_updated=list(set(jobs_to_process_updated))
	jobs_processing_updated=list(set(jobs_processing_updated))
	jobs_completed_updated=list(set(jobs_completed_updated))
	jobs_with_errors_updated=list(set(jobs_with_errors_updated))
	jobs_aborted_forever_updated=list(set(jobs_aborted_forever_updated))
		# don't mess with order of jobs_aborted_to_restart_updated or
			# the associated jobs_aborted_newtimes_updated and
			# jobs_aborted_newmems_updated

	# update trackfile
	jobs_to_process_updated_string = ';'.join(str(x) for x in
		jobs_to_process_updated)
	jobs_processing_updated_string = ';'.join(str(x) for x in
		jobs_processing_updated)
	jobs_completed_updated_string = ';'.join(str(x) for x in
		jobs_completed_updated)
	jobs_with_errors_updated_string = ';'.join(str(x) for x in
		jobs_with_errors_updated)
	jobs_aborted_forever_updated_string = ';'.join(str(x) for x in
		jobs_aborted_forever_updated)
	jobs_aborted_to_restart_updated_string = ';'.join(str(x) for x in
		jobs_aborted_to_restart_updated)
	jobs_aborted_newtimes_updated_string = ';'.join(str(x) for x in
		jobs_aborted_newtimes_updated)
	jobs_aborted_newmems_updated_string = ';'.join(str(x) for x in
		jobs_aborted_newmems_updated)
	
	trackfile_contents[1][1] = \
		jobs_to_process_updated_string
	trackfile_contents[2][1] = \
		len(jobs_to_process_updated)
	trackfile_contents[1][2] = \
		jobs_processing_updated_string
	trackfile_contents[2][2] = \
		len(jobs_processing_updated)
	trackfile_contents[1][3] = \
		jobs_completed_updated_string
	trackfile_contents[2][3] = \
		len(jobs_completed_updated)
	trackfile_contents[1][4] = \
		jobs_with_errors_updated_string
	trackfile_contents[2][4] = \
		len(jobs_with_errors_updated)
	trackfile_contents[1][5] = \
		jobs_aborted_forever_updated_string
	trackfile_contents[2][5] = \
		len(jobs_aborted_forever_updated)
	trackfile_contents[1][6] = \
		jobs_aborted_to_restart_updated_string
	trackfile_contents[2][6] = \
		len(jobs_aborted_to_restart_updated)
	trackfile_contents[1][7] = \
		jobs_aborted_newtimes_updated_string
	trackfile_contents[2][7] = ''
	trackfile_contents[1][8] = \
		jobs_aborted_newmems_updated_string
	trackfile_contents[2][8] = ''	
	with open(current_trackfile, 'w') as csvfile:
		trackfile_writer = csv.writer(csvfile)
		for trackfile_row in trackfile_contents:
			trackfile_writer.writerow(trackfile_row[:])

	return [new_jobs_to_submit,aborted_jobs_to_submit,aborted_newtimes_to_submit,\
		aborted_newmems_to_submit]

def Consecutive_Parser(data, stepsize = 1):
	# Splits list of integers into list of lists of consecutive numbers
	# Returns the list of consecutive integer lists, as well as a list
		# of lists of the indices of each member of each consecutive
		# integer lists in the original data
	np_data = numpy.array(data)

	sorted_indices = numpy.argsort(np_data)
	sorted_data = np_data[sorted_indices]

	# Find position of indices where the number is higher than previous
		# number + stepsize
	split_indices = numpy.where(numpy.diff(sorted_data) != stepsize)[0]+1

	# Split data at split_indices
	split_data = numpy.split(sorted_data,split_indices)
	# Create list of lists of indices from original data corresponding
		# to every list of values in split_data
	split_sorted_indices = numpy.split(sorted_indices,split_indices)

	return([split_data,split_sorted_indices])

def Job_List_String_Converter(job_list):
	# Converts a list of job numbers into a job submission string for
		# SLURM

	# Split job_list into lists of consecutive jobs that can be
		# submitted in one group to SLURM
	# Consecutive_Parser takes a list of integers, returns numpy
		# array of numpy arrays
	[split_job_list,_] = Consecutive_Parser(map(int,job_list))
	job_string = ''

	# Denote consecutive job sublists by dashes between the min
		# and max of the sublist; separate job sublists by commas
	for current_job_sublist in split_job_list:
		current_job_num = len(current_job_sublist)
		if current_job_num==1:
			current_job_string = (str(current_job_sublist[0])+',')
		else:
			current_job_string = (str(min(current_job_sublist))+'-'+
				str(max(current_job_sublist))+',')
		job_string = job_string + current_job_string

	return(job_string)

def Job_List_Parser(job_list,max_char_num=4096):
	# Taking in a list of jobs, reorganize them so that consecutive jobs
		# can be called as intervals (e.g. '5-8' for '5,6,7,8', and find
		# the optimal arrangement that maximized job number so that the
		# character count for the list of jobs doesn't exceed max_char_num
	# Assumes that max_char_num is sufficiently large
		# i.e. >(2*(max number of digits in a job number)+2)

	# Split job_list into lists of consecutive jobs that can be
		# submitted in one group to SLURM
	# Consecutive_Parser takes a list of integers, returns numpy
		# array of numpy arrays
	[split_job_list,split_index_list] = Consecutive_Parser(map(int,job_list))

	job_group_number = len(split_job_list)
	#string_job_list = []
	char_count_list = numpy.empty(job_group_number)
	job_number_list = numpy.empty(job_group_number)

	# For each set of consecutive jobs in split_job_list, create a
		# string that can be used to submit those jobs to SLURM,
		# and count the number of characters in that string
	# Keep track of the job submission strings (string_job_list),
		# character counts (char_count_list), and number of jobs
		# (job_number_list) in each consecutive list of jobs
	for counter,current_job_sublist in enumerate(split_job_list):
		current_jobs_string = Job_List_String_Converter(current_job_sublist)
		current_job_num = len(current_job_sublist)
		current_char_count = len(current_jobs_string)
		char_count_list[counter] = current_char_count
		job_number_list[counter] = current_job_num
		#string_job_list.append(current_jobs_string)

	# sort character_count_list based on job number submitted in
		# descending order
	order_by_job_num = numpy.argsort(-job_number_list)
	char_count_list_sorted = char_count_list[order_by_job_num]

	# maximize the number of jobs that can be submitted without
		# exceeding max_char_num
	jobs_submitted = []
	jobs_submitted_indices = []
	submission_string = ''
	additional_jobs_possible = True
	remaining_chars = max_char_num-1
	used_positions = numpy.array([False]*job_group_number)
	unused_char_count_list_sorted = char_count_list_sorted[:]

	while additional_jobs_possible:
		# One set at a time, submit sets of jobs that would not equal
			# max_char_num when submitted
		# find the number of characters that will be incurred for every n if n
			# sets of jobs from string_job_list are submitted
		unused_char_count_list_sorted[used_positions] = max_char_num+1
			# eliminate used positions from following rounds of job
				# set selection
		allowed_for_submission = unused_char_count_list_sorted <= remaining_chars
		allowed_indices = numpy.where(allowed_for_submission)[0] # need [0] since numpy.where returns a tuple
		number_submittable_job_lists = sum(allowed_for_submission)
		if number_submittable_job_lists > 0:
			current_sorted_index_to_use = allowed_indices[0]
			current_original_index_to_use = order_by_job_num[current_sorted_index_to_use]
			current_sublist = split_job_list[current_original_index_to_use]
			current_index_sublist = split_index_list[current_original_index_to_use]
			jobs_submitted.extend(current_index_sublist)
			jobs_submitted_indices.extend(current_index_sublist)
			#submission_string = submission_string+string_job_list[current_original_index_to_use]
			remaining_chars = remaining_chars-char_count_list[current_original_index_to_use]
			used_positions[current_sorted_index_to_use] = True
		else:
			additional_jobs_possible = False

	#return([jobs_submitted, jobs_submitted_indices, submission_string])
	return(jobs_submitted_indices)

def Jobs_to_Run_Selector(username,list_of_job_property_lists,space_in_queue,
	max_sbatches_in_one_run = float('Inf')):
	# Based on the amount of space available in the queue, splits jobs
		# into ones that can be run now vs those that have to be run
		# later

	#################
	# SLURM-based cluster doesn't appear to slow down as a penalty for
		# submitting jobs one-by-one so don't reset
		# max_sbatches_in_one_run if submitting jobs one-by-one

	#if len(list_of_job_property_lists) > 1:
	#	max_sbatches_in_one_run = 25
			# if more than ~25 jobs are submitted in a row, submission
				# slows down very significantly, at least on PBS
			# Not sure this is true on slurm, but keep this for now
	#################

	# list_of_job_property_lists contains 1 or more lists, where the
		# first is a list of jobs and the latter are properties of each
		# of the job (e.g. runtime) that need to be retained
	initial_job_list = list_of_job_property_lists[0]

	list_number = len(list_of_job_property_lists)

	# Take the maximum number of jobs that can be submitted from
		# initial_job_list, and parse them to ensure that
		#	1. they're listed in the most efficient SLURM-readable
		#		format
		#	2. this list doesn't exceed the max number of chars
		#		allowed by SLURM
	num_jobs_to_start = min(len(initial_job_list),space_in_queue,max_sbatches_in_one_run)
	initial_truncated_job_list = initial_job_list[0:num_jobs_to_start]
	#[jobs_submitted, jobs_submitted_indices, submission_string] = Job_List_Parser(initial_truncated_job_list)
	if num_jobs_to_start > 0:
		jobs_submitted_indices = Job_List_Parser(initial_truncated_job_list)
	else:
		jobs_submitted_indices = []

	list_of_job_property_lists_to_submit = []
	list_of_job_property_lists_for_later = []

	# Make list of indices that aren't being submitted in this round
	indices_for_later = numpy.setdiff1d(numpy.array(range(0,len(initial_job_list))),
		jobs_submitted_indices)

	for sublist_idx,current_sublist in enumerate(list_of_job_property_lists):
		list_of_job_property_lists_to_submit.append([])
		list_of_job_property_lists_for_later.append([])
		for current_idx_to_submit in jobs_submitted_indices:
			list_of_job_property_lists_to_submit[sublist_idx].append(current_sublist[current_idx_to_submit])
		for current_idx_for_later in indices_for_later:
			list_of_job_property_lists_for_later[sublist_idx].append(current_sublist[current_idx_for_later])

	return([list_of_job_property_lists_to_submit,list_of_job_property_lists_for_later])

def Create_Slurm_Job(single_job_time,single_job_mem,current_sbatch_filename,
	parallel_processors,code_dir,job_list_string,user_email,module_to_use,
	code_run_string,current_job_name,additional_beginning_lines_in_sbatch,
	additional_end_lines_in_sbatch):
	# Writes .q files to submit to cluster queue

	# convert single_job_time from minutes to hh:mm:ss	
	single_job_hours = int(math.floor(float(single_job_time)/60))
	if single_job_hours > 0:
		single_job_minutes = int(float(single_job_time) % (float(single_job_hours)*60))
	else:
		single_job_minutes = int(float(single_job_time))
	single_job_time_string = str(single_job_hours).zfill(2)+':'+ \
		str(single_job_minutes).zfill(2)+':00'

	# convert single_job_mem from # of Mb to memory string
	# cluster only cares about units of GB
	single_job_mem_GB = math.ceil(float(single_job_mem)/1024)
	single_job_mem_string = str(int(single_job_mem_GB))+'GB'

	with open(current_sbatch_filename,'w') as sbatch_job_file:
		sbatch_job_file.write('#!/bin/bash\n')

		sbatch_job_file.write('#SBATCH --job-name='+current_job_name+'\n')
			# name of current job
		sbatch_job_file.write('#SBATCH --output='+current_job_name+'.o%A-%a\n')
		sbatch_job_file.write('#SBATCH --error='+current_job_name+'.e%A-%a\n')
		sbatch_job_file.write('#SBATCH --time='+single_job_time_string+'\n')
			# amount of time allowed per job
		sbatch_job_file.write('#SBATCH --mem='+single_job_mem_string+'\n')
			# memory allocated to each job
		sbatch_job_file.write('#SBATCH --nodes=1\n')
		sbatch_job_file.write('#SBATCH --cpus-per-task='+str(parallel_processors)+'\n')
			# nodes and processors used per job
		sbatch_job_file.write('#SBATCH --array='+job_list_string+'\n')
			# list of job IDs to be submitted
		sbatch_job_file.write('#SBATCH --mail-type=FAIL\n')
			# only sends email if job array fails
		sbatch_job_file.write('#SBATCH --mail-user='+user_email+'\n')
			# email that gets notification about aborted job

		if additional_beginning_lines_in_sbatch:
			for additional_sbatch_beginning_row in additional_beginning_lines_in_sbatch:
				sbatch_job_file.write(additional_sbatch_beginning_row+'\n')

		sbatch_job_file.write('cd '+code_dir+'\n')
		sbatch_job_file.write('module purge\n')
			# cd into code directory
		if module_to_use == 'matlab':
			sbatch_job_file.write('module load matlab/2017a\n')
				# load matlab module
			#sbatch_job_file.write('matlab -nodisplay -nosplash -nodesktop -r '+
			sbatch_job_file.write('matlab -nodisplay -r '+
				code_run_string+'\n')
				# write appropriate code-running line
		elif module_to_use == 'R':
			sbatch_job_file.write('module load r/intel/3.3.2\n')
				# load R module
			sbatch_job_file.write('R CMD BATCH --vanilla '+code_run_string+'\n')

		if additional_end_lines_in_sbatch:
			for additional_sbatch_end_row in additional_end_lines_in_sbatch:
				sbatch_job_file.write(additional_sbatch_end_row+'\n')

		sbatch_job_file.write('\n\n')
			# need additional returns at end of shell scripts

	return None

def Job_List_Submitter(new_jobs_to_submit,aborted_jobs_to_submit,aborted_newtimes_to_submit,
	aborted_newmems_to_submit,current_sbatch_filename,parallel_processors,
	code_dir,user_email,module_to_use,code_run_string,current_job_name,
	additional_beginning_lines_in_sbatch,additional_end_lines_in_sbatch,
	default_single_job_time,default_single_job_mem):

	if aborted_jobs_to_submit:

		# submit aborted jobs one at a time, with appropriate memory and time requirements
		for current_aborted_job_index, current_aborted_job in enumerate(aborted_jobs_to_submit):
			current_single_job_time =  aborted_newtimes_to_submit[current_aborted_job_index] # in minutes
			current_single_job_mem = aborted_newmems_to_submit[current_aborted_job_index] # in MB
			job_list_string = str(current_aborted_job)

			############### ??? TO DO ??? ###############
			# Need to update treatment of aborted jobs here
			############### ??? TO DO ??? ###############

			Create_Slurm_Job(current_single_job_time,current_single_job_mem,current_sbatch_filename,
				parallel_processors,code_dir,job_list_string,user_email,
				module_to_use,code_run_string,current_job_name,
				additional_beginning_lines_in_sbatch,additional_end_lines_in_sbatch)

			# cd into sbatch directory
			os.chdir(slurm_folder)			
			# Run .q file for sim
			subprocess.call('sbatch '+current_sbatch_filename,shell=True)

	if new_jobs_to_submit:

		#job_list_string = ','.join(str(x) for x in new_jobs_to_submit)
		job_list_string = Job_List_String_Converter(new_jobs_to_submit)

		Create_Slurm_Job(default_single_job_time,default_single_job_mem,current_sbatch_filename,
			parallel_processors,code_dir,job_list_string,user_email,
			module_to_use,code_run_string,current_job_name,
			additional_beginning_lines_in_sbatch,additional_end_lines_in_sbatch)

		# cd into sbatch directory
		os.chdir(slurm_folder)			
		# Run .q file for sim

		subprocess.call('sbatch '+current_sbatch_filename,shell=True)



	return None