#!/usr/bin/python

# Estimates mutational effects and proportion petites for each strain
	# based on differences in average growth rates between mean test
	# and reference strain GRs in a group within a field

import sys, os, getopt
import csv
import subprocess
import math
import numpy
from scipy.stats import chi2
import re
from shutil import copyfile
import pandas

usage = '''

Usage:

    MLE_runner_slurm_GRdiff.py -s setup_file

Defaulting to ~/setup_data_MLE.csv

'''

# code below necessary to avoid errors in reading csv files with large fields
	# from user1251007 on stackoverflow
maxInt = sys.maxsize
decrement = True

while decrement:
    # decrease the maxInt value by factor 10 
    # as long as the OverflowError occurs.

    decrement = False
    try:
        csv.field_size_limit(maxInt)
    except OverflowError:
        maxInt = int(maxInt/10)
        decrement = True

def Check_Input():
	# gets input passed to program, which should contain filename of
		# setup file
	if len(sys.argv) <1 :
		print >> sys.stderr, 'Setup file needed.'
		print >> sys.stderr, usage
	else:
		opts, args = getopt.gnu_getopt(sys.argv[1:], "s:h", ["--help"])

	setup_file = '~/setup_data.csv'
	try:
		for o, a in opts:
			if o == "-h" or o == "--help":
				print usage
				exit()
			if o == '-s':
				setup_file = a
	except:
		print >> sys.stderr, 'Setup file needed.'
		print >> sys.stderr, usage
		exit()

	return (setup_file)

def Read_Setup_File(file_name):
	# Processes file_name to extract info about directories in which to
		# find code and data, etc

	with open(file_name, 'rU') as data_file:
		#setup_data = list(csv.reader(data_file, dialect='excel', quotechar='"'))
		setup_data = list(csv.reader(data_file))

#	data_file = open(file_name,'rU')
#	setup_data = list(csv.reader(data_file, dialect='excel', quotechar='"'))

	username = setup_data[1][1]
	home_dir = setup_data[2][1]
	composite_data_dir = setup_data[3][1]
	backup_dir = setup_data[4][1]
	code_dir = setup_data[5][1]
	temp_storage_dir = setup_data[6][1] # /scratch/yp19/mut_effect_estimation
	starting_L = float(setup_data[7][1])
		# original length of x-axis used for pdf in MLE
	starting_gridpower = int(setup_data[8][1])
		# log2(original number of x-axis gridpoints used in MLE)
	sim_repeats = int(setup_data[9][1])
		# number of repeats of simulation
	profile_points_by_mode_string = str(setup_data[10][1])
		# of points to evaluate in likelihood profile for each parameter
	MLE_memory = float(setup_data[11][1])
		# amount of Mb required to run MLE
	MLE_time = float(setup_data[12][1])
		# amount of time required to run MLE in minutes
	user_email = setup_data[13][1]
		# email of user who gets email notifications from these jobs
	strain_number = int(setup_data[14][1])
		# number of strains phenotyped
		# this number will also be simulated
	mode_list = setup_data[15][1].split('|') # which modes should be run

	# Each mode has an associated list of fixed_parameters to try,
		# as well as minimum and maximum fixed values to use for each
		# of these parameters
	# modes are separated by vertical bars while parameters are
		# separated by semicolons
	parameters_by_mode_string = setup_data[16][1]
	parameter_min_by_mode_string = setup_data[17][1]
	parameter_max_by_mode_string = setup_data[18][1]
	parameter_start_by_mode_string = setup_data[19][1]
		# list of parameters corresponding to each mode in mode_list;
		# this is meant to produce a list of lists; in
			# parameters_by_mode_string, modes are separated by
			# vertical bars, whereas parameter names are separated by
			# semicolons

	pdf_points = int(setup_data[20][1])
		# number of gridpoints in empirical pdf returned for each mode
	parallel_processors_max = int(setup_data[21][1])
		# max number of parallel processors to use

	lambda_SNM = float(setup_data[22][1])
		# empirically determined mean number of SNPs per strain

	max_memory = float(setup_data[23][1])
		# max amount of Mb allowed for one job
	max_time = float(setup_data[24][1])
		# max amount of mins allowed for one job
	ms_positions = int(setup_data[25][1])

	parameter_profile_lb_by_mode_string = setup_data[26][1]
	parameter_profile_ub_by_mode_string = setup_data[27][1]

	profile_points_by_mode = [map(int,current_param_string.split(';')) \
		for current_param_string in profile_points_by_mode_string.split('|')]

	parameters_by_mode = [current_param_string.split(';') \
		for current_param_string in parameters_by_mode_string.split('|')]

#	parameter_max_by_mode_list = [current_param_max_string.split(';') \
#		for current_param_max_string in \
#		parameter_max_by_mode_string.split('|')]

#	parameter_max_by_mode = [[map(float,current_val) for current_val in \
#		parameter_max_by_mode_list] for parameter_max_by_mode_sublist in \
#		parameter_max_by_mode_list]

	parameter_max_by_mode = [map(float,current_param_max_string.split(';')) \
		for current_param_max_string in \
		parameter_max_by_mode_string.split('|')]

	parameter_min_by_mode = [map(float,current_param_min_string.split(';')) \
		for current_param_min_string in \
		parameter_min_by_mode_string.split('|')]

	parameter_start_by_mode = [map(float,current_param_start_string.split(';')) \
		for current_param_start_string in \
		parameter_start_by_mode_string.split('|')]

	parameter_profile_ub_by_mode = [map(float,current_param_max_string.split(';')) \
		for current_param_max_string in \
		parameter_profile_ub_by_mode_string.split('|')]

	parameter_profile_lb_by_mode = [map(float,current_param_min_string.split(';')) \
		for current_param_min_string in \
		parameter_profile_lb_by_mode_string.split('|')]

	return (username,home_dir,composite_data_dir,backup_dir,code_dir,
		temp_storage_dir,starting_L,starting_gridpower,sim_repeats,
		profile_points_by_mode,MLE_memory,MLE_time,user_email,strain_number,
		mode_list,parameters_by_mode,parameter_max_by_mode,
		parameter_min_by_mode,parameter_start_by_mode,pdf_points,parallel_processors_max,lambda_SNM,
		max_memory,max_time,ms_positions,parameter_profile_ub_by_mode,parameter_profile_lb_by_mode)

def Setup_Directories(current_folder,temp_storage_dir,output_dir_name,username):
	trackfile_folder = current_folder + '/trackfiles'
	completefile_folder = current_folder + '/completefiles'
	slurm_folder = temp_storage_dir + '/' + output_dir_name + '/slurm_folder'
	sim_output_folder = temp_storage_dir +  '/' + output_dir_name \
		+ '/simulated_phenotypes'
	MLE_output_folder = temp_storage_dir +  '/' + output_dir_name \
		+ '/MLE_output'
	LL_profile_folder = current_folder + '/LL_profiles'
#	epilogue_folder = '/home/' + username + '/mut_effect_epilogue_files/' \
#		+ output_dir_name
	MLE_combined_sim_folder = current_folder + '/MLE_sim_outputs'

#	new_directory_list = (trackfile_folder,completefile_folder,slurm_folder,
#		sim_output_folder,MLE_output_folder,LL_profile_folder,epilogue_folder,
#		MLE_combined_sim_folder)
	new_directory_list = (trackfile_folder,completefile_folder,slurm_folder,
		sim_output_folder,MLE_output_folder,LL_profile_folder,
		MLE_combined_sim_folder)

	setup_complete_file = completefile_folder + '/folder_setup_complete.txt'
	if not os.path.isfile(setup_complete_file):
		for current_new_directory in new_directory_list:
			if not os.path.isdir(current_new_directory):
				os.makedirs(current_new_directory)
		open(setup_complete_file,'a').close()

	return (new_directory_list)

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
#	max_allowed_jobs = round(max_allowed_jobs*.9)

	# find max amount of jobs that can be added to queue without
		# making it overflow
	space_in_queue = max_allowed_jobs-jobs_in_queue

	return(space_in_queue)

#def CI_position_calculator(LL_a,LL_b,param_a,param_b,CI_boundary):
#	# calculates slope and intercept between points [param_a,LL_a] and
#		# [param_b,LL_b], and then calculates the parameter value
#		# (param_c) that would intersect CI_boundary with those slope
#		# and intercept
#
#	slope = (LL_b-LL_a)/(param_b-param_a)
#	intercept = LL_b-slope*param_b
#	param_c = (CI_boundary-intercept)/slope
#
#	return(param_c)

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
	# convert all lists to sets and back to ensure you have unique values
#	jobs_to_process_updated = list(sorted(set(jobs_to_process_updated)))
#	jobs_processing_updated = list(sorted(set(jobs_processing_updated)))
#	jobs_completed_updated = list(sorted(set(jobs_completed_updated)))
#	jobs_with_errors_updated = list(sorted(set(jobs_with_errors_updated)))

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

def Simulation_Setup(sim_repeats,current_mode,
	username,output_dir_name,slurm_folder,user_email,strain_number,
	MLE_output_file,growth_condition,completefile_folder,trackfile_folder,
	sim_output_folder,original_phenotype_file,epilogue_folder,max_memory,
	max_time,lambda_SNM):
	# Sets up the simulation of each type of data based on the ML value

	sim_completefile = completefile_folder + '/' + current_mode + '_' + growth_condition + '_sim_completefile.txt'

	if not os.path.isfile(sim_completefile):
		# Update the trackfile; get list of jobs, if any, for which
			# simulation still needs to be run
		sim_trackfile = trackfile_folder + '/' + current_mode + '_' + growth_condition + '_sim_trackfile.csv'

		current_job_name = output_dir_name + '-sim-' + current_mode + '_'\
			+ growth_condition
		current_output_filename = 'phenotype_table_' + current_mode + '_'\
			 + growth_condition
		current_output_extension = 'csv'

		default_single_job_time = 3 # in minutes
		default_single_job_mem = 300 # in MB

		[new_jobs_to_submit,aborted_jobs_to_submit,aborted_newtimes_to_submit,\
			aborted_newmems_to_submit] = Trackfile_Processor(sim_repeats,
				current_mode,sim_trackfile,sim_completefile,username,
				current_job_name,slurm_folder,sim_output_folder,
				current_output_filename,current_output_extension,
				epilogue_folder,output_dir_name,max_memory,max_time,
				default_single_job_mem,default_single_job_time)

		current_sbatch_filename = slurm_folder + '/' + current_job_name + '.q'
		parallel_processors = 1
		code_run_string = ('\"MC_strain_simulator_' + current_mode + '(\'' 
			+ MLE_output_file + '\',\'' + original_phenotype_file + '\',\'' + sim_output_folder
			+ '\',\'' + growth_condition + '\',' + str(strain_number) + ',' + str(lambda_SNM)
			+ ',\"$PBS_ARRAYID\");exit\"')
		additional_beginning_lines_in_sbatch = []
		additional_end_lines_in_sbatch = []
		module_to_use = 'matlab'

		# submit aborted jobs one-by-one, new_jobs_to_submit as a batch job
		Job_List_Submitter(new_jobs_to_submit,aborted_jobs_to_submit,aborted_newtimes_to_submit,
			aborted_newmems_to_submit,current_sbatch_filename,
			parallel_processors,code_dir,user_email,module_to_use,
			code_run_string,current_job_name,
			additional_beginning_lines_in_sbatch,additional_end_lines_in_sbatch,
			default_single_job_time,default_single_job_mem)

	return sim_completefile

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

	current_MLE_completefile = completefile_folder + '/MLE_' + growth_condition \
		+ '_' + current_mode + '_' + current_fixed_parameter \
		+ '_' + rep + '_completefile.txt'

	if not os.path.isfile(current_MLE_completefile):
		# Update the trackfile; get list of jobs, if any, for which
			# simulation still needs to be run

		# For 'sim', mat_phenotype_file needs to include the job iteration
		if rep == 'sim':
			#phenotype_table_lambdafit_',growth_condition,'_',num2str(ID),'.csv'
			mat_phenotype_file = phenotype_file + '${SLURM_ARRAY_TASK_ID}.csv'
		else:
			mat_phenotype_file = phenotype_file
#			# MLE needs to be run only once in 'unfixed' condition,
#				# so 'profile_points' needs to be 1
#			if current_fixed_parameter == 'unfixed':
#				profile_points = 1

		current_MLE_trackfile = trackfile_folder + '/MLE_' + current_mode \
			+ '_' + current_fixed_parameter + '_' + growth_condition + '_' \
			+ rep + '_trackfile.csv'
		current_job_name = output_dir_name + '-MLE-' + current_mode + '_rep' + rep + '_'\
			+ current_fixed_parameter + '_' + growth_condition
		
		current_output_filename = growth_condition + '_data_' + current_mode + '_' + current_fixed_parameter
		current_output_extension = 'csv'
		folder_with_csv_output_files = MLE_output_folder + '/csv_output'

#		if current_mode =='strain_pairwise_diff':
			#parallel_processors = min(20,strain_number)
		if current_fixed_parameter=='unfixed':
			fitted_dimensions = 4
		else:
			fitted_dimensions = 3
		parallel_processors = min(parallel_processors_max,ms_positions**fitted_dimensions)
#		else:
#			parallel_processors = 1
#			additional_beginning_lines_in_sbatch = []
#			additional_end_lines_in_sbatch = []

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

def MLE_combiner(completefile_folder,growth_condition,current_mode,rep,
	current_folder,MLE_summary_file,profile_points_to_loop_over,fixed_parameter_list,
	MLE_output_folder,LL_profile_folder,lambda_SNM,total_parameter_list):
	# combines output of of ML_Estimator and writes, for each 
		# growth_condition and mode:
		# 1. LL_file, which contains all MLE results for each fixed parameter
		# 2. MLE_summary_file, which contains 
	current_ML_combination_completefile = completefile_folder \
		+ '/ML_combination_' + growth_condition + '_' + current_mode + '_' \
		+ rep + '_completefile.txt'

	if not os.path.isfile(current_ML_combination_completefile):

		total_parameter_number = len(total_parameter_list)
		max_profile_points = numpy.max(profile_points_to_loop_over)

		MLE_summary_matrix = numpy.empty(((total_parameter_number+3),3))
		MLE_summary_matrix[:] = numpy.NAN


		unfixed_MLE_output = (MLE_output_folder + '/csv_output/'
				+ growth_condition + '_data_' + current_mode + '_unfixed_1.csv')

		MLE_output_param_and_LL_indices = [0]+range(2,total_parameter_number+2)

		warning_line_list = ['']*(total_parameter_number+3)

		# include highest and lowest values with calculated LL above
			# and below CI, respectively
		highest_calc_below_CI = numpy.empty((total_parameter_number+3))
		lowest_calc_above_CI = numpy.empty((total_parameter_number+3))
		highest_calc_below_CI[:] = numpy.nan
		lowest_calc_above_CI[:] = numpy.nan

		runtime_list = numpy.array([])

		if os.path.isfile(unfixed_MLE_output):
			with open(unfixed_MLE_output, 'rU') as unfixed_MLE_output_contents:
				unfixed_MLE_output_data = list(csv.reader(unfixed_MLE_output_contents))

			unfixed_MLE_np = numpy.array(map(float,unfixed_MLE_output_data[0]))
			current_time = unfixed_MLE_np[1]
			runtime_list = numpy.append(runtime_list,current_time)
				# warning - doesn't return list in python 3
			current_warning_line = ''
		else:
			current_warning_line = 'warning! no unfixed condition file found;'
			unfixed_MLE_np = numpy.array([float('NaN')]*(total_parameter_number+2))

		max_likelihood_params = unfixed_MLE_np

		warning_line_list[1] = current_warning_line

		# make an array that will hold all likelihood profile results

		LL_combined_array = numpy.empty((max_profile_points+1,total_parameter_number+2,len(fixed_parameter_list)))
		LL_combined_array[:] = numpy.NAN # avoid problems with numpy array
		# Populate each LL_array
		for current_fixed_parameter_index, current_fixed_parameter in enumerate(fixed_parameter_list):

			current_index_in_full_param_list = total_parameter_list.index(current_fixed_parameter)
			current_profile_point_num = profile_points_to_loop_over[(current_fixed_parameter_index+1)]

			current_warning_line = warning_line_list[(current_index_in_full_param_list+2)]

			# Loop through all profile_points_to_loop_over of current_fixed_parameter and
				# populate LL_array with results
			LL_array = numpy.empty((max_profile_points+1,total_parameter_number+2))
			LL_array[:] = numpy.NAN # avoid problems with numpy array

			# leave the first row to be filled later with max_likelihood_params
				# the other rows should be the profile points

			for current_profile_point in range(1,current_profile_point_num+1):
				current_MLE_output = (MLE_output_folder + '/csv_output/'
					+ growth_condition + '_data_' + current_mode + '_' + current_fixed_parameter
					+ '_' + str(current_profile_point) + '.csv')
				if os.path.isfile(current_MLE_output):
					with open(current_MLE_output, 'rU') as current_MLE_output_contents:
						current_MLE_output_data = list(csv.reader(current_MLE_output_contents))

					current_MLE_output_np = numpy.array(map(float,current_MLE_output_data[0]))
					current_time = current_MLE_output_np[1]
					runtime_list = numpy.append(runtime_list,current_time)
						# warning - doesn't return list in python 3

					LL_array[current_profile_point,:] = current_MLE_output_np

					# in case 'unfixed' search didn't find actual MLE, update max_likelihood_params
					if numpy.isnan(max_likelihood_params[0]) or current_MLE_output_np[0] > max_likelihood_params[0]:
						max_likelihood_params = current_MLE_output_np
						current_warning_line = current_warning_line + \
							'warning! unfixed condition did not find true ML! ML found in fixed ' + \
							current_fixed_parameter + ' search;'

			LL_combined_array[:,:,current_fixed_parameter_index] = LL_array

			warning_line_list[(current_index_in_full_param_list+2)] = current_warning_line

		# Fill first position in each LL array with max_likelihood_params
		LL_combined_array[0,:,:] = numpy.transpose(numpy.tile( \
			max_likelihood_params,(len(fixed_parameter_list),1)))

		# Loop through parameters again, sorting and saving LL arrays for each parameter
		for current_fixed_parameter_index, current_fixed_parameter in enumerate(fixed_parameter_list):

			current_index_in_full_param_list = total_parameter_list.index(current_fixed_parameter)
			current_LL_file = (LL_profile_folder + '/' + growth_condition + '_'
				+ current_mode + '_' + current_fixed_parameter + '_' + rep
				+ '_LL_file.csv')
			# sort LL_array by the fixed_parameter
				# drop all-NaN rows later

			LL_array = LL_combined_array[:,:,current_fixed_parameter_index]
			LL_array_sorted = LL_array[LL_array[:,current_index_in_full_param_list+2].argsort()]
			numpy.sort(LL_array_sorted)
#			numpy.savetxt(current_LL_file, LL_array_sorted, delimiter=",")
			LL_df = pandas.DataFrame(
				columns = (['LL','runtime_in_mins'] + total_parameter_list),
				data = LL_array_sorted)
			LL_df.to_csv(path_or_buf=current_LL_file,index=False)

			LL_combined_array[0:numpy.shape(LL_array_sorted)[0], \
				0:numpy.shape(LL_array_sorted)[1], \
				current_fixed_parameter_index] = LL_array_sorted

		MLE_summary_matrix[1:-1,0] = max_likelihood_params[MLE_output_param_and_LL_indices]
		max_likelihood = max_likelihood_params[0]

		conf_int_cutoff = max_likelihood-chi2.ppf(.975,1)/2

		runtime_mean = numpy.mean(runtime_list)
		runtime_se = numpy.std(runtime_list)/numpy.sqrt(len(runtime_list))
		runtime_params = numpy.array([runtime_mean,(runtime_mean-1.96*runtime_se),(runtime_mean+1.96*runtime_se)])/3600

		MLE_summary_matrix[[1],[1,2]] = conf_int_cutoff
		MLE_summary_matrix[0,0:3] = runtime_params
		MLE_summary_matrix[-1,0:3] = lambda_SNM

		# Loop through parameters again, this time looking for points
			# bounding conf_int_cutoff; if first point is within C.I.,
			# set lower bound to -Inf; if last point is within C.I., set
			# upper bound to Inf
		for current_fixed_parameter_index, current_fixed_parameter in enumerate(fixed_parameter_list):

			current_index_in_full_param_list = total_parameter_list.index(current_fixed_parameter)

			current_warning_line = warning_line_list[(current_index_in_full_param_list+2)]
			current_LL_array = numpy.squeeze(LL_combined_array[:,:,current_fixed_parameter_index])
			# drop rows containing all NaN
			all_nan_rows = numpy.all(numpy.isnan(current_LL_array),axis=1)
			current_LL_array_clean = current_LL_array[~all_nan_rows,:]

			real_profile_point_number = (current_LL_array_clean).shape[0]

			if real_profile_point_number > 0:
				# important to treat current_CI_lower_bound and current_CI_upper_bound
					# as lists here, so that you can check for emptiness
				if current_LL_array_clean[real_profile_point_number-1,0] > conf_int_cutoff:
					current_CI_upper_bound_list = numpy.array([float('Inf')])
				else:
					current_CI_upper_bound_list = numpy.array([])

				# if point of lowest param value is inside CI, set CI bound to -Inf
				if current_LL_array_clean[0,0] > conf_int_cutoff:
					current_CI_lower_bound_list = numpy.array([float('-Inf')])
				else:
					current_CI_lower_bound_list = numpy.array([])

				# Find positions where likelihood profile intersects conf_int_cutoff
				y_vals = current_LL_array_clean[:,0]
				x_vals = current_LL_array_clean[:,current_index_in_full_param_list+2]

				# rounding errors can cause incorrect behavior of the above,
					# leading to x-positions following a very large negative
					# number to be incorrectly identified as being a
					# confidence interval intersection point; thus, remove
					# positions with large absolute values
				x_vals[numpy.absolute(y_vals) > pow(2,51)] = numpy.nan
				y_vals[numpy.absolute(y_vals) > pow(2,51)] = numpy.nan

				if not all(numpy.isnan(y_vals)):
					current_param_MLE_val = max_likelihood_params[current_index_in_full_param_list+2]
					# Look for predicted intersections based on slope and
						# intercept of each interpoint line segment, which are
						# within that line segment
					y_diffs = numpy.diff(y_vals)
					x_diffs = numpy.diff(x_vals)
					slopes = y_diffs/x_diffs
					intersection_predictions = (conf_int_cutoff-y_vals[0:-1])*(1/slopes)+x_vals[0:-1]
					# Filter predicted intersections of conf_int_cutoff to
						# include only ones within the segment used for prediction
					intersection_positions = intersection_predictions[(intersection_predictions >= \
						x_vals[0:-1])*(intersection_predictions < x_vals[1:])]

					current_CI_lower_bound_list = numpy.append(current_CI_lower_bound_list, \
						intersection_positions[intersection_positions <= current_param_MLE_val])
					current_CI_upper_bound_list = numpy.append(current_CI_upper_bound_list, \
						intersection_positions[intersection_positions > current_param_MLE_val])

					# In simple and accurately estimated LL landscape, LL
						# profile expected to increase monotonically up until
						# max LL, then decrease monotonically after; if this
						# isn't the case, throw a warning
					increasing_LL_section = y_diffs[x_vals[:-1] < current_param_MLE_val]
					decreasing_LL_section = y_diffs[x_vals[:-1] >= current_param_MLE_val]

						# remember indexing is different for y_diffs and y_vals,
							# it's correct here
					monotonicity_state = numpy.all(increasing_LL_section >= 0) and numpy.all(decreasing_LL_section <= 0)
					if not monotonicity_state:
						current_warning_line = current_warning_line + 'warning! non-monotonic LL profile before or after ML;'

				if len(current_CI_lower_bound_list) > 0:
					current_CI_lower_bound = numpy.min(current_CI_lower_bound_list)
					if current_CI_lower_bound < numpy.min(x_vals):
						current_highest_calc_below_CI = numpy.nan
					else:
						current_highest_calc_below_CI = numpy.max(x_vals[x_vals < current_CI_lower_bound])
				else:
					current_CI_lower_bound = numpy.nan
					current_highest_calc_below_CI = numpy.nan

				if len(current_CI_upper_bound_list) > 0:
					current_CI_upper_bound = numpy.min(current_CI_upper_bound_list)
					if current_CI_upper_bound > numpy.max(x_vals):
						current_lowest_calc_above_CI = numpy.nan
					else:
						current_lowest_calc_above_CI = numpy.min(x_vals[x_vals > current_CI_upper_bound])
				else:
					current_CI_upper_bound = numpy.nan
					current_lowest_calc_above_CI = numpy.nan

			else:
				current_CI_lower_bound = numpy.nan
				current_CI_upper_bound_list = numpy.nan
				current_highest_calc_below_CI = numpy.nan
				current_lowest_calc_above_CI = numpy.nan

			MLE_summary_matrix[current_index_in_full_param_list+2,1] = current_CI_lower_bound
			MLE_summary_matrix[current_index_in_full_param_list+2,2] = current_CI_upper_bound
			warning_line_list[(current_index_in_full_param_list+2)] = current_warning_line
			highest_calc_below_CI[current_index_in_full_param_list+2] = current_highest_calc_below_CI
			lowest_calc_above_CI[current_index_in_full_param_list+2] = current_lowest_calc_above_CI

		# save output by creating a pandas df
		output_df = pandas.DataFrame(columns = ['Parameter','LL_val', \
			'highest_calculated_LL_val_below_CI','lowest_calculated_LL_val_above_CI', \
			'CI_lower_bound','CI_upper_bound','Warnings'])
		output_df.Parameter = ['average_runtime_in_hrs','Overall_Max_LL'] + total_parameter_list + ['global_fixed_lambda_SNM']
		output_df.LL_val = MLE_summary_matrix[:,0]
		output_df.CI_lower_bound = MLE_summary_matrix[:,1]
		output_df.CI_upper_bound = MLE_summary_matrix[:,2]
		output_df.highest_calculated_LL_val_below_CI = highest_calc_below_CI
		output_df.lowest_calculated_LL_val_above_CI = lowest_calc_above_CI
		output_df.Warnings = warning_line_list
		output_df.to_csv(path_or_buf=MLE_summary_file,index=False)

	#		numpy.savetxt(MLE_summary_file, MLE_summary_matrix, delimiter=",")
	#		with open(MLE_summary_file, "a") as MLE_summary_openfile:
	#			MLE_summary_openfile.write(warning_line)
		open(current_ML_combination_completefile,'a').close()		

		############### ??? TO DO ??? ###############
		# at some point need to get rid of illegitimate values (ones where parameter has bumped into limits)
		############### ??? TO DO ??? ###############

	return(current_ML_combination_completefile)

def pdf_Plotter(completefile_folder,trackfile_folder,pdf_points,
	current_mode,username,slurm_folder,current_folder,
	growth_condition,current_L,current_gridpower,MLE_summary_file,
	phenotype_file,rep,lambda_SNM,epilogue_folder,output_dir_name,
	max_memory,max_time,use_single_mutation):
	# runs MLE

	if use_single_mutation:
		single_mutation_bool = 'true'
		pdf_type = 'single-mut'
	else:
		single_mutation_bool = 'false'
		pdf_type = 'multiple-mut'

	current_pdf_plotter_completefile = completefile_folder + '/pdf_plotter_' + growth_condition \
		+ '_' + current_mode + '_' + rep + '_' + pdf_type + '_completefile.txt'

	if not os.path.isfile(current_pdf_plotter_completefile):
		# Update the trackfile; get list of jobs, if any, for which
			# simulation still needs to be run

		current_pdf_plotter_trackfile = trackfile_folder + '/pdf_plotter_' + current_mode \
			+ '_' + growth_condition + '_' \
			+ rep + '_' + pdf_type + '_trackfile.csv'
		current_job_name = output_dir_name + '-pdf_plotter-' + current_mode + '_rep' + rep + '_' + pdf_type \
			+ '_' + growth_condition

		default_single_job_time = 30 # in minutes
		default_single_job_mem = 1024*10 # in MB
		
		current_output_filename = current_mode + '_' + growth_condition + '_' + pdf_type + '_distribution'
		current_output_extension = 'csv'
		[new_jobs_to_submit,aborted_jobs_to_submit,aborted_newtimes_to_submit,aborted_newmems_to_submit] = Trackfile_Processor(1,current_mode,
			current_pdf_plotter_trackfile,current_pdf_plotter_completefile,username,
			current_job_name,
			slurm_folder,current_folder,current_output_filename,current_output_extension,epilogue_folder,output_dir_name,max_memory,max_time,default_single_job_mem,default_single_job_time)

		if new_jobs_to_submit or aborted_jobs_to_submit:
			# write a sbatch file to submit a batch process
			current_sbatch_filename = slurm_folder + '/' + current_job_name + '.q'
			# when doing MLE in 'knownmuts' mode, LL for each strain
				# has to be calculated separately, so MLE code uses
				# parfor loop, making code run faster on multiple
				# processors
			parallel_processors = 1
			additional_beginning_lines_in_sbatch = []
			additional_end_lines_in_sbatch = []
			module_to_use = 'matlab'
			code_run_string = ('\'Distribution_Plotter(\'' + '\"\''
				+ current_folder + '\"\"\',\"\"\'' + growth_condition
				+ '\"\"\',\"\"\'' + current_mode + '\"\"\',\"\"'
				+ str(current_L) + '\"\",\"\"' + str(int(current_gridpower))
				+ '\"\",\"\"' + str(lambda_SNM)+ '\"\"\',\"\"\'' 
				+ phenotype_file + '\"\"\',\"\"' + str(pdf_points)
				+ '\"\",\"\"\'' + MLE_summary_file
				+ '\"\"\',\"\"\'' + rep + '\"\"\',\"\"' + single_mutation_bool + '\"\");exit\"')

			# submit aborted jobs one-by-one, new_jobs_to_submit as a batch job
			Job_List_Submitter(new_jobs_to_submit,aborted_jobs_to_submit,aborted_newtimes_to_submit,
				aborted_newmems_to_submit,current_sbatch_filename,
				parallel_processors,code_dir,user_email,module_to_use,
				code_run_string,current_job_name,
				additional_beginning_lines_in_sbatch,additional_end_lines_in_sbatch,
				default_single_job_time,default_single_job_mem)

	return current_pdf_plotter_completefile

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

def Simulation_Processor(completefile_folder,current_mode,rep,
	growth_condition,true_MLE_output_file,sim_repeats,fixed_parameter_list,
	MLE_combined_sim_folder,sim_MLE_output_folder,current_fixed_parameter):
	# combines output of of ML_Estimator and writes, for each 
		# growth_condition and mode:
		# 1. LL_file, which contains all MLE results for each fixed parameter
		# 2. MLE_summary_file, which contains 
	current_ML_combination_completefile = completefile_folder \
		+ '/ML_combination_' + growth_condition + '_' + current_mode + '_' \
		+ rep + '_completefile.txt'

	if not os.path.isfile(current_ML_combination_completefile):

		percentile_summary_matrix = numpy.empty((1,2))
		percentile_summary_matrix[:] = numpy.NAN
		percentile_summary_file = (MLE_combined_sim_folder + '/' + growth_condition + '_'
			+ current_mode + '_sim_MLE_percentile.csv')

		with open(true_MLE_output_file, 'rU') as true_MLE_output_contents:
			true_MLE_output_data = list(csv.reader(true_MLE_output_contents))

		true_LL = float(true_MLE_output_data[0][0])

		percentile_summary_matrix[0][0] = true_LL

		current_sim_combo_file = (MLE_combined_sim_folder + '/' + growth_condition + '_'
			+ current_mode + '_sim_MLE_output_file.csv')

		# Loop through all simulations of current_mode and
			# populate sim_output_array with results
		sim_output_array = numpy.empty((sim_repeats,len(fixed_parameter_list)+3))
		sim_output_array[:] = numpy.NAN # avoid problems with numpy array

		for current_sim in range(1,sim_repeats+1):
			current_MLE_output = (sim_MLE_output_folder + '/csv_output/'
				+ growth_condition + '_data_' + current_mode + '_' + current_fixed_parameter
				+ '_' + str(current_sim) + '.csv')
			if os.path.isfile(current_MLE_output):
				with open(current_MLE_output, 'rU') as current_MLE_output_contents:
					current_MLE_output_data = list(csv.reader(current_MLE_output_contents))

				current_MLE_output_np = numpy.array(map(float,current_MLE_output_data[0]))
					# warning - doesn't return list in python 3

				sim_output_array[current_sim-1,:] = current_MLE_output_np


		# sort sim_output_array by the LL and save
		sim_output_array_sorted = sim_output_array[sim_output_array[:,0].argsort()]
		numpy.sort(sim_output_array_sorted)
		numpy.savetxt(current_sim_combo_file, sim_output_array_sorted, delimiter=",")

		# figure out percentile of true_LL
		all_nan_rows = numpy.all(numpy.isnan(sim_output_array),axis=1)
		sim_output_array_clean = sim_output_array[~all_nan_rows,:]

		sim_output_array_clean_shape = numpy.shape(sim_output_array_clean)
		total_successful_sim_rows = sim_output_array_clean_shape[0]
		sim_rows_below_true_LL = sum(sim_output_array_clean[:,0]<=true_LL)
		true_LL_p_val = float(sim_rows_below_true_LL)/total_successful_sim_rows
		percentile_summary_matrix[0][1] = true_LL_p_val
		numpy.savetxt(percentile_summary_file, percentile_summary_matrix, delimiter=",")

		open(current_ML_combination_completefile,'a').close()		

		############### ??? TO DO ??? ###############
		# at some point need to get rid of illegitimate values (ones where parameter has bumped into limits)
		############### ??? TO DO ??? ###############

	return(current_ML_combination_completefile)

##############################################################################
# main

setup_file = Check_Input()
	# setup_file is given as input when running the code, and contains
		# info about directories containing code and data

(username,home_dir,composite_data_dir,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_) = \
	Read_Setup_File(setup_file)
	# get info needed to loop through folders; the rest of setup info
		# is read from the setup_file in each individual folder, in
		# case the general setup_file changes as an experiment is
		# processing

# path of file that determines whether MLE_runner already running
currently_running_checkfile = '/' + home_dir + '/' + username + '/MLE_running.txt'

# Loop through all directories in composite_data_dir
for output_dir_name in os.walk(composite_data_dir).next()[1]:

	print(output_dir_name)

	current_folder = composite_data_dir + '/' + output_dir_name
	complete_checkfile = current_folder + '/processing_complete.txt'

	# Only run MLE_runner if it's not already running and if this
		# folder hasn't been processed yet
	if not os.path.isfile(currently_running_checkfile) and \
		not os.path.isfile(complete_checkfile):

		# If setup_file doesn't exist in current directory, copy it to
			# current_folder
		local_setup_file = current_folder+'/setup_file.csv'
		if not os.path.isfile(local_setup_file):
			copyfile(setup_file, local_setup_file)

		# Read setup_file from current_folder
		(_,_,_,backup_dir,code_dir,temp_storage_dir,
			starting_L,starting_gridpower,sim_repeats,profile_points_by_mode,MLE_memory,
			MLE_time,user_email,strain_number,mode_list,parameters_by_mode,
			parameter_max_by_mode,parameter_min_by_mode,parameter_start_by_mode,
			pdf_points,parallel_processors_max,lambda_SNM,
			max_memory,max_time,ms_positions,parameter_profile_ub_by_mode,
			parameter_profile_lb_by_mode) = Read_Setup_File(local_setup_file)
		# get general info necessary to run the rest of code

		# create MLE_running.txt so no new instances of MLE_runner.py run
		open(currently_running_checkfile,'a').close()

		# Get names of necessary folders, and, if this hasn't been done
			# before, create those folders

	#	(trackfile_folder,completefile_folder,slurm_folder,sim_output_folder,
	#		MLE_output_folder,LL_profile_folder,epilogue_folder,
	#		MLE_combined_sim_folder) = \
		(trackfile_folder,completefile_folder,slurm_folder,sim_output_folder,
			MLE_output_folder,LL_profile_folder,MLE_combined_sim_folder) = \
			Setup_Directories(current_folder,temp_storage_dir,output_dir_name,username)

		##########
		# Run data analysis

		try:
			phenotype_file_list = subprocess.check_output('ls -lrt '
				+ current_folder + '/' + '*_phenotype_file.csv',shell=True)
		except subprocess.CalledProcessError:
			phenotype_file_list = ''
		growth_condition_list = re.findall(current_folder + '/(.+?)_phenotype_file.csv',
			phenotype_file_list,re.MULTILINE)
#		growth_condition_list = ['SC']
		############### ??? TO DO ??? ###############
		# somehow get growth_condition_list, possibly as output of data analysis?
#		rep_list = ['1','2']
		rep_list = ['1']
		# somehow create starting_parameter_val_file
		# as output of data analysis, create phenotype_file
		############### ??? TO DO ??? ###############

		# Run all steps for one phenotype (growth condition) at a time,
			# but loop through modes separately for each step
		for current_growth_condition in growth_condition_list:

			############### ??? TO DO ??? ###############
			# Check for data analysis complete_file before running anything in this loop
			############### ??? TO DO ??? ###############

			current_phenotype_file = current_folder + '/' \
				+ current_growth_condition + '_phenotype_file.csv'
			current_petite_file = current_folder + '/' \
				+ current_growth_condition + '_petite_GR_data.csv'
				# contains 4 columns, with data for each individual
					# measured strain in each row.
					# 1: observed phenotype;
					# 2: s.e. of observed phenotype;
					# 3. known SNP number in strain;
					# 4. proportion of genome sequenced to satisfactory
						# depth to make mutation calls
		
			##########
			# Run MLE - fixed and unfixed
			# Run 2 reps - one 'fast' with standard grid_power and L,
				# and a second one to ensure same results when
				# grid_power and L are increased

			for rep_index, rep in enumerate(rep_list):
				rep_float = float(rep)
				current_L = starting_L*pow(1.5,rep_float-1)
					# returns starting_L for rep1 and starting_L*1.5 for rep2
				current_gridpower = starting_gridpower+(rep_float-1)
				current_MLE_time = MLE_time*rep_float
				current_MLE_memory = MLE_memory*rep_float

				current_MLE_output_subfolder = MLE_output_folder + '/rep_' + rep
				if not os.path.isdir(current_MLE_output_subfolder):
					os.makedirs(current_MLE_output_subfolder)

				for current_mode_index, current_mode in enumerate(mode_list):

					print('##################################################################################')
					print(current_mode)

					# Use parameter lists appropriate to current_mode, but
						# adding 'unfixed' parameter regime to the list
					current_parameter_list = parameters_by_mode[current_mode_index]
					current_parameter_max_list = numpy.array(
						parameter_max_by_mode[current_mode_index])
					current_parameter_min_list = numpy.array(
							parameter_min_by_mode[current_mode_index])
					current_start_values = numpy.array(
							parameter_start_by_mode[current_mode_index])
					current_parameter_profile_ub_list = numpy.array(
						parameter_profile_ub_by_mode[current_mode_index])
					current_parameter_profile_lb_list = numpy.array(
							parameter_profile_lb_by_mode[current_mode_index])

					# find the total number of parameters, including fixed ones
					total_parameter_number = len(current_parameter_list)

					current_profile_points_list = numpy.array(
						profile_points_by_mode[current_mode_index])
					if len(current_profile_points_list)==1:
						current_profile_points_list = \
							numpy.repeat(current_profile_points_list,total_parameter_number)

					# Determine which parameters are always fixed
						# (the ones where min and max values are the same)
					pernamently_fixed_parameter_bool = \
						current_parameter_max_list==current_parameter_min_list

					parameters_to_loop_over_bool = numpy.invert(pernamently_fixed_parameter_bool)* \
						numpy.invert((current_profile_points_list < 1))

					non_permafixed_parameters = [item for (item,bool_val) in \
						zip(current_parameter_list,parameters_to_loop_over_bool) \
						if bool_val]
					parameters_to_loop_over = ['unfixed'] + non_permafixed_parameters
					point_numbers_to_loop_over = numpy.append([1], \
						current_profile_points_list[ \
							numpy.invert(pernamently_fixed_parameter_bool)])

	#				starting_parameter_val_file = (current_folder + '/'
	#					+ current_growth_condition + '_' + current_mode
	#					+ '_MLE_starting_vals.csv')

					MLE_summary_file = (current_folder + '/' + current_mode + '_'
						+ current_growth_condition + '_' + rep + '_MLE_file.csv')

					MLE_completefile_list = []

					for current_fixed_parameter in parameters_to_loop_over:

						current_fixed_parameter_bool_list = numpy.empty_like(pernamently_fixed_parameter_bool)
						numpy.copyto(current_fixed_parameter_bool_list,pernamently_fixed_parameter_bool)


						if current_fixed_parameter == 'unfixed':
							current_profile_point_num = 1
						else:
							# find index of current_fixed_parameter in
								# parameter list, and add it to bool
								# list of fixed parameters
							current_fixed_parameter_index = current_parameter_list.index(current_fixed_parameter)
							current_fixed_parameter_bool_list[current_fixed_parameter_index] = True
							current_profile_point_num = current_profile_points_list[current_fixed_parameter_index]

						#current_fixed_param_min = \
						#	current_parameter_min_list[current_fixed_parameter_index]
						#current_fixed_param_max = \
						#	current_parameter_max_list[current_fixed_parameter_index]

						print(current_fixed_parameter)

						current_MLE_completefile = ML_Estimator(
							completefile_folder,trackfile_folder,
							current_profile_point_num,current_mode,mode_list,username,
							current_MLE_time,current_MLE_memory,slurm_folder,
							current_fixed_parameter,
							current_fixed_parameter_bool_list,
							current_MLE_output_subfolder,
							current_growth_condition,current_L,current_gridpower,
							current_phenotype_file,current_petite_file,rep,strain_number,
							parallel_processors_max,output_dir_name,max_memory,
							max_time,current_start_values,
							current_parameter_max_list,
							current_parameter_min_list,
							pernamently_fixed_parameter_bool,
							current_parameter_list,total_parameter_number,
							ms_positions,current_parameter_profile_ub_list,
							current_parameter_profile_lb_list)

						MLE_completefile_list.append(current_MLE_completefile)

						print('done')


					# Check that MLE is complete for all parameters in this folder
					MLE_completeness_test_list = []
					for MLE_completefile_to_test in MLE_completefile_list:
						if os.path.isfile(MLE_completefile_to_test):
							MLE_completeness_test_list.append(1)
					print('completeness done')

					if sum(MLE_completeness_test_list) == len(parameters_to_loop_over):
					#if sum(MLE_completeness_test_list) >= 24:
						# run MLE on every non-fixed parameter and in unfixed mode
						
						fixed_parameter_list = parameters_to_loop_over[1:]
						#fixed_parameter_list = parameters_to_loop_over[1:24]
							# doesn't include 'unfixed'

						current_ML_combination_completefile = MLE_combiner(
							completefile_folder,current_growth_condition,
							current_mode,rep,current_folder,MLE_summary_file,
							point_numbers_to_loop_over,fixed_parameter_list,
							current_MLE_output_subfolder,LL_profile_folder,lambda_SNM,current_parameter_list)







#						# plot pdf for multiple mutations
#						use_single_mutation = False
#						if not current_mode == 'known_muts':							
#							current_pdf_completefile = pdf_Plotter(completefile_folder,trackfile_folder,
#								pdf_points,current_mode,username,
#								slurm_folder,current_folder,current_growth_condition,
#								current_L,current_gridpower,MLE_summary_file,
#								current_phenotype_file,rep,lambda_SNM,epilogue_folder,output_dir_name,max_memory,max_time,use_single_mutation)
#						# plot pdf for single mutational effect
#						use_single_mutation = True
#						single_mut_L = 'NaN'
#						single_mut_gridpower = current_gridpower+8
#						single_mut_pdf_points = pdf_points*5
#						current_pdf_completefile = pdf_Plotter(completefile_folder,trackfile_folder,
#							single_mut_pdf_points,current_mode,username,
#							slurm_folder,current_folder,current_growth_condition,
#							single_mut_L,single_mut_gridpower,MLE_summary_file,
#							current_phenotype_file,rep,lambda_SNM,epilogue_folder,output_dir_name,max_memory,max_time,use_single_mutation)
#
#						if current_mode == 'mixedgauss':
#							# also plot 'combined' single pdf whose
#								# poisson mixture is equivalent to
#								# mixedgauss distribution
#							use_single_mutation = True
#							single_mut_L = 'NaN'
#							single_mut_gridpower = current_gridpower+8
#							single_mut_pdf_points = pdf_points*5
#							current_pdf_completefile = pdf_Plotter(completefile_folder,trackfile_folder,
#								single_mut_pdf_points,'combo_poissonmut_mixedgauss',username,
#								slurm_folder,current_folder,current_growth_condition,
#								single_mut_L,single_mut_gridpower,MLE_summary_file,
#								current_phenotype_file,rep,lambda_SNM,epilogue_folder,output_dir_name,max_memory,max_time,use_single_mutation)
#						elif current_mode == 'mixed':
#							# also plot 'combined' single pdf whose
#								# poisson mixture is equivalent to
#								# mixedgauss distribution
#							use_single_mutation = True
#							single_mut_L = 'NaN'
#							single_mut_gridpower = current_gridpower+8
#							single_mut_pdf_points = pdf_points*5
#							current_pdf_completefile = pdf_Plotter(completefile_folder,trackfile_folder,
#								single_mut_pdf_points,'combo_poissonmut_mixed',username,
#								slurm_folder,current_folder,current_growth_condition,
#								single_mut_L,single_mut_gridpower,MLE_summary_file,
#								current_phenotype_file,rep,lambda_SNM,epilogue_folder,output_dir_name,max_memory,max_time,use_single_mutation)
#
#			#########
#			# Compare MLE results between the two reps?
#
#			#########
#			# Run LR: unfixed (and fixed?) AFTER MLE COMPLETE FOR ALL PARAMETERS
#
#			# Run Simulation to create simulated phenotype files
#			rep = 'sim'
#			# doesn't make sense to run simulation on knownmuts mode?
#			mode_list_random_muts = mode_list[:]
#			mode_list_random_muts.remove('knownmuts')
#
#			parameters_by_mode_random_muts = parameters_by_mode[:]
#			del parameters_by_mode_random_muts[mode_list.index('knownmuts')]
#
#			for current_mode_index, current_mode in enumerate(mode_list_random_muts):
#
#				# Determine which parameters are always fixed
#					# (the ones where min and max values are the same)
#				current_parameter_max_list = numpy.array(parameter_max_by_mode[current_mode_index])
#				current_parameter_min_list = numpy.array(parameter_min_by_mode[current_mode_index])
#				pernamently_fixed_parameter_bool = \
#					current_parameter_max_list==current_parameter_min_list
#
#				fixed_parameter_list = parameters_by_mode_random_muts[current_mode_index]
#
#				MLE_summary_file = (current_folder + '/' + current_mode + '_'
#							+ current_growth_condition + '_1_MLE_file.csv')
#
#				if os.path.isfile(MLE_summary_file):
#
#					simulation_setup_completefile = Simulation_Setup(sim_repeats,current_mode,
#						username,output_dir_name,slurm_folder,user_email,strain_number,
#						MLE_summary_file,current_growth_condition,completefile_folder,trackfile_folder,
#						sim_output_folder,current_phenotype_file,epilogue_folder,max_memory,
#						max_time,lambda_SNM)
#
#		#			for current_mode in mode_list:
#
#		#				Simulation_Setup(sim_repeats,current_mode,
#		#					mode_list,username,output_dir_name,slurm_folder,user_email,
#		#					strain_number,MLE_output_file,current_growth_condition,
#		#					completefile_folder,trackfile_folder,sim_output_folder)
#					if os.path.isfile(simulation_setup_completefile):
#						current_fixed_parameter = 'unfixed'
#						current_fixed_param_min = 0;
#						current_fixed_param_max = 0;
#
#						current_MLE_output_subfolder = MLE_output_folder + '/rep_' + rep
#						if not os.path.isdir(current_MLE_output_subfolder):
#							os.makedirs(current_MLE_output_subfolder)
#
#						starting_parameter_val_file = MLE_summary_file
#							# start parameter search at MLE values being simulated
#
#						current_MLE_phenotype_file = sim_output_folder + '/phenotype_table_'\
#							 + current_mode + '_' + current_growth_condition + '_'
#							# ML_Estimator will add id# and .csv
#
#						sim_MLE_completefile = ML_Estimator(
#										completefile_folder,trackfile_folder,
#										sim_repeats,current_mode,mode_list,username,
#										MLE_time,MLE_memory,slurm_folder,
#										current_fixed_parameter,current_fixed_param_min,
#										current_fixed_param_max,current_MLE_output_subfolder,
#										current_growth_condition,starting_L,
#										starting_gridpower,starting_parameter_val_file,
#										current_MLE_phenotype_file,rep,strain_number,
#										parallel_processors_max,epilogue_folder,output_dir_name,max_memory,max_time,pernamently_fixed_parameter_bool)
#
#						if os.path.isfile(sim_MLE_completefile):
#
#							sim_process_completefile = Simulation_Processor(
#								completefile_folder,current_mode,rep,
#								current_growth_condition,MLE_summary_file,sim_repeats,
#								fixed_parameter_list,MLE_combined_sim_folder,
#								current_MLE_output_subfolder,current_fixed_parameter)

						############### ??? TO DO ??? ###############
						# Currently, Simulation_MLE only runs an 'unfixed' MLE on simulated data - correct if necessary
						# Simulation_MLE is NOT ACTUALLY FINISHED
						############### ??? TO DO ??? ###############
		#				Simulation_MLE()




				

			# Compile LR results and calculate p-value

############### ??? TO DO ??? ###############
# Remove temporary folder on scratch?
# Remove epilogue folder
# Back up current_folder on archive
# !!! write complete_checkfile !!!
############### ??? TO DO ??? ###############

		# remove MLE_running.txt so MLE_runner.py can run again
		os.remove(currently_running_checkfile)




