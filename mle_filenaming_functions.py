#!/usr/bin/python

# Contains functions necessary for naming files in MLE

import os

def generate_filename(output_file_path, profile_point_as_str, \
	output_id, parameter_name, output_prename):
	# creates a filename to which output file of MLE will be written,
		# to be read by LLProfile
	output_file_label = generate_file_label(output_prename, output_id, \
		parameter_name)
	output_file = '_'.join([output_file_label, profile_point_as_str]) + '.csv'
	output_filename = os.path.join(output_file_path,output_file)
	return(output_filename)

def generate_file_label(output_prename, output_id, parameter_name):
	# creates a file label that can be used by Cluster_Functions to
		# track the completeness of the file
	# the file label doesn't include the path, profile point, or file
		# extension
	output_file_label = '_'.join([output_prename, output_id, parameter_name])
	return(output_file_label)