"""
Allows tracking of job submission to computer or server, including
tracking error status of jobs that did not run successfully.
Currently supports batch job submission on SLURM or one-by-one
submission of jobs on MacOSX.
Takes an input setup_file.csv (see example), which can be processed by
the cluster_functions.parse_setup_file() function; all job submission
and job status verification is handled by cluster_functions.job_flow_handler
"""