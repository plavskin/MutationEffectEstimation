----------------------------
# ML parameter estimation code
----------------------------
This repository contains general code for running maximal likelihood parameter estimation on a cluster. It contains python code for managing cluster job submissions, as well as matlab code to handle the MLE.

## Components and pipeline
--------------------------
This reporistory contains code to run maximum likelihood estimation (MLE) pipeline. This includes:

    * The [cluster_wrangler package](cluster_wrangler), for assigning, running, and tracking batch jobs on a cluster (currently runs on SLURM and, with highly limited functionality, on a personal computer with a unix-based OS such as a Mac)
    * The [mlestimator package](mlestimator) for running MLE components, including initial estimation, log-likelihood profile estimation, data simulation, and confidence interval estimation from log-likelihood profiles-based and simulation-based likelihood ratio tests.
    * [Test code](exponential_fam_LL_and_sim) for implementing the maximum likelihood estimation pipeline for numbers drawn from a normal, exponential, or gamma distribution.

## Function
-----------
In short, this repository can be used for ML parameter estimation as long as a likelihood function is provided as matlab code; the inputs of the code should include a vector of parameter values, and two map objects: `input_value_dict`, which can contain parameters specified in the setup file (although these can also be modified by upstream code, e.g. names of parameters, for an example see [setup_file_gauss_unix.csv](setup_file_gauss_unix.csv)), and `pre_MLE_output_dict`, which should contain any objects set at the beginning of the run (e.g. the loaded in object that contains the data).

The optimization steps themselves are run through the [`MLE_finder`](mle_finder/MLE_finder.m) function; the optimization can optionally include constraint functions, functions that specify the gradient of the likelihood, functions that pre- or post-process data and/or results, and functions that simulate data from a set of parameters to allow for simulation-based confidence interval estimation.

## Implementation
-----------------
While mlestimator and cluster_wrangler are written in python 2.7, all maximum likelihood estimation functions themselves run in matlab, and require the Optimization toolkit.

Estimating mutational effects, and especially confidence intervals on mutational effect parameters, requires extensive computational time and is best run on a high-performance cluster (currently, only clusters running SLURM are supported by cluster_wrangler). However, a quick test of this system can be performed using the included [exponential_fam_LL_and_sim](exponential_fam_LL_and_sim) code, which allows mlestimator to be run on data containing normally-, exponentially-, or gamma-distributed random numbers, and can be run on a unix-based OS via commandline.

To run maximum likelhood estimation using the test data:

1. Select the appropriate setup file (setup_file_gauss_slurm.csv if running code on SLURM, setup_file_gauss_unix.csv if running code on unix-based OS) and edit the rows corresponding to **username**, **home_folder**, **composite_data_folder**, **backup_folder**, **code_folder**, and **temp_storage_folder**. **composite_data_folder** should be the path to the directory that *contains* [gaussian_test_folder](data/gaussian_test_folder).
1. *Optional:* change the data in **gaussian_test_folder/phenotype_file.csv**; this is the data that the maximum likelihood parameters will be estimated for, and it can contain any number of integers or floats under 'data'.
1. Run *MLE_runner.py*, selecting the appropriate setup_file under the '-s' option; e.g.:
    ```
    MLE_runner.py -s setup_file_gauss_unix.csv
    ```
    On a cluster running SLURM, running *MLE_runner.py* will result in the submission of a batch job; on a personal computer with a unix-based OS, it results in the necessary jobs being submitted to matlab consecutively through the active terminal session.
1. Run *MLE_runner.py* again as above; this will:
    1. Update the status of the current jobs in **gaussian_test_folder/trackfiles/**; these can be most conveniently viewed in the *summary* files in this folder.
    1. Submit any jobs not submitted in the first run
    1. Resubmit any jobs that ran out of time/memory in the previous run (on SLURM)
1. Run *MLE_runner.py* several times after all jobs have been completed to compile log likelihood profiles (in **gaussian_test_folder/LLprofiles/**), identify confidence interval bounds (this will also submit jobs running matlab scripts), and write and update an **MLE_output** file in **gaussian_test_folder** that will contain the maximum likelihood estimates and confidence intervals for each parameter, as well as info on the runtime for each maximum likelihood estimation.
	* Note: it may be useful to set MLE_runner.py running on *cron*; for this, modify the MLE_runner_shell.sh file as necessary, and run this file on the cron.



