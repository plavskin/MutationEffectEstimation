----------------------------
# Mutation Effect Estimation
----------------------------
This repository contains code for estimating mutational effects from growth rate data for individual yeast colonies

## Components and pipeline
--------------------------
The code to perform analysis of mutational effects from colony growth data consists of the following components, some of which are still under development:
1. [ ] Image analysis and colony tracking code **'PIE'**, described in a [preprint by Li, Plavskin *et al.*](http://dx.doi.org/10.1101/253724).
1. [ ] Code to process the PIE output into data required by the mutation effect maximum likelihood estimation pipeline. This includes:
    1. For petite strains, a *.csv* file containing growth rates for each colony and information on the strain identities and the subexperiment in which each datapoint was collected.
    1. For strains of interest (i.e. MA strains and reference controls), a *.csv* file containing difference within random pairs of colonies of the reference strain and a strain of interest, and information on the strain identities and the subexperiment in which each datapoint was collected.
1. [x] Code to run maximum likelihood estimation (MLE) pipeline. This includes:
    * The [cluster_wrangler package](cluster_wrangler), for assigning, running, and tracking batch jobs on a cluster (currently runs on SLURM and, with highly limited functionality, on a Mac)
    * The [mlestimator package](mlestimator) for running MLE components, including initial estimation, log-likelihood profile estimation, data simulation, and confidence interval estimation from log-likelihood profiles-based and simulation-based likelihood ratio tests.
    * [Test code](gaussian_mle) for implementing the maximum likelihood estimation pipeline for numbers drawn from a normal distribution.
1. [ ] Code to idenitify the maximum likelihood parameter values for growth rate data, including ML estimates for each strain's individual effect, and for the overall distribution(s) of strain mutational effects

## Implementation
-----------------
While mlestimator and cluster_wrangler are written in python 2.7, all maximum likelihood estimation functions themselves run in matlab, and require the Optimization toolkit.
Estimating mutational effects, and especially confidence intervals on mutational effect parameters, requires extensive computational time and is best run on a high-performance cluster (currently, only clusters running SLURM are supported by cluster_wrangler). However, a quick test of this system can be performed using the included [gaussian_mle](gaussian_mle) code, which allows mlestimator to be run on data containing normally-distributed random numbers, and can be run on a Mac via commandline.
To run maximum likelhood estimation using the test data:
1. Select the appropriate setup file (setup_file_gauss_slurm.csv if running code on SLURM, setup_file_gauss_mac.csv if running code on Mac) and edit the rows corresponding to **username**, **home_folder**, **composite_data_folder**, **backup_folder**, **code_folder**, and **temp_storage_folder**. **composite_data_folder** should be the path to the directory that *contains* [gaussian_test_folder](data/gaussian_test_folder).
1. *Optional:* change the data in **gaussian_test_folder/SC_phenotype_file.csv**; this is the data that the maximum likelihood parameters will be estimated for, and it can contain any number of integers or floats under 'data'.
1. Run *MLE_runner.py*, selecting the appropriate setup_file under the '-s' option; e.g.:
    ```
    MLE_runner.py -s setup_file_gauss_mac.csv
    ```
    On a cluster running SLURM, running *MLE_runner.py* will result in the submission of a batch job; on a Mac computer, it currently results in the necessary jobs being submitted to matlab consecutively through the active terminal session.
1. Run *MLE_runner.py* again as above; this will:
    1. Update the status of the current jobs in **gaussian_test_folder/trackfiles/**; these can be most conveniently viewed in the *summary* files in this folder.
    1. Submit any jobs not submitted in the first run
    1. Resubmit any jobs that ran out of time/memory in the previous run (on SLURM)
1. Run *MLE_runner.py* several times after all jobs have been completed to compile log likelihood profiles (in **gaussian_test_folder/LLprofiles/**), identify confidence interval bounds (this will also submit jobs running matlab scripts), and write and update an **MLE_output** file in **gaussian_test_folder** that will contain the maximum likelihood estimates and confidence intervals for each parameter, as well as info on the runtime for each maximum likelihood estimation.
	* Note: it may be useful to set MLE_runner.py running on *cron*; for this, modify the MLE_runner_shell.sh file as necessary, and run this file on the cron.



