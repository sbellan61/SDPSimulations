#!/bin/bash
#
# Simple SLURM script for submitting multiple serial
# jobs (e.g. parametric studies) using a script wrapper
# to launch the jobs.
#
# To use, build the launcher executable and your
# serial application(s) and place them in your WORKDIR
# directory.  Then, edit the CONTROL_FILE to specify 
# each executable per process.
#-------------------------------------------------------
#-------------------------------------------------------
# 
#------------------Scheduler Options--------------------
#SBATCH -J Disc	          # Job name
#SBATCH -N 1                 # Total number of nodes (24 cores/node)
#SBATCH -n 48                # Total number of tasks
#SBATCH -p normal          # Queue name
#SBATCH -o Disco.o%j      # Name of stdout output file (%j expands to jobid)
#SBATCH -t 24:00:00            # Run time (hh:mm:ss)
#SBATCH --mail-user=steve.bellan@gmail.com
#SBATCH --mail-type=all  # email me all messages
#      <------------ Account String ------------>
# <--- (Use this ONLY if you have MULTIPLE accounts) --->
#------------------------------------------------------
#
# Usage:
#	#$ -pe <parallel environment> <number of slots> 
#	#$ -l h_rt=hours:minutes:seconds to specify run time limit
# 	#$ -N <job name>
# 	#$ -q <queue name>
# 	#$ -o <job output file>
#	   NOTE: The env variable $JOB_ID contains the job id. 
#
#------------------------------------------------------

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE="HetCounterFactualAcute.txt"

ml launcher
$LAUNCHER_DIR/paramrun
