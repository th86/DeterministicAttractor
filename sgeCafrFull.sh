#!/bin/sh
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -t 1-100
export INPUTFILE=$1
export JOBNAME=$2
Rscript parallel.sgeFull.R $INPUTFILE $SGE_TASK_ID $SGE_TASK_LAST $JOBNAME
