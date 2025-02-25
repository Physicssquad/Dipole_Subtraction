#!/bin/bash

############# Bash Submission Script for c5pc00  ############
### Amol Rahane: Tue Jul 20 18:40:18 IST 2010
##               Wed Jul 28 09:57:00 IST 2010
#################################################

############ Enter no. of nodes (N) and total no. of processors (n) required (Maximum 4 Processors per Node). ##########
####### Do Not Remove comment(#)
#SBATCH  -N 5 -n 20

############ Export working directory (Do not modify) ###########
export work_dir=`pwd`
cd $work_dir
echo $work_dir

############ create mpd ring, DO NOT modify this section unless you really know ############
srun hostname >hosts
mpdboot -n $SLURM_NNODES -v -f hosts

############ ONLY change the path of the executable (program) to your own application and necessary arguments  #########

##### VASP Gamma Point Calculations
mpiexec -l -machinefile hosts -n $SLURM_NPROCS $work_dir/cteq.x 

##### VASP BAND Calculations
##mpiexec -l -machinefile hosts -n $SLURM_NPROCS /c5scratch/mrinal/bin/vasp.4621.band.par.x

############ exit the mpd ring and clean off the nodes (Do not change) ###################

mpdallexit
mpdcleanup

`rm hosts`

exit

#################################################
######## Submit this script using "sbatch script.sh"
