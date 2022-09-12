#!/bin/bash
############################## SET YOUR PARAMETERS!!!!! ##################################################
############################## SET YOUR PARAMETERS!!!!! ##################################################
############################## SET YOUR PARAMETERS!!!!! ##################################################
# Set the number of MPI tasks you want to run
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=tsubasa
#SBATCH --time=1:00:00
#SBATCH --exclusive
#SBATCH --job-name=mare2dem-job
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

############################## SET YOUR PARAMETERS!!!!! ##################################################
############################## SET YOUR PARAMETERS!!!!! ##################################################
############################## SET YOUR PARAMETERS!!!!! ################################################## 
# IT MAY NOT WORK WITHOUT ABSOLUTE PATHS
git_repo="/home/users/dancarastan/mare2dem-omeep-felix/mare2dem-omeep-newest"
#input_dir="/home/users/dancarastan/mare2dem-omeep-newest/Examples/mare2dem_examples_v4.1/Demo/inversion_CSEM/"
input_dir="/home/users/dancarastan/mare2dem-omeep-felix/mare2dem-omeep-newest/Examples/Linha_04Tx13/"
MACHINEFILE="/home/users/dancarastan/mare2dem-omeep-felix/mare2dem-omeep-newest/nodes.$SLURM_JOB_ID"
run_args=04Tx013.0.resistivity
source /home/users/dancarastan/.bashrc

export SCOREP_ENABLE_PROFILING=false
export SCOREP_ENABLE_TRACING=true
export SCOREP_VERBOSE=false
export SCOREP_TOTAL_MEMORY=4G
export SCOREP_PAGE_SIZE=8K
export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=100
export SCOREP_MPI_ENABLE_GROUPS=all
#export M2D_NPROCS=$((SLURM_NTASKS+1))
export M2D_NPROCS=$((SLURM_NTASKS+0))

#export SCOREP_SELECTIVE_CONFIG_FILE=/home/users/dancarastan/mare2dem-omeep-newest/configScorep/regionsSelected.cfg
##########################################################################################################

#MACHINEFILE="/home/users/dancarastan/mare2dem-omeep-newest/nodes.$SLURM_JOB_ID"
srun -l hostname | sort -n | awk '{print $2}' > $MACHINEFILE

source /home/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64

cd $git_repo
make clean > compile.out 2> compile.err
make INCLUDE=./include/omeep_nec.inc NEC_PATH=/opt/nec TRACING=0 ITERATION=1 >> compile.out 2>> compile.err
#. ./configScorep/scorep.cfg 
#make INCLUDE=./include/omeep.inc PREP=scorep INTEL_PATH=/home/intel/ >> compile.out 2>> compile.err

#dumping old/untracked scorep directories into scorep-dmp                                                                                                                               
#mv scorep-* dmp-scorep/

# Need to execute from the input directory
bin_app=$input_dir
cd $bin_app
m2d_dir=$git_repo
app=$m2d_dir/MARE2DEM
arg=$run_args

timestamp=`date +"%Y%m%d_%H%M"`
app_short_hash=`git rev-parse --short HEAD`

exp="EXP_${SLURM_JOB_ID}_${timestamp}_${app_short_hash}"
mkdir -p "$exp"
filedat="METADATA_${exp}.org"

#echo "$exp"
#exit 1

running="mpirun -n $M2D_NPROCS -machinefile $MACHINEFILE $app $arg"
# SLURM_NTASKS is an environment variable set by slurm accordingly with the number of MPI tasks
#running="mpirun -n $M2D_NPROCS $app $arg"
#running="ls -la"

$m2d_dir/machine_info.sh $filedat

echo "############################################" >> $filedat
echo "* APPLICATION INFO:" >> $filedat

app_hash=`git show --quiet HEAD`
echo "** REPOSITORY INFORMATION:" >> $filedat
echo "$app_hash" >> $filedat 

echo "############################################" >> $filedat
echo "** PROGRAM BINARY LIBRARIES:" >> $filedat
ldd $app >> $filedat

echo "############################################" >> $filedat
echo "Execute..." >> $filedat

# Writing results
rm -f stdout.out
rm -f stderr.out

# Executing and taking the time
time1=$(date +%s.%N)
set +e # In order to detect and print execution errors
eval $running 1> stdout.out 2> stderr.out
set -e
time2=$(date +%s.%N)
echo "* ELAPSED TIME:" >> $filedat
echo "Elapsed:    $(echo "$time2 - $time1"|bc ) seconds" >> $filedat

# Checking for errors
echo "* STDERR OUTPUT:" >> $filedat
cat stderr.out >> $filedat
if [ ! -s stdout.out ]; then
      echo "ERROR DURING THE EXECUTION!!!" >> stdout.out
fi

echo "* STDOUT OUTPUT:" >> $filedat
(echo -n "Makespan (in ms): "; cat stdout.out) >> $filedat
if [[ $verbose == 1 ]]; then
      cat stderr.out
          echo -n "Makespan (in ms): "
              cat stdout.out
fi
rm -f stdout.out
rm -f stderr.out

#moving scorep directory to the exp directory
scorep_out_dir=`ls $bin_app | grep scorep`
#mv ${bin_app}/scorep-* ${exp}/scorep-${exp}
mv $scorep_out_dir ${exp}/scorep-${exp}

#mv $m2d_dir/compile.out $exp/
mv $filedat $exp/
