source /home/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64

# Recompiling mare2dem (and redirecting the output to a file to avoid garbage)
cd ../mare2dem-omeep-newest/
make clean >> compile.out 2<&1
make INCLUDE=./include/omeep.inc INTEL_PATH=/home/intel TRACING=1 >> compile.out 2<&1

#dumping old/untracked scorep directories into scorep-dmp                                                                                                                               
#mv scorep-* dmp-scorep/

# SET YOUR PARAMETERS
# Set scorep config file
m2d_dir=~/mare2dem-omeep-newest
cd $m2d_dir/Examples/configScorep
set -a
. scorep.cfg
set +a
# Need to execute from the input directory
bin_app="/home/users/fpmjunior/mare2dem-omeep-newest/Examples/mare2dem_examples_v4.1/Demo/inversion_CSEM"
cd $bin_app
app=$m2d_dir/MARE2DEM
arg=Demo.0.resistivity

timestamp=`date +"%Y%m%d_%I%M"`
app_short_hash=`git rev-parse --short HEAD`

exp="EXP_${app_short_hash}_${timestamp}" # TODO nome dos experimentos
mkdir -p "$exp"
filedat="METADATA_${exp}.org"

#echo "$exp"
#exit 1

running="mpirun -n 34 $app $arg"
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

# CREATE TRACING PLOT
cd $exp/scorep-${exp}

~/otf2utils/otf22csv traces.otf2 | grep -E 'occamCompute|compute|FourierDeriv' > traces.csv

Rscript $m2d_dir/Examples/configScorep/trace.r 'traces1.pdf'

~/otf2utils/otf22csv traces.otf2 | grep -E 'jacobianCompute|smoothingOccam|compute|FourierDeriv' > traces.csv

Rscript $m2d_dir/Examples/configScorep/trace.r 'traces2.pdf'

~/otf2utils/otf22csv traces.otf2 | grep -E 'occamCompute|primal_solve|error_estimate|derivs_comp_adj|lhs_gen|em2dkx_localRefinement|FourierDeriv' > traces.csv

Rscript $m2d_dir/Examples/configScorep/trace.r 'traces3.pdf'

~/otf2utils/otf22csv traces.otf2 | grep -E 'jacobianCompute|smoothingOccam|primal_solve|error_estimate|derivs_comp_adj|lhs_gen|em2dkx_localRefinement|FourierDeriv' > traces.csv

Rscript $m2d_dir/Examples/configScorep/trace.r 'traces4.pdf'

