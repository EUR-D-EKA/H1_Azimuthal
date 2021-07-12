#!/bin/bash
#
# ------------------------------------------------ #
#
# Submit the entire analysis to the bird HTC cluster
# 
# 1) generate a new directory to 
#    + preserve the code/steerings
#    + hold the output files
# 2) run all jobs, by executing at a bird node: run_HTC.sh
#
# ------------------------------------------------ #

# Check if the input is valid
if [[ -z $2 ]]; then
    echo "Usage: "
    echo "$0 <jobname> <Steering File>"
    echo ""
    echo "Example:"
    echo "submit_HTC.sh jobname name.steer"
    exit 0
fi

# Check the steering file
STEERFILE=$2
if [[ -f $STEERFILE ]] ; then
    echo "OK! Steering file found. It is:"
    head -2 $STEERFILE
    echo "[...]"
    tail -2 $STEERFILE
else
    echo "ERROR! Steering file NOT found!"
    exit 1
fi

# making job directory
JOBNAME=$1-`date +%Y-%m-%d-%H:%M`
echo $JOBNAME
if [[ -d $JOBNAME ]]; then
    echo "Output directory for the given job name already exists. Exiting"
    exit 0
fi
mkdir $JOBNAME


# compile executable and library
# make all
# status=$?
# if [ $status -eq 0 ]; then
#     echo "Compilation completed."
# else
#     echo ""
#     echo "ERROR! make failed. Please check source code and run $0 again."
#     exit 1;
# fi


# # copy all files to job directory
# echo "Now copying files to work directory $JOBNAME"
# ln -s ../LumiFiles .
# mkdir -p $JOBNAME/lib
# cp *.C *.h *.steer $JOBNAME/.
# cp -r ../lib/$PLATFORM/ $JOBNAME/lib/.
# cp ../bin/$PLATFORM/EventShapes $JOBNAME/.
# cp -r Steering $JOBNAME/.
# cp QUEUE* $JOBNAME/.
cp run_HTC.sh $JOBNAME/.

# mkdir -p $JOBNAME/output
# mkdir -p $JOBNAME/log

# go to job directory and make CONDOR submit file:
cd $JOBNAME
cat <<EOF > HTC.condor
executable               = run_HTC.sh
transfer_executable    = True
universe               = vanilla
arguments                = -f ../$2
output                   = job.out
error                    = job.error
log                      = job.log
RequestMemory            = 2048
# should_transfer_files  = Yes
getenv                 = True
when_to_transfer_output  = ON_EXIT
requirements             = OpSysAndVer == "CentOS7"
# +RequestRuntime        = 3500
+RequestRuntime        = 86399
# +RequestRuntime        = 604800
# +RequestRuntime        = 259200
# +MyProject             = "h1"
queue 1
EOF

echo "Submitting jobs for runperiod: all"
condor_submit HTC.condor 
echo "Check job status with '$ condor_q'"
echo "Check if option -g (gen level only) was selected in run_HTC.sh"
