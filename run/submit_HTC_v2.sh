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
    echo "$0 <jobname> <QUEUE-file>"
    echo ""
    echo "Example:"
    echo "submit_HTC.sh jobname QUEUE_all.txt"
    exit 0
fi

# Check the steering file
QUEUEFILE=$2
if [[ -f $QUEUEFILE ]] ; then
    echo "OK! QUEUE-file found. QUEUE-file is:"
    head -2 $QUEUEFILE
    echo "[...]"
    tail -2 $QUEUEFILE
else
    echo "ERROR! QUEUE-file NOT found!"
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

cp run_HTC_v2.sh $JOBNAME/.
cp Queue.txt $JOBNAME/.
cp SelectedRuns.root FidVolCut.steer $JOBNAME/.

# go to job directory and make CONDOR submit file:
cd $JOBNAME
cat <<EOF > HTC.condor
executable               = run_HTC.sh
transfer_executable      = True
universe                 = vanilla
arguments                = \$(STEERING) \$(CHAIN) \$(OUT) ${PWD}
output                   = \$(CHAIN).out
error                    = \$(CHAIN).error
log                      = \$(CHAIN).log
RequestMemory            = 2048
# should_transfer_files    = Yes
getenv                   = True
when_to_transfer_output  = ON_EXIT
requirements             = OpSysAndVer == "CentOS7"
+RequestRuntime          = 3500
# +RequestRuntime         = 86399
#+RequestRuntime         = 604800
#+RequestRuntime         = 259200
queue STEERING CHAIN OUT from ${QUEUEFILE}
EOF

echo "Submitting jobs for runperiod: all"
condor_submit HTC.condor 
echo "Check job status with '$ condor_q'"
echo "Check if option -g (gen level only) was selected in run_HTC.sh"
