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

cat <<EOF > HTC.condor
executable               = run_HTC_mini.sh
transfer_executable      = True
universe                 = vanilla
output                   = mini.out
error                    = mini.error
log                      = mini.log
RequestMemory            = 64000
# should_transfer_files    = Yes
getenv                   = True
when_to_transfer_output  = ON_EXIT
requirements             = OpSysAndVer == "CentOS7"
# +RequestRuntime          = 3500
+RequestRuntime         = 86399
#+RequestRuntime         = 604800
#+RequestRuntime         = 259200
queue 1
EOF

echo "Submitting jobs for runperiod: all"
condor_submit HTC.condor 
echo "Check job status with '$ condor_q'"
echo "Check if option -g (gen level only) was selected in run_HTC.sh"
