#!/bin/zsh

# echo "Running EventShape analysis"
echo "Repeat execution, by running:"
echo " ------------------------------------------------------------------------ "
echo "run_HTC.sh $@"
echo " ------------------------------------------------------------------------ "
echo " "
echo "Setting up environment"

echo " ------------------------------------------------------------------------ "
echo "This node:"
uname -a
echo "HOST=$HOST"
echo "HOSTNAME=$HOSTNAME"
echo " ------------------------------------------------------------------------ "

THISDIR=$PWD

echo ""
echo $THISDIR
echo "Now running analysis"
../../bin/x86_64-centos7-gcc9-opt/AzimuthalAnalyzer $1 $2

