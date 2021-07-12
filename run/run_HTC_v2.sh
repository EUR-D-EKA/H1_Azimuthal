#!/bin/zsh

echo "Running EventShape analysis"
echo "Repeat execution, by running:"
echo " ------------------------------------------------------------------------ "
echo "run_HTC.sh  $@"
echo " ------------------------------------------------------------------------ "
echo "  STEERING:   $1"
echo "  CHAIN:      $2"
echo "  OUTPUTDIR:  $3  (optional - default: 'output')"
echo "  PWD:        $4  (optional - default '.')"
echo " ------------------------------------------------------------------------ "
echo " "
echo "Setting up environment"

echo " ------------------------------------------------------------------------ "
echo "This node:"
uname -a
echo "HOST=$HOST"
echo "HOSTNAME=$HOSTNAME"
echo " ------------------------------------------------------------------------ "
echo " ------------------------------------------------------------------------ "

THISDIR=$PWD
cd $THISDIR

OUTDIR=output

echo ""
echo "Now running analysis"
../../bin/x86_64-centos7-gcc9-opt/AzimuthalAnalyzer -f $1 -o $2.root
