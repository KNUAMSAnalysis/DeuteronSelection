#!/bin/bash

# PURPOSE: This script helps to submit massive amount of jobs
# AUTHOR: Wooyoung Jang
# DATE: 2014. 12. 22

date -u  # print current time

USRNAME=`whoami`
CMDEXE="/afs/cern.ch/work/${USRNAME:0:1}/$USRNAME/Pass6/bin/main"
EXEC="/afs/cern.ch/work/${USRNAME:0:1}/$USRNAME/Pass6/exec.sh"
OUTDIR="/afs/cern.ch/work/${USRNAME:0:1}/$USRNAME/Pass6/Output"
LOGDIR=$OUTDIR/$1/log
ERRDIR=$OUTDIR/$1/err
LISTDIR="/afs/cern.ch/work/${USRNAME:0:1}/$USRNAME/Pass6/FileList"

if [ -e $CMDEXE ]; then
  echo "\$CMDEXE found .. $CMDEXE"
else
  echo "Can not find $CMDEXE"
fi

if [ -e $EXEC ]; then
  echo "\$EXEC found ... $EXEC"
else
  echo "Can not find $EXEC"
fi

if [ -d $OUTDIR ]; then
  echo "Output files will be located at $OUTDIR"
else
  mkdir $OUTDIR
  echo "Warning: $OUTDIR is made to store output files."
fi

if [ -d $OUTDIR/$1 ]; then
  echo "Output files will be located at subdirectory $OUTDIR/$1"
else
  mkdir $OUTDIR/$1
  echo "Warning: $OUTDIR/$1 subdirectory is made to store output files."
fi

if [ -d $LOGDIR ]; then
  echo "Log files will be stored at $LOGDIR"
else
  mkdir $LOGDIR
  echo "Warning: $LOGDIR is made to store log files."
fi

if [ -d $ERRDIR ]; then
  echo "Error files will be stored at $ERRDIR"
else
  mkdir $ERRDIR
  echo "Warning: $ERRDIR is made to store error files."
fi

ls -1 $LISTDIR/$1 > templist_$1

while read irun
do
  echo "/usr/bin/bsub -J $irun -o $LOGDIR/$irun.log -e $ERRDIR/$irun.err -q 8nh -n 1 $EXEC $LISTDIR/$irun $OUTDIR/$1/$irun.root"
  /usr/bin/bsub -J $irun -o $LOGDIR/$irun.log -e $ERRDIR/$irun.err -q 8nh -n 1 $EXEC $LISTDIR/$1/$irun $OUTDIR/$1/$irun.root
  sleep 5
done < "templist_$1"

# Delete the temporary file list
rm templist_$1
echo "templist_$1 is deleted."
