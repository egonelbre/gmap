#!/bin/bash

mkdir -p _test_
cd _test_

DATASET="../data/test"
PHENO="../data/pheno.txt.pset"
LOG=../test_gmap.log
BIN=../../bin

rm -f $LOG

echo test started `date` >> $LOG
echo >> $LOG

echo `date`
echo

function timeit(){
  echo `date +"%T"` : $*
  echo -n "$*" >> $LOG
  time ( $* &> /dev/null ) 2>> $LOG
  echo 1>> $LOG
}

echo -e "gmap converting\n" >> $LOG

timeit $BIN/gmapconvert -i $DATASET test
timeit $BIN/gmappack -i test ptest

echo -e "gmap (unpacked)\n" >> $LOG

timeit $BIN/gmapfreq test
timeit $BIN/gmaphardyweinberg test
timeit $BIN/gmapassoc -g $PHENO test

echo -e "gmap (packed)\n" >> $LOG

timeit $BIN/gmapfreq ptest
timeit $BIN/gmaphardyweinberg ptest
timeit $BIN/gmapassoc -g $PHENO ptest
