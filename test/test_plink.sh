#!/bin/bash

mkdir -p _test_
cd _test_

PLINK="../plink"
DATASET="../data/test"
PHENO="../data/pheno.txt"
LOG=../test_plink.log

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

echo -e "plink converting\n" >> $LOG

timeit $PLINK --file $DATASET --out bedtest --make-bed
timeit $PLINK --file $DATASET --out tpedtest --recode --transpose

echo -e "plink (ped)\n" >> $LOG

timeit $PLINK --file $DATASET --hardy2 --nonfounders
timeit $PLINK --file $DATASET --freq
timeit $PLINK --file $DATASET --assoc --pheno $PHENO --nonfounders --allow-no-sex

echo -e "plink (tped)\n" >> $LOG

timeit $PLINK --tfile tpedtest --hardy2 --nonfounders
timeit $PLINK --tfile tpedtest --freq
timeit $PLINK --tfile tpedtest --assoc --pheno pheno.txt --nonfounders --allow-no-sex

echo -e "plink (bed)\n" >> $LOG

timeit $PLINK --bfile bedtest --hardy2 --nonfounders
timeit $PLINK --bfile bedtest --freq
timeit $PLINK --bfile bedtest --assoc --pheno pheno.txt --nonfounders --allow-no-sex
