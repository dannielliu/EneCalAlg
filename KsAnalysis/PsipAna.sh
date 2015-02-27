#!/bin/bash
#datadir="./data"
#dataset="xyz4230"
#echo "there are $# input arguements"
#if (( $# >= 1 ));then
#  dataset=$1
#fi
#datadir="/Volumes/data/besIIIwork/uidata/liud/${dataset}_Ks_2"



dirc_exist(){
  if [ -d $1 ] ; then
     echo "The directory : $1 exist, please remove it"
     exit 1
  fi
}

date
#for datadir in `ls -d /Volumes/data/besIIIwork/Rvalue_fpipill5/combinedroot`;do
for datadir in `ls -d /Volumes/data2/Rvalue`;do
  dataset=`echo $datadir | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}'`
  
  outdir="$PWD/Psi_${dataset}_cut1p_20range"
  dirc_exist $outdir
  mkdir -p $outdir

  for data in `ls $datadir/*_fpipill_*.root`;do
    echo $data
    ene=`echo $data | awk -F "." '{print $1}' | awk -F "/" '{print $NF}'` # | awk -F "_" '{print $NF}'
    mkdir -p "$outdir/$ene"
    echo "../checkf $data $outdir/$ene > $outdir/$ene/log"
    ../analysis $data $outdir/$ene > $outdir/$ene/log 
  date
  done
done
date
