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
#for datadir in `ls -d /Volumes/data/besIIIwork/uidata/liud/Rvalue_Ks_2`;do
for datadir in `ls -d /Volumes/data2/newFiles`;do
  #dataset=`echo $datadir | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}'`
  dataset="Rvalue"
  outdir="$PWD/Ks_${dataset}_cut2p_20range1_0"
  dirc_exist $outdir
  mkdir -p $outdir

  for data in `ls $datadir/Rvalue_Ks_all.root`;do
    echo $data
    ene=`echo $data | awk -F "." '{print $1}' | awk -F "/" '{print $NF}'` # | awk -F "_" '{print $NF}'
    mkdir -p "$outdir/$ene"
    echo "../analysis $data $outdir/$ene > $outdir/$ene/log"
    ../analysis $data $outdir/$ene > $outdir/$ene/log 
  done
done
date
