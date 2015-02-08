#!/bin/bash
#datadir="./data"
datadir="$datadir/xyz4260_Ks_2/combinedroot"
outdir="$PWD/Ks_4260"
mkdir -p $outdir

for data in `ls $datadir/*.root`;do
  ene=`echo $data | awk -F "." '{print $1}' | awk -F "/" '{print $NF}'` # | awk -F "_" '{print $NF}'
  mkdir -p "$outdir/$ene"
  ./analysis $data $outdir/$ene > $outdir/$ene/log &
done

 
