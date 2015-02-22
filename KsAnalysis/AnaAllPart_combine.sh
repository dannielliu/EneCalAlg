#!/bin/bash
datadir="."
filelist=""
echo $#
if (( $# >= 1 ));then
  datadir=$1
fi
echo "data directory is $datadir"
for dir in `ls -d $datadir/*`;do
  echo "in $dir"
  if [ -d $dir ];then
    filelist="$filelist ${dir}/parspur.txt"
    #echo "ana $dir"
    #echo "./AnaSinglePart ${dir}/parspur.txt $dir"
    #./AnaSinglePart ${dir}/parspur.txt $dir
  fi
done
../AnaSomePart $filelist
mv combinedpar.txt $datadir
echo "moving Ks_factors_sum.pdf to $datadir"
mv Ks_factors_sum.pdf $datadir
