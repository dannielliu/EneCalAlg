#!/bin/bash
datadir="."
echo $#
if (( $# >= 1 ));then
  datadir=$1
fi
echo "data directory is $datadir"
for dir in `ls -d $datadir/*`;do
  echo "in $dir"
  if [ -d $dir ];then
    #echo "ana $dir"
    #echo "./AnaSinglePart ${dir}/parspur.txt $dir"
    suffix=`echo $dir | awk -F "/" '{print $NF}' | awk -F "_" '{print $NF}'`
    echo $suffix
    mv ${dir}/Ks_factors.pdf $dir/../Ks_factors_part$suffix.pdf
  fi
done
