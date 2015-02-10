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
    echo "ana $dir"
    echo "./AnaSinglePart ${dir}/parspur.txt $dir"
    ./AnaSinglePart ${dir}/parspur.txt $dir
  fi
done
