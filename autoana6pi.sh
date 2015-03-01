#!/bin/bash
#datadir="./data"
datadir="/Volumes/data/besIIIwork/"
#mv par.txt par.txt.old
#mv analog analog.old
#mv detail.txt detail.txt.old

for algdir in `ls $datadir`;do
  if test -d "$datadir/$algdir";then
    if [ $algdir == "fast6pi" ];then
	  for datafile in `ls ${datadir}/${algdir}/*.root`;do
	    ene=`echo $datafile | awk -F "." '{print $1}' | awk -F "_" '{print $NF}'`
		mkdir -p "./graphs/$ene"
#		echo $ene>>"par.txt"
#		echo $ene>>"detail.txt"
        echo "analysis $ene ..."
        ./analysis ${datafile} "./graphs/$ene" >> analog
	  done
    fi
  fi
done
