#!/bin/bash
#datadir="./data"
datadir="/Volumes/data/besIIIwork/"
#mv parpipill.txt parpipill.txt.old
#mv analog analog.old
#mv detail.txt detail.txt.old

for algdir in `ls $datadir`;do
  if test -d "$datadir/$algdir";then
    if [ $algdir == "fastpipill" ];then
	  for datafile in `ls ${datadir}/${algdir}`;do
	    ene=`echo $datafile | awk -F "." '{print $1}' | awk -F "_" '{print $NF}'`
		mkdir -p "./graphs/$ene"
		echo $ene>>"parpipill.txt"
		echo $ene>>"detail.txt"
		echo "analysis $ene"
        ./analysis ${datadir}/${algdir}/${datafile} "./graphs/$ene" >> analog
	  done
    fi
  fi
done