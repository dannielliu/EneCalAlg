#!/bin/bash
#datadir="./data"
datadir="/Volumes/data/besIIIwork/"
parfile="parf6pi"
mv parf6pi.txt parf6pi.txt.old
#mv analog analog.old
#mv detailkpi.txt detailkpi.txt.old
echo -e "energy\t\tpeak_nofit\terror\tpeak_fit\terror">$parfile

for algdir in `ls $datadir`;do
  if test -d "$datadir/$algdir";then
    if [ $algdir == "fast6pi" ];then
      for datafile in `ls ${datadir}/${algdir}`;do
        ene=`echo $datafile | awk -F "." '{print $1}' | awk -F "_" '{print $NF}'`
        mkdir -p "./graphs/$ene"
        echo $ene>>"parf6pi.txt"
        echo $ene>>"detailf6pi.txt"
        echo "analysis $ene"
        ./analysis ${datadir}/${algdir}/${datafile} "./graphs/$ene" >> anaf6pilog
        echo -n -e "$ene\t\t">>$parfile
        cat par>>$parfile
        echo  >>$parfile
      done
    fi
  fi
done
