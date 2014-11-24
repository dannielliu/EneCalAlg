#!/bin/bash
#datadir="./data"
datadir="/Volumes/data/besIIIwork/"
parfile="parkpi"
mv parkpi.txt parkpi.txt.old
#mv analog analog.old
mv detailkpi.txt detailkpi.txt.old
echo -e "energy\t\tfactor pi\terror">$parfile

for algdir in `ls $datadir`;do
  if test -d "$datadir/$algdir";then
    if [ $algdir == "Rvalue_kpi" ];then
      for datafile in `ls ${datadir}/${algdir}`;do
        ene=`echo $datafile | awk -F "." '{print $1}' | awk -F "_" '{print $NF}'`
        mkdir -p "./graphs/$ene"
        echo $ene>>"parkpi.txt"
        echo $ene>>"detailkpi.txt"
        echo "analysis $ene"
        ./analysis ${datadir}/${algdir}/${datafile} "./graphs/$ene" >> anakpilog
        echo -n -e "$ene\t\t">>$parfile
        cat par>>$parfile
        echo  >>$parfile
        done
    fi
  fi
done
