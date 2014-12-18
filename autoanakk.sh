#!/bin/bash
#datadir="./data"
datadir="/Volumes/data/besIIIwork/"
parfile="parkk"
mv parkk.txt parkk.txt.old
#mv analog analog.old
#mv detail.txt detail.txt.old
echo -e "energy\t\tfactor k\terror">$parfile

for algdir in `ls $datadir`;do
  if test -d "$datadir/$algdir";then
    if [ $algdir == "Rvalue_kk" ];then
      for datafile in `ls ${datadir}/${algdir}`;do
        ene=`echo $datafile | awk -F "." '{print $1}' | awk -F "_" '{print $NF}'`
        mkdir -p "./graphs/$ene"
        echo $ene>>"parkk.txt"
        echo $ene>>"detailkk.txt"
        echo "analysis $ene"
        ./analysis ${datadir}/${algdir}/${datafile} "./graphs/$ene" >> anakklog
        echo -n -e "$ene\t\t">>$parfile
        cat par>>$parfile
        echo  >>$parfile
      done
    fi
  fi
done
