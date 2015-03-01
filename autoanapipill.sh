#!/bin/bash
#datadir="./data"
datadir="/Volumes/data/besIIIwork"
#datadir="/Volumes/data2"
parfile="parpipill"
#mv parkk.txt parkk.txt.old
#mv analog analog.old
#mv detail.txt detail.txt.old
echo -e "energy\t\tfactor pi\terror">$parfile

for algdir in `ls $datadir`;do
  if test -d "$datadir/$algdir";then
    #if [ $algdir == "CombinedRootPipillXYZ4260" ];then
    if [ $algdir == "Rvalue_fastpipill" ];then
      for datafile in `ls ${datadir}/${algdir}`;do
        ene=`echo $datafile | awk -F "." '{print $1}' | awk -F "_" '{print $NF}'`
        mkdir -p "./graphs2/$ene"
        echo $ene>>"parpipill.txt"
        echo $ene>>"detailpipill.txt"
        echo "analysis $ene"
        ./analysis ${datadir}/${algdir}/${datafile} "./graphs2/$ene" >> anapipilllog
        echo -n -e "$ene\t">>$parfile
        cat par>>$parfile
        echo  >>$parfile
      done
    fi
  fi
done
