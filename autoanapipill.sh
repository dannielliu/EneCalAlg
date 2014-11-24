#!/bin/bash
#datadir="./data"
datadir="/Volumes/data/besIIIwork/"
parfile="parpipill"
mv parpipille.txt parpipille.txt.old
mv parpipillmu.txt parpipillmu.txt.old
mv parpipillpi.txt parpipillpi.txt.old
#mv analog analog.old
#mv detail.txt detail.txt.old
echo -e "energy\t\tfactor e\terror\t\tfactor mu\terror\t\tfactor pi\terror">$parfile

for algdir in `ls $datadir`;do
  if test -d "$datadir/$algdir";then
    if [ $algdir == "fastpipill" ];then
      for datafile in `ls ${datadir}/${algdir}`;do
        ene=`echo $datafile | awk -F "." '{print $1}' | awk -F "_" '{print $NF}'`
        mkdir -p "./graphs/$ene"
        echo $ene>>"parpipille.txt"
        echo $ene>>"parpipillmu.txt"
        echo $ene>>"parpipillpi.txt"
        echo $ene>>"detail.txt"
        echo "analysis $ene"
        ./analysis ${datadir}/${algdir}/${datafile} "./graphs/$ene" >> analog
        echo -n -e "$ene\t\t">>$parfile
        cat par>>$parfile
        echo  >>$parfile
        done
    fi
  fi
done

