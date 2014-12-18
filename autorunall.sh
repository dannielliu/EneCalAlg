#!/bin/bash
#datadir="./data"
#datadir="/Volumes/data/besIIIwork/"
#mv parpipill.txt parpipill.txt.old
#mv analog analog.old
#mv detail.txt detail.txt.old
make clean
make
./autoanapipill.sh
cp ./src/gepep_kpi.C.vkpi ./src/gepep_kpi.C
make
./autoanakpi.sh
cp ./src/gepep_kpi.C.vkpipi ./src/gepep_kpi.C
make
./autoanakpipi.sh
root -l weightedpar.C
#./autoanaf4pi.sh
./autoanaf6pi.sh
echo "done"



