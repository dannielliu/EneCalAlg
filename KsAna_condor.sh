#!/bin/bash
#datadir="./data"
datadir="$datadir/xyz4600_Ksto2pi_3/xyz4600_Ksto2pi_3"
outdir="$PWD/jobs"

for data in `ls $datadir/*.root`;do
  ene=`echo $data | awk -F "." '{print $1}' | awk -F "/" '{print $NF}'` # | awk -F "_" '{print $NF}'
  mkdir -p "$outdir/$ene"
  cat > ${outdir}/$ene/${ene}.cmd <<EOF
  Universe             = vanilla
  Notification         = Never
  GetEnv               = True
  Executable           = $PWD/run.sh
  Arguments            = $data ${outdir}/$ene
  Output               = $outdir/$ene/${ene}.out
  Error                = $outdir/$ene/${ene}.err
  Log                  = $outdir/$ene/${ene}.log
  +Group               = "BESIII"
  #should_transfer_files= yes
  requirements         = (substr(Machine,0,4)!="bl-0"&&ARCH=="X86_64")&& (machine != "bl-3-15.hep.ustc.edu.cn") && (machine != "bl-3-16.hep.ustc.edu.cn") && (machine != "bl-3-06.hep.ustc.edu.cn") &&(machine != "bl-2-15.hep.ustc.edu.cn")
  WhenToTransferOutput = ON_EXIT_OR_EVICT 
  OnExitRemove         = TRUE 
  Queue
EOF

done

 for cmd in `ls $outdir`;do
   echo $cmd
   if [ -d $cmd ];then
     echo $outdir/$cmd/$cmd.cmd
     condor_submit $outdir/$cmd/$cmd.cmd
   fi
 done
