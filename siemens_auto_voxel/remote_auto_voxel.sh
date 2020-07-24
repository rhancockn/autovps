#!/bin/bash
SUBJ=$1
SESSION=$2
SERIES=$3
VOXEL=$4
cd ~/auto_voxels
outdir=${VOXEL}-${SUBJ}-${SESSION}-${SERIES}
mkdir $outdir
ssh -t auto_voxel@psypacs.psy.uconn.edu "source /usr/local/anaconda3/bin/activate && /data1/auto_voxel/scripts/auto_voxel.sh $SUBJ $SESSION $SERIES $VOXEL" |tee $outdir/log
scp -r auto_voxel@psypacs.psy.uconn.edu:/tmp/${outdir} .
open ./${outdir}/vois_rendered.png
cd ./${outdir}
cat voi_*.txt > coordinates.txt

cat coordinates.txt | lpr
cat coordinates.txt


ssh -t auto_voxel@psypacs.psy.uconn.edu "rm -rf /tmp/${outdir}"

