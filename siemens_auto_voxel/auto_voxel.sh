#!/bin/bash
source /usr/local/anaconda3/bin/activate
export FSLDIR=/data1/fsl
export PATH=$PATH:/data1/fsl/bin
source ${FSLDIR}/etc/fslconf/fsl.sh
export FSLOUTPUTTYPE=NIFTI
export AUTOVOXELDIR=/data1/auto_voxel
export PATH=$PATH:$AUTOVOXELDIR/bin
export PATH=$PATH:$AUTOVOXELDIR/bin/ROBEX
export NIDBDIR=/data1/nidb/archive

SUBJ=$1
SESSION=$2
SERIES=$3
VOXEL=$4

#tdir=`mktemp -d`
tdir=/tmp/${VOXEL}-${SUBJ}-${SESSION}-${SERIES}
mkdir $tdir

echo "Converting..."
dcm2niix -o $tdir -x n -f t1 ${NIDBDIR}/${SUBJ}/${SESSION}/${SERIES}

echo "Skull strip..."
runROBEX.sh ${tdir}/t1.nii ${tdir}/t1_brain.nii

echo "Normalization...."
flirt -in ${tdir}/t1_brain -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -out ${tdir}/t1_brain_mni -omat ${tdir}/t1_brain_mni.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear

for voi in ${AUTOVOXELDIR}/voxels/${VOXEL}/*.mat; do
  voi_name=`basename $voi`
  voi_name=${voi_name%.mat}
  echo "======== $voi_name ========" > ${tdir}/voi_${voi_name}.txt
  $AUTOVOXELDIR/scripts/siemens_online_voxel.py ${tdir}/t1.nii $voi ${tdir}/t1_brain_mni.mat ${tdir}/voi_${voi_name} >> ${tdir}/voi_${voi_name}.txt
  echo "================" >> ${tdir}/voi_${voi_name}.txt
#  echo ""

done

fslmerge -t ${tdir}/vois ${tdir}/voi_*
fslmaths ${tdir}/vois -Tmax ${tdir}/vois_combined
overlay 1 1 ${tdir}/t1 -a ${tdir}/vois_combined .5 1 ${tdir}/vois_rendered
slicer ${tdir}/vois_rendered -L -l /data1/fsl/etc/luts/renderhot.lut -c -S 5 2048 ${tdir}/vois_rendered.png

cat ${tdir}/voi_*.txt

echo $tdir
#rm -rf $tdir


