#!/bin/bash
DCM_DIR=$1
OUT_DIR=`mktemp -d`
if [[ "$#" -eq 2 ]]; then
  OUT_DIR=$2
fi
mkdir $OUT_DIR

echo "Converting..."
dcm2niix -o $OUT_DIR -x n -f t1 $DCM_DIR

echo "Skull strip..."
runROBEX.sh ${OUT_DIR}/t1.nii ${OUT_DIR}/t1_brain.nii

echo "Normalization...."
flirt -in ${OUT_DIR}/t1_brain \
  -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain \
  -out ${OUT_DIR}/t1_brain_mni \
  -omat ${OUT_DIR}/t1_brain_mni.mat \
  -bins 256 -cost corratio \
  -searchrx -90 90 -searchry -90 90 -searchrz -90 90 \
  -dof 12 -interp trilinear

for voi in /siemens_auto_voxel/vois/*.mat; do
  voi_name=`basename $voi`
  voi_name=${voi_name%.mat}
  echo "======== $voi_name ========" > ${OUT_DIR}/voi_${voi_name}.txt
  /siemens_auto_voxel/siemens_online_voxel.py ${OUT_DIR}/t1.nii $voi ${OUT_DIR}/t1_brain_mni.mat ${OUT_DIR}/voi_${voi_name} >> ${OUT_DIR}/voi_${voi_name}.txt
  echo "================" >> ${OUT_DIR}/voi_${voi_name}.txt
#  echo ""

done
cat ${OUT_DIR}/voi_*.txt

fslmerge -t ${OUT_DIR}/vois ${OUT_DIR}/voi_*.nii.gz
fslmaths ${OUT_DIR}/vois -Tmax ${OUT_DIR}/vois_combined
overlay 1 1 ${OUT_DIR}/t1 -a ${OUT_DIR}/vois_combined .5 1 ${OUT_DIR}/vois_rendered
slicer ${OUT_DIR}/vois_rendered -L -l /data1/fsl/etc/luts/renderhot.lut -c -S 5 2048 ${OUT_DIR}/vois_rendered.png




