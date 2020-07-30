from autovps import make_voi
import nibabel as nib
import numpy as np
from autovps.transform import Transform

#from method
##$PVM_VoxelGeoCub=( 1 )
(((0.994521895375187 - 0.104528463201876 0
0.104528463201876 0.994521895375187 0
0 0 1, 55 - 32 - 8), 20 30 25,
# Multiply translation by [-1 1 -1]

R = np.array([[0.994521895375187, -0.104528463201876, 0, -55],
                [0.104528463201876, 0.994521895375187, 0, -32],
                [0, 0, 1, 8],
                [0,0,0,1]])

tform = Transform(R)
tform.scale([20, 30, 25])

t1_path = '214016/dicom/BRAIN/11_HASKI/31/1/Brain--02_AxialScout_NS70--RM--2031617.nii.gz'

t1 = nib.load(t1_path)
img = make_voi.make_voi(t1, tform)
nib.save(img, 'mask.nii.gz')




# convert to dicoms
bruker2dicom convert 20180428_092608_Haskins_1_11_214016/ dicom/

# align axial t2 to t1
/usr/local/fsl/bin/flirt \
-in 214016/dicom/BRAIN/11_HASKI/31/1/Brain--02_AxialScout_NS70--RM--2031617.nii.gz \
-ref /Users/rhancock/scratch/a214/sub-pa0579_ses-01_T1w.nii.gz \
-out /Users/rhancock/scratch/a214/t2_to_t1 \
-omat /Users/rhancock/scratch/a214/t2_to_t1.mat \
-bins 256 -cost normcorr -searchrx -90 90 \
-searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear


# apply
/usr/local/fsl/bin/flirt -in /Users/rhancock/scratch/a214/mask.nii.gz -applyxfm -init /Users/rhancock/scratch/a214/t2_to_t1.mat -out /Users/rhancock/scratch/a214/mask_t1.nii.
