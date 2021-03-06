#!/usr/bin/env python

#puts MNI voxel into T1 and Siemens space
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import numpy as np
import os
import nibabel as nib
import scipy.spatial as spatial
from autovps.transform import Transform
from autovps.dataset import siemens
from autovps.transform import Transform
from autovps.dataset import siemens

import argparse


parser = argparse.ArgumentParser(description="Generate ROI voxels from MRS affine and T1 files.")
parser.add_argument('t1_nifti', help="Path to T1 nifti, converted from dicoms.")
parser.add_argument('svs_transform', help="A file containing the affine transform of the MRS voxel")
parser.add_argument('native2mni', help="A file containing the native to mni transform")

parser.add_argument('roi_prefix', help="Prefix for output files")

args = parser.parse_args()

# create a labelled position string
def position_string(pos):
    directions_pos = ['R', 'A', 'H']
    directions_neg = ['L', 'P', 'F']
    direction_string = []
    for i in range(3):
        if pos[i] < 0:
            d = directions_neg[i]
        else:
            d = directions_pos[i]

        direction_string.append('%s%0.0f' % (d, np.abs(pos[i])))

    return(' '.join(direction_string))


# load the MNI template
FSLDIR = os.getenv('FSLDIR', '/usr/local/fsl')
mni_nifti = nib.load(os.path.join(FSLDIR, 'data/standard/MNI152_T1_1mm.nii.gz'))
mni_aff = Transform(mni_nifti.get_qform())

# load T1
t1_nifti = nib.load(args.t1_nifti)
t1_aff = t1_nifti.get_qform()
t1_inv = np.linalg.pinv(t1_aff)

# load native2mni FLIRT transform
native2mni = np.loadtxt(args.native2mni)

# load the template voxel spec, in mm
mrs_aff_orig = np.loadtxt(args.svs_transform)
mrs_mm_tform = Transform(mrs_aff_orig)
# get scaling factor of native2mni transform and scale the voxel to preserve size
# this is mm, so do before voxel conversion

mrs_mm_tform.scale(Transform(native2mni).get_scale())
# change to voxel coordinates
#mrs_aff_vox = np.dot(np.linalg.pinv(mni_aff), t_mrs_aff_orig.get_matrix())
mrs_aff_vox = mni_aff.get_inverse() * mrs_mm_tform


## Do the transform

#see https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;4b44c252.0907

# flip matrix for reference (T1)
Wref = np.eye(4)
if np.linalg.det(t1_aff) > 0:
    Wref[0,0]=-1 #TODO: make Transform subscriptable
    Wref[0,3]=t1_nifti.get_shape()[0]-1

#scale matrix for reference (T1)
t1_dims = [np.sqrt(np.dot(t1_aff[:, i].T.tolist(), t1_aff[:, i].tolist())) for i in range(3)]
Sref = np.diag([t1_dims[0], t1_dims[1], t1_dims[2], 1])

#scale matrix for input (MNI)
Sin = np.eye(4)

# flip matrix for input (MNI)
Win = np.eye(4)
Win[0,0]=-1
Win[0,3]=mni_nifti.get_shape()[0]-1

composed_affine = np.dot(Wref, np.linalg.pinv(Sref))
composed_affine = np.dot(composed_affine, np.linalg.pinv(native2mni))
composed_affine = np.dot(composed_affine, Sin)
composed_affine = np.dot(composed_affine, Win)

#the template voxel in T1 space (vox coordinates)
composed_affine = np.dot(composed_affine, mrs_aff_vox.get_matrix())
#print(composed_affine)
#template voxel in mm (scanner) coordinates
composed_affine_mm=np.array(np.dot(t1_aff, composed_affine))

vox_tform = Transform(composed_affine_mm)

ori = vox_tform.siemens_orientation()
pos = vox_tform.get_position()
# vox_size = [np.sqrt(np.dot(composed_affine_mm[:, i].T.tolist(), composed_affine_mm[:, i].tolist())) for i in range(3)]
# use the original affine to avoid rounding errors
vox_size = [np.sqrt(np.dot(mrs_aff_orig[:, i].T.tolist(), mrs_aff_orig[:, i].tolist())) for i in range(3)]

ip_rot = -1.0*float(ori[1])

if ip_rot > 90.0:
    ip_rot = 180.0 - ip_rot
if ip_rot < -90.0:
    ip_rot = 180 + ip_rot

print('Orientation: %s' % ori[0])
print('Rotation: %0.1f deg' %  ip_rot)
print('Position: %0.0f %0.0f %0.0f mm' % (-1*pos[0], -1*pos[1], pos[2]))
print('          %s' % position_string(pos))
print('VOI: R>>L%0.0f A>>P%0.0f F>>H%0.0f mm' % (vox_size[0], vox_size[1], vox_size[2]))

# ### make a volume in T1 space
mrs_corners = [[-0.5, -0.5, -0.5, 1], #0
               [-0.5, -0.5,  0.5, 1], #1
               [-0.5,  0.5, -0.5, 1], #2
               [-0.5,  0.5,  0.5, 1], #3
               [ 0.5, -0.5, -0.5, 1], #4
               [ 0.5, -0.5,  0.5, 1], #5
               [ 0.5,  0.5, -0.5, 1], #6
               [ 0.5,  0.5,  0.5, 1]] #7


#don't round off here
t1_corners = np.array([(np.dot(composed_affine, c)) for c in mrs_corners]).squeeze()

# resample to match voxel grid
#for i in range(3):
#    t1_corners[:,i] = t1_corners[:,i]/t1_dims[i]

t1_data = t1_nifti.get_data().squeeze()
mrs_roi = np.ones_like(t1_data) * 0

#exhaustive search for points in the voxel
tri = spatial.Delaunay(t1_corners[:,0:3])
for i in range(np.floor(min(t1_corners[:,0])).astype(int),np.ceil(max(t1_corners[:,0])).astype(int)):
    for j in range(np.floor(min(t1_corners[:,1])).astype(int),np.ceil(max(t1_corners[:,1])).astype(int)):
        for k in range(np.floor(min(t1_corners[:,2])).astype(int),np.ceil(max(t1_corners[:,2])).astype(int)):
            mrs_roi[i,j,k] = tri.find_simplex([i,j,k]) >=0

img=nib.Nifti1Image(mrs_roi, t1_aff)
nib.save(img, args.roi_prefix + '.nii.gz')
