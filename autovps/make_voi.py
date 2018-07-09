import numpy as np
from autovps.transform import Transform
import nibabel as nib
import scipy.spatial as spatial

def make_voi(t1_nifti, voi_tform):
	""" Creates a volumetric ROI for SVS data.

	Parameters:
		t1_nifti (Nifti1Image): A Nifti1Image with an affine transform
		voi_tform (Transform): A Transform with the affine for the voxel

	Returns:
		img (Nifti1Image): A binary image of the voxel in T1 space

	"""
	
	t1_aff = Transform(t1_nifti.get_qform())
	t1_aff_inv = t1_aff.get_inverse()
	composed_affine = t1_aff_inv * voi_tform


	mrs_corners = [[-0.5, -0.5, -0.5, 1], #0
               [-0.5, -0.5,  0.5, 1], #1
               [-0.5,  0.5, -0.5, 1], #2
               [-0.5,  0.5,  0.5, 1], #3
               [ 0.5, -0.5, -0.5, 1], #4
               [ 0.5, -0.5,  0.5, 1], #5
               [ 0.5,  0.5, -0.5, 1], #6
               [ 0.5,  0.5,  0.5, 1]] #7

	#don't round off here
	t1_corners = np.array([(np.dot(composed_affine.get_matrix(), c)) for c in mrs_corners]).squeeze()
	t1_data = t1_nifti.get_data().squeeze()
	mrs_roi = np.ones_like(t1_data) * 0

	#exhaustive search for points in the voxel
	tri = spatial.Delaunay(t1_corners[:,0:3])
	for i in range(np.floor(min(t1_corners[:,0])).astype(int),np.ceil(max(t1_corners[:,0])).astype(int)):
	    for j in range(np.floor(min(t1_corners[:,1])).astype(int),np.ceil(max(t1_corners[:,1])).astype(int)):
	        for k in range(np.floor(min(t1_corners[:,2])).astype(int),np.ceil(max(t1_corners[:,2])).astype(int)):
	            mrs_roi[i,j,k] = tri.find_simplex([i,j,k]) >=0

	img = nib.Nifti1Image(mrs_roi, t1_aff.get_matrix())

	return(img)

