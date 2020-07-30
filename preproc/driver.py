%load_ext autoreload
%autoreload 2

import autovps.dataset.siemens as siemens
from autovps.dataset import svsdata
from autovps import make_voi
import nibabel as nib

base_path = '/Users/roh17004/Downloads/104DI'
#press
data = siemens.Siemens(base_path + '/12')
svs = data.get_svsdata()
svs.save_fida(base_path + '/press.mat')

data = siemens.Siemens(base_path + '/15')
svs = data.get_svsdata()
svs.save_fida(base_path + '/MEGA1.mat')

data = siemens.Siemens(base_path + '/18')
svs = data.get_svsdata()
svs.save_fida(base_path + '/MEGA2.mat')

data = siemens.Siemens(base_path + '/21')
svs = data.get_svsdata()
svs.save_fida(base_path + '/MEGA3.mat')

tform= data.calculate_transform()
t1=nib.load(base_path + '/3/104DI-T1w_MPR.nii.gz')

img = make_voi.make_voi(t1, tform)
nib.save(img, base_path+'/roi.nii.gz')

x=data.get_sidecar()

with open('test.json', 'w+') as fp:
    json.dump(x, fp, indent=3)

# save the transform
np.savetxt('%s_from-Device_to-orig_mode-image_xfm.mat' % source_bids, data.qform.get_matrix())
