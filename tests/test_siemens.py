import pytest
import autovps.dataset.siemens as siemens
import nibabel.nicom.dicomwrappers as dcmwrapper
import numpy as np
import struct

DIM_TOL = .1
def test_orientations():
	siemens_orientations = {'3':'Tra', '4':'Tra', '5':'Tra', '6':'Sag', '7':'Cor', 
	'8':'Sag', '9':'Cor', '10':'Tra>Sag 30.0', '11':'Tra>Cor 30.0', 
	'12':'Tra>Cor 30.0 >Sag 15.0', '13': 'Tra>Sag 30.0 >Cor 15.0',
	'14': 'Tra>Sag -30.0', '15': 'Tra>Cor -30.0',
	'16': 'Tra>Cor 30.0 >Sag -15.0', '17': 'Tra>Sag 30.0 >Cor -15.0',
	'18': 'Tra>Sag 30.0', '19':'Tra>Cor 30.0', 
	'20': 'Tra>Cor 30.0 >Sag 15.0', '21': 'Tra>Sag 30.0 >Cor 15.0' }
	for idx,ori in siemens_orientations.iteritems():
		data = siemens.Siemens('tests/data/siemens/%s' % idx)
		data.calculate_transform()
		
		#orientation string is correct
		assert ori==data.T.siemens_orientation()
		size = data.T.get_scale()
		real_size = [20,25,30]
		position = data.T.get_position() 
		for i in range(3):
			assert abs(size[i]-real_size[i]) < DIM_TOL
			assert abs(position[i]) < DIM_TOL
		

def test_read_single():
	test_file = 'tests/data/siemens/3/S8457LTU_2_3_00001_00001_173218510000__1263534865.dcm'
	data_file = siemens.Siemens(test_file)
	data_file.read_data()

	data_dir = siemens.Siemens('tests/data/siemens/3')
	data_dir.read_data()

	#is it the right size?
	assert data_file.data.shape == data_dir.data.shape == (1,1,2048)

	dcm = dcmwrapper.wrapper_from_file(test_file)
	packed=dcm.get((0x7fe1, 0x1010)).value
	data = struct.unpack("<%df" % (len(packed) / 4), packed)
	cmplx = [data[i]+data[i+1]*1j for i in range(0,len(data),2)]

	assert (data_file.data==np.conj(cmplx)).all()
	assert (data_file.data==data_dir.data).all()




