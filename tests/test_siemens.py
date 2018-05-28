import struct
import numpy as np
import six  # for 2/3 compatibility
import nibabel.nicom.dicomwrappers as dcmwrapper
from pytest import approx

import autovps.dataset.siemens as siemens
from autovps.dataset import svsdata


DIM_TOL = .1


def test_orientations():
    """ Test that the axis string is correctly recovered
    """
    siemens_orientations = {'3': 'Tra', '4': 'Tra', '5': 'Tra', '6': 'Sag',
                            '7': 'Cor', '8': 'Sag', '9': 'Cor', '10': 'Tra>Sag 30.0',
                            '11': 'Tra>Cor 30.0', '12': 'Tra>Cor 30.0 >Sag 15.0',
                            '13': 'Tra>Sag 30.0 >Cor 15.0', '14': 'Tra>Sag -30.0',
                            '15': 'Tra>Cor -30.0', '16': 'Tra>Cor 30.0 >Sag -15.0',
                            '17': 'Tra>Sag 30.0 >Cor -15.0', '18': 'Tra>Sag 30.0',
                            '19': 'Tra>Cor 30.0', '20': 'Tra>Cor 30.0 >Sag 15.0',
                            '21': 'Tra>Sag 30.0 >Cor 15.0'}

    siemens_inplane_orientations = {'3': 0.0, '4': 30.0, '5': -30.0, '6': 0.0,
                                    '7': 0.0, '8': 30.0, '9': 30.0, '10': 0.0,
                                    '11': 0.0, '12': 0.0, '13': 0.0, '14': 0.0,
                                    '15': 0.0, '16': 0.0, '17': 0.0, '18': 20.0,
                                    '19': 20.0, '20': 20.0, '21': 20.0}
    for idx, ori in six.iteritems(siemens_orientations):
        data = siemens.Siemens('tests/data/siemens/%s' % idx)
        data.calculate_transform()
        R = data.qform.get_rotation()

        print('%s: %f' % (idx, np.degrees(np.arctan2(-R[2,1], -R[2,0]))))

        # orientation string is correct
        assert ori == data.tform.siemens_orientation()
        size = data.tform.get_scale()
        real_size = [20, 25, 30]
        position = data.tform.get_position()
        for i in range(3):
            assert abs(size[i]-real_size[i]) < DIM_TOL
            assert abs(position[i]) < DIM_TOL

        #in plane rotation
        assert approx(np.degrees(data.meta['VoiInPlaneRotation'])) \
            == siemens_inplane_orientations[idx]


def test_read_single():
    """Test reading a single combined DICOM
    """
    test_file = 'tests/data/siemens/3/S8457LTU_2_3_00001_00001_173218510000__1263534865.dcm'
    data_file = siemens.Siemens(test_file)
    svs_file = data_file.get_svsdata()

    data_dir = siemens.Siemens('tests/data/siemens/3')
    svs_dir = data_dir.get_svsdata()

    # is it the right size?
    assert data_file.data.shape == data_dir.data.shape == (1, 1, 1, 1, 2048)

    dcm = dcmwrapper.wrapper_from_file(test_file)
    packed = dcm.get((0x7fe1, 0x1010)).value
    data = struct.unpack("<%df" % (len(packed) / 4), packed)
    cmplx = [data[i]+data[i+1]*1j for i in range(0, len(data), 2)]

    assert (data_file.data == np.conj(cmplx)).all()
    assert (data_file.data == data_dir.data).all()



def test_read_uncombined():
    """Test reading a directory of dicoms
    """
    data = siemens.Siemens('tests/data/siemens/eja_svs_press_uncombined')
    data.get_svsdata()

    assert data.data.shape == (52, 32, 1, 1, 2048)
    for i in range(data.data.shape[0]):
        for j in range(data.data.shape[1]):
            assert data.data[i, j, 0, 0, -1] != 0



def test_read_unaveraged():
    """Test reading a directory of dicoms
    """
    data = siemens.Siemens('tests/data/siemens/eja_svs_press_combined_noave')
    data.get_svsdata()

    assert data.data.shape == (1, 8, 1, 1, 2048)
    for i in range(data.data.shape[0]):
        for j in range(data.data.shape[1]):
            assert data.data[i, j, 0, 0, -1] != 0
