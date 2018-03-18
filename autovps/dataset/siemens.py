"""
Utility class for manipulating Siemens SVS DICOM files
"""


import os
import struct
import re
import warnings
import nibabel.nicom.csareader as csareader
import nibabel.nicom.dicomwrappers as dcmwrapper
import numpy as np
from autovps.transform import Transform


class Siemens(object):
    """Utility class for manipulating Siemens SVS DICOM files.

    Attributes:
        path (str): Full path to the originating DICOM
        csa (obj): CSA header dictionary
        meta (obj): Selected CSA fields
        T (:obj:`Transform`): NIfTI compatible transform
    """

    csa = None
    T = None
    path = None
    meta = {}
    data = None

    def __init__(self, path):
        """Initialize a Siemens object from a DICOM file or directory.
        If a directory is given, read the first DICOM file in the directory.


        Args:
            path (str): The path to a Siemens SVS DICOM file or directory of 
                DICOM files.

        """
        self.path = os.path.realpath(path)
        dcmfile = None
        if os.path.isdir(self.path):
            # find the first dicom file to read meta data from
            files = os.listdir(self.path)
            # exclude any non dicoms
            files = [f for f in files if (os.path.isfile(os.path.join(self.path, f)) and f.endswith(('.DCM', '.dcm')))]
            files = sorted(files)
            dcmfile = os.path.join(self.path, files[0])

            # a single file, change the path
            if len(files) == 1:
                path = os.path.join(self.path, files[0])
                self.path = path

            if not(dcmfile):
                raise FileNotFoundError('No DICOM files found')

        else:
            dcmfile = path

        dcm = dcmwrapper.wrapper_from_file(dcmfile)
        self.csa = csareader.get_csa_header(dcm.dcm_data)

        scalar_fields = ['MagneticFieldStrength', 'ImagingFrequency',
                         'MixingTime', 'EchoTime', 'RepetitionTime', 'ImaCoilString',
                         'SequenceName', 'VoiReadoutFoV', 'VoiPhaseFoV', 'VoiThickness',
                         'VoiInPlaneRotation', 'DataPointColumns', 'RealDwellTime',
                         'PixelBandwidth', 'ImagedNucleus']

        for k in scalar_fields:
            self.meta[k] = csareader.get_scalar(self.csa, k)

        self.meta['ImagePositionPatient'] = \
            csareader.get_vector(self.csa, 'ImagePositionPatient', 3)
        self.meta['VoiPosition'] = csareader.get_vector(self.csa, 'VoiPosition', 3)
        self.meta['VoiOrientation'] = csareader.get_vector(self.csa, 'VoiOrientation', 3)
        self.meta['ImageOrientationPatient'] = \
            csareader.get_vector(self.csa, 'ImageOrientationPatient', 6)

    def get_parameter(self, k):
        """Get the value for the specified CSA key.

        Args:
            k (str): a selected CSA key

        Returns:
            The value of the key
        """

        return(self.meta[k])

    def calculate_transform(self):
        """Form the linear transform matrix.

        Returns:
            Transform: a Transform

        """

        image_orientation = np.array(self.get_parameter('ImageOrientationPatient'))
        flip = np.array([-1, -1, 1])
        R = np.vstack((image_orientation[0:3]*flip, image_orientation[3:6]*flip))
        R = np.vstack((R, -np.cross(R[0, :], R[1, :]))).T

        T_matrix = np.eye(4)
        T_matrix[0:3, 0:3] = R

        T_matrix = np.dot(T_matrix, np.diag([self.get_parameter('VoiReadoutFoV'),
            self.get_parameter('VoiPhaseFoV'), self.get_parameter('VoiThickness'), 1]))
        T_matrix[0:3, 3] = np.array(self.get_parameter('VoiPosition'))*flip

        self.T_matrix = T_matrix
        self.T = Transform(T_matrix)

        return(self.T)

    def read_data(self, conj=True):
        """Read the associated fids. 
        If the instance was initialized with a directory, the directory is assumed
        to contain a single series in order.

        """

        if os.path.isfile(self.path):
            dcm = dcmwrapper.wrapper_from_file(self.path)
            data = np.array(_read_fid(dcm), ndmin=3)

        elif os.path.isdir(self.path):
            # read a directory of DICOMS, each containing one fid
            # the directory must contain more than one file, as this is checked
            # in the class constructor

            channels = []

            files = os.listdir(self.path)
            # exclude any non dicoms
            files = [f for f in files
                     if (os.path.isfile(os.path.join(self.path, f))
                     and f.endswith(('.DCM', '.dcm')))]
            files = sorted(files)

            # find the instance number of the last dicom to calculate data size
            # this assumes different interleaved acquisitions have the same 
            # (0020, 0012) Acquisition Number 
            # true for eja sequences
            lastdcmfile = os.path.join(self.path, files[-1])
            lastdcm = dcmwrapper.wrapper_from_file(lastdcmfile)
            lastinstance = int(lastdcm.get((0x0020, 0x0013)).value)

            # figure out which channels are on
            csa_series = csareader.get_csa_header(lastdcm.dcm_data, 'series')
            csa_image = csareader.get_csa_header(lastdcm.dcm_data, 'image')
            siemens_hdr = csa_series['tags']['MrPhoenixProtocol']['items'][0]
            m = re.findall(r"""sCoilSelectMeas.aRxCoilSelectData\[0\].asList\[(?P<coilnum>\d+)\].sCoilElementID.tElement\t = \t""(?P<coilname>[HENC]+\d+)""""", siemens_hdr)
            channels = dict(m)
            channels = dict(zip(channels.values(), channels.keys()))
            n_channels = len(channels)


            # is the data combined over channels?
            # mri_probedicom reports ucUncombImages, but where is this in the CSA?
            

            # the first two instances of uncombined eja sequences are single channels
            # TODO: figure out what they are
            # Assume the first match the last two, which are missing the same channel
            n_reps = lastinstance - 2

            is_combined = False
            # TODO: handle channel combined data
            if len(files) != (n_reps*(n_channels + 1)):
                # not enough files for uncombined data
                if len(files) == lastinstance:
                    warnings.warn('Assuming channels are combined')
                    n_reps = lastinstance
                    n_channels = 1
                    is_combined = True

                else:
                    raise Exception('Expected n_reps[%d] * (n_channels[%d] + 1 files' 
                                    % (n_reps, n_channels))

            data = np.zeros((n_reps, n_channels, int(csareader.get_scalar(csa_image, 'DataPointColumns'))), dtype=complex)

            for fi in range(len(files)):
                dcmfile = os.path.join(self.path, files[fi])
                dcm = dcmwrapper.wrapper_from_file(dcmfile)
                csa = csareader.get_csa_header(dcm.dcm_data)
                fid = np.array(_read_fid(dcm), ndmin=3)
                channel = csareader.get_scalar(csa, 'ImaCoilString')
                inst = int(dcm.get((0x0020, 0x0013)).value)

                if not is_combined:
                    if inst <= 2:
                        ri = n_reps - inst
                    else:
                        ri = inst - 2 - 1
                else:
                    ri = inst - 1

                # there are combined coils (HC1-7) in the dicoms?
                # make sure this channel is one that is turned on
                if not is_combined:
                    if channel in channels.keys():
                        ci = int(channels[channel])
                        data[ri, ci, :] = fid
                else:
                    data[ri, 0, :] = fid
        
        # take the complex conjugate, which Tarquin seems to expect            
        if conj:
            data = np.conj(data)

        self.data = data
        return(data)

    def to_svsdata(self):
        data = dataset.SVSData()
        data.sequence_name = self.get_parameter('SequenceName')
        data.sequence_type = 'UNKNOWN'

        if data.sequence_name in ['eja_svs_press', 'svs_se']:
            data.sequence_type = 'PRESS'
        elif data.sequence_name in ['eja_svs_steam', 'svs_st']:
            data.sequence_type = 'STEAM'
        elif data.sequence_name in ['eja_svs_mpress']:
            data.sequence_type = 'MEGAPRESS'
        else:
            warnings.warn('Unknown or unsupported sequence')



def _read_fid(dcm):
    """Read a fid from a PyDicom wrapper.

    This converts the data in (0x7fe1, 0x1010) from a series of
    little endian 4 byte real-imaginary pairs to a list of complex values.

    Args:
        dcm (:obj:`SiemensWrapper`): a Siemens SVS DICOM object

    Returns:
        fid (obj): the complex-valued fid
    """
    TAG = (0x7fe1, 0x1010)

    packed = dcm.get(TAG).value
    data = struct.unpack("<%df" % (len(packed) / 4), packed)
    cmplx = [data[i]+data[i+1]*1j for i in range(0, len(data), 2)]

    return(cmplx)







