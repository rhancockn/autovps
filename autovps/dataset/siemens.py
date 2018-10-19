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
import json
from autovps.transform import Transform
from autovps.dataset import svsdata


class Siemens(object):
    """Utility class for manipulating Siemens SVS DICOM files.

    Attributes:
        path (str): Full path to the originating DICOM
        csa (obj): CSA header dictionary
        meta (obj): Selected CSA fields
        T (:obj:`Transform`): NIfTI compatible transform
    """

    csa = None
    qform = None
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
        """Get the value for the specified CSA key, with attempts to cast

        Args:
            k (str): a selected CSA key

        Returns:
            The value of the key
        """

        if k not in self.meta.keys():
            raise KeyError()


        # try to cast

        v = self.meta[k]

        # return non scalar values directly
        if isinstance(v, (list, tuple, np.ndarray)):
            return(v)

        # int
        iv = None
        try:
            iv = int(v)
        except ValueError:
            pass

        # float
        fv = None
        try:
            fv = float(v)
        except ValueError:
            pass

        if (iv is not None) and (fv is not None):
            if iv == fv:
                return(iv)
            else:
                return(fv)
        else:
            if iv is not None:
                return(iv)
            if fv is not None:
                return(fv)

        # a string, if not leave the exception unhandled
        return(str(v))



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
        self.qform = Transform(T_matrix)

        return(self.qform)

    def read_data(self, conj=True):
        """Read the associated fids.
        If the instance was initialized with a directory, the directory is assumed
        to contain a single series in order.

        Dimensions are channel x rep x mega x isis x t
        """

        print('Reading data...')
        if os.path.isfile(self.path):
            dcm = dcmwrapper.wrapper_from_file(self.path)
            data = np.array(_read_fid(dcm), ndmin=5)
            self.dcm_data = dcm.dcm_data

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
            self.dcm_data = lastdcm.dcm_data
            lastinstance = int(lastdcm.get((0x0020, 0x0013)).value)

            # figure out which channels are on
            csa_series = csareader.get_csa_header(lastdcm.dcm_data, 'series')
            csa_image = csareader.get_csa_header(lastdcm.dcm_data, 'image')
            siemens_hdr = csa_series['tags']['MrPhoenixProtocol']['items'][0]
            m = re.findall(r"""sCoilSelectMeas.aRxCoilSelectData\[0\].asList\[(?P<coilnum>\d+)\].sCoilElementID.tElement\t = \t""(?P<coilname>[HENC]+\d+)""""", siemens_hdr)
            channels = dict(m)
            channels = dict(zip(channels.values(), channels.keys()))
            n_channels = len(self.channels)
            self.channels = channels

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

            data = np.zeros((n_channels, n_reps, 1, 1, int(csareader.get_scalar(csa_image, 'DataPointColumns'))), dtype=complex)

            for fi in range(len(files)):
                dcmfile = os.path.join(self.path, files[fi])
                dcm = dcmwrapper.wrapper_from_file(dcmfile)
                csa = csareader.get_csa_header(dcm.dcm_data)
                fid = np.array(_read_fid(dcm), ndmin=5)
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
                        data[ci, ri, 0, 0, :] = fid
                else:
                    data[0, ri, 0, 0, :] = fid

            print('Read %d acquisitions from %d channels' % (n_reps, n_channels))
            print(channels.keys())
        # take the complex conjugate, which Tarquin seems to expect
        if conj:
            data = np.conj(data)


        # permute data for SPECIAL and MEGA

        if self.get_parameter('SequenceName') in ['eja_svs_mpress', 'eja_svs_mslaser']:
            data_on = data[:, 0:int(n_reps/2.), 0, ::]
            data_off = data[:, int(n_reps/2.):n_reps, 0, ::]
            data = np.stack((data_off, data_on), 2)
            #data = np.reshape(data, (n_channels, int(n_reps/2.), 2, 1, self.get_parameter('DataPointColumns')))

        self.data = data
        return(data)

    def get_svsdata(self):
        if self.data is None:
            self.read_data()

        if self.qform is None:
            self.calculate_transform()

        data = svsdata.SVSData()
        data.fid = self.data
        data.transform = self.qform
        data.sequence_name = self.get_parameter('SequenceName')
        data.tr = self.get_parameter('RepetitionTime')/1000.0
        data.te = self.get_parameter('EchoTime')/1000.0
        data.tm = self.get_parameter('MixingTime')/1000.0
        data.nucleus = self.get_parameter('ImagedNucleus')
        data.larmor = self.get_parameter('ImagingFrequency')*1.0e6
        dt = self.get_parameter('RealDwellTime')*1.0e-9
        sw = 1.0/dt
        data.sw = sw
        npts = self.get_parameter('DataPointColumns')
        data.f = np.arange(-sw/2.0+sw/(2*npts), sw/2, sw/npts)

        data.t = np.arange(0, npts*dt, dt)

        return(data)

    def get_sidecar(self):
        """Returns a BIDS-style sidecar
        """

        json_dict = _csa_to_dict(self.csa)
        # data has to be read in to count channels
        if self.data is None:
            self.read_data()
        json_dict['ReceiveCoilActiveElements'] = list(self.channels.keys())
        json_dict.pop('ImaCoilString', None)
        if self.qform is None:
            self.calculate_transform()
        json_dict['Transform'] = np.squeeze(np.asarray(self.qform.get_matrix())).tolist()

        dicom_tags = [((0x0008, 0x0060), 'Modality'),
                       ((0x0008, 0x0008), 'ImageType'),
                       ((0x0008, 0x0070), 'Manufacturer'),
                       ((0x0008, 0x0080), 'InstitutionName'),
                       ((0x0008, 0x0081), 'InstitutionAddress'),
                       ((0x0008, 0x1010), 'StationName'),
                       ((0x0008, 0x1030), 'StudyDescription'),
                       ((0x0008, 0x103e), 'SeriesDescription'),
                       ((0x0008, 0x1090), 'ManufacturersModelName'),
                       ((0x0018, 0x0015), 'BodyPartExamined'),
                       ((0x0018, 0x1000), 'DeviceSerialNumber'),
                       ((0x0018, 0x1020), 'SoftwareVersions'),
                       ((0x0018, 0x1030), 'ProtocolName'),
                       ((0x0018, 0x5100), 'PatientPosition'),
                       ((0x0020, 0x0011), 'SeriesNumber'),
                       ((0x0020, 0x0012), 'AcquisitionNumber'),
                       ((0x0020, 0x4000), 'ImageComments')
        ]

        for tag in dicom_tags:
            json_dict[tag[1]] = self.dcm_data.get(tag[0]).value

        json_dict.pop('ReferencedImageSequence', None)
        json_dict.pop('ScanningSequence')

        json_dict.update({'ImageType': list(json_dict['ImageType'])})
        json_dict.update({'SeriesNumber': int(json_dict['SeriesNumber'])})
        json_dict.update({'AcquisitionNumber': int(json_dict['AcquisitionNumber'])})
        return(json_dict)

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

def _csa_to_dict(csa_info):
    """Convert a csa header to a dictionary

    Arguments:
        csa_header (dict): the CSA record from `csareader.get_csa_header`

    Returns:
        json_dict: a dictionary with non-empty field value pairs
    """
    csa_dict = csa_info.get('tags')
    json_dict = dict()
    for k in csa_dict.keys():
        if csa_dict[k]['n_items'] != 0:
            item = csa_dict[k]['items']
            if len(item) > 1:
                json_dict[k] = item
            else:
                json_dict[k] = item[0]

    return(json_dict)
