import numpy as np
import scipy.io
import nibabel as nib
import warnings
import autovps.util as util

class SVSData(object):
    """Class for SVS spectroscopy data

    Attributes:
        tr (float): TR (sec)
        te (float): echo time (sec)
        tm (float): mixing time (sec)
        f (list): frequency axis in Hz
        ppm (list): frequency axis in ppm
        larmor (float): Larmor frequency in Hz
        sw (float): spectral width
        fid (:obj:`numpy.array`): complex repetition x acquisition x channel x time
        spec (:obj:`numpy.array`): complex repetition x acquisition x channel x time
        isis_dim (int): index of ISIS dimension or None
        mega_dim (int): index of on/off dimension or None
        sequence_name (str): name of the acquisition sequence
        sequence_type (str): type of sequence: STEAM, PRESS, MEGAPRESS, SPECIAL, MEGASPECIAL, LASER, sLASER, UNKNOWN

    """
    def __init__(self):
        self._tr = None
        self._te = None
        self._tm = None
        self._fid = None
        self._spec = None
        self._f = None
        self._ppm = None
        self.larmor = None
        self._sequence_type = 'UNKNOWN'
        self._sequence_name = 'UNKNOWN'
        self.nucleus = '1H'
        self.isis_dim = 3
        self.mega_dim = 2
        self._ppmshift = 4.65
        self.sw = None
        self.t = None
        self.transform = None
        self.is_mega = False
        self.is_special = False
        self.dimnames = {'channel': 0, 'rep': 1, 'mega': 2, 'isis': 3, 't': 4}


    @property
    def tr(self):
        return(self._tr)

    @tr.setter
    def tr(self, tr):
        self._tr = tr

    @property
    def te(self):
        return(self._te)

    @te.setter
    def te(self, te):
        self._te = te


    @property
    def tm(self):
        return(self._tm)

    @tm.setter
    def tm(self, tm):
        self._tm = tm

    @property
    def fid(self):
        """Get the fid.

        If there is a spec, but no fid, the spec is iffted and returned.
        The frequency axis is unchanged since this may have been shifted.
        """

        if self._fid is not None:
            return(self._fid)

        if self._spec is not None:
            self._fid = np.fft.ifft(np.fft.ifftshift(self._spec))
            return(self._fid)

        return(None)

    @fid.setter
    def fid(self, fid):
        self._fid = fid
        self._spec = None

    @property
    def spec(self):
        """Get the spectrum.

        If there is a fid, but no spec, the fid is ffted and returned.
        The frequency axis is unchanged since this may have been shifted
        """
        if self._spec is not None:
            return(self._spec)

        if self._fid is not None:
            self._spec = np.fft.fftshift(np.fft.fft(self._fid))
            return(self._spec)

        return(None)

    @spec.setter
    def spec(self, spec):
        self._spec = spec
        self._fid = None

    @property
    def ppm(self):
        """Return the ppm axis.

        If there is no stored ppm and the Larmor frequency and f axis are defined,
        the ppm axis will be calculated and stored

        """
        if self._ppm is not None:
            return(self._ppm)
        if self.larmor is None:
            return(None)
        if self._f is not None:
            self._ppm = self._f/(self.larmor)*1.0e6 + self._ppmshift
            return(self._ppm)

        return(None)

    @ppm.setter
    def ppm(self, ppm):
        """Set the ppm axis.

        This will invalidate the frequency axis to ensure consistency
        """
        self._ppm = ppm
        self._f = None

    @property
    def f(self):
        """Return the frequency axis in Hz.

        If there is no stored f and the Larmor frequency and ppm axis are defined,
        the f axis will be calculated and stored

        """

        if self._f is not None:
            return(self._f)
        if self.larmor is None:
            return(None)
        if self._ppm is not None:
            self._f = (self._ppm - self._ppmshift) * self.larmor * 1e-6;
            return(self._f)
        return(None)

    @f.setter
    def f(self, f):
        """Set the frequency axis.

        This will invalidate the ppm axis to ensure consistency
        """
        self._ppm = None
        self._f = f

    @property
    def sequence_type(self):
        return (self._sequence_type)
    @sequence_type.setter
    def sequence_type(self, sequence_type):
        known_types = {'STEAM', 'PRESS', 'LASER', 'sLASER', 'MEGAPRESS', 'MEGASPECIAL', 'UNKNOWN'}

        if sequence_type in known_types:
            self._sequence_type = sequence_type
        else:
            warnings.warn('%s is not a known sequence type' % sequence_type)
            self._sequence_type = 'UNKNOWN'

        if sequence_type in {'MEGAPRESS', 'MEGASPECIAL'}:
            self.is_mega = True

        if sequence_type in {'SPECIAL', 'MEGASPECIAL'}:
            self.is_special = True


    @property
    def sequence_name(self):
        return (self._sequence_name)

    @sequence_name.setter
    def sequence_name(self, sequence_name):
        self._sequence_name = sequence_name

        if sequence_name in {'eja_svs_press', 'svs_se'}:
            self.sequence_type = 'PRESS'
        elif sequence_name in {'eja_svs_mpress'}:
            self.sequence_type = 'MEGAPRESS'
        elif sequence_name in {'eja_svs_mpress'}:
            self.sequence_type = 'MEGAPRESS'
        elif sequence_name in {'eja_svs_steam', 'svs_st'}:
            self.sequence_type = 'STEAM'
        elif sequence_name in {'eja_svs_slaser'}:
            self.sequence_type = 'sLASER'
        elif sequence_name in {'eja_svs_laser'}:
            self.sequence_type = 'LASER'
        elif sequence_name in {'eja_svs_mslaser'}:
            self.sequence_type = 'MEGAsLASER'
        else:
            self.sequence_type = 'UNKNOWN'


    def has_isis(self):
        return(self.fid.shape[self.isis_dim] > 1)

    def has_mega(self):
        return(self.fid.shape[self.mega_dim] > 1)


    def save_fida(self, path):
        """ Convert the dataset to a FID-A compatible Matlab structure
        """
        # input is Dimensions are channel x rep x mega x isis x t
        # FID-A seems to accomodate only 4: t x chan x rep x subSpecs
        # TODO: see if ISIS and MEGA subspecs are differentiated

        # permute the axes to t x chan x rep x mega x isis
        fids = np.transpose(self.fid, (4, 0,1,2,3))
        specs = np.transpose(self.spec, (4, 0,1,2,3))
        # reshape to combine subspecs
        dims = list(fids.shape[0:-2])
        dims.append(-1)
        fids = np.reshape(fids, tuple(dims))
        specs = np.reshape(specs, tuple(dims))

        # remove last dimensi if there are no subSpecs
        fids = np.squeeze(fids)
        specs = np.squeeze(specs)

        # fp to avoid int64 errors
        dim_dict = {'t': 1.0, 'coils': 2.0, 'averages': 3.0, 'subSpecs': 0.0, 'extras': 0.0}

        # there are still subSpectra
        if fids.ndim == 4:
            subspecs = fids.shape[-1]
            rawSubspecs = fids.shape[-1]
            dim_dict['subSpecs'] = 4.0
        else:
            subspecs = 0
            rawSubspecs = 0

        if self.fid.shape[0]==1:
            addedrcvrs = 1
        else:
            addedrcvrs = 0

        B0 = self.larmor/util.GYROMAGNETIC_RATIO[self.nucleus]

        n_averages = float(self.fid.shape[self.dimnames['rep']])
        # fids - time domain MRS data.
        # specs - frequency domain MRS data.
        # t - vector of time values for plotting in the time domain [s]
        # ppm - vector of frequency values for plotting in the frequency domain
        # [ppm]
        # sz - size of the fids and specs arrays
        # date - date that the data was acquired or simulated
        # averages - number of averages in the dataset (possibly altered by
        # processing)
        # rawAverages - number of averages in the original dataset (not altered by
        # processing).
        # subspecs - number of subspectra (ISIS, edit on/off, etc) in the dataset
        # (possibly altered by processing).
        # rawSubspecs - number of subspectra (ISIS, edit on/off, etc) in the original
        # dataset (not altered by processing). Bo - magnetic field strength [Tesla]
        # txfrq - Centre frequnecy [MHz];
        # linewidth - linewidth of data (only used for simulated data) [Hz]
        # n - number of spectral points
        # dwelltime - dwell time of the data in the time domain [s] (dwelltime =
        # 1/spectralwidth)
        # sim - type of simulation (ideal vs. shaped pulses), only used for
        #               simulated data.
        # te seq dims
        # - echo time of acquisition [ms], only used for simulated data - type of sequence used (only used for simulated data).
        # - structure specifying which data dimensions are stored along
        # which dimensions of the fids/specs arrays. Fields include:
        # t - time/frequency dimension (usually this is 1, the first
        # dimension of the fids/specs array).
        # coils - for multiple receiver array, this is the dimension of
        # the arrayed receiver data (can be 2, 3 or 4). averages - for multiple averages, this is the dimension of the
        # averages (can be 2, 3 or 4).
        # subSpecs - in the case of subtraction data (ISIS, MEGA-PRESS), this
        # is the dimension of the subSpectra (can be 2, 3 or 4).


        mdict = {'fids': fids, 'specs': specs, 't': self.t,
                 'ppm': self.ppm, 'sz': np.float_(fids.shape), 'date': '',
                 'averages': n_averages, 'rawAverages': n_averages,
                 'subspecs': float(subspecs), 'rawSubspecs': float(rawSubspecs), 'Bo': B0,
                 'txfrq': self.larmor, 'dwelltime': 1.0/self.sw,
                 'spectralwidth': self.sw, 'seq': self._sequence_name,
                 'dims': dim_dict, 'te': self.te * 1e3, 'tr': self.tr * 1e3,
                 'pointsToLeftshift': 0}

        # writtentostruct
        # gotparams
        # filtered
        # zeropadded
        # freqcorrected
        # phasecorrected
        # averaged
        # addedrcvrs
        # Subtracted
        # Writtentotext
        # Downsampled
        # avgNormalized
        # isISIS
        # - Has the dataset been written to a structure (1 or 0)
        # - Have the parameters been retrieved from the dataset (1 or 0)
        # - Has the dataset been filtered (1 or 0)
        # - Has the dataset been zeropadded (1 or 0)
        # - Has the dataset been frequency corrected (1 or 0) - Has the dataset been phase corrected (1 or 0)
        # - Have the averages been combined (1 or 0)
        # - Have the rcvr channels been combined (1 or 0).
        # - Have the subspecs been subtracted (1 or 0)
        # - Has the data been written to text file (1 or 0) - has the data been resampled to a different
        # spectral resolution (1 or 0)
        # - Has the data been amplitude scaled following
        # combination of the averages (1 or 0)
        # - Does the dataset contain ISIS subspectra (1 or 0)

        flags = {'writtentostruct': 1, 'gotparams': 1, 'filtered': 0,
                 'zeropadded': 0, 'freqcorrected': 0, 'phasecorrected': 0,
                 'averaged': int(n_averages == 1), 'addedrcvrs': addedrcvrs,
                 'Subtracted': 0, 'Writtentotext': 0, 'Downsampled': 0,
                 'avgNormalized': 0, 'isISIS': int(self.is_special),
                 'leftshifted': 0}

        if self.sequence_type == 'STEAM':
            mdict['tm'] = self.tm

        mdict['flags'] = flags
        scipy.io.savemat(path, {'svs': mdict}, format='5', long_field_names=True)

    def save_nifti(self, path):
        """ Convert the dataset to a NIFTI
        """
        meta = {'te': self.te, 'tr': self.tr, 'sw': self.sw}
        if self.sequence_type == 'STEAM':
            meta['tm'] = self.tm

        # store real and imaginary components in last 2 dims
        component_fid = np.stack((np.real(self.fid),np.imag(self.fid)), -2)
        nifti = nib.Nifti2Image(component_fid, self.transform.get_matrix(), extra=meta)
        nib.save(nifti, path)
