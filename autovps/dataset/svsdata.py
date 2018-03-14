import numpy as np

class svsdata(object):
"""Class for SVS spectroscopy data

Attributes:
	tr (float): TR (sec)
	te (float): echo time (sec)
	tm (float): mixing time (sec)
	fid (:obj:`numpy.array`): complex repetition x acquisition x channel x time
	spec (:obj:`numpy.array`): complex repetition x acquisition x channel x time
	f (list): frequency axis in Hz
	ppm (list): frequency axis in ppm
	larmor (float): Larmor frequency in Hz
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
		self._txfrq = None
		self.sequence_type = 'UNKNOWN'
		self.sequence_name = 'UNKNOWN'
		self.nucleus = '1H'


	@property
	def tr(self):
		return(self._tr)

	@tr.setter
	def tr(self, tr):
		self._tr = tr

	@property
	def fid(self):
		"""Get the fid.

		If there is a spec, but no fid, the spec is iffted and returned.
		The frequency axis is unchanged since this may have been shifted.
		"""

		if self._fid:
			return(self._fid)
		elif self._spec:
			self._fid = np.fft.ifft(np.fft.ifftshift(self._spec))
			return(self._fid)
		else:
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
		if self._spec:
			return(self._spec)
		elif self._fid:
			self._spec = np.fft.fftshift(np.fft.fft(self._fid))
			return(self._spec)
		else:
			return(None)

	@spec.setter
	def spec(self, spec):
		self._spec = spec
		self._fid = None



	@property
	def ppm(self):
		if self._ppm:
			return(self._ppm)
		elif self._f:
			self._ppm = self._f/(self.B0*util.gyromagnetic_ratio[self.nucleus])
		else:
			return(None)

	def has_isis(self):
		if not isis_dim:
			return(False)
		else:
			return(True)

	def is_mega(self):
		if not mega_dim:
			return(False)
		else:
			return(True)


	def to_fida(self, path=None):
		""" Convert the dataset to a FID-A compatible Matlab structure
		"""

	def to_nifti(self, path=None):
		""" Convert the dataset to a NIFTI
		"""


