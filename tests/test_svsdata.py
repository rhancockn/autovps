import numpy as np
from pytest import approx, warns

import autovps.dataset.siemens as siemens
from autovps.dataset import svsdata


def test_axis():
    siemens_data = siemens.Siemens('tests/data/siemens/eja_svs_press_combined_noave')
    data = siemens_data.get_svsdata()

    sw = data.sw
    npts = data.fid.shape[-1]

    f = np.arange(-sw/2.0+sw/(2*npts), sw/2, sw/npts)

    # is the frequnecy axis correct?
    orig_f = data.f
    assert approx(f) == orig_f
    assert not (approx(data.f) != orig_f)

    # shift the ppm axis
    ppm_shifted = data.ppm - data._ppmshift
    data.ppm = ppm_shifted
    shifted_f = data.f
    assert approx(shifted_f) != orig_f


def test_nifti():
    with warns(UserWarning, match='Assuming channels are combined'):
        siemens_data = siemens.Siemens(
            'tests/data/siemens/eja_svs_press_combined_noave')
        data = siemens_data.get_svsdata()
        data.save_nifti('test.nii')
