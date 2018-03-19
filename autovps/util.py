GYROMAGNETIC_RATIO = {'1H': 42577478, '13C': 10705000, '31P': 17235000}

PPM_OFFSET = 4.65

def f_to_ppm(f, larmor):
    """ Converts frequency (Hz) to ppm

    Parameters:
        f (float): frequency or array in Hz
        larmor (float): Larmor frequency (in Hz)

    Returns:
        ppm (float)
    """

    ppm = f/larmor*1.0e6 + PPM_OFFSET
    return(ppm)


def ppm_to_f(ppm):
    f = ppm - PPM_OFFSET * larmor * 1e-6;
    return(f)
