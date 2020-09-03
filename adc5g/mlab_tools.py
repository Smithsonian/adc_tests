from matplotlib.mlab import psd, detrend_mean
from .tools import get_snapshot


def get_psd(roach, snap_name, samp_freq, bitwidth=8, nfft=256):
    """
    Reads data off a given channel on a ROACH and calculates
    the power spectral density of the time-series.
    """

    data = get_snapshot(roach, snap_name, bitwidth)
    power, freqs = psd(data, nfft, Fs=samp_freq, detrend=detrend_mean, scale_by_freq=True)

    return power, freqs
