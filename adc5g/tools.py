from math import pi
from time import sleep
from struct import unpack
from matplotlib.mlab import (
    psd, detrend_mean,
)
from numpy import (
    array, sinc, sum, max, log10,
)


CTRL_FMT = "snap_{chan}_ctrl"
BRAM_FMT = "snap_{chan}_bram_{lmsb}"


def interleave(lsb, msb, bitwidth=8):
    """
    Re-interleaves the snapshots from the 32-bit LSB 
    and MSB block memories to form the correctly
    ordered original snapshot.
    """

    out = []
    demux = int(32/bitwidth)

    for i in range(0, len(lsb), demux):
        out.extend(msb[i:i+demux])
        out.extend(lsb[i:i+demux])

    return array(out)


def get_snapshot(roach, chan, length=2**14, 
                 bitwidth=8, ctrl_fmt=CTRL_FMT, bram_fmt=BRAM_FMT, 
                 force=True):
    """
    Reads a one-channel, 64-bit snapshot off the given 
    ROACH and returns the time-ordered samples.
    """

    roach.write_int(ctrl_fmt.format(chan=chan), 0)
    roach.write_int(ctrl_fmt.format(chan=chan), 1)

    if force:
        roach.write_int('force_capture', 0)
        roach.write_int('force_capture', 1)

    lsb_bram_name = bram_fmt.format(chan=chan, lmsb='lsb')
    lsb_raw = roach.read(lsb_bram_name, length)
    lsb = unpack('%iB'%length, lsb_raw)

    msb_bram_name = bram_fmt.format(chan=chan, lmsb='msb')
    msb_raw = roach.read(msb_bram_name, length)
    msb = unpack('%iB'%length, msb_raw)

    return interleave(lsb, msb, bitwidth)


def get_psd(roach, chan, samp_freq, length=2**14,
            bitwidth=8, ctrl_fmt=CTRL_FMT, bram_fmt=BRAM_FMT,
            force=True, NFFT=256):
    """
    Reads data off a given channel on a ROACH and calculates
    the power spectral density of the time-series.
    """

    data = get_snapshot(roach, chan, length, bitwidth, ctrl_fmt, bram_fmt, force)
    power, freqs = psd(data, NFFT=NFFT, Fs=samp_freq, detrend=detrend_mean, scale_by_freq=True)

    return power, freqs


def get_sinad(power, freqs=None, peak_width=10):
    """ 
    Computes the SIgnal-to-Noise-And-Distortion ratio.
    """

    total_power = power.sum()
    peak_argfreq = power.argmax()
    peak_power = sum(power[peak_argfreq-peak_width/2:peak_argfreq+peak_width/2])
    sinad = 10*log10(total_power/(total_power-peak_power))

    return sinad


def go(freq):
    q = adc5g.get_snapshot(r, 'q')
    q1 = ceil(fitfunc(leastsq(errfunc, [6., 2500./freq, 0., 8.], args=(Tx, q))[0], Tx))
    figure(1).clear(); plot(q, '.'); plot(q1, '-'); plot(q-q1); xlim(0, 100); figure(2).clear(); psd(q); psd(q-q1, linestyle='--');
    sinad = 10*log10(sum(q)/sum(q-q1))
    enob = (sinad-1.76)/6.02
    return sinad, enob
