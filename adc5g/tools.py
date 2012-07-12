from struct import pack, unpack
from numpy import array
from matplotlib.mlab import (
    psd, detrend_mean,
)
from opb import (
    OPB_CONTROLLER,
    OPB_DATA_FMT,
    inc_mmcm_phase,
    )
from spi import (
    get_spi_control,
    set_spi_control,
    )


def get_snapshot(roach, snap_name, bitwidth=8, man_trig=True, wait_period=2):
    """
    Reads a one-channel snapshot off the given 
    ROACH and returns the time-ordered samples.
    """

    grab = roach.snapshot_get(snap_name, man_trig=man_trig, wait_period=2)
    data = unpack('%iB' %grab['length'], grab['data'])

    return array(data)


def get_test_vector(roach, snap_name, bitwidth=8, man_trig=True, wait_period=2):
    """
    Sets the ADC to output a test ramp and reads off the ramp,
    one per core. This should allow a calibration of the MMCM
    phase parameter to reduce bit errors.

    core_a, core_c, core_b, core_d = get_test_vector(roach, snap_name)

    NOTE: This function requires the ADC to be in "test" mode, please use 
    set_spi_control(roach, zdok_n, test=1) before-hand to be in the correct 
    mode.
    """
    data = get_snapshot(roach, snap_name, bitwidth, man_trig=man_trig, wait_period=2)
    data_bin = (data>>1) ^ data
    return data_bin[0::4], data_bin[1::4], data_bin[2::4], data_bin[3::4]


def calibrate_mmcm_phase(roach, zdok_n, snap_name, bitwidth=8, man_trig=True, wait_period=2):
    """
    This function steps through all 56 steps of the MMCM clk-to-out 
    phase and finds total number of glitchss in the test vector ramp 
    per core. It then finds the least glitchy phase step and sets it.
    """
    ramp_max = 2**bitwidth - 1
    def total_glitches(core):
        glitches = 0
        for i in range(len(core)-1):
            diff = core[i+1] - core[i]
            if (diff!=1) and (diff!=-ramp_max):
                glitches += 1
        return glitches
    orig_control = get_spi_control(roach, zdok_n)
    set_spi_control(roach, zdok_n, test=1)
    glitches_per_ps = []
    for ps in range(56):
        core_a, core_c, core_b, core_d = get_test_vector(roach, snap_name, man_trig=man_trig, wait_period=2)
        glitches = total_glitches(core_a) + total_glitches(core_c) + \
            total_glitches(core_b) + total_glitches(core_d)
        glitches_per_ps.append(glitches)
        inc_mmcm_phase(roach, zdok_n)
    zero_glitches = [gl==0 for gl in glitches_per_ps]
    n_zero = 0
    longest_min = None
    while True:
        try:
            rising  = zero_glitches.index(True, n_zero)
            n_zero  = rising + 1
            falling = zero_glitches.index(False, n_zero)
            n_zero  = falling + 1
            min_len = falling - rising
            if min_len > longest_min:
                longest_min = min_len
                optimal_ps = rising + int((falling-rising)/2)
        except ValueError:
            break
    set_spi_control(roach, zdok_n, **orig_control)
    if longest_min==None:
        #raise ValueError("No optimal MMCM phase found!")
        return None, glitches_per_ps
    else:
        for ps in range(optimal_ps):
            inc_mmcm_phase(roach, zdok_n)
        return optimal_ps, glitches_per_ps


def get_psd(roach, snap_name, samp_freq, bitwidth=8, nfft=256):
    """
    Reads data off a given channel on a ROACH and calculates
    the power spectral density of the time-series.
    """

    data = get_snapshot(roach, snap_name, bitwidth)
    power, freqs = psd(data, nfft, Fs=samp_freq, detrend=detrend_mean, scale_by_freq=True)

    return power, freqs
