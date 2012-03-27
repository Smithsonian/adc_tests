from math import pi, floor
from time import sleep
from struct import pack, unpack
from matplotlib.mlab import (
    psd, detrend_mean,
)
from numpy import (
    array, zeros, sinc, sum, max, log10
)


SPI_DATA_FMT = '>H2B'
CONTROL_REG_ADDR = 0x01 + 0x80
CHANSEL_REG_ADDR = 0x0f + 0x80
EXTOFFS_REG_ADDR = 0x20 + 0x80
EXTGAIN_REG_ADDR = 0x22 + 0x80
EXTPHAS_REG_ADDR = 0x24 + 0x80
FIRST_EXTINL_REG_ADDR = 0x30 + 0x80
CALCTRL_REG_ADDR = 0x10 + 0x80
ADC_CONTROLLER = 'adc5g_controller'


def set_spi_register(roach, zdok_n, reg_addr, reg_val):
    """
    Sets the value of an ADC's register over SPI
    """
    spi_data = pack(SPI_DATA_FMT, reg_val, reg_addr, 0x01)
    roach.blindwrite(ADC_CONTROLLER, spi_data, offset=0x4+zdok_n*0x4)


def get_spi_register(roach, zdok_n, reg_addr):
    """
    Gets the value of an ADC's register over SPI
    """
    spi_data = pack(SPI_DATA_FMT, 0x0, reg_addr, 0x01)
    roach.blindwrite(ADC_CONTROLLER, spi_data, offset=0x4+zdok_n*0x4)
    raw = roach.read(ADC_CONTROLLER, 0x4, offset=0x4+zdok_n*0x4)
    #print raw.encode('string_escape')
    reg_val, old_reg_addr, config_done = unpack(SPI_DATA_FMT, raw)
    if old_reg_addr is not reg_addr:
        raise ValueError("Could not read SPI register!")
    else:
        return reg_val


def set_spi_control(roach, zdok_n, adcmode=8, stdby=0, dmux=1, bg=1, bdw=3, fs=0, test=0):
    """
    Sets the control register of an ADC over SPI.
    
    Default mode is DMUX=1:1, gray-code, and channel A only.

    See: http://www.e2v.com/e2v/assets/File/documents/broadband-data-converters/doc0846I.pdf
         specifically section 8.7.3 for more details and possible values.
    """
    reg_val = adcmode + (stdby<<4) + (dmux<<6) + (bg<<7) + (bdw<<8) + (fs<<10) + (test<<12)
    set_spi_register(roach, zdok_n, CONTROL_REG_ADDR, reg_val)


def get_spi_control(roach, zdok_n):
    """
    Gets the current value of the control register of an ADC over SPI.
    
    See: http://www.e2v.com/e2v/assets/File/documents/broadband-data-converters/doc0846I.pdf
         specifically section 8.7.3 for more details and an explanation of returned values.
    """
    #reg_val = adcmode + (stdby<<4) + (dmux<<6) + (bg<<7) + (bdw<<8) + (fs<<10) + (test<<12)
    reg_val = get_spi_register(roach, zdok_n, CONTROL_REG_ADDR-0x80)
    return {'adcmode' : reg_val & 0xf,
            'stdby'   : (reg_val>>4) & 0x3,
            'dmux'    : (reg_val>>6) & 0x1,
            'bg'      : (reg_val>>7) & 0x1,
            'bdw'     : (reg_val>>8) & 0x3,
            'fs'      : (reg_val>>10) & 0x1,
            'test'    : (reg_val>>12) & 0x1}


def set_spi_offset(roach, zdok_n, chan, offset):
    """
    Sets the offset value of one of the four channels on an ADC over SPI.
    
    The offset is a float ranging from -50 mV to +50 mV (with a resolution of 0.4 mV).

    See: http://www.e2v.com/e2v/assets/File/documents/broadband-data-converters/doc0846I.pdf
         specifically section 8.7.14 for more details and possible values.
    """
    reg_val = floor(0.5 + offset*(255/100.)) + 0x80
    set_spi_register(roach, zdok_n, CHANSEL_REG_ADDR, chan)
    set_spi_register(roach, zdok_n, EXTOFFS_REG_ADDR, reg_val)
    set_spi_register(roach, zdok_n, CALCTRL_REG_ADDR, 2<<2)


def get_spi_offset(roach, zdok_n, chan):
    """
    Sets the offset value of one of the four channels on an ADC over SPI.
    
    The offset is a float ranging from -50 mV to +50 mV (with a resolution of 0.4 mV).

    See: http://www.e2v.com/e2v/assets/File/documents/broadband-data-converters/doc0846I.pdf
         specifically section 8.7.14 for more details and possible values.
    """
    set_spi_register(roach, zdok_n, CHANSEL_REG_ADDR, chan)
    reg_val = get_spi_register(roach, zdok_n, EXTOFFS_REG_ADDR-0x80)
    return (reg_val - 0x80)*(100./255)


def set_spi_gain(roach, zdok_n, chan, gain):
    """
    Sets the gain value of one of the four channels on an ADC over SPI.
    
    The gain is a float ranging from -18% to +18% (with a resolution of 0.14%).

    See: http://www.e2v.com/e2v/assets/File/documents/broadband-data-converters/doc0846I.pdf
         specifically section 8.7.16 for more details and possible values.
    """
    reg_val = floor(0.5+gain*(255/36.)) + 0x80
    set_spi_register(roach, zdok_n, CHANSEL_REG_ADDR, chan)
    set_spi_register(roach, zdok_n, EXTGAIN_REG_ADDR, reg_val)
    set_spi_register(roach, zdok_n, CALCTRL_REG_ADDR, 2<<4)


def get_spi_gain(roach, zdok_n, chan):
    """
    Gets the gain value of one of the four channels on an ADC over SPI.
    
    The gain is a float ranging from -18% to +18% (with a resolution of 0.14%).

    See: http://www.e2v.com/e2v/assets/File/documents/broadband-data-converters/doc0846I.pdf
         specifically section 8.7.16 for more details and possible values.
    """
    set_spi_register(roach, zdok_n, CHANSEL_REG_ADDR, chan)
    reg_val = get_spi_register(roach, zdok_n, EXTGAIN_REG_ADDR-0x80)
    return (reg_val - 0x80)*(36./255)


def set_spi_phase(roach, zdok_n, chan, phase):
    """
    Sets the phase value of one of the four channels on an ADC over SPI.
    
    The phase is a float ranging from -14 ps to +14 ps (with a resolution of 110 fs).

    See: http://www.e2v.com/e2v/assets/File/documents/broadband-data-converters/doc0846I.pdf
         specifically section 8.7.18 for more details and possible values.
    """
    reg_val = floor(0.5+phase*(255/28.)) + 0x80
    set_spi_register(roach, zdok_n, CHANSEL_REG_ADDR, chan)
    set_spi_register(roach, zdok_n, EXTPHAS_REG_ADDR, reg_val)
    set_spi_register(roach, zdok_n, CALCTRL_REG_ADDR, 2<<6)


def get_spi_phase(roach, zdok_n, chan):
    """
    Gets the phase value of one of the four channels on an ADC over SPI.
    
    The phase is a float ranging from -14 ps to +14 ps (with a resolution of 110 fs).

    See: http://www.e2v.com/e2v/assets/File/documents/broadband-data-converters/doc0846I.pdf
         specifically section 8.7.18 for more details and possible values.
    """
    set_spi_register(roach, zdok_n, CHANSEL_REG_ADDR, chan)
    reg_val = get_spi_register(roach, zdok_n, EXTPHAS_REG_ADDR-0x80)
    return (reg_val - 0x80)*(28./255)

def set_inl_registers(roach, zdok_n, chan, offs):
    """
    Sets the Integral NonLinearity bits of one of the four channels on
    an ADC over SPI.
    
    The bits are packed into six 16-bit registers on the adc in a way
    that must make sense to the hardware designer. This subroutine takes
    its arguments in a way that is easier to explain

    The argument offs should be a list or array of 17 floats containing
    the fraction of an lsb to offset the adc's reference ladder at 0,
    16, ... 240, 255.  The possible offsets are 0, +-0.15, +-0.3, +-0.45
    and +-0.0.  The values given will be rounded to the nearest of these
    values and converted to the bits in the hardware registerd.

    See: http://www.e2v.com/e2v/assets/File/documents/broadband-data-converters/doc0846I.pdf,
     specifically section 8.7.19 through 8.8, for more details.
    """
    level_to_bits = array([5,4,6,1,0,2,9,8,10])

    # create a array of 6 ints to hold the values for the registers
    regs = zeros((6), dtype='int32')
    r = 2	# r is the relative register number.  R = 2 for 0x32 and 0x35
    regbit = 8 #  regbit is the bit in the register
    for level in range(17):	# n is the bit number in the incoming bits aray
        n = int(floor(0.5 + offs[level]/0.15))
        if n > 4:
            n = 4
        if n < -4:
            n = -4
	i = level_to_bits[4-n]
	regs[r] |= ((i >>2) & 3) << regbit
	regs[r + 3] |= (i & 3)<< regbit
	if regbit == 14:
	    r -= 1
	    regbit = 0
	else:
	    regbit += 2

    set_spi_register(roach, zdok_n, CHANSEL_REG_ADDR, chan)
    for n in range(6):
	reg_val = float(regs[n])
        set_spi_register(roach, zdok_n, FIRST_EXTINL_REG_ADDR+n, reg_val)
    set_spi_register(roach, zdok_n, CALCTRL_REG_ADDR, 2)

def get_inl_registers(roach, zdok_n, chan):
    bits_to_off = array([0,1,-1,0,3,4,2,0,-3,-2,-4])
    offs = zeros((17), dtype = float)
    regs = zeros((6), dtype='int32')

    set_spi_register(roach, zdok_n, CHANSEL_REG_ADDR, chan)
    for n in range(6):
        regs[n] = get_spi_register(roach, zdok_n, FIRST_EXTINL_REG_ADDR-0x80+n)

    r = 2	# r is the relative register number.  R = 2 for 0x32 and 0x35
    regbit = 8	#  regbit is the bit in the register
    for level in range(17):	# n is the bit number in the incoming bits aray
        bits = 0xc & ((regs[r]>>regbit)<<2) | 3 & (regs[r+3]>>regbit)
	offs[level] = 0.15 * bits_to_off[bits]
	if regbit == 14:
	    regbit = 0
	    r -= 1
	else:
	    regbit += 2
    return offs

def inc_mmcm_phase(roach, zdok_n, inc=1):
    """
    This increments (or decrements) the MMCM clk-to-data phase relationship by 
    (1/56) * Pvco, where VCO is depends on the MMCM configuration.

    inc_mmcm_phase(roach, zdok_n)        # default increments
    inc_mmcm_phase(roach, zdok_n, inc=0) # set inc=0 to decrement
    """
    reg_val = pack(SPI_DATA_FMT, (1<<(zdok_n*4)) + (inc<<(1+zdok_n*4)), 0x0, 0x0)
    #print reg_val.encode('string_escape')
    roach.blindwrite(ADC_CONTROLLER, reg_val, offset=0x0)


def get_snapshot(roach, snap_name, bitwidth=8):
    """
    Reads a one-channel snapshot off the given 
    ROACH and returns the time-ordered samples.
    """

    grab = roach.snapshot_get(snap_name)
    data = unpack('%iB' %grab['length'], grab['data'])

    return array(data)


def get_test_vector(roach, zdok_n, snap_name, bitwidth=8):
    """
    Sets the ADC to output a test ramp and reads off the ramp,
    one per core. This should allow a calibration of the MMCM
    phase parameter to reduce bit errors.

    core_a, core_c, core_b, core_d = get_test_vector(roach, zdok_n, snap_name)

    NOTE: This function requires the ADC to be in "test" mode, please use 
    set_spi_control(roach, zdok_n, test=1) before-hand to be in the correct 
    mode.
    """
    data = get_snapshot(roach, snap_name, bitwidth)
    data_bin = (data>>1) ^ data
    return data_bin[0::4], data_bin[1::4], data_bin[2::4], data_bin[3::4]


def calibrate_mmcm_phase(roach, zdok_n, snap_name, bitwidth=8):
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
        core_a, core_c, core_b, core_d = get_test_vector(roach, zdok_n, snap_name)
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
        raise ValueError("No optimal MMCM phase found!")
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
