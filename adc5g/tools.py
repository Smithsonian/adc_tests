from math import pi
from time import sleep
from struct import pack, unpack
from matplotlib.mlab import (
    psd, detrend_mean,
)
from numpy import (
    array, sinc, sum, max, log10, floor,
)


SPI_DATA_FMT = '>H2B'
CONTROL_REG_ADDR = 0x01 + 0x80
CHANSEL_REG_ADDR = 0x0f + 0x80
EXTOFFS_REG_ADDR = 0x20 + 0x80
EXTGAIN_REG_ADDR = 0x22 + 0x80
EXTPHAS_REG_ADDR = 0x24 + 0x80
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
    print raw.encode('string_escape')
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
    reg_val = floor(offset*(255/100.)) + 0x80
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
    reg_val = floor(gain*(255/36.)) + 0x80
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
    reg_val = floor(phase*(255/28.)) + 0x80
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


def get_snapshot(roach, snap_name, bitwidth=8):
    """
    Reads a one-channel snapshot off the given 
    ROACH and returns the time-ordered samples.
    """

    grab = roach.snapshot_get(snap_name)
    data = unpack('%iB' %grab['length'], grab['data'])

    return array(data)


def get_psd(roach, snap_name, samp_freq, bitwidth=8, NFFT=256):
    """
    Reads data off a given channel on a ROACH and calculates
    the power spectral density of the time-series.
    """

    data = get_snapshot(roach, snap_name, bitwidth)
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
