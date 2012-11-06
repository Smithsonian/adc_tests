from math import floor
from struct import pack, unpack
from numpy import array, zeros
from opb import (
    OPB_CONTROLLER,
    OPB_DATA_FMT,
    )

CONTROL_REG_ADDR = 0x01 + 0x80
CHANSEL_REG_ADDR = 0x0f + 0x80
EXTOFFS_REG_ADDR = 0x20 + 0x80
EXTGAIN_REG_ADDR = 0x22 + 0x80
EXTPHAS_REG_ADDR = 0x24 + 0x80
CALCTRL_REG_ADDR = 0x10 + 0x80
FIRST_EXTINL_REG_ADDR = 0x30 + 0x80


def set_spi_register(roach, zdok_n, reg_addr, reg_val):
    """
    Sets the value of an ADC's register over SPI
    """
    spi_data = pack(OPB_DATA_FMT, reg_val, reg_addr, 0x01)
    roach.blindwrite(OPB_CONTROLLER, spi_data, offset=0x4+zdok_n*0x4)


def get_spi_register(roach, zdok_n, reg_addr):
    """
    Gets the value of an ADC's register over SPI
    """
    spi_data = pack(OPB_DATA_FMT, 0x0, reg_addr, 0x01)
    roach.blindwrite(OPB_CONTROLLER, spi_data, offset=0x4+zdok_n*0x4)
    raw = roach.read(OPB_CONTROLLER, 0x4, offset=0x4+zdok_n*0x4)
    reg_val, old_reg_addr, config_done = unpack(OPB_DATA_FMT, raw)
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
