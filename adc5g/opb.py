from struct import pack


OPB_CONTROLLER = 'adc5g_controller'
OPB_DATA_FMT = '>H2B'


def inc_mmcm_phase(roach, zdok_n, inc=1):
    """
    This increments (or decrements) the MMCM clk-to-data phase relationship by 
    (1/56) * Pvco, where VCO is depends on the MMCM configuration.

    inc_mmcm_phase(roach, zdok_n)        # default increments
    inc_mmcm_phase(roach, zdok_n, inc=0) # set inc=0 to decrement
    """
    reg_val = pack(OPB_DATA_FMT, (1<<(zdok_n*4)) + (inc<<(1+zdok_n*4)), 0x0, 0x0)
    roach.blindwrite(OPB_CONTROLLER, reg_val, offset=0x0)
