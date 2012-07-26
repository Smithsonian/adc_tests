import unittest
from optparse import OptionParser

import adc5g
try:
    from corr import katcp_wrapper
    REMOTE_POSSIBLE = True
except ImportError:
    REMOTE_POSSIBLE = False


class BofList(list):

    def __repr__(self):
        size = self.__len__()
        return "%d available BOF files" % size


class DevList(list):

    def __repr__(self):
        size = self.__len__()
        return "%d software accessible devices" % size


class ADC5GTestResult(unittest.TextTestResult):

    def getDescription(self, test):
        return test.shortDescription()


class TestBase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        global roach, boffile, zdok_n, clk_rate
        cls._roach = roach
        cls._dut = boffile
        cls._zdok_n = zdok_n
        cls._clk_rate = clk_rate


class TestSetup(TestBase):

    def test_connected(self):
        "testing roach connectivity"
        self.assertTrue(self._roach.is_connected())

    def test_ping(self):
        "testing roach pingability"
        self.assertTrue(self._roach.ping())

    def test_listbof(self):
        "checking if requested bof is available"
        bofs = BofList(self._roach.listbof())
        self.assertIn(self._dut, bofs)


class TestProgramming(TestBase):

    def test_progdev(self):
        "programming the requested bof"
        ret = self._roach.progdev(self._dut)
        self.assertEqual(ret, "ok")


class TestBasics(TestBase):

    @classmethod
    def setUpClass(cls):
        TestBase.setUpClass()
        cls._devices = DevList(cls._roach.listdev())

    def test_clk_rate(self):
        "estimated clock rate should be within 1 MHz of expected"
        rate = self._roach.est_brd_clk()
        self.assertLess(rate, self._clk_rate/8. + 1.0)
        self.assertGreater(rate, self._clk_rate/8. - 1.0)

    def test_has_adc_controller(self):
        "confirm the design has the ADC SPI controller"
        self.assertIn('adc5g_controller', self._devices)

    def test_has_scope(self):
        "confirm the design has the needed scope"
        self.assertIn('raw_%d_bram' % self._zdok_n, self._devices)
        self.assertIn('raw_%d_ctrl' % self._zdok_n, self._devices)
        self.assertIn('raw_%d_status' % self._zdok_n, self._devices)

        
class TestCalibration(TestBase):

    @classmethod
    def setUpClass(cls):
        TestBase.setUpClass()
        cls._optimal_phase, cls._glitches = adc5g.calibrate_mmcm_phase(
            cls._roach, cls._zdok_n, 'raw_%d' % cls._zdok_n)

    def test_optimal_solution_found(self):
        "testing if calibration found optimal MMCM phase"
        self.assertIsNotNone(self._optimal_phase)


class TestInitialSPIControl(TestBase):

    @classmethod
    def setUpClass(cls):
        TestBase.setUpClass()
        cls._control_dict = adc5g.get_spi_control(cls._roach, cls._zdok_n)

    def assertControlParameterIs(self, param, value, msg=None):
        if self._control_dict[param] != value:
            standardMsg = "Control parameter '%s' is not %r" % (param, value)
            self.fail(self._formatMessage(msg, standardMsg))

    def test_adc_mode(self):
        "mode must be single-input A"
        self.assertControlParameterIs('adcmode', 8)

    def test_bandwidth(self):
        "bandwidth should be set to full 2 GHz"
        self.assertControlParameterIs('bdw', 3)

    def test_gray_code(self):
        "gray code should be enabled"
        self.assertControlParameterIs('bg', 1)

    def test_demux(self):
        "demux should be 1:1"
        self.assertControlParameterIs('dmux', 1)

    def test_full_scale(self):
        "full scale should be set to 500 mVpp"
        self.assertControlParameterIs('fs', 0)

    def test_standby(self):
        "board should be in full active mode"
        self.assertControlParameterIs('stdby', 0)

    def test_ramp_disabled(self):
        "ramp mode is off"
        self.assertControlParameterIs('test', 0)


def run_tests(verbosity):
    loader = unittest.TestLoader()
    full_suite = unittest.TestSuite([
            loader.loadTestsFromTestCase(TestSetup),
            loader.loadTestsFromTestCase(TestProgramming),
            loader.loadTestsFromTestCase(TestBasics),
            loader.loadTestsFromTestCase(TestCalibration),
            loader.loadTestsFromTestCase(TestInitialSPIControl),
            ])
    runner = unittest.TextTestRunner(
        verbosity=verbosity, failfast=True, resultclass=ADC5GTestResult)
    runner.run(full_suite)


def main():
    global roach, boffile, zdok_n, clk_rate
    parser = OptionParser()
    parser.add_option("-v", action="store_true", dest="verbose",
                      help="use verbose output while testing")
    parser.add_option("-r", "--remote",
                      dest="remote", metavar="HOST:PORT",
                      help="run tests remotely over katcp using HOST and PORT")
    parser.add_option("-b", "--boffile",
                      dest="boffile", metavar="BOFFILE", default="digicom_r2_2500MHz.bof",
                      help="test using the BOFFILE bitcode")
    parser.add_option("-z", "--zdok",
                      dest="zdok_n", metavar="ZDOK", type='int', default=0,
                      help="test the ADC in the ZDOK port")
    parser.add_option("-c", "--clk-rate",
                      dest="clk_rate", metavar="CLK_MHZ", type='float', default=2500.0,
                      help="specify the input clock frequency in MHz")
    (options, args) = parser.parse_args()
    if options.verbose:
        verbosity = 2
    else:
        verbosity = 1
    if options.remote:
        if REMOTE_POSSIBLE:
            roach = katcp_wrapper.FpgaClient(*options.remote.split(':'))
            roach.wait_connected(1)
        else:
            raise ImportError("corr package was not found, "
                              "remote operation not possible!")
    else:
        roach = adc5g.LocalRoachClient()
    boffile = options.boffile
    zdok_n = options.zdok_n
    clk_rate = options.clk_rate
    run_tests(verbosity)


if __name__ == "__main__":
    main()
