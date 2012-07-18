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


class ADC5GTestResult(unittest.TextTestResult):

    def getDescription(self, test):
        return test.shortDescription()


class TestBase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        global roach, boffile, zdok_n
        cls._roach = roach
        cls._dut = boffile
        cls._zdok_n = zdok_n


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


class TestCalibration(TestBase):

    @classmethod
    def setUpClass(cls):
        TestBase.setUpClass()
        cls._optimal_phase, cls._glitches = adc5g.calibrate_mmcm_phase(
            cls._roach, cls._zdok_n, 'raw_%d' % cls._zdok_n)

    def test_optimal_solution_found(self):
        "testing if calibration found optimal MMCM phase"
        self.assertIsNotNone(self._optimal_phase)


class TestBasics(TestBase):
    pass


def run_tests(verbosity):
    loader = unittest.TestLoader()
    full_suite = unittest.TestSuite([
            loader.loadTestsFromTestCase(TestSetup),
            loader.loadTestsFromTestCase(TestProgramming),
            loader.loadTestsFromTestCase(TestCalibration),
            ])
    runner = unittest.TextTestRunner(
        verbosity=verbosity, failfast=True, resultclass=ADC5GTestResult)
    runner.run(full_suite)


def main():
    global roach, boffile, zdok_n
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
    run_tests(verbosity)


if __name__ == "__main__":
    main()
