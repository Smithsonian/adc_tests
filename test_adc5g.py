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


class TestBasics(unittest.TestCase):
    
    def setUp(self):
        if not hasattr(self, '_roach'):
            global roach, boffile
            self._roach = roach
            self._dut = boffile

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

    
def run_tests(verbosity):
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBasics)
    unittest.TextTestRunner(verbosity=verbosity).run(suite)


def main():
    global roach, boffile
    parser = OptionParser()
    parser.add_option("-v", action="store_true", dest="verbose",
                      help="use verbose output while testing")
    parser.add_option("-r", "--remote",
                      dest="remote", metavar="HOST:PORT",
                      help="run tests remotely over katcp using HOST and PORT")
    parser.add_option("-b", "--boffile",
                      dest="boffile", metavar="BOFFILE", default="adc5g_r2_snap.bof",
                      help="test using the BOFFILE bitcode")
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
    run_tests(verbosity)


if __name__ == "__main__":
    main()
