import unittest
from optparse import OptionParser

import adc5g
try:
    from corr import katcp_wrapper
    REMOTE_POSSIBLE = True
except ImportError:
    REMOTE_POSSIBLE = False


class TestBasics(unittest.TestCase):
    
    def setUp(self):
        if not hasattr(self, '_roach'):
            global roach
            self._roach = roach

    def test_connected(self):
        "testing katcp connectivity"
        self.assertTrue(self._roach.is_connected())

    def test_ping(self):
        "testing katcp pingability"
        self.assertTrue(self._roach.ping())

    
def run_tests(verbosity):
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBasics)
    unittest.TextTestRunner(verbosity=verbosity).run(suite)


def main():
    global roach
    parser = OptionParser()
    parser.add_option("-v", action="store_true", dest="verbose",
                      help="use verbose output while testing")
    parser.add_option("-r", "--remote",
                      dest="remote", metavar="HOST:PORT",
                      help="run tests remotely over katcp using HOST and PORT")
    (options, args) = parser.parse_args()
    if options.verbose:
        verbosity = 2
    else:
        verbosity = 1
    if options.remote:
        if REMOTE_POSSIBLE:
            roach = katcp_wrapper.FpgaClient(*options.remote.split(':'))
        else:
            raise ImportError("corr package was not found, "
                              "remote operation not possible!")
    else:
        roach = adc5g.LocalRoachClient()
    run_tests(verbosity)


if __name__ == "__main__":
    main()
