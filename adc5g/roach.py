"""
This implements some useful classes that mimic corr.katcp_wrapper.FpgaClient

Note: some of this code is copied directly from Jason Manley's Corr package!
"""
import os
from time import time, sleep
from struct import pack, unpack
from subprocess import Popen, PIPE
from signal import SIGKILL


class DummyRoachClient:

    def __init__(self):
        """ 
        A dummy ROACH client for testing.
        """
        self._boffiles = ['fake.bof', 'also_fake.bof']
        self._devices = {"fake_reg": pack('>I', 0), "fake_bram": pack('>1024I', *list([0,]*1024))}

    def blindwrite(self, device_name, data, offset=0):
        """
        Fake Unchecked data write.
        """
        # The following code is taken directly from
        # corr.katcp_wrapper.FpgaClient.blindwrite
        # --->
        assert (type(data)==str) , 'You need to supply binary packed string data!'
        assert (len(data)%4) ==0 , 'You must write 32bit-bounded words!'
        assert ((offset%4) ==0) , 'You must write 32bit-bounded words!'
        # <---
        if device_name not in self._devices:
            raise RuntimeError("Request write failed, i.e. device not present!")
        else:
            old_val = self._devices[device_name]
            old_size = len(self._devices[device_name])
            new_val = old_val[:offset] + data + old_val[offset+len(data):]
            self._devices[device_name] = new_val[:old_size]

    def est_brd_clk(self):
        """
        Returns a fake clock frequency.
        """
        sleep(2)
        return 133.7

    def is_connected(self):
        """
        Always True!
        """
        return True

    def listbof(self):
        """
        Returns a fake list of bof files.
        """
        return self._boffiles

    def listdev(self):
        """
        Returns a fake list of devices.
        """
        return list(self._devices.keys())

    def progdev(self, boffile):
        """
        Waits 5 seconds, then returns 'ok'.
        """
        sleep(5)
        if boffile in self._boffiles:
            return 'ok'
        else:
            raise RuntimeError("Request progdev failed.")

    def ping(self):
        """
        Always True!
        """
        return True

    def read(self, device_name, size, offset=0):
        """
        Read fake data from a device.
        """
        if offset+size > len(self._devices[device_name]):
            raise RuntimeError("Request read failed, offset+size goes beyond device!")
        try:
            return self._devices[device_name][offset:offset+size]
        except KeyError:
            raise RuntimeError("Request read failed, device not presnet!")

    def read_uint(self, device_name, offset=0):
        """
        Read an unsigned integer from a device.
        """
        return unpack('>I', self.read(device_name, 4, offset=offset))[0]

    def write_uint(self, device_name, uint, offset=0):
        """
        Blindly write an unsigned integer to a device.
        """
        self.blindwrite(device_name, pack('>I', uint))


class LocalRoachClient(DummyRoachClient):

    def __init__(self):
        """
        This is a simplified mimicry of corr.katcp_wrapper.FpgaClient
        that runs locally on a ROACH or ROACH2 PowerPC without the need
        for tcpborphserver.
        """
        self._find_proc()
        if os.path.exists(os.path.join(os.getcwd(), 'boffiles')):
            self._bofpath = os.path.join(os.getcwd(), 'boffiles')
        elif os.path.exists(os.path.join(os.path.sep, 'boffiles')):
            self._bofpath = os.path.join(os.path.sep, 'boffiles')
        self.listbof()

    def _find_proc(self):
        """
        Find a running hardware process and set up devices.
        """
        found_proc = False
        for path, dirs, files in os.walk(os.path.join(os.path.sep, 'proc')):
            split_path = path.split(os.path.sep)
            if split_path[-1] == 'ioreg':
                self._pid = int(split_path[2])
                self._devices_path = path
                self._devices = dict((k, '') for k in files)
                found_proc = True
        if not found_proc:
            self._pid = None
            self._devices = {}
            self._devices_path = ''

    def blindwrite(self, device_name, data, offset=0):
        """
        Blindly write data to the device's file handle.
        """
        DummyRoachClient.blindwrite(self, device_name, data, offset)
        device_file = os.path.join(self._devices_path, device_name)
        with open(device_file, 'wb') as file_:
            file_.seek(offset)
            file_.write(data)

    # The following code is basically identical to
    # corr.katcp_wrapper.FpgaClient.est_brd_clk
    # --->
    def est_brd_clk(self):
        """
        Estimate the clock rate of the FPGA fabric.
        """
        first_pass = self.read_uint('sys_clkcounter')
        sleep(2)
        second_pass = self.read_uint('sys_clkcounter')
        if first_pass > second_pass:
            second_pass = second_pass + 2**32
        return (second_pass - first_pass)/2000000.

    def listbof(self):
        """
        Returns the local list of bof files.
        """
        self._boffiles = os.listdir(self._bofpath)
        return DummyRoachClient.listbof(self)

    def progdev(self, boffile):
        """
        Programs the FPGA with the specified bitcode.
        """
        retcode = DummyRoachClient.progdev(self, boffile)
        if self._pid:
            os.kill(self._pid, SIGKILL)
        proc = Popen([os.path.join(self._bofpath, boffile)])
        self._find_proc()
        return retcode

    def read(self, device_name, size, offset=0):
        """
        Read data from a device.
        """
        device_file = os.path.join(self._devices_path, device_name)
        with open(device_file, 'rb') as file_:
            file_.seek(offset)
            data = file_.read(size)
        return data

    # The following code is a modified form of
    # corr.katcp_wrapper.FpgaClient.snapshot_get
    # --->
    def snapshot_get(self, device_name, man_trig=False, man_valid=False, wait_period=1, 
                     offset=-1, circular_capture=False, get_extra_val=False):
        """
        Fake capture of snapshot data.
        """
        self.write_uint(device_name+'_ctrl', (0 + (man_trig<<1) + (man_valid<<2) + (circular_capture<<3)))
        self.write_uint(device_name+'_ctrl', (1 + (man_trig<<1) + (man_valid<<2) + (circular_capture<<3)))

        done = False
        start = time()
        while not done and ((time()-start)<wait_period or (wait_period<0)):
            addr = self.read_uint(device_name + '_status')
            done = not bool(addr & 0x80000000)
            sleep(0.05)

        bram_size = addr & 0x7fffffff
        bram_dmp = dict()
        bram_dmp['length'] = bram_size
        if (bram_size != self.read_uint(device_name+'_status')&0x7fffffff) or bram_size==0:
            raise RuntimeError("A snap block logic error occurred or it didn't finish "
                               "capturing in the allotted %2.2f seconds. Reported %i "
                               "bytes captured."%(wait_period, bram_size))

        if bram_size == 0:
            bram_dmp['data'] = []
        else:
            bram_dmp['data'] = self.read(device_name+'_bram', bram_size)
        return bram_dmp
    # <---
