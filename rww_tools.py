# run 'ipython --pylab -i'
# and then in ipython: '%run -i rww_tools.py'
# the -i there causes rww_tools.py to be run in the interactive name space
import sys
import os
import time
from corr import katcp_wrapper
from corr.snap import snapshots_get
from struct import pack, unpack
from scipy.signal import correlate
from scipy.fftpack import fft, rfft
#import katcp_wrapper
roach2=katcp_wrapper.FpgaClient('roach2-00.cfa.harvard.edu')
zdok=0
import adc5g
import matplotlib.pyplot as plt
from matplotlib import mlab
import numpy as np
from numpy import math
import fit_cores
lanio = "lanio 131.142.9.146 "
DEFAULT_BITCODE = 'sma_corr_2014_Apr_21_1603.bof.gz'

freq = 10.070801
pwr = 1.0
numpoints=16384
#samp_freq = 5000.0
#snap_name = "scope_raw_0_snap"
#prog_name = 'adc5g_test.bof'

def dosnap(fr=0, name=None, rpt = 1, donot_clear=False, plot=True):
  """
  Takes a snapshot and uses fit_cores to fit a sine function to each
  core separately assuming a CW signal is connected to the input.  The
  offset, gain and phase differences are reoprted for each core as
  well as the average of all four.

  The parameters are:
    fr   The frequency of the signal generator.  It will default to the last
         frequency set by set_freq()
    name the name of the file into which the snapshot is written.  5 other
         files are written.  Name.c1 .. name.c4 contain themeasurements from
	 cores a, b, c and d.  Note that data is taken from cores in the order
	 a, c, b, d.  A line is appended to the file name.fit containing
	 signal freq, average zero, average amplitude followed by triplets
	 of zero, amplitude and phase differences for cores a, b, c and d
	 name defaults to if0 or if1, depending on the current zdok

    rpt  The number of repeats.  Defaults to 1.  The c1 .. c4 files mentioned
         above are overwritten with each repeat, but new rows of data are added
	 to the .fit file for each pass.
  """
  global freq
  if name == None:
    name = "if%d" % (zdok)
  avg_pwr_sinad = 0
  if fr == 0:
    fr = freq
  for i in range(rpt):
    snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
    if(plot):
      plt.clf()
      plt.plot(snap)
      plt.show(block = False)
      rmsSnap = np.std(snap)
      loadingFactor = -20.0*math.log10(128/rmsSnap)
      print "Rms = %f, loading factor = %f" % (rmsSnap,loadingFactor)
    if i == rpt-1:
      np.savetxt(name, snap,fmt='%d')
    ogp, pwr_sinad = fit_cores.fit_snap(snap, fr, samp_freq, name,\
       clear_avgs = ((i == 0) and not donot_clear), prnt = (i == rpt-1))
    avg_pwr_sinad += pwr_sinad
  return ogp, avg_pwr_sinad/rpt

def get_rms():
  global zdok

  save_zdok = zdok
  set_zdok(0)
  snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
  rmsSnap = np.std(snap)
  loadingFactor = -20.0*math.log10(128/rmsSnap)
  print "If0 Rms = %f, loading factor = %f" % (rmsSnap,loadingFactor)
  set_zdok(1)
  snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
  rmsSnap = np.std(snap)
  loadingFactor = -20.0*math.log10(128/rmsSnap)
  print "If1 Rms = %f, loading factor = %f" % (rmsSnap,loadingFactor)
  if zdok != save_zdok:
    set_zdok(save_zdok)

def dosim(freq=10.070801, name="sim", rpt = 1, exact=True):
  """
  Do the same analysis as dosnap, but on simulated data.
  The arguments have the same meaning as dosnap except that freq is not
  coupled to the global variable.  A random phase is generated for each pass
  """

  for i in range(rpt):
    snap=get_sim_data(freq, exact)
    if i == rpt-1:
      np.savetxt(name, snap,fmt='%d')
    fit_cores.fit_snap(snap, freq, samp_freq, name, i == 0)

def simpsd(freq=318.0, rpt = 1, exact=True):
  """
  Make a simulated snapshot and do the psd analysis on it.  The
  sine wave will have a randon start phase.
  """
  for i in range(rpt):
    data = get_sim_data(freq, exact)
    power, freqs = mlab.psd(data, numpoints, Fs=samp_freq*1e6, \
        detrend=mlab.detrend_mean, scale_by_freq=True)
    plt.clf()
    if i == 0:
      sp = power
    else:
      sp += power
  sp /= rpt
  print "about to plot", len(freqs)
  plt.step(freqs, 10*np.log10(sp))
  plt.show(block = False)
  fd = open("sim.psd", 'w')
  for i in range(len(sp)):
    print >>fd, "%7.2f %6.1f" % (freqs[i]/1e6, 10*np.log10(sp[i]))

def get_sim_data(freq, exact=True):
  """
  Make a simulated snapshot of data
  """
  if exact:
    offs = [0,0,0,0]
    gains = [1,1,1,1]
  else:
    offs = [.2, .3, -.2, -.1]
    gains = [1.001, .9984, .999, 1.002]
  del_phi = 2 * math.pi * freq / samp_freq
  data = np.empty((numpoints), dtype='int32')
  phase = 2*math.pi * np.random.uniform()
  for n in range(numpoints):
    core = n&3
    data[n] = (math.floor(0.5 + 119.0 * math.sin(del_phi * n + phase) + \
        offs[core]))*gains[core]
  return data

def dotest(plotcore = 1):
  """
  Put the adc in test mode and get a sample of the test vector.  Plot core 1
  by default.
  """
  global snap_name
  adc5g.set_spi_control(roach2, zdok, test=1)
  cores = (corea, corec, coreb, cored) = adc5g.get_test_vector(roach2, [snap_name,])
  if plotcore == 2:
    plotcore = 3
  elif plotcore == 3:
    plotcore = 2
  plt.clf()
  plt.plot(cores[plotcore])
  plt.show(block=False)
  adc5g.set_spi_control(roach2, zdok)

def dopsd(nfft = None, rpt = 10, plotdB=True):
  """
  Takes a snapshot, then computes, plots and writes out the Power Spectral
  Density functions.  The psd function is written into a file named "psd".
  This file will be overwritten with each call.  Arguments:

  nfft The number of points in the psd function.  Defaults to the number of
       points in a snapshot, the maximum which should be used.

  rpt  The numper of mesurements to be averaged for the plot and output file. 
       Defaults to 10.

  plotdB controls whether the plot is linear in power or in dB
  """
  if nfft == None:
    nfft = numpoints
  for i in range(rpt):
    power, freqs = adc5g.get_psd(roach2, snap_name, samp_freq*1e6, 8, nfft)
    if i == 0:
      sp = power
    else:
      sp += power
  sp /= rpt
  if plotdB:
    plt.step(freqs, 10*np.log10(sp))
  else:
    plt.step(freqs, (sp))
  plt.show(block = False)
  data = np.column_stack((freqs/1e6, 10*np.log10(sp)))
  np.savetxt("psd", data, fmt=('%7.2f', '%6.1f'))

def dopsdcores(nfft = None, rpt = 10):
  if nfft == None:
    nfft = numpoints/4
  for i in range(rpt):
    snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
    power, freqs = mlab.psd(snap, nfft*4, Fs=samp_freq*1e6, detrend=mlab.detrend_mean, \
        scale_by_freq=False)
    if i == 0:
      psd_all = power[:1+nfft/2]
    else:
      psd_all += power[:1+nfft/2]
    power, freqs = mlab.psd(snap[0:: 4], nfft, Fs=samp_freq*.25e6, detrend=mlab.detrend_mean,\
        scale_by_freq=True)
    if i == 0:
      psd1 = power
    else:
      psd1 += power
    power, freqs = mlab.psd(snap[1:: 4], nfft, Fs=samp_freq*.25e6, detrend=mlab.detrend_mean,\
        scale_by_freq=True)
    if i == 0:
      psd2 = power
    else:
      psd2 += power
    power, freqs = mlab.psd(snap[2:: 4], nfft, Fs=samp_freq*.25e6, detrend=mlab.detrend_mean,\
        scale_by_freq=True)
    if i == 0:
      psd3 = power
    else:
      psd3 += power
    power, freqs = mlab.psd(snap[3:: 4], nfft, Fs=samp_freq*.25e6, detrend=mlab.detrend_mean,\
        scale_by_freq=True)
    if i == 0:
      psd4 = power
    else:
      psd4 += power
  data = np.column_stack((freqs*1e-6, 10*np.log10(psd_all/rpt), \
      10*np.log10(psd1/rpt), \
      10*np.log10(psd2/rpt), 10*np.log10(psd3/rpt), 10*np.log10(psd4/rpt)))
  np.savetxt("psd_cores", data, fmt=('%7.2f'))

def hist_from_snapshots(rpt = 10):
#  hist_all = np.zeros(256,dtype=int)
  hist1 = np.zeros(256,dtype=int)
  hist2 = np.zeros(256,dtype=int)
  hist3 = np.zeros(256,dtype=int)
  hist4 = np.zeros(256,dtype=int)
  for i in range(rpt):
    snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
    snap = 128 + np.array(snap)
#    hist = np.bincount(snap, minlength=256)
#    hist_all += hist
    hist = np.bincount(snap[0:: 4], minlength=256)
    hist1 += hist
    hist = np.bincount(snap[1:: 4], minlength=256)
    hist2 += hist
    hist = np.bincount(snap[2:: 4], minlength=256)
    hist3 += hist
    hist = np.bincount(snap[3:: 4], minlength=256)
    hist4 += hist
  data=np.column_stack((np.arange(-128., 128, dtype=int), hist1, hist2,
      hist3, hist4))
  np.savetxt("hist_cores", data, fmt=("%d"))
#  print "all ",np.sum(hist_all[0:128]), np.sum(hist_all[128:256])
  print "core a  ",np.sum(hist1[0:128]), np.sum(hist1[129:256])
  print "core b  ",np.sum(hist3[0:128]), np.sum(hist3[129:256])
  print "core c  ",np.sum(hist2[0:128]), np.sum(hist2[129:256])
  print "core d  ",np.sum(hist4[0:128]), np.sum(hist4[129:256])

# For now get_histogram has cores b and c reversed.
def get_hist(fname="hist_cores"):
  data = np.empty(shape=(256,5), dtype=int)
  for c in range(4):
    data[:, c+1] = adc5g.get_histogram(roach2, zdok, "acbd"[c])
  data[:,0] = range(-128, 128)
  np.savetxt(fname, data, fmt=("%d"))

def multifreq(start=100, end=560, step=50, repeat=5, do_sfdr=False):
  """
  Calls dosnap for a range of frequenciesi in MHz.  The actual frequencies are
  picked to have an odd number of cycles in the 16384 point snapshot.
  """
  global ogp_name

  name = "if%d" % (zdok)
  sfd = open('sinad', 'a')
  f = samp_freq / numpoints
  nstart = int(0.5+start/f)
  nend = int(0.5+end/f)
  nstep = int(0.5+step/f)
  for n in range(nstart, nend, nstep):
    freq = f*n
    set_freq(freq)
#    ogp, avg_pwr_sinad = dosnap(rpt=repeat, donot_clear = False)
    ogp, avg_pwr_sinad = dosnap(rpt=repeat, name=name,\
          donot_clear = n!=nstart, plot=False)
    sinad = 10.0*np.log10(avg_pwr_sinad)
    print >>sfd, "%8.3f %7.2f" % (freq, sinad)
    if do_sfdr:
      dopsd(rpt=3)
      fit_cores.dosfdr(freq)
  np.savetxt(ogp_name+".meas", ogp[3:], fmt="%8.4f")
  fit_cores.fit_inl(name+'.res')

def freqResp(start=100, end=2400, delta=50, repeat=10,powerlevel=7):
  """
  Runs adc5g.get_psd for a range of frequenciesi in MHz.  
  The actual frequencies are picked to have an odd number of 
  cycles in the 16384 point snapshot. (This part is copied from multifreq()).
  Writes out freq and max() of spectrum to freqResponse.dat file.
  """
  set_pwr(powerlevel)
  frfile = open('freqResponse.dat', 'a')
  f = samp_freq / numpoints
  nstart = int(0.5+start/f)
  nend = int(0.5+end/f)
  nstep = int(0.5+delta/f)
  for n in range(nstart, nend, nstep):
    freq = f*n
    set_freq(freq)
#    dopsd(rpt=3)
    for i in range(repeat):
      power, freqs = adc5g.get_psd(roach2, snap_name, samp_freq*1e6, 8, numpoints)
      if i == 0:
        sp = power
      else:
        sp += power
      sp /= repeat
    power=10*log10(sp)
#    step(freqs, power)
    peakpower=max(power)
    print freq,peakpower
    output="%f %f\n" % (freq,peakpower)
    frfile.write(output)
  frfile.close()


def multipwr(start = 1, end = -40, step = -3, repeat=10):
  """
  Calls dosnap for a range of powers
  """
  for n in range(start, end, step):
    set_pwr(n)
    dosnap(rpt=repeat, plot=False)

def update_ogp(fname = None, set=True):
  """
  Retreive the ogp data from the ADC and add in the corrections from
  the measured ogp (in ogp.meas).  Store in the standard ogp file
  in/instance/configFiles for the roach.
  """
  global ogp_name

  if fname == None:
    fname = ogp_name+".meas"
  cur_ogp = get_ogp_array()
  meas_ogp = np.genfromtxt(fname)
  # Correct for the ~1.4X larger effect of the phase registers than expected
  for i in (2,5,8,11):
    meas_ogp[i] *= 0.65
  np.savetxt(ogp_name, cur_ogp+meas_ogp, fmt="%8.4f")
  if set:
    set_ogp()

def update_inl(fname = None, set=True):
  """
  Retreive the INL data from the ADC and add in the corrections from
  the measured inl (in inl.meas).  Store in the file 'inl'
  """
  global inl_name

  if fname == None:
    fname = inl_name+".meas"
  cur_inl = get_inl_array()
  meas_inl = np.genfromtxt(fname)
  for level in range(17):
    cur_inl[level][1:] += meas_inl[level][1:]
  np.savetxt(inl_name, cur_inl, fmt=('%3d','%7.4f','%7.4f','%7.4f','%7.4f'))
  if set:
    set_inl()

def program():
  """
  Program the roach2 with the standard program.  After this, calibrate()
  should be called
  """
  global prog_name, roach2, samp_freq

  roach2.progdev(prog_name)
  adc5g.set_spi_control(roach2, zdok)
  if prog_name[:10] != 'adc5g_test':
#    print "set up for correlator code"
    roach2.write_int('source_ctrl', 18)
    roach2.write_int('scope_ctrl', 1536)
    samp_freq = 2288.
    numpoints = 32768
  set_zdok(zdok)

def calibrate(verbose=False):
  """
  Call Rurik's routine to calibrate the time delay at the adc interface.
  """
  global zdok

  adc5g.set_test_mode(roach2, 0)
  adc5g.set_test_mode(roach2, 1)
  adc5g.sync_adc(roach2)
  save_zdok = zdok
  set_zdok(0)
  opt0, glitches0 = adc5g.calibrate_mmcm_phase(roach2, 0, \
      [snap_name,])
  if verbose or (opt0 == None):
    print "zodk0 ", opt0, glitches0
  else:
    print "zodk0", opt0
  set_zdok(1)
  opt1, glitches1 = adc5g.calibrate_mmcm_phase(roach2, 1, \
      [snap_name,])
  if verbose or (opt1 == None):
    print "zodk1 ", opt1, glitches1
  else:
    print "zodk1", opt1
  set_zdok(save_zdok)
  adc5g.unset_test_mode(roach2, 0)
  adc5g.unset_test_mode(roach2, 1)

#  t = adc5g.calibrate_mmcm_phase(roach2, zdok, [snap_name,], bitwidth=8)
#  print t

def clear_ogp():
  """
  Clear all of the Offset, Gain and Phase corrections registers on the adc.
  """
  for core in range(1,5):
    adc5g.set_spi_gain(roach2,zdok, core, 0)
    adc5g.set_spi_offset(roach2,zdok, core, 0)
    adc5g.set_spi_phase(roach2,zdok, core, 0)

def get_ogp():
  """
  Use get_ogp_array to get the Offset, Gain and Phase corrections
  and print them.
  """
  ogp = get_ogp_array()
  print "zero(mV) amp(%%)  dly(ps) (adj by .4, .14, .11)"
  print "core A  %7.4f %7.4f %8.4f" %  (ogp[0], ogp[1], ogp[2])
  print "core B  %7.4f %7.4f %8.4f" %  (ogp[3], ogp[4], ogp[5])
  print "core C  %7.4f %7.4f %8.4f" %  (ogp[6], ogp[7], ogp[8])
  print "core D  %7.4f %7.4f %8.4f" %  (ogp[9], ogp[10], ogp[11])

def set_ogp(fname = None):
  """
  Clear the control register and then load the offset, gain and phase
  registers for each core.  These values are hard coded for now.
  fname defaults to the standard ogp name in /instance/configFiles on the roach.
  """
  global ogp_name

  if fname == None:
    fname = ogp_name
  adc5g.set_spi_control(roach2, zdok)
  t = np.genfromtxt(fname)
  set_offs(t[0], t[3], t[6], t[9])
  set_gains(t[1], t[4], t[7], t[10])
  set_phase(t[2], t[5], t[8], t[11])

def clear_inl():
  """
  Clear the INL registers on teh ADC
  """
  offs = [0.0]*17
  for chan in range(1,5):
    adc5g.set_inl_registers(roach2, zdok, chan, offs)

def get_inl():
  """
  Use get_inl_array to get the INL registers from the ADC and then
  print them.
  """
  a = get_inl_array()
  print "lvl  A     B     C     D"
  for level in range(17):
    print "%3d %5.2f %5.2f %5.2f %5.2f" % tuple(a[level])
    

def set_inl(fname = None):
  """
  Set the INL registers for all four cores from a file containing 17 rows
  of 5 columns.  The first column contains the level and is ignored.
  Columns 2-5 contain the inl correction for cores a-d
  fname defaults to the standard name in /instance/configFiles on the roach
  """
  global inl_name

  if fname == None:
    fname = inl_name
  c = np.genfromtxt(fname, usecols=(1,2,3,4), unpack=True)
  adc5g.set_inl_registers(roach2,zdok,1,c[0])
  adc5g.set_inl_registers(roach2,zdok,2,c[1])
  adc5g.set_inl_registers(roach2,zdok,3,c[2])
  adc5g.set_inl_registers(roach2,zdok,4,c[3])

def set_freq(fr, centered = True, prnt=True):
  """
  Set the synthesizer frequency (MHz) and save the value for use by dosnap(),
  etc.  If centered is True, pick the closest frequency in the center of a
  channel with an odd number of cycles in a snapshot.
  """
  global freq
  if centered:
    base_freq = samp_freq / numpoints
    n = 2*int(fr/(2.0*base_freq))+1
    freq = base_freq*n
    print "n, freq =  ", n, freq
  else:
    freq=fr
  os.system(lanio + "\":FREQ " + str(freq) + " MHz\"")
  if prnt:
    print "%.6f" % (freq)
  time.sleep(0.5)
  
def get_freq():
  """
  Retreive the frequency from the Agilent Synthesizer and print it (in Hz).
  """
  print os.system(lanio + "\"FREQ?\"")

def set_pwr(p):
  """
  Set the synthesizer power and save the value for use by dosnap(), etc.
  """
  global pwr
  pwr = p
  os.system(lanio + "\":POW " + str(p) + " dBm\"")
  os.system(lanio + "\":OUTP 1\"")
  
def get_pwr():
  """
  Retreive the power level from the Agilent Synthesizer and print it.
  """
  print os.system(lanio + "\"POW?\"")

def set_offs(o1, o2, o3, o4):
  """
  Set the offsets for each core in the order a, b, c, d.
  """
  t = float(o1)
  print math.floor(.5+t*255/100.)+0x80,
  adc5g.set_spi_offset(roach2,zdok, 1, t)
  t = float(o2)
  print math.floor(.5+t*255/100.)+0x80,
  adc5g.set_spi_offset(roach2,zdok, 2, t)
  t = float(o3)
  print math.floor(.5+t*255/100.)+0x80,
  adc5g.set_spi_offset(roach2,zdok, 3, t)
  t = float(o4)
  print math.floor(.5+t*255/100.)+0x80
  adc5g.set_spi_offset(roach2,zdok, 4, t)
def get_offs():
  """
  Get and print the offsets for the four cores of the ADC.
  """
  for i in range(1,5):
    print "%.3f " % adc5g.get_spi_offset(roach2,zdok,i),
  print

def set_gains(g1, g2, g3, g4):
  """
  Set the gains for each core in the order a, b, c, d.
  """
  t = float(g1)
  print math.floor(.5+t*255/36.)+0x80,
  adc5g.set_spi_gain(roach2,zdok, 1, t)
  t = float(g2)
  print math.floor(.5+t*255/36.)+0x80,
  adc5g.set_spi_gain(roach2,zdok, 2, t)
  t = float(g3)
  print math.floor(.5+t*255/36.)+0x80,
  adc5g.set_spi_gain(roach2,zdok, 3, t)
  t = float(g4)
  print math.floor(.5+t*255/36.)+0x80
  adc5g.set_spi_gain(roach2,zdok, 4, t)
def get_gains():
  """
  Get and print the gains for the four cores of the ADC.
  """
  for i in range(1,5):
    print "%.3f " % adc5g.get_spi_gain(roach2,zdok,i),
  print
def set_phase(p1, p2, p3, p4):
  """
  Set the phases (delays) for each core in the order a, b, c, d.
  """
  t = float(p1)
  print math.floor(.5+t*255/28.)+0x80,
  adc5g.set_spi_phase(roach2,zdok, 1, t)
  t = float(p2)
  print math.floor(.5+t*255/28.)+0x80,
  adc5g.set_spi_phase(roach2,zdok, 2, t)
  t = float(p3)
  print math.floor(.5+t*255/28.)+0x80,
  adc5g.set_spi_phase(roach2,zdok, 3, t)
  t = float(p4)
  print math.floor(.5+t*255/28.)+0x80
  adc5g.set_spi_phase(roach2,zdok, 4, t)
def get_phase():
  """
  Get and print the delays for the four cores of the ADC.
  """
  for i in range(1,5):
    print "%.3f " % adc5g.get_spi_phase(roach2,zdok,i),
  print

def get_inl_array():
  """
  Read the INL corrections from the adc and put in an array
  """
  inl = np.zeros((5,17), dtype='float')
  for chan in range(1,5):
    inl[chan] = adc5g.get_inl_registers(roach2, zdok, chan)
  inl[0] = range(0, 257,16)
  return inl.transpose()

def get_ogp_array():
  """
  Read  the Offset, Gain and Phase corrections for each core from the ADC
  and return in a 1D array
  """
  ogp = np.zeros((12), dtype='float')
  indx = 0
  for chan in range(1,5):
    ogp[indx] = adc5g.get_spi_offset(roach2,zdok,chan)
    indx += 1
    ogp[indx] = adc5g.get_spi_gain(roach2,zdok,chan)
    indx += 1
    ogp[indx] = adc5g.get_spi_phase(roach2,zdok,chan)
    indx += 1
  return ogp

def set_zdok(zd):
  global zdok, snap_name, prog_name, ogp_name, inl_name
  if prog_name[:8] == 'sma_corr':
    snap_name = "scope_snap%d" % (zd)
  else:
    snap_name = "scope_raw_%d_snap" % (zd)
  zdok = zd
  d = roach_name[-1] if roach_name[-2] == '0' else roach_name[-2:]
  ogp_name = inst_name+"/configFiles/ogp_if%d" % (zd)
  inl_name = inst_name+"/configFiles/inl_if%d" % (zd)

def get_zdok():
  print "zdok %d, snapshot %s" % (zdok, snap_name)

def setup(r_name='roach2-00', prg_nam=DEFAULT_BITCODE):
  global roach_name, inst_name, roach2, prog_name, zdok, samp_freq, numpoints

  roach2=katcp_wrapper.FpgaClient(r_name)
  connected = roach2.wait_connected(timeout=2)
  if connected == False:
    raise RuntimeError("Unable to connect to %s" %(r_name))
#    print "Unable to connect to %s" %(r_name)
#    return connected
  prog_name = prg_nam
  if prog_name[:8] == 'sma_corr':
#    print "this is correlator code"
    samp_freq = 2288.
    numpoints = 32768
  elif prog_name[:10] == 'adc5g_test':
    print "this is adc test code"
    samp_freq = 5000.0
    numpoints = 16384
  else:
    print "I do not recognize the bit code name, treating it as correlator code"
    samp_freq = 2288.
    numpoints = 32768
  roach_name = r_name
  d = r_name[-1] if r_name[-2] == '0' else r_name[-2:]
  inst_name = "/otherInstances/roach2/"+d
  os.chdir(inst_name+"/adcTests")
  set_zdok(zdok)
#  return connected

def og_from_noise(fname=None, rpt=100, printEach=False):
  """
  Take a number of snapshots of noise.  Analyze for offset and gain
  for each core separately.
  """
  global ogp_name

  if fname == None:
    fname = ogp_name+".meas"
  sum_result = np.zeros((15), dtype=float)
  sum_cnt = 0
  for n in range(rpt):
    result = np.zeros((15), dtype=float)
    snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
    if(rpt == 1):
      np.savetxt("t.og_noise", snap,fmt='%d')
    l=float(len(snap))
    snap_off=np.sum(snap)/l
    snap_amp=np.sum(abs(snap-snap_off))/l
    result[0]=snap_off*(-500.0/256.0)
    result[1]=snap_amp
    for core in range(4):
      # This will actually sample the cores in the order A,C,B,D
      # index will fix this up when data is put in the result array
      index=(3,9,6,12)[core]
      c=snap[core::4]
      l=float(len(c))
      off=np.sum(c)/l
      result[index] = off*(-500.0/256.0)
      amp=np.sum(abs(c-off))/l
      result[index+1]= 100.0*(snap_amp-amp)/snap_amp
    sum_result += result
    sum_cnt += 1
    if printEach:
      print "%.4f "*15 % tuple(result)
  sum_result /= sum_cnt
  print "%.4f "*15 % tuple(sum_result)
  np.savetxt(fname, sum_result[3:], fmt="%8.4f")

def phase_curve():
  f= 28/255.
  ofd = open('phasecurve', 'w')
  p = {}
  for i in range(1,5):
     p[i] = adc5g.get_spi_phase(roach2,zdok,i)
  for i in range(-10,11):
    set_phase(p[1]+ f*i,p[2] -f*i,p[3],p[4])
    ogp, gar = dosnap(rpt=5, plot=False)
    print >>ofd, "%.3f %.3f %.3f" % (f*i, ogp[5], ogp[8])
  set_phase(p[1],p[2],p[3],p[4])

def dohist(base_name='hist', type='sin', gethist=True, plt=True):
  hc_name=base_name+'_cores'
  if gethist:
    get_hist(fname=hc_name)
  res = np.empty([5, 256], dtype=float)
  res[0] = np.arange(256, dtype=float)
  z_fact = 500.0/256.0
  (a1,z1), res[1] =fit_cores.fit_hist(1,type, hc_name)
  (a2,z2), res[2] =fit_cores.fit_hist(2,type, hc_name)
  (a3,z3), res[3] =fit_cores.fit_hist(3,type, hc_name)
  (a4,z4), res[4] =fit_cores.fit_hist(4,type, hc_name)
  avamp = (a1+a2+a3+a4)/4.0
  # Reverse the amplitude and zero differences so they can be applied to the
  # offset and gain registers directly.  The phase registers don't need the
  # reversal
  a1p = 100*(avamp -a1)/avamp
  a2p = 100*(avamp -a2)/avamp
  a3p = 100*(avamp -a3)/avamp
  a4p = 100*(avamp -a4)/avamp
  ogp=np.array([z_fact*z1, a1p, 0, z_fact*z2, a2p, 0, z_fact*z3, a3p, 0, \
      z_fact*z4, a4p, 0])
  avz=(z1+z2+z3+z4)*z_fact/4.0
  print "#avg    %7.4f %7.4f %8.4f" %  (ogp[1], avamp, 0)
  print "core A  %7.4f %7.4f %8.4f" %  tuple(ogp[0:3])
  print "core B  %7.4f %7.4f %8.4f" %  tuple(ogp[3:6])
  print "core C  %7.4f %7.4f %8.4f" %  tuple(ogp[6:9])
  print "core D  %7.4f %7.4f %8.4f" %  tuple(ogp[9:12])
  np.savetxt(base_name+"_ogp.meas", ogp, fmt= "%8.4f")
  r_name=base_name+'.res'
  np.savetxt(r_name, np.transpose(res), fmt='%3i %6.3f %6.3f %6.3f %6.3f')
  fit_cores.fit_inl(fname=r_name)
  if plt:
    plotres(r_name)

def plotres(fname="hist.res",title=""):
  
  res = np.genfromtxt(fname, unpack=True)
  plt.clf()
  plt.plot(res[0][1:-1], res[1][1:-1], label='core a')
  plt.plot(res[0][1:-1], res[2][1:-1], label='core b')
  plt.plot(res[0][1:-1], res[3][1:-1], label='core c')
  plt.plot(res[0][1:-1], res[4][1:-1], label='core d')
  plt.legend(loc=0)
  if title != "":
    plt.title(title)
  plt.show(block=False)

def get_snaps():
  global raw0, raw1
  data = snapshots_get([roach2, roach2], ['scope_raw_0_snap', \
      'scope_raw_1_snap'])
  raw0 = np.array(unpack('%ib' % data['lengths'][0], data['data'][0]), \
      dtype=float)
  raw1 = np.array(unpack('%ib' % data['lengths'][1], data['data'][1]), \
      dtype=float)

def set_adc_delay(cntr_chan = 0):
  del0 = 0
  del1 = 0
  if(cntr_chan > 0):
    del1 = cntr_chan
  else:
    del0 = -cntr_chan
  roach2.write('cdelay_ctrl', pack('>I', (1<<31) + (1<<30) + \
       (del1<<15) + del0),0)

def xcorr_snaps(set_delay = True):
  global raw0, raw1, xcn

  if(set_delay):
    set_adc_delay(0)
  get_snaps()
  xcn = xcorr(raw0, raw1, maxlags=100, normed=True)
  maxchan = np.argmax(xcn[1])
  cntr_chan = xcn[0][maxchan]
  if(set_delay):
    set_adc_delay(cntr_chan)
  print  xcn[1][maxchan], cntr_chan

def fx_snaps(n = 10):
  global raw0, raw1, xc, xn
  
  for cnt in range(n):
    get_snaps()
    ft0=fft(raw0)[:8192]
    ft1=fft(raw1)[:8192]
    if cnt == 0:
      xc = ft0*conjugate(ft1)
      ac0 = ft0*conjugate(ft0)
      ac1 = ft1*conjugate(ft1)
#      ac0 = abs(ft0)*abs(ft0)
#      ac1 = abs(ft1)*abs(ft1)
    else:
      xc += ft0*conjugate(ft1)
      ac0 += ft0*conjugate(ft0)
      ac1 += ft1*conjugate(ft1)
#      ac0 += abs(ft0)*abs(ft0)
#      ac1 += abs(ft1)*abs(ft1)
  xn = xc/(sqrt(ac0*ac1))
  np.savetxt('xn.txt', xn.view(float).reshape(-1, 2))
  print np.mean(abs(xn)), "+-", np.std(abs(xn))

#setup()
#print name = "__name__"

if __name__ == "__main__" and len(sys.argv) > 2:

  command = sys.argv[1]

  if command == "update":
    r_list = []
    for s in sys.argv[2:]:
      for a in s.split(','):
        l = a.split('..')
        if len(l) == 1:
          r_list.append(int(l[0]))
        else:
          for i in range(int(l[0]), int(l[-1])+1):
            r_list.append(i)
        print "Roaches to be set up", r_list
    
    for n in r_list:
      roach2_host = 'roach2-%02d' % (int(n))
      setup(r_name = roach2_host)
      for zdok in [0, 1]:
        print "Running og_from_noise for %s:zdok=%d" % (roach2_host, zdok)
        set_zdok(zdok)
        set_ogp()
        og_from_noise()
        update_ogp()
  elif command == 'setup':
    setup(sys.argv[2])
  else:
    print "Please run setup('roach2-nn') before trying anything else"

