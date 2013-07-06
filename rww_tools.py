# run ipython rww_tools -pylab -i
import sys
import os
import time
from corr import katcp_wrapper
#import katcp_wrapper
roach2=katcp_wrapper.FpgaClient('roach2-00.cfa.harvard.edu')
zdok=0
import adc5g
#import matplotlib.pyplot as plt
import numpy
import fit_cores
lanio = "lanio 131.142.9.146 "

freq = 10.070801
pwr = 1.0
numpoints=16384
samp_freq = 5000.0
snap_name = "scope_raw_0_snap"

def dosnap(fr=0, name="t", rpt = 1, donot_clear=False):
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

    rpt  The number of repeats.  Defaults to 1.  The c1 .. c4 files mentioned
         above are overwritten with each repeat, but new rows of data are added
	 to the .fit file for each pass.
  """
  from numpy import savetxt
  global freq
  avg_pwr_sinad = 0
  if fr == 0:
    fr = freq
  for i in range(rpt):
    snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
    savetxt(name, snap,fmt='%d')
    ogp, pwr_sinad = fit_cores.fit_snap(fr, samp_freq, name,\
       clear_avgs = i == 0 and not donot_clear, prnt = i == rpt-1)
    avg_pwr_sinad += pwr_sinad
  return ogp, avg_pwr_sinad/rpt

def dosim(freq=10.070801, name="sim", rpt = 1, exact=True):
  """
  Do the same analysis as dosnap, but on simulated data.
  The arguments have the same meaning as dosnap except that freq is not
  coupled to the global variable.  A random phase is generated for each pass
  """

  from numpy import savetxt
  for i in range(rpt):
    snap=get_sim_data(freq, exact)
    savetxt(name, snap,fmt='%d')
    fit_cores.fit_snap(freq, samp_freq, name, i == 0)

def simpsd(freq=318.0, rpt = 1, exact=True):
  """
  Make a simulated snapshot and do the psd analysis on it.  The
  sine wave will have a randon start phase.
  """
  import matplotlib.pyplot as plt
  from numpy import log10
  from pylab import psd, detrend_mean
  for i in range(rpt):
    data = get_sim_data(freq, exact)
    power, freqs = psd(data, numpoints, Fs=samp_freq*1e6, \
        detrend=detrend_mean, scale_by_freq=True)
    plt.clf()
    if i == 0:
      sp = power
    else:
      sp += power
  sp /= rpt
  print "about to plot", len(freqs)
  plt.step(freqs, 10*log10(sp))
  fd = open("sim.psd", 'w')
  for i in range(len(sp)):
    print >>fd, "%7.2f %6.1f" % (freqs[i]/1e6, 10*log10(sp[i]))

def get_sim_data(freq, exact=True):
  """
  Make a simulated snapshot of data
  """
  from numpy import math
  if exact:
    offs = [0,0,0,0]
    gains = [1,1,1,1]
  else:
    offs = [.2, .3, -.2, -.1]
    gains = [1.001, .9984, .999, 1.002]
  del_phi = 2 * math.pi * freq / samp_freq
  data = numpy.empty((numpoints), dtype='int32')
  phase = 2*math.pi * numpy.random.uniform()
  for n in range(numpoints):
    core = n&3
    data[n] = (128 + math.floor(0.5 + 119.0 * math.sin(del_phi * n + phase) + \
        offs[core]))*gains[core]
  return data

def dotest(plotcore = 1):
  """
  Put the adc in test mode and get a sample of the test vector.  Plot core 1
  by default.
  """
  import matplotlib.pyplot as plt
  global snap_name
  adc5g.set_spi_control(roach2, zdok, test=1)
  cores = (corea, corec, coreb, cored) = adc5g.get_test_vector(roach2, [snap_name,])
  if plotcore == 2:
    plotcore = 3
  elif plotcore == 3:
    plotcore = 2
  plt.plot(cores[plotcore])
  adc5g.set_spi_control(roach2, zdok)

def dopsd(nfft = numpoints, rpt = 10):
  """
  Takes a snapshot, then computes, plots and writes out the Power Spectral
  Density functions.  The psd function is written into a file named "psd".
  This file will be overwritten with each call.  Arguments:

  nfft The number of points in the psd function.  Defaults to 16384.  Since
       a snapshot has 16384 points, this is the maximum which should be used
  rpt  The numper of mesurements to be averaged for the plot and output file. 
  """
  from numpy import savetxt, log10, column_stack
  import matplotlib.pyplot as plt
  for i in range(rpt):
    power, freqs = adc5g.get_psd(roach2, snap_name, samp_freq*1e6, 8, nfft)
    if i == 0:
      sp = power
    else:
      sp += power
  sp /= rpt
  plt.step(freqs, 10*log10(sp))
  data = column_stack((freqs/1e6, 10*log10(sp)))
  savetxt("psd", data, fmt=('%7.2f', '%6.1f'))

def dopsdcores(nfft = numpoints/4, rpt = 10):
  from numpy import savetxt, log10, column_stack
  from matplotlib.mlab import psd, detrend_mean
  for i in range(rpt):
    snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
    power, freqs = psd(snap, nfft*4, Fs=samp_freq*1e6, detrend=detrend_mean, \
        scale_by_freq=False)
    if i == 0:
      psd_all = power[:1+nfft/2]
    else:
      psd_all += power[:1+nfft/2]
    power, freqs = psd(snap[0:: 4], nfft, Fs=samp_freq*.25e6, detrend=detrend_mean,\
        scale_by_freq=True)
    if i == 0:
      psd1 = power
    else:
      psd1 += power
    power, freqs = psd(snap[1:: 4], nfft, Fs=samp_freq*.25e6, detrend=detrend_mean,\
        scale_by_freq=True)
    if i == 0:
      psd2 = power
    else:
      psd2 += power
    power, freqs = psd(snap[2:: 4], nfft, Fs=samp_freq*.25e6, detrend=detrend_mean,\
        scale_by_freq=True)
    if i == 0:
      psd3 = power
    else:
      psd3 += power
    power, freqs = psd(snap[3:: 4], nfft, Fs=samp_freq*.25e6, detrend=detrend_mean,\
        scale_by_freq=True)
    if i == 0:
      psd4 = power
    else:
      psd4 += power
  data = column_stack((freqs*1e-6, 10*log10(psd_all/rpt), 10*log10(psd1/rpt), \
      10*log10(psd2/rpt), 10*log10(psd3/rpt), 10*log10(psd4/rpt)))
  savetxt("psd_cores", data, fmt=('%7.2f'))

def hist_from_snapshots(rpt = 10):
  from numpy import savetxt, column_stack, bincount, zeros, arange, sum, array
#  hist_all = zeros(256,dtype=int)
  hist1 = zeros(256,dtype=int)
  hist2 = zeros(256,dtype=int)
  hist3 = zeros(256,dtype=int)
  hist4 = zeros(256,dtype=int)
  for i in range(rpt):
    snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
    snap = 128 + array(snap)
#    hist = bincount(snap, minlength=256)
#    hist_all += hist
    hist = bincount(snap[0:: 4], minlength=256)
    hist1 += hist
    hist = bincount(snap[1:: 4], minlength=256)
    hist2 += hist
    hist = bincount(snap[2:: 4], minlength=256)
    hist3 += hist
    hist = bincount(snap[3:: 4], minlength=256)
    hist4 += hist
  data=column_stack((arange(-128., 128, dtype=int), hist1, hist2,
      hist3, hist4))
  savetxt("hist_cores", data, fmt=("%d"))
#  print "all ",sum(hist_all[0:128]), sum(hist_all[128:256])
  print "c1  ",sum(hist1[0:128]), sum(hist1[129:256])
  print "c2  ",sum(hist2[0:128]), sum(hist2[129:256])
  print "c3  ",sum(hist3[0:128]), sum(hist3[129:256])
  print "c4  ",sum(hist4[0:128]), sum(hist4[129:256])

def get_hist(fname="hist_cores"):
  from numpy import empty, sum, savetxt
  data = empty(shape=(256,5), dtype=int)
  for c in range(4):
    data[:, c+1] = adc5g.get_histogram(roach2, zdok, "abcd"[c])
  data[:,0] = range(-128, 128)
  savetxt(fname, data, fmt=("%d"))

def multifreq(start=100, end=560, step=50, repeat=10, do_sfdr=False):
  """
  Calls dosnap for a range of frequenciesi in MHz.  The actual frequencies are
  picked to have an odd number of cycles in the 16384 point snapshot.
  """
  from numpy import savetxt, log10
  sfd = open('sinad', 'a')
  f = samp_freq / numpoints
  nstart = int(0.5+start/f)
  nend = int(0.5+end/f)
  nstep = int(0.5+step/f)
  for n in range(nstart, nend, nstep):
    freq = f*n
    set_freq(freq)
#    ogp, avg_pwr_sinad = dosnap(rpt=repeat, donot_clear = False)
    ogp, avg_pwr_sinad = dosnap(rpt=repeat, donot_clear = n!=nstart)
    sinad = 10.0*log10(avg_pwr_sinad)
    print >>sfd, "%8.3f %7.2f" % (freq, sinad)
    if do_sfdr:
      dopsd(rpt=3)
      fit_cores.dosfdr(freq)
  savetxt("ogp.meas", ogp[3:], fmt="%8.4f")
  fit_cores.fit_inl()

def freqResp(start=100, end=2400, delta=50, repeat=10,powerlevel=7):
  """
  Runs adc5g.get_psd for a range of frequenciesi in MHz.  
  The actual frequencies are picked to have an odd number of 
  cycles in the 16384 point snapshot. (This part is copied from multifreq()).
  Writes out freq and max() of spectrum to freqResponse.dat file.
  """
  from numpy import log10
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
    dosnap(rpt=repeat)

def update_ogp(fname = 'ogp', set=True):
  """
  Retreive the ogp data from the ADC and add in the corrections from
  the measured ogp (in ogp.meas).  Store in the file 'ogp'
  """
  from numpy import savetxt, genfromtxt
  cur_ogp = get_ogp_array()
  meas_ogp = genfromtxt("ogp.meas")
  # Correct for the ~1.4X larger effect of the phase registers than expected
  for i in (2,5,8,11):
    meas_ogp[i] *= 0.65
  savetxt(fname, cur_ogp+meas_ogp, fmt="%8.4f")
  if set:
    set_ogp()

def update_inl(fname = 'inl.meas', set=True):
  """
  Retreive the INL data from the ADC and add in the corrections from
  the measured inl (in inl.meas).  Store in the file 'inl'
  """
  from numpy import savetxt, genfromtxt
  cur_inl = get_inl_array()
  meas_inl = genfromtxt(fname)
  for level in range(17):
    cur_inl[level][1:] += meas_inl[level][1:]
  savetxt("inl", cur_inl, fmt=('%3d','%7.4f','%7.4f','%7.4f','%7.4f'))
  if set:
    set_inl()

def program(progname='adc5g_test.bof'):
  """
  Program the roach2 with the standard program.  After this, calibrate()
  should be called
  """

  roach2.progdev(progname)
  adc5g.set_spi_control(roach2, zdok)

def calibrate(verbose=False):
  """
  Call Rurik's routine to calibrate the time delay at the adc interface.
  """
  adc5g.set_test_mode(roach2, 0)
  adc5g.set_test_mode(roach2, 1)
  adc5g.sync_adc(roach2)
  opt0, glitches0 = adc5g.calibrate_mmcm_phase(roach2, 0, \
      ['scope_raw_0_snap',])
  if verbose or (opt0 == None):
    print "zodk0 ", opt0, glitches0
  else:
    print "zodk0", opt0
  opt1, glitches1 = adc5g.calibrate_mmcm_phase(roach2, 1, \
      ['scope_raw_1_snap',])
  if verbose or (opt1 == None):
    print "zodk1 ", opt1, glitches1
  else:
    print "zodk1", opt1
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

def set_ogp(fname = 'ogp'):
  """
  Clear the control register and then load the offset, gain and phase
  registers for each core.  These values are hard coded for now.
  """
  from numpy import genfromtxt
  adc5g.set_spi_control(roach2, zdok)
  t = genfromtxt(fname)
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
    

def set_inl(fname = 'inl'):
  """
  Set the INL registers for all four cores from a file containing 17 rows
  of 5 columns.  The first column contains the level and is ignored.
  Columns 2-5 contain the inl correction for cores a-d
  """
  from numpy import genfromtxt
  c = genfromtxt(fname, usecols=(1,2,3,4), unpack=True)
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
  from numpy import math
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
  from numpy import math
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
  from numpy import math
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
  from numpy import zeros
  inl = zeros((5,17), dtype='float')
  for chan in range(1,5):
    inl[chan] = adc5g.get_inl_registers(roach2, zdok, chan)
  inl[0] = range(0, 257,16)
  return inl.transpose()

def get_ogp_array():
  """
  Read  the Offset, Gain and Phase corrections for each core from the ADC
  and return in a 1D array
  """
  from numpy import zeros
  ogp = zeros((12), dtype='float')
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
  global zdok, snap_name
  snap_name = "scope_raw_%d_snap" % (zd)
  zdok = zd

def get_zdok():
  print "zdok %d, snapshot %s" % (zdok, snap_name)

def og_from_noise(fname="ogp.noise", rpt=100):
  """
  Take a number of snapshots of noise.  Analyze for offset and gain
  for each core separately.
  """
  from numpy import savetxt, zeros
  sum_result = zeros((15), dtype=float)
  sum_cnt = 0
  for n in range(rpt):
    result = zeros((15), dtype=float)
    snap=adc5g.get_snapshot(roach2, snap_name, man_trig=True, wait_period=2)
    l=float(len(snap))
    snap_off=sum(snap)/l
    snap_amp=sum(abs(snap-snap_off))/l
    result[0]=snap_off
    result[1]=snap_amp
    for core in range(4):
      # This will actually sample the cores in the order A,C,B,D
      # index will fix this up when data is put in the result array
      index=(3,9,6,12)[core]
      c=snap[core::4]
      l=float(len(c))
      off=sum(c)/l
      result[index] = (snap_off-off)*500.0/256.0
      amp=sum(abs(c-off))/l
      result[index+1]= 100.0*(snap_amp-amp)/snap_amp
    sum_result += result
    sum_cnt += 1
  sum_result /= sum_cnt
  print "%.4f "*15 % tuple(sum_result)
  savetxt(fname, sum_result[3:], fmt="%8.4f")

def phase_curve():
  f= 28/255.
  ofd = open('phasecurve', 'w')
  p = {}
  for i in range(1,5):
     p[i] = adc5g.get_spi_phase(roach2,zdok,i)
  for i in range(-10,11):
    set_phase(p[1]+ f*i,p[2] -f*i,p[3],p[4])
    ogp, gar = dosnap(rpt=5)
    print >>ofd, "%.3f %.3f %.3f" % (f*i, ogp[5], ogp[8])
  set_phase(p[1],p[2],p[3],p[4])
