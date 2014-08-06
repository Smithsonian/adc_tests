#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
# Program to fit sine waves to each core of the adc in adc snapshots
#  and calculate the offsets, gains, etc. of  the individual cores
import sys
import os
import math
#from scipy import *
#import adc5g
from numpy import array, zeros, savetxt, genfromtxt, shape, size
from numpy import sum, cumsum, genfromtxt, max, min
from scipy.optimize import leastsq
from numpy import arccos, pi, empty, arange, array, absolute
from matplotlib.pyplot import plot
from scipy.special import erfc

sum_result = zeros((15), dtype=float)
result_cnt = 0
code_errors = zeros((256,4), dtype='float')
ce_counts = zeros((256, 4), dtype='int32')

def fitsin(p, s, c):
  return p[0] +  p[1] * s + p[2] * c

def sin_residuals(p, s, c, adc):
  res = adc - fitsin(p, s, c)
  for i in range(adc.size):
    if adc[i] == -128 or adc[i] == 127:
      res[i] = 0
  return res

def fit_snap(adc, sig_freq, samp_freq, fname, clear_avgs=True, prnt=True):
  """
  Given a file containing a snapshot of data, separate the data from the
  4 cores and fit a separate sine wave to each.  From the dc offset, gain
  and phase of the four fits, report the average and the difference of each
  core from the average.  Write a line in fname.fit giving these values.
  Compute the average difference between the fitted value and measured value
  for each level of each core averaged over the samples (the raw data for
  INL corrections) and write to fname.res.

  Internally, cores 1-4 are in time sequence, but when the data is
  written out, write in the sequence 1324 for cores abcd.
  """
  global sum_result, result_cnt, code_errors, ce_counts
  p0 = [128.0, 90.0, 90.0]
  ogp = ()
  
  # Create lists to hold parameters at sample rate
#  adc = []	# adc value
  s = []		# sin of signal freq
  c = []		# cos of signal freq
  
  del_phi = 2 * math.pi * sig_freq / samp_freq
#  ifd=open(fname, 'r')
  if prnt:
    outfile = fname + ".a"
    cfd1 = open(outfile, 'w')
    outfile = fname + ".b"
    cfd2 = open(outfile, 'w')
    outfile = fname + ".c"
    cfd3 = open(outfile, 'w')
    outfile = fname + ".d"
    cfd4 = open(outfile, 'w')
#  data_cnt = 0
#  for line in ifd:
#    if line[0] == "#":
#      continue
  for data_cnt in range(len(adc)):
#    adc += [int(line)]
    s += [math.sin(del_phi * data_cnt)]
    c += [math.cos(del_phi * data_cnt)]
#    data_cnt += 1
  
  core1 = adc[0:: 4]
  core2 = adc[1:: 4]
  core3 = adc[2:: 4]
  core4 = adc[3:: 4]
  s1 = s[0:: 4]
  s2 = s[1:: 4]
  s3 = s[2:: 4]
  s4 = s[3:: 4]
  c1 = c[0:: 4]
  c2 = c[1:: 4]
  c3 = c[2:: 4]
  c4 = c[3:: 4]
  
#  ifd.close()


# express offsets as mV.  1 lsb = 500mV/256. z_fact converts from lsb to mV
# negate z_fact for negative feedback
  z_fact = -500.0/256.0
  true_zero = 0.0 * z_fact
#  z_fact = 1.0
# Express delay in ps.  d_fact converts from angle at sig_freq(MHz) to ps
  d_fact = 1e12/(2*math.pi*sig_freq*1e6)
# This d_fact converts dly from angle at sig_freq to fraction of sample period.
#  d_fact = samp_freq/(2*math.pi*sig_freq)
#  d_fact = 1

  args0 = (array(s), array(c), array(adc))
  plsq0 = leastsq(sin_residuals, p0, args0)
  #if plsq1[1] != 1:
  #  print "Fit failed to converge"
  z0 = z_fact * plsq0[0][0]
  s0a = plsq0[0][1]
  c0a = plsq0[0][2]
  amp0 = math.sqrt(s0a**2 + c0a**2)
  dly0 = d_fact*math.atan2(s0a, c0a)
  Fit0 = fitsin(plsq0[0], args0[0], args0[1])
  savetxt(fname+".fit", Fit0)
  ssq0 = 0.0
  for i in range(data_cnt):
    ssq0 += (adc[i] - Fit0[i])**2
  pwr_sinad = (amp0**2)/(2*ssq0/data_cnt)

  args1 = (array(s1), array(c1), array(core1))
  plsq1 = leastsq(sin_residuals, p0, args1)
  #if plsq1[1] != 1:
  #  print "Fit failed to converge"
  z1 = z_fact * plsq1[0][0]
  s1a = plsq1[0][1]
  c1a = plsq1[0][2]
  amp1 = math.sqrt(s1a**2 + c1a**2)
  dly1 = d_fact*math.atan2(s1a, c1a)
  
  args2 = (array(s2), array(c2), array(core2))
  plsq2 = leastsq(sin_residuals, p0, args2)
  z2 = z_fact * plsq2[0][0]
  s2a = plsq2[0][1]
  c2a = plsq2[0][2]
  amp2 = math.sqrt(s2a**2 + c2a**2)
  dly2 = d_fact*math.atan2(s2a, c2a)
  
  args3 = (array(s3), array(c3), array(core3))
  plsq3 = leastsq(sin_residuals, p0, args3)
  z3 = z_fact * plsq3[0][0]
  s3a = plsq3[0][1]
  c3a = plsq3[0][2]
  amp3 = math.sqrt(s3a**2 + c3a**2)
  dly3 = d_fact*math.atan2(s3a, c3a)
  
  args4 = (array(s4), array(c4), array(core4))
  plsq4 = leastsq(sin_residuals, p0, args4)
  z4 = z_fact * plsq4[0][0]
  s4a = plsq4[0][1]
  c4a = plsq4[0][2]
  amp4 = math.sqrt(s4a**2 + c4a**2)
  dly4 = d_fact*math.atan2(s4a, c4a)

  avz = (z1+z2+z3+z4)/4.0
  avamp = (amp1+amp2+amp3+amp4)/4.0
  # Reverse the amplitude and zero differences so they can be applied to the
  # offset and gain registers directly.  The phase registers don't need the
  # reversal
  a1p = 100*(avamp -amp1)/avamp
  a2p = 100*(avamp -amp2)/avamp
  a3p = 100*(avamp -amp3)/avamp
  a4p = 100*(avamp -amp4)/avamp
  avdly = (dly1+dly2+dly3+dly4)/4.0
#  print s1a, c1a, s2a, c2a, s3a, c3a, s4a, c4a
  if prnt:
    print "#%6.2f  zero(mV) amp(%%)  dly(ps) (adj by .4, .14, .11)" % (sig_freq)
    print "#avg    %7.4f %7.4f %8.4f" %  (avz, avamp, avdly)
    print "core A  %7.4f %7.4f %8.4f" %  (z1 -true_zero, a1p, dly1-avdly)
    print "core B  %7.4f %7.4f %8.4f" %  (z3 -true_zero, a3p, dly3-avdly)
    print "core C  %7.4f %7.4f %8.4f" %  (z2 -true_zero, a2p, dly2-avdly)
    print "core D  %7.4f %7.4f %8.4f" %  (z4 -true_zero, a4p, dly4-avdly)
    print "\nsinad = %.2f" % (10.0*math.log10(pwr_sinad))

  if clear_avgs:
    sum_result = zeros((15), dtype=float)
    result_cnt = 0
    code_errors = zeros((256,4), dtype='float')
    ce_counts = zeros((256, 4), dtype='int32')

  result = (sig_freq, avz, avamp,\
      z1-true_zero, a1p, dly1-avdly, z3-true_zero, a3p, dly3-avdly, \
      z2-true_zero, a2p, dly2-avdly, z4-true_zero, a4p, dly4-avdly)
  result_fmt = "%8.4f "*15
  sum_result += array(result)
  result_cnt += 1
#  print "%.3f %.3f %.3f %.3f %d" % (sum_result[5]/result_cnt,\
#     sum_result[8]/result_cnt,\
#     sum_result[11]/result_cnt, sum_result[14]/result_cnt, result_cnt)
  if prnt and result_cnt > 1:
    ofd = open(fname + ".ogp", 'a')
    avg_result = sum_result/result_cnt
    avg_result[0] = sig_freq
    ogp = tuple(avg_result)
    print >>ofd, result_cnt,
    print >>ofd,  result_fmt % ogp
    ofd.close()
    print "average of %d measurements" % (result_cnt)
    print "#avg    %7.4f %7.4f %8.4f" %  (ogp[1], ogp[2], 0)
    print "core A  %7.4f %7.4f %8.4f" %  ogp[3:6]
    print "core B  %7.4f %7.4f %8.4f" %  ogp[6:9]
    print "core C  %7.4f %7.4f %8.4f" %  ogp[9:12]
    print "core D  %7.4f %7.4f %8.4f" %  ogp[12:15]
    print

  # for each core (n), accumulate the sum of the residuals at each output code
  # in code_errors[code][n]
  # and the count of residuals added in ce_counts[code][n]
  Fit1 = fitsin(plsq1[0], args1[0], args1[1])
  for i in range(data_cnt/4):
    code = core1[i]
    if prnt:
      print >>cfd1, "%d %d %.2f" % (4 * i, code, Fit1[i])
    code_errors[code+128][0] += code - Fit1[i]
    ce_counts[code+128][0] += 1
  Fit2 = fitsin(plsq2[0], args2[0], args2[1])
  for i in range(data_cnt/4):
    code = core2[i]
    if prnt:
      print >>cfd3, "%d %d %.2f" % (4 * i + 1, code, Fit2[i])
    code_errors[code+128][1] += code - Fit2[i]
    ce_counts[code+128][1] += 1
  Fit3 = fitsin(plsq3[0], args3[0], args3[1])
  for i in range(data_cnt/4):
    code = core3[i]
    if prnt:
      print >>cfd2, "%d %d %.2f" % (4 * i + 2, code, Fit3[i])
    code_errors[code+128][2] += code - Fit3[i]
    ce_counts[code+128][2] += 1
  Fit4 = fitsin(plsq4[0], args4[0], args4[1])
  for i in range(data_cnt/4):
    code = core4[i]
    if prnt:
      print >>cfd4, "%d %d %.2f" % (4 * i + 3, code, Fit4[i])
    code_errors[code+128][3] += code - Fit4[i]
    ce_counts[code+128][3] += 1
  if prnt:
    rfd = open(fname + '.res', "w")
    # Since the INL registers are addressed as offset binary, generate the
    # .res file that way
    for code in range(256):
      if ce_counts[code].min() > 1:
        e = code_errors[code]/ce_counts[code]
        print >>rfd, "%3d %5.3f %5.3f %5.3f %5.3f" % \
	    (code, e[0], e[2], e[1], e[3])
      else:
        print >>rfd, "%3d %5.3f %5.3f %5.3f %5.3f" % (code,0,0,0,0)
  if prnt:
    cfd1.close()
    cfd2.close()
    cfd3.close()
    cfd4.close()
    rfd.close()
  return ogp, pwr_sinad

def fit_inl(fname='t.res'):
  """
  Read the raw residuals from fname.res and compute the INL corrections
  Assume that the residuals file in in core order ie. a,b,c,d.
  """
  
  corrections = zeros((17,5), dtype='float')
  wts = array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15,16,\
        15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.])

  data = genfromtxt(fname, unpack=True)
  start_data = int(data[0][0])
  file_limit = len(data[0])
  data_limit = start_data + file_limit
  print start_data, data_limit
  if data[0][data_limit - start_data - 1] != data_limit - 1:
    raise RuntimeError("there are holes in the data file")
#    print "there are holes in the data file"
#    return
  for corr_level in range(17):
    # a and b are positions in the file
    a = corr_level*16 - 15 - start_data
    b = a + 31
    if a < 0:
      a = 0
    if a == 0 and start_data == 0:
      a = 1
    if b > file_limit:
      b = file_limit
    if b == file_limit and data_limit == 256:
      b -= 1
    if a > b:
      continue
    wt_a = a - corr_level*16 + 15 + start_data
    wt_b = wt_a -a + b
    wt = sum(wts[wt_a:wt_b])
    av1 = sum(data[1][a:b]*wts[wt_a:wt_b])/wt
    av2 = sum(data[2][a:b]*wts[wt_a:wt_b])/wt
    av3 = sum(data[3][a:b]*wts[wt_a:wt_b])/wt
    av4 = sum(data[4][a:b]*wts[wt_a:wt_b])/wt
    print "%d %7.5f %7.5f %7.5f %7.5f" %  (16*corr_level,av1,av2,av3,av4)
    corrections[corr_level][0] = 16*corr_level
    corrections[corr_level][1] = av1
    corrections[corr_level][2] = av2
    corrections[corr_level][3] = av3
    corrections[corr_level][4] = av4
  if fname[:4] == 'hist':
    outname = 'hist_inl.meas'
  else:
    outname = 'inl.meas'
  savetxt(outname, corrections, fmt=('%3d','%7.4f','%7.4f','%7.4f','%7.4f'))


def dosfdr(sig_freq, fname = 'psd'):
  """
  Read the psd data from a file and calculate the SFDR and SINAD.  Write the
  results in a file named sfdr
  """
  
  tot_pwr = 0.0
  in_peak = False
  spur_pwr = 0.0
  for line in open(fname, 'r'):
    f, d = line.split()
    freq = float(f)
    if abs(freq - sig_freq) < 4:
      test = -70
    else:
      test = -90
    db = float(d)
    pwr = 10**(float(db)/10.)
    tot_pwr += pwr
    if in_peak:
      if db < test:
        in_peak = False
        if abs(peak_freq - sig_freq) < 1:
          sig_pwr = pwr_in_peak
	  sig_db = peak_db
	  peak_sig_freq = peak_freq
        else:
  	  if pwr_in_peak > spur_pwr:
  	    spur_pwr = pwr_in_peak
            spur_db = peak_db
  	    spur_freq = peak_freq
      else:
        pwr_in_peak += 10**(float(db)/10.)
        if db > peak_db:
          peak_db = db
  	peak_freq = freq
    elif db > test:
      pwr_in_peak = 10**(float(db)/10.)
      peak_freq = freq
      peak_db = db
      in_peak = True
  outfd = open('sfdr', 'a')
  sfdr = 10.0*math.log10(sig_pwr / spur_pwr)
  sinad = 10.0*math.log10(sig_pwr/(tot_pwr - sig_pwr))
  print >> outfd, "%8.3f %6.2f %6.2f %6.2f %7.2f" %\
      (sig_freq, sig_db, sfdr, sinad, spur_freq)

# Start of code for histograms
def cumsin(p, codes):
  amp = p[0]
  arg=(codes-127.5+p[1])/amp
  for i in range(size(arg)):
    a = arg[i]
    if a > 1.0:
      arg[i] = 1.0
    elif a < -1.0:
      arg[i] = 0.0
    else:
      arg[i] = 1-arccos(a)/pi
  return arg
#  return array([(sign(a)+1)/2 if abs(a) > 1.0 else 1-arccos(a)/pi for a in arg])

def pltcumsin(p):
  codes = array(range(0,256), dtype=float)
  cum = cumsin(p, codes)
  plot(codes,cum)
  return cum

def cumgaussian(p, codes):
  amp = p[0]
  arg=(codes-127.5+p[1])/amp
  return (1-erfc(arg)/2.0)

def pltcumgaussian(p):
  codes = array(range(0,256), dtype=float)
  cum = cumgaussian(p, codes)
  plot(codes,cum)
  return cum

def hist_residuals(p, codes, cumhist, fit_function):
  return cumhist - fit_function(p, codes)

def fit_hist(core=1, type='sin', fname='hist_cores'):
  global cumhist, hist, plsq, cumresid, extended_fit

  coderesid=empty(256,dtype=float)
  if type == "sin":
    fit_function = cumsin
  else:
    fit_function = cumgaussian
  codes = arange(0,256, dtype=float)
  # get the data as a 6x256 array, ie. hist[0] == arange(256)
  hist = genfromtxt(fname, dtype=float, unpack=True)
  t=float(sum(hist[core]))
  cumhist=cumsum(hist[core])/t
  args = (codes[0:255], cumhist[0:255], fit_function)
  plsq = leastsq(hist_residuals, [135,0], args)
  extended_fit = empty(258, dtype=float)
  cumresid = hist_residuals(plsq[0],codes, cumhist, fit_function)
  extended_fit[1:257] = cumhist-cumresid
  extended_fit[0] = 2 * extended_fit[1] - extended_fit[2]
  extended_fit[257] = 2 * extended_fit[256] - extended_fit[255]
  # invert the sign of cumresid so inl corrections will be correct
  for i in range(256):
    if(cumresid[i] > 0 and extended_fit[i+2] > extended_fit[i+1]):
      coderesid[i] = -cumresid[i] / (extended_fit[i+2] - extended_fit[i+1])
    elif(extended_fit[i+1] > extended_fit[i]):
      coderesid[i] = -cumresid[i] / (extended_fit[i+1] - extended_fit[i])
    else:
      coderesid[i]=0
  return plsq[0], coderesid
