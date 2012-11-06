#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
# Program to fit sine waves to each core of the adc in adc snapshots
#  and calculate the offsets, gains, etc. of  the individual cores
import sys
import os
import math
#from scipy import *
import adc5g
from numpy import array, zeros, savetxt, genfromtxt, shape
from scipy.optimize import leastsq

sum_result = zeros((15), dtype=float)
result_cnt = 0
code_errors = zeros((256,4), dtype='float')
ce_counts = zeros((256, 4), dtype='int32')

def fitval(p, s, c):
  return p[0] +  p[1] * s + p[2] * c

#This is an array function if the arguments are of type array
def residuals(p, s, c, adc):
  return adc - fitval(p, s, c)

def fit_snap(sig_freq, samp_freq, df_name, clear_avgs=True, prnt=True):
  """
  Given a file containing a snapshot of data, separate the data from the
  4 cores and fit a separate sine wave to each.  From the dc offset, gain
  and phase of the four fits report the average and the difference of each
  core from the average.  Write a line in fname.fit giving these values.
  Compute the average difference between the fitted value and measured value
  for each level of each core averaged over the samples (the raw data for
  INL corrections) and write to fname.res.
  """
  global sum_result, result_cnt, code_errors, ce_counts
  p0 = [128.0, 90.0, 90.0]
  ogp = ()
  
  # Create lists to hold parameters at sample rate
  adc = []	# adc value
  s = []		# sin of signal freq
  c = []		# cos of signal freq
  
  del_phi = 2 * math.pi * sig_freq / samp_freq
  ifd=open(df_name, 'r')
  ofd = open(df_name + ".ogp", 'a')
  if prnt:
    outfile = df_name + ".c1"
    cfd1 = open(outfile, 'w')
    outfile = df_name + ".c2"
    cfd2 = open(outfile, 'w')
    outfile = df_name + ".c3"
    cfd3 = open(outfile, 'w')
    outfile = df_name + ".c4"
    cfd4 = open(outfile, 'w')
  data_cnt = 0
  for line in ifd:
    if line[0] == "#":
      continue
    adc += [int(line)]
    s += [math.sin(del_phi * data_cnt)]
    c += [math.cos(del_phi * data_cnt)]
    data_cnt += 1
  
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
  
  ifd.close()


# express offsets as mV.  1 lsb = 500mV/256. z_fact conf=verts from lsb to mV
  z_fact = 500.0/256.0
  true_zero = 127 * z_fact
#  z_fact = 1.0
# Express delay in ps.  d_fact converts from angle at sig_freq(MHz) to ps
  d_fact = 1e12/(2*math.pi*sig_freq*1e6)
# This d_fact converts dly from angle at sig_freq to fraction of sample period.
#  d_fact = samp_freq/(2*math.pi*sig_freq)
#  d_fact = 1

  args0 = (array(s), array(c), array(adc))
  plsq0 = leastsq(residuals, p0, args0)
  #if plsq1[1] != 1:
  #  print "Fit failed to converge"
  z0 = z_fact * plsq0[0][0]
  s0a = plsq0[0][1]
  c0a = plsq0[0][2]
  amp0 = math.sqrt(s0a**2 + c0a**2)
  dly0 = d_fact*math.atan2(s0a, c0a)
  Fit0 = fitval(plsq0[0], args0[0], args0[1])
  ssq0 = 0.0
  for i in range(data_cnt):
    ssq0 += (adc[i] - Fit0[i])**2
  pwr_sinad = (amp0**2)/(2*ssq0/data_cnt)

  args1 = (array(s1), array(c1), array(core1))
  plsq1 = leastsq(residuals, p0, args1)
  #if plsq1[1] != 1:
  #  print "Fit failed to converge"
  z1 = z_fact * plsq1[0][0]
  s1a = plsq1[0][1]
  c1a = plsq1[0][2]
  amp1 = math.sqrt(s1a**2 + c1a**2)
  dly1 = d_fact*math.atan2(s1a, c1a)
  
  args2 = (array(s2), array(c2), array(core2))
  plsq2 = leastsq(residuals, p0, args2)
  z2 = z_fact * plsq2[0][0]
  s2a = plsq2[0][1]
  c2a = plsq2[0][2]
  amp2 = math.sqrt(s2a**2 + c2a**2)
  dly2 = d_fact*math.atan2(s2a, c2a)
  
  args3 = (array(s3), array(c3), array(core3))
  plsq3 = leastsq(residuals, p0, args3)
  z3 = z_fact * plsq3[0][0]
  s3a = plsq3[0][1]
  c3a = plsq3[0][2]
  amp3 = math.sqrt(s3a**2 + c3a**2)
  dly3 = d_fact*math.atan2(s3a, c3a)
  
  args4 = (array(s4), array(c4), array(core4))
  plsq4 = leastsq(residuals, p0, args4)
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
  if prnt:
    print "#%6.2f  zero(mV) amp(%%)  dly(ps) (adj by .4, .14, .11)" % (sig_freq)
    print "#avg    %7.4f %7.4f %8.4f" %  (avz, avamp, avdly)
    print "core A  %7.4f %7.4f %8.4f" %  (true_zero-z1, a1p, dly1-avdly)
    print "core B  %7.4f %7.4f %8.4f" %  (true_zero-z3, a3p, dly3-avdly)
    print "core C  %7.4f %7.4f %8.4f" %  (true_zero-z2, a2p, dly2-avdly)
    print "core D  %7.4f %7.4f %8.4f" %  (true_zero-z4, a4p, dly4-avdly)
    print "\nsinad = %.2f" % (10.0*math.log10(pwr_sinad))

  if clear_avgs:
    sum_result = zeros((15), dtype=float)
    result_cnt = 0
    code_errors = zeros((256,4), dtype='float')
    ce_counts = zeros((256, 4), dtype='int32')

  result = (sig_freq, avz, avamp,\
      true_zero-z1, a1p, dly1-avdly, true_zero-z3, a3p, dly3-avdly, \
      true_zero-z2, a2p, dly2-avdly, true_zero-z4, a4p, dly4-avdly)
  result_fmt = "%8.4f "*15
  sum_result += array(result)
  result_cnt += 1
#  print "%.3f %.3f %.3f %.3f %d" % (sum_result[5]/result_cnt,\
#     sum_result[8]/result_cnt,\
#     sum_result[11]/result_cnt, sum_result[14]/result_cnt, result_cnt)
  if prnt and result_cnt > 1:
    avg_result = sum_result/result_cnt
    avg_result[0] = sig_freq
    ogp = tuple(avg_result)
    print >>ofd, result_cnt,
    print >>ofd,  result_fmt % ogp
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
  Fit1 = fitval(plsq1[0], args1[0], args1[1])
  for i in range(data_cnt/4):
    code = core1[i]
    if prnt:
      print >>cfd1, "%d %d %.2f" % (4 * i, code, Fit1[i])
    code_errors[code][0] += code - Fit1[i]
    ce_counts[code][0] += 1
  Fit2 = fitval(plsq2[0], args2[0], args2[1])
  for i in range(data_cnt/4):
    code = core2[i]
    if prnt:
      print >>cfd3, "%d %d %.2f" % (4 * i + 1, code, Fit2[i])
    code_errors[code][1] += code - Fit2[i]
    ce_counts[code][1] += 1
  Fit3 = fitval(plsq3[0], args3[0], args3[1])
  for i in range(data_cnt/4):
    code = core3[i]
    if prnt:
      print >>cfd2, "%d %d %.2f" % (4 * i + 2, code, Fit3[i])
    code_errors[code][2] += code - Fit3[i]
    ce_counts[code][2] += 1
  Fit4 = fitval(plsq4[0], args4[0], args4[1])
  for i in range(data_cnt/4):
    code = core4[i]
    if prnt:
      print >>cfd4, "%d %d %.2f" % (4 * i + 3, code, Fit4[i])
    code_errors[code][3] += code - Fit4[i]
    ce_counts[code][3] += 1
  if prnt:
    rfd = open(df_name + '.res', "w")
    for code in range(256):
      if ce_counts[code].min() > 4:
        e = code_errors[code]/ce_counts[code]
        print >>rfd, "%3d %5.3f %5.3f %5.3f %5.3f" % \
	    (code, e[0], e[1], e[2], e[3])
  return ogp, pwr_sinad

def fit_inl(df_name='t.res'):
  """
  Read the raw residuals from fname.res and compute the INL corrections
  """
  
  corrections = zeros((17,5), dtype='float')
  wts = array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15,16,\
        15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.])

  data = genfromtxt(df_name, unpack=True)
  start_data = int(data[0][0])
  file_limit = len(data[0])
  data_limit = start_data + file_limit
  print start_data, data_limit
  if data[0][data_limit - start_data - 1] != data_limit - 1:
    print "there are holes in the data file"
    return
  for corr_level in range(17):
    # a and b are positions in the file
    a = corr_level*16 - 15 - start_data
    b = a + 31
    if a < 0:
      a = 0
    if b > file_limit:
      b = file_limit
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
    corrections[corr_level][2] = av3
    corrections[corr_level][3] = av2
    corrections[corr_level][4] = av4
  savetxt("inl.meas", corrections, fmt=('%3d','%7.4f','%7.4f','%7.4f','%7.4f'))


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

#def dosfdr(sig_freq, fname = 'psd'):
#  sig_peak = -100.0
#  spur_peak = -100.0
#  spur_freq = 0.0
#  outfd = open('sfdr', 'a')
#  for line in open(fname, 'r'):
#    freq, db = line.split()
#    db = float(db)
#    freq = float(freq)
#    if abs(freq - sig_freq) < 1.5:
#      if db > sig_peak:
#        sig_peak = db
#    elif db > spur_peak:
#        spur_peak = db
#	spur_freq = freq
#  print >> outfd, "%8.3f %6.2f %6.2f %7.2f" %\
#      (sig_freq, sig_peak, sig_peak - spur_peak, spur_freq)
#
#def get_inl_array(roach, zdok_n):
#  """
#  Read the INL corrections from the adc and put in an array
#  """
#  inl = zeros((5,17), dtype='float')
#  for chan in range(1,5):
#    inl[chan] = adc5g.get_inl_registers(roach, zdok_n, chan)
#  inl[0] = range(0, 257,16)
#  return inl.transpose()
#
#def get_ogp_array(roach, zdok_n):
#  """
#  Read  the Offset, Gain and Phase corrections for each core from the ADC
#  and return in a 1D array
#  """
#  ogp = zeros((12), dtype='float')
#  indx = 0
#  for chan in range(1,5):
#    ogp[indx] = adc5g.get_spi_offset(roach,zdok_n,chan)
#    indx += 1
#    ogp[indx] = adc5g.get_spi_gain(roach,zdok_n,chan)
#    indx += 1
#    ogp[indx] = adc5g.get_spi_phase(roach,zdok_n,chan)
#    indx += 1
#  return ogp
