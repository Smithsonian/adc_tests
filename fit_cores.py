#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
# Program to fit sine waves to each core of the adc in adc snapshots
#  and calculate the offsets, gains, etc. of  the individual cores
import sys
import os
import math
from scipy import *
from scipy.optimize import leastsq

def fitval(p, s, c):
  return p[0] +  p[1] * s + p[2] * c

#This is an array function if the arguments are of type array
def residuals(p, s, c, adc):
  return adc - fitval(p, s, c)

def fitc(sig_freq, samp_freq, df_name):
  p0 = [7.5, 5.0, 5.0]
  
  # Create lists to hold parameters at sample rate
  adc = []	# adc value
  s = []		# sin of signal freq
  c = []		# cos of signal freq
  
  del_phi = 2 * math.pi * sig_freq / samp_freq
  ifd=open(df_name, 'r')
  ffd = open(df_name + ".fit", 'a')
  outfile = df_name + ".c1"
  ffd1 = open(outfile, 'w')
#  print >>ffd, 'f1 = "%s"' % (outfile)
  outfile = df_name + ".c2"
  ffd2 = open(outfile, 'w')
#  print >>ffd, 'f2 = "%s"' % (outfile)
  outfile = df_name + ".c3"
  ffd3 = open(outfile, 'w')
#  print >>ffd, 'f3 = "%s"' % (outfile)
  outfile = df_name + ".c4"
  ffd4 = open(outfile, 'w')
#  print >>ffd, 'f4 = "%s"' % (outfile)
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
#  z_fact = 1.0
# Express delay in ps.  d_fact converts from angle at sig_freq(MHz) to ps
  d_fact = 1e12/(2*math.pi*sig_freq*1e6)
# This d_fact converts dly from angle at sig_freq to fraction of sample period.
#  d_fact = samp_freq/(2*math.pi*sig_freq)
#  d_fact = 1

  args1 = (array(s1), array(c1), array(core1))
  plsq1 = leastsq(residuals, p0, args1)
  #if plsq1[1] != 1:
  #  print "Fit failed to converge"
  z1 = z_fact * plsq1[0][0]
  s1a = plsq1[0][1]
  c1a = plsq1[0][2]
  amp1 = sqrt(s1a**2 + c1a**2)
  dly1 = d_fact*math.atan2(s1a, c1a)
  
  args2 = (array(s2), array(c2), array(core2))
  plsq2 = leastsq(residuals, p0, args2)
  z2 = z_fact * plsq2[0][0]
  s2a = plsq2[0][1]
  c2a = plsq2[0][2]
  amp2 = sqrt(s2a**2 + c2a**2)
  dly2 = d_fact*math.atan2(s2a, c2a)
  
  args3 = (array(s3), array(c3), array(core3))
  plsq3 = leastsq(residuals, p0, args3)
  z3 = z_fact * plsq3[0][0]
  s3a = plsq3[0][1]
  c3a = plsq3[0][2]
  amp3 = sqrt(s3a**2 + c3a**2)
  dly3 = d_fact*math.atan2(s3a, c3a)
  
  args4 = (array(s4), array(c4), array(core4))
  plsq4 = leastsq(residuals, p0, args4)
  z4 = z_fact * plsq4[0][0]
  s4a = plsq4[0][1]
  c4a = plsq4[0][2]
  amp4 = sqrt(s4a**2 + c4a**2)
  dly4 = d_fact*math.atan2(s4a, c4a)
  
  avz = (z1+z2+z3+z4)/4.0
  avamp = (amp1+amp2+amp3+amp4)/4.0
  a1p = 100*(amp1 - avamp)/avamp
  a2p = 100*(amp2 - avamp)/avamp
  a3p = 100*(amp3 - avamp)/avamp
  a4p = 100*(amp4 - avamp)/avamp
  avdly = (dly1+dly2+dly3+dly4)/4.0
  print "#%6.2f  zero(mV) amp(%%)  dly(ps)" % (sig_freq)
  print "#avg    %7.4f %7.4f %8.4f" %  (avz, avamp, avdly)
  print "core A  %7.4f %7.4f %8.4f" %  (z1-avz, a1p, dly1-avdly)
  print "core B  %7.4f %7.4f %8.4f" %  (z3-avz, a3p, dly3-avdly)
  print "core C  %7.4f %7.4f %8.4f" %  (z2-avz, a2p, dly2-avdly)
  print "core D  %7.4f %7.4f %8.4f" %  (z4-avz, a4p, dly4-avdly)
  print

  print >> ffd, "%9.4f %9.4f %9.4f %7.4f %7.4f %8.4f %7.4f %7.4f %8.4f %7.4f %7.4f %8.4f %7.4f %7.4f %8.4f" % (sig_freq, avz, avamp, z1-avz, a1p, dly1-avdly, z3-avz, a3p, dly3-avdly, z2-avz, a2p, dly2-avdly, z4-avz, a4p, dly4-avdly)
  
  Fit1 = fitval(plsq1[0], args1[0], args1[1])
  for i in range(data_cnt/4):
    print >>ffd1, "%d %d %.2f" % (4 * i, core1[i], Fit1[i])
  Fit2 = fitval(plsq2[0], args2[0], args2[1])
  for i in range(data_cnt/4):
    print >>ffd3, "%d %d %.2f" % (4 * i + 1, core2[i], Fit2[i])
  Fit3 = fitval(plsq3[0], args3[0], args3[1])
  for i in range(data_cnt/4):
    print >>ffd2, "%d %d %.2f" % (4 * i + 2, core3[i], Fit3[i])
  Fit4 = fitval(plsq4[0], args4[0], args4[1])
  for i in range(data_cnt/4):
    print >>ffd4, "%d %d %.2f" % (4 * i + 3, core4[i], Fit4[i])
