import adc5g

roach = adc5g.LocalRoachClient()
print "Estimated FPGA clock rate: %.2f MHz" % round(roach.est_brd_clk(), 2)

print "Devices found:"
print '\n'.join(roach.listdev())
