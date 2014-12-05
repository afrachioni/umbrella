#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import numpy

for i in [30]:
	print "Reading from out/out_" + str(i) + ".txt"
	#data = numpy.loadtxt("out/out_" + str(i) + ".txt")
	#plt.plot(data[10000:,1])
	data = numpy.loadtxt("out_short_" + str(i) + ".txt")
	plt.plot(data[:,1])

	plt.grid(True)
	pylab.savefig("plots/scatter_" + str(i) + ".png")
	plt.clf()
