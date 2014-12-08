#!/usr/bin/env python

print "importing things..."
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import re
import numpy

min = 0
max = -1

i = 0
for i in range(40):
	print "Working on log of index " + str(i)
	fname = "logs/log_" + str(i) + ".txt"
	data_file = open (fname, "r")
	for line in data_file:
		if re.match ("#Step", line):
			break
	data_file.close()
	headers = line.split()
	data = numpy.loadtxt(fname)

	for col in [1]: #range(len(data[0])):
		plt.plot(data[min:max,col])

		plt.grid(True)
		#pylab.savefig("plots/" + headers[col] + "_" + str(i) + ".png")
		#plt.clf()

pylab.savefig("ag_scatter.png")
