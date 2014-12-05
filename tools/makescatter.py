#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import re
import numpy

i = 0
fname = "logs/log_" + str(i) + ".txt"
data_file = open (fname, "r")
for line in data_file:
	if re.match ("#Step", line):
		break
data_file.close()
headers = line.split()
data = numpy.loadtxt(fname)

for col in range(len(data[0])):
	plt.plot(data[1000:,col])

	plt.grid(True)
	pylab.savefig("plots/" + headers[col] + "_" + str(i) + ".png")
	plt.clf()
