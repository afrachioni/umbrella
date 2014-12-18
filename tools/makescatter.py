#!/usr/bin/env python

print "importing things..."
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import re
import numpy
import os

min = 0
max = -1

dirname = "log_sub/" # must end in path separator
dirname = "logs/" # must end in path separator
filenames = os.listdir(dirname)
filenames.sort(key=lambda x:int(re.search("\d+", x).group(0)))


for filename in filenames:
	if re.search(".lammps", filename) or re.search(".hist", filename):
		continue
	print "Working on log: " + filename
	#data_file = open (dirname + filename, "r")
	#for line in data_file:
	#	if re.match ("#Step", line):
	#		break
	#data_file.close()
	#headers = line.split()
	data = numpy.loadtxt(dirname + filename)

	for col in [1]: #range(len(data[0])):
		plt.plot(data[min:max,col])

		plt.grid(True)
		#pylab.savefig("plots/" + headers[col] + "_" + str(i) + ".png")
		#plt.clf()

plt.xlabel("Step index")
plt.ylabel("n")
pylab.savefig("ag_scatter.png")
