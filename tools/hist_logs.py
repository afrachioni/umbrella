#!/usr/bin/env python
#
# I am a little script for drawing a histogram of some quantity
#  recorded in each per-window log.
#

import numpy as np
import os
import re
import matplotlib
matplotlib.use('Agg')
import pylab
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


perfile = False
column = 2
dirmode = 1
dirname = "logs/" # must end in path separator
single_filename = "qdotq.txt"


# I'm tired of not being able to easily sort a bunch of plots by timestep
if dirmode:
	filenames = []
	for filename in os.listdir(dirname):
		filenames.append(filename)
	filenames.sort(key=lambda x:int(re.search("\d+", x).group(0)))
else:
	dirname = ""
	filenames = [single_filename]

for filename in filenames:
	if re.search(".lammps", filename):
		continue
	infile = open (dirname + filename, "r")
	line = ""
	while not re.search("#Step", line):
		line = infile.readline()

	header = line.split()[ column ]
	print "Working on histogram of " + header + " for " + filename

	data = []
	for line in infile:
		data.append(float( line.split()[column] ))
	infile.close()
	del infile

	weights = np.ones_like(data) / float(len(data))

	n, bins, patches = plt.hist(data, bins=map(lambda x:x-0.5, np.arange(0, 15, 0.1)), facecolor='green', alpha=0.75, align='mid')
	del data[:]

	plt.xlabel('Connections / atom')
	plt.ylabel('Frequency')
	plt.title(filename)
	plt.axis([0, 15, 0, 10000])
	plt.grid(True)

	if perfile:
		if dirmode:
			if not os.path.isdir("plots"):
				os.mkdir("plots")
			plotfname = "plots/" + filename + "_hist.png"
		else:
			plotfname = filename + "_hist.png"
		pylab.savefig(plotfname)
		plt.clf()

if not perfile:
	pylab.savefig("ag_hist.png")
