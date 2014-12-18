#!/usr/bin/env python
#
# I am a little script for drawing a histogram of some quantity
#  recorded in each per-window log.
#

print "Importing things..."

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

colors = [ 'red', 'blue', 'green' ]
color_index = -1
def getColor():
	global color_index, colors
	if color_index < len(colors) - 1:
		color_index = color_index + 1
		return colors[color_index]
	else:
		color_index = -1
		return getColor()

max_pop = -float("inf")
min_Q = float("inf")
max_Q = -float("inf")

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
# TODO read hist files
	if re.search(".lammps", filename) or re.search(".hist", filename):
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

	n, bins, patches = plt.hist(data, bins=map(lambda x:x-0.5, np.arange(0, 15, 0.01)), facecolor=getColor(), alpha=0.75, align='mid')
	if max(n) > max_pop:
		max_pop = max(n)

	for i in range(len(n)):
		if n[i] > 0 and bins[i] < min_Q:
			min_Q = bins[i]
		if n[i] > 0 and bins[i + 1] > max_Q:
			max_Q = bins[i + 1]

	del data[:]

	plt.xlabel(header)
	plt.ylabel('Frequency')
	#plt.title(filename)
	#plt.axis([1, 3, 0, 10000])
	plt.axis([min_Q, max_Q, 0, max_pop])
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