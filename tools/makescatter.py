#!/usr/bin/env python

print "importing things..."
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import re
import os

min = 0
max = -1

#dirname = "ag_logs_interval/" # must end in path separator
dirname = "logs/" # must end in path separator
filenames = os.listdir(dirname)
filenames.sort(key=lambda x:int(re.search("\d+", x).group(0)))

if not os.path.exists("plots"):
	os.makedirs("plots")


for filename in filenames:
	if re.search(".lammps", filename) or re.search(".hist", filename) or filename.endswith(".png") or filename.startswith("."):
		continue
	print "Working on log: " + filename
	data_file = open (dirname + filename, "r")
	for line in data_file:
		if re.match ("#Step", line):
			break
	headers = line.split()
	header_indices = range(len(headers))
	data = [[] for i in header_indices]

	for line in data_file:
		tokens = line.split()
		#if int(tokens[0]) % 1000:
			#continue
		for i in header_indices:
			data[i].append(float(tokens[i]))

	for col in range(len(headers)):
		plt.plot(data[0], data[col])

		plt.grid(True)
		plt.xlabel("Step index")
		plt.ylabel(headers[col])
		pylab.savefig("plots/" + headers[col] + "_" + filename + ".png")
		plt.clf()

#plt.xlabel("Step index")
#plt.ylabel("n")
#pylab.savefig("ag_scatter.png")
