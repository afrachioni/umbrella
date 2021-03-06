#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import pylab
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import re

import sys

arg = sys.argv[1]

f = open(arg,"r")
#f = open("logs/log_7.txt","r")

mylist = []
zeroes = []
minus = []

column = 1

line = "#"
while (re.match ("^#", line)):
	line = f.readline()

old_sample = line.split()[column]
for line in f:
	if (re.match ("^#", line)): # the carat may be unnecessary
		continue
	try:
		vals = line.split()
		lineval = vals[column]
		index = int(vals[0])
	except:
		print "Parse error! Line: " + str(i)
		print line
		exit()
	if (index % 1 == 0):
		sample = lineval
		zeroes.append(sample)
		minus.append(old_sample)
		old_sample = sample


plt.plot(zeroes, minus, 'b.')
plt.xlabel('value of sample')
plt.ylabel('value of previous sample')
#plt.plot(mylist[:,0], mylist[:,1])
pylab.savefig(arg + ".png")
