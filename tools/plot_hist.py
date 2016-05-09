#!/usr/bin/env python
import numpy
import re, os
import matplotlib.pyplot as plt
import pylab

a = 0.8
plt.rc('axes', color_cycle=[[0,0,a], [0,a,0], [a,0,0]])

files = os.listdir('logs')
lookup = sorted([[int(re.search('\d+', elem).group(0)), elem]
    for elem in files], key=lambda x:x[0])

for n, f in lookup:
    if not f.endswith('.hist'): continue;
    if n % 10 > 0: continue;
    x = numpy.loadtxt('logs/' + f, ndmin=2)
    if not x.size: continue;
    #plt.plot(x[:,2], numpy.log10(x[:,1]))
    plt.plot(x[:,2], x[:,1])

pylab.savefig('plot.png')
