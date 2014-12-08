#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import pylab
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy

def rescale(arr):
	min = 1e20
	for x in arr:
		if x < min:
			min = x
	for i in range(len(arr)):
		arr[i] = (arr[i] - min) # / something


def flatten(arr):
	for elem in arr:
		elem[1] -= 0
		#elem[1] -= lj(elem[0])
#############################

kt = 8.617e-4
#f = lambda x : 1 * kt * numpy.log(4 * numpy.pi * x**2)
f = lambda x : 0

sigma = 1
epsilon = 10 * kt

trunc = [0, 1]

lj = lambda x : 4 * epsilon * ((x / sigma)**(-12) - (x / sigma)**(-6)) + epsilon
x = numpy.arange(1, 3, 0.01)
y = map(lj, x)
plt.plot(x, y, color='black')

#################
## Trapezoids
#data = numpy.loadtxt("A_trap.txt")[(2*trunc[0]):(-2*trunc[1])]
#flatten(data)
#corrected = data[:,1] + map(f, data[:,0])
#rescale(corrected)
## line up tails
##corrected = map(lambda x : x + lj(data[-1,0]) - corrected[-1], corrected)
#plt.plot(data[:,0], corrected, 'b+')
#
#
#################
## Simpson's rule
#data = numpy.loadtxt("A_simps.txt")[(trunc[0]):(-trunc[1])]
#flatten(data)
#corrected = data[:,1] + map(f, data[:,0])
#rescale(corrected)
## line up tails
##corrected = map(lambda x : x + lj(data[-1,0]) - corrected[-1], corrected)
#plt.plot(data[:,0], corrected, 'rx')



data = numpy.loadtxt("wham_free.txt")[8:-8]
flatten(data)
corrected = data[:,1] + map(f, data[:,0])
rescale(corrected)
plt.plot(data[:,0], corrected, 'g-')

#x = numpy.arange(-1.5, 1.5, 0.05)
#y = x**2
#plt.plot(x, y)

plt.xlabel('x')
plt.ylabel('A')
plt.grid(True)
plotfname = "pmf.png"
pylab.savefig(plotfname)
