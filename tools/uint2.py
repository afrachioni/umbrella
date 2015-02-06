#!/usr/bin/env python
import numpy
import scipy.integrate


dadq = numpy.loadtxt("pmf.txt")
bin_centers = dadq[:,0]
derivative = dadq[:,1]

d_Q = bin_centers[1] - bin_centers[0]

A = scipy.integrate.cumtrapz(derivative, x=bin_centers)

out = []
out_centers = []
for i in range( (len(bin_centers)-1)/2 ):
	left = 2 * i
	center = 2 * i + 1
	right = 2 * i + 2

	width = bin_centers[right] - bin_centers[left]
	out.append( width/6 * (derivative[left] + 4*derivative[center] + derivative[right]) )
	out_centers.append(bin_centers[center])

aa = numpy.cumsum(out)
min = 1e20
for x in aa:
	if x < min:
		min = x
for i in range(len(aa)):
	aa[i] = aa[i] - min

simps_file = open("A_simps.txt", "w")
for i in range(len(aa)):
	simps_file.write(str(out_centers[i] + d_Q) + "\t" + str(aa[i]) + "\n")
simps_file.close()


min = 1e20
for x in A:
	if x < min:
		min = x
for i in range(len(A)):
	A[i] = A[i] - min

trap_file = open("A_trap.txt", "w")
for i in range(len(bin_centers) - 1):
	trap_file.write( str(bin_centers[i] + 1*d_Q) + "\t" + str(A[i]) + "\n")
trap_file.close()
