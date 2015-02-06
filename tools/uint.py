#!/usr/bin/env python
import re
import numpy
import scipy.integrate # Simpson's rule implementation for integrating dA/dQ
import os
import pickle

kb = 1
log_dir = 'ag_logs_interval_target/'

kt = 0.741
all_springs = 0.03*kt
min_Q = 0
max_Q = 390
num = max_Q

# WindowStats class to hold window data and do some math
class WindowStats:
	def __init__(self, N, center, spring, mean, sd):
		self.N = N
		self.center = center
		self.spring = spring
		self.mean = mean
		self.sd = sd
	def P_i(self, Q): # Numerator of (8)
		return numpy.exp(-0.5*((Q-self.mean)/self.sd)**2)/ \
			(self.sd*numpy.sqrt(2*numpy.pi))
	def dA_idQ(self, Q): # (6)
		return kt * (Q-self.mean)/self.sd**2 - kb*self.spring*(Q - self.center)

def comb (min_val, max_val, num, index):
	return min_val + index * (max_val - min_val) / float(num - 1)

# Do things
log_files = os.listdir(log_dir)
window_data = []


window_data_fname = "window_data.pickle"
dQ = (max_Q - min_Q) / float(num)

if os.path.isfile(window_data_fname):
	print "Found " + window_data_fname + "..."
	window_data_file = open (window_data_fname, "r")
	window_data = pickle.load(window_data_file)
else:
	print "No window data file, parsing log files..."
	for log_fname in log_files:
		log_file = open(log_dir + log_fname, 'r')
		list = []
		for line in log_file:
			if (re.match (".*Target", line)):
				window_center = float(line.split()[-1])
			if (re.match ("^#", line)): # the carat may be unnecessary
				continue

			tokens = line.split()
			try:
				list.append(float(tokens[1]))
			except:
				print "Error with line: " + line

		print (log_fname + "\t" + str(len(list)) + "\t" + \
			str(window_center) + "\t" + str(numpy.mean(list)) + "\t" + \
			str(numpy.std(list)))
		
		window_data.append(WindowStats( \
			len(list), window_center, all_springs, \
			numpy.mean(list), numpy.std(list)))

	window_data_file = open (window_data_fname, "w")
	pickle.dump(window_data, window_data_file, 2)
window_data_file.close()


dAdQ = []
bin_centers = []
for i in range(num):
	bin_center = comb(min_Q, max_Q, num, i)
	bin_centers.append(bin_center)

	weight_denom = 0 # Denominator of (8)
	for w in window_data:
		w.weight = w.N * w.P_i(bin_center)
		weight_denom += w.weight
	if (weight_denom == 0):
		dAdQ.append(0)
		continue
	
	#print(weight_denom)
	dAdQ_here = 0 # (6)
	for w in window_data:
		dAdQ_here += w.weight/weight_denom * w.dA_idQ(bin_center)
	dAdQ.append(dAdQ_here)

print("integrating...")
#A = scipy.integrate.simps(dAdQ, x=bin_centers, even='avg')
A = scipy.integrate.cumtrapz(dAdQ, x=bin_centers)

print("Writing file...")
out_file = open("pmf.txt", "w")
out_file.write("#Q\t\t\t\tdA/dQ\n")
for i in range(len(dAdQ)):
	out_file.write(str(bin_centers[i]) + "\t" + str(dAdQ[i]) + "\n")
out_file.close()

# Set A = 0 at global minimum
min = 1e20 # big
for x in A:
	if x < min:
		min = x
for i in range(len(A)):
	A[i] = A[i] - min


out_file = open ("A.txt", "w")
out_file.write("#Q\t\t\t\tA\n")
for i in range(len(A)):
	out_file.write(str(bin_centers[i] + 1*dQ/2) + "\t" + str(A[i]) + "\n")
out_file.close()
