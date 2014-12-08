#!/usr/bin/env python
meta_file = open ("wham_meta_40.txt", "w")

for i in range(40):
	meta_file.write ("logs/log_" + str(i) + ".hist\t" + str(1+0.05*i) + "\t1\n")

meta_file.close()
