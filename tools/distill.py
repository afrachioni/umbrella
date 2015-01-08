#!/usr/bin/env python
import os, subprocess
stride = 1000

if not os.path.exists("logs_distilled"):
	os.makedirs("logs_distilled")

logs = os.listdir("logs")
for fname in logs:
	in_f = "logs/" + fname
	out_f = "logs_distilled/" + fname
	print in_f + " -> " + out_f
	subprocess.call("~/driver/tools/log_distill " 
			+ in_f + " " + out_f + " " + str(stride), shell = True)
