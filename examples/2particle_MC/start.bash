#!/bin/bash

mpirun -n 40 ~/umbrella/driver_* --windows 40 --script 2part_MC.lu --count 100000 -l
