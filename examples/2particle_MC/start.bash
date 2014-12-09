#!/bin/bash

 mpirun -n 40 ~/driver/driver_white --windows 40 --script 2part_MC.lu --count 100000 -l
