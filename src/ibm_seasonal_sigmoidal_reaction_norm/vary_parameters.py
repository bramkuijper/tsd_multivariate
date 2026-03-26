#!/usr/bin/env python3

import numpy as np
import socket
import datetime
import sys
import copy

# name of the executable
exe_name = "./tsd_seasonal.exe"

# get the timestamp to create output file names
# which will be passed as a command-line argument
date = datetime.datetime.now()
base_name = "sim_seasonal_tsd_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

# auxiliary variable counting the number of simulations
# that need to be run 
ctr = 0

#  maximum duration the simulation is running
#max_simulation_time = int(20000000)
max_simulation_time = int(1_000_000)

# annual survival probability of female and male parents
sm = [ 0 ]
sf = [ 0 ]

# optimal temperature of females and males
t_opt = [[0.5,0.8],[0.8,0.5]]

# width of the selection function
omega = 0.01

# number of simulation replicates
nrep = 5

# how much temperature should change in a particular generation
temp_change = [0.6] #np.linspace(start = 0, stop=1.0,num=60)

# mutation rate of the threshold at which an individual should breed
mu_t = [0.02]

mu_tb = [0.0,0.02]

init_threshold = None

# mutation rate of the threshold 
mu_threshold = None
temp_error =[0.0] #[0,0.01,0.05,0.5,1.0]
cue_error = [0.0] #[0,0.01,0.05,0.5,1.0]

# aim for roughly 25000 lines
skip = max([int(round(max_simulation_time / 25000, 0)),1])

simulation_time_change = max_simulation_time - 1000

for mu_ti in mu_t:
    for mu_tbi in mu_tb:

        if mu_tbi > 0:
            init_threshold = 0
            mu_threshold_i = 0
        else:
            init_threshold = 5
            mu_threshold_i = 0.1

        for change_i in temp_change:
            for rep_i in range(0,nrep):
                for sm_i in sm:
                    for sf_i in sf:
                        if sf_i != sm_i:
                            continue

                        for t_opt_i in t_opt:
                            t_opt_fi = t_opt_i[0]
                            t_opt_mi = t_opt_i[1]

                            for temp_error_i in temp_error:
                                for cue_error_i in cue_error:

                                    # get the basename
                                    ctr += 1
                                    base_name_i = base_name + "_" + str(ctr)

                                    print(f"{exe_name} {base_name_i} " \
                                            f"{max_simulation_time} " \
                                            f"{simulation_time_change} " \
                                            f"{sm_i} " \
                                            f"{sf_i} " \
                                            f"{omega} " \
                                            f"{omega} " \
                                            f"{t_opt_mi} " \
                                            f"{t_opt_fi} " \
                                            f"{skip} " \
                                            f"{change_i} " \
                                            f"{mu_ti} " \
                                            f"{mu_tbi} " \
                                            f"{mu_threshold_i} " \
                                            f"{temp_error_i} " \
                                            f"{cue_error_i} " \
                                            f"{init_threshold} " \
                                          )
