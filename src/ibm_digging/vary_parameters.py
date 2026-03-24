#!/usr/bin/env python3

import numpy as np
import socket
import datetime
import sys
import copy

exe_name = "./tsd_seasonal.exe"
date = datetime.datetime.now()
base_name = "sim_seasonal_tsd_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

ctr = 0
max_simulation_time = int(1e07)
max_simulation_time = int(1000000)
simulation_time_change = max_simulation_time - 1000
interval_data = 1000

sm = [ 0 ]
sf = [ 0 ]

t_opt = [[0.5,0.8],[0.8,0.5]]

omega = 0.05

nrep = 10

# how much temperature should change in a particular generation
temp_change = [0.25, 0.5] #np.linspace(start = 0, stop=1.0,num=60)
mu_t = [0.1]
mu_d = [0.01]
mu_slope = [0,0.01]

init_time_threshold = 6

# aim for roughly 25000 lines
skip = 1#int(round(max_simulation_time / 25000, 0))

if skip < 1:
    skip = 1

for mu_d_slope_i in mu_slope:
    for mu_di in mu_d:
        for mu_ti in mu_t:
            for change_i in temp_change:
                for rep_i in range(0,nrep):
                    for sm_i in sm:
                        for sf_i in sf:
                            if sf_i != sm_i:
                                continue

                            for t_opt_i in t_opt:
                                t_opt_f = t_opt_i[0]
                                t_opt_m = t_opt_i[1]

                                # get the basename
                                ctr += 1
                                base_name_i = base_name + "_" + str(ctr)

                                print(f"{exe_name} " \
                                        f"{base_name_i} " \
                                        f"{max_simulation_time} " \
                                        f"{sm_i} " \
                                        f"{sf_i} " \
                                        f"{omega} " \
                                        f"{omega} " \
                                        f"{t_opt_m} " \
                                        f"{t_opt_f} " \
                                        f"{skip} " \
                                        f"{change_i} " \
                                        f"{simulation_time_change} " \
                                        f"{interval_data} " \
                                        f"{mu_ti} " \
                                        f"{mu_di} " \
                                        f"{mu_d_slope_i} " \
                                        f"{init_time_threshold} " \
                                        )
