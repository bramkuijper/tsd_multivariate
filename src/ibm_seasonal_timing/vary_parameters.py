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

sm = [ 0.5]
sf = [ 0.5]

t_opt_f = [ 0.5 ]
t_opt_m = 0.8

omega = 0.05

nrep = 3

temp_change = [0.1]
mu_t = [0,0.05]

# aim for roughly 25000 lines
skip = int(round(max_simulation_time / 25000, 0))

if skip < 1:
    skip = 1

for mu_ti in mu_t:
    for change_i in temp_change:
        for rep_i in range(0,nrep):
            for sm_i in sm:
                for sf_i in sf:
                    if sf_i != sm_i:
                        continue

                    for t_opt_fi in t_opt_f:

                        # get the basename
                        ctr += 1
                        base_name_i = base_name + "_" + str(ctr)

                        print(f"{exe_name} {base_name_i} " \
                                f"{max_simulation_time} " \
                                f"{sm_i} " \
                                f"{sf_i} " \
                                f"{omega} " \
                                f"{omega} " \
                                f"{t_opt_m} " \
                                f"{t_opt_fi} " \
                                f"{skip} " \
                                f"{change_i} " \
                                f"{mu_ti}")
