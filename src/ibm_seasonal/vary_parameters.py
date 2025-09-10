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
max_simulation_time = 150000

sm = [ 0.95]
sf = [ 0.95]

t_opt_f = [ 0.1 ]
t_opt_m = 0.25

freq = [ 1, 0.1, 0.01 ]

omega = 0.05

nrep = 3

for rep_i in range(0,nrep):
    for sm_i in sm:
        for sf_i in sf:
            if sf_i != sm_i:
                continue

            for t_opt_fi in t_opt_f:
                for freq_i in freq:

                    # get the basename
                    ctr += 1
                    base_name_i = base_name + "_" + str(ctr)

                    print(f"{exe_name} {base_name_i} " \
                            f"{max_simulation_time} " \
                            f"{sm_i} " \
                            f"{sf_i} " \
                            f"{omega} " \
                            f"{omega} " \
                            f"{t_opt_fi} " \
                            f"{t_opt_m} " \
                            f"{freq_i}")



