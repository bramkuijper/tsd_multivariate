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
max_simulation_time = 50000

sm = [ 0.95]
sf = [ 0.95]

t_opt_f = [ -0.1, -0.5 ]

freq = [ 1, 0.1, 0.01 ]

for sm_i in sm:
    for sf_i in sf:
        if sf_i != sm_i:
            continue

        for t_opt_fi in t_opt_f:
            for freq_i in freq:

                # get the basename
                ctr += 1
                base_name_i = base_name + "_" + str(ctr)

                print(f"{exe_name} {base_name_i} {max_simulation_time} {sm_i} {sf_i} {t_opt_fi} {freq_i}")



