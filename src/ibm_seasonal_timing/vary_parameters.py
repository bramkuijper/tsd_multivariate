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

ctr = 0
max_simulation_time = int(20000000)

sm = [ 0 ]
sf = [ 0 ]

t_opt_f = [0.5]
t_opt_m = [0.5,0.8]

# width of the selection function
omega = 0.01

nrep = 3

temp_change = [0.5]
mu_t = [0.05]

temp_error = [0,0.01,0.05,0.5,1.0]

# aim for roughly 25000 lines
skip = max([int(round(max_simulation_time / 25000, 0)),1])

for change_i in temp_change:
    for rep_i in range(0,nrep):
        for sm_i in sm:
            for sf_i in sf:
                if sf_i != sm_i:
                    continue

                for t_opt_fi in t_opt_f:
                    for t_opt_mi in t_opt_m:
                        for temp_error_i in temp_error:

                            # get the basename
                            ctr += 1
                            base_name_i = base_name + "_" + str(ctr)

                            print(f"{exe_name} {base_name_i} " \
                                    f"{max_simulation_time} " \
                                    f"{sm_i} " \
                                    f"{sf_i} " \
                                    f"{omega} " \
                                    f"{omega} " \
                                    f"{t_opt_mi} " \
                                    f"{t_opt_fi} " \
                                    f"{skip} " \
                                    f"{change_i} " \
                                    f"{mu_t[0]} " \
                                    f"{temp_error_i} " \
                                  )
