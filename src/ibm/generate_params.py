#!/usr/bin/env python3

# ptyhon script that 
# generates a batch file with all the different 
# simulations that need to run for all parameter
# combinations

import numpy as np
import socket
import datetime
import sys
import copy

nf_nm = [[10,10]]

init_d = [[0.2,0.2]]

init_b = [0.0]

init_sr = [[0.5,0.5]]

vf = [[0.5,1.0]]

vm = [[1.0,1.0]]

s = [[0.5,0.5]]

spatial = [1]

clutch_max = 25

mu_sr = 0.01

mu_b = [0, 0.01]

mu_d = [[0,0]]

sdmu = 0.01

burrow_mod_survival = [0,1]

date = datetime.datetime.now()
base_name = "sim_tsd_multivariate_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

nrep = 5

ctr = 0

exe = "./tsd_multivariate.exe"

max_generations = 50000

for replicate_i in range(0,nrep):
    for nf_nm_i in nf_nm:

        n_patches_i = 5000 / nf_nm_i[0]

        for init_d_i in init_d:
            for init_b_i in init_b:
                for init_sr_i in init_sr:
                    for vf_i in vf:
                        for vm_i in vm:
                            for s_i in s:
                                for mu_b_i in mu_b:
                                    for mu_d_i in mu_d:
                                        for spatial_i in spatial:
                                            for burrow_mod_survival_i in burrow_mod_survival:
                                                base_name_i = base_name + "_" + str(ctr)

                                                print("echo " + str(ctr))

                                                ctr+=1

                                                print(exe + " \t" +
                                                        f"{nf_nm_i[0]} " + 
                                                        f"{nf_nm_i[1]} " + 
                                                        f"{n_patches_i:n} " +
                                                        f"{clutch_max:n} " +
                                                        f"{init_d_i[0]} " +
                                                        f"{init_d_i[1]} " +
                                                        f"{init_b_i} " +
                                                        f"{init_sr_i[0]} " +
                                                        f"{init_sr_i[1]} " +
                                                        f"{vf_i[0]} " +
                                                        f"{vf_i[1]} " +
                                                        f"{vm_i[0]} " +
                                                        f"{vm_i[1]} " +
                                                        f"{burrow_mod_survival_i} " +
                                                        f"{s_i[0]} " +
                                                        f"{s_i[1]} " +
                                                        f"{spatial_i} " +
                                                        f"{mu_sr} " +
                                                        f"{mu_b_i} " +
                                                        f"{mu_d_i[0]} " +
                                                        f"{mu_d_i[1]} " +
                                                        f"{sdmu} " +
                                                        f"{max_generations} " +
                                                        base_name_i
                                                        )

                                






