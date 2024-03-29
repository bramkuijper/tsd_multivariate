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

# function to calculate a risk value
# and autocorrelation into frequencies
# risk is given by freq of envt 2, s12/(s21 + s12)
# autocorr is given by 1.0 - s12 - s21

def riskAutocorr2freq(risk, autocorr):

    s12 = risk * (1.0 - autocorr)
    s21 = 1.0 - autocorr - s12

    return([s12,s21])

nf_nm = [[5,5]]

init_d = [[1.0,1.0],[0.2,0.2]]

init_b = [0.0]

init_sr = [[0.5,0.5]]

vf = [[1.0,0.9]]

vm = [[1.0,1.0]]

d_cue_is_envt = 0

vf_pert = [[1.0,0.3]]
vm_pert = [[1.0,1.0]]

#r = list(np.linspace(0.01,0.99,10))
#p1 = list(np.linspace(0.01,0.99,10))
#
#s = []
#
#for r_i in r:
#    for p1_i in p1:
#        s += [[r_i/(2 * p1_i),r_i/(2 * (1.0 - p1_i))]]

s_orig = [riskAutocorr2freq(0.875, 0.75),riskAutocorr2freq(0.1, 0.75)]

s_pert = [riskAutocorr2freq(0.1, 0.75),riskAutocorr2freq(0.875, 0.75)]

spatial = [1]

clutch_max = 25

mu_sr = 0.01

mu_b = [0.0,0.1]

mu_d = [[0,0],[0.01,0.01]]

sdmu = 0.01

burrow_mod_survival = [1]

date = datetime.datetime.now()
base_name = "sim_tsd_multivariate_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

nrep = 5

ctr = 0

# whether jobs should be run in the background
run_in_background = False

# never run background jobs on cluster
hostname = socket.gethostname()
if "athena" in hostname:
    run_in_background = False

exe = "./tsd_multivariate.exe"

max_generations = 50000
time_perturb = 20000#max_generations/2.0

for replicate_i in range(0,nrep):
    for nf_nm_i in nf_nm:

        n_patches_i = 5000 / nf_nm_i[0]

        for init_d_i in init_d:
            for init_b_i in init_b:
                for init_sr_i in init_sr:
                    for vf_i in vf:
                        for vm_i in vm:
                            for vf_pert_i in vf_pert:
                                for vm_pert_i in vm_pert:
                                    for s_i in s_orig:
                                        for s_pert_i in s_pert:
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
                                                                    f"{vf_pert_i[0]} " +
                                                                    f"{vf_pert_i[1]} " +
                                                                    f"{vm_pert_i[0]} " +
                                                                    f"{vm_pert_i[1]} " +
                                                                    f"{burrow_mod_survival_i} " +
                                                                    f"{s_i[0]} " +
                                                                    f"{s_i[1]} " +
                                                                    f"{s_pert_i[0]} " +
                                                                    f"{s_pert_i[1]} " +
                                                                    f"{spatial_i} " +
                                                                    f"{d_cue_is_envt} " +
                                                                    f"{mu_sr} " +
                                                                    f"{mu_b_i} " +
                                                                    f"{mu_d_i[0]} " +
                                                                    f"{mu_d_i[1]} " +
                                                                    f"{sdmu} " +
                                                                    f"{max_generations} " +
                                                                    f"{time_perturb:n} " +
                                                                    base_name_i + " " +
                                                                    ("&" if run_in_background else "")
                                                                    )

                                            






