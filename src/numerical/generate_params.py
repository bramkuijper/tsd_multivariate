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
init_dfdm = [[0.5,0.5]]
init_b = [0]
init_sr = 0.5
init_vf = 1.0
init_vm = 1.0
init_lambda = 1.0

# sigma1->2, sigma2->1
sigma = []

# vary risk of having the cold environment
risk = list(np.arange(0.01,0.99,0.01))

# environmental autocorrelation (within a patch)
autocorr = list(np.arange(0.01,0.99,0.01))

for risk_i in risk:
    for autocorr_i in autocorr:
        s12 = (1 - risk_i) * (1 - autocorr_i)
        s21 = risk_i * (1 - autocorr_i)

        sigma += [[s12, s21]]

# [eul_d, eul_sr, eul_b]
#euls = [[0,0.01,0],[0.01,0.01,0],[0.01,0.01,0.01],[0,0.01,0.01]]
euls = [[0,0.01,0.01]]


survival_vf2 = [0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]

survival_values = []

for surv_i in survival_vf2:
    survival_values += [[[1.0,1.0],[1.0,surv_i]]]

burrow_mod_survival = [1]

date = datetime.datetime.now()
base_name = "iter_tsd_multivariate_" + f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

ctr = 0

exe = "./tsd_multivariate.exe"

for nf_nm_i in nf_nm:
    for init_dfdm_i in init_dfdm:
        for init_b_i in init_b:
            for sigma_i in sigma:
                for burrow_mod_survival_i in burrow_mod_survival:
                    for survival_val_i in survival_values:
                        for eul_i in euls:
                            base_name_i = base_name + "_" + str(ctr)

                            print("echo " + str(ctr))

                            ctr+=1

                            print(exe + " \t" +
                                    f"{survival_val_i[0][0]} " + 
                                    f"{survival_val_i[0][1]} " + 
                                    f"{survival_val_i[1][0]} " + 
                                    f"{survival_val_i[1][1]} " + 
                                    f"{init_dfdm_i[0]} " +
                                    f"{init_dfdm_i[1]} " +
                                    f"{init_b_i} " +
                                    f"{init_sr} " +
                                    f"{init_sr} " +
                                    f"{sigma_i[0]} " +
                                    f"{sigma_i[1]} " +
                                    f"{init_lambda} " +
                                    f"{burrow_mod_survival_i} " +
                                    f"{nf_nm_i[0]} " +
                                    f"{nf_nm_i[1]} " +
                                    f"{eul_i[0]} " +
                                    f"{eul_i[1]} " +
                                    f"{eul_i[2]} " +
                                    base_name_i
                                    )

            






