#!/usr/bin/env python3

# plot the division of labor simulations

import pandas as pd
import itertools
import subprocess 
import math, string
import argparse
import numpy as np
import sys, re, os.path, pprint
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm

# import stats for kernel density estimation
from scipy import stats


#########################################
#           check where data ends
#########################################
def find_parameter_linenum(filename):

    # get all the lines
    f = open(filename)
    fl = f.readlines()
    f.close()

    # find the first line from below that
    # starts with a number. That is the last
    # line of the data
    frange = list(range(0,len(fl)))
    frange.reverse()

    # ok loop over all lines (backwards)
    for line_num in frange:
        if re.match("^\d.*",fl[line_num]) is not None:
            return(line_num+1)

    return(len(fl))


#########################################
#           read in the data
#########################################

if len(sys.argv) < 2:
    print("please provide a filename")
    sys.exit()

line_num_params = find_parameter_linenum(sys.argv[1])

dat = pd.read_csv(sys.argv[1],
        nrows=line_num_params-1,
        sep=';',
        index_col=False)

#########################################

#           make figure

#########################################

def make_panel(title, gs, row_ctr, panel_dict, data):
    
    panel_title = title

    ax = plt.subplot(gs[row_ctr,0])

    if type(panel_dict["x"]) != type(""):
        raise "X axis column name " + pprint.pformat(panel_dict["x"]) + " is not a string."

    if panel_dict["x"] not in data.columns.values:
        raise "Cannot find column " + panel_dict["x"] + " in dataset"

    # see whether there are labels to make a legend
    legend_labels = False

    if len(panel_dict["y"]) > 1 and "legend_labels" in panel_dict.keys():

        # check whether the list of similar length as the list of 
        # y variables
        if len(panel_dict["labels"]) != len(panel_dict["y"]):
            raise "Length of list of labels " + \
                    pprint.pformat(panel_dict["labels"]) + \
                    " differs from length of list of y axis variables " + \
                    ppprint.pformat(panel_dict["y"])

        legend_labels = True

    if type(panel_dict["y"]) == type(""):
        panel_dict["y"] = [panel_dict["y"]]

    # multiple y values?
    if type(panel_dict["y"]) == type([]):

        # go through each of the y values
        for i, y_val in enumerate(panel_dict["y"]):

            ax.plot(
                    data[panel_dict["x"]]
                    ,data[y_val]
                    ,label=panel_dict["legend_labels"][i] if \
                            legend_labels else y_val
                    )
    else: # just 
        raise "Y variable container for " + title + " is not a string or list, " + \
                pprint.format(panel_dict["y"])

    # should we add a legend
    if "legend" in panel_dict.keys() and type(panel_dict["legend"]) == \
            type(True) and panel_dict["legend"]:
        ax.legend()

    # only print tick labels in case it is explicitly
    # stated
    if "xticks" not in panel_dict.keys() or not panel_dict["xticks"]:
        ax.set_xticklabels(labels=[])
    else:
        ax.set_xlabel(panel_dict["xlabel"])

    if "xlim" in panel_dict.keys() and \
        type(panel_dict["xlim"]) == type([]):
        ax.set_xlim(panel_dict["xlim"])
        
    if "ylim" in panel_dict.keys() and \
        type(panel_dict["ylim"]) == type([]):
        ax.set_ylim(panel_dict["ylim"])
        
    ax.set_ylabel(panel_title)
        
    


def plot_from_dict(plot_dict, data, sim_file):

    # initialize the figure
    fig = plt.figure(figsize=(5,5 * len(plot_dict)))
    
    # generate the grid of the graph
    # see: 
    widths = [ 1 ]
    heights = [ 1 for x in range(0, len(plot_dict))]
    numrows = len(heights)
    numcols  = len(widths)

    # make the grid
    gs = gridspec.GridSpec(
            nrows=numrows,
            ncols=numcols,
            width_ratios=widths,
            height_ratios=heights)

    # go through the rows
    for iter_i, (key, panel_dict) in enumerate(plot_dict.items()):
        make_panel(title = key
                ,gs = gs
                ,row_ctr = iter_i
                ,panel_dict = panel_dict
                ,data = data
                )

    format = "png"

    filename = os.path.join(
            os.path.dirname(sim_file),
            "graph_" + os.path.basename(sim_file) + "." + format
            )

    plt.savefig(
            fname=filename
            ,format=format
            ,dpi="figure"
            ,bbox_inches="tight")


traitdict = {
        "Sex ratio": {
            "y" : ["sr1","sr2"],
            "x" : "generation",
            "legend" : True,
            "ylim" : [-.05,1.05]
        },
        "Variance in sex ratio": {
            "y" : ["varsr1","varsr2"],
            "x" : "generation",
            "legend" : True
            },
        "Dispersal" : {
            "y" : ["df","dm"],
            "x" : "generation",
            "legend" : True,
            "ylim" : [-.05,1.05]
            },
        "Digging" : {
            "x" : "generation",
            "y" : "b",
            "ylim" : [-0.05,1.05],
            "xticks" : True,
            "xlabel" : "generation"
            }
}




plot_from_dict(
        plot_dict = traitdict
        ,data = dat
        ,sim_file = sys.argv[1])
