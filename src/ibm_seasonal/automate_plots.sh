#!/usr/bin/env bash

find . -iname "sim_seasonal_*" -print0 | xargs -0 -I % -P6 ./plot_output.r %
