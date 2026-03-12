#!/usr/bin/env bash

arg="$1"
if [ -z "$1" ]; then
    arg="sim_seasonal_tsd*"
fi

find . -iname "$arg" -print0 | xargs -0 -I % -P6 ./plot_output.r %

rm -rf Rplots.pdf
