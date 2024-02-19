#!/bin/bash
g++ $(root-config --cflags) -o plasma_deceleration_plots plasma_deceleration_plots.C $(root-config --libs)
./plasma_deceleration_plots
