# coupled-systems-generator
Toolbox of functions that generate time-series of coupled systems of your choice.
This version is designed for Matlab/Octave.

Completed:
* adjacencygen.m: quickly build the adjacency matrix for networks of interest (linear, parallel, wheatstone, distance-3)
* coupledlogistic.m:  build the time-series for a system of coupled Logistic maps

Under development:
* coupledlogisterr.m: build the time-series for a system of coupled Logistic maps with a stochastic dynamical error

Note: for all the systems generated the user also gets the weighted laplacian and its pseudo-inverse.

For specific document of each item type in the Matlab/Octave command window: help [function-name]

--------------------------------
(C) Arthur Valencio* and Dr Murilo S. Baptista

ICSMB, University of Aberdeen    

(AV thanks to CNPq (Brazil) for a scholarship)

This toolbox is available free and with no warranty. Use it at your own risk.
