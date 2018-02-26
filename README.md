# coupled-logistic-maps
Functions for generating the time-series from networks of coupled logistic systems

This version is written for Matlab/Octave.

--------------------------------
Functions:

* adjacencygen.m: quickly build the adjacency matrix for networks of interest (linear, parallel, wheatstone, distance-3)
* coupledlogistic.m:  build the time-series for a system of coupled Logistic maps
* noisycoupledlogistic.m: build the time-series for a system of coupled Logistic maps with a stochastic noise in the dynamical process

Coupling options:

* diffuse: linear diffusive coupling with the adjacent node (i.e., of the form (x(j)-x(i), j representing the neighbour node),
* kaneko: Kaneko coupling with the adjacent node (i.e., of the form f(x(j)), j representing the neighbour node, f(x) the logistic function)

---------------------------------

For specific document of each item type in the Matlab/Octave command window: help [function-name]

--------------------------------
(C) Arthur Valencio* and Dr Murilo S. Baptista

ICSMB, University of Aberdeen    

* Support: CNPq (Brazil)

This toolbox is available as is, without any warranty. Use it at your own risk.

--------------------------------

If useful, please cite:

Arthur Valencio and Murilo S. Baptista. Coupled Logistic Maps: functions for generating the time-series from networks of coupled logistic systems. Open source codes for Matlab. 2018. Available at: https://github.com/artvalencio/coupled-logistic-maps.

Bibtex entry:
@misc{cami, author={Valencio, Arthur and Baptista, Murilo da Silva}, title={Coupled Logistic Maps: functions for generating the time-series from networks of coupled logistic systems}, note={Open source codes for Matlab. Available at \url{https://github.com/artvalencio/coupled-logistic-maps}.}, year={2018} }
