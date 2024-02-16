# T4.1-Tools-for-Ancillary-Services-Procurement-in-day-ahead-operation-planning-of-Distrib.-Networks
ATTEST-WP4-T4.1-Tools for Ancillary Services Procurement in day-ahead operation planning of distribution networks

This folder contains the tractable Sequential Linearization Algorithm (the GEM!), implemented in Julia, for flexibility procurement in active
distribution systems.

Steps to run the Tractable Tool

1. Open 'main_sc_milp.jl' file. This file contains all the packages, functions and main algorithm. 
2. Select the network input file as well the output file for the respective network. 
3. Select the flexible options that user wants to include in the optimization model
4. Run the file. 

The description of algorihtm is presented in the following paper. Please cite this paper if you use this algorithm for your work. Thanks!

M. Usman and F. Capitanescu, "A Novel Tractable Methodology to Stochastic Multi-Period AC OPF in Active Distribution Systems Using Sequential Linearization Algorithm," in IEEE Transactions on Power Systems, vol. 38, no. 4, pp. 3869-3883, July 2023, doi: 10.1109/TPWRS.2022.3197884
