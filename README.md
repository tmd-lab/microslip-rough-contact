# microslip-rough-contact
An open source physics-based friction model for the dynamics of jointed structures.

This code is provided to aid in research, no guarantee is made for the accuracy of the results when applied to other structures. If this code is useful, please cite the journal paper:

**Insert Reference Info for Final Paper Here**

## Running the code

1. Download the needed *.mat files from **[insert url]**.
2. Run the top level function from the list below.



## .mat Files to Download

* file1.mat to FOLDER/ % Description

## Top Level Scripts

### SURFACE
SURFACE folder contains script for generating surface parameters for use in 
models. Some surface parameter files are provided, so this is not needed to 
rerun existing dynamic simulations. 

* SURFACE/asp_id_ltw2.m % Processing of surface scan data into surface model     

### FJSIMS
FJSIMS folder contains RQNMA simulations. 

* FJSIMS/main_rqnm.m % script for running RQNMA simulations
* FJSIMS/rqnm_pars.m % script for setting parameters used for RQNMA simulations, this file is not directly run

### EPMC_SIMS
EPMC_SIMS folder contains EPMC simulations. 

* EPMC_SIMS/main_epmc.m % main file to run to do an EPMC simulation. Pass in argument for which parameter number (default 0)
* EPMC_SIMS/epmc_pars.m % parameters for running specific EPMC simulations, this file is not directly run.
* EPMC_SIMS/submitX.sh % used to submit runs to slurm system with "source submitX.sh"
* EPMC_SIMS/run_epmcX.slurm % file for running an individual slurm job

### TESTS
Tests are provided for several individual functions and 
verification of derivatives. Multiple cases for
derivative checking have to be run manually (see
comments). In some cases, derivatives do not match
numerical differentiation. Comments note the subtle
differences, frequently due to taking the derivative on
only one side of slope discontinuities. Some tests
produce plots relevant to the paper. 

* TESTS/test_elliptical_contact.m
* TESTS/test_asp_fun.m % Derivative verification for the elastic ellipsoid models and some verification plots
* TESTS/test_plastic_contact.m 
* TESTS/test_plastic_contact_unload.m % Plots related to code development and unloading models
* TESTS/verification_normal_plasticity.m % Figures can be compared to unloading paper to verify normal model implementation.
* TESTS/EPMC/test_harmonic_asp_deriv.m % Harmonic Balance Derivatives checking
* TESTS/EPMC/test_harmonic_rc_traction.m % Harmonic Balance Derivatives checking
* TESTS/test_asp_fit.m % Tests for asperity fitting routine
* TESTS/test_asp_id.m % Quick test of watershed routine to verify that it looks roughtly as expected.    

### PLOTS
PLOTS folder contains script for plotting final results. For backbones in the
paper, not all data is provided, but scripts show an example plot. 

* PLOTS/mindlin_compare_plot.m % Plots of tangential friction response versus Mindlin theory.  
* PLOTS/elastic_bb.m % Backbone plotting - not all results are provided for this. File includes plastic backbone plotting. Data is provided for plotting case 6 to serve as an example of the output plots. This example includes both an RQNMA and EPMC simulation.
* PLOTS/epmc_hyst_plots.m % Hysteresis Loop Plotting - requires manually running loops over subsets of indices to prevent overloading graphics card. Running the full script on linux may run some plot commands out of order and not give the expected results. Data is provided for the elastic beam example for plot_set=2
* PLOTS/plot_mesh.m % Plots Mesh
* PLOTS/plot_surf_scan.m % Plot Surface Scan data. 'red_plot=true' only plots a reduced set to check the figure before plotting everything.
* PLOTS/asperity_fit_plot.m % Plot of ellipsoid fitting an asperity

### EXPERIMENTAL_DATA
EXPERIMENTAL_DATA folder - files are placed here for experimental data.
 


## Acknowledgments

This material is based upon work
	supported by the U.S. Department of Energy, Office of
	Science, Office of Advanced Scientific Computing
	Research, Department of Energy Computational Science
	Graduate Fellowship under Award Number(s) DE-SC0021110.
This work was supported in part by the Big-Data Private-Cloud Research Cyberinfrastructure MRI-award funded by NSF under grant CNS-1338099 and by Rice University. %This is required for using NOTS: https://docs.rice.edu/confluence/display/CD/NOTS+Overview
The authors are thankful for the support of the National
Science Foundation under Grant Number 1847130.

This report was prepared as an account of
work sponsored by an agency of the United States
Government. Neither the United States Government nor
any agency thereof, nor any of their employees, makes
any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness,
or usefulness of any information, apparatus, product, or
process disclosed, or represents that its use would not
infringe privately owned rights. Reference herein to any
specific commercial product, process, or service by trade
name, trademark, manufacturer, or otherwise does not
necessarily constitute or imply its endorsement,
recommendation, or favoring by the United States
Government or any agency thereof. The views and
opinions of authors expressed herein do not necessarily
state or reflect those of the United States Government or
any agency thereof.

