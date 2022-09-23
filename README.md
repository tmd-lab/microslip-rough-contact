# microslip-rough-contact
An open source physics-based friction model for the dynamics of jointed structures.

This code is provided to aid in research, no guarantee is made for the accuracy of the results when applied to other structures. If this code is useful, please cite the journal paper:


@article{porterTowardsAPredictive,
	title = {Towards a Predictive, Physics-Based Friction Model for the Dynamics of Jointed Structures},
	journal = {Mechanical Systems and Signal Processing},
	author = {Porter, J.~H. and Brake, M.~R.~W.},
	year = {Under Review},
}


**Insert Reference Info for Final Paper Here**

## Running the code

All computations for the project were run with MATLAB/2020a.

1. Download the needed *.mat files from **[insert url]**.
2. Most of the time, just run the top level function from the list below. Additional details are noted below.

### SURFACE

1. Open 'SURFACE/asp_id_ltw2.m'
2. Set the scans to use on lines 15 and 16. 
3. Settings can be changed under "Process Settings" for a variety of the steps. Code comments describe the different settings. These settings are as used in the paper. The only setting changed for the paper was "settings.cleanLength" which was varied from 2 to 75 to test different erosion sizes. Most simulations used a size of 2.
4. Run the full code. 
5. A temporary file is written out on line 149 in case one wants to restart after that point. If restarting, the code will need to be modified to load the saved data.
6. Make sure to save the final outputs. Save lines are currently commented out at the bottom of the script. 

### RQNMA

1. RQNMA simulations are run with "FJSIMS/main_rqnm.m", which contains the function "main_rqnm()". Passing a number into the "main_rqnm" function sets the choice of parameters. The different parameter cases are defined in "FJSIMS/rqnm_pars.m"
2. To create a new case, copy an existing case such as case 4 that defines all of the parameters needed for the input. If only modifying a few parameters from an existing case, you can also reference a previous case as shown in some examples (e.g., case 12). parameters are:
   1. mu - friction coefficient
   2. ElasticModulus - elastic modulus (Pa)
   3. nu - Poisson's ratio
   4. mesoscale - true uses the defined mesoscale topology in the surface processing. false uses a flat interface.
   5. asp_fun - function handle for the asperity contact model. "ELLIPSOID_IWAN_DECOUPLE" is the decoupled Mindlin-Iwan Constant model from the paper. "ELLIPSOID_IWAN_FIT_(DE)COUPLE" refers to the (de)coupled Mindlin-Iwan Fit model from the paper.
   6. Nqp_radius - number of quadrature points in the radial direction for the asperity contact. Commented out lines have both this and the asperity function with the values used in the paper. 
   7. alpha - Asperity rotation angle (see Section 2.3.1 of paper). 
   8. useSphere - true sets the contact models to be spheres, false uses ellipsoids. 
   9. Nqp_heights - number of quadrature points for integrating over initial asperity gaps. 
   10. meshName - "UROM232" refers to the mesh used in the paper and the file provided. "zte152" refers to a smaller mesh that is also provided and can be used for debugging. 
   11. meshNameRelative - not relevant, just leave as ''. 
   12. output_name - Name for file output. 
   13. load_initial - this should not be relevant anymore. 
   14. repeatLoop - number of times to repeat the cycle from +q to -q to +q in RQNMA solver. 
   15. Nhp - number of quadrature points per segment in RQNMA solver. 
3. Set "nots_min_cores" at the top of "main_rqnm.m". If there are fewer than this number of cores, it assumes it is running locally and does not need to set variables for running on a cluster. If there are more cores, additional steps are taken. 
4. If running on a cluster (with more cores than "nots_min_cores"), set the work directory on line 40 of "main_rqnm.m". This should be the directory where the code is saved and where it is being run. 
5. Define the amplitude levels of interest with "Qamps" on line 431. These are modal amplitudes as defined by RQNMA that the solution will be found at. 
6. Other properties for RQNMA can be set at around line 435-444, but these should not need to be changed. 
7. Save the output. Currently, it is set to only save on line 456 if running on a cluster. Make sure to put a break point or change this otherwise the function call may exit without saving the output data. 
8. Run the code. An example of slurm file is provided only for EPMC if running on a cluster, but could easily be modified for RQNMA runs. 

### EPMC

The procedure for running EPMC is very similar to that of RQNMA. 

1. EPMC simulations are run with "EPMC_SIMS/main_epmc.m", which contains function "main_epmc()". Passing a number into the function sets the choice of parameters as defined in "EPMC_SIMS/epmc_pars.m". 
2. To create a new case, copy an existing case such as case 11 that contains all of the parameters already. The parameters are: 
   1. Sys - Yield Strength in Pa
   2. mu - friction coefficient (set to a large value if using CEB friction coefficient). 
   3. ElasticModulus - elastic modulus of the asperities. This is just used for the asperities for the provided meshes. The provided value was also used in the model creation for the provided models. 
   4. nu - Poisson's ratio
   5. sliptype - slip type to be used for the model. 1 is for just using friction coefficient mu, 3 is CEB friction and limited by mu. 
   6. mesoscale - true to use the measured mesoscale topology, false to use a flat interface.
   7. unloadModel - setting previously used to change between different unloading models for spheres after plastic deformation. Code is only provided for the option 'brake'
   8. useSphere - calculates some parameters for contact as if the asperities were spheres instead of ellipsoids. EPMC implementation with plasticity only supports spheres.
   9. Et - tangent modulus after yielding of the asperities, units of Pa. 
   10. meshName - 'UROM232' was used for the paper. 'zte152' is also provided and can be used as a smaller model for debugging. 
   11. output_name - output filename to save after every completed solution step of continuation. Results are automatically put in a folder called 'EPMC_SIMS/Results'
   12. load_initial - use False to start a new simulation, set to true along with 'prevInd' and 'loadName' to restart a simulation.
   13. prevInd - index to use from a loaded file when restarting. 
   14. loadName - file to load to restart, needs to start with 'Results/' if loading a file that was saved in the Results folder. 
   15. arcSettings.da - initial step size to take in continuation. 
   16. arcSettings.dsmax - maximum step size to allow in continuation. 
   17. arcSettings.dsmin - minimum step size to allow in continuation. 
3. Set "nots_min_cores" at the top of "main_epmc.m". If there are fewer than this number of cores, it assumes it is running locally and does not need to set variables for running on a cluster. If there are more cores, additional steps are taken. 
4. If running on a cluster (with more cores than "nots_min_cores"), set the work directory on line 40 of "main_epmc.m". This should be the directory where the code is saved and where it is being run. 
5. Define the amplitude levels of interest with As (start amplitude on a log scale) and Ae (end amplitude on a log scale). 
6. Additional settings for continuation and the solver can be found on line 638, but these should not need to be changed. 
7. Run the code. An example of slurm file is provided as 'EPMC_SIMS/run_epmcX.slum' which runs the epmc simulation with parameters set by the environment variable X. Calling 'source submitX.sh' will submit EPMC runs based on the do loop in the file 'submitX.sh'.

### Key Plots

PLOTS/elastic_bb.m produces the backbone plots for the entire paper:

1. Set the appropriate level of viscous damping for the modes of interest on line 15. 
2. Select the plot set of interest with the variable 'plot_set' (line 17). 
3. Copy an existing case. Case 6 provides an example for both RQNMA and EPMC.
4. For RQNMA simulations only one filename is needed for the data. For EPMC simulations, two files are provided, one for the data saved before the simulation starts and another for the iteration results.
5. All of the other parameters control formatting the plot. 
6. Data is provided for plot_set = 6 as an example.

PLOTS/epmc_hyst_plots.m produces the plots of local hysteresis loops for the paper:

1. Set plot_set on line 35 to select which case of loaded data to use. Data is provided for plot_set = 2. 
2. Set the run_name for saving output. 
3. Load the vector of harmonic displacements of interest into the variable Uwxa. This vector is a stack of the harmonic displacements, the frequency, the damping from EPMC, and the amplitude. 
4. plot_indices sets which elements to show plots of. 
5. paper_plot controls the output indices for the paper format of figures. 
6. plot_letters on line 203 controls the letters used for labeling the mesh. 

## .mat Files to Download

| Filename | Destination Folder | Description|
|:---           |:---                            |:---              |
| R05A_Before_R1.mat | SURFACE/SCANS/LongTermWear/ | Scan data for side A |
| R05B_Before_R1.mat | SURFACE/SCANS/LongTermWear/ | Scan data for side B |
| combined_14sept21_R1.mat.  | SURFACE/OUT/ | Baseline processed scans for most simulations|
| combined_2dec21_R1_erode($X).mat.  | SURFACE/OUT/ | Processed scan properties for different erosion sizes. Replace ($X) with 25, 50, or 75. |
| ROM_PD_152ELS.mat | FJSIMS/ROMS/ | Reduced order model with 152 ZTEs. Not used for paper, but is a smaller model that can be helpful for debugging. | 
| ROM_U_232ELS.mat  | FJSIMS/ROMS/ | Reduced order model with 232 ZTEs. This version is used in the paper. |
|rqnm_paper18.mat | FJSIMS/Results/PAPER/alphaRuns/ | Numerical Results for elastic model with spheres used in plotting examples |
|epmc_paper_run10_v2_pre.mat | EPMC_SIMS/Results/PAPER/Run10/ | Numerical results for initial state of solver for EPMC with elastic spheres used in plotting examples | 
|epmc_paper_run10_v2_iter.mat | EPMC_SIMS/Results/PAPER/Run10/ | Numerical results for EPMC solver with elastic spheres used in plotting examples |
|LTW_Mode1Trimmed_27Sept2021 |EXPERIMENTAL_DATA/ | Experimental backbone data for comparison to numerical results, mode 1 | 
|LTW_Mode2Trimmed_27Sept2021 |EXPERIMENTAL_DATA/ | Experimental backbone data for comparison to numerical results, mode 2 | 

Notes:

1. Additional scans of the same set of interfaces ending in 'R2' are also provided and should be in a format that can be used with the surface processing scripts. These were not used in this project.


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
* EPMC_SIMS/submitX.sh % used to submit runs to slurm system with 'source submitX.sh'
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


## References

Please see references in the journal paper above. 

Additional Notes:

Many nonlinear solver routines are used from: https://github.com/Nidish96/octave-jim

Distinguishable colors for plots are generated use code from: https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors


## Acknowledgments

The authors would like to thank Nidish Narayanaa Balaji, Iyabo Lawal, and Scott A. Smith for conducting experiments and sharing data for the Brake-Reuss Beam. 

Funding: 
This material is based upon work
	supported by the U.S. Department of Energy, Office of
	Science, Office of Advanced Scientific Computing
	Research, Department of Energy Computational Science
	Graduate Fellowship under Award Number(s) DE-SC0021110.
  This work was supported in part by the Big-Data Private-Cloud Research
  Cyberinfrastructure MRI-award funded by NSF under grant CNS-1338099 
  and by Rice University.
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
