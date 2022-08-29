% Script for Plotting Elastic Backbones and Comparing to Experiments

clear;

% set(0, 'defaultaxesinterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');


addpath('../ROUTINES/')


% Viscous damping to add to all sims
z_visc_opts = [0.087e-2, 0.034e-2]; % Mode 1,2

plot_set = 12;

%% Plot Data to Load etc

show_z0 = true;
SuppressLegend = 'on';
extra_info = '';
eliminateInds = [];
leg_down = 0; % Amount to move legend down. 

exp_mode = 1; % Default 1st mode

switch plot_set
    case 0
        %%%%%%%%% Just Plot the Experimental data
        NAMES = {};
        RQNMA_EPMC = [];
        
        Rindex = 6;
        
        legend_names = {};

        
        plot_markers = {'v-','h--', 'o:', '^-', 's--', 'd:', 'p-', '<--', '>:'};
        color_inds = [8, 15, 16, 7, 17, 18, 9, 19, 20];
        marker_size = [6, 6, 6, 6, 6, 6, 6, 6, 6];
        
        line_width = 3;

%         numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [1e-1, 1e2];
%         Wlim = [155, 180];
%         Zlim = [2e-4 1e-1];
        
        numTrim = [0, 3]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [2.8e-1, 1.35e1];
        Wlim = [175.5, 180];
        Zlim = [5e-4 10e-3];
        
        show_z0 = false;
        
    case 1
        %%%%%%%%% Model Convergence
        NAMES = {'../FJSIMS/Results/PAPER/rqnm_paper2.mat', '../FJSIMS/Results/PAPER/rqnm_paper3.mat', ...
                '../FJSIMS/Results/PAPER/rqnm_paper1.mat', ...
                {'../EPMC_SIMS/Results/PAPER/Initial/epmc_paper7_pre_initial.mat', '../EPMC_SIMS/Results/PAPER/Initial/epmc_paper7_iter70_initial.mat'}};
        RQNMA_EPMC = [1, 1, 1, 2];
        FIG_TITLE = 'Convergence: $\mu = 0.02$';
        
        Rindex = 6;
        
        legend_names = {'122 ZTE (RQNMA)', '232 ZTE (RQNMA)', '588 ZTE (RQNMA)', '232 ZTE (EPMC)'};

        
        plot_markers = {'v-','h--', 'o:', '^-', 's--', 'd:', 'p-', '<--', '>:'};
        color_inds = [8, 15, 16, 7, 17, 18, 9, 19, 20];
        marker_size = [6, 6, 6, 6, 6, 6, 6, 6, 6];
        
        line_width = 3;

%         numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [1e-1, 1e2];
%         Wlim = [155, 180];
%         Zlim = [2e-4 1e-1];
        
        numTrim = [0, 3]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [1e-1, 2e1];
        Wlim = [170, 180];
        Zlim = [2e-4 2e-2];
    
    case 2
        %%%%%%%%%% PAPER: Friction Coefficient
        NAMES = {'../FJSIMS/Results/PAPER/MuRuns/rqnm_paper4.mat', ...
                 '../FJSIMS/Results/PAPER/MuRuns/rqnm_paper5.mat', ...
                 '../FJSIMS/Results/PAPER/MuRuns/rqnm_paper6.mat', ...
                 '../FJSIMS/Results/PAPER/MuRuns/rqnm_paper7.mat', ...
                 '../FJSIMS/Results/PAPER/MuRuns/rqnm_paper8.mat', ...
                 '../FJSIMS/Results/PAPER/MuRuns/rqnm_paper9.mat', ...
                 '../FJSIMS/Results/PAPER/MuRuns/rqnm_paper10.mat'};
             
        RQNMA_EPMC = [1, 1, 1, 1, 1, 1, 1];
        FIG_TITLE = 'Varying $\mu$';
        
        Rindex = 6;
        
        legend_names = {'$\mu$=0.01', '$\mu$=0.02', '$\mu$=0.03', '$\mu$=0.05', ...
                        '$\mu$=0.10', '$\mu$=0.50', '$\mu$=0.90'};

        
        plot_markers = {'v-','h--', 'o:', '^-', 's--', 'd-', 'p:', '<--', '>:'};
        color_inds = [1:3, 5:8];
        marker_size = [6, 7, 7, 6, 7, 7, 7, 6, 6];
        
        line_width = 3;

%         numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [1e-1, 1e2];
%         Wlim = [155, 180];
%         Zlim = [2e-4 1e-1];
        
        numTrim = [0, 3]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [1e-1, 2.1e1];
        Wlim = [162, 180];
        Zlim = [0.7e-3 4e-2];

    case 3
        %%%%%%%%%% PAPER: Asperity Function Differences
        NAMES = {'../FJSIMS/Results/PAPER/MuRuns/rqnm_paper6.mat', ...
                 '../FJSIMS/Results/PAPER/ASPFUN/rqnm_paper11.mat', ...
                 '../FJSIMS/Results/PAPER/ASPFUN/rqnm_paper12.mat', ...
                 '../FJSIMS/Results/PAPER/ASPFUN/rqnm_paper13.mat'};
             
        RQNMA_EPMC = [1, 1, 1, 1, 1];
        FIG_TITLE = 'Varying Asperity Function, $\alpha=\pi/2$';
        
        Rindex = 6;
        
        legend_names = {'Ellipsoid Tangent', 'Ellipsoid MIC', ...
                        'Ellipsoid MIF', 'Ellipsoid MIF Coupled'};

        
        % First one for run 6 copied from previous figure.
        plot_markers = {'o:', 'p--', '<-', '>:', 'v-', 'h--', 'o-', '^--', 's--', 'd:', };
        color_inds = [3, 9, 11, 14:20];
        marker_size = [7, 7, 7, 7, 7, 7, 6, 6, 6, 6];
        
        line_width = 3;

%         numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [1e-1, 1e2];
%         Wlim = [155, 180];
%         Zlim = [2e-4 1e-1];
        
        numTrim = [0, 3]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [1e-1, 2.1e1];
        Wlim = [171, 180];
        Zlim = [7e-4 2e-2];
        
    case 4
        %%%%%%%%%% PAPER: Asperity Angle
        % Using Orthotropic Ellipsoid Tangent Model - not isotropic causes
        % mu to varying with rotation. 
%         NAMES = {'../FJSIMS/Results/PAPER/alphaRuns/rqnm_paper14.mat', ...
%                  '../FJSIMS/Results/PAPER/alphaRuns/rqnm_paper15.mat', ...
%                  '../FJSIMS/Results/PAPER/alphaRuns/rqnm_paper16.mat', ...
%                  '../FJSIMS/Results/PAPER/alphaRuns/rqnm_paper17.mat', ...
%                  '../FJSIMS/Results/PAPER/MuRuns/rqnm_paper6.mat'};
             
        % Using Isotropic Ellipsoid Fit Model Instead
        NAMES = {'../FJSIMS/Results/PAPER/alphaRuns/rqnm_paper18.mat' ... % Sphere
                 '../FJSIMS/Results/PAPER/alphaRuns/rqnm_paper19.mat' ... % Ellipsoid, alpha=0
                 '../FJSIMS/Results/PAPER/alphaRuns/rqnm_paper20.mat' ... % Ellipsoid alpha=pi/4
                 '../FJSIMS/Results/PAPER/alphaRuns/rqnm_paper21.mat' ... % Ellipsoid alpha=-pi/4
                 '../FJSIMS/Results/PAPER/ASPFUN/rqnm_paper13.mat'}; % pi/2
        
        RQNMA_EPMC = [1, 1, 1, 1, 1];
        FIG_TITLE = 'Varying $\alpha$';
        
        Rindex = 6;
        
        legend_names = {'Sphere', '$\alpha$=0', '$\alpha$=$\pi$/4', '$\alpha$=$-\pi$/4', '$\alpha$=$\pi$/2'};

        % Entry 5 copied to match previous figure
        plot_markers = {'s-', 'h--', 'o--', '^:', '>:', '>:'}; %'d:', 'p-', '<--', '>:'
        color_inds = [19, 16, 17, 18, 14];
        marker_size = [10, 7, 7, 7, 7];
        
        line_width = 3;
        
%         numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [1e-1, 1e2];
%         Wlim = [155, 180];
%         Zlim = [2e-4 1e-1];
        
        numTrim = [0, 3]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [1e-1, 2.1e1];
        Wlim = [171, 180];
        Zlim = [7e-4 2e-2];
        
    case 5
        %%%%%%%%%% IMAC 2022 Mesoscale
        NAMES = {
                 '../FJSIMS/Results/IMAC2022/Mesoscale/rqnm_paper24.mat', ...
                 '../FJSIMS/Results/IMAC2022/Mesoscale/rqnm_paper23.mat', ...
                 '../FJSIMS/Results/PAPER/MuRuns/rqnm_paper6.mat', ...
                 '../FJSIMS/Results/IMAC2022/Mesoscale/rqnm_paper22.mat', ...
                 '../FJSIMS/Results/IMAC2022/Mesoscale/rqnm_paper27.mat', ...
                 '../FJSIMS/Results/IMAC2022/Mesoscale/rqnm_paper25.mat'};%, ...
                 %'../FJSIMS/Results/IMAC2022/Mesoscale/rqnm_paper26.mat'};
             
        RQNMA_EPMC = [1, 1, 1, 1, 1, 1, 1];
        FIG_TITLE = 'Mesoscale';
        
        Rindex = 6;
        
        legend_names = {'Center Gap - 250 $\mu$m', 'Center Gap - 100 $\mu$m', ...
                        'Measured Mesoscale', 'Center Gap - 50 $\mu$m', ...
                        'Flat', ...
                        'Edge Gap - 50 $\mu$m'};%,    'Edge Gap - 250 $\mu$m'};

        
        plot_markers = {'v--','h-', 'o:', '^--', 's-', 'd-', 'p:', '<--', '>:'};
        color_inds = [1:3, 5:8];
        marker_size = [6, 7, 7, 6, 7, 7, 7, 6, 6];
        
        line_width = 3;

%         numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [1e-1, 1e2];
%         Wlim = [155, 180];
%         Zlim = [2e-4 1e-1];
        
        numTrim = [0, 3]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [.75e-1, 2.1e1];
%         Wlim = [142, 184]; % If including 250 edge gap
        Wlim = [153, 184];
        Zlim = [7e-4 4e-2];
    
    case 6
        %%%%%%%%%% IMAC 2022 Physics
        NAMES = {'../FJSIMS/Results/IMAC2022/Mesoscale/rqnm_paper27.mat', ...
                 '../FJSIMS/Results/PAPER/MuRuns/rqnm_paper6.mat', ...
                 {'../FJSIMS/Results/IMAC2022/Physics/epmc_imac2021_run9_pre.mat', ...
                            '../FJSIMS/Results/IMAC2022/Physics/epmc_imac2021_run9_iter16.mat '}...
                 {'../FJSIMS/Results/IMAC2022/Physics/epmc_imac2021_run8_pre.mat', ...
                            '../FJSIMS/Results/IMAC2022/Physics/epmc_imac2021_run8_iter12.mat ', ...
                            '../FJSIMS/Results/IMAC2022/Physics/epmc_imac2021_run8_v2_iter18.mat '}...
                 '../FJSIMS/Results/TMP/rqnm_run0_adhesion_.mat'
                };
             
        RQNMA_EPMC = [1, 1, 2, 2, 1, 1, 1];
        FIG_TITLE = 'Mesoscale';
        
        Rindex = 6;
        
        legend_names = {'Elastic,  $\mu$=0.03, Flat (RQNMA)', 'Elastic,  $\mu$=0.03, Measured Mesoscale (RQNMA)', ...
                        'Plastic, $\mu$=0.03, Measured Mesoscale (EPMC)', 'Plastic, $\mu$=CEB, Measured Mesoscale (EPMC)', ...
                        'Elastic $\mu$=Adhesion'};%,    'Edge Gap - 250 $\mu$m'};

        
        plot_markers = {'v--','h-', 'o:', '^--', 's-', 'd-', 'p:', '<--', '>:'};
        color_inds = [1:3, 5:10];
        marker_size = [6, 7, 7, 6, 7, 7, 7, 6, 6];
        
        line_width = 3;

%         numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [1e-1, 1e2];
%         Wlim = [155, 180];
%         Zlim = [2e-4 1e-1];
        
        numTrim = [0, 3]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [.75e-1, 2.1e1];
%         Wlim = [142, 184]; % If including 250 edge gap
        Wlim = [140, 184];
        Zlim = [2e-4 10e-2];
        
    case 7
        %%%%%%%%%%%%% PAPER: EPMC Physics
        NAMES = {'../FJSIMS/Results/PAPER/alphaRuns/rqnm_paper18.mat', ... % Sphere
                 {'../EPMC_SIMS/Results/PAPER/Run10/epmc_paper_run10_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run10/epmc_paper_run10_v2_iter.mat'},... % Elastic EPMC
                 {'../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_iter.mat'},... % 'Plastic 1, $\mu$=CEB'
                 {'../EPMC_SIMS/Results/PAPER/Run12/epmc_paper_run12_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run12/epmc_paper_run12_v2_iter.mat'},... % Plastic 2, $\mu$=CEB
                 {'../EPMC_SIMS/Results/PAPER/Run13/epmc_paper_run13_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run13/epmc_paper_run13_v2_iter.mat'},... % Plastic 1, $\mu$=0.02
                 {'../EPMC_SIMS/Results/PAPER/Run14/epmc_paper_run14_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run14/epmc_paper_run14_v2_iter.mat'},... % Plastic 1, $\mu$=0.03
                 {'../EPMC_SIMS/Results/PAPER/Run17/epmc_paper_run17_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run17/epmc_paper_run17_iter.mat'},... % Plastic 1, $\mu$=0.04
                 {'../EPMC_SIMS/Results/PAPER/Run15/epmc_paper_run15_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run15/epmc_paper_run15_iter.mat'},... % Plastic 1, Flat, $\mu$=0.03
                 {'../EPMC_SIMS/Results/PAPER/Run16/epmc_paper_run16_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run16/epmc_paper_run16_iter.mat'},... % Plastic 1, Flat, $\mu$=CEB
                 };
             
        RQNMA_EPMC = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
        FIG_TITLE = 'EPMC_Physics';
        
        Rindex = 6;
        
        legend_names = {'Elastic, RQNMA, $\mu$=0.03', ... % rqnm 18
                        'Elastic, EPMC, $\mu$=0.03', ... %10
                        'Plastic 1, $\mu$=CEB', ... %11
                        'Plastic 2, $\mu$=CEB', ... %12
                        'Plastic 1, $\mu$=0.02', ... %13
                        'Plastic 1, $\mu$=0.03', ... %14
                        'Plastic 1, $\mu$=0.04', ... %17
                        'Plastic 1, Flat, $\mu$=0.03', ... %15
                        'Plastic 1, Flat, $\mu$=CEB', ... %16
                        };


        plot_markers = {'s-', '^--','h-', 'o--', 'v-', '>--', 'd:', 'p-', '<--', '>:'};
        color_inds = [19, 13, 12, 15, 21, 22, 24:30];
        marker_size = [10, 7, 7, 6, 7, 7, 7, 6, 6];
        
        line_width = 3;

        %%%%%% FULL PLOT
        numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [1e-1, 1e2];
        Wlim = [154.5, 185];
        Zlim = [7e-4 1e-1];
        leg_down = .05;

                
        eliminateInds = {{[]}, {[12]}, {[]}, {[11]}, {[13:15]}, {[14]}, {[]}, {[]}, {[]}};
        
%         %%%% ZOOM IN PLOT
%         numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [1e-1, 1e2];
%         Wlim = [171.3, 184];
%         Zlim = [7e-4 1.1e-2];
%         SuppressLegend = 'off';
%         extra_info = '_Zoom';
        
    case 8
        %%%%%%%%%%%%% PAPER: EPMC Surface Settings
        NAMES = {{'../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_iter.mat'},...
                 {'../EPMC_SIMS/Results/PAPER/Run20/epmc_paper_run20_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run20/epmc_paper_run20_iter.mat'},...
                 {'../EPMC_SIMS/Results/PAPER/Run21/epmc_paper_run21_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run21/epmc_paper_run21_iter.mat'},...
                 {'../EPMC_SIMS/Results/PAPER/Run19/epmc_paper_run19_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run19/epmc_paper_run19_iter.mat'},...
                 };
             
        RQNMA_EPMC = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
        FIG_TITLE = 'EPMC_Physics';
        
        Rindex = 6;
        
        legend_names = {'Plastic 1, $\mu$=CEB, erode 2', ... %11
                        'Plastic 1, $\mu$=CEB, erode 25', ... %20
                        'Plastic 1, $\mu$=CEB, erode 50', ... %21
                        'Plastic 1, $\mu$=CEB, erode 75', ... %19
                        };%,    'Edge Gap - 250 $\mu$m'};

        
        eliminateInds = {{[]}, {[9:12, 19:22, 27:28]}, {[1:14]}, {[1:13]}};
        
        plot_markers = {'h-','d:', 'o-', 'v--', 's-', 'd-', 'p:', '<--', '>:'};
        color_inds = [12, 27, 30, 33:40];
        marker_size = [6, 7, 7, 6, 7, 7, 7, 6, 6];
        
        line_width = 3;

        numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [1e-1, 1e2];
        Wlim = [171.3, 181];
        Zlim = [7e-4 3e-2];
        
        leg_down = .15;
        
%         numTrim = [0, 3]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [.75e-1, 2.1e1];
%         Wlim = [140, 184];
%         Zlim = [2e-4 10e-2];
        
        
    case 9
        %%%%%%%%%%% IMAC 2022: Summary Figure
        NAMES = {{'../EPMC_SIMS/Results/PAPER/Run10/epmc_paper_run10_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run10/epmc_paper_run10_v2_iter.mat'},... % Elastic EPMC
                 {'../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_iter.mat'},... % 'Plastic 1, $\mu$=CEB'
                 {'../EPMC_SIMS/Results/PAPER/Run14/epmc_paper_run14_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run14/epmc_paper_run14_v2_iter.mat'},... % Plastic 1, $\mu$=0.03
                 {'../EPMC_SIMS/Results/PAPER/Run16/epmc_paper_run16_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run16/epmc_paper_run16_iter.mat'},... % Plastic 1, Flat, $\mu$=CEB
                 };
             
        RQNMA_EPMC = [2, 2,2, 2, 2, 2, 2, 2];
        FIG_TITLE = 'EPMC_Physics';
        
        Rindex = 6;
        
        legend_names = {'Elastic, $\mu$=0.03', ... %10
                        'Plastic 1, $\mu$=CEB', ... %11
                        'Plastic 1, $\mu$=0.03', ... %14
                        'Plastic 1, Flat, $\mu$=CEB', ... %16
                        };

        
        plot_markers = {'^--','h-', 'd:', '<--', '>:'};
        color_inds = [13, 12, 24 26:30];
        marker_size = [7, 7, 7, 7, 6, 6];
        
        line_width = 3;

        %%%%%% FULL PLOT
        numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [1e-1, 1e2];
        Wlim = [159.5, 182];
        Zlim = [7e-4 5e-2];
%         leg_down = .15;

                
        eliminateInds = {{[12]}, {[]}, {[14]}, {[]}};
        
    case 10
        %%%%%%%%%%% IMAC 2022: Physics Slides
        
        NAMES = {{'../EPMC_SIMS/Results/PAPER/Run10/epmc_paper_run10_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run10/epmc_paper_run10_v2_iter.mat'},... % Elastic EPMC
                 {'../EPMC_SIMS/Results/PAPER/Run13/epmc_paper_run13_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run13/epmc_paper_run13_v2_iter.mat'},... % Plastic 1, $\mu$=0.02
                 {'../EPMC_SIMS/Results/PAPER/Run14/epmc_paper_run14_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run14/epmc_paper_run14_v2_iter.mat'},... % Plastic 1, $\mu$=0.03
                 {'../EPMC_SIMS/Results/PAPER/Run17/epmc_paper_run17_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run17/epmc_paper_run17_iter.mat'},... % Plastic 1, $\mu$=0.04
                 {'../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_iter.mat'},... % 'Plastic 1, $\mu$=CEB'
                 {'../EPMC_SIMS/Results/PAPER/Run12/epmc_paper_run12_v2_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run12/epmc_paper_run12_v2_iter.mat'},... % Plastic 1, $\mu$=CEB
                 {'../EPMC_SIMS/Results/PAPER/Run15/epmc_paper_run15_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run15/epmc_paper_run15_iter.mat'},... % Plastic 1, Flat, $\mu$=0.03
                 {'../EPMC_SIMS/Results/PAPER/Run16/epmc_paper_run16_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Run16/epmc_paper_run16_iter.mat'},... % Plastic 1, Flat, $\mu$=CEB
                 };
             
        RQNMA_EPMC = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
        FIG_TITLE = 'EPMC_Physics';
        
        Rindex = 6;
        
        legend_names = {'Elastic, $\mu$=0.03', ... %10
                        'Plastic 1, $\mu$=0.02', ... %13
                        'Plastic 1, $\mu$=0.03', ... %14
                        'Plastic 1, $\mu$=0.04', ... %17
                        'Plastic 1, $\mu$=CEB', ... %11
                        'Plastic 2, $\mu$=CEB', ... %12
                        'Plastic 1, Flat, $\mu$=0.03', ... %15
                        'Plastic 1, Flat, $\mu$=CEB', ... %16
                        };

        
        plot_markers = {'^--', 'v-', '>--', 'd:','h-', 'o--', 'p-', '<--', '>:'};
        color_inds = [13, 21, 22, 24, 12, 15, 25:30];
        marker_size = [7, 7, 6, 7, 7, 7, 6, 6];
        
        line_width = 3;

        %%%%%% FULL PLOT
        numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [1e-1, 1e2];
        Wlim = [154.5, 185];
        Zlim = [7e-4 6e-2];
        leg_down = .05;
        
        max_plot = 8; % Run with 1,4,5,6,8
        
        extra_info = sprintf('_slide%u', max_plot);
        
        NAMES = NAMES(1:max_plot);
        legend_names = legend_names(1:max_plot);

                
        eliminateInds = {{[12]}, {[13:15]}, {[14]}, {[]}, {[]}, {[11]}, {[]}, {[]}, {[]}};
        
%         %%%% ZOOM IN PLOT
%         numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
%         Qlim = [1e-1, 1e2];
%         Wlim = [171.3, 184];
%         Zlim = [7e-4 1.1e-2];
%         SuppressLegend = 'off';
%         extra_info = '_Zoom';
        
    case 11
        % PAPER, Mode 2
        exp_mode = 2;
        
        %%%%%%%%%%%%% PAPER: EPMC Mode 2
        NAMES = {{'../EPMC_SIMS/Results/PAPER/Mode2/Run22/epmc_paper_run22_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Mode2/Run22/epmc_paper_run22_iter.mat'},... % Elastic EPMC
                 {'../EPMC_SIMS/Results/PAPER/Mode2/Run24/epmc_paper_run24_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Mode2/Run24/epmc_paper_run24_iter.mat'},... % 'Plastic 1, $\mu$=CEB'
                 {'../EPMC_SIMS/Results/PAPER/Mode2/Run23/epmc_paper_run23_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Mode2/Run23/epmc_paper_run23_iter.mat'},... % Plastic 1, $\mu$=0.03
                 {'../EPMC_SIMS/Results/PAPER/Mode2/Run25/epmc_paper_run25_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Mode2/Run25/epmc_paper_run25_iter.mat'},... % Plastic 1, Flat, $\mu$=0.03
                 {'../EPMC_SIMS/Results/PAPER/Mode2/Run26/epmc_paper_run26_pre.mat'...
                            '../EPMC_SIMS/Results/PAPER/Mode2/Run26/epmc_paper_run26_iter.mat'},... % Plastic 1, Flat, $\mu$=CEB
                 };
             
        RQNMA_EPMC = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
        FIG_TITLE = 'EPMC_Physics';
        
        Rindex = 6;
        
        legend_names = {'Elastic, $\mu$=0.03', ... %10 -> 22
                        'Plastic 1, $\mu$=CEB', ... %11 -> 24
                        'Plastic 1, $\mu$=0.03', ... %14 -> 23
                        'Plastic 1, Flat, $\mu$=0.03', ... %15 -> 25
                        'Plastic 1, $\mu$=CEB, erode 75', ... %19 -> 26
                        };

        
                       % 10,    11,     14,       15,      19
        plot_markers = {'^--',  'h-',   '>--',    'p-',    'v--',};
        color_inds =   [13,      12,    22,       25,      33];
        marker_size =  [7,       7,     7,        6,       6];
        
        line_width = 3;

        %%%%%% FULL PLOT
        numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [0.2, 1.2e3];
        Wlim = [562, 597];
        Zlim = [2e-4 0.017];
        leg_down = .05;

                
        eliminateInds = {{[]}, {[]}, {[]}, {[]}, {[]}, {[]}, {[]}, {[]}, {[]}};
        
    case 12
        %%%%%%%%% Just Plot the Experimental data, mode 2
        exp_mode = 2;
        NAMES = {};
        RQNMA_EPMC = [];
        
        Rindex = 6;
        
        legend_names = {};

        
        plot_markers = {'v-','h--', 'o:', '^-', 's--', 'd:', 'p-', '<--', '>:'};
        color_inds = [8, 15, 16, 7, 17, 18, 9, 19, 20];
        marker_size = [6, 6, 6, 6, 6, 6, 6, 6, 6];
        
        line_width = 3;

        %%%%%% FULL PLOT
        numTrim = [0, 0]; % Number of amplitudes to trim off [start, end] of rqnma
        Qlim = [0.6, 45];
        Wlim = [591.8, 595];
        Zlim = [2e-4 1e-2];
        
        show_z0 = false;
    otherwise
        error('Plot set not defined.');
end

color_plot = DISTINGUISHABLE_COLORS(max(color_inds), 'w');

%% Set Marker Sizes Consistently

markerSizesDictionary = {'v', 'h', 'o', '^', 's', 'd', 'p', '<', '>';
                          7,   7,   7.5,   7,   10,  7,   7,   7,   7};
            
for ii = 1:length(plot_markers)
    
    inds = cellfun(@(c)isequal(c, plot_markers{ii}(1)), markerSizesDictionary(1, :));
    
    marker_size(ii) = markerSizesDictionary{2, inds}(1);
    
end

%% Update Eliminate Indices if needed

if(isempty(eliminateInds))    
    eliminateInds = repmat({{[]}}, size(NAMES));
end

%% Load Data (Numerical Backbones)


% Viscous damping
z_visc = z_visc_opts(exp_mode);

BB = cell(size(NAMES));

for ii = 1:length(NAMES)
    if(RQNMA_EPMC(ii) == 1)
        [BB{ii}] = load_rqma(NAMES{ii}, Rindex, z_visc, numTrim);
    elseif(RQNMA_EPMC(ii) == 2)
        [BB{ii}] = load_epmc(NAMES{ii}, Rindex, z_visc, eliminateInds{ii});
    else
        error('Load type undefined');
    end
end


%% Load Experimental Data

if(exp_mode == 1)
    [BB_exp, BB_exp_low, BB_exp_high] = load_exp2('../EXPERIMENTAL_DATA/LTW_Mode1Trimmed_27Sept2021');
elseif(exp_mode == 2)
    [BB_exp, BB_exp_low, BB_exp_high] = load_exp2('../EXPERIMENTAL_DATA/LTW_Mode2Trimmed_27Sept2021');
else
    error('Not processed that experimental mode.');
end

legend_names = [{'Experimental Mean', 'Experimental Bounds'}, legend_names];


%% Plot Data

font_size = 14;

figure('Position', [200, 150, 900, 600]); % [left, bottom, width, height]
hold on;
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);


% %% Add Assumed Viscous Damping
subplot(2,1,2);
hold on;
if(show_z0)
    loglog(Qlim, [z_visc, z_visc], 'k--', 'LineWidth', line_width);
end
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlim(Qlim);
ylim(Zlim);
hold on;

% Plot Experiments
plot_bb(BB_exp, '-', 1, [0, 0, 0], line_width);
plot_bb(BB_exp_low, ':', 1, [0, 0, 0], line_width);
plot_bb(BB_exp_high, ':', 1, [0, 0, 0], line_width, 'off');

% %% Plot Numerical Results
for ii = 1:length(NAMES)
    plot_bb(BB{ii}, plot_markers{ii}, marker_size(ii), color_plot(color_inds(ii), :), line_width);
end

% %% Reset To Correct limits for W
subplot(2,1,1);
xlim(Qlim);
ylim(Wlim);
% xlabel('Amplitude [m/s^2]');
ylabel('Frequency [Hz]');

% %% Add Legend
subplot(2,1,2);
if(show_z0)
    legend_names = horzcat('Assumed Viscous Damping', legend_names);
end
leg = legend(legend_names, 'Location', 'nw');
xlabel('Amplitude [m/s$^2$]');
ylabel('Damping Factor');

% %% Size the Subplots as desired
ax1.XTickLabel = '';

height = 0.4;
bottom = 0.1;
gap = 0.03;
left = 0.1;
width = 0.75;

ax1.Position = [left, bottom+height+gap, width, height]; %[left, bottom, width, height]
ax2.Position = [left, bottom, width, height]; %[left, bottom, width, height]


set(ax1, 'FontSize', font_size);
set(ax2, 'FontSize', font_size);

legPos = leg.Position;
if(show_z0)
    leg.Position = [left+0.02, 0.5-legPos(4)/2-leg_down, legPos(3:end)];%[left, bottom, width, height]
end


ax1.Box = 'on';
ax2.Box = 'on';

leg.Visible = SuppressLegend;

set(gcf, 'Renderer', 'painters');
drawnow;

% print(sprintf('./OUTPUTS/Elastic_BB_Set%u%s.eps', plot_set, extra_info), '-depsc', '-r400');
% print(sprintf('./OUTPUTS/Elastic_BB_Set%u%s.svg', plot_set, extra_info), '-dsvg', '-r400');

% % For IMAC, Mesoscale Plot
% leg.Position = [0.12    0.25    0.2787    0.3157];


%% Functions for Plotting Stuff
function plot_bb(BB, marker, marker_size, color, line_width, handlevis)

    if(nargin < 6)
        handlevis = 'on';
    end
    
    subplot(2,1,1)
    h = semilogx(BB.accelQ, BB.W, marker, ...
                    'LineWidth', line_width, 'MarkerSize', marker_size, ...
                    'Color', color, 'MarkerFaceColor', color, ...
                    'HandleVisibility', handlevis); hold on
                
    set(h, 'MarkerFaceColor', get(h,'Color')); 
    
    
    
    subplot(2,1,2)
    h = loglog(BB.accelQ, BB.Z, marker, ...
                    'LineWidth', line_width, 'MarkerSize', marker_size, ...
                    'Color', color, 'MarkerFaceColor', color, ...
                    'HandleVisibility', handlevis); hold on
                
    set(h, 'MarkerFaceColor', get(h,'Color')); 

end

%% Load RQNMA Function
function [BB] = load_rqma(name, Rindex, z_visc, numTrim)

    tmp = load(name, 'BB', 'MESH', 'R', 'pars');
            
    % Displacement Amplitude [m]
    dispQ = abs(tmp.R(Rindex, :) * tmp.BB.Xavg(1:end-1, :));
    
    accelQ = dispQ.*(tmp.BB.W.^2);
    
    BB.dispQ = dispQ(1+numTrim(1):end-numTrim(2));
    BB.accelQ = accelQ(1+numTrim(1):end-numTrim(2));
    BB.W = tmp.BB.W(1+numTrim(1):end-numTrim(2))/2/pi;
    BB.Z = tmp.BB.Z(1+numTrim(1):end-numTrim(2))+z_visc;
    
    fprintf('The maximum RQNMA residual was %s\n', max(abs(tmp.BB.RHist), [], 'all'));
end

%% Load EPMC Function
function [BB] = load_epmc(name, Rindex, z_visc, eliminateInds)

    %%%%%% Assuming that names has two entries for now. If one finishes
    %%%%%% without needing pre/iter, then add a case here.

    tmp_pre = load(name{1}, 'MESH', 'R', 'pars');
    
    for ii = 2:length(name)
        tmp_bb = load(name{ii}, 'U');
        Ulc = tmp_bb.U;
        
        Ulc(:, eliminateInds{ii-1}) = [];

        if(ii > 2)
            Ulc = [Ulc_prev, Ulc];
        end
        
        W = Ulc(end-2, :);

        % Displacement Amplitude [m]
        maxind = size(Ulc, 2);
        Q1 = abs(tmp_pre.R(Rindex,:)*Ulc( (size(tmp_pre.R, 2)+1):(2*size(tmp_pre.R, 2)), 1:maxind)).*(10.^Ulc(end, 1:maxind));
        Q2 = abs(tmp_pre.R(Rindex,:)*Ulc( (2*size(tmp_pre.R, 2)+1):(3*size(tmp_pre.R, 2)), 1:maxind)).*(10.^Ulc(end, 1:maxind));
        dispQ = sqrt(Q1.^2 + Q2.^2);

        % Eliminate discarded indices here.
        
%         [~, maxind] = max(dispQ);
%         Ulc_prev = Ulc(:, 1:maxind);
    end
    
    accelQ = dispQ.*(W.^2);
    
    BB.dispQ = dispQ(1:maxind);
    BB.accelQ = accelQ(1:maxind);
    BB.W = W(1:maxind)/2/pi;
    BB.Z = Ulc(end-1, 1:maxind)/2./W(1:maxind);
    
end


%% Load Experimental Data
function [BB, BB_low, BB_high] = load_exp2(NAME)

    tmp = load(NAME, 'BB');
    
    BB.dispQ = tmp.BB.Qaccel./(2*pi*tmp.BB.Wmean).^2;
    BB.accelQ = tmp.BB.Qaccel;
    BB.W = tmp.BB.Wmean;
    BB.Z = tmp.BB.Zmean;
    
    BB_low.accelQ = BB.accelQ;
    BB_low.W = tmp.BB.Wmin;
    BB_low.Z = tmp.BB.Zmin;
% 
    BB_high.accelQ = BB.accelQ;
    BB_high.W = tmp.BB.Wmax;
    BB_high.Z = tmp.BB.Zmax;
end



