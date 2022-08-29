function inputpars = epmc_pars(pars_num)
% Function for defining parameter sets used in EPMC simulations

inputpars.mode = 1; % Default First Mode

switch pars_num
    case 0
        warning('Using default parameters - random for testing.');
        
        inputpars.Sys = 331.7e6; %330e6; %Pa
        inputpars.mu = 0.9;
        inputpars.ElasticModulus = 192.31e9; % 192.31e9 % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %mu + CEB
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake'; %'brake', 'etsion', 'elastic_ellip_tan';
        inputpars.Et = 620e6;
        inputpars.useSphere = true;
        
%         inputpars.alpha = pi/2;
%         inputpars.area_density = 5.485e6;
%         inputpars.Rprime  = 5.46e-3/2;
%         inputpars.Rpprime = 5.46e-3/2;
   
        inputpars.meshName = 'zte152'; % 'zte152', 'zte588', 'zte588Jul21', 'zte588convgX' (X = 7, 8, 9, 10), 'UROM122', 'UROM232'
        inputpars.meshNameRelative = ''; % '' or '_NOREL'
        
        inputpars.output_name = 'epmc_delete';
        inputpars.load_initial = false;
%         inputpars.prevInd = 
%         inputpars.loadName = 

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.2;
        inputpars.arcSettings.dsmin = 0.001;

        
    case 1
        
        % Full / Best Model, TMD Abstract
        
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %mu + CEB
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'zte152'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = 'epmc_tmd_absV2_run1';
        inputpars.load_initial = false;
%         inputpars.prevInd = 
%         inputpars.loadName = 'Results/

        inputpars.arcSettings.da = 1e3;
        inputpars.arcSettings.dsmax = 80e5;
        inputpars.arcSettings.dsmin = 1e3;

        
    case 2
        
        % Without Mesoscale Model, TMD Abstract
        
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %mu + CEB
        inputpars.mesoscale = false; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'zte152'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = 'epmc_tmd_absV2_run2';
        inputpars.load_initial = false;
%         inputpars.prevInd = 
%         inputpars.loadName = 'Results/

        inputpars.arcSettings.da = .5e3;
        inputpars.arcSettings.dsmax = 80e5;
        inputpars.arcSettings.dsmin = 1e2;

        
    case 3
        
        % Etsion Unloading, TMD Abstract
        
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %mu + CEB
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'ETSION';
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'zte152'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = 'epmc_tmd_absV2_run3';
        inputpars.load_initial = false;
%         inputpars.prevInd = 
%         inputpars.loadName = 'Results/

        inputpars.arcSettings.da = 1e3;
        inputpars.arcSettings.dsmax = 80e5;
        inputpars.arcSettings.dsmin = 1e3;

        
    case 4
        
        % Elastic Model, TMD Abstract
        
        inputpars.Sys = Inf; %Pa
        inputpars.mu = 0.02;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 1; %mu 
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'zte152'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = 'epmc_tmd_absV2_run4d';
        inputpars.load_initial = true;
        inputpars.prevInd = 8;
        inputpars.loadName = 'Results/epmc_tmd_absV2_run4c_iter8.mat';

        inputpars.arcSettings.da = 4e6;
        inputpars.arcSettings.dsmax = 5e7;
        inputpars.arcSettings.dsmin = 5e5;

        
        
        
    case 5
        
        % Full + mu, TMD Abstract
        
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 0.02;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %mu + CEB
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'zte152'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = 'epmc_tmd_absV2_run5d';
        inputpars.load_initial = true;
        inputpars.prevInd = 78;
        inputpars.loadName = 'Results/epmc_tmd_absV2_run5_iter78.mat';

        inputpars.arcSettings.da = 5e6;
        inputpars.arcSettings.dsmax = 5e7;
        inputpars.arcSettings.dsmin = 5e5;

        
    case 6
        
        % Elastic + mu, Mesoscale = flat TMD Presentation
        
        inputpars.Sys = Inf; %Pa
        inputpars.mu = 0.02;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 1; %mu + CEB
        inputpars.mesoscale = false; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'zte152'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = 'epmc_tmd_absV2_run6d';
        inputpars.load_initial = false;
        inputpars.prevInd = 78;
        inputpars.loadName = 'Results/epmc_tmd_absV2_run6d_iter*.mat';

        inputpars.arcSettings.da = 5e6;
        inputpars.arcSettings.dsmax = 5e7;
        inputpars.arcSettings.dsmin = 5e5;
        
    case 7
        
        % Elastic Model, Comparison to RQNMA, Paper
        
        inputpars.Sys = Inf; %Pa
        inputpars.mu = 0.02;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 1; %mu 
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.Et = 620e6;
        
        inputpars.useSphere = true;

        inputpars.meshName = 'UROM232'; % 'zte152', 'zte588', 'zte588Jul21', 'zte588convgX' (X = 7, 8, 9, 10), 'UROM122', 'UROM232'
        inputpars.meshNameRelative = ''; % '' or '_NOREL'
        
        inputpars.output_name = 'epmc_paper7';
        inputpars.load_initial = false;
%         inputpars.prevInd = 8;
%         inputpars.loadName = 'Results/epmc_paper7.mat';

        inputpars.arcSettings.da = 4e6;
        inputpars.arcSettings.dsmax = 5e7;
        inputpars.arcSettings.dsmin = 5e5;

    case 8
        
        % Full / Best Model, IMAC 2021 Abstract
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %mu + CEB
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_imac2021_run%u_v3', pars_num);
        inputpars.load_initial = true;
        inputpars.prevInd = 18;
        inputpars.loadName = 'Results/epmc_imac2021_run8_v2_iter18';

        inputpars.arcSettings.da = 4e5;
        inputpars.arcSettings.dsmax = 64e6;
        inputpars.arcSettings.dsmin = 1e3;

    case 9
        
        % Just mu, IMAC 2021 Abstract
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 0.03;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 1; %mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_imac2021_run%u_v3', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = 
%         inputpars.loadName = 'Results/

        inputpars.arcSettings.da = 4e7;
        inputpars.arcSettings.dsmax = 4e8;
        inputpars.arcSettings.dsmin = 1e3;
    
    case 10
        % Paper: Backbone 1 - Fully Elastic
        inputpars.Sys = Inf; %Pa
        inputpars.mu = 0.03;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 1; %mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 192.31e9;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_paper_run%u_v2', pars_num);
        inputpars.load_initial = true;
        inputpars.prevInd = 15;
        inputpars.loadName = 'Results/epmc_paper_run10_v2_iter';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.1;
        inputpars.arcSettings.dsmin = 0.001;
    
    case 11
        % Paper: Backbone 2 - SNL Properties
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %CEB+mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_paper_run%u_v2', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = 
%         inputpars.loadName = 'Results/

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.2;
        inputpars.arcSettings.dsmin = 0.001;
    
    case 12
        % Paper: Backbone 3 - Qu et al. properties
        inputpars.Sys = 150e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %CEB+mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 3.5e9; % Pa
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_paper_run%u_v2', pars_num);
        inputpars.load_initial = true;
        inputpars.prevInd = 10;
        inputpars.loadName = 'Results/epmc_paper_run12_v2_iter.mat';

        inputpars.arcSettings.da = 0.1; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.2;
        inputpars.arcSettings.dsmin = 0.001;
    
    
    case 13
        % Paper: Backbone 4 - SNL Properties, low mu
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 0.02;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 1; %mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_paper_run%u_v2', pars_num);
        inputpars.load_initial = true;
        inputpars.prevInd = 13;
        inputpars.loadName = 'Results/epmc_paper_run13_v2_iter.mat';
        
        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.1;
        inputpars.arcSettings.dsmin = 0.001;
    
    
    case 14
        % Paper: Backbone 5 - SNL Properties, upper mu
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 0.03;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 1; %mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_paper_run%u_v2', pars_num);
        inputpars.load_initial = true;
        inputpars.prevInd = 14;
        inputpars.loadName = 'Results/epmc_paper_run14_v2_iter.mat';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.2;
        inputpars.arcSettings.dsmin = 0.001;
    
    
    case 15
        % Paper: Backbone 6 - SNL Properties, w/o Mesoscale + mu=0.03
        
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 0.03;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 1; %mu only
        inputpars.mesoscale = false; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = 
%         inputpars.loadName = 'Results/

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.1;
        inputpars.arcSettings.dsmin = 0.001;
    
    
    case 16
        % Paper: Backbone 7 - SNL Properties, w/o Mesoscale + mu=CEB
        
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %CEB
        inputpars.mesoscale = false; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = 
%         inputpars.loadName = 'Results/

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.1;
        inputpars.arcSettings.dsmin = 0.001;
    
    
    case 17
        % Paper: Backbone 8 - SNL Properties, upper mu=0.04
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 0.04;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 1; %mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = 25;
%         inputpars.loadName = 'Results/epmc_paper_run17_iter.mat';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.1;
        inputpars.arcSettings.dsmin = 0.001;
    
    case 18
        % Paper: Backbone 9 - SNL Properties
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 5; %CEB+mu only w/ UTS CEB
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        inputpars.UTS = 713.5e6; % Pa
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = 
%         inputpars.loadName = 'Results/

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.2;
        inputpars.arcSettings.dsmin = 0.001;
    
        
    
    case 19
        % Paper: Backbone 10 - SNL Properties + erroded asperities
        
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %CEB+mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.surfaceName = 'combined_2dec21_R1_erode75.mat'; % Do not include relative path to SURFACE/OUT
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
        inputpars.prevInd = 9;
        inputpars.loadName = 'Results/epmc_paper_run19_iter';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.12;
        inputpars.arcSettings.dsmin = 0.001;

    
    
    case 20
        % Paper: Backbone 10 - SNL Properties + erroded asperities
        
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %CEB+mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.surfaceName = 'combined_2dec21_R1_erode25.mat'; % Do not include relative path to SURFACE/OUT
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = true;
        inputpars.prevInd = 26;
        inputpars.loadName = 'Results/epmc_paper_run20_iter';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.12;
        inputpars.arcSettings.dsmin = 0.001;

    
    
    case 21
        % Paper: Backbone 10 - SNL Properties + erroded asperities
        
        inputpars.Sys = 330e6; %Pa
        inputpars.mu = 1e10;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        inputpars.sliptype = 3; %CEB+mu only
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        inputpars.unloadModel = 'brake';
        inputpars.useSphere = true;
        inputpars.Et = 620e6;
        
        inputpars.meshName = 'UROM232'; %'zte152', 'zte588', 'zte588Jul21'
        
        inputpars.surfaceName = 'combined_2dec21_R1_erode50.mat'; % Do not include relative path to SURFACE/OUT
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
        inputpars.prevInd = 9;
        inputpars.loadName = 'Results/epmc_paper_run21_iter';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.12;
        inputpars.arcSettings.dsmin = 0.001;

            
    case 22
        
        % Elastic Spheres, Mode 2
        inputpars = epmc_pars(10);
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = X;
%         inputpars.loadName = 'Results/epmc_paper_runXX_iter';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.12;
        inputpars.arcSettings.dsmin = 0.001;
        
        inputpars.mode = 2; % Second Bending Mode
   
    case 23
        
        % Elastic Spheres, Mode 2
        inputpars = epmc_pars(14);
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = X;
%         inputpars.loadName = 'Results/epmc_paper_runXX_iter';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.12;
        inputpars.arcSettings.dsmin = 0.001;
        
        inputpars.mode = 2; % Second Bending Mode
        
        
    case 24
        
        % Elastic Spheres, Mode 2
        inputpars = epmc_pars(11);
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = X;
%         inputpars.loadName = 'Results/epmc_paper_runXX_iter';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.12;
        inputpars.arcSettings.dsmin = 0.001;
        
        inputpars.mode = 2; % Second Bending Mode
        
    case 25
        
        % Elastic Spheres, Mode 2
        inputpars = epmc_pars(15);
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = X;
%         inputpars.loadName = 'Results/epmc_paper_runXX_iter';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.12;
        inputpars.arcSettings.dsmin = 0.001;
        
        inputpars.mode = 2; % Second Bending Mode
        
    case 26
        
        % Elastic Spheres, Mode 2
        inputpars = epmc_pars(19);
        
        inputpars.output_name = sprintf('epmc_paper_run%u', pars_num);
        inputpars.load_initial = false;
%         inputpars.prevInd = X;
%         inputpars.loadName = 'Results/epmc_paper_runXX_iter';

        inputpars.arcSettings.da = 0.04; % If displacements are near constant, then steps are 0.02 to start (since dA is scaled by b=0.5).
        inputpars.arcSettings.dsmax = 0.12;
        inputpars.arcSettings.dsmin = 0.001;
        
        inputpars.mode = 2; % Second Bending Mode
        
        inputpars.ITMAX = 6;
        
    otherwise
        error('Did not supply valid parameter input');
end



end
