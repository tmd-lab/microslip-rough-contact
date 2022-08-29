function inputpars = rqnm_pars(pars_num)
% Function for defining parameter sets used in RQNM simulations

switch pars_num
    case 0
        warning('Using default parameters - random for testing.');
        
        inputpars.mu = 0.05;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        
        % Aperity Function
%         inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE; inputpars.Nqp_radius = 1; 
        inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE_ADHESION; inputpars.Nqp_radius = 1; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_COUPLE; inputpars.Nqp_radius = 100; 
        
        inputpars.alpha = 0; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        inputpars.Nqp_heights = 100;
        
        inputpars.meshName = 'zte152'; % 'zte152', 'zte588', 'zte588Jul21', 'zte588convgX' (X = 7, 8, 9)
        inputpars.meshNameRelative = ''; % '' or '_NOREL'

        inputpars.output_name = 'rqnm_run0_adhesion_';
        inputpars.load_initial = false;
        
        % RQNMA Solver Settings
        inputpars.repeatLoop = 4;
        inputpars.Nhp = 9;
        
        inputpars.DeltaGamma = 5; % J/m^2
        
    case 1 
        inputpars.mu = 0.02;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        
        % Aperity Function
        inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE; inputpars.Nqp_radius = 1; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_COUPLE; inputpars.Nqp_radius = 100; 
        
        inputpars.alpha = 0; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = true; % Convert ellipsoid properties to sphere
        inputpars.Nqp_heights = 100;
        
        inputpars.meshName = 'zte588convg10'; % 'zte152', 'zte588', 'zte588Jul21', 'zte588convgX' (X = 7, 8, 9, 10), 'UROM122', 'UROM232'
        inputpars.meshNameRelative = '_NOREL'; % '' or '_NOREL'

        inputpars.output_name = 'Results/rqnm_paper1';
        inputpars.load_initial = false;
        
        % RQNMA Solver Settings
        inputpars.repeatLoop = 4;
        inputpars.Nhp = 9;
        
    case 2 
        inputpars.mu = 0.02;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        
        % Aperity Function
        inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE; inputpars.Nqp_radius = 1; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_COUPLE; inputpars.Nqp_radius = 100; 
        
        inputpars.alpha = 0; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = true; % Convert ellipsoid properties to sphere
        inputpars.Nqp_heights = 100;
        
        inputpars.meshName = 'UROM122'; % 'zte152', 'zte588', 'zte588Jul21', 'zte588convgX' (X = 7, 8, 9, 10), 'UROM122', 'UROM232'
        inputpars.meshNameRelative = ''; % '' or '_NOREL'

        inputpars.output_name = 'Results/rqnm_paper2';
        inputpars.load_initial = false;
        
        % RQNMA Solver Settings
        inputpars.repeatLoop = 4;
        inputpars.Nhp = 9;
         
    case 3 
        inputpars.mu = 0.02;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        
        % Aperity Function
        inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE; inputpars.Nqp_radius = 1; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_COUPLE; inputpars.Nqp_radius = 100; 
        
        inputpars.alpha = 0; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = true; % Convert ellipsoid properties to sphere
        inputpars.Nqp_heights = 100;
        
        inputpars.meshName = 'UROM232'; % 'zte152', 'zte588', 'zte588Jul21', 'zte588convgX' (X = 7, 8, 9, 10), 'UROM122', 'UROM232'
        inputpars.meshNameRelative = ''; % '' or '_NOREL'

        inputpars.output_name = 'Results/rqnm_paper3';
        inputpars.load_initial = false;
        
        % RQNMA Solver Settings
        inputpars.repeatLoop = 4;
        inputpars.Nhp = 9;
    
    case 4 
        inputpars.mu = 0.01;
        inputpars.ElasticModulus = 192.31e9; % Pa
        inputpars.nu = 0.3;
        
        inputpars.mesoscale = true; % Multiplies ZTE heights;
        
        % Aperity Function
        inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE; inputpars.Nqp_radius = 1; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_DECOUPLE; inputpars.Nqp_radius = 100; 
%         inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_COUPLE; inputpars.Nqp_radius = 100; 
        
        inputpars.alpha = pi/2; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        inputpars.Nqp_heights = 100;
        
        inputpars.meshName = 'UROM232'; % 'zte152', 'zte588', 'zte588Jul21', 'zte588convgX' (X = 7, 8, 9, 10), 'UROM122', 'UROM232'
        inputpars.meshNameRelative = ''; % '' or '_NOREL'

        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
        inputpars.load_initial = false;
        
        % RQNMA Solver Settings
        inputpars.repeatLoop = 4;
        inputpars.Nhp = 9;
    
    case 5 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.02;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    case 6 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    case 7 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.05;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    case 8 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.1;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    case 9 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.5;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    case 10 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.9;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
        
    
    case 11 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        
        inputpars.asp_fun = @ELLIPSOID_IWAN_DECOUPLE; inputpars.Nqp_radius = 100; 
        inputpars.alpha = pi/2; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    
    case 12 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        
        inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_DECOUPLE; inputpars.Nqp_radius = 100; 
        inputpars.alpha = pi/2; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    
    case 13 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        
        inputpars.asp_fun = @ELLIPSOID_IWAN_FIT_COUPLE; inputpars.Nqp_radius = 100; 
        inputpars.alpha = pi/2; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    
    case 14 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        
        inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE; inputpars.Nqp_radius = 1; 
        inputpars.alpha = 0; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = true; % Convert ellipsoid properties to sphere
        
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    
    case 15 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        
        inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE; inputpars.Nqp_radius = 1; 
        inputpars.alpha = 0; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    
    case 16 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        
        inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE; inputpars.Nqp_radius = 1; 
        inputpars.alpha = pi/4; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    
    case 17 
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        
        inputpars.asp_fun = @ELLIPSOID_TAN_DECOUPLE; inputpars.Nqp_radius = 1; 
        inputpars.alpha = -pi/4; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    
    case 18 
        inputpars = rqnm_pars(13);
        
        inputpars.alpha = 0; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = true; % Convert ellipsoid properties to sphere
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
       
    case 19
        inputpars = rqnm_pars(13);
        
        inputpars.alpha = 0; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
       
    case 20
        inputpars = rqnm_pars(13);
        
        inputpars.alpha = pi/4; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
       
    case 21
        inputpars = rqnm_pars(13);
        
        inputpars.alpha = -pi/4; % Rotation from x to small ellipsoid axis
        inputpars.useSphere = false; % Convert ellipsoid properties to sphere
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
    
    case 22
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
        
        %%%% Center Gap
        l = 0.121412; % length of interface
        h = 50e-6; % height of shape
        R = h/2 + l^2/8/h;
        inputpars.overrideMesoscale = @(XY) sqrt(R^2 - XY(:, 1).^2);
    
    case 23
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
        
        %%%% Center Gap
        l = 0.121412; % length of interface
        h = 100e-6; % height of shape
        R = h/2 + l^2/8/h;
        inputpars.overrideMesoscale = @(XY) sqrt(R^2 - XY(:, 1).^2);
        
        
    case 24
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
        
        %%%% Center Gap
        l = 0.121412; % length of interface
        h = 250e-6; % height of shape
        R = h/2 + l^2/8/h;
        inputpars.overrideMesoscale = @(XY) sqrt(R^2 - XY(:, 1).^2);
        
        
    case 25
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
        
        %%%% Center Gap
        l = 0.121412; % length of interface
        h = 50e-6; % height of shape
        R = h/2 + l^2/8/h;
        inputpars.overrideMesoscale = @(XY) -sqrt(R^2 - XY(:, 1).^2);
        
        
    case 26
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
        
        %%%% Center Gap
        l = 0.121412; % length of interface
        h = 250e-6; % height of shape
        R = h/2 + l^2/8/h;
        inputpars.overrideMesoscale = @(XY) -sqrt(R^2 - XY(:, 1).^2);
       
        
    case 27
        % Flat interface for comparison.
        inputpars = rqnm_pars(4);
        inputpars.mu = 0.03;
        inputpars.output_name = sprintf('Results/rqnm_paper%u', pars_num);
        
        inputpars.mesoscale = false; % Multiplies ZTE heights;
        
        
    otherwise
        error('Did not supply valid parameter input');
end



end
