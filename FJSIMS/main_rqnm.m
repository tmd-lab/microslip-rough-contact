function [BB] = main_rqnm(varargin)

nots_min_cores = 10;

if(nargin == 0)
    warning('Defaulting to parameter set 0.');
    inputpars = rqnm_pars(0);
    runnum = 0;
else
    fprintf('Loading parameters number %u. \n', varargin{1});
    inputpars = rqnm_pars(varargin{1});
    runnum = varargin{1};
end

if(nargin > 1)
    nots_min_cores = varargin{2};
end

% %% Configure parallel and NOTS

cores_avail = feature('numcores');

if(cores_avail < nots_min_cores)
    fprintf('There are %u cores, running on Local.\n', cores_avail);
    
    p = gcp('nocreate');
    
    if isempty(p)
        myCluster = parpool('local');
    end
    
    onnots = false;

else
    fprintf('There are %u cores, running on NOTS.\n', cores_avail);

    myCluster = parpool('threads');

    %set working directory:
    workdir = [getenv('SHARED_SCRATCH') filesep getenv('USER') filesep 'microslip-rough-contact' filesep 'FJSIMS'];
    cd(workdir);
    
    mkdir('Results'); % Issues a warning if exists, prevents fatal save error.
    
    onnots = true;
end


addpath('../ROUTINES/FRIC_MODELS')
addpath('../ROUTINES/FRIC_MODELS/ASP_FUN/')
% addpath('../ROUTINES/FRIC_MODELS/PLASTIC/')
addpath('../ROUTINES/FRIC_MODELS/MINDLIN/')
addpath('../ROUTINES/RQNMA')
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')


output_name = inputpars.output_name;

%% Load ROM

if(isequal(inputpars.meshName, 'zte152'))
    % 152 ZTE Model
    load('../FJSIMS/ROMS/ROM_PD_152ELS.mat', 'M', 'K', 'R', 'Fv', 'L', 'TFMfhcb', 'MESH');
    
elseif(isequal(inputpars.meshName, 'zte588'))
    % Full order (588 ZTE) model.
    load('../FJSIMS/ROMS/5_SET_NULLRED.mat', 'M', 'K', 'R', 'Fv', 'L', 'MESH');
end

if( isequal(inputpars.meshName, 'zte588') || isequal(inputpars.meshName, 'zte152') )
    % Rescale density and stiffness for the new beams
    Ecurr = inputpars.ElasticModulus; %Pa
    M = M*(7784.3/7857.8); %Densities in kg/m^3
    K = K*(Ecurr / 1.9285e11); % Elastic modulus in Pa. 

    % Accelerometer Mass
    maccel = 0.8e-3; % 0.8 grams
    M = M + (R'*R)*maccel;
    
end

%% Load MESH

if( isequal(inputpars.meshName, 'zte588Jul21') )
    
    load('../FJSIMS/ROMS/MATS_30Jun21/MATRICES_NR.mat', 'M', 'K', 'R', 'Fv', 'L');

    Nds = dlmread('../FJSIMS/ROMS/MATS_30Jun21/Nodes.dat');
    Quad = dlmread('../FJSIMS/ROMS/MATS_30Jun21/Elements.dat');

    Nq = 1;
    MESH = MESH2D(Nds, 3, [], Quad, Nq);

    % Accelerometer Mass
    maccel = 0.8e-3; % 0.8 grams
    M = M + (R(1:6, :)'*R(1:6, :))*maccel;

end


%% MESH Convergence Checking

if( length(inputpars.meshName) >= 11 && isequal(inputpars.meshName(1:11), 'zte588convg') )
    
    load(sprintf('../FJSIMS/ROMS/CONVERGENCE/%s_SET_NULLRED%s', inputpars.meshName(12:end), inputpars.meshNameRelative), 'M', 'K', 'R', 'Fv', 'L');

    Nds = dlmread('../FJSIMS/ROMS/MATS_30Jun21/Nodes.dat');
    Quad = dlmread('../FJSIMS/ROMS/MATS_30Jun21/Elements.dat');

    Nq = 1;
    MESH = MESH2D(Nds, 3, [], Quad, Nq);

    % Accelerometer Mass
    maccel = 0.8e-3; % 0.8 grams
    M = M + (R(4:9, :)'*R(4:9, :))*maccel;
    
    fprintf('There are %u Fixed Interface Modes.\n', mod(size(M, 1), 3*MESH.Nn));

end

%% Uniform ROM

if( isequal(inputpars.meshName(1:4), 'UROM') )
    
    load(sprintf('../FJSIMS/ROMS/ROM_U_%sELS', inputpars.meshName(5:end)), 'M', 'K', 'R', 'Fv', 'L', 'MESH');

    % Accelerometer Mass
    maccel = 0.8e-3; % 0.8 grams
    M = M + (R(4:9, :)'*R(4:9, :))*maccel;
    
    fprintf('There are %u Fixed Interface Modes.\n', mod(size(M, 1), 3*MESH.Nn));

end



%% Rename Prestress Force Vector

% This is just for the name making sense to prevent several bugs where the
% naming has confused me.
Fvnorm = Fv;


%% Model Parameters

Prestress = (12002+12075+12670)/3;%12845;
mds = [1 3 5];
mode_ind = 1; %set to 1,2,3 for first three bending modes
mi = mds(mode_ind);  % Mode of interest


%% Mapping Matrices
MESH.dpn = 3;
[Q1, T1] = ZTE_ND2QP(MESH, 1);

% %% ZTE Quadrature Matrices
No   = 1;                  % Number of GLQ integration points in each direction
[Qm,Tm] = ZTE_ND2QP(MESH,No);
QuadMats.Q = Qm;
QuadMats.T = Tm;%*Qm*inv(Qm'*Qm)*Qm';
SCP = (Qm'*Qm)\Qm';
 
assert(abs(sum(Tm(:)) - 0.002933013353519)/0.002933013353519 < 0.0051, 'Area of interface has changed');
%WJROM Patch Areas - 0.002918218473670
%588 DOF - 0.002918218473670
%152 DOF - 0.002933013353519

%% Verify the Recovery Nodes before prestress

ktkn = 1e10;

Kstuck = zeros(size(L,1));
Kstuck(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = Tm*ktkn*Qm*0; %Frictinless prestress
Kstuck(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = Tm*ktkn*Qm*0; %frictionless prestress
Kstuck(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = Tm*ktkn*Qm;

K0 = K+L'*Kstuck*L;

[Vst, Wst] = eigs(K0, M, 10, 'SM');
[Wst, si] = sort(sqrt(diag(Wst)));
Vst = Vst(:, si);  Vst = Vst./sqrt(diag(Vst'*M*Vst)');

% Order is [X1 Y1 Z1, X2 Y2 Z2, . . . ];
responseXYZ = R*Vst(:, 1);

responseXYZ(3:3:end)



%% Physics Parameters

%Initialize the sphere parameters
E = inputpars.ElasticModulus; %304L from current (e.g., Nidish updating for Long Term wear).
nu = inputpars.nu; %304L taken from Traction paper
G = E / 2/ (1 + nu);%304L approx Shear modulus

% pars.mu     = 0.075; % set with friction model, what is reasonable changes
pars.E      = E;
pars.Estar  = E/2/(1 - nu^2);
pars.nu     = nu;
pars.G      = G;
pars.Gstar  = G/2/(2 - nu);

pars.mu = inputpars.mu;

if(isfield(inputpars, 'DeltaGamma'))
    pars.DeltaGamma = inputpars.DeltaGamma;
end

%% Set Discretization + Contact Function

% Number of Quadrature Points
Nqp_radius = inputpars.Nqp_radius; %overridden to 1 if using a tangent or secant asperity model below.
Nqp_heights = inputpars.Nqp_heights;


ASP_FUN = inputpars.asp_fun;
ASP_FUN_PRE = @ELLIPSOID_PRE;
ELEM_TRAC = @RC_TRACTION;

if(isfield(inputpars, 'alpha') && inputpars.alpha ~=0)
        
    UN_ROT_FUN = ASP_FUN;
    
    ASP_FUN = @(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq)ASP_ROTATION(UN_ROT_FUN, pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq);
    
    pars.alpha = inputpars.alpha;
end

%% Processed surface load

matpars = pars;

%%%% Load Surface Processed Info

%%%%% Set 8 - Initial Paper
surf = load('../SURFACE/OUT/combined_14sep21_R1.mat'); 
combined_surf = surf.combinedSurf;
Qps = QuadMats.Q*MESH.Nds;

% Calculate gaps based on meshing the interface at nodes
if(isfield(inputpars, 'overrideMesoscale'))
    nd_gaps = inputpars.overrideMesoscale(MESH.Nds);
else
    nd_gaps = -surf.meso_surf{1}(MESH.Nds) - surf.meso_surf{2}(MESH.Nds);
end

nd_gaps = nd_gaps - min(nd_gaps);
zte_gaps = inputpars.mesoscale*QuadMats.Q*nd_gaps;


if(~onnots)
    
    figure;
    scatter3(Qps(:, 1), Qps(:, 2), zte_gaps, 1000, zte_gaps, '.'); %9 is the point size
    hold on;
    xlabel('x');
    ylabel('y');
    zlabel('z');

    colormap(jet); colorbar;

    aspectRatios = [range(MESH.Nds), 0.5*range(MESH.Nds(:, 2))];
    xlim([min(MESH.Nds(:, 1)), max(MESH.Nds(:, 1))]);
    ylim([min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))]);

    pbaspect(aspectRatios);
    
    if( isequal(inputpars.meshName, 'zte588Jul21') )

        figure('Position', [600 400 1200 400]);
        MESH.SHOWFIELD3D(nd_gaps)
        c = colorbar;
        c.Label.String = 'Gap [m]';

        xrange = [min(MESH.Nds(:, 1)), max(MESH.Nds(:, 1))];
        yrange = [min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))];
        zrange = [0, 1.5e-4];
        
        xlim(xrange);
        ylim(yrange);
        zlim(zrange);
        view([-5.5448   36.0170]);
        
        pbaspect([range(xrange), range(yrange), 0.5*range(yrange)])
        
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('Gap [m]');
        
        set(gca, 'fontsize', 12);
    
        
        set(gcf, 'Renderer', 'painters');

        % print('./tmd_ltw_mesoscale.svg', '-dsvg', '-r400');
    end
    
    
end

%%%%% End Load Surface Process

%Update parameters with the surface scan.
[pars, PZFUN, CFUN, CFUN_PRE, area_density, zmax, zmin] = add_surf_pars(matpars, ...
    combined_surf, ASP_FUN, ASP_FUN_PRE, ELEM_TRAC, Nqp_heights, ...
    Nqp_radius, zte_gaps, inputpars.useSphere);

fprintf('R''=%.3s, R"=%.3s, Re=%.3s, zmax=%.3s \n', pars.Rprime, pars.Rpprime, pars.Re, zmax-zmin)

clear surf combined_surf

%% Stuck interface initial guess

% uxyn_est = [0, 0, abs(zmax)/2/4/2];
% uxyn_est = [0, 0, 4e-5]; %45 iterations
uxyn_est = [0, 0, .1e-5];

%initialize prev.
clear prev;

prev.tx0 = zeros(Nqp_heights, Nqp_radius); %traction for each slider of each asperity height - X
prev.ty0 = zeros(Nqp_heights, Nqp_radius); %traction for each slider of each asperity height - X
prev.rq0 = ones(Nqp_heights, 1)*linspace(0, 1, Nqp_radius); %radii for each slip radius used
prev.uxyw0 = ones(Nqp_heights, 1)*[0, 0, 0]; %displacement at previous step. 
prev.deltam = zeros(Nqp_heights, 1); %Max previous normal displacement
prev.Fm = zeros(Nqp_heights, 1);%Max previous normal force
prev.am = zeros(Nqp_heights, 1);%reversal contact radius


[txyn, dtxynduxyn, prev] = ELEM_TRAC(pars, uxyn_est, prev, ASP_FUN, PZFUN, Nqp_heights, Nqp_radius, zmin, zmax, area_density);


Kstuck = zeros(size(L,1));
Kstuck(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = Tm*dtxynduxyn(1,1)*Qm*0; %Frictinless prestress
Kstuck(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = Tm*dtxynduxyn(2,2)*Qm*0; %frictionless prestress
Kstuck(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = Tm*dtxynduxyn(3,3)*Qm;

K0 = K+L'*Kstuck*L;
X0 = K0\(Fvnorm*Prestress);




%% Simple Prestress

opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true, 'MaxIterations', 100);

% History Initialize
prevS1.tx0 = zeros(Nqp_heights, Nqp_radius); %traction for each slider of each asperity height - X
prevS1.ty0 = zeros(Nqp_heights, Nqp_radius); %traction for each slider of each asperity height - X
prevS1.rq0 = ones(Nqp_heights, 1)*linspace(0, 1, Nqp_radius); %radii for each slip radius used
prevS1.uxyw0 = ones(Nqp_heights, 1)*[0, 0, 0]; %displacement at previous step. 
prevS1.deltam = zeros(Nqp_heights, 1); %Max previous normal displacement
prevS1.Fm = zeros(Nqp_heights, 1);%Max previous normal force
prevS1.am = zeros(Nqp_heights, 1);%contact radius at reversal
prevS = repmat({prevS1}, size(QuadMats.Q, 1), 1);

% Only use Coulomb friction limit, need nonzero stiffness at prestressed
% state. 
pars_pre = pars;
pars_pre.sliptype = 1;

[Xstat, eflg] = fsolve(@(X) NLRES_TMP(CFUN_PRE, [X;0], K, L, Prestress*Fvnorm, Fvnorm*0, prevS, pars_pre, 0, QuadMats, MESH), X0, opt);

if eflg<0
    error('Non-convergent prestress!');
end

uxynStat = QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xstat, 3, [])';

% Update previous information
[tx, ty, tn, ~, ~, ~, ~, ~, ~, ~, dtxdp, dtydp, dtndp, prevS] = ...
    CFUN_PRE(uxynStat, 0, pars_pre, prevS);


% Reference quantities
[tx, ty, tn, dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun, dtxdp, dtydp, dtndp, prevS] = ...
    CFUN(QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xstat, 3, [])', 0, pars_pre, prevS);


% Updated stiffness matrix at prestressed state
[~,dR0,~,~,~] = NLRES_TMP(CFUN, [Xstat;0], K, L, Prestress*Fvnorm, Fvnorm*0, prevS, pars_pre, 0, QuadMats, MESH, ones(size(X0)));

% % Look at the normal displacements
if(~onnots)
    uxyn_stat = QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xstat, 3, [])';
    
    figure;
    semilogy(abs(uxyn_stat(:, 3)), 'LineWidth', 3);
    hold on;
    semilogy((uxyn_stat(:, 3)- zte_gaps), '--', 'LineWidth', 3);
%     semilogy(zte_gaps, 'LineWidth', 3);
    hold on;
    ylabel('un');
    legend('un', 'un - gap');
    set(gcf, 'Renderer', 'painters');


end

%% Eigen Analysis after Prestress

[Vst, Wst] = eigs(dR0, M, 10, 'SM');
[Wst, si] = sort(sqrt(diag(Wst)));
Vst = Vst(:, si);  Vst = Vst./sqrt(diag(Vst'*M*Vst)');

%% Select mode of interest

resp_amp = abs(R(3,:)*Vst);
[val, mis] = maxk(resp_amp, 3);
mis = sort(mis);

mi = mis(mode_ind); %index of the Vst matrix for the mode of interest
assert(abs(val(mode_ind)-1) < 0.1, 'Normalization of experimental data is bad')

Wst([1, 3])/2/pi

if(~onnots)
    keyboard; %dbcont
%     assert(false, 'Stopping'); % Just want the frequencies
end

%% Call RQNMA Solution

% Amplitude Range
As = -7; %log10(startAmp) 
Ae = -4.2; % log10(endAmp)
Qamps = 10.^linspace(As, Ae, 15);
% Qamps = Qamps([end]);

% RQNMA Settings
Nhp = inputpars.Nhp;  % DONT DO EVEN NUMBERS - You land at 0 on the hysteretic paths in that case
optFJ.viscousDamping    = 0;
optFJ.repeatLoop        = inputpars.repeatLoop; % Max number of loop iterations including final one
optFJ.zRelTol           = 1e-5; 
optFJ.calcQuads         = true; % Calculate all quad points on initial iterations, needed if using a relative tolerance or a penalty type model
optFJ.RC_prev           = 2;
optFJ.Nqp_heights       = Nqp_heights; 
optFJ.Nqp_radius        = Nqp_radius; 
optFJ.frictionlessPrestress = true; 
optFJ.CFUN_PRE          = CFUN_PRE; 



% NEED TO UPDATE TO A CORRECT RQNMA BACKBONE MODEL FUNCTION
BB = FJMODEL_BBFUN_LOOPXN(pars, mi, Qamps, K, M, Xstat, Prestress*Fvnorm, L, QuadMats, MESH, CFUN, Nhp, opt, optFJ);

% Save the git hash with the data. 
[~,git_hash_string] = system('git rev-parse HEAD');

% Save results on NOTS before exiting
if(onnots)
    save(inputpars.output_name);
end


% if(~onnots)
%     % If on local do not exit in case you want to look at something
%     keyboard %use "dbcont" to continue
% end



end


%% Function for initializing pars
function [pars, PZFUN, CFUN, CFUN_PRE, area_density, zmax, zmin] ...
    = add_surf_pars(matpars, curr_combine, ASP_FUN, ASP_FUN_PRE, ELEM_TRAC, ...
    Nqp_heights, Nqp_radius, zte_gaps, useSphere)
% Function for updating all of the parameters to the current surface scan
% Inputs:
%   matpars - material parameters
%   curr_combine - combined surface processed data. 
%
%
% NOTES:
%   1. Assumes that the min value of z is moved to zero.

    pars = matpars;

    if(useSphere)
        pars.Rprime = sqrt(curr_combine.Rprime*curr_combine.Rpprime);
        pars.Rpprime = pars.Rprime;
    else
        pars.Rprime = curr_combine.Rprime;
        pars.Rpprime = curr_combine.Rpprime;
    end
    
    pars = INITIALIZE_ECCENTRIC(pars);
    pars.R = pars.Re;


    area_density = curr_combine.area_density;    
    zmax = curr_combine.z_max;
    zmin = 0;
    
    PZFUN = @(z) interp1(curr_combine.normzinterp, curr_combine.pzinterp, z/zmax);
    
    CFUN = @(uxyn,pn0,pars,prev) RC_FULL_TR(uxyn,pn0,pars, ASP_FUN, PZFUN, ELEM_TRAC, Nqp_heights, Nqp_radius, zmin, zmax, area_density, prev, zte_gaps);
    CFUN_PRE = @(uxyn,pn0,pars,prev) RC_FULL_TR(uxyn,pn0,pars, ASP_FUN_PRE, PZFUN, ELEM_TRAC, Nqp_heights, Nqp_radius, zmin, zmax, area_density, prev, zte_gaps);
    
%     % Reduced arguments for element traction
%     ELEM_TRAC_RED = @(uxyn_elem, pars, prev)ELEM_TRAC(pars, uxyn_elem, prev, ASP_FUN, PZFUN, Nqp_heights, Nqp_radius, zmin, zmax, area_density);

end

