function [] = main_epmc(varargin)

nots_min_cores = 10;

if(nargin == 0)
    warning('Defaulting to parameter set 0.');
    inputpars = epmc_pars(0);
    runnum = 0;
else
    fprintf('Loading parameters number %u. \n', varargin{1});
    inputpars = epmc_pars(varargin{1});
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
    workdir = [getenv('SHARED_SCRATCH') filesep getenv('USER') filesep 'microslip-rough-contact' filesep 'EPMC_SIMS'];
    cd(workdir);
    
    mkdir('Results'); % Issues a warning if exists, prevents fatal save error.
    
    onnots = true;
end

% %% Routines locations

addpath('../ROUTINES/')

% Friction Model
addpath('../ROUTINES/FRIC_MODELS')
addpath('../ROUTINES/FRIC_MODELS/ASP_FUN/')
addpath('../ROUTINES/FRIC_MODELS/PLASTIC/')
addpath('../ROUTINES/FRIC_MODELS/MINDLIN/')

% EPMC Specific
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

% ROM/FEM Info
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/RQNMA/') %For the shakedown prestress


% output_name = 'epmc_tmd_abs_run4b';
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

%% Newer FEM Model

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

    % Switch to MESH Class for sake of plots
    Nq = 1;
    MESH = MESH2D(MESH.Nds, 3, [], MESH.Quad, Nq);

    
    % Accelerometer Mass
    maccel = 0.8e-3; % 0.8 grams
    M = M + (R(4:9, :)'*R(4:9, :))*maccel;
    
    fprintf('There are %u Fixed Interface Modes.\n', mod(size(M, 1), 3*MESH.Nn));

end

%% Model Parameters

Prestress = (12002+12075+12670)/3;%12845
mds = [1 3 5];

if(isfield(inputpars, 'mode'))
    mode_ind = inputpars.mode;
else
    mode_ind = 1; %set to 1,2,3 for first three bending modes
end

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


%% Set Discretization

% Number of Quadrature Points
Nqp_radius = 100; %overridden to 1 if using a tangent or secant asperity model below.
Nqp_heights = 100;


%% Input deck style:

if(isequal('BRAKE', upper(inputpars.unloadModel)))
    
    % Plastic Models - Brake unloading
    ASP_FUN = @SPHERE_PL_TAN_DECOUPLE_PL_NORM; modeloutname = 'sphere_pltan_de_plnorm'; Nqp_radius = 1;
    ASP_FUN_PRE = @SPHERE_PLASTIC_PRE;
    
    ELEM_TRAC = @RC_TRACTION_PLASTIC;

elseif(isequal('ELASTIC_ELLIP_TAN', upper(inputpars.unloadModel)))
    
    ASP_FUN = @ELLIPSOID_TAN_DECOUPLE; modeloutname = 'ellip_tan_de'; Nqp_radius = 1;
    
    ASP_FUN_PRE = @ELLIPSOID_PRE;
    ELEM_TRAC = @RC_TRACTION;

else
    
    error('No unloading model defined');

end

if(isfield(inputpars, 'alpha') && inputpars.alpha ~=0)
    
    assert(contains(upper(inputpars.unloadModel), "ELASTIC"), 'Can only apply rotations to elastic models.')
    
    UN_ROT_FUN = ASP_FUN;
    
    ASP_FUN = @(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq)ASP_ROTATION(UN_ROT_FUN, pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq);
    
    pars.alpha = inputpars.alpha;
end



pars.mu = inputpars.mu;
pars.sliptype = inputpars.sliptype;

% Adhesion Energy
if(isfield(inputpars, 'DeltaGamma'))
    pars.DeltaGamma = inputpars.DeltaGamma;
end

%% Processed surface load

matpars = pars;

%%%% Load Surface Processed Info

%%%%% Set 8 - Initial Paper
if(isfield(inputpars, 'surfaceName'))
    surf = load(sprintf('../SURFACE/OUT/%s', inputpars.surfaceName));
else
    surf = load('../SURFACE/OUT/combined_14sep21_R1.mat'); 
end
combined_surf = surf.combinedSurf;
Qps = QuadMats.Q*MESH.Nds;

% Calculate gaps based on meshing the interface at nodes
nd_gaps = -surf.meso_surf{1}(MESH.Nds) - surf.meso_surf{2}(MESH.Nds);
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
    drawnow;
    
    if( isequal(class(MESH), 'MESH2D') )

        set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaultTextInterpreter','latex');

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
        
        xlabel('$x$ [m]');
        ylabel('$y$ [m]');
        zlabel('Gap [m]');
        
        set(gca, 'fontsize', 14);
        c.FontSize = 14;
    
        
        set(gcf, 'Renderer', 'painters');
        drawnow;
%         print('./ltw_mesoscale.svg', '-dsvg', '-r400');
%         print(sprintf('./ltw_mesoscale_%uzte.png', MESH.Ne), '-dpng', '-r600'); %png doesn't have pdf rendering issues in my LaTeX.

        
    end
    
    
end


% combined_surf.z_max = 1.1528e-04; warning('Changing zmax');

%%%%% Potentially Modify Surface Parameters
if(isfield(inputpars, 'Rprime'))
    combined_surf.Rprime  = inputpars.Rprime;
    combined_surf.Rpprime = inputpars.Rpprime;
end

if(isfield(inputpars, 'area_density'))
    combined_surf.area_density = inputpars.area_density;
end


%%%%% End Load Surface Process

%Update parameters with the surface scan.
[pars, PZFUN, CFUN, CFUN_PRE, area_density, zmax, zmin] = add_surf_pars(matpars, ...
    combined_surf, ASP_FUN, ASP_FUN_PRE, ELEM_TRAC, Nqp_heights, ...
    Nqp_radius, zte_gaps, inputpars.useSphere);

clear surf combined_surf
%% Yield Parameters

% Sys = 200e6; %Pa
% Sys = 230e6; %Pa

% Sys = 330e6; %Pa % From current test reference
% Sys = Inf; % Fully elastic

Sys = inputpars.Sys;

C = 1.295*exp(0.736*pars.nu);

delta_y1s = (pi*C*Sys/(2*(2*pars.Estar)))^2*(2*pars.R); %displacement of one sphere against a rigid flat to cause yielding.

% Store all of these results as system parameters to prevent recalculation:
pars.delta_y = delta_y1s*2;
pars.Sys = Sys;

% warning('Fully Elastic-ish');
% pars.delta_y = 100000;
% pars.Sys = 100000*Sys;

pars.Et = inputpars.Et;

if(isfield(inputpars, 'UTS'))
    pars.UTS = inputpars.UTS;
end

%% Harmonic Force Functions

% Just need the 3rd Output to update prev
HARM_INIT = @(uxyn_init, uxyn0, h, t, prev, quadz) RC_TRACTION_PL_HARMONIC_INIT(uxyn_init, uxyn0, h, t, ASP_FUN_PRE, ...
                                        PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density, quadz);
      
% Outputs: [txyn, dtxynduxynh, prev]
HARM_FUN = @(uxyn, h, t, cst, prev, quadz) RC_TRACTION_PL_HARMONIC(uxyn, h, t, ASP_FUN, ...
                    PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density, quadz, cst);

                

%% Setup model / Cleanup

clear C; % Previous C was for plasticity not the damping matrix that this wants.

%% Stuck interface initial guess

uxyn_est = [0, 0, abs(zmax)/2/4/2];
uxyn_est = [0, 0, .1e-5]; % 3e-4=37 588 iters,8e-5=37. 4e-5

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
X0 = K0\(Fv*Prestress);

%% Simple Prestress

opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true, 'MaxIterations', 200); %, 'StepTolerance', 1e-12

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

[Xstat, eflg] = fsolve(@(X) NLRES_TMP(CFUN_PRE, [X;0], K, L, Prestress*Fv, Fv*0, prevS, pars_pre, 0, QuadMats, MESH), X0, opt);

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
[~,dR0,~,~,~] = NLRES_TMP(CFUN, [Xstat;0], K, L, Prestress*Fv, Fv*0, prevS, pars_pre, 0, QuadMats, MESH, ones(size(X0)));

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

    
    if(  isequal(class(MESH), 'MESH2D') )

        figure('Position', [600 400 1200 400]);
        MESH.SHOWFIELD2D(QuadMats.T*tn./sum(QuadMats.T, 2))
        c = colorbar;
        c.Label.String = 'Contact Pressure [Pa]';

        xrange = [min(MESH.Nds(:, 1)), max(MESH.Nds(:, 1))];
        yrange = [min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))];
%         zrange = [0, 1.5e-4];
        
        xlim(xrange);
        ylim(yrange);
%         zlim(zrange);
%         view([-5.5448   36.0170]);
        
        pbaspect([range(xrange), range(yrange), 0.5*range(yrange)])
        
        xlabel('x [m]');
        ylabel('y [m]');
%         zlabel('Gap [m]');
        
        set(gca, 'fontsize', 12);
    
        
        set(gcf, 'Renderer', 'painters');

        % print('./tmd_ltw_mesoscale.svg', '-dsvg', '-r400');

        
    end
    
    

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
    fprintf('%.4f & %.4f & %.4f\n', Wst([1, 3, 5])/2/pi);
%     assert(false, 'Stopping'); % Just want the frequencies
    keyboard %dbcont
end

% Check that if doing second bending mode, it is correctly pulling the
% third system mode.
if(mode_ind == 2)
    assert(mi == 3, '2nd bending mode not using third system mode as expected');
end


%% Damping Factors / Matrix
Zetas = [0.087e-2; 0.034e-2]; 

ab = [1./(2*Wst([1, 3])) Wst([1, 3])/2]\Zetas;
C = ab(1)*M + ab(2)*dR0;

Zetas_wst = 1./(2*Wst)*ab(1) + Wst/2*ab(2);


%% Setup model 

GM = MDOFGEN(M, K, C, L);

clear prev0;

prev0.tx0 = zeros(Nqp_heights, Nqp_radius); %traction for each slider of each asperity height - X
prev0.ty0 = zeros(Nqp_heights, Nqp_radius); %traction for each slider of each asperity height - X
prev0.rq0 = ones(Nqp_heights, 1)*linspace(0, 1, Nqp_radius); %radii for each slip radius used
prev0.uxyw0 = ones(Nqp_heights, 1)*[0, 0, 0]; %displacement at previous step. 
prev0.deltam = zeros(Nqp_heights, 1); %Max previous normal displacement
prev0.Fm = zeros(Nqp_heights, 1);%Max previous normal force
prev0.am = zeros(Nqp_heights, 1);%reversal contact radius


% GM = GM.SETNLFUN(2+5, QuadMats.Q, CFUN, QuadMats.T, prev, HARM_INIT);

% Say we have,
% Cell (Npoints,1) Ls where, Ls{i} is (3 x GM.Ndofs) : Go from global DoFs to local Dofs of nonlinearity "i"
% Cell (Npoints,1) Lf where, Lf{i} is (GM.Ndofs x 3) : Go from nonlinear forces from nonlinearity "i" to global DoFs

Npoints = MESH.Ne*No;

Ls = cell(Npoints, 1);
Lf = cell(Npoints, 1);

% QL = [kron(QuadMats.Q, eye(3)), zeros(size(QuadMats.Q,1)*3, size(L,1)-size(QuadMats.Q,2)*3)]*L;
QL = kron(QuadMats.Q, eye(3))*L(1:3*MESH.Nn, :);
LTT = L(1:3*MESH.Nn, :)'*kron(QuadMats.T, eye(3));

for i=1:Npoints % MESH.Ne*No - > "No" points per element

    Ls{i} = QL(((i-1)*3+1):((i-1)*3+3), :);
    Lf{i} = LTT(:, ((i-1)*3+1):((i-1)*3+3));
    GM = GM.SETNLFUN(2+5, Ls{i}, ...
        @(uxyn, h, t, cst, prev)HARM_FUN(uxyn, h, t, cst, prev, zte_gaps(i)), ...
        Lf{i}, prev0, ...
        @(uxyn_init, uxyn0, h, t, prev)HARM_INIT(uxyn_init, uxyn0, h, t, prev, zte_gaps(i)));
end




%% Backbone
% h = 0:5; %Harmonics 
h = 0:3; %Harmonics - converges
Nhc = sum((h==0)+2*(h~=0)); %Number of Harmonic Coordinates

Nt = 2^7; %Number of time points per AFT, recommended 2^7

Fl = kron([0; 1; 0; zeros(Nhc-3,1)], L(end,:)'); %Pull out the point where the 90deg phase constraint

if(mode_ind == 1)
    As = -7; %log10(startAmp) Use -7 for plotting?
    Ae = -4.2; % log10(endAmp)
elseif(mode_ind == 2)
    As = -7.7; %log10(startAmp) Use -7 for plotting?
    Ae = -4.2; % log10(endAmp)
end

da = inputpars.arcSettings.da; 


% Simulation for mode mi
Uwx0 = [kron([0; 0; 1; zeros(Nhc-3,1)], Vst(:,mi)) + kron([1; zeros(Nhc-1,1)], Xstat); ...
        Wst(mi); 2*Zetas_wst(mi)*Wst(mi)]; %This is: [U0; U1cos; U1sin; . . . ; W; Xi] with each U being the full DOF vector

Dscale = [mean(abs(Xstat))*ones(size(Xstat));
            mean(abs(Vst(:, mi)))*ones(length(Vst(:,mi))*(Nhc-1),1); 
            Wst(mi); 
            100*2*Zetas_wst(mi)*Wst(mi); 
            1];

% K0norm = speye(length(Dscale) - 1);
K0norm = sparse(diag(1./Dscale(1:end-1).^2));
normc = (Uwx0)'*K0norm*(Uwx0);

Copt = struct('Nmax', 1000, 'dsmax', inputpars.arcSettings.dsmax, ...
                'dsmin', inputpars.arcSettings.dsmin, 'DynDscale', 0, 'Dscale', Dscale, ...
                'Display', true, ...
                'itDisplay', true, 'ITMAX', 12, 'crit', 11, ...% NL Solver settings at each step
                'itopt', 7, ... % Optimal number of iterations used to adjust step size
                'lsrch', 0, 'lsit', 2, ... % Line search settings true/false, iters
                'solverchoice', 2, 'arclengthparm', 'K0NormalizedArcLength', ... % Solver choices
                'arcsettings', struct('b', 0.5, 'normc', normc, 'K0', K0norm, 'normb', Dscale(end)^2)); 

if(isfield(inputpars, 'ITMAX'))
    Copt.ITMAX = inputpars.ITMAX;
end

if(onnots)
    itersave_name = ['Results/' output_name '_iter'];
    Copt.callbackfun = @(U, J0, dUdlam) CONTINUE_SAVE(U, J0, dUdlam, itersave_name, [itersave_name '_backup']);
else
    keyboard;
end
            
            
optEPMC = struct('tol', 1e-6, 'maxAFT', 3);


Fl = Fl + kron([1; zeros(Nhc-1,1)], Prestress*Fv); % Add static forces


if(inputpars.load_initial)
    
    prevInd = inputpars.prevInd;

    tmp = load(inputpars.loadName);
    
    Uwx0 = tmp.U(1:length(Uwx0), prevInd);

    As = tmp.U(end, prevInd);

end


% Save the git hash with the data. 
[~,git_hash_string] = system('git rev-parse HEAD');

if(onnots)
    save(['Results/' output_name '_pre']);
end


[UwxC, dUwxC] = CONTINUE(@(uwxa) GM.EPMCRESFUN(uwxa, Fl, h, Nt, optEPMC), Uwx0, As, Ae, da, Copt);

% Nlvib Continuation (try if PRECOCONT fails) (from Malte Krack & Johann Gross)
% Sopt = struct('jac', 'full', 'dynamicDscale', 1);
% UwxCs = solve_and_continue(Uwx0, @(Uwxa) GM.EPMCRESFUN(Uwxa, Fl, h, Nt, 1e-6, 1), As, Ae, da, Sopt);


if(onnots)
    save(['Results/' output_name]);
end


assert(~onnots, 'Stopping here since this is on NOTS\n');


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

