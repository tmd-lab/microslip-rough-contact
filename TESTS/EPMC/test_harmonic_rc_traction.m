% This script is for basic debugging and writing of the element wise force
% function for evaluation of forces at ROM Quadrature points. 
clear;

addpath('../../ROUTINES/')
addpath('../../ROUTINES/FRIC_MODELS/')
addpath('../../ROUTINES/FRIC_MODELS/ASP_FUN/')
addpath('../../ROUTINES/FRIC_MODELS/PLASTIC/')

addpath('../../ROUTINES/HARMONIC/')
addpath('../../ROUTINES/SOLVERS/')

%% Probability Distribution, max height of 0, nominal gap of 0 for highest points. 

% Asperity Heights of z
zmax = 2e-5;
zmin = 0;

z = linspace(zmin, zmax, 1000);

betaA = 4;
betaB = 4;
% betaA = 8;
% betaB = 1;

pz = betapdf((z - zmin)/(zmax - zmin), betaA, betaB);

figure; 
hold on;
xlabel('Asperity Height (m)');
ylabel('Probability');

plot(z, pz, 'LineWidth', 2);

%probability function
p_fun = @(z) betapdf((z - zmin)/(zmax - zmin), betaA, betaB);

%% Parameters

% Number of Quadrature Points
Nqp_radius = 100;
Nqp_heights = 10; %3 tests a single asperity since it goes to zeros on the ends, 10 is probably reasonable for computation times


%% Friction Model


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brake unloading equation
% ASP_FUN = @SPHERE_PLASTIC_PRE;
ASP_FUN = @SPHERE_PL_TAN_DECOUPLE_PL_NORM; Nqp_radius = 1;

ASP_FUN_PRE = @SPHERE_PLASTIC_PRE;

pars.sliptype = 3; % Must be 3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Alternative unloading equation
% % ASP_FUN = @SPHERE_PLASTIC_PRE2;
% ASP_FUN = @SPHERE_PL_TAN_DECOUPLE_PL_NORM2; Nqp_radius = 1;
% 
% ASP_FUN_PRE = @SPHERE_PLASTIC_PRE2;
% 
% pars.sliptype = 3; % Setting divisible by 2 will consider plasticity and Coulomb., Divisible by 3 will do min(others, CEB)

%% Quadrature info

% Radius Quadrature
rq = linspace(0, 1, Nqp_radius);
wq_r = ones(size(rq));
wq_r(2:end-1) = 2;
wq_r = wq_r/sum(wq_r);

% Asperity Height Quadrature
zq = linspace(zmin, zmax, Nqp_heights);
wq_z = ones(size(zq));
wq_z(2:end-1) = 2;
wq_z = wq_z/sum(wq_z);

%% Initialize Parameters ect.

%Initialize the sphere parameters
E = 192.85e9; %304L from Traction paper
nu = 0.30; %304L taken from Traction paper
G = E / 2/ (1 + nu);%304L approx Shear modulus

% pars.mu     = 0.075; % set with friction model, what is reasonable changes
pars.E      = E;
pars.Estar  = E/2/(1 - nu^2);
pars.nu     = nu;
pars.G      = G;
pars.Gstar  = G/2/(2 - nu);


%Ellipsoid initialization, comment out if Rprime=Rpprime.
pars.Rprime = 0.0080; %(1/0.018453533867678 + 1/0.019130030937404)^(-1);
pars.Rpprime = 6.0630e-05; %(1/0.000172616666847 + 1/0.000199388648556)^(-1);
pars = INITIALIZE_ECCENTRIC(pars);
pars.R = pars.Re;


%Area density of asperities
area_density = 0.2*(8.0597e+06);

%elastic dry friction normal tractions on the order of 2e7, tangential 1e7
%% Yield Parameters

Sys = 200e6; %Pa

C = 1.295*exp(0.736*pars.nu);

delta_y1s = (pi*C*Sys/(2*(2*pars.Estar)))^2*(2*pars.R); %displacement of one sphere against a rigid flat to cause yielding.

% Store all of these results as system parameters to prevent recalculation:
pars.delta_y = delta_y1s*2;
pars.Sys = Sys;

%tangent stiffness after yield.
pars.Et = 0.015*pars.E;

pars.mu = 0.6;

%% TEST RC Friction Function


%Previous info
prev.tx0 = zeros(Nqp_heights, Nqp_radius); %traction for each slider of each asperity height - X
prev.ty0 = zeros(Nqp_heights, Nqp_radius); %traction for each slider of each asperity height - X
prev.rq0 = ones(Nqp_heights, 1)*linspace(0, 1, Nqp_radius); %radii for each slip radius used
prev.uxyw0 = ones(Nqp_heights, 1)*[0, 0, 0]; %displacement at previous step. 
prev.deltam = zeros(Nqp_heights, 1); %Max previous normal displacement
prev.Fm = zeros(Nqp_heights, 1);%Max previous normal force
prev.am = zeros(Nqp_heights, 1);%reversal contact radius


prev0 = prev;

PZFUN = p_fun;

%% Test just a single set of displacements

prev = prev0;

Nt = 2^9;
h = [0, 1, 2]';
uxyn0 = [0, .5e-5, 1e-5];
uxyn1c = [0.1e-5, 0.2e-5, 0.1e-5];
uxyn1s = [0.3e-5, 0.2e-5, 0.1e-5];
uxyn2c = [0.1e-5, 0.1e-5, 0.1e-5];

uxynharmonics = [uxyn0; uxyn1c; uxyn1s; uxyn2c; zeros(1, 3)];

%%%%%%%%% Block from EPMCRESFUN
uxyn_t = TIMESERIES_DERIV(Nt, h, uxynharmonics, 0);  % Nt x Ndnl
Nhc = sum((h==0)+2*(h~=0));
t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);  
% sct = w*TIMESERIES_DERIV(Nt, h, eye(Nhc), 1);
%%%%%%%%%%


figure; 
hold on;
plot(uxyn_t(:, 1));
plot(uxyn_t(:, 2));
plot(uxyn_t(:, 3));

legend('x', 'y', 'z');

uxyn_init = uxyn_t(1, :);
[uxyn_init(3), unmax_ind] = max(uxyn_t(:, 3));
uxyn0 = uxyn_t(1, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Call the initialization of the normal history %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% derivative of maximum normal w.r.t. the harmonic coefficients
ddeltamduxynh = reshape(cst(unmax_ind, :), 1, 1, []);
warning('THIS IS ONLY PART OF THE DERIVATIVE, need to also account for the time of the maximum point shifting');

prev.ddeltamduxynh = ddeltamduxynh;
prev.duxyn0duxynh = reshape(cst(1, :), 1, 1, []).*eye(3);

quadz = 0;

[txyn, dtxynduxynh, prev] = RC_TRACTION_PL_HARMONIC_INIT(uxyn_init, uxyn0, h, t(unmax_ind), ...
                            ASP_FUN_PRE, PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density, quadz);

ii = 1;

[txyn, dtxynduxynh] = RC_TRACTION_PL_HARMONIC(uxyn_t(ii, :), h, t(ii), ASP_FUN, ...
    PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density, quadz, cst(ii, :));


%% Test Derivatives


prev = prev0;

derivative_test = 1; %1-6+?
%6 is very slow / high resolution in time

delta = [1e-8, 1e-9, 1e-10, 1e-12]; %displacement variation (m)
% delta = [1e-8, 1e-9, 1e-10]; %displacement variation (m)

area_density = 1; %easier to debug.

switch derivative_test    
    case 1
        
        % Random case
        Nt = 2^5;
        h = [0, 1, 2]';
        uxyn0 = [0, .5e-5, 2e-5];
        uxyn1c = [0.1e-5, 0.2e-5, 0.05e-5];
        uxyn1s = [0.3e-5, 0.2e-5, 0.05e-5];
        uxyn2c = [0.1e-5, 0.1e-5, 0.02e-5];

        uxynharmonics = [uxyn0; uxyn1c; uxyn1s; uxyn2c; zeros(1, 3)];
        
    case 2
        % Separation in the middle
        
        Nt = 2^5;
        h = [0, 1, 2]';
        uxyn0 = [0, .5e-5, 1e-5];
        uxyn1c = [0.1e-5, 0.2e-5, 0.05e-5];
        uxyn1s = [0.3e-5, 0.2e-5, 0.05e-5];
        uxyn2c = [0.1e-5, 0.1e-5, 0.02e-5];

        uxynharmonics = [uxyn0; uxyn1c; uxyn1s; uxyn2c; zeros(1, 3)];
        
            
    case 3
        % Low Amplitude had a bug originally, fixed now
        
        Nt = 2^3;
        h = [0, 1, 2]';
        uxyn0 = [0, .5e-5, 1e-5];
        uxyn1c = [0.1e-5, 0.2e-5, 0.05e-5];
        uxyn1s = [0.3e-5, 0.2e-5, 0.05e-5];
        uxyn2c = [0.1e-5, 0.1e-5, 0.02e-5];

        uxynharmonics = [uxyn0; uxyn1c; uxyn1s; uxyn2c; zeros(1, 3)];
        
        
    case 4 
        % Minimal contact/separation throughout
        
        Nt = 2^3;
        h = [0, 1, 2]';
        uxyn0 = [0, .5e-5, .1e-5];
        uxyn1c = [0.1e-5, 0.2e-5, 0.1e-5];
        uxyn1s = [0.3e-5, 0.2e-5, 0.1e-5];
        uxyn2c = [0.1e-5, 0.1e-5, 0.05e-5];

        uxynharmonics = [uxyn0; uxyn1c; uxyn1s; uxyn2c; zeros(1, 3)];
        
    case 5 
        % Late peak contact
        Nt = 2^7;
        h = [0, 1, 2]';
        uxyn0 = [0, .5e-5, -.1e-5];
        uxyn1c = [0.1e-5, 0.2e-5, -0.1e-5];
        uxyn1s = [0.3e-5, 0.2e-5, -0.1e-5];
        uxyn2c = [0.1e-5, 0.1e-5, -0.5e-5];

        uxynharmonics = [uxyn0; uxyn1c; uxyn1s; uxyn2c; zeros(1, 3)];
        
    otherwise
        % Very high resolution test 5 to see if there is notable error for
        % the time of max normal varying with coefficients
        
        % Late peak contact
        Nt = 2^10;
        h = [0, 1, 2]';
        uxyn0 = [0, .5e-5, -.1e-5];
        uxyn1c = [0.1e-5, 0.2e-5, -0.1e-5];
        uxyn1s = [0.3e-5, 0.2e-5, -0.1e-5];
        uxyn2c = [0.1e-5, 0.1e-5, -0.5e-5];

        uxynharmonics = [uxyn0; uxyn1c; uxyn1s; uxyn2c; zeros(1, 3)];
        
        
        
end


% Generation function
GEN_FUN = @(uxynharmonics) generate_time_series(uxynharmonics, Nt, h, ASP_FUN, ...
            ASP_FUN_PRE, PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density);


% Reference
[txyn_t, dtxynduxynh_t, t, uxyn_t] = GEN_FUN(uxynharmonics);


% Checking that reference makes useful sense
figure; 
hold on;
plot(uxyn_t(:, 1));
plot(uxyn_t(:, 2));
plot(uxyn_t(:, 3));

legend('x', 'y', 'z');

% txyn_t

dtxynduxynh_t_num = zeros(size(dtxynduxynh_t));

errors = zeros(size(delta));

for ii = 1:length(delta)
    for dir = 1:3
        for harm = 1:size(uxynharmonics,1)
            
            uxynharmonics_high = uxynharmonics;
            uxynharmonics_low  = uxynharmonics;
            
            uxynharmonics_high(harm, dir) = uxynharmonics(harm, dir) + delta(ii);
            uxynharmonics_low(harm, dir)  = uxynharmonics(harm, dir) - delta(ii);
            
            [txyn_t_high, ~, ~, ~] = GEN_FUN(uxynharmonics_high);
            [txyn_t_low,  ~, ~, ~] = GEN_FUN(uxynharmonics_low);
            
            dtxynduxynh_t_num(:, :, dir, harm) = (txyn_t_high - txyn_t_low)/(2*delta(ii));

        end
    end
        
    errors(ii) = norm(abs(dtxynduxynh_t_num(:) - dtxynduxynh_t(:)))/(norm(dtxynduxynh_t(:)) + (norm(dtxynduxynh_t(:)) == 0));
    
end

log10error = log10(errors)

%% Look at matrices side by side

for dir = 1:size(dtxynduxynh_t, 3)
    for harm = 1:size(dtxynduxynh_t, 4)
        
        fprintf('Varying displacement direction %d, harmonic component %d\n', dir, harm);
        
        [dtxynduxynh_t(:, :, dir, harm), NaN(size(dtxynduxynh_t(:, 1,1,1))) dtxynduxynh_t_num(:, :, dir, harm)]
        
        3
    end
end

disp('Make sure the both for loops cover all entries');

%% Function to generate time series
function [txyn_t, dtxynduxynh_t, t, uxyn_t] = generate_time_series(uxynharmonics, Nt, h, ASP_FUN, ASP_FUN_PRE, PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density)


    %%%%%%%%% Block from EPMCRESFUN
    uxyn_t = TIMESERIES_DERIV(Nt, h, uxynharmonics, 0);  % Nt x Ndnl
    Nhc = sum((h==0)+2*(h~=0));
    t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);  
    % sct = w*TIMESERIES_DERIV(Nt, h, eye(Nhc), 1);
    %%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Call the initialization of the normal history %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Find maximum displacement
    uxyn_init = uxynharmonics(1, :); % Start Tangential Models with Zeroth Harmonic displacements (phase invariance if does not slip).
    [uxyn_init(3), unmax_ind] = max(uxyn_t(:, 3));
%     uxyn0 = uxyn_t(1, :);
    uxyn0 = uxyn_init;

    % derivative of maximum normal w.r.t. the harmonic coefficients
    ddeltamduxynh = reshape(cst(unmax_ind, :), 1, 1, []);

    prev.ddeltamduxynh = ddeltamduxynh;
%     prev.duxyn0duxynh = reshape(cst(1, :), 1, 1, []).*eye(3);
    prev.duxyn0duxynh = zeros(3,3,size(cst, 2));
    prev.duxyn0duxynh(:, :, 1) = eye(3);
    prev.duxyn0duxynh(3, 3, :) = prev.ddeltamduxynh;
    

    quadz = 0;
    
    [~, ~, prev] = RC_TRACTION_PL_HARMONIC_INIT(uxyn_init, uxyn0, h, t(unmax_ind), ...
                                ASP_FUN_PRE, PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density, quadz);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Construct Time history %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialize time history
    txyn_t = zeros(length(t), 3);
    dtxynduxynh_t = zeros(length(t), 3, 3, Nhc);
                            
    % Construct the full time history
    for ii = 1:length(t)
        
        [txyn_t(ii, :), dtxynduxynh_t(ii, :, :, :), prev] = RC_TRACTION_PL_HARMONIC(uxyn_t(ii, :), h, t(ii), ASP_FUN, ...
            PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density, quadz, cst(ii, :));

    end

end
