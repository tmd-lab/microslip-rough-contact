% This script is for testing the force function for a single asperity
% against different cases. 

clear;

addpath('../ROUTINES/FRIC_MODELS/ASP_FUN/')
addpath('../ROUTINES/FRIC_MODELS/MINDLIN/')

%% Define Parameters 
% Parameters are shared between all models 

%Initialize the sphere parameters
E = 192.85e9; %304L from Traction paper
nu = 0.3; %304L taken from Traction paper
G = E / 2/ (1 + nu);%304L approx Shear modulus

%Mean of top and bottom
Rx = mean([0.000273924650289, 0.000314503650665]); 
Ry = mean([0.008183455622513, 0.008521098239706]);

pars.mu     = 0.3; %may want to do a function of displacement - see Eriten's work
pars.R      = (1/Rx + 1/Ry)^(-1);%, 1e-4 is about a maximum
pars.E      = E;
pars.Estar  = E/2/(1 - nu^2);
pars.nu     = nu;
pars.G      = G;
pars.Gstar  = G/2/(2 - nu);



pars.Rprime = (1/0.008183455622513 + 1/0.008521098239706)^(-1);
pars.Rpprime = (1/0.000273924650289 + 1/0.000314503650665)^(-1);

% warning('Testing spheres');
% pars.Rpprime = pars.Rprime; % Test Spheres

pars = INITIALIZE_ECCENTRIC(pars);

pars.R = pars.Re;


%% Choose the Asperity Contact Model To Test

Nqp_radius = 100; %for TAN/SEC anything >1 is ignored and is just wasted memory. 

%Ellipsoids:

% ASP_FUN = @ELLIPSOID_PRE; %prestress, normal only 
% ASP_FUN = @ELLIPSOID_TAN_DECOUPLE; Nqp_radius = 1;
% ASP_FUN = @ELLIPSOID_IWAN_DECOUPLE;
% ASP_FUN = @ELLIPSOID_IWAN_FIT_DECOUPLE;
ASP_FUN = @ELLIPSOID_IWAN_FIT_COUPLE;
% ASP_FUN = @ELLIPSOID_TAN_DECOUPLE_ADHESION; Nqp_radius = 1; pars.DeltaGamma = 5; % J/m^2

%Analytical Function to Compare to:
ANALYTICAL_FUN = @MINDLIN_MONOTONIC;


%Chose direction to test, 1 = x, 2 = y
xynum = 1;


%Potential Functions of the form:
%(SPHERE or ELLIPSOID)_(TAN, SEC, or IWAN)_(COUPLE or DECOUPLE)_(optional PRE)

%Functions to make comparison graphs
ASP_FUN_LIST = {@ELLIPSOID_TAN_DECOUPLE, @ELLIPSOID_TAN_DECOUPLE, ...
                @ELLIPSOID_IWAN_DECOUPLE, @ELLIPSOID_IWAN_DECOUPLE, ...
                @ELLIPSOID_IWAN_FIT_DECOUPLE, @ELLIPSOID_IWAN_FIT_DECOUPLE};
            
ASP_NAMES_LIST = {  'Ellipsoid Tanget x', 'Ellipsoid Tanget y', ...
                    'Ellipsoid Iwan x', 'Ellipsoid Iwan y', ...
                    'Ellipsoid Fit Iwan x', 'Ellipsoid Fit Iwan y'};
                
ASP_NQP_RAD_LIST = [1, 1, 100, 100, 100, 100];
ASP_LIST_XYINDEX = [1, 2, 1, 2, 1, 2];

SPHERE_TESTS = [false, false, false, false, false, false]; %Which listed items are spherical contact v. elliptical


%% Initialize General Info

tx0 = zeros(1, Nqp_radius); 
ty0 = zeros(1, Nqp_radius); 


rq = linspace(0, 1, Nqp_radius); 
wq = ones(size(rq));
wq(2:end-1) = 2;
wq = wq/sum(wq);

rq0 = linspace(0, 1, Nqp_radius); 


uxyn0 = [0, 0, 0]; %displacement at previous step. 

%% Simple Test at one Point


uxyn0 = [0, 0, 0];
% uxyn = [0.5*2e-6, 0.5e-6, 2e-5];
uxyn = [200e-6, 100e-6, 2e-5]; %sqrt(sum(fxyn(1:2).^2))/fxyn(3)


[fxyn, dfxynduxyn, ~, ~, ~] = ASP_FUN(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq);


% Verification of initial implementation, the old code is no longer kept up
% to date, so do not use the following.
%
% % 10/20/2020 - dividing by two since Mindlin displacement may be only one
% % sphere relative to the contact patch.
% kt_fun = @(rq, w, N, a, pars) ones(size(rq)).*8.*pars.Gstar/pi/a/2;
% dktdun_fun = @(rq, w, N, a, pars) -ones(size(rq)).*4.*pars.Gstar.*pars.R/pi/(w*pars.R).^(1.5)/2; 
% 
% [f, ~, ~, dfdx, dfdn] = MINDLIN_IWAN(rq, wq, fxyn(3), uxyn(1), uxyn(3), uxyn0(1), rq0, tx0, kt_fun, dktdun_fun, pars);
% 
% [(f - fxyn(1))./fxyn(1), f, fxyn]

%% Quick Timing Test

N_run = 10000;

%Fully stuck + some displacement test
uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
uxyn2 = [0.5e-6, 0.5e-6, 2e-5];

[~, ~, rq1, tx1, ty1] = ASP_FUN(pars, uxyn1, uxyn0, rq0, tx0, ty0, rq, wq);

tic;
for ii = 1:N_run
    [fxyn, dfxynduxyn, ~, ~, ~] = ASP_FUN(pars, uxyn2, uxyn1, rq1, tx1, ty1, rq, wq);
end

time = toc;

fprintf('Average time is % .5s \n', time/N_run);

%Ellipsoid Iwan Decouple (old derivatives): N_run = 10000;
% Nqp_radius = 100; time: 1.11392e-04 
% Nqp_radius = 300; time: 3.54488e-04 
% Nqp_radius = 1000; time: 2.02723e-03 

%Ellipsoid Iwan Decouple (new derivatives): N_run = 10000;
% Nqp_radius = 100; time: 2.51009e-04
% Nqp_radius = 300; time: 2.88068e-04
% Nqp_radius = 1000; time: 4.04922e-04 

%% Derivative Test at a point

%Initialize at 0, then go to 1, then test derivative at uxyn2

uxyn0 = [0, 0, 0];

derivative_test = 1; %1-8

switch derivative_test
    case 1
        %Test 1 that was used for the sphere models
        uxyn1 = [0, 0, 1.9e-5]; %z must be greater than zero. 
        uxyn2 = [1e-6, 0.5e-6, 2e-5];

    case 2
        %Fully slipped test
        uxyn1 = [0, 0, 1.9e-5]; %z must be greater than zero. 
        uxyn2 = [2e-6, 5e-6, 2e-5];

    case 3
        %Fully stuck test
        uxyn1 = [0, 0, 1.9e-5]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 2e-5];

    case 4 
        %Fully stuck + some displacement test
        uxyn1 = [0, 0, 1.9e-5]; %z must be greater than zero. 
        uxyn2 = [0.5e-6, 0.5e-6, 2e-5];

    case 5
        %Fully stuck test
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [0e-6, 0e-6, 2e-5];

    case 6
        %Fully stuck + some displacement test
        uxyn1 = [0, 0, 2e-5]; %z must be greater than zero. 
        uxyn2 = [0.5e-6, 0.5e-6, 2e-5];

    case 7
        %Re-establishment of contact case. - microslip
        uxyn1 = [0, 0, -2e-5]; %z must be greater than zero. 
        uxyn2 = [0.1e-6, 0.1e-6, 2e-5];

    otherwise
        %Re-establishment of contact case. - macroslip
        uxyn1 = [0, 0, -2e-5]; %z must be greater than zero. 
        uxyn2 = [100e-6, 100e-6, 2e-5];
end

if(uxyn1(3) > 0)
    %Initialize
    [~, ~, rq1, tx1, ty1] = ASP_FUN(pars, uxyn1, uxyn0, rq0, tx0, ty0, rq, wq);
else
    rq1 = rq;
    tx1 = zeros(size(rq));
    ty1 = zeros(size(rq));
end


%Baseline/analytical solution
[fxyn, dfxynduxyn, ~, ~, ~] = ASP_FUN(pars, uxyn2, uxyn1, rq1, tx1, ty1, rq, wq);

% [fxyn, dfxynduxyn, ~, ~, ~, ~, duxyn1duxyn2] = ASP_FUN(pars, uxyn2, uxyn1, rq1, tx1, ty1, rq, wq);
% num_duxynprev = zeros(size(duxyn1duxyn2));



% delta = 1e-8; %displacement variation (m)
delta = [1e-8, 1e-9, 1e-10, 1e-12]; %displacement variation (m)

errorK = zeros(size(delta));

numK = zeros(size(dfxynduxyn));

for ii = 1:length(delta)
    
    for jj = 1:3
        
        uxyn_low = uxyn2;
        uxyn_high = uxyn2;
        
        uxyn_low(jj) = uxyn_low(jj) - delta(ii);
        uxyn_high(jj) = uxyn_high(jj) + delta(ii);
        
        [fxyn_high, ~, ~, ~, ~] = ASP_FUN(pars, uxyn_high, uxyn1, rq1, tx1, ty1, rq, wq);
        [fxyn_low, ~, ~, ~, ~] = ASP_FUN(pars, uxyn_low, uxyn1, rq1, tx1, ty1, rq, wq);
        
%         [fxyn_high, ~, ~, ~, ~, uxynprev_high] = ASP_FUN(pars, uxyn_high, uxyn1, rq1, tx1, ty1, rq, wq);
%         [fxyn_low,  ~, ~, ~, ~, uxynprev_low]  = ASP_FUN(pars, uxyn_low, uxyn1, rq1, tx1, ty1, rq, wq);
        
        
        numK(:, jj) = (fxyn_high' - fxyn_low')/(2*delta(ii));
        
%         num_duxynprev(:, jj) = (uxynprev_high' - uxynprev_low')/(2*delta(ii));

    end
    
%     num_duxynprev
    
    errorK(ii) = norm(dfxynduxyn - numK)/norm(dfxynduxyn);
    
    
end

errorK

% %% Debugging more derivatives

% [~, ~, ~, ~, ~, ~, dadun, ~, dtsdun, ~, dktdun] = ASP_FUN(pars, uxyn2, uxyn1, rq1, tx1, ty1, rq, wq);
% [~, ~, ~, ~, ~, a_low, ~, ts_low, ~, kt_low, ~] = ASP_FUN(pars, uxyn_low, uxyn1, rq1, tx1, ty1, rq, wq);
% [~, ~, ~, ~, ~, a_high, ~, ts_high, ~, kt_high, ~] = ASP_FUN(pars, uxyn_high, uxyn1, rq1, tx1, ty1, rq, wq);
% 
% num_dadun = (a_high - a_low) / (2*delta(ii));
% num_dtsdun = (ts_high - ts_low) / (2*delta(ii));
% num_dktdun = (kt_high - kt_low) / (2*delta(ii));
% 
% [num_dtsdun; dtsdun]
% 
% a_error = (num_dadun- dadun)./dadun
% ts_error = sum(num_dtsdun - dtsdun)/sum(num_dtsdun)
% kt_error = (num_dktdun - dktdun)/dktdun


%% Normal Loading Curve verification

un = linspace(0, 1e-5, 100);

fn_ASP = zeros(size(un));

uxyn0 = [0, 0, 0];

prev0.uxnyn  = [0, 0, 0, 0]; %From the previous step
prev0.fxnyn  = [0, 0, 0, 0];
prev0.uxnyn0 = [0, 0, 0, 0]; %From previous reversal point
prev0.fxnyn0 = [0, 0, 0, 0];
prev0.uxynOrigin = [0, 0, 0]; %Origin for where to consider displacement. 


%Mindlin normal force
fxyn_MINDLIN = MINDLIN_MONOTONIC([zeros(length(un), 2), un'], pars);
fn_MINDLIN = fxyn_MINDLIN(:, 3)';

% Elliptical contact normal force
[fxyn_ELLIPSOID, uxy_slip] = ELLIPSOID_MONOTONIC([zeros(length(un), 2), un'], pars);
fn_ELLIPSOID = fxyn_ELLIPSOID(:, 3)';

for ii = 1:length(un)

    uxyn_tmp(3) = un(ii);
    
    %Asperity function
    [fxyn, ~, ~, ~, ~] = ASP_FUN(pars, uxyn_tmp, uxyn0, rq0, tx0, ty0, rq, wq);
    
    fn_ASP(ii) = fxyn(3);

end


figure;
hold on;
xlabel('un');
ylabel('fn');

plot(un, fn_MINDLIN, '-', 'LineWidth', 2);
plot(un, fn_ELLIPSOID, '--', 'LineWidth', 2);
plot(un, fn_ASP, ':', 'LineWidth', 2);

legend('Hertz Sphere', 'Hertz Ellipsoid', 'Asperity');


%% Complicated Tests + Comparisons to Mindlin Theory

hyst_amp = 1; 
if(hyst_amp > 1)
    warning('If the system has entered macroslip, the current transition method is not correct');
end

r0 = 0; %Gap between the asperities
N_test = 1000;
area_density = 0.2*(8.0597e+06);



% %% First Set of Hysteresis Loop
un = 1e-5;

[fxyn, u_slip] = ANALYTICAL_FUN([0, 0, un], pars);
if(length(u_slip) > 1)
    u_slip = u_slip(xynum); %always use x direction slip for cases of ellipsoid.
end
N1 = fxyn(1, 3);

ux_test{1} = [linspace(0, hyst_amp*u_slip, N_test), ...
          linspace(hyst_amp*u_slip, -hyst_amp*u_slip, N_test), ...
          linspace(-hyst_amp*u_slip, hyst_amp*u_slip, N_test)];
      
un_test{1} = un*ones(size(ux_test{1}));



%Reversal force
[fxyn, ~] = ANALYTICAL_FUN([hyst_amp*u_slip, 0, un], pars);
fx_star = fxyn(1, 1);

% % %% Second Set of Hysteresis Loop
un = 1.2e-5;
uxyn_tmp(1, 1:2) = 0;
uxyn_tmp(1, 3)  = un;

[fxyn, ~] = ANALYTICAL_FUN([0, 0, un], pars);
N2 = fxyn(1,3);

%Displacement to match (fx_star + mu N)
a = sqrt((uxyn_tmp(1,3)-r0).*pars.R); 
muDeltaN = pars.mu*(N2 - N1);
u_max = 3*pars.mu*fxyn(1, 3)/16/pars.Gstar/a*(1 - (1 - (fx_star + muDeltaN)/pars.mu/fxyn(1,3))^(2/3));

ux_test{2} = [linspace(0, u_max, N_test), ...
          linspace(u_max, -u_max, N_test), ...
          linspace(-u_max, u_max, N_test)];
      
un_test{2} = un*ones(size(ux_test{1}));

% %% Third Set of Hysteresis Loop
un = 1.4e-5;
uxyn_tmp(1, 1:2) = 0;
uxyn_tmp(1, 3)  = un;

[fxyn, ~] = ANALYTICAL_FUN([0, 0, un], pars);
N3 = fxyn(1,3);
    
% u_slip = 3*pars.mu*fxyn(1, 3)/16/pars.Gstar/sqrt((uxyn_tmp(1, 3)-r0)*pars.R/2);
%Displacement to match (fx_star + mu N)
a = sqrt((uxyn_tmp(1,3)-r0).*pars.R); 
muDeltaN = pars.mu*(N3 - N1);
u_max = 3*pars.mu*fxyn(1, 3)/16/pars.Gstar/a*(1 - (1 - (fx_star + muDeltaN)/pars.mu/fxyn(1,3))^(2/3));

ux_test{3} = [linspace(0, u_max, N_test), ...
          linspace(u_max, -u_max, N_test), ...
          linspace(-u_max, u_max, N_test)];
      
un_test{3} = un*ones(size(ux_test{1}));


% %% Fourth Hysteresis loop with variable Normal load

ind_set1 = 1:(N_test*1.5);
ind_set2 = N_test*1.5 + (1:N_test);
ind_set3 = N_test*2.5 + (1:(0.5*N_test));

ux_test{4} = [ux_test{1}(ind_set1), ux_test{2}(ind_set2), ux_test{3}(ind_set3)];
un_test{4} = [un_test{1}(ind_set1), un_test{2}(ind_set2), un_test{3}(ind_set3)];


% %% Fifth Hysteresis loop with decreasing normal loads

ind_set1 = 1:(N_test*1.5);
ind_set2 = N_test*1.5 + (1:N_test);
ind_set3 = N_test*2.5 + (1:(0.5*N_test));

ux_test{5} = [ux_test{3}(ind_set1), ux_test{2}(ind_set2), ux_test{1}(ind_set3)];
un_test{5} = [un_test{3}(ind_set1), un_test{2}(ind_set2), un_test{1}(ind_set3)];


% %% End Making Hysteresis loop inputs

symbol_set(1:3) = {'-'};
symbol_set(4:5) = {'--'};


fx_test = cell(size(ux_test));
fx_test_ASP = fx_test;

for ii = 1:length(ux_test)
    fx_test{ii} = zeros(size(ux_test{ii}));
    fx_test_ASP{ii} = zeros(size(ux_test{ii}));
end

% Generate the fx_test for the first three outputs with the masing
% hypotheses
for ii = 1:3
    
    [fxyn, ~] = ANALYTICAL_FUN([ux_test{ii}(1:N_test)', zeros(N_test, 1), un_test{ii}(1:N_test)'], pars);

    %monotonic loading
    fx_test{ii}(1:N_test) = fxyn(:, 1)';
    
    %hysteretic unloading
    fx_test{ii}((N_test+1):(2*N_test)) = fxyn(N_test, 1)' - 2*fxyn(:, 1)';
    
    
    %hysteretic unloading
    fx_test{ii}((2*N_test+1):(3*N_test)) = -fxyn(N_test, 1)' + 2*fxyn(:, 1)';
    
end

%take the appropriate parts of the Mindlin loading
for ii = 4:length(fx_test)
    
    fx_test{ii} = (un_test{ii} == un_test{1}).*fx_test{1} ... 
                    + (un_test{ii} == un_test{2}).*fx_test{2} ... 
                    + (un_test{ii} == un_test{3}).*fx_test{3};
    
end


t0 = zeros(size(rq));


for ii = 1:length(ux_test)
    
    uxyn_tmp(1, 1:2) = 0;
    uxyn_tmp(1, 3)  = un_test{ii}(1);
    
    %Initialize to the first normal load for the ASP_FUN as well 
    uxyn0_ASP = [0,0,0];
    uxyn_ASP = uxyn_tmp;
    rq0_ASP = rq;
    tx0_ASP = zeros(size(rq));
    ty0_ASP = zeros(size(rq));
    
    [fxyn, dfxynduxyn, rq0_ASP, tx0_ASP, ty0_ASP] = ASP_FUN(pars, uxyn_ASP, uxyn0_ASP, rq0_ASP, tx0_ASP, ty0_ASP, rq, wq);
    
    uxyn0_ASP = uxyn_ASP;

    
    for jj = 1:length(ux_test{ii})
        
        %Asperity function test.        
        uxyn_ASP(1,xynum) = ux_test{ii}(jj);
        uxyn_ASP(1, 3)    = un_test{ii}(jj);
        
        [fxyn, ~, rq0_ASP, tx0_ASP, ty0_ASP] = ASP_FUN(pars, uxyn_ASP, uxyn0_ASP, rq0_ASP, tx0_ASP, ty0_ASP, rq, wq);

        uxyn0_ASP = uxyn_ASP;
        
        fx_test_ASP{ii}(jj) = fxyn(xynum);
        
    end
end



figure; 
hold on;
xlabel('Displacement');
ylabel('Approx Traction');

for ii = 1:length(ux_test)
    
    plot(ux_test{ii}, fx_test{ii}*area_density, '-', 'LineWidth', 2);
    
    plot(ux_test{ii}, fx_test_ASP{ii}*area_density, '--', 'LineWidth', 2);
    
end

legend('Analytical Solution', 'Iwan Representation');


for ii = 1:length(ux_test)
    
    figure; 
    hold on;
    xlabel('Displacement');
    ylabel('Force');
    title('Compare Models');

    plot(ux_test{ii}, fx_test{ii}, '-', 'LineWidth', 2);
    plot(ux_test{ii}, fx_test_ASP{ii}, '--', 'LineWidth', 2);
    
    legend('Mindlin Solution (w/0 tranisitions)', 'Asperity Function');
end

%% Test all of the functions in the asperity function list

% ASP_FUN_LIST 
% ASP_NAMES_LIST 
% ASP_NQP_RAD_LIST 

hyst_amp = 1.5; 

% %% First Set of Hysteresis Loop
un = 1e-5;

[fxyn, u_slip] = MINDLIN_MONOTONIC([0, 0, un], pars);
N1 = fxyn(1,3);  

ux_test = [linspace(0, hyst_amp*u_slip, N_test), ...
          linspace(hyst_amp*u_slip, -hyst_amp*u_slip, N_test), ...
          linspace(-hyst_amp*u_slip, hyst_amp*u_slip, N_test)];
      
un_test = un*ones(size(ux_test));

% Initialize force storage
fx_test = zeros(size(ux_test));
fx_test_ASP = cell(size(ASP_FUN_LIST));

tx0 = cell(size(ASP_FUN_LIST));
ty0 = cell(size(ASP_FUN_LIST));
rq = cell(size(ASP_FUN_LIST));
rq0 = cell(size(ASP_FUN_LIST));
wq = cell(size(ASP_FUN_LIST));
uxyn0_ASP = cell(size(ASP_FUN_LIST));

for ii = 1:length(ASP_FUN_LIST)
    fx_test_ASP{ii} = zeros(size(ux_test));
    
    tx0{ii} = zeros(1, ASP_NQP_RAD_LIST(ii));
    ty0{ii} = zeros(1, ASP_NQP_RAD_LIST(ii));
    rq{ii} = linspace(0, 1, ASP_NQP_RAD_LIST(ii)); 
    rq0{ii} = linspace(0, 1, ASP_NQP_RAD_LIST(ii)); 
    
    
    wq{ii} = ones(size(rq{ii}));
    wq{ii}(2:end-1) = 2;
    wq{ii} = wq{ii}/sum(wq{ii});
    
    uxyn0_ASP{ii} = [0,0,0];
end

%%%% Generate Mindlin Solution data w/Masing assumptions
[fxyn, ~] = MINDLIN_MONOTONIC([ux_test(1:N_test)', zeros(N_test, 1), un_test(1:N_test)'], pars);

%monotonic loading
fx_test(1:N_test) = fxyn(:, 1)';

%hysteretic unloading
fx_test((N_test+1):(2*N_test)) = fxyn(N_test, 1)' - 2*fxyn(:, 1)';

%hysteretic unloading
fx_test((2*N_test+1):(3*N_test)) = -fxyn(N_test, 1)' + 2*fxyn(:, 1)';

%%%%%%%%%%%%LOOP TO UPDATE GIVEN LIST OF ASPERITY MODELS

for ii = 1%:length(ux_test)
    
    %Initialize to the first normal load for the ASP_FUN as well 
    for kk = 1:length(ASP_FUN_LIST)
        
        uxyn_ASP = [0, 0, un_test(1)];
        
        [fxyn, dfxynduxyn, rq0{kk}, tx0{kk}, ty0{kk}] = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk});

        uxyn0_ASP{kk} = uxyn_ASP;
        
    end
    

    
    for jj = 1:length(ux_test)

        
        %Loop over asperity functions
        for kk = 1:length(ASP_FUN_LIST)
        
            xynum_asp = ASP_LIST_XYINDEX(kk);
            
            %Asperity function test. 
            uxyn_ASP = [0, 0, 0];
            uxyn_ASP(1,xynum_asp) = ux_test(jj);
            uxyn_ASP(1, 3)    = un_test(jj);
            
            [fxyn, dfxynduxyn, rq0{kk}, tx0{kk}, ty0{kk}] = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk});

            uxyn0_ASP{kk} = uxyn_ASP;
            
            fx_test_ASP{kk}(jj) = fxyn(xynum_asp);
%             jj

        end
        
        
    end
end


% Plot the comparison
figure; 
hold on;
xlabel('ux/u slip');
ylabel('Fx/mu N');

xlim(hyst_amp*[-1.1, 1.1]);
ylim([-1.1, 1.1]);

legend_names = cell(length(fx_test_ASP)+1, 1);

plot(ux_test/u_slip, fx_test/pars.mu/N1, '-', 'LineWidth', 2.5);
legend_names{1} = 'Mindlin Sphere';

symbol_set = {'--', '-.', ':','--', '-.', ':','--', '-.', ':', '--', '-.', ':'};

for kk = 1:length(fx_test_ASP)
    
    plot(ux_test/u_slip, fx_test_ASP{kk}/pars.mu/N1, symbol_set{kk}, 'LineWidth', 2.5);
    
    legend_names{kk+1} = ASP_NAMES_LIST{kk};
    
end

legend(legend_names, 'Location', 'se');


%%%%%% Look at dissipation with the first N_test entries and Masing
%%%%%% assumptions

D_mindlin = cumtrapz(ux_test(1:N_test)/u_slip, fx_test(1:N_test)/pars.mu/N1);
D_mindlin = 8*D_mindlin - 4 * ux_test(1:N_test)/u_slip.*fx_test(1:N_test)/pars.mu/N1;


D_ASP = cell(size(fx_test_ASP));

for kk = 1:length(fx_test_ASP)
    
    D_ASP{kk} = cumtrapz(ux_test(1:N_test)/u_slip, fx_test_ASP{kk}(1:N_test)/pars.mu/N1);
    D_ASP{kk} = 8*D_ASP{kk} - 4 * ux_test(1:N_test)/u_slip.*fx_test_ASP{kk}(1:N_test)/pars.mu/N1;
    
end


% Dissipation v. displacement
figure; 
hold on;
xlabel('ux/u slip');
ylabel('Dissipation [non-dimensional]');

xlim(hyst_amp*[0, 1.1]);

plot(ux_test(1:N_test)/u_slip, D_mindlin, '-', 'LineWidth', 2.5);

for kk = find(SPHERE_TESTS) %1:3%length(fx_test_ASP)
    
    plot(ux_test(1:N_test)/u_slip, D_ASP{kk}, symbol_set{kk}, 'LineWidth', 2.5);
        
end

legend(legend_names, 'Location', 'se');



% Dissipation v. force
figure; 
hold on;
xlabel('fx/f slip');
ylabel('Dissipation [non-dimensional]');

xlim(hyst_amp*[0, 1.1]);

plot(fx_test(1:N_test)/pars.mu/N1, D_mindlin, '-', 'LineWidth', 2.5);

for kk = find(SPHERE_TESTS) %1:3%length(fx_test_ASP)
    
    plot(fx_test_ASP{kk}(1:N_test)/pars.mu/N1, D_ASP{kk}, symbol_set{kk}, 'LineWidth', 2.5);
        
end

legend(legend_names, 'Location', 'se');


% Error in dissipation calculation. 
figure; 
hold on;
xlabel('ux/u slip');
ylabel('Dissipation [non-dimensional]');

xlim(hyst_amp*[0, 1.1]);

% plot(ux_test(1:N_test)/u_slip, D_mindlin, '-', 'LineWidth', 2.5);

for kk = find(SPHERE_TESTS) %1:3%length(fx_test_ASP)
    
    plot(ux_test(1:N_test)/u_slip, (D_ASP{kk}-D_mindlin)./D_mindlin, symbol_set{kk}, 'LineWidth', 2.5);
        
end

legend(legend_names(2:(sum(SPHERE_TESTS)+1)), 'Location', 'se');


% Error in Secant stiffness/Secant Stiffness
figure; 
hold on;
xlabel('un/u slip');
ylabel('Secant Stiffness [non-dimensional]');

xlim(hyst_amp*[0, 1.1]);

secStiff = (fx_test(1:N_test)/pars.mu/N1)./(ux_test(1:N_test)/u_slip);

plot(ux_test(1:N_test)/u_slip, secStiff, '-', 'LineWidth', 2.5);

for kk = find(SPHERE_TESTS) %1:3%length(fx_test_ASP)
    
    secStiff = (fx_test_ASP{kk}(1:N_test)/pars.mu/N1)./(ux_test(1:N_test)/u_slip);
    
    plot(ux_test(1:N_test)/u_slip, secStiff, symbol_set{kk}, 'LineWidth', 2.5);
        
end

legend(legend_names(2:(sum(SPHERE_TESTS)+1)), 'Location', 'ne');



% Error in Secant stiffness/Secant Stiffness
figure; 
hold on;
xlabel('un/u slip');
ylabel('Frac Error in Secant Stiffness [non-dimensional]');

xlim(hyst_amp*[0, 1.1]);

secStiff_Mindlin = (fx_test(1:N_test)/pars.mu/N1)./(ux_test(1:N_test)/u_slip);

for kk = find(SPHERE_TESTS) %1:3%length(fx_test_ASP)
    
    secStiff = (fx_test_ASP{kk}(1:N_test)/pars.mu/N1)./(ux_test(1:N_test)/u_slip);
    
    plot(ux_test(1:N_test)/u_slip, (secStiff-secStiff_Mindlin)./secStiff_Mindlin, symbol_set{kk}, 'LineWidth', 2.5);
        
end

legend(legend_names(2:(sum(SPHERE_TESTS)+1)), 'Location', 'ne');


%% Normalized Comparison for IMAC paper

% set(0, 'defaultaxesinterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

addpath('../ROUTINES/')
color_plot = DISTINGUISHABLE_COLORS(5, 'w');

color_plot(4, :) = color_plot(end, :);

% SPHERE_TESTS(4) = false;

% Plot the comparison
figure; 
hold on;
xlabel('$u_t/u_{slip}$');
ylabel('$f_t/(\mu N)$');

xlim(hyst_amp*[-1.1, 1.1]);
ylim([-1.1, 1.1]);

legend_names = cell(length(fx_test_ASP)+1, 1);

plot(ux_test/u_slip, fx_test/pars.mu/N1, 'k-', 'LineWidth', 3);
legend_names{1} = 'Mindlin Sphere';

symbol_set = {'m--', 'r-.', 'b:','--', '-.', ':','--', '-.', ':', '--', '-.', ':'};
symbol_set = {':', '--', '-.','--', '-.', ':','--', '-.', ':', '--', '-.', ':'};

legend_ind = 2;

for kk = [2, 3, 1, 4] %include 4 only if want fitted Iwan model
    
    plot(ux_test/u_slip, fx_test_ASP{kk}/pars.mu/N1, symbol_set{kk}, 'LineWidth', 3, 'color', color_plot(kk, :));
    
    legend_names{legend_ind} = ASP_NAMES_LIST{kk};
    
    legend_ind = legend_ind + 1;
    
end

% legend(legend_names([1, find(SPHERE_TESTS)+1]), 'Location', 'nw');
legend('Analytical', 'Tangent', 'Secant', 'Iwan', 'Iwan Fit', 'Location', 'nw');



set(gcf, 'Renderer', 'painters');

% %get figure details
% fig = gcf;
% fig.Position = [250 250 500 500];
    
% print('AspHysteresisLoops.eps', '-depsc')

% set(0, 'defaultaxesinterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
set(groot, 'defaultLegendInterpreter','default');
set(groot, 'defaultTextInterpreter','default');


% print('../FJSIMS/Results/IMAC_Slides/ModelHyst_1.svg', '-dsvg')
% print('../FJSIMS/Results/IMAC_Slides/ModelHyst_2.svg', '-dsvg')

% SPHERE_TESTS(4) = true;


%% Ellipsoid Contact and Partial Slip Solution

hyst_amp = 1.5; 

% %% First Set of Hysteresis Loop
un = 1e-5;

[fxyn, uxy_slip] = ELLIPSOID_MONOTONIC([0, 0, un], pars);
u_slip = uxy_slip(1);
N1 = fxyn(1,3);  

ux_test = [linspace(0, hyst_amp*u_slip, N_test), ...
          linspace(hyst_amp*u_slip, -hyst_amp*u_slip, N_test), ...
          linspace(-hyst_amp*u_slip, hyst_amp*u_slip, N_test)];
      
un_test = un*ones(size(ux_test));

% Initialize storage

% Initialize force storage
fx_test = zeros(size(ux_test));
fx_test_ASP = cell(size(ASP_FUN_LIST));

tx0 = cell(size(ASP_FUN_LIST));
ty0 = cell(size(ASP_FUN_LIST));
rq = cell(size(ASP_FUN_LIST));
rq0 = cell(size(ASP_FUN_LIST));
wq = cell(size(ASP_FUN_LIST));
uxyn0_ASP = cell(size(ASP_FUN_LIST));

for ii = 1:length(ASP_FUN_LIST)
    fx_test_ASP{ii} = zeros(size(ux_test));
    
    tx0{ii} = zeros(1, ASP_NQP_RAD_LIST(ii));
    ty0{ii} = zeros(1, ASP_NQP_RAD_LIST(ii));
    rq{ii} = linspace(0, 1, ASP_NQP_RAD_LIST(ii)); 
    rq0{ii} = linspace(0, 1, ASP_NQP_RAD_LIST(ii)); 
    
    
    wq{ii} = ones(size(rq{ii}));
    wq{ii}(2:end-1) = 2;
    wq{ii} = wq{ii}/sum(wq{ii});
    
    uxyn0_ASP{ii} = [0,0,0];
end


%%%% Generate Ellipsoid Solution data w/Masing assumptions
[fxyn, uxy_slip] = ELLIPSOID_MONOTONIC([ux_test(1:N_test)', ux_test(1:N_test)', un_test(1:N_test)'], pars);

%monotonic loading
fx_test_ell(1:N_test) = fxyn(:, 1)';
fy_test_ell(1:N_test) = fxyn(:, 2)';

%hysteretic unloading
fx_test_ell((N_test+1):(2*N_test)) = fxyn(N_test, 1)' - 2*fxyn(:, 1)';
fy_test_ell((N_test+1):(2*N_test)) = fxyn(N_test, 2)' - 2*fxyn(:, 2)';

%hysteretic unloading
fx_test_ell((2*N_test+1):(3*N_test)) = -fxyn(N_test, 1)' + 2*fxyn(:, 1)';
fy_test_ell((2*N_test+1):(3*N_test)) = -fxyn(N_test, 2)' + 2*fxyn(:, 2)';



%%%%%%%%%%%%%%%%%%%%
% Now run asperity models

%Initialize to the first normal load for the ASP_FUN as well 
for kk = find(~SPHERE_TESTS) %1:length(ASP_FUN_LIST)

    uxyn_ASP = [0, 0, un_test(1)];

    [fxyn, dfxynduxyn, rq0{kk}, tx0{kk}, ty0{kk}] = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk});

    uxyn0_ASP{kk} = uxyn_ASP;

end



for jj = 1:length(ux_test)

    %Loop over asperity functions
    for kk = find(~SPHERE_TESTS) %1:length(ASP_FUN_LIST)

        xynum_asp = ASP_LIST_XYINDEX(kk);

        %Asperity function test. 
        uxyn_ASP = [0, 0, 0];
        uxyn_ASP(1,xynum_asp) = ux_test(jj);
        uxyn_ASP(1, 3)    = un_test(jj);

        [fxyn, dfxynduxyn, rq0{kk}, tx0{kk}, ty0{kk}] = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk});

        uxyn0_ASP{kk} = uxyn_ASP;

        fx_test_ASP{kk}(jj) = fxyn(xynum_asp);
%             jj

    end
end

figure;
hold on;
xlabel('Displacement');
ylabel('Force');

plot(ux_test, fx_test_ell, '-', 'LineWidth', 2.5)
plot(ux_test, fy_test_ell, '-', 'LineWidth', 2.5)

symbol_set = {'--', '-.', ':','--', '-.', ':','--', '-.', ':', '--', '-.', ':'};

legend_names = {'Analytical x', 'Analytical y'};
legend_ind = 2;

for kk = find(~SPHERE_TESTS) %1:length(fx_test_ASP)
    
    plot(ux_test, fx_test_ASP{kk}, symbol_set{kk}, 'LineWidth', 2.5);
    
    legend_names{legend_ind+1} = ASP_NAMES_LIST{kk};
    legend_ind = legend_ind + 1;
    
end

legend(legend_names);

%%%%%% Look at dissipation with the first N_test entries and Masing
%%%%%% assumptions

D_ell_x = cumtrapz(ux_test(1:N_test), fx_test_ell(1:N_test));
D_ell_x = 8*D_ell_x - 4 * ux_test(1:N_test).*fx_test_ell(1:N_test);

D_ell_y = cumtrapz(ux_test(1:N_test), fy_test_ell(1:N_test));
D_ell_y = 8*D_ell_y - 4 * ux_test(1:N_test).*fy_test_ell(1:N_test);

D_ASP = cell(size(fx_test_ASP));

for kk = find(~SPHERE_TESTS) %1:length(fx_test_ASP)
    
    D_ASP{kk} = cumtrapz(ux_test(1:N_test), fx_test_ASP{kk}(1:N_test));
    D_ASP{kk} = 8*D_ASP{kk} - 4 * ux_test(1:N_test).*fx_test_ASP{kk}(1:N_test);
    
end

figure; 
hold on;
xlabel('uxy');
ylabel('Dissipation');

plot(ux_test(1:N_test), D_ell_x, '-', 'LineWidth', 2.5);
plot(ux_test(1:N_test), D_ell_y, '-', 'LineWidth', 2.5);

for kk = find(~SPHERE_TESTS) %1:3%length(fx_test_ASP)
    
    plot(ux_test(1:N_test), D_ASP{kk}, symbol_set{kk}, 'LineWidth', 2.5);
        
end

legend(legend_names, 'Location', 'se');


% Error in dissipation calculation. 
figure; 
hold on;
xlabel('ux');
ylabel('Dissipation [non-dimensional]');


for kk = find(~SPHERE_TESTS) %1:3%length(fx_test_ASP)
    
    if(ASP_LIST_XYINDEX(kk) == 1)
        plot(ux_test(1:N_test), (D_ASP{kk}-D_ell_x)./D_ell_x, symbol_set{kk}, 'LineWidth', 2.5);
    else
        plot(ux_test(1:N_test), (D_ASP{kk}-D_ell_y)./D_ell_y, symbol_set{kk}, 'LineWidth', 2.5);
    end
end

legend(legend_names(3:end), 'Location', 'se');



% Error in Secant stiffness/Secant Stiffness
figure; 
hold on;
xlabel('uxy');
ylabel('Secant Stiffness');


plot(ux_test(1:N_test), fx_test_ell(1:N_test)./ux_test(1:N_test), '-', 'LineWidth', 2.5)
plot(ux_test(1:N_test), fy_test_ell(1:N_test)./ux_test(1:N_test), '-', 'LineWidth', 2.5)

for kk = find(~SPHERE_TESTS) %1:3%length(fx_test_ASP)
    
    secStiff = (fx_test_ASP{kk}(1:N_test))./(ux_test(1:N_test));
    
    plot(ux_test(1:N_test), secStiff, symbol_set{kk}, 'LineWidth', 2.5);
        
end

legend(legend_names, 'Location', 'ne');


% Error in Secant stiffness/Secant Stiffness
figure; 
hold on;
xlabel('uxy/u slip');
ylabel('Frac Error in Secant Stiffness');

secStiff_Ellx = (fx_test_ell(1:N_test))./(ux_test(1:N_test));
secStiff_Elly = (fy_test_ell(1:N_test))./(ux_test(1:N_test));

for kk = find(~SPHERE_TESTS) %1:3%length(fx_test_ASP)
    
    secStiff = (fx_test_ASP{kk}(1:N_test))./(ux_test(1:N_test));
        
    if(ASP_LIST_XYINDEX(kk) == 1)
        plot(ux_test(1:N_test), (secStiff-secStiff_Ellx)./secStiff_Ellx, symbol_set{kk}, 'LineWidth', 2.5);
    else
        plot(ux_test(1:N_test), (secStiff-secStiff_Elly)./secStiff_Elly, symbol_set{kk}, 'LineWidth', 2.5);
    end
        
end

legend(legend_names(3:end), 'Location', 'ne');

