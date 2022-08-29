% This script tests the harmonic derivatives of an asperity function for
% use with EPMC. Written specifically for the plastic asperity models, so
% may require adjustments if using with the elastic asperity models

clear;

addpath('../')
addpath('../../ROUTINES/FRIC_MODELS')
addpath('../../ROUTINES/FRIC_MODELS/ASP_FUN/')
addpath('../../ROUTINES/FRIC_MODELS/MINDLIN/')
addpath('../../ROUTINES/FRIC_MODELS/PLASTIC/')


%% Baseline Properties and Yield

%Initialize the sphere parameters
E = 192.85e9; %304L from Traction paper
nu = 0.30; %304L taken from Traction paper
G = E / 2/ (1 + nu);%304L approx Shear modulus

pars.mu     = 0.05; %may want to do a function of displacement - see Eriten's work
% pars.R      = (1/Rx + 1/Ry)^(-1);%, 1e-4 is about a maximum
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

%% Yield Parameters

Sys = 200e6; %Pa

C = 1.295*exp(0.736*pars.nu);

delta_y1s = (pi*C*Sys/(2*(2*pars.Estar)))^2*(2*pars.R); %displacement of one sphere against a rigid flat to cause yielding.

% Store all of these results as system parameters to prevent recalculation:
pars.delta_y = delta_y1s*2;
pars.Sys = Sys;

%tangent stiffness after yield.
pars.Et = 0.015*pars.E;


%% Functions

Nqp_radius = 1;

testvarargout = false;

%%%%%%%%%%%%%%%%%%%
% Brake Unloading
% ASP_FUN = @SPHERE_PLASTIC_PRE;
ASP_FUN = @SPHERE_PL_TAN_DECOUPLE_PL_NORM; Nqp_radius = 1;

ASP_FUN_PRE = @SPHERE_PLASTIC_PRE;

pars.sliptype = 3; % Setting divisible by 2 will consider plasticity and Coulomb., Divisible by 3 will do min(others, CEB)

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
uxyn = [1e-6, 0.5e-6, 1e-5];
% uxyn = [200e-6, 100e-6, 2e-5]; %sqrt(sum(fxyn(1:2).^2))/fxyn(3)

deltam = 0;
Fm = 0;
am = 0;

if(testvarargout)
    [fxyn, dfxynduxyn, ~, ~, ~, ~, ~, ~, a, dadun, Sys, dSysdun] = ASP_FUN(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
else
    [fxyn, dfxynduxyn, ~, ~, ~, ~, ~, ~] = ASP_FUN(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
end

%% Test With an Initial Load and then Retest with tangent

uxyn0 = [0, 0, 0];
uxyn1 = [0, 0, 2e-5]; %ONLY NORMAL IS NONZERO!
% uxyn2 = [0.01e-6, 0.01e-6, 1.9e-5];
uxyn2 = [0.0e-6, 0.0e-6, 1.9e-5];

deltam = 0;
Fm = 0;
am = 0;

if(testvarargout)
    [fxyn, dfxynduxyn, ~, ~, ~, deltam1, Fm1, am1, a, dadun, Sys, dSysdun] = ASP_FUN(pars, uxyn1, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
else
    [fxyn, dfxynduxyn, ~, ~, ~, deltam1, Fm1, am1] = ASP_FUN(pars, uxyn1, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
end


if(testvarargout)
    [fxyn, dfxynduxyn, ~, ~, ~, deltam2, Fm2, am2, a, dadun, Sys, dSysdun] = ASP_FUN(pars, uxyn2, uxyn0, rq0, tx0, ty0, rq, wq, deltam1, Fm1, am1);
else
    [fxyn, dfxynduxyn, ~, ~, ~, deltam2, Fm2, am2] = ASP_FUN(pars, uxyn2, uxyn0, rq0, tx0, ty0, rq, wq, deltam1, Fm1, am1);
end

% fxyn
dfxynduxyn


%% Test Derivatives

uxyn0 = [0, 0, 0];

derivative_test = 16; %1-40+?


% These tests need to be verified that they are reasonable, but this is a
% good list of tests for general asperities
[uxyn2, uxyn3, pars] = TEST_DISPLACEMENTS(derivative_test, pars);

% Initialization
[fxyn2, fxyn3, dfxyn3duxyn3, dfxyn3duxyn2, dfxyn3dfxy2] = calc_force_test(ASP_FUN_PRE, ASP_FUN, pars, uxyn2, uxyn3, rq0, tx0, ty0, rq, wq, testvarargout);


normalize_error2 = norm(dfxyn3duxyn2)~=0; 
normalize_error3 = norm(dfxyn3duxyn3)~=0; 

% If want to manually set normalization:
% normalize_error2 = false; %False is generally desired when the denominator may be zero frequently.
% normalize_error3 = true; %False is generally desired when the denominator may be zero frequently.




delta = [1e-8, 1e-9, 1e-10, 1e-12]; %displacement variation (m)

error_dfxyn3duxyn3 = zeros(size(delta));
error_dfxyn3duxyn2 = zeros(size(delta));

for ii = 1:length(delta)
    
    num_dfxyn3duxyn3 = zeros(size(dfxyn3duxyn3));
    num_dfxyn3duxyn2 = zeros(size(dfxyn3duxyn2));
    
    for dim = 1:3
        
        uxyn2_low  = uxyn2;
        uxyn2_high = uxyn2;
        
        uxyn3_low  = uxyn3;
        uxyn3_high = uxyn3;
        
        uxyn2_low(dim)  = uxyn2_low(dim)  - delta(ii);
        uxyn2_high(dim) = uxyn2_high(dim) + delta(ii);
        
        uxyn3_low(dim)  = uxyn3_low(dim)  - delta(ii);
        uxyn3_high(dim) = uxyn3_high(dim) + delta(ii);
        
        % Calculate perturbed forces
        [~, fxyn3_low2,  ~, ~, ~] = calc_force_test(ASP_FUN_PRE, ASP_FUN, pars, uxyn2_low, uxyn3, rq0, tx0, ty0, rq, wq, testvarargout);
        [~, fxyn3_high2, ~, ~, ~] = calc_force_test(ASP_FUN_PRE, ASP_FUN, pars, uxyn2_high, uxyn3, rq0, tx0, ty0, rq, wq, testvarargout);
        
        [~, fxyn3_low3,  ~, ~, ~] = calc_force_test(ASP_FUN_PRE, ASP_FUN, pars, uxyn2, uxyn3_low, rq0, tx0, ty0, rq, wq, testvarargout);
        [~, fxyn3_high3, ~, ~, ~] = calc_force_test(ASP_FUN_PRE, ASP_FUN, pars, uxyn2, uxyn3_high, rq0, tx0, ty0, rq, wq, testvarargout);
        
        % Stiffness matrix
        num_dfxyn3duxyn3(:, dim) = (fxyn3_high3' - fxyn3_low3')/(2*delta(ii));
        num_dfxyn3duxyn2(:, dim) = (fxyn3_high2' - fxyn3_low2')/(2*delta(ii));
        
    end
    
    error_dfxyn3duxyn3(ii) = norm(dfxyn3duxyn3 - num_dfxyn3duxyn3) / (norm(dfxyn3duxyn3).*normalize_error3 + ~normalize_error3);
    error_dfxyn3duxyn2(ii) = norm(dfxyn3duxyn2 - num_dfxyn3duxyn2) / (norm(dfxyn3duxyn2).*normalize_error2 + ~normalize_error2);
    
end

error_dfxyn3duxyn2
error_dfxyn3duxyn3



%% Explicitly perturb force / all variables and only evaluate at uxyn3

uxyn0 = [0, 0, 0];
deltam0 = 0;
Fm0 = 0;
am0 = 0;


derivative_test = 16; %1-40+?

% WARNING: the derivative w.r.t. uxyn2(3) uses a central difference derivative
% which can sometimes allow for inputs of uxyn2(3) > deltam. This is ignored
% since deltam is treated as independent of uxyn2(3) and thus derivatives
% are not related. This could lead to this test being wrong for some cases
% not yet considered.

% These errors have to do with deltam = uxyn3(3) so overlapped meanings: I
% believe they work out okay.
% For Tests: 24, 28, 36, 40: errors_dfddeltam - shows errors, but I think
% may be okay.
% 24, 28, 36, 40: errors_dfxyn3duxyn3 - shows errors


% These tests need to be verified that they are reasonable, but this is a
% good list of tests for general asperities
[uxyn2, uxyn3, pars] = TEST_DISPLACEMENTS(derivative_test, pars);


normalize_deltam = true; %False is generally desired when the denominator may be zero frequently.
normalize_fxy2 = false;
normalize_uxyn2 = false;
normalize_uxyn3 = true;

% Get the initial state
[~, ~, ~, ~, ~, inputs3] = calc_force_test(ASP_FUN_PRE, ASP_FUN, pars, uxyn2, uxyn3, rq0, tx0, ty0, rq, wq, testvarargout);

% Apply uxyn3
if(testvarargout)
    [fxyn3, dfxyn3duxyn3, rq3, tx3, ty3, deltam3, Fm3, am3, dfxyn3duxyn2, ...
        dfxyn3dfxy2, dfddeltam, a, dadun, Sys, dSysdun, daddeltam, dSysddeltam] = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
        inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
        inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam, testvarargout);
else
    [fxyn3, dfxyn3duxyn3, rq3, tx3, ty3, deltam3, Fm3, am3, dfxyn3duxyn2, ...
        dfxyn3dfxy2, dfddeltam] = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
        inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
        inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam, testvarargout);
    
end


delta = [1e-8, 1e-9, 1e-10, 1e-12]; %displacement variation (m)
deltaF = 100*delta/1e-5; 

errors_dfxyn3duxyn3 = zeros(size(delta));
errors_dfxyn3duxyn2 = zeros(size(delta));
errors_dfxyn3dfxy2 = zeros(size(delta));
errors_dfddeltam = zeros(size(delta));

errors_daddeltam = zeros(size(delta));
errors_dSysddeltam = zeros(size(delta));

errors_dadun = zeros(size(delta));
errors_dSysdun = zeros(size(delta));

% Errors for checking when some derivatives are ill defined
errors_dfxyn3duxyn3_nodun = zeros(size(delta));
errors_dfxyn3duxyn3_dun = zeros(size(delta));
errors_dfnddeltam = zeros(size(delta));


for ii = 1:length(delta)
    
    
    % Perturb deltam (also requires pertubing the other am/Fm quantities a
    % corresponding amount.
    % Assymetric derivative so deltam is always deltam and uxyn3 never
    % becomes deltam
    deltam_low  = inputs3.deltam2;
    deltam_high = inputs3.deltam2 + delta(ii);
    
    [~, ~, ~, ~, ~, deltam1_low, Fm1_low, am1_low] ...
        = ASP_FUN_PRE(pars, [0,0,deltam_low], uxyn0, rq0, tx0, ty0, rq, wq, deltam0, Fm0, am0);
    
    [~, ~, ~, ~, ~, deltam1_high, Fm1_high, am1_high] ...
        = ASP_FUN_PRE(pars, [0,0,deltam_high], uxyn0, rq0, tx0, ty0, rq, wq, deltam0, Fm0, am0);
    
    if(testvarargout)
        [fxyn3_low_deltam, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, a_low, ~, Sys_low, ...
            ~, ~, ~] = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
            inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
            deltam_low, Fm1_low, am1_low, inputs3.dFmddeltam, inputs3.damddeltam);

        [fxyn3_high_deltam, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, a_high, ~, Sys_high, ...
            ~, ~, ~] = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
            inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
            deltam_high, Fm1_high, am1_high, inputs3.dFmddeltam, inputs3.damddeltam);
    else
        [fxyn3_low_deltam] = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
            inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
            deltam_low, Fm1_low, am1_low, inputs3.dFmddeltam, inputs3.damddeltam);

        [fxyn3_high_deltam] = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
            inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
            deltam_high, Fm1_high, am1_high, inputs3.dFmddeltam, inputs3.damddeltam);
    end
    
    dfddeltam_num = (fxyn3_high_deltam' - fxyn3_low_deltam')/(delta(ii));
    
    errors_dfddeltam(ii) = norm(dfddeltam_num - dfddeltam) / (norm(dfddeltam) + (norm(dfddeltam)==0));
    errors_dfnddeltam(ii) = abs(dfddeltam_num(3) - dfddeltam(3)) / (abs(dfddeltam(3)) + (dfddeltam(3)==0));
    
    if(testvarargout)
        daddeltam_num = (a_high - a_low)/(delta(ii));
        dSysddeltam_num = (Sys_high - Sys_low)/(delta(ii));
        
        errors_daddeltam(ii) = abs(daddeltam - daddeltam_num)/abs(daddeltam+(daddeltam == 0));
        errors_dSysddeltam(ii) = abs(dSysddeltam - dSysddeltam_num)/abs(dSysddeltam + (dSysddeltam == 0));
    end
    
    % Perturb fx2
    fx2_high = inputs3.tx2 + deltaF(ii);
    fx2_low  = inputs3.tx2 - deltaF(ii);

    [fxyn3_high_fx2] = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
        inputs3.rq2, fx2_high, inputs3.ty2, inputs3.rq, inputs3.wq, ...
        inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);

    [fxyn3_low_fx2]  = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
        inputs3.rq2, fx2_low, inputs3.ty2, inputs3.rq, inputs3.wq, ...
        inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);
    
    dfxyn3dfxy2_num = zeros(3,2);
    dfxyn3dfxy2_num(:,1) = (fxyn3_high_fx2' - fxyn3_low_fx2')/(2*deltaF(ii));
    
    % Perturb fy0
    fy2_high = inputs3.ty2 + deltaF(ii);
    fy2_low  = inputs3.ty2 - deltaF(ii);
    
    [fxyn3_high_fy2] = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
        inputs3.rq2, inputs3.tx2, fy2_high, inputs3.rq, inputs3.wq, ...
        inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);

    [fxyn3_low_fy2]  = ASP_FUN(pars, inputs3.uxyn3, inputs3.uxyn2, ...
        inputs3.rq2, inputs3.tx2, fy2_low, inputs3.rq, inputs3.wq, ...
        inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);
    
    dfxyn3dfxy2_num(:,2) = (fxyn3_high_fy2' - fxyn3_low_fy2')/(2*deltaF(ii));
    
    errors_dfxyn3dfxy2(ii) = norm(dfxyn3dfxy2_num - dfxyn3dfxy2) / (norm(dfxyn3dfxy2).*normalize_fxy2 + ~normalize_fxy2);
    
    % Perturb uxyn2 / uxyn3
    dfxyn3duxyn2_num = zeros(3,3);
    dfxyn3duxyn3_num = zeros(3,3);
    
    for dim = 1:3
        
        uxyn2_low  = uxyn2;
        uxyn2_high = uxyn2;
        
        uxyn3_low  = uxyn3;
        uxyn3_high = uxyn3;
        
        uxyn2_low(dim)  = uxyn2_low(dim) - delta(ii);
        uxyn2_high(dim) = uxyn2_high(dim) + delta(ii);
        
        % Assymetric derivative so that uxyn3(3) is always <= deltam
        % This also causes other problems since Sys has it's derivative
        % asigned towards loading while dFdun has its derivative assigned
        % towards unloading with uxyn(3) = deltam
        uxyn3_low(dim)  = uxyn3_low(dim)-delta(ii);
        uxyn3_high(dim) = uxyn3_high(dim);
        
        [fxyn3_low2]  = ASP_FUN(pars, inputs3.uxyn3, uxyn2_low, ...
            inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
            inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);
        
        [fxyn3_high2] = ASP_FUN(pars, inputs3.uxyn3, uxyn2_high, ...
            inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
            inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);
        
        if(testvarargout)
            [fxyn3_low3, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, a_low, ~, Sys_low, ...
                ~, ~, ~]  = ASP_FUN(pars, uxyn3_low, inputs3.uxyn2, ...
                inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
                inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);

            [fxyn3_high3, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, a_high, ~, Sys_high, ...
                ~, ~, ~] = ASP_FUN(pars, uxyn3_high, inputs3.uxyn2, ...
                inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
                inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);
        else
            [fxyn3_low3]  = ASP_FUN(pars, uxyn3_low, inputs3.uxyn2, ...
                inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
                inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);

            [fxyn3_high3] = ASP_FUN(pars, uxyn3_high, inputs3.uxyn2, ...
                inputs3.rq2, inputs3.tx2, inputs3.ty2, inputs3.rq, inputs3.wq, ...
                inputs3.deltam2, inputs3.Fm2, inputs3.am2, inputs3.dFmddeltam, inputs3.damddeltam);
        end
        
        dfxyn3duxyn2_num(:, dim) = (fxyn3_high2' - fxyn3_low2')/(2*delta(ii));
        dfxyn3duxyn3_num(:, dim) = (fxyn3_high3' - fxyn3_low3')/(delta(ii));
        
        if(testvarargout)
            dadun_num = (a_high - a_low)/(delta(ii));
            dSysdun_num = (Sys_high - Sys_low)/(delta(ii));
            
            errors_dadun(ii) = abs(dadun_num - dadun)/(abs(dadun) + (dadun == 0 ));
            errors_dSysdun(ii) = abs(dSysdun_num - dSysdun)/(abs(dSysdun) + (dSysdun == 0 ));
            
        end

    end
    
    errors_dfxyn3duxyn2(ii) = norm(dfxyn3duxyn2_num - dfxyn3duxyn2) / (norm(dfxyn3duxyn2).*normalize_uxyn2 + ~normalize_uxyn2);
    errors_dfxyn3duxyn3(ii) = norm(dfxyn3duxyn3_num - dfxyn3duxyn3) / (norm(dfxyn3duxyn3).*normalize_uxyn3 + ~normalize_uxyn3);
    
    dfxyn3duxyn3_num_nodun = dfxyn3duxyn3_num.*[1, 1, 0; 1, 1, 0; 1, 1, 1];
    dfxyn3duxyn3_nodun = dfxyn3duxyn3.*[1, 1, 0; 1, 1, 0; 1, 1, 1];
    
    
    errors_dfxyn3duxyn3_nodun(ii) = norm(dfxyn3duxyn3_num_nodun - dfxyn3duxyn3_nodun) / (norm(dfxyn3duxyn3_nodun) + (norm(dfxyn3duxyn3_nodun)==0));
    errors_dfxyn3duxyn3_dun(ii) = norm(dfxyn3duxyn3_num(1:2, 3) - dfxyn3duxyn3(1:2, 3)) / (norm(dfxyn3duxyn3(1:2, 3)) + (norm(dfxyn3duxyn3(1:2, 3))==0));
    
end


disp('Go through all errors.');
errors_dfddeltam
errors_dfnddeltam
errors_dfxyn3dfxy2
errors_dfxyn3duxyn2
errors_dfxyn3duxyn3

errors_dfxyn3duxyn3_nodun
% errors_dfxyn3duxyn3_dun

if(testvarargout)
    errors_daddeltam
    errors_dSysddeltam
    errors_dadun
    errors_dSysdun
end


%%
% 
% a = linspace(0, am3, 1000);
% 
% Sys_a = Fm3./(pi.*a.^2).*(1 - (1 - (a./am3).^2).^1.5);
% 
% figure; 
% hold on;
% xlabel('a');
% ylabel('Sys');
% 
% plot(a, Sys_a);
% 
% (Sys_a(end-1)-Sys_a(end))/(a(end-1)-a(end))



%% Function to Calculate Final State after 3 loading steps
function [fxyn2, fxyn3, dfxyn3duxyn3, dfxyn3duxyn2, dfxyn3dfxy2, inputs3, extraDerivatives] = calc_force_test(ASP_FUN_PRE, ASP_FUN, pars, uxyn2, uxyn3, rq0, tx0, ty0, rq, wq, testvarargout)

uxyn0 = [0, 0, 0];

deltam0 = 0;
Fm0 = 0;
am0 = 0;


uxyn1 = [0, 0, max(uxyn2(3), uxyn3(3))];

% Apply uxyn1

% if(uxyn1(3) > 0)
    %Initialize
    [fxyn, dfxynduxyn, rq1, tx1, ty1, deltam1, Fm1, am1, dfxyn1duxyn0, dfxyn1dfxy0, ...
        dfddeltam, a, dadun, Sys, dSysdun, daddeltam, dSysddeltam] ...
        = ASP_FUN_PRE(pars, uxyn1, uxyn0, rq0, tx0, ty0, rq, wq, deltam0, Fm0, am0);
    
    % Since this state is by construction the loading to deltam, set these
    % variables here.
    dFmddeltam = dfxynduxyn(3,3);
    damddeltam = dadun;
% else
%     rq1 = rq;
%     tx1 = zeros(size(rq));
%     ty1 = zeros(size(rq));
%     
%     deltam1 = 0;
%     Fm1 = 0;
%     am1 = 0;
% end

% Apply uxyn2
[fxyn2, dfxyn2duxyn2, rq2, tx2, ty2, deltam2, Fm2, am2] ...
    = ASP_FUN(pars, uxyn2, uxyn1, rq1, tx1, ty1, rq, wq, deltam1, Fm1, am1, dFmddeltam, damddeltam);

if(~testvarargout)
    % Apply uxyn3
    [fxyn3, dfxyn3duxyn3, rq3, tx3, ty3, deltam3, Fm3, am3, dfxyn3duxyn2, ...
        dfxyn3dfxy2, dfddeltam] = ASP_FUN(pars, uxyn3, uxyn2, rq2, tx2, ty2, rq, wq, deltam2, Fm2, am2, dFmddeltam, damddeltam);
    
    extraDerivatives = 0;
    
else
    [fxyn3, dfxyn3duxyn3, rq3, tx3, ty3, deltam3, Fm3, am3, dfxyn3duxyn2, ...
        dfxyn3dfxy2, dfddeltam, a, dadun, Sys, dSysdun, daddeltam, dSysddeltam] = ASP_FUN(pars, uxyn3, uxyn2, rq2, tx2, ty2, rq, wq, deltam2, Fm2, am2, dFmddeltam, damddeltam);
    
    extraDerivatives.dadun = dadun;
    extraDerivatives.dSysdun = dSysdun;
    extraDerivatives.daddeltam = daddeltam;
    extraDerivatives.dSysddeltam = dSysddeltam;
end

inputs3.uxyn3 = uxyn3;
inputs3.uxyn2 = uxyn2;
inputs3.rq2 = rq2;
inputs3.tx2 = tx2;
inputs3.ty2 = ty2;
inputs3.rq = rq;
inputs3.wq = wq;
inputs3.deltam2 = deltam2;
inputs3.Fm2 = Fm2;
inputs3.am2 = am2;
inputs3.dFmddeltam = dFmddeltam;
inputs3.damddeltam = damddeltam;


% Collect and return derivatives - These should be double checked.

%d(currentForce)/d(uxyn3) including contributions from deltam
dfxyn3duxyn3 = dfxyn3duxyn3 + dfddeltam*[0, 0, 1]*(uxyn1(3) == uxyn3(3));

%d(currentForce)/d(uxyn2) including contributions from deltam
dfxyn3duxyn2 = dfxyn3duxyn2 + dfddeltam*[0, 0, 1]*(uxyn1(3) == uxyn2(3)) + dfxyn3dfxy2*dfxyn2duxyn2(1:2, :);

end



