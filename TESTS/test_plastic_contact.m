%% Tests of Plastic Contact Eqns

clear;

addpath('../ROUTINES/FRIC_MODELS')
addpath('../ROUTINES/FRIC_MODELS/ASP_FUN/')
addpath('../ROUTINES/FRIC_MODELS/MINDLIN/')
addpath('../ROUTINES/FRIC_MODELS/PLASTIC/')

addpath('../ROUTINES/')



%% Baseline Properties and Yield

%Initialize the sphere parameters
E = 1.9231e+11; %304L from Traction paper
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
pars.Rprime = 0.003166873052212; %(1/0.018453533867678 + 1/0.019130030937404)^(-1);
pars.Rpprime = 5.882905464446780e-04; %(1/0.000172616666847 + 1/0.000199388648556)^(-1);
pars = INITIALIZE_ECCENTRIC(pars);
pars.R = pars.Re;

%% Yield Parameters

Sys = 330e6; %Pa

C = 1.295*exp(0.736*pars.nu);

delta_y1s = (pi*C*Sys/(2*(2*pars.Estar)))^2*(2*pars.R); %displacement of one sphere against a rigid flat to cause yielding.

% Store all of these results as system parameters to prevent recalculation:
pars.delta_y = delta_y1s*2;
pars.Sys = Sys;

%tangent stiffness after yield.
pars.Et = 620e6;

pars.UTS = 713.5e6; % Pa


%% Functions

Nqp_radius = 1;

testvarargout = false;

% ASP_FUN = @SPHERE_PLASTIC_PRE;
ASP_FUN = @SPHERE_PL_TAN_DECOUPLE_PL_NORM; Nqp_radius = 1;

% Under development:
% ASP_FUN = @ELLIPSOID_PLASTIC_PRE;

pars.sliptype = 5; % (2=surface yielding) (3=CEB Model) (5=CEB based on UTS)

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
uxyn = [1e-6, 0.5e-6, 1e-7];
% uxyn = [200e-6, 100e-6, 2e-5]; %sqrt(sum(fxyn(1:2).^2))/fxyn(3)

deltam = 0;
Fm = 0;
am = 0;

if(testvarargout)
    [fxyn, dfxynduxyn, ~, ~, ~, ~, ~, ~, ~, ~, ~, a, dadun, Sys, dSysdun] = ASP_FUN(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
else
    [fxyn, dfxynduxyn, ~, ~, ~, ~, ~, ~] = ASP_FUN(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
end

%% Test With an Initial Load and then Retest with tangent

uxyn0 = [0, 0, 0];
uxyn1 = [0, 0, 2e-5]; %ONLY NORMAL IS NONZERO!
% uxyn2 = [0.01e-6, 0.01e-6, 1.9e-5];
% uxyn2 = [0.0e-6, 0.0e-6, 1.9e-5];
uxyn2 = [1.0e-6, 1.0e-6, 2e-5-5.7e-6];


% un_ref = [0.1115e-3, 0.8612e-4, 0.7920e-4]-15e-6; % Max displacements from a test
% uxyn1 = [0, 0, un_ref(1)]; %ONLY NORMAL IS NONZERO!
% uxyn2 = [0.0e-6, 0.0e-6, un_ref(1)];


deltam = 0;
Fm = 0;
am = 0;

if(testvarargout)
    [fxyn, dfxynduxyn, ~, ~, ~, deltam1, Fm1, am1, ~, ~, ~, a, dadun, Sys, dSysdun] = ASP_FUN(pars, uxyn1, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
else
    [fxyn, dfxynduxyn, ~, ~, ~, deltam1, Fm1, am1] = ASP_FUN(pars, uxyn1, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
end


if(testvarargout)
    [fxyn, dfxynduxyn, ~, ~, ~, deltam2, Fm2, am2, ~, ~, ~, a, dadun, Sys, dSysdun] = ASP_FUN(pars, uxyn2, uxyn0, rq0, tx0, ty0, rq, wq, deltam1, Fm1, am1);
else
    [fxyn, dfxynduxyn, ~, ~, ~, deltam2, Fm2, am2] = ASP_FUN(pars, uxyn2, uxyn0, rq0, tx0, ty0, rq, wq, deltam1, Fm1, am1);
end

% fxyn
dfxynduxyn

%% Derivative Test at a point (run for many values)


pars.mu = 0.05; %the default that most tests were written for.
pars.sliptype = 2; % Default to be used. 

% testvarargout = true;

%Initialize at 0, then go to 1, then test derivative at uxyn2

uxyn0 = [0, 0, 0];

% Tests 1-15 handle plasticity and Coulomb friction
% Tests 16-? Are for checking the CEB model
derivative_test = 42; %1-40+, start at 16 for CEB models (e.g. with Brake unloading)

[uxyn1, uxyn2, pars] = TEST_DISPLACEMENTS(derivative_test, pars);

%Tests with repeated normal displacements in the yield regime are expected
%to fail. - E.G. Case 9 is expected to fail on dfxynuxyn(3,3)



deltam0 = 0;
Fm0 = 0;
am0 = 0;

if(uxyn1(3) > 0)
    %Initialize
    [~, ~, rq1, tx1, ty1, deltam1, Fm1, am1] = ASP_FUN(pars, uxyn1, uxyn0, rq0, tx0, ty0, rq, wq, deltam0, Fm0, am0);
else
    rq1 = rq;
    tx1 = zeros(size(rq));
    ty1 = zeros(size(rq));
    
    
    deltam1 = 0;
    Fm1 = 0;
    am1 = 0;

end


%Baseline/analytical solution
if(testvarargout)
    [fxyn, dfxynduxyn, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, dadun, ~, dSysdun] = ASP_FUN(pars, uxyn2, uxyn1, rq1, tx1, ty1, rq, wq, deltam1, Fm1, am1); %extra derivative check
else
    [fxyn, dfxynduxyn, ~, ~, ~, ~, ~, ~] = ASP_FUN(pars, uxyn2, uxyn1, rq1, tx1, ty1, rq, wq, deltam1, Fm1, am1); %w/o extra check
end


% [fxyn, dfxynduxyn, ~, ~, ~, ~, duxyn1duxyn2] = ASP_FUN(pars, uxyn2, uxyn1, rq1, tx1, ty1, rq, wq);
% num_duxynprev = zeros(size(duxyn1duxyn2));



% delta = 1e-8; %displacement variation (m)
delta = [1e-8, 1e-9, 1e-10, 1e-12]; %displacement variation (m)

errorK = zeros(size(delta));
errorA = zeros(size(delta));
errorSys = zeros(size(delta));

numK = zeros(size(dfxynduxyn));

for ii = 1:length(delta)
    
    for jj = 1:3
        
        uxyn_low = uxyn2;
        uxyn_high = uxyn2;
        
        uxyn_low(jj) = uxyn_low(jj) - delta(ii);
        uxyn_high(jj) = uxyn_high(jj) + delta(ii);
        
        if(testvarargout)
            [fxyn_high, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, a_high, ~, Sys_high, ~] = ASP_FUN(pars, uxyn_high, uxyn1, rq1, tx1, ty1, rq, wq, deltam1, Fm1, am1);
            [fxyn_low, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, a_low, ~, Sys_low, ~] = ASP_FUN(pars, uxyn_low, uxyn1, rq1, tx1, ty1, rq, wq, deltam1, Fm1, am1);
            
            if(jj == 3)
                dadun_num = (a_high - a_low)/(2*delta(ii));
                dSysdun_num = (Sys_high - Sys_low)/(2*delta(ii));
            end
            
        else
            [fxyn_high, ~, ~, ~, ~, ~] = ASP_FUN(pars, uxyn_high, uxyn1, rq1, tx1, ty1, rq, wq, deltam1, Fm1, am1);
            [fxyn_low, ~, ~, ~, ~, ~] = ASP_FUN(pars, uxyn_low, uxyn1, rq1, tx1, ty1, rq, wq, deltam1, Fm1, am1);
        end
        
        numK(:, jj) = (fxyn_high' - fxyn_low')/(2*delta(ii));
        
    end
    
%     num_duxynprev
    
    errorK(ii) = norm(dfxynduxyn - numK)/(norm(dfxynduxyn) + (norm(dfxynduxyn) == 0) );
    
    if(testvarargout)        
        errorA(ii) = abs(dadun_num - dadun) / abs(dadun);
        errorSys(ii) = abs(dSysdun_num - dSysdun)/abs(dSysdun);
    end
    
end

if(testvarargout)
    errorA
    errorSys %Probably always will be inf since I always chose the unloading direction for the derivative.
end

errorK

%% Normal Loading Curve verification - Near Yield

% THIS BLOCK IS TO VERIFY THAT IT IS CONTINUOUS AROUND YIELD. DO NOT
% INCREASE unMax.
% If you want to plot normal force for larger range use the next block!!!!

unMax = 4*pars.delta_y;

un = [linspace(0, unMax, 400), linspace(unMax, 0, 100)];

fn_ASP = zeros(size(un));
ahist = zeros(size(un));
Syshist = zeros(size(un));

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


deltam = 0;
Fm = 0;
am = 0;


for ii = 1:length(un)

    uxyn_tmp(3) = un(ii);
    
%     if(ii > 405)
%         3
%     end
    
    %Asperity function
    if(testvarargout)
        [fxyn, ~, ~, ~, ~, deltam, Fm, am, ~, ~, ~, ahist(ii), ~, Syshist(ii), ~] = ASP_FUN(pars, uxyn_tmp, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
    else
        [fxyn, ~, ~, ~, ~, deltam, Fm, am] = ASP_FUN(pars, uxyn_tmp, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
    end
    
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

if(testvarargout)
    
    figure;
    hold on;
    xlabel('un');
    ylabel('Contact Radius a');

    plot(un, ahist, 'LineWidth', 2);

    yyaxis right
    plot(un, Syshist, 'LineWidth', 2);
    hold on;
    ylabel('Sys');

    legend('Contact Radius', 'Yield Stress');
end


%% Normal Loading Curve verification - Whatever desired amplitude

% unMax = 1.5e-5; %generally what I have been plotting.
% unMax = .5e-5; %generally what I have been plotting.
% unMax = 3e-5; %
unMax = .8e-6; %Clear elastic + plastic transition - what is used in paper


loadingPoints = 400;
un = [linspace(0, unMax, loadingPoints), linspace(unMax, 0, 100)];

fn_ASP = zeros(size(un));
ahist = zeros(size(un));
Syshist = zeros(size(un));

uxyn0 = [0, 0, 0];

prev0.uxnyn  = [0, 0, 0, 0]; %From the previous step
prev0.fxnyn  = [0, 0, 0, 0];
prev0.uxnyn0 = [0, 0, 0, 0]; %From previous reversal point
prev0.fxnyn0 = [0, 0, 0, 0];
prev0.uxynOrigin = [0, 0, 0]; %Origin for where to consider displacement. 


%Mindlin normal force
fxyn_MINDLIN = MINDLIN_MONOTONIC([zeros(length(un), 2), un'], pars);
fn_MINDLIN = fxyn_MINDLIN(:, 3)';

% % Elliptical contact normal force
[fxyn_ELLIPSOID, uxy_slip] = ELLIPSOID_MONOTONIC([zeros(length(un), 2), un'], pars);
fn_ELLIPSOID = fxyn_ELLIPSOID(:, 3)';


deltam = 0;
Fm = 0;
am = 0;


for ii = 1:length(un)

    uxyn_tmp(3) = un(ii);
    
    %Asperity function
    if(testvarargout)
        [fxyn, ~, ~, ~, ~, deltam, Fm, am, ~, ~, ~, ahist(ii), ~, Syshist(ii), ~] = ASP_FUN(pars, uxyn_tmp, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
    else
        [fxyn, ~, ~, ~, ~, deltam, Fm, am] = ASP_FUN(pars, uxyn_tmp, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
    end
    
    fn_ASP(ii) = fxyn(3);

end

close all;




colors_list = DISTINGUISHABLE_COLORS(13, 'w');
colors_list = colors_list([1:3, 5:end], :);
% set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
set(groot, 'defaultAxesTickLabelInterpreter','latex');  %Tex
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');


figure;
hold on;
xlabel('Normal Displacement [m]');
ylabel('Normal Force [N]');

linewidth = 3;

plot(un, fn_MINDLIN, '-', 'LineWidth', linewidth, 'Color', colors_list(1, :));
plot(un(1:loadingPoints), fn_ELLIPSOID(1:loadingPoints), '-.', 'LineWidth', linewidth, 'Color', colors_list(3, :));
plot(un(1:loadingPoints), fn_ASP(1:loadingPoints), '--', 'LineWidth', linewidth, 'Color', colors_list(2, :));
plot(un(loadingPoints+1:end), fn_ASP(loadingPoints+1:end), ':', 'LineWidth', linewidth, 'Color', colors_list(6, :));

% legend('Hertz Sphere', 'Hertz Ellipsoid', 'Asperity');
legend('Fully Elastic Sphere', 'Fully Elastic Ellipsoid', 'Elastic-Plastic Loading (Sphere)', 'Elastic Unloading (Sphere)', 'Location', 'nw');


set(gca,'FontSize',14)

%%%%%%%%%%%%
% Manually add some tick information:

% Find deltabar (approximately)
f_unload = fn_ASP(loadingPoints+1:end);
u_unload = un(loadingPoints+1:end);

f_unload = f_unload(1:find(f_unload == 0, 1));
u_unload = u_unload(1:length(f_unload));
deltabar = interp1(f_unload, u_unload, 0);

ax = gca;

labels = [ax.XTickLabel; {'$\delta_y$'}; {'$\begin{array}{c} \\ 1.9\delta_y\end{array}$'}; {'$\bar{\delta}$'}];
vals = [ax.XTick, pars.delta_y, 1.9*pars.delta_y, deltabar];

[vals, inds] = sort(vals);
labels = labels(inds);

ax.XTickLabel = labels;
ax.XTick = vals;

% ax.XTickLabelRotation = 60;

%%%%%%%%%%%%
% Ytick information
[fxyn_y] = ASP_FUN(pars, [0, 0, pars.delta_y], [0, 0, 0], rq0, tx0, ty0, rq, wq, 0, 0, 0);
[fxyn_19y] = ASP_FUN(pars, [0, 0, 1.9*pars.delta_y], [0, 0, 0], rq0, tx0, ty0, rq, wq, 0, 0, 0);

labels = [ax.YTickLabel; {'$f_{n,\delta y} \ \ $'}; {'$f_{n, 1.9\delta y}$'}];
vals = [ax.YTick, fxyn_y(3), fxyn_19y(3)];

[vals, inds] = sort(vals);
labels = labels(inds);

ax.YTickLabel = labels;
ax.YTick = vals;

ax.TickLength = [0.02, 0.05];

%%%%%%%%%%%%


set(gcf, 'Renderer', 'painters');
drawnow;

% print('../TMD2021_EXTABS/Figures/normal_asperity_force.eps', '-depsc', '-r400');
% print('../TMD2021_EXTABS/Figures/normal_asperity_force.svg', '-dsvg', '-r400');
% print('../PLOTS/OUTPUTS/normal_asperity_force.eps', '-depsc', '-r400');


if(testvarargout)
    
    figure;
    hold on;
    xlabel('un');
    ylabel('Contact Radius a');

    plot(un, ahist, 'LineWidth', 2);

    yyaxis right
    plot(un, Syshist, 'LineWidth', 2);
    hold on;
    ylabel('Sys');

    legend('Contact Radius', 'Yield Stress');
end


