%% Tests of Plastic Contact Eqns

clear;

addpath('../ROUTINES/FRIC_MODELS')
addpath('../ROUTINES/FRIC_MODELS/ASP_FUN/')
addpath('../ROUTINES/FRIC_MODELS/MINDLIN/')
addpath('../ROUTINES/FRIC_MODELS/PLASTIC/')


%% Baseline Properties and Yield

%Initialize the sphere parameters
E = 192.85e9; %304L from Traction paper
nu = 0.29; %304L taken from Traction paper
G = E / 2/ (1 + nu);%304L approx Shear modulus

pars.mu     = 0.05; %may want to do a function of displacement - see Eriten's work
% pars.R      = (1/Rx + 1/Ry)^(-1);%, 1e-4 is about a maximum
pars.E      = E;
pars.Estar  = E/2/(1 - nu^2);
pars.nu     = nu;
pars.G      = G;
pars.Gstar  = G/2/(2 - nu);

%Ellipsoid initialization, comment out if Rprime=Rpprime.
pars.Rprime = (1/0.018453533867678 + 1/0.019130030937404)^(-1);
pars.Rpprime = (1/0.000172616666847 + 1/0.000199388648556)^(-1);
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

testvarargout = true;

ASP_FUN = @SPHERE_PLASTIC_PRE; pars.sliptype = 3;


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
uxyn = [0.5*2e-6, 0.5e-6, 1e-5];
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
uxyn2 = [1e-6, 0.5e-6, 1.9e-5];

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

fxyn
dfxynduxyn


%% Normal Loading Curve verification

% un = [linspace(0, 3e-5, 100), linspace(3e-5, 0, 100)];
unLoad   = linspace(0, 1.5e-5, 100);
unUnload = linspace(1.5e-5, 0, 100);


fn_ASP = zeros(size(unLoad));
ahist = zeros(size(unLoad));
Syshist = zeros(size(unLoad));

uxyn0 = [0, 0, 0];

prev0.uxnyn  = [0, 0, 0, 0]; %From the previous step
prev0.fxnyn  = [0, 0, 0, 0];
prev0.uxnyn0 = [0, 0, 0, 0]; %From previous reversal point
prev0.fxnyn0 = [0, 0, 0, 0];
prev0.uxynOrigin = [0, 0, 0]; %Origin for where to consider displacement. 


%Mindlin normal force
fxyn_MINDLIN = MINDLIN_MONOTONIC([zeros(length(unLoad), 2), unLoad'], pars);
fn_MINDLIN = fxyn_MINDLIN(:, 3)';

% Elliptical contact normal force
[fxyn_ELLIPSOID, uxy_slip] = ELLIPSOID_MONOTONIC([zeros(length(unLoad), 2), unLoad'], pars);
fn_ELLIPSOID = fxyn_ELLIPSOID(:, 3)';


deltam = 0;
Fm = 0;
am = 0;

% Forward Loading Theories
for ii = 1:length(unLoad)

    uxyn_tmp(3) = unLoad(ii);
    
    %Asperity function
    [fxyn, ~, ~, ~, ~, deltam, Fm, am, ~, ~, ~, ahist(ii), ~, Syshist(ii), ~] = ASP_FUN(pars, uxyn_tmp, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
    
    fn_ASP(ii) = fxyn(3);

end

% Elastic Unloading Theory
fn_unload = zeros(3, length(unUnload));
a_unload  = zeros(3, length(unUnload));
unload_legend = cell(3, 1);

am = ahist(end);

% Brake's Unloading:
unload_legend{1} = 'Brake: Analtyical... strain hardening...';

deltabar = deltam*(1 - Fm/(4/3*pars.Estar*sqrt(pars.R)*deltam^(3/2)));
Rebar = Fm^2/( (4/3*pars.Estar)^2 * (deltam - deltabar)^3);
a_unload(1, :) = sqrt(Rebar.*(unUnload - deltabar));
fn_unload(1, :) = 4/3*pars.Estar*sqrt(Rebar)*(unUnload - deltabar).^(3/2).*((unUnload - deltabar) >= 0);
        
% Etsion et al. - From Jackson 2010 Prediction Coef of resitution
unload_legend{2} = 'Etsion et al.: based on am';

deltabar = deltam*(1- 3*Fm/(4*pars.Estar*am*deltam));
Rebar = 4*pars.Estar*(am^3)/(3*Fm);
fn_unload(2, :) = 4/3*pars.Estar*sqrt(Rebar)*(unUnload - deltabar).^(3/2).*((unUnload - deltabar) >= 0);
a_unload(2, :) = sqrt( Rebar * (unUnload - deltabar) );

% Jackson, Chusoipin, and Green
unload_legend{3} = 'Jackson, Chusoipin, and Green: based on FEM fit un res';

deltabar = deltam*1.02*(1 - ( (deltam/pars.delta_y + 5.9)/6.9  )^(-0.54));
Rebar = 1/(deltam - deltabar)^3*( (3*Fm)/(4*pars.Estar) )^2;
fn_unload(3, :) = 4/3*pars.Estar*sqrt(Rebar)*(unUnload - deltabar).^(3/2).*((unUnload - deltabar) >= 0);
a_unload(3, :) = sqrt( Rebar * (unUnload - deltabar) );

figure;
hold on;
xlabel('un');
ylabel('fn');

plot(unLoad, fn_MINDLIN, '-', 'LineWidth', 2);
plot(unLoad, fn_ELLIPSOID, '--', 'LineWidth', 2);
plot(unLoad, fn_ASP, ':', 'LineWidth', 2);
plot(unUnload, real(fn_unload(1, :)), '-', 'LineWidth', 2);
plot(unUnload, real(fn_unload(2, :)), '--', 'LineWidth', 2);
plot(unUnload, real(fn_unload(3, :)), ':', 'LineWidth', 2);

legend('Hertz Sphere', 'Hertz Ellipsoid', 'Asperity', unload_legend{1}, unload_legend{2}, unload_legend{3});


figure;
hold on;
xlabel('un');
ylabel('Contact Radius a');

plot(unLoad, ahist, 'LineWidth', 2);
plot(unUnload, real(a_unload(1, :)), 'LineWidth', 2);
plot(unUnload, real(a_unload(2, :)), '--', 'LineWidth', 2);
plot(unUnload, real(a_unload(3, :)), ':', 'LineWidth', 2);

legend('Loading', unload_legend{1}, unload_legend{2}, unload_legend{3});


%% Test Sys behavior

% unmax = 0.7e-7; %behavior because the strain hardening center is the Sys until it expands?
unmax = 1.5e-5;

un_test = [linspace(0, unmax, 100), linspace(unmax, 0, 100)];
uxyn_tmp = [0, 0, 0];


Syshist = zeros(size(un_test));
fn_ASP  = zeros(size(un_test));

deltam = 0;
Fm = 0;
am = 0;

% Forward Loading Theories
for ii = 1:length(un_test)

    uxyn_tmp(3) = un_test(ii);
    
    %Asperity function
    [fxyn, ~, ~, ~, ~, deltam, Fm, am, ~, ~, ~, ahist(ii), ~, Syshist(ii), ~] = ASP_FUN(pars, uxyn_tmp, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am);
    
    fn_ASP(ii) = fxyn(3);

end


figure;
plot(un_test, fn_ASP, 'LineWidth', 2);
hold on;
ylabel('Normal Force');

set(gcf, 'Renderer', 'painters');


figure;
plot(un_test, Syshist, 'LineWidth', 2);
hold on;
ylabel('Sys');

set(gcf, 'Renderer', 'painters');



