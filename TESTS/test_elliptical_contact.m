% Script tests the determination of elliptical contact parameters. This
% script only considers the normal contact force and pressure distribution
% for now. 

clear;

addpath('../ROUTINES/FRIC_MODELS/ASP_FUN/')

%% Initialize Parameters

Rprime = 30;
Rpprime = 1;

aperb0 = Rprime/Rpprime;

[R, dRdaperb] = CONTACT_ECCENTRICITY_RES(aperb0, Rprime, Rpprime);

%% Test the derivative


delta = aperb0*[0.1, 0.05, 0.01, 0.005, 0.001];
aperb_test = aperb0*[0.5, 1, 1.5];

error = zeros(length(delta), length(aperb_test), 4);

for ii = 1:length(aperb_test)
    
        [~, dRdaperb, ~, ~, ~, dKdaperb, dEdaperb, dedaperb] = CONTACT_ECCENTRICITY_RES(aperb_test(ii), Rprime, Rpprime);

    for jj = 1:length(delta)
        
        [R_high, ~, K_high, E_high, e_high, ~, ~, ~] = CONTACT_ECCENTRICITY_RES(aperb_test(ii)+delta(jj), Rprime, Rpprime);
        [R_low, ~, K_low, E_low, e_low, ~, ~, ~] = CONTACT_ECCENTRICITY_RES(aperb_test(ii)-delta(jj), Rprime, Rpprime);
        
        dRdaperb_Num = (R_high - R_low)/2/delta(jj);
        dKdaperb_Num = (K_high - K_low)/2/delta(jj);
        dEdaperb_Num = (E_high - E_low)/2/delta(jj);
        dedaperb_Num = (e_high - e_low)/2/delta(jj);
        
        error(jj, ii, 1) = (dRdaperb_Num - dRdaperb)/dRdaperb;
        error(jj, ii, 2) = (dKdaperb_Num - dKdaperb)/dKdaperb;
        error(jj, ii, 3) = (dEdaperb_Num - dEdaperb)/dEdaperb;
        error(jj, ii, 4) = (dedaperb_Num - dedaperb)/dedaperb;
        
    end
    
end

format long
% Rows should have decreasing error compared to previous row
error
format


%% Numerically solve the problem

Rprime = 130;
Rpprime = 1;

aperb0 = Rprime/Rpprime;

fopts = optimoptions('fsolve','Display','iter',...
                            'FunctionTolerance', 1e-20, ...
                            'MaxIterations', 400, ...
                            'MaxFunctionEvaluations', 400, ...
                            'StepTolerance', 1e-20);

aperb_sol = fsolve(@(x)CONTACT_ECCENTRICITY_RES(x, Rprime, Rpprime), aperb0, fopts);


[R, dRdaperb] = CONTACT_ECCENTRICITY_RES(aperb_sol, Rprime, Rpprime);


%% Plot the ellipse of contact v. the ellipse of constant heights.

%Rprime associated with x

theta = linspace(0, 2*pi, 10000);


%%%% Geometric ellipse
% 3D shape is approximated as z = 1/(2*Rprime1)*x^2 + 1/(2*Rpprime1)*y^2

Rprime1 = Rprime*2;
Rpprime1 = Rpprime*2;

zunity_y = 1/(2*Rprime1)*aperb_sol^2;
zunity_y = 1/(2*Rpprime1)*1^2;

major_geo = sqrt(2*Rprime1*zunity_y);
minor_geo = sqrt(2*Rpprime1*zunity_y);

r_geometry = sqrt(major_geo.^2.*minor_geo.^2 ./ (minor_geo.^2.*cos(theta).^2 + major_geo.^2.*sin(theta).^2));

x_geo = r_geometry.*cos(theta);
y_geo = r_geometry.*sin(theta);

%%%% Contact Ellipse
a = aperb_sol;
b = 1;
r_cont = sqrt(a.^2.*b.^2 ./ (b.^2.*cos(theta).^2 + a.^2.*sin(theta).^2));

x_cont = r_cont.*cos(theta);
y_cont = r_cont.*sin(theta);


f = figure;
hold on;
xlabel('x');
ylabel('y');

plot(x_geo, y_geo, 'LineWidth', 2)
plot(x_cont, y_cont, 'LineWidth', 2)


legend('Geometric Ellipse', 'Contact Ellipse');


xrange = xlim;
yrange = ylim;

pbaspect([xrange(2)-xrange(1), yrange(2)-yrange(1), 1]);

f.Position = [461 678 1460 420];

fprintf('From Johnson pg96, the contact ellipse should be more slender than the geometric ellipse\n');

%% Normal Contact Tests - Initialize Parameters

%Initialize the sphere parameters
E = 192.85e9; %304L from Traction paper
nu = 0.29; %304L taken from Traction paper
G = E / 2/ (1 + nu);%304L approx Shear modulus

%Mean of top and bottom
Rx = mean([0.000273924650289, 0.000314503650665]); 
Ry = mean([0.008183455622513, 0.008521098239706]);

pars.mu     = 0.3; %may want to do a function of displacement - see Eriten's work
% pars.R      = (1/Rx + 1/Ry)^(-1);%, 1e-4 is about a maximum %R is not used when doing ellipsoid modeling

pars.Rprime = (1/0.008183455622513 + 1/0.008521098239706)^(-1);
pars.Rpprime = (1/0.000273924650289 + 1/0.000314503650665)^(-1);
% pars.Rpprime = 1/1*pars.Rprime ;

pars.E      = E;
pars.Estar  = E/2/(1 - nu^2);
pars.nu     = nu;
pars.G      = G;
pars.Gstar  = G/2/(2 - nu);

pars = INITIALIZE_ECCENTRIC(pars);

%% Initialize General Info

Nqp_radius = 100;

tx0 = zeros(1, Nqp_radius); 
ty0 = zeros(1, Nqp_radius); 


rq = linspace(0, 1, Nqp_radius); 
wq = ones(size(rq));
wq(2:end-1) = 2;
wq = wq/sum(wq);

rq0 = linspace(0, 1, Nqp_radius); 


uxyn0 = [0, 0, 0]; %displacement at previous step. 

%% Normal Contact Function

un = linspace(0, 1e-5, 100);

fn = zeros(size(un));

for ii = 1:length(un)
    
    uxyn = uxyn0;
    uxyn(3)= un(ii);
    
    [fxyn, dfxynduxyn, rq0, tx0, ty0, b, dbdun, a, dadun, p0, dp0dun] = ELLIPSOID_PRE(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq);
    
    fn(ii) = fxyn(3);
end

figure;
hold on;
xlabel('un');
ylabel('fn');

plot(un, fn, '-', 'LineWidth', 2.5);

%% Test Derivatives of the Normal contact model

delta = [1e-8, 1e-9, 1e-10, 1e-11]; %displacement variation (m)

un_test = [1e-6, 2.5e-5, 1e-5];

error_Nbap0 = zeros(length(delta), length(un_test), 4);

uxyn0 = [0, 0, 0]; %displacement at previous step. 



for jj = 1:length(un_test)

    uxyn = uxyn0;
    uxyn(3) = un_test(jj);
    
    [fxyn, dfxynduxyn, rq0, tx0, ty0, b, dbdun, a, dadun, p0, dp0dun] = ELLIPSOID_PRE(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq);

    
    for ii = 1:length(delta)


        uxyn_low = uxyn;
        uxyn_high = uxyn;

        uxyn_low(3) = uxyn_low(3) - delta(ii);
        uxyn_high(3) = uxyn_high(3) + delta(ii);
        
        [fxyn_high, ~, rq0, tx0, ty0, b_high, ~, a_high, ~, p0_high, ~] = ELLIPSOID_PRE(pars, uxyn_high, uxyn0, rq0, tx0, ty0, rq, wq);
        [fxyn_low, ~, rq0, tx0, ty0, b_low, ~, a_low, ~, p0_low, ~] = ELLIPSOID_PRE(pars, uxyn_low, uxyn0, rq0, tx0, ty0, rq, wq);

        dfxynduxyn_num = (fxyn_high(3) - fxyn_low(3))/2/delta(ii);
        dbdun_num = (b_high - b_low)/2/delta(ii);
        dadun_num = (a_high - a_low)/2/delta(ii);
        dp0dun_num = (p0_high - p0_low)/2/delta(ii);


        error_Nbap0(ii, jj, 1) = norm(dfxynduxyn(3,3) - dfxynduxyn_num)/norm(dfxynduxyn);
        error_Nbap0(ii, jj, 2) = norm(dbdun - dbdun_num)/norm(dbdun);
        error_Nbap0(ii, jj, 3) = norm(dadun - dadun_num)/norm(dadun);
        error_Nbap0(ii, jj, 4) = norm(dp0dun - dp0dun_num)/norm(dp0dun);

    end
end


format long
% Rows should have decreasing error compared to previous row
error
format

%% Test Integration Method over an ellipse.

un_test = 1e-5;
uxyn = [0, 0, un_test];

[fxyn, dfxynduxyn, rq0, tx0, ty0, b, dbdun, a, dadun, p0, dp0dun] = ELLIPSOID_PRE(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq);


%Long dimension is the yaxis, corresponds to a
Nqp_radius = 1000;

tx0 = zeros(1, Nqp_radius); 
ty0 = zeros(1, Nqp_radius); 

rq = linspace(0, 1, Nqp_radius); 
wq = ones(size(rq));
wq(2:end-1) = 2;
wq = wq/sum(wq);

%assume that rq in [0, 1]

% pressure distribution w.r.t. rho. - rq is rho at quadrature points when rq
% in [0, 1]
p_rho = p0*sqrt(1-rq.^2);

%Normal force approximation:
N = 2*pi*a*b*(p_rho.*rq)*wq';

% [N, fxyn(3)]
error_N = abs(N-fxyn(3))/fxyn(3)

%% Verify the Implementation of Elliptical Contact Functions

clear;

Rpprime = 1;
% Rprime = linspace(2, 900, 500);
Rprime = logspace(log10(1), log10(900), 500);

f1exact = zeros(size(Rprime));
f2exact = zeros(size(Rprime));

eccentricity = zeros(size(Rprime));
aperb = zeros(size(Rprime));

f1approx = zeros(size(Rprime));
f2approx = zeros(size(Rprime));


for ii = 1:length(Rprime)

    %%% F1 and F2 from Johnson's book (exact)

    pars.Rprime = Rprime(ii);
    pars.Rpprime = Rpprime;
    pars = INITIALIZE_ECCENTRIC(pars);
    pars.R = pars.Re;

    f1exact(ii) = pars.F1;
    f2exact(ii) = pars.F2;
    
    eccentricity(ii) = pars.e;
    aperb(ii) = pars.aperb;

    %%% F1 and F2 from Becker Paper (Approximate)

    Rprime1 = Rprime(ii)*2;
    Rpprime1 = Rpprime*2;

    Rprime2 = Rprime(ii)*2;
    Rpprime2 = Rpprime*2;

%     mat = [1, 1; -1, 1];
%     rhs = [0.5*(1/Rprime1 + 1/Rprime2 + 1/Rpprime1 + 1/Rpprime2);
%            0.5*( (1/Rprime1 - 1/Rpprime1)^2 + (1/Rprime2 - 1/Rpprime2)^2 ...
%                  + 2*(1/Rprime1 - 1/Rpprime1)*(1/Rprime2 - 1/Rpprime2)*cos(2*0))];
% 
%     AB = mat \ rhs;

    f1approx(ii) = 1 - ( ( Rprime(ii) / Rpprime )^0.0602 - 1 )^1.456;
    f2approx(ii) = 1 - ( ( Rprime(ii) / Rpprime )^0.0684 - 1 )^1.531;
    
% Rprime(ii) / Rpprime
% Rpprime / Rprime(ii)

end

linewidth = 4;

figure;
hold on;
xlabel('sqrt(Rprime / Rpprime)');
ylabel('f1, f2')
title('Compare to Figure 4.4 (pg97) in Johnson');

plot(sqrt(Rprime./Rpprime), f1exact, 'k-', 'LineWidth', linewidth);
plot(sqrt(Rprime./Rpprime), f1approx, 'r--', 'LineWidth', linewidth);


plot(sqrt(Rprime./Rpprime), f2exact, 'g:', 'LineWidth', linewidth);
plot(sqrt(Rprime./Rpprime), f2approx, 'c-.', 'LineWidth', linewidth);


plot(sqrt(Rprime./Rpprime), eccentricity, 'b-', 'LineWidth', linewidth);
plot(sqrt(Rprime./Rpprime), sqrt(Rprime./Rpprime)./aperb, 'm-', 'LineWidth', linewidth);


legend('F1 Exact', 'F1 Approx', 'F2 Exact', 'F2 Approx', 'Eccentricity', 'b/a*sqrt(R''/R")');

ylim([0.4, 1]);
grid on;

set(gca, 'xscale', 'log');

set(gcf, 'Renderer', 'painters');

