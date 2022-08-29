% This script seeks to recreate the graphs from:
%
% Strain Hardening From Elastic-Perfectly Plastic to Perfectly Elastic
% Flattening Single Asperity Contact
% by H. Ghaednia, M.R.W. Brake, M. Berryhill, R.L. Jackson, 2019, Journal of
% Tribology
%
% to verify that the normal plasticity model is correctly implemented.
%
% This uses 1 function to implement the equations as given in the paper.
% Then at the end, my sphere on sphere equations/corrections are verified
% against the implementation that is compared to the paper. 

clear;

addpath('../ROUTINES/FRIC_MODELS')
addpath('../ROUTINES/FRIC_MODELS/ASP_FUN/')
addpath('../ROUTINES/FRIC_MODELS/MINDLIN/')
addpath('../ROUTINES/FRIC_MODELS/PLASTIC/')

ASP_FUN = @GHAEDNIA_ET_AL_FLATTENING;
ASP_FUN_SPHERES = @SPHERE_PLASTIC_PRE;

%% Initialize Parameters

R = 1e-3; % See start of Section 2.2. 


% Pars 1
Es = 200e9;
Ef = 200e9;
Et = 0.04*Es;
Sys = 500e6;
nus = 0.33;
nuf = 0.33;

pars1 = struct('Es', Es, 'Ef', Ef, 'Et', Et, 'nus', nus, 'nuf', nuf, 'Sys', Sys, 'R', R);

% Pars 2
Es = 200e9;
Ef = 200e9;
Et = 0.02*Es;
Sys = 750e6;
nus = 0.3;
nuf = 0.3;

pars2 = struct('Es', Es, 'Ef', Ef, 'Et', Et, 'nus', nus, 'nuf', nuf, 'Sys', Sys, 'R', R);

% Pars 3
Es = 71e9;
Ef = 200e9;
Et = 0.005*Es;
Sys = 200e6;
nus = 0.33;
nuf = 0.33;

pars3 = struct('Es', Es, 'Ef', Ef, 'Et', Et, 'nus', nus, 'nuf', nuf, 'Sys', Sys, 'R', R);

% Pars 4
Es = 71e9;
Ef = 200e9;
Et = 0.1*Es;
Sys = 500e6;
nus = 0.33;
nuf = 0.3;

pars4 = struct('Es', Es, 'Ef', Ef, 'Et', Et, 'nus', nus, 'nuf', nuf, 'Sys', Sys, 'R', R);

% Pars 5
Es = 200e9;
Ef = 200e9;
Et = 0.02*Es;
Sys = 300e6;
nus = 0.33;
nuf = 0.3;


pars5 = struct('Es', Es, 'Ef', Ef, 'Et', Et, 'nus', nus, 'nuf', nuf, 'Sys', Sys, 'R', R);


%% Figure 8, page 7
pars = pars1;
fignum = 8;

Delta_perR = linspace(0, 0.05, 100);

xlims = [0. 0.05];
ylimsdeltas = [0.4, 1];
ylimsa = [0, 4e-4];
ylimsF = [0, 2000];
give_ylims = true;

create_plots(pars, fignum, Delta_perR, xlims, ylimsdeltas, ylimsa, ylimsF, give_ylims);


%% Figure 9, page 8

pars = pars2;
fignum = 9;

Delta_perR = linspace(0, 0.05, 100);


xlims = [0. 0.05];
ylimsdeltas = [0.4, 1];
ylimsa = [0, 4e-4];
ylimsF = [0, 2000];
give_ylims = true;

create_plots(pars, fignum, Delta_perR, xlims, ylimsdeltas, ylimsa, ylimsF, give_ylims);

%% Figure 10, page 8

pars = pars3;
fignum = 10;

Delta_perR = linspace(0, 0.05, 100);


xlims = [0. 0.05];
ylimsdeltas = [0.5, 1];
ylimsa = [0, 4e-4];
ylimsF = [0, 1000];
give_ylims = true;

create_plots(pars, fignum, Delta_perR, xlims, ylimsdeltas, ylimsa, ylimsF, give_ylims);



%% Figure 11, page 8

pars = pars4;
fignum = 11;

Delta_perR = linspace(0, 0.05, 100);


xlims = [0. 0.05];
ylimsdeltas = [0.5, 1];
ylimsa = [0, 4e-4];
ylimsF = [0, 1000];
give_ylims = true;

create_plots(pars, fignum, Delta_perR, xlims, ylimsdeltas, ylimsa, ylimsF, give_ylims);



%% Figure 12, page 9

pars = pars5;
fignum = 12;

Delta_perR = linspace(0, 0.05, 100);

xlims = [0. 0.05];
ylimsdeltas = [0.4, 1];
ylimsa = [0, 4e-4];
ylimsF = [0, 2000];
give_ylims = true;

create_plots(pars, fignum, Delta_perR, xlims, ylimsdeltas, ylimsa, ylimsF, give_ylims);



%% Personal Test of the model - Elastic Transition

pars = pars5;
fignum = -100;

Delta_perR = linspace(0, 10e-5, 1000);

xlims = [0. Delta_perR(end)];
ylimsa = []; %[0, 4e-4];
ylimsF = []; %[0, 2000];
give_ylims = false;

create_plots(pars, fignum, Delta_perR, xlims, ylimsdeltas, ylimsa, ylimsF, give_ylims);

%% Cross Comparison with my Sphere on Sphere Implementation

fignum = -200;

Es = 192.85e9; %304L from Traction paper
Ef = Inf;
Et = 0.02*Es;
Sys = 200e6; %Pa
nus = 0.29; %304L taken from Traction paper
nuf = 0.29; %Shouldn't matter
R = 7e-4;

pars6 = struct('Es', Es, 'Ef', Ef, 'Et', Et, 'nus', nus, 'nuf', nuf, 'Sys', Sys, 'R', R);

pars6s = initialize_pars_spheres(R, Es, Et, nus, Sys);

Delta = [linspace(0, 2e-8, 100), linspace(2e-8, 2e-5, 100)];
Delta = [linspace(0, 2e-8, 100), linspace(2e-8, pars6s.delta_y*4e3, 100)];
% Delta = linspace(0, 2e-8, 100);
% Delta = 1e-8;

give_ylims = false;

[fnhist_flat, ahist_flat] = create_plots(pars6, fignum, Delta/R, [0 (Delta(end)/R)], [], [], [], give_ylims);


figtitle = 'Comparison of my implementations';

pars6s.sliptype = 1;

create_plots_spheres(ASP_FUN_SPHERES, pars6s, figtitle, Delta/R, fnhist_flat, ahist_flat);



%% Function for initializing parameters
function [pars] = initialize_pars_spheres(R, E, Et, nu, Sys)

    pars.R = R/2;
    pars.Re = R/2;
    
    G = E / 2/ (1 + nu);%304L approx Shear modulus
    
    % pars.R      = (1/Rx + 1/Ry)^(-1);%, 1e-4 is about a maximum
    pars.E      = E;
    pars.Estar  = E/2/(1 - nu^2);
    pars.nu     = nu;
    pars.G      = G;
    pars.Gstar  = G/2/(2 - nu);

    % %% Yield Parameters
    C = 1.295*exp(0.736*pars.nu);

    delta_y1s = (pi*C*Sys/(2*(2*pars.Estar)))^2*(2*pars.R); %displacement of one sphere against a rigid flat to cause yielding.

    % Store all of these results as system parameters to prevent recalculation:
    pars.delta_y = delta_y1s*2;
    pars.Sys = Sys;

    %tangent stiffness after yield.
    pars.Et = Et;
    
end

%% Function for figure creation
function [fnhist, ahist] = create_plots(pars, fignum, Delta_perR, xlims, ylimsdeltas, ylimsa, ylimsF,give_ylims)


    Delta = Delta_perR.*pars.R;

    ahist = zeros(size(Delta));
    fnhist = zeros(size(Delta));
    
    aehist = zeros(size(Delta));
    fehist = zeros(size(Delta));

    aphist = zeros(size(Delta));
    fphist = zeros(size(Delta));
    
    deltas_star = zeros(size(Delta));
    
    Upsilon = zeros(size(Delta));

    for ii = 1:length(Delta)

        %Asperity function
        [fnhist(ii), ahist(ii), fehist(ii), aehist(ii), fphist(ii), aphist(ii), deltas_star(ii), Upsilon(ii)] = GHAEDNIA_ET_AL_FLATTENING(pars.Es, pars.nus, pars.Ef, pars.nuf, pars.Et, pars.Sys, pars.R, Delta(ii));


    end

    linewidth = 2.5;
    
    figure;
    subplot(3,1,1);
    hold on;
    title(sprintf('Figure %u', fignum));
    grid on;
    xlabel('\Delta / R');
    ylabel('\delta_s^*');

    plot(Delta_perR, deltas_star, 'LineWidth', linewidth);
    xlim(xlims);
    if(give_ylims)
        ylim(ylimsdeltas);
    end
        
    
    subplot(3,1,2);
%     title(sprintf('Figure %u', fignum));
    hold on;
    grid on;
    xlabel('\Delta / R');
    ylabel('a (m)');

    plot(Delta_perR, aehist, '--', 'LineWidth', linewidth);
    plot(Delta_perR, aphist, '-', 'LineWidth', 1.5*linewidth);
    plot(Delta_perR, ahist, 'LineWidth', linewidth);
    
    xlim(xlims);
    if(give_ylims)
        ylim(ylimsa);
    end
    
    legend('Hertzian', 'JG Model', 'Ghaednia et al', 'Location', 'nw');

    subplot(3,1,3);
    hold on;
    grid on;
    xlabel('\Delta / R');
    ylabel('Force (N)');

    plot(Delta_perR, fehist, '--', 'LineWidth', linewidth);
    plot(Delta_perR, fphist, '-', 'LineWidth', 1.5*linewidth);
    plot(Delta_perR, fnhist, 'LineWidth', linewidth);
    xlim(xlims);
    
    if(give_ylims)
        ylim(ylimsF);
    end
    
    legend('Hertzian', 'JG Model', 'Ghaednia et al', 'Location', 'nw');
    set(gcf, 'Renderer', 'painters');
    
    
    fig = gcf;
    first_pos = fig.Position;
    
    newWidth = 830;
    newHeight = 1000;
    
    fig.Position = [first_pos(1) + (first_pos(3) - newWidth), first_pos(2) + (first_pos(4) - newHeight), newWidth, newHeight];
    
    figure; 
    title(sprintf('Figure %u', fignum));
    hold on;
    xlabel('Delta/R');
    ylabel('Upsilon');
    plot(Delta_perR, Upsilon, 'LineWidth', linewidth);

end

%% Function for testing sphere on sphere model
function [] = create_plots_spheres(ASP_FUN, pars, figtitle, Delta_perR, fnhist_flat, ahist_flat)


    un = 2*Delta_perR.*(2*pars.R);

    % My calculations:
    deltam0 = 0;
    Fm0 = 0;
    am0 = 0;
    uxyn0 = [0, 0, 0];
    rq0 = 0;
    tx0 = 0;
    ty0 = 0;
    rq = 0;
    wq = 1;

    ahist = zeros(size(un));
    fnhist = zeros(size(un));


    for ii = 1:length(un)

        uxyn_tmp(3) = un(ii);

        %Asperity function
        [fxyn, ~, ~, ~, ~, deltam, Fm, am, ~, ~, ~, ahist(ii), ~, ~, ~] = ASP_FUN(pars, uxyn_tmp, uxyn0, rq0, tx0, ty0, rq, wq, deltam0, Fm0, am0);

        fnhist(ii) = fxyn(3);

    end
    
    
    linewidth = 2.5;

    figure;
    subplot(2,1,1);
    title(figtitle);
    hold on;
    grid on;
    xlabel('\Delta / R');
    ylabel('a (m)');

    plot(Delta_perR, ahist, 'LineWidth', linewidth);
    plot(Delta_perR, ahist_flat, '--', 'LineWidth', linewidth);
    
    legend('Sphere Flat Implementation', 'Sphere-Sphere Implementation', 'Location', 'se');
    
    set(gcf, 'Renderer', 'painters');

    subplot(2,1,2);
    hold on;
    grid on;
    xlabel('\Delta / R');
    ylabel('Force (N)');

    plot(Delta_perR, fnhist, 'LineWidth', linewidth);
    plot(Delta_perR, fnhist_flat, '--', 'LineWidth', linewidth);
    
    legend('Sphere Flat Implementation', 'Sphere-Sphere Implementation', 'Location', 'nw');
    
    set(gcf, 'Renderer', 'painters');
    
    
    fprintf('The maximum difference between model forces is %s \n', max(abs(fnhist-fnhist_flat)));
    fprintf('The maximum difference between model contact radii is: %s \n', max(abs(fnhist-fnhist_flat)));

end
