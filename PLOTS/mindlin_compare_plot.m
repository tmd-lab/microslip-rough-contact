% Script for plotting comparisons between asperity contact models and
% Mindlin Theory. 

clear;

addpath('../ROUTINES/FRIC_MODELS/ASP_FUN/')
addpath('../ROUTINES/FRIC_MODELS/PLASTIC/')
addpath('../ROUTINES/FRIC_MODELS/MINDLIN/')



% set(0, 'defaultaxesinterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

addpath('../ROUTINES/')

%% Define Parameters 
% Parameters are shared between all models 

%Initialize the sphere parameters
E = 192.85e9; %304L from Traction paper
nu = 0.29; %304L taken from Traction paper
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

pars.Rpprime = pars.Rprime; % Set to spherical contact for plots

pars = INITIALIZE_ECCENTRIC(pars);

pars.R = pars.Re;

pars.mu = 0.2;

% %% Plasticity Parameters 
pars.sliptype = 1;

Sys = Inf;
C = 1.295*exp(0.736*pars.nu);
delta_y1s = (pi*C*Sys/(2*(2*pars.Estar)))^2*(2*pars.R); %displacement of one sphere against a rigid flat to cause yielding.

pars.delta_y = delta_y1s*2;
pars.Sys = Sys;
pars.Et = 0.01*pars.E;

%% Asperity Model Lists + Details

%Functions to make comparison graphs
ASP_FUN_LIST = {@ELLIPSOID_TAN_DECOUPLE,...
                @ELLIPSOID_IWAN_DECOUPLE, ...
                @ELLIPSOID_IWAN_FIT_DECOUPLE, ...
                @ELLIPSOID_IWAN_FIT_COUPLE, ...
                @SPHERE_PL_TAN_DECOUPLE_PL_NORM};
            
            
ASP_NAMES_LIST = {'Tangent', 'Mindlin-Iwan Constant (MIC)', 'Mindlin-Iwan Fit (MIF)', 'Iwan Fit Coupled', 'Plastic (elastic regime) Sphere Tangent'};
                
% ASP_NQP_RAD_LIST = [1, 10000, 10000, 10000];
ASP_NQP_RAD_LIST = [1, 100, 100, 100, 1];

ASP_LIST_XYINDEX = [1, 1, 1, 1, 1]; % Test x
% ASP_LIST_XYINDEX = [2, 2, 2, 2, 2]; % Test y

EXTRA_ARGS = [false, false, false, false, true];

UY = 0.1e-5; % If set to a constant this will change the behavior of the coupled model
NUM_PLOT = 3; % Set to 3 for paper plots
% NUM_PLOT = length(ASP_FUN_LIST); % Set to 3 for paper plots


color_plot = DISTINGUISHABLE_COLORS(6, 'w');

color_plot(4:5, :) = color_plot(5:6, :);

%% Compare all of the Initial Stiffness Jacobians

un = 1e-5;
ref_ind = 1; % which asperity function to compare against

% Initialize storage

tx0 = cell(size(ASP_FUN_LIST));
ty0 = cell(size(ASP_FUN_LIST));
rq = cell(size(ASP_FUN_LIST));
rq0 = cell(size(ASP_FUN_LIST));
wq = cell(size(ASP_FUN_LIST));
uxyn0_ASP = cell(size(ASP_FUN_LIST));

for ii = 1:length(ASP_FUN_LIST)
    
    tx0{ii} = zeros(1, ASP_NQP_RAD_LIST(ii));
    ty0{ii} = zeros(1, ASP_NQP_RAD_LIST(ii));
    rq{ii} = linspace(0, 1, ASP_NQP_RAD_LIST(ii)); 
    rq0{ii} = linspace(0, 1, ASP_NQP_RAD_LIST(ii)); 
    
    
    wq{ii} = ones(size(rq{ii}));
    wq{ii}(2:end-1) = 2;
    wq{ii} = wq{ii}/sum(wq{ii});
    
    uxyn0_ASP{ii} = [0,0,0];
end


%%%%%%%%%%%% Get Jacobian stiffnesses at initial normal displacement

dfxynduxyn = cell(size(ASP_FUN_LIST));

%Initialize to the first normal load for the ASP_FUN as well 
for kk = 1:length(ASP_FUN_LIST)

    uxyn_ASP = [0, 0, un];

    if(EXTRA_ARGS(kk) )
        [~, dfxynduxyn{kk}, rq0{kk}, tx0{kk}, ty0{kk}] ...
            = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk}, 0, 0, 0);
    else
        [~, dfxynduxyn{kk}, rq0{kk}, tx0{kk}, ty0{kk}] ...
            = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk});
    end

    uxyn0_ASP{kk} = uxyn_ASP;
    

end

for kk = 1:length(ASP_FUN_LIST)

    fprintf('Jacobian for function %u differed from %u by %s \n', kk, ...
        ref_ind, norm(dfxynduxyn{kk} - dfxynduxyn{ref_ind}) / norm(dfxynduxyn{ref_ind}));
end

fprintf('Note that for Nq = 100, errors on the order of 1e-2 are expected, errors decrease for more integration radii \n');


%% Comparison Plot 1 - Constant Normal Load

if( sum(ASP_NQP_RAD_LIST > 100) )
    warning('Excessive number of quadrature radii used from here on.');
end

hyst_amp = 1.5; 
un = 1e-5;
N_test = 1000;

[fxyn, u_slip] = MINDLIN_MONOTONIC([0, 0, un], pars);
N1 = fxyn(1,3);  

ux_test = [linspace(0, hyst_amp*u_slip, N_test), ...
          linspace(hyst_amp*u_slip, -hyst_amp*u_slip, N_test), ...
          linspace(-hyst_amp*u_slip, hyst_amp*u_slip, N_test)];
      
un_test = un*ones(size(ux_test));

%%%% Generate Mindlin Solution data w/Masing assumptions
[fxyn, ~] = MINDLIN_MONOTONIC([ux_test(1:N_test)', zeros(N_test, 1), un_test(1:N_test)'], pars);

%monotonic loading
fx_test(1:N_test) = fxyn(:, 1)';

%hysteretic unloading
fx_test((N_test+1):(2*N_test)) = fxyn(N_test, 1)' - 2*fxyn(:, 1)';

%hysteretic unloading
fx_test((2*N_test+1):(3*N_test)) = -fxyn(N_test, 1)' + 2*fxyn(:, 1)';

[fx_test_ASP] = generate_data(ux_test, un_test, ASP_FUN_LIST, ...
                    ASP_NQP_RAD_LIST, ASP_LIST_XYINDEX, EXTRA_ARGS, UY, pars);

%% Plot the comparison 1 

ulimits = hyst_amp*[-1.1, 1.1];
flimits = [-1.1, 1.1];

plot_comparison('ConstN', pars, N1, ASP_NAMES_LIST, ux_test, fx_test, ...
    fx_test_ASP, u_slip, color_plot, NUM_PLOT, ulimits, flimits)

plot_comparison('', pars, N1, ASP_NAMES_LIST, ux_test, fx_test, ...
    fx_test_ASP, u_slip, color_plot, length(ASP_NAMES_LIST), ulimits, flimits)

%% Comparison Plot 2 - Decrease Normal Load In Middle

hyst_amp = 1.2; 
un1 = 1e-5;
un2 = 0.6^(2/3) * un1; % 60% Force
N_test = 1000;



[fxyn, u_slip1] = MINDLIN_MONOTONIC([0, 0, un1], pars);
N1 = fxyn(1,3);  

[fxyn, u_slip2] = MINDLIN_MONOTONIC([0, 0, un2], pars);
N2 = fxyn(1,3);  


% Histories to consider
un_test = [un1*ones(1, N_test), un2*ones(1, N_test)];

ux_trans = 0.5*u_slip2;

ux_test = [linspace(0, ux_trans, N_test), ...
            linspace(ux_trans, hyst_amp*u_slip2, N_test)];
        
% Mindlin Solution with displacement control
fx_test = zeros(size(ux_test));

% Load at high normal load
[fxyn, ~] = MINDLIN_MONOTONIC([ux_test(1:N_test)', zeros(N_test, 1), un_test(1:N_test)'], pars);
fx_test(1:N_test) = fxyn(:, 1)';


% Integration of increasing T, Decreasing N transition
dNdT = -1/pars.mu * 10;
Trange = [fx_test(N_test), fx_test(N_test)+(N2 - N1)/dNdT];

dNUx_dT = @(T, NUx)MINDLIN_INCREASET_DECREASEN_DOT(T, NUx, pars, dNdT);
NUx0 = [N1; ux_trans];
options = odeset('MaxStep',abs(range(Trange))/100);
[T, NUx] = ode45(dNUx_dT, Trange, NUx0, options);

% Cut the Transition Curves to be appropriately sized
%%%%%%%FILL
%%%%%%%%%%%%%%%%%%%%%%%%%%

UN_insert = (9.*NUx(:, 1)'.^2/16/pars.R/pars.Estar^2).^(1/3);

% Second portion - constant at N2
ux_test(N_test+1:2*N_test) = linspace(NUx(end, 2), hyst_amp*u_slip2, N_test);

[fxyn, ~] = MINDLIN_MONOTONIC([ux_test(N_test+1:2*N_test)', zeros(N_test, 1), un_test(N_test+1:2*N_test)'], pars);
fx_test(N_test+1:2*N_test) = fxyn(:, 1)';


%%%%%%%FILL
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Insert the integrated curve
un_test = [un_test(1:N_test), UN_insert, un_test(N_test+1:end)];
ux_test = [ux_test(1:N_test), NUx(:, 2)', ux_test(N_test+1:end)];
fx_test = [fx_test(1:N_test), T', fx_test(N_test+1:end)];
fn_test = [N1*ones(1, N_test), NUx(:, 1)', N2*ones(1, N_test)];


[fx_test_ASP] = generate_data(ux_test, un_test, ASP_FUN_LIST, ...
                    ASP_NQP_RAD_LIST, ASP_LIST_XYINDEX, EXTRA_ARGS, UY, pars);

%% Plot the comparison 2 

ulimits = hyst_amp*[0, 1]*u_slip2/u_slip1;
flimits = [0, 1.2]*N2/N1;
Nlimits = [0.55, 1.05];

% fn_test = N1*(un_test == un1) + N2 * (un_test == un2);

plot_comparison('DecN', pars, N1, ASP_NAMES_LIST, ux_test, fx_test, ...
    fx_test_ASP, u_slip1, color_plot, NUM_PLOT, ulimits, flimits, ...
    N2, u_slip2, fn_test, Nlimits)

plot_comparison('', pars, N1, ASP_NAMES_LIST, ux_test, fx_test, ...
    fx_test_ASP, u_slip1, color_plot, length(ASP_NAMES_LIST), ulimits, flimits, ...
    N2, u_slip2, fn_test, Nlimits)

%% Comparison 3 - Increasing N and T (Mindlin and Deresiewicz 1953, Section 7)

hyst_amp = 1.2; 
un1 = 1e-5;
un2 = 1.2^(2/3) * un1; % 120% Force
N_test = 1000;

[fxyn, u_slip1] = MINDLIN_MONOTONIC([0, 0, un1], pars);
N1 = fxyn(1,3);  

[fxyn, u_slip2] = MINDLIN_MONOTONIC([0, 0, un2], pars);
N2 = fxyn(1,3);  


% Histories to consider
un_test = [un1*ones(1, N_test), un2*ones(1, N_test)];

ux_trans = 0.6*u_slip1;

ux_test = [linspace(0, ux_trans, N_test), ...
            linspace(ux_trans, hyst_amp*u_slip2, N_test)];
        
% Mindlin Solution with displacement control
fx_test = zeros(size(ux_test));

% Load at high normal load
[fxyn, ~] = MINDLIN_MONOTONIC([ux_test(1:N_test)', zeros(N_test, 1), un_test(1:N_test)'], pars);
fx_test(1:N_test) = fxyn(:, 1)';

% Drop to low normal load curve
[fxyn, ~] = MINDLIN_MONOTONIC([ux_test(N_test+1:2*N_test)', zeros(N_test, 1), un_test(N_test+1:2*N_test)'], pars);
fx_test(N_test+1:2*N_test) = fxyn(:, 1)';

% Apply correction for the transition region
trans_mask = (fx_test > fx_test(N_test)) & (fx_test < fx_test(N_test) + pars.mu*(N2 - N1));

fx_test(trans_mask) = NaN;
trans_offset = sum(trans_mask)+1;

% Integration of increasing transition
Trange = [fx_test(N_test), fx_test(N_test)+pars.mu*(N2 - N1)];
dNUx_dT = @(T, NUx)MINDLIN_INCREASE_DOT(T, NUx, pars);
NUx0 = [N1; ux_trans];
options = odeset('MaxStep',pars.mu*(N2 - N1)/100);
[T, NUx] = ode45(dNUx_dT, Trange, NUx0, options);
UN_insert = (9.*NUx(:, 1)'.^2/16/pars.R/pars.Estar^2).^(1/3);

% Insert the integrated curve
un_test = [un_test(1:N_test), UN_insert, un_test(N_test+trans_offset:end)];
ux_test = [ux_test(1:N_test), NUx(:, 2)', ux_test(N_test+trans_offset:end)];
fx_test = [fx_test(1:N_test), T', fx_test(N_test+trans_offset:end)];
fn_test = [N1*ones(1, N_test), NUx(:, 1)', N2*ones(1, N_test-trans_offset+1)];

% Generate the asperity data
[fx_test_ASP] = generate_data(ux_test, un_test, ASP_FUN_LIST, ...
                    ASP_NQP_RAD_LIST, ASP_LIST_XYINDEX, EXTRA_ARGS, UY, pars);

%% Plot the comparison 3 

ulimits = hyst_amp*[0, 1]*u_slip2/u_slip1;
flimits = [0, 1.05]*N2/N1;
Nlimits = [0.95, 1.25];

% fn_test = N1*(un_test == un1) + N2 * (un_test == un2);

plot_comparison('IncN', pars, N1, ASP_NAMES_LIST, ux_test, fx_test, ...
    fx_test_ASP, u_slip1, color_plot, NUM_PLOT, ulimits, flimits, ...
    N2, u_slip2, fn_test, Nlimits)

plot_comparison('', pars, N1, ASP_NAMES_LIST, ux_test, fx_test, ...
    fx_test_ASP, u_slip1, color_plot, length(ASP_NAMES_LIST), ulimits, flimits, ...
    N2, u_slip2, fn_test, Nlimits)



%% Comparison Plot 1 - Constant Normal Load (again)

if( sum(ASP_NQP_RAD_LIST > 100) )
    warning('Excessive number of quadrature radii used from here on.');
end

hyst_amp = 1.5; 
un = 1e-5;
N_test = 5001;

[fxyn, u_slip] = MINDLIN_MONOTONIC([0, 0, un], pars);
N1 = fxyn(1,3);  

ux_test = linspace(0, hyst_amp*u_slip, N_test);
      
un_test = un*ones(size(ux_test));

%%%% Generate Mindlin Solution data w/Masing assumptions
[fxyn, ~] = MINDLIN_MONOTONIC([ux_test(1:N_test)', zeros(N_test, 1), un_test(1:N_test)'], pars);

%monotonic loading
fx_test = fxyn(:, 1)';


[fx_test_ASP] = generate_data(ux_test, un_test, ASP_FUN_LIST, ...
                    ASP_NQP_RAD_LIST, ASP_LIST_XYINDEX, EXTRA_ARGS, UY, pars);

% %% Plot Secant Stiffness and Damping For Comparison 1

sec_k_asp = cell(size(fx_test_ASP));
diss_asp = cell(size(fx_test_ASP));


D_mindlin = cumtrapz(ux_test, fx_test);
D_mindlin = 8*D_mindlin - 4 * ux_test.*fx_test;


for ii = 1:length(fx_test_ASP)
    
    sec_k_asp{ii} = fx_test_ASP{ii}(1:5:end)./fx_test(1:5:end);
    
    tmp = cumtrapz(ux_test, fx_test_ASP{ii});
    tmp = 8*tmp - 4*ux_test.*fx_test_ASP{ii};
    
    diss_asp{ii} = tmp(1:5:end)./D_mindlin(1:5:end);
    
end

%% Plot Secant Stiffness and Dissipation

NUM_PLOT = 3;

% Secant Plot
plot_comp_secant(ux_test(1:5:end)/u_slip, sec_k_asp, ...
                '$f_t / f_{t, Mindlin}$', [0, hyst_amp], [0.95, 1.25], ...
                color_plot, ASP_NAMES_LIST, NUM_PLOT, 'Secant');


% Dissipation Plot
plot_comp_secant(ux_test(1:5:end)/u_slip, diss_asp, ...
                '$ D / D_{Mindlin}$', [0, hyst_amp], [-0.0, 1.9], ...
                color_plot, ASP_NAMES_LIST, NUM_PLOT, 'Dissipation');


%% Function for plotting secant stiffness / dissipation
function plot_comp_secant(xplot, yplot, ytext, ulimits, flimits, ...
                            color_plot, ASP_NAMES_LIST, NUM_PLOT, printName)
    font_size = 14;
    line_width = 2.5;

    figWidth = 1200;

    figure('Position', [50, 50, figWidth, 700]); % Left bottom width height 
    
    bot_height = 0.4;
    label_height = 0.1;
    
    axPosUF =   [0.13*600/figWidth, label_height,  0.74*600/figWidth, bot_height]; % Left bottom width height 
    
%     t = tiledlayout(1,1);
    ax1 = axes('Position', axPosUF);
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    
    hold on;
    
    xlabel('$u_t/u_{t,slip}$');
    ylabel(ytext);
    
    set(ax1, 'FontSize', font_size);

    xlim(ax1, ulimits);
    ylim(ax1, flimits);
    box on;

    legend_names = cell(length(yplot), 1);

    symbol_set = {'--', '-.', ':','--', ':', ':','--', '-.', ':', '--', '-.', ':'};

    for kk = 1:NUM_PLOT

        plot(xplot, yplot{kk}, symbol_set{kk}, 'LineWidth', line_width, 'Color', color_plot(kk, :));

        legend_names{kk} = ASP_NAMES_LIST{kk};

    end

    h = legend(legend_names(1:NUM_PLOT), 'NumColumns', 3);
    legPos = h.Position;
    h.Position = [0.2, ... % Left 
                    axPosUF(2)+axPosUF(4)+ 0.2, ... %bottom
                    legPos(3:4) ];%width height 
                
    set(gcf, 'Renderer', 'painters');

    drawnow;
    
    if(~isempty(printName))
        if(~isequal(printName, 'ConstN'))
            h.Visible = 'off'; % Turn off legend for saving
        end
%         print(sprintf('./OUTPUTS/MindlinTest_%s.eps', printName), '-depsc', '-r400');
    end
    
    %%% Save Legend
    h.Visible = 'on';
    ax = gca;
    ax.Position = [-1, -1, 0, 0];
    
    drawnow;
    
%     print(sprintf('./OUTPUTS/MindlinTest_FlatLegend.eps', printName), '-depsc', '-r400');

    ax.Position = axPosUF;


end

%% Function For plotting each time

function plot_comparison(printName, pars, N1, ASP_NAMES_LIST, ux_test, ...
        fx_test, fx_test_ASP, u_slip1, color_plot, NUM_PLOT, ulimits, flimits, ...
        N2, u_slip2, fn_test, Nlimits)

    font_size = 14;
    line_width = 2.5;
    
    figWidth = 1200;

    figure('Position', [50, 50, figWidth, 700]); % Left bottom width height 
    
    bot_height = 0.4;
    label_height = 0.1;
    
    axPosTime = [0.13*600/figWidth, bot_height+3*label_height, 0.74*600/figWidth, 1-bot_height-4*label_height]; % Left bottom width height 
    axPosUF =   [0.13*600/figWidth, label_height,              0.74*600/figWidth, bot_height]; % Left bottom width height 
    
%     t = tiledlayout(1,1);
    ax1 = axes('Position', axPosUF);
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    
    hold on;
    if(nargin > 12)
        xlabel('$u_t/u_{t,slip, 1}$');
        ylabel('$f_t/\mu f_{n,1}$');
    else
        xlabel('$u_t/u_{t,slip}$');
        ylabel('$f_t/\mu f_{n}$');
    end

    set(ax1, 'FontSize', font_size);

    xlim(ax1, ulimits);
    ylim(ax1, flimits);
    box on;

    legend_names = cell(length(fx_test_ASP)+1, 1);

    plot(ux_test/u_slip1, fx_test/pars.mu/N1, 'k-', 'LineWidth', 2.5);
    legend_names{1} = 'Analytical';

    symbol_set = {'--', '-.', ':','--', ':', ':','--', '-.', ':', '--', '-.', ':'};

    for kk = 1:NUM_PLOT

        plot(ux_test/u_slip1, fx_test_ASP{kk}/pars.mu/N1, symbol_set{kk}, 'LineWidth', line_width, 'Color', color_plot(kk, :));

        legend_names{kk+1} = ASP_NAMES_LIST{kk};

    end

    h = legend(legend_names(1:NUM_PLOT+1));
    legPos = h.Position;
    h.Position = [2*axPosUF(1)+axPosUF(3), ... % Left 
                    axPosUF(2)+axPosUF(4)/2-legPos(4)/2, ... %bottom
                    legPos(3:4) ];%width height 

    if(nargin > 12)
        
        ax2 = axes('Position', ax1.Position);
        
        ax2.XAxisLocation = 'top';
        ax2.YAxisLocation = 'right';
        ax2.Color = 'none';
%         ax2.Xcolor = 'k';
%         ax2.Ycolor = 'k';
        
        xlabel(ax2, '$u_t/u_{t,slip,2}$');
        ylabel(ax2, '$f_t/\mu f_{n,2}$');

        ax1.Box = 'off';
        ax2.Box = 'off';
        
        
        xlim(ax2, ulimits*u_slip1/u_slip2);
        ylim(ax2, flimits*N1/N2);
        
        set(ax2, 'FontSize', font_size);
    end

    if(nargin > 14)
        %%%%%%%%% Plot Some History for visibility
        dux = [0, abs(ux_test(2:end)-ux_test(1:end-1))];
        cum_ux_test = cumsum(dux);

        ax3 = axes('Position', axPosTime);
        ax3.XColor = 'k';
        ax3.YColor = 'k';

        hold on;
        xlabel(ax3, '$u_t/u_{t,slip, 1}$');
        ylabel(ax3, '$f_{n}/\mu f_{n,1}$');
        xlim(ax3, ulimits);
        ylim(ax3, Nlimits);

        plot(cum_ux_test/u_slip1, fn_test/N1, 'k', 'LineWidth', line_width);

        set(ax3, 'FontSize', font_size);
        
        
        ax4 = axes('Position', ax3.Position);
        ax4.XAxisLocation = 'top';
        ax4.YAxisLocation = 'right';
        ax4.Color = 'none';
        xlabel(ax4, '$u_t/u_{t,slip,2}$');
        ylabel(ax4, '$f_{n}/\mu f_{n,2}$');
        ax3.Box = 'off';
        ax4.Box = 'off';
        
        xlim(ax4, ulimits*u_slip1/u_slip2);
        ylim(ax4, Nlimits*N1/N2);
        
        set(ax4, 'FontSize', font_size);
    end
    
    set(gcf, 'Renderer', 'painters');

    drawnow;
    
    if(~isempty(printName))
        if(~isequal(printName, 'ConstN'))
            h.Visible = 'off'; % Turn off legend for saving
        end
%         print(sprintf('./OUTPUTS/MindlinTest_%s.eps', printName), '-depsc', '-r400');
    end

end

%% Function for tracing history 

function [fx_test_ASP] = generate_data(ux_test, un_test, ASP_FUN_LIST, ...
                ASP_NQP_RAD_LIST, ASP_LIST_XYINDEX, EXTRA_ARGS, UY, pars)


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


%%%%%%%%%%%%LOOP TO UPDATE GIVEN LIST OF ASPERITY MODELS

%Initialize to the first normal load for the ASP_FUN as well 
for kk = 1:length(ASP_FUN_LIST)

    uxyn_ASP = [0, 0, un_test(1)];

    if(EXTRA_ARGS(kk) )
        [fxyn, dfxynduxyn, rq0{kk}, tx0{kk}, ty0{kk}] ...
            = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk}, 0, 0, 0);
    else
        [fxyn, dfxynduxyn, rq0{kk}, tx0{kk}, ty0{kk}] ...
            = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk});
    end
    
    uxyn0_ASP{kk} = uxyn_ASP;

end

for jj = 1:length(ux_test)

    %Loop over asperity functions
    for kk = 1:length(ASP_FUN_LIST)

        xynum_asp = ASP_LIST_XYINDEX(kk);

        %Asperity function test. 
        uxyn_ASP = [0, 0, 0];
        uxyn_ASP(1,xynum_asp) = ux_test(jj);
        uxyn_ASP(1,3-xynum_asp) = UY;
        uxyn_ASP(1, 3)    = un_test(jj);

        
        if(EXTRA_ARGS(kk) )
            [fxyn, dfxynduxyn, rq0{kk}, tx0{kk}, ty0{kk}] ...
                = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk}, 0, 0, 0);
        else
            [fxyn, dfxynduxyn, rq0{kk}, tx0{kk}, ty0{kk}] ...
                = ASP_FUN_LIST{kk}(pars, uxyn_ASP, uxyn0_ASP{kk}, rq0{kk}, tx0{kk}, ty0{kk}, rq{kk}, wq{kk});
        end
        
        uxyn0_ASP{kk} = uxyn_ASP;

        fx_test_ASP{kk}(jj) = fxyn(xynum_asp);

    end

end


end

%% Function for integrating increasing T and N
function [dNUx_dT] = MINDLIN_INCREASE_DOT(T, NUx, pars)
% Function to use for ode45 integration of increasing T and N
% Time Variable (pseudo):
%   T - Tangential Force
% States:
%   N  - normal load
%   Ux - Tangential Displacement
    
    dNUx_dT = zeros(size(NUx));
    
    dNUx_dT(1) = 1/pars.mu; % The limitting case where the compliance is valid
    
    a = (3*NUx(1)*pars.R/4/pars.Estar)^(1/3); % pg93, 4.22 from Johnson, 1985
    
    dNUx_dT(2) = 1/(8*pars.Gstar*a); % Mindlin and Deresiewicz, 1953
    
end


%% Function for integrating increasing T and decreasing N
function [dNUx_dT] = MINDLIN_INCREASET_DECREASEN_DOT(T, NUx, pars, dNdT)
% Function to use for ode45 integration of increasing T and N
% Time Variable (pseudo):
%   T - Tangential Force
% States:
%   N  - normal load
%   Ux - Tangential Displacement
    
    dNUx_dT = zeros(size(NUx));
    
%     dNUx_dT(1) = -1/pars.mu; 
    dNUx_dT(1) = dNdT; 
    
    a = (3*NUx(1)*pars.R/4/pars.Estar)^(1/3); % pg93, 4.22 from Johnson, 1985
    
    dNUx_dT(2) = 1/(8*pars.Gstar*a)*(-pars.mu* abs(dNUx_dT(1)) ...
        + (1 + pars.mu* abs(dNUx_dT(1)))*(1 - T / (pars.mu * NUx(1)))^(-1/3) ); % Mindlin and Deresiewicz, 1953 - Eqn 30
    
end