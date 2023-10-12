%% Save out test ouputs to be used to verify a python implementation 
% asperity level

clear; 


addpath('../ROUTINES/FRIC_MODELS')
addpath('../ROUTINES/FRIC_MODELS/ASP_FUN/')
addpath('../ROUTINES/FRIC_MODELS/MINDLIN/')
addpath('../ROUTINES/FRIC_MODELS/PLASTIC/')

addpath('../ROUTINES/')

ASP_FUN_PRE = @SPHERE_PLASTIC_PRE; Nqp_radius = 1;
ASP_FUN_TAN = @SPHERE_PL_TAN_DECOUPLE_PL_NORM; Nqp_radius = 1;


%% Baseline Parameters That Should Be Tested

initial_pars.E = 1.92e+11; %304L from Traction paper
initial_pars.nu = 0.30; %304L taken from Traction paper
initial_pars.R = 1.4e-3;

initial_pars.Et = 620e6;
initial_pars.Sys = 330e6; %Pa

%% Alternative Parameters to Capture Full Changes

alt_pars.E = 2.5e+11; %304L from Traction paper
alt_pars.nu = 0.2; %304L taken from Traction paper
alt_pars.R = 5e-3;

alt_pars.Et = 800e6;
alt_pars.Sys = 200e6; %Pa

%% Another Alternative pars 
% Et = 0 has the potential to break gradients, so try it here.

alt2_pars = initial_pars;
alt2_pars.Et = 0;

%% Full Pars

full_pars = initialize_pars(initial_pars);

full_pars_alt = initialize_pars(alt_pars);

full_pars_alt2 = initialize_pars(alt2_pars);


% Save everything into a nice format for looping
pars_list = {full_pars, full_pars, full_pars_alt, full_pars_alt2};

input_list = {initial_pars, initial_pars, alt_pars, alt2_pars};

un_shift = {0, 5*full_pars.delta_y, 0, 0};

%% Calculate Results

res_list = cell(size(pars_list));

for i = 1:length(pars_list)
    
    [normal_disp, a_list, fn_list] = produce_results(pars_list{i}, ASP_FUN_PRE, un_shift{i});
    
    results.normal_disp = normal_disp;
    results.contact_radius = a_list;
    results.normal_force = fn_list;
    
    res_list{i} = results;
    
end


%% Write Out a Yaml File

fileID = fopen('normal_asperity.yaml','w');


% fprintf(fileID,'%6s %12s\n','x','exp(x)');
% fprintf(fileID,'%6.2f %12.8f\n',A);


% Start with the Input parameters
fn = fieldnames(input_list{1});
%loop through the fields
for i = 1:numel(fn)
    field_vals = zeros(size(input_list));
    for j = 1:length(input_list)
        field_vals(j) = input_list{j}.(fn{i});
    end
    
    fprintf(fileID, '%s : %s \n', fn{i}, mat2str_comma(field_vals));
end

% Start with the Output Results
fn = fieldnames(res_list{1});
%loop through the fields
for i = 1:numel(fn)
    fprintf(fileID, '%s : [', fn{i});
    field_vals = res_list{1}.(fn{i});
    fprintf(fileID, '%s', mat2str_comma(field_vals));
    for j = 2:length(res_list)
        field_vals = res_list{j}.(fn{i});
        fprintf(fileID, ', %s', mat2str_comma(field_vals));
    end
    
    fprintf(fileID, '] \n');
end

fclose(fileID);

%% Matrix to Str with comma: 
% n = [12345 6789 10234 3452]*1e-3;
% 
% str = mat2str_comma(n)
    
function str = mat2str_comma(mat)
    
    str = sprintf('%.12e,' , mat);
    str = str(1:end-1);% strip final comma

    str = "[" + str + "]";
end
%% Produce Test Results Function

function [normal_disp, a_list, fn_list] = produce_results(full_pars, ASP_FUN_PRE, un_shift)

    un_max = 3.0*full_pars.delta_y;
    uxyn = [0, 0, un_max];
    uxyn0 = [0, 0, 0];
    rq0 = 0;
    fx0 = 0;
    fy0 = 0;
    rq = 0;
    wq = 1;
    deltam = 0;
    Fm = 0; 
    am = 0;


    % Determine the unloading point
    [fxyn, dfxynduxyn, ~, ~, ~, deltam, Fm, am, ...
        dfduxyn0, dfdfxy0, dfddeltam, a, dadun, Sys, ...
        dSysdun, daddeltam, dSysddeltam, deltabar, ddeltabar_ddeltam, ...
        Rebar, dRebar_ddeltam] ...
            = ASP_FUN_PRE(full_pars, uxyn, uxyn0, rq0, fx0, fy0, rq, wq, deltam, Fm, am);


    [fxyn, dfxynduxyn, ~, ~, ~, deltam, Fm, am, ...
        dfduxyn0, dfdfxy0, dfddeltam, a, dadun, Sys, ...
        dSysdun, daddeltam, dSysddeltam, deltabar, ddeltabar_ddeltam, ...
        Rebar, dRebar_ddeltam] ...
            = ASP_FUN_PRE(full_pars, uxyn, uxyn, rq0, fx0, fy0, rq, wq, deltam, Fm, am);

    normal_disp = [-full_pars.delta_y, 0.5*full_pars.delta_y, ...
                   1.0*full_pars.delta_y, 0.5*full_pars.delta_y, ...
                   un_max, deltabar + 0.5*(un_max - deltabar), ...
                   0.5*deltabar];

    normal_disp = normal_disp + un_shift;
    
    a_list = zeros(size(normal_disp));
    fn_list = zeros(size(normal_disp));


    deltam = 0;
    Fm = 0; 
    am = 0;
    uxyn0 = [0, 0, 0];

    for i = 1:length(normal_disp)

        uxyn = [0, 0, normal_disp(i)];

        % Determine the unloading point
        [fxyn, dfxynduxyn, ~, ~, ~, deltam, Fm, am, ...
            dfduxyn0, dfdfxy0, dfddeltam, a, dadun, Sys, ...
            dSysdun, daddeltam, dSysddeltam, deltabar, ddeltabar_ddeltam, ...
            Rebar, dRebar_ddeltam] ...
                = ASP_FUN_PRE(full_pars, uxyn, uxyn0, rq0, fx0, fy0, rq, wq, deltam, Fm, am); 

        fn_list(i) = fxyn(3);
        if fxyn(3) ~= 0
            a_list(i) = a;
        end

        uxyn0 = uxyn;
    end

end

%% Function for initializing properties

function full_pars = initialize_pars(initial_pars)

    nu = initial_pars.nu;
    E = initial_pars.E;

    G = E / 2/ (1 + nu);%304L approx Shear modulus

    full_pars.mu     = 0.05; %may want to do a function of displacement - see Eriten's work
    % pars.R      = (1/Rx + 1/Ry)^(-1);%, 1e-4 is about a maximum
    full_pars.E      = E;
    full_pars.Estar  = E/2/(1 - nu^2);
    full_pars.nu     = nu;
    full_pars.G      = G;
    full_pars.Gstar  = G/2/(2 - nu);

    full_pars.R = initial_pars.R;
    full_pars.Re = initial_pars.R;

    Sys = initial_pars.Sys; %Pa
    Et = initial_pars.Et; 

    C = 1.295*exp(0.736*full_pars.nu);

    delta_y1s = (pi*C*Sys/(2*(2*full_pars.Estar)))^2*(2*full_pars.R); %displacement of one sphere against a rigid flat to cause yielding.

    % Store all of these results as system parameters to prevent recalculation:
    full_pars.delta_y = delta_y1s*2;
    full_pars.Sys = Sys;

    %tangent stiffness after yield.
    full_pars.Et = Et;
    
    full_pars.sliptype = 1;

end










