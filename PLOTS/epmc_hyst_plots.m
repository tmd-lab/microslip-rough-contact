% Script for plotting hysteresis loops generated with epmc

clear;

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


set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');


colors_plots = DISTINGUISHABLE_COLORS(13, 'w');

%% Load File information

plot_indices = 1;% [1, 5, 32];

plot_set = 2;

desired_freq_shift = 3.5; % Hz - for choosing points consistently
clims_input = [0, 0.08];


switch plot_set
    

    case 2
        %%%%%%% PAPER: Elastic, mu = 0.03
        run_name = sprintf('plots%u', plot_set);
                        
        tmp = load('../EPMC_SIMS/Results/PAPER/Run10/epmc_paper_run10_v2_iter.mat');

        load('../EPMC_SIMS/Results/PAPER/Run10/epmc_paper_run10_v2_pre.mat');

        % % 232 ZTE Model        
        load(sprintf('../FJSIMS/ROMS/ROM_U_%uELS', 232), 'M', 'K', 'R', 'Fv', 'L', 'MESH');

        [~, ind] = min(abs(abs(tmp.U(end-2, :)/2/pi - tmp.U(end-2, 1)/2/pi) - desired_freq_shift));
        
%         % Switch to MESH Class for sake of plots
%         Nq = 1;
%         MESH = MESH2D(MESH.Nds, 3, [], MESH.Quad, Nq);

        Uwxa = tmp.U(:, ind); 

        plot_indices = [1, 9, 10, 91, 75, 51, 47, 60, 52, 22, 16, 8]; 
        paper_plot = [1, 4, 7]; % Indices of plot set to save for use in paper. 

    case 3
        %%%%%%% PAPER: Plastic, mu = 0.03
        run_name = sprintf('plots%u', plot_set);
                        
        tmp = load('../EPMC_SIMS/Results/PAPER/Run14/epmc_paper_run14_v2_iter.mat');

        load('../EPMC_SIMS/Results/PAPER/Run14/epmc_paper_run14_v2_pre.mat');

        % % 232 ZTE Model        
        load(sprintf('../FJSIMS/ROMS/ROM_U_%uELS', 232), 'M', 'K', 'R', 'Fv', 'L', 'MESH');

        [~, ind] = min(abs(abs(tmp.U(end-2, :)/2/pi - tmp.U(end-2, 1)/2/pi) - desired_freq_shift));
        
%         ind = size(tmp.U, 2);

        
%         % Switch to MESH Class for sake of plots
%         Nq = 1;
%         MESH = MESH2D(MESH.Nds, 3, [], MESH.Quad, Nq);

        Uwxa = tmp.U(:, ind); 


        plot_indices = [5, 9, 10, 91, 75, 51, 47, 60, 52, 22, 16, 8];    
        paper_plot = [1, 4, 11]; % Indices of plot set to save for use in paper. 

    case 4
        %%%%%%% PAPER: Plastic, mu = CEB
        run_name = sprintf('plots%u', plot_set);
                        
        tmp = load('../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_iter.mat');
        
        load('../EPMC_SIMS/Results/PAPER/Run11/epmc_paper_run11_v2_pre.mat');

        % % 232 ZTE Model        
        load(sprintf('../FJSIMS/ROMS/ROM_U_%uELS', 232), 'M', 'K', 'R', 'Fv', 'L', 'MESH');

        [~, ind] = min(abs(abs(tmp.U(end-2, :)/2/pi - tmp.U(end-2, 1)/2/pi) - desired_freq_shift));
        
%         % Switch to MESH Class for sake of plots
%         Nq = 1;
%         MESH = MESH2D(MESH.Nds, 3, [], MESH.Quad, Nq);

        Uwxa = tmp.U(:, ind); 


%         plot_indices = [9, 10, 91, 75, 51, 60, 16, 15, 37, 205, 209, 212];
        plot_indices = [9, 209, 205, 91, 75, 51, 212, 60, 37, 16, 15, 10]; %Sorted
        paper_plot = [1, 2, 4, 5, 7, 9, 10, 11, 12]; % Indices of plot set to save for use in paper. 

        clims_input = [0, 2.75];
    
        
    case 5
        %%%%%%% PAPER: Plastic, mu = 0.03, Flat Interface
        run_name = sprintf('plots%u', plot_set);
                        
        tmp = load('../EPMC_SIMS/Results/PAPER/Run15/epmc_paper_run15_iter.mat');
                        
        load('../EPMC_SIMS/Results/PAPER/Run15/epmc_paper_run15_pre.mat');

        % % 232 ZTE Model        
        load(sprintf('../FJSIMS/ROMS/ROM_U_%uELS', 232), 'M', 'K', 'R', 'Fv', 'L', 'MESH');

        [~, ind] = min(abs(abs(tmp.U(end-2, :)/2/pi - tmp.U(end-2, 1)/2/pi) - desired_freq_shift));
                
        
%         ind = size(tmp.U, 2);
        
%         % Switch to MESH Class for sake of plots
%         Nq = 1;
%         MESH = MESH2D(MESH.Nds, 3, [], MESH.Quad, Nq);

        Uwxa = tmp.U(:, ind); 

        plot_indices = [9, 17, 26, 91, 75, 51, 50, 47, 22, 23, 5, 11];    
        paper_plot = [11, 1, 4]; %[1, 2, 4, 6, 10, 11]; % Indices of plot set to save for use in paper. 

        clims_input = [0, 0.024];

    otherwise
        error('Undefined plot set');

end

ActualShift = (tmp.U(end-2, ind) - tmp.U(end-2, 1))/2/pi



%% Produce Hysteresis Loops:

[unlt_quad, fnlt_quad] = GM.EPMC_HYST(Uwxa, Fl, h, Nt, 1e-6);

%% Calculate Energy flux

diss_flux = zeros(size(unlt_quad));

for ii = 1:length(diss_flux)
    
    ux = [unlt_quad{ii}(:, 1); unlt_quad{ii}(1, 1)];
    uy = [unlt_quad{ii}(:, 2); unlt_quad{ii}(1, 2)];
    tx = [fnlt_quad{ii}(:, 1); fnlt_quad{ii}(1, 1)];
    ty = [fnlt_quad{ii}(:, 2); fnlt_quad{ii}(1, 2)];
    
    diss_flux(ii) = trapz(ux, tx) + trapz(uy, ty);
    
end

%% Reorder for Paper


% %% ZTE Quadrature Matrices
No   = 1;                  % Number of GLQ integration points in each direction
[Qm,Tm] = ZTE_ND2QP(MESH,No);

Nds_coords = Qm(plot_indices(paper_plot), :)*MESH.Nds(:, 1:2);
Nds_coords(:, 2) = -Nds_coords(:, 2); % Start numbering in the top left
Nds_coords = round(Nds_coords, 4);
[~, inds] = sortrows(Nds_coords);

paper_plot = paper_plot(inds);

%% Plot Loops

line_width = 3;
font_size = 14;


disp_size = [100, 700, 560, 420]; %[x, y, width, height]
loop_size = [100, 100, 560, 420];
disp_size_subs = [700, 700, 560, 600]; %[x, y, width, height]

xyn_colors = colors_plots([1,2,3], :);
xyn_style = {'.-', '--', '-'};

% Generally should be set in load block
% plot_indices = [1, 8, 9, 10, 16, 28, 22, 60, 52, 47, 51, 75];
plot_letters = {'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p'};

xyhystloops = cell(length(plot_indices), 2);

for pind = plot_indices
% for pind = 1:length(unlt_quad)
    
    % Plot what the displacements are doing
    figure('Position', disp_size);
    hold on;
    xlabel('$t/T$');
    ylabel('$u_i$ [m]');
    
    t_T = linspace(0, 1, length(unlt_quad{pind}(:, 1)));% Normalized time
    
    plot(t_T, unlt_quad{pind}(:, 1), xyn_style{1}, 'Color', xyn_colors(1, :), 'LineWidth', line_width);
    plot(t_T, unlt_quad{pind}(:, 2), xyn_style{2}, 'Color', xyn_colors(2, :), 'LineWidth', line_width);
    plot(t_T, unlt_quad{pind}(:, 3), xyn_style{3}, 'Color', xyn_colors(3, :), 'LineWidth', line_width);
    
    legend('$u_x$', '$u_y$', '$u_n$', 'Location', 'ne');
    
    set(gcf, 'Renderer', 'painters');

    set(gca,'FontSize',font_size)

    drawnow;
    
%     print(sprintf('./OUTPUTS/disp_hist_DOF%u_%s.eps', pind, run_name), '-depsc', '-r400');
%     print(sprintf('./OUTPUTS/disp_hist_DOF%u_%s.svg', pind, run_name), '-dsvg', '-r400');
    
    title(sprintf('DOF=%u', pind));
    
    drawnow; 
    
    % Plot what the displacements are doing
    figure('Position', disp_size_subs);
    ax = subplot(3,1,1);

    for dir = 1:3
        subplot(3,1,dir);
        hold on;
        ylabel(sprintf('$u_%s$ [m]', 'x'*(dir == 1) + 'y'*(dir == 2) + 'n'*(dir == 3)));

        t_T = linspace(0, 1, length(unlt_quad{pind}(:, 1)));% Normalized time

        plot(t_T, unlt_quad{pind}(:, dir), xyn_style{dir}, 'Color', xyn_colors(dir, :), 'LineWidth', line_width);

        set(gcf, 'Renderer', 'painters');

        set(gca,'FontSize',font_size)
        subheight = 0.18;
        set(gca, 'Position', [.2, 0.05+(4-dir)*0.1 + (subheight*(dir<3) + subheight*(dir < 2)), 0.7, subheight]); %[left bottom width height]
    end
    
    xlabel('$t/T$');
    
    title(ax, sprintf('DOF=%u', pind));

    drawnow;
    
%     print(sprintf('./OUTPUTS/TMD/disp_hist_split_DOF%u_%s.eps', pind, run_name), '-depsc', '-r400');
%     print(sprintf('./OUTPUTS/TMD/disp_hist_split_DOF%u_%s.svg', pind, run_name), '-dsvg', '-r400');
%     subplot(3,1,1);
%     title(sprintf('DOF=%u', pind));
    drawnow;

    
    % Plot the Effect on Forces
    for pxy = 1:2

        if(pxy == 1)
            xysym = 'x';
        elseif(pxy == 2)
            xysym = 'y';
        elseif(pxy == 3)
            xysym = 'n';
        else
            xysym = '';
        end
        
        drawnow;
        
        figure('Position', loop_size + loop_size(3)*1.5*[1, 0, 0, 0]*(pxy==2));
        hold on;
        xlabel(sprintf('$u_%s$ [m]', xysym));
        ylabel(sprintf('$t_%s$ [Pa]', xysym));


        plot(unlt_quad{pind}(:, pxy), fnlt_quad{pind}(:, pxy), xyn_style{pxy}, 'Color', xyn_colors(pxy, :), 'LineWidth', line_width);

%         legend('$x$', '$y$', '$n$', 'Location', 'ne');

        set(gcf, 'Renderer', 'painters');
        
        loopfig = gcf;

        set(gca,'FontSize',font_size)

        ax = gca;
        xyhystloops{pind, pxy} = ax;

        drawnow;
        
        ax.YAxis.Exponent = max(floor(log10(abs(ax.YLim))));

        drawnow;
        
        plot_num = find(plot_indices == pind);
        if(pxy == 1 && ismember(plot_num, paper_plot))
            % Save the output figure
            letter_ind = find(paper_plot == plot_num);
            
            figname = sprintf('OUTPUTS/LoopSet%u_Sub%s.eps', plot_set, plot_letters{letter_ind});
            
%             print(loopfig, figname, '-depsc', '-r400');

        end
        
        title(ax, sprintf('DOF=%u, Direction=%s', pind, xysym));
        drawnow;
        
        if(pxy == 1)
%             print(sprintf('./OUTPUTS/TMD/hyst_Dir%s_DOF%u_%s.eps', xysym, pind, run_name), '-depsc', '-r400');
%             print(sprintf('./OUTPUTS/TMD/hyst_Dir%s_DOF%u_%s.svg', xysym, pind, run_name), '-dsvg', '-r400');
        end
    
    end
    
    
end




%% Plot of the Locations of Quadrature Points

% %% ZTE Quadrature Matrices
No   = 1;                  % Number of GLQ integration points in each direction
[Qm,Tm] = ZTE_ND2QP(MESH,No);

font_size = 20;


figure('Position', [200, 400, 1000, 400]);

SHOW2DMESH(MESH.Nds,MESH.Tri,MESH.Quad,'w',-1,-100)  % plot mesh without node numbers
hold on; 

drawnow;

xoffset = 0.0005;
yoffset = 0.0005;

for pind = plot_indices
    
    
    plot_num = find(plot_indices == pind);
    
    if(ismember(plot_num, paper_plot))

        letter_ind = find(paper_plot == plot_num);
        
    
        % Letter label
        text(Qm(pind,:)*MESH.Nds(:,1) + xoffset, Qm(pind,:)*MESH.Nds(:,2) + yoffset, plot_letters{letter_ind}, 'FontSize', font_size)%, 'BackgroundColor', [0.9, 0.9,0.9])
        
        plot(Qm(pind,:)*MESH.Nds(:,1), Qm(pind,:)*MESH.Nds(:,2), 'r.', 'MarkerSize', 25)
    end
    
    % Number labels (all)
%     plot(Qm(pind,:)*MESH.Nds(:,1), Qm(pind,:)*MESH.Nds(:,2), 'r.', 'MarkerSize', 25)
%     text(Qm(pind,:)*MESH.Nds(:,1) + xoffset, Qm(pind,:)*MESH.Nds(:,2) + yoffset, num2str(pind), 'FontSize', font_size)
    
end
% Qm(1,:)

ax = gca;


set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('$x$');
ylabel('$y$');
ylim(ax, [min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))]);

set(gca,'FontSize',font_size)


% spy(Qm)
% spy(Tm)
% clf()
% SHOW2DMESH(MESH.Nds,MESH.Tri,MESH.Quad,'w',-1,-100)  % plot mesh without node numbers
% SHOW2DMESH(MESH.Nds,MESH.Tri,MESH.Quad,'w',-1,1)  % plot mesh with node numbers
% SHOW2DMESH(MESH.Nds,[],MESH.Quad,'b',-1,-100)  % plot mesh without node numbers

axis equal;

xlim(ax, [min(MESH.Nds(:, 1)), max(MESH.Nds(:, 1))]);
ylim(ax, [min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))]);

set(gcf, 'Renderer', 'painters');




drawnow;

% print(sprintf('OUTPUTS/LoopSet%u_locs.eps', plot_set), '-depsc', '-r400');


% Rerun y limit before saving
% print('./OUTPUTS/TMD/quad_labeled_mesh152.svg', '-dsvg', '-r400');

%% Mesh Plot without any details

%% Plot of the Locations of Quadrature Points

% %% ZTE Quadrature Matrices
No   = 1;                  % Number of GLQ integration points in each direction
[Qm,Tm] = ZTE_ND2QP(MESH,No);

%%%%%%%%%%%
% Reorder Dissipate Fluxes Based on Qm
% Probably unnecessary, but Just in case Qm is out of order for some
% reason, this fixes it (at least I double checked for quad only meshes)
%%%%%%%%%%
elem_order_diss_flux = zeros(size(diss_flux));
mapping_inds = zeros(size(diss_flux));

QuadNds_sorted = [MESH.Quad(:, 1), sort(MESH.Quad(:, 2:end), 2)];
TriNds_sorted = [MESH.Tri(:, 1), sort(MESH.Tri(:, 2:end), 2)];

for ii = 1:length(diss_flux)
    
    Nds = find(Qm(ii, :));
    
    if(length(Nds) == 3)
        [~, loc] = ismember(TriNds_sorted(:, 2:end), Nds, 'rows');
        elem_num = QuadNds_sorted(loc==1, 1);
    elseif(length(Nds) == 4)
        [~, loc] = ismember(QuadNds_sorted(:, 2:end), Nds, 'rows');  
        elem_num = QuadNds_sorted(loc==1, 1);      
    end
    
    elem_order_diss_flux(elem_num) = diss_flux(ii);
    mapping_inds(ii) = elem_num;
    
end



%%%%%%%%%%%
% Actually Plot Mesh
%%%%%%%%%%%

font_size = 20;


figure('Position', [200, 400, 1000, 500]);


colormap jet;


% Cface = 'w'; % White mesh;
% Cface = ones(MESH.Ne, 1)*[0, 1, 0]; % Solid Mesh
Cface = elem_order_diss_flux; % color map properties
checkfin = 0;

SHOW2DMESH(MESH.Nds,MESH.Tri,MESH.Quad,Cface,-1,-100, checkfin)  % plot mesh without node numbers
hold on; 

color_bar = colorbar('northoutside');
color_bar.Label.String = 'Dissipative Flux [J/m$^2$]';
color_bar.Label.Interpreter = 'latex';
color_bar.FontSize = font_size;

% caxis(clims_input);


drawnow;

% Qm(1,:)

ax = gca;

xlim(ax, [min(MESH.Nds(:, 1)), max(MESH.Nds(:, 1))]);
ylim(ax, [min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))]);


set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('$x$');
ylabel('$y$');
ylim(ax, [min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))]);

set(gca,'FontSize',font_size)

meshplot = ax;


axis equal;
set(gcf, 'Renderer', 'painters');

xlim(ax, [min(MESH.Nds(:, 1)), max(MESH.Nds(:, 1))]);
ylim(ax, [min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))]);

drawnow;

% print(sprintf('OUTPUTS/LoopSet%u_flux.eps', plot_set), '-depsc', '-r400');


%% Create the Mesh surrounded by figures

% 3 x 5 grid of hysteresis loops with inner 1 x 3 taken by mesh

xgap = 0.05; % Side to side
ygap = 0.09; % Up and down

font_size = 12;


subheight = (1 - ygap*4)/3;
subwidth = (1 - xgap*6)/5;

figHeight_Width = subwidth / subheight;

monitor_size = get(0, 'MonitorPositions');

figwidth = 0.9*min(monitor_size(3), monitor_size(4)/figHeight_Width);
figheight = figwidth*figHeight_Width;

% Set size if the monitor is large enough
if(figwidth > 1500)
    figwidth = 1500;
    figheight = round(figwidth*figHeight_Width);
end

figure('Position', [(monitor_size(3) - figwidth )/ 2, (monitor_size(4) - figheight )/ 2, figwidth, figheight]); %[left bottom width height]
set(gcf, 'Renderer', 'painters');


pxy = 1; % 1 for x, 2 for y



plot_indices = plot_indices(1:12);  % Manual reordering for nicely setup

lineCoords = zeros(12, 2);

% Start at the top left and go clockwise

for ii = 1:12
    
    % Zero indexed row and column positions of figures
    subrow = 2.*(ii<=5) + 1.*(mod(ii, 6) == 0);%+0.*(6 < ii < 12)
    
    subcol = (ii-1).*(subrow == 2) ... % Top row
                + 4.*(ii == 6) ... % Right plot
                + 0.*(ii == 12) ... % Left plot
                + (11-ii).*(subrow == 0); % Bottom
    
    subpos = [subcol*(subwidth + xgap) + xgap, ...% Left
              subrow*(subheight + ygap) + ygap, ... % Bottom
              subwidth, subheight];
            
    curraxis = subplot('Position', subpos); % [left bottom width height]. 
    hold on;
%     title(sprintf('Subplot %u', ii));
    xlabel('$u_x$ [m]');
    ylabel('$t_x$ [Pa]');
    
    ax1Chil = xyhystloops{plot_indices(ii), pxy}.Children; 
    copyobj(ax1Chil, curraxis)
    
    curraxis.YAxis.Exponent = max(floor(log10(abs(ylim))));
    
    set(gca,'FontSize',font_size)

    % Set Coordinates for lines to the mesh
    lineCoords(ii, 1) = subpos(1) + 0.5*subwidth + 0.4*subwidth*(subcol == 0) - 0.4*subwidth*(subcol == 4);
    lineCoords(ii, 2) = subpos(2) + 0.5*subheight + 0.4*subheight*(subrow == 0) - 0.4*subheight*(subrow == 2);
    
end

% Insert the Mesh


meshHeight = subheight;
meshWidth = 3*subwidth + 2*xgap;

% Correct Aspect Ratio
aspectRatio = range(MESH.Nds(:, 2))/range(MESH.Nds(:, 1));
meshHeight = min(meshHeight*figheight, meshWidth*aspectRatio*figwidth)/figheight;
meshWidth = (meshHeight*figheight)/aspectRatio/figwidth;

subpos = [0.5 - meshWidth/2, 0.5 - meshHeight/2, meshWidth, meshHeight];

curraxis = subplot('Position', subpos); % [left bottom width height]. 

ax1Chil = meshplot.Children; 
copyobj(ax1Chil, curraxis)

%%%%%%%%
% Color bar
colormap jet
color_bar = colorbar('northoutside');
color_bar.Label.String = 'Dissipative Flux [J/m$^2$]';
color_bar.Label.Interpreter = 'latex';
color_bar.FontSize = font_size;

% caxis(clims_input);

curraxis.Position = subpos + [0, -0.04, 0, 0];
%%%%%%%

xlim([min(MESH.Nds(:, 1)), max(MESH.Nds(:, 1))]);
ylim([min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))]);
hold on;
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('$x$');
ylabel('$y$');
set(gca,'FontSize',font_size)


marker_size = 15;
line_width = 2;

for ii = 1:length(plot_indices)
    
    pind = plot_indices(ii);
    
    xyCoords = Qm(pind,:)*MESH.Nds;
    
    [xf, yf] = DS2NFU(xyCoords(1), xyCoords(2));
    
    plot(xyCoords(1), xyCoords(2), 'r.', 'MarkerSize', marker_size)
    
    connectline = annotation('line',[xf, lineCoords(ii, 1)],[yf, lineCoords(ii, 2)],'color','r', 'LineWidth', line_width);

    uistack(connectline, 'bottom');
    
end


% How to draw lines between subplots: https://stackoverflow.com/questions/3635411/draw-line-between-two-subplots

% annotation('line',[0.25, 0.75],[0.5, 0.5],'color','k', 'LineWidth', 3);



