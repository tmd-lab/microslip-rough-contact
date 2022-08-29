% Script For Plotting Meshes

clear;

set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

addpath('../ROUTINES/FEM/'); % ROM/FEM Info


%% Uniform Remeshed Interface

load(sprintf('../FJSIMS/ROMS/ROM_U_%uELS', 232), 'M', 'K', 'R', 'Fv', 'L', 'MESH');

% Switch to MESH Class for sake of plots
Nq = 1;
MESH = MESH2D(MESH.Nds, 3, [], MESH.Quad, Nq);

%% Plot Mesh

figure;

MESH.SHOWFIELD2D(zeros(MESH.Nn, 2))  % plot mesh without node numbers
hold on;

xRange = [min(MESH.Nds(:, 1)), max(MESH.Nds(:, 1))];
yRange = [min(MESH.Nds(:, 2)), max(MESH.Nds(:, 2))];

axis equal;

xlim(xRange);
ylim(yRange);

ax = gca;
ax.XTick = [];
ax.YTick = [];

set(gcf, 'Renderer', 'painters');

drawnow;

% print(sprintf('./OUTPUTS/urom%u_mesh.eps', MESH.Ne), '-depsc', '-r600');

