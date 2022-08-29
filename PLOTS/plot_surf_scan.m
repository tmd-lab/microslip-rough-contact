% Plot an example of the surface scan

clear;

set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');


%% Surface To Plot

load('../SURFACE/SCANS/LongTermWear/R05A_Before_R1.mat');


%% Loop over patches and create one set of data

patch_data = cell(size(RawDatabyElem, 1), 1);

for patch = 1:size(RawDatabyElem, 1)
    
    patch_data{patch} = cell2mat(cellfun(@(c) c, RawDatabyElem(patch, :), 'UniformOutput', false)');
    
end

%% Plot Patch Data

figure('Position', [250, 250, 600, 400]);

patch_data(7:end) = patch_data(end:-1:7); 

red_plot = true; % Plot only a reduced set of points to check before doing everything.

zunit_convert = 1e6;
xyunit_convert = 1e3;

for patch = 1:length(patch_data)
    
    if(red_plot)
        Nplot = 100;
        step = round(size(patch_data{patch}, 1)/Nplot);
        scatter3(patch_data{patch}(1:step:end, 1)*xyunit_convert, ...
            patch_data{patch}(1:step:end, 2)*xyunit_convert, ... 
            patch_data{patch}(1:step:end, 3)*zunit_convert);
    else
        scatter3(patch_data{patch}(:, 1)*xyunit_convert, ...
                    patch_data{patch}(:, 2)*xyunit_convert, ...
                    patch_data{patch}(:, 3)*zunit_convert);
    end
    hold on;
end



all_points = cell2mat(cellfun(@(c) c, patch_data', 'UniformOutput', false)');

xRange = xyunit_convert*[min(all_points(:, 1)), max(all_points(:, 1))];
yRange = xyunit_convert*[min(all_points(:, 2)), max(all_points(:, 2))];
zRange = zunit_convert*[min(all_points(:, 3)), max(all_points(:, 3))];

xlim(xRange);
ylim(yRange);
zlim(zRange);

pbaspect([range(xRange), range(yRange), 0.8*range(yRange)]);

view([-72, 30]);

xlabel('$x$ [mm]');
ylabel('$y$ [mm]');
zlabel('$z$ [$\mu$m]');

set(gca, 'FontSize', 12);
ax = gca;


set(gcf, 'Renderer', 'painters');
drawnow;


% assert(~red_plot); print('./OUTPUTS/paper_rawscan.eps', '-depsc', '-r400');

