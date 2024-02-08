%% Filenames for Input and Output

input_name = '../SURFACE/OUT/combined_14sep21_R1.mat';
output_name = 'brb_surface_data.mat';
output_dir = '../SURFACE/OUT/for_py';

%% Pull out the parameters to be saved

surf = load(input_name); 

combined = surf.combinedSurf;

Re = sqrt(combined.Rprime*combined.Rpprime);
area_density = combined.area_density;
z_max = combined.z_max;
normzinterp = combined.normzinterp;
pzinterp = combined.pzinterp;

%% Calculate Mesoscale Gap at Some Points

approx_meso_grid = 0.5e-3; %  0.5 mm

% Expand out a single point in each direction
meso_x_set = linspace(surf.mesoSettings.xrange(1)-approx_meso_grid, ...
                      surf.mesoSettings.xrange(2)+approx_meso_grid, ...
                      ceil(range(surf.mesoSettings.xrange)/approx_meso_grid)+3);

meso_y_set = linspace(surf.mesoSettings.yrange(1)-approx_meso_grid, ...
                      surf.mesoSettings.yrange(2)+approx_meso_grid, ...
                      ceil(range(surf.mesoSettings.yrange)/approx_meso_grid)+3);

[X,Y] = meshgrid(meso_x_set, meso_y_set);

mesoscale_xy = [X(:), Y(:)];

gaps = -surf.meso_surf{1}(mesoscale_xy) - surf.meso_surf{2}(mesoscale_xy);
gaps = gaps - min(gaps); % move to start at 0.0

mesoscale_xygap = [mesoscale_xy, gaps(:)];

%% Resave mat File

mkdir(output_dir)

output_loc = fullfile(output_dir, output_name);

save(output_loc, 'Re', 'area_density', 'z_max', 'normzinterp', 'pzinterp', ...
     'mesoscale_xygap', '-v7')


