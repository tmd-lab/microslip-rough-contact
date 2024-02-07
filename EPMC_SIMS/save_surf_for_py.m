%% Filenames for Input and Output

% input_name = '../SURFACE/OUT/combined_14sep21_R1.mat';
% output_name = 'combined_14sep21_R1_4py.mat';
% output_dir = '../SURFACE/OUT/for_py';

input_name = '../SURFACE/OUT/combined_2dec21_R1_erode2.mat';
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

meso_xyh_top = xyh_meso_points{1};
meso_xyh_bot = xyh_meso_points{2};

%% Resave mat File

mkdir(output_dir)

output_loc = fullfile(output_dir, output_name);

save(output_loc, 'Re', 'area_density', 'z_max', 'normzinterp', 'pzinterp', ...
     'meso_xyh_top', 'meso_xyh_bot', '-v7')


