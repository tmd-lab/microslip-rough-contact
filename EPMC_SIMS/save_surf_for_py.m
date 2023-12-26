%% Filenames for Input and Output

input_name = '../SURFACE/OUT/combined_14sep21_R1.mat';
output_name = 'combined_14sep21_R1_4py.mat';
output_dir = '../SURFACE/OUT/for_py';

%% Pull out the parameters to be saved

surf = load(input_name); 

combined = surf.combinedSurf;

Re = sqrt(combined.Rprime*combined.Rpprime);
area_density = combined.area_density;
z_max = combined.z_max;
normzinterp = combined.normzinterp;
pzinterp = combined.pzinterp;


%% Resave mat File

mkdir(output_dir)

output_loc = fullfile(output_dir, output_name);

save(output_loc, 'Re', 'area_density', 'z_max', 'normzinterp', 'pzinterp', '-v7')


