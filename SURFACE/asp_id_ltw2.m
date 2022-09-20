% Script for processing long term wear study surface data
% This version was used after TMD 2021. 

clear;
clc;

addpath('./ASP_ID/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/')

tic

%% Load Data

top = load('SCANS/LongTermWear/R05A_Before_R1.mat');
bot = load('SCANS/LongTermWear/R05B_Before_R1.mat');
scanNum = 1;


% top = load('SCANS/LongTermWear/R05A_Before_R2.mat');
% bot = load('SCANS/LongTermWear/R05B_Before_R2.mat');
% scanNum = 2;


%% Process Settings

%%%%%%%%%%%%% Initial Grid Step
xhalf = 0.060706; 
yhalf = 0.0127;
trimsize = 0.5e-3; % How much to increase holes from tolerance / trim edges
gridSettings.xrange = [-xhalf + trimsize, xhalf - trimsize];
gridSettings.yrange = [-yhalf + trimsize, yhalf - trimsize];


gridSettings.xyres = 23.6e-6; % Resolution in xy directions 
gridSettings.overlap = 10*gridSettings.xyres; % number of overlapped points. 

gridSettings.xyholes = [-1.188*.0254, 0;
           0,            0;
            1.188*.0254, 0]; %+/-0.005

gridSettings.hole_rad = 5e-3*ones(3, 1); %~= (0.332/2+0.0075)*0.0254*ones(3, 1) + trimsize; %+0.075 is for the tolerance of hole + location. Third addition is to try to get rid of drop off near holes for mesoscale
gridSettings.levelData = true; % Level data before putting on grid?

%%%%%%%%%%%%% Mesoscale Settings
mesoSettings.rotateMat = eye(3); 
ind = 1; %[~, ind] = max(range(xyzpoints_top));
mesoSettings.rotateMat(ind, ind) = -1;
mesoSettings.translateXYZ = [0         0         0]; %range(xyzpoints_top).*((1:3) == ind)*combineSettings.toMConvert;
mesoSettings.meshAlign = [0 0  0]; %-range(xyzpoints_top)*combineSettings.toMConvert/2.*[1,1,0];
mesoSettings.xrange = [-xhalf, xhalf]; % Full area so can interpolate onto mesh - extrapolation in combination is nearest.
mesoSettings.yrange = [-yhalf, yhalf];
mesoSettings.xyres = 23.6e-6; % Resolution in xy directions 
mesoSettings.mesoFilterSigma = round(5.08e-3/gridSettings.xyres); %Number of data points for the standard deviation of the filter
mesoSettings.filterSupportSigma = 3; % number of standard deviations in the support of the filter. 2 is MATLAB default, 3 looks better
mesoSettings.drawPlots = true;


%%%%%%%%%%%%% Asperity ID Step

settings.cleanLength = 2; % Default should be 2.
settings.partialArea = true;

%%%%%%%%%%%%% Interface Combination Step

combineSettings.IQR_Factor = 3;
combineSettings.NumInterpPoints = 1001;
combineSettings.drawPlots = true;
combineSettings.maxGaps = 1e5; % sort and then reduce so MATLAB doesn't die. - have about 400e6 gaps possible
combineSettings.maxHistogram = 1e6; % Maximum points to put in the histogram so it doesn't die


%% Clean input data:

scandat = {top, bot};

griddat = cell(size(scandat));

for topbot = 1:length(scandat)
    
    [xyzPatches, vertsPatches, xyzPatches_ungridded, gridX, gridY] ...
        = LEVEL_GRID_TRIM(scandat{topbot}.RawDatabyElem, ...
        gridSettings.xrange, gridSettings.yrange, gridSettings.overlap, ...
        gridSettings.xyres, gridSettings.xyholes, gridSettings.hole_rad, ...
        gridSettings.levelData);
    
    tmp.xyzPatches = xyzPatches;
    tmp.vertsPatches = vertsPatches;
    tmp.X = gridX;
    tmp.Y = gridY;
    
    griddat{topbot} = tmp;

end

clear tmp xyzPatches vertsPatches xyzPatches_ungridded

%% Mesoscale Fit

meso_surf = cell(2, 1);
xyh_points = cell(2, 1);

for topbot = 1:length(scandat)
    
    mesoSettings.flip = (topbot == 2);
    
    [meso_surf{topbot}, xyh_points{topbot}] = CALC_MESOSCALE(griddat{topbot}, mesoSettings);
end

% New Vertices are the grid range vertices
grid_verts = [gridSettings.xrange(2), gridSettings.yrange(2);
                gridSettings.xrange(2), -gridSettings.yrange(2);
                -gridSettings.xrange(2), -gridSettings.yrange(2);
                -gridSettings.xrange(2), gridSettings.yrange(2)];

%% Process asperities

surfdat = cell(size(griddat));
Xmat = cell(size(griddat));
Ymat = cell(size(griddat));
Zmat = cell(size(griddat));
peakMat = cell(size(griddat));
Holes = cell(size(griddat));

for topbot = 1:length(griddat)
            
    % Find Asperities
    [tmp, Xmat{topbot}, Ymat{topbot}, Zmat{topbot}, ...
        peakMat{topbot}, Holes{topbot}] = PROFILE_SURF(xyh_points{topbot}, settings, grid_verts);

    % Eliminate Asperities that are not in the patch area
    in = inpolygon(tmp.xyzloc(:, 1), tmp.xyzloc(:, 2), ...
                    grid_verts(:, 1), grid_verts(:, 2));

    tmp.xyzloc   = tmp.xyzloc(in, :);
    tmp.valid    = tmp.valid(:, in');
    tmp.radii    = tmp.radii(:, in');
    tmp.alpha    = tmp.alpha(:, in');
    tmp.area     = tmp.area;

    surfdat{topbot} = tmp;
    
end

% Save the git hash with the data. 
[~,git_hash_string] = system('git rev-parse HEAD');

mkdir('OUT'); % Issues a warning if exists, prevents fatal save error.
save('./OUT/TMP.mat'); % Always have a save just in case MATLAB crashes


% %%%%%%%%%%%% Save Intermediate Data
date = '2dec21';
% save(sprintf('./OUT/LongTermWear/asperities_%s_R%u_erode%u.mat', date, scanNum, settings.cleanLength), ...
%         'surfdat', 'meso_surf', ...
%         'combineSettings', 'gridSettings', 'mesoSettings', 'settings', ...
%         'git_hash_string', ...
%         'Xmat', 'Ymat', 'Zmat', 'peakMat', 'Holes');

surfs_time = toc

%% Combine two surfaces with meso-scale analysis

% %%%%%%%%%%%% Load Interface Information
% scanNum = 1;
% date = '14sep21';
% load(sprintf('./OUT/LongTermWear/asperities_%s_R%u.mat', date, scanNum));

tic;

% Apply combination
combinedSurf = COMBINE_SURFS2(surfdat{1}, surfdat{2}, combineSettings);

combine_time = toc


% Save the git hash with the data. 
[~,git_hash_string] = system('git rev-parse HEAD');

% %%%%%%%%%%%% Save Results
% save(sprintf('./OUT/LongTermWear/combined_%s_R%u.mat', date, scanNum), ...
%                 'combinedSurf', 'meso_surf', ...
%                 'combineSettings', 'gridSettings', 'mesoSettings', 'settings', ...
%                 'git_hash_string');


% save(sprintf('./OUT/LongTermWear/combined_%s_R%u_erode%u.mat', date, scanNum, settings.cleanLength), ...
%                 'combinedSurf', 'meso_surf', ...
%                 'combineSettings', 'gridSettings', 'mesoSettings', 'settings', ...
%                 'git_hash_string');
 

