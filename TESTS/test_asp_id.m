% Script to Test functions used in asperity identification with the
% Watershed algorithm:


clear;

addpath('../SURFACE/ASP_ID/')

%% Generate a Surface to ID

Npoints = 30;
freqx = 3.75;
freqy = 2;

[X, Y] = meshgrid(linspace(0, 1, Npoints), linspace(0, 1, Npoints));

Z = sin(2*pi*X*freqx).*sin(2*pi*Y*freqy);

Z = (Z - min(Z(:)))/range(Z(:));

%% Hole Locations

Holes = false(size(Z));

%% Add some noise:

mean = 0;
var = 0.05^2;

rng('default');
Z = imnoise(Z, 'gaussian', mean, var);

%% Plot Original Surface:

figure;
surf(X, Y, Z, Z, 'EdgeColor', 'none');
hold on;
title('Height Colored');


settings.cleanLength = 3;

segmented = FIND_PEAKS_WATERSHED(Z, Holes, settings);

figure;
surf(X, Y, Z, segmented, 'EdgeColor', 'none');
hold on;
title('Segment Colored');

