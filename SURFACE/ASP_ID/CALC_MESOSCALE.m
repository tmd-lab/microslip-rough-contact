function [mesoscale_fit, xyh_points] = CALC_MESOSCALE(griddat, settings)
% Function calculates the mesoscale fit of the surface data in griddat
%
% Inputs:
%   griddat  - data on a regular grid from LEVEL_GRID_TRIM
%   settings - has fields:
%       .flip - boolean, true means that the surface will be flipped. If
%               the joints are scanned in the same orientation, then one
%               should be flipped and the other not. 
%       .drawPlots - flag for drawing mesoscale plots
%       .rotateMat - rotation matrix for flipping data to align
%       .translateXYZ - Translation for after rotation to get the data
%                       aligned with the mesh after flipping
%       .meshAlign - alignment added to the coordinates after rotation to
%                           align meaured data with mesh coordinates
%       .xrange    - range of x points for the Mesh
%       .yrange    - range of y points for the Mesh
%       .xyres     - resolution to use for the interpolation grid
%       .mesoFilterSigma - Number of data points for the standard deviation
%                               of the filter
%       .filterSupportSigma - number of standard deviations to be used in
%                              the gaussian filter support (uses
%                              filterSupportSigma*ceil(2*sigma) + 1 points)
% 
% Outputs:
%   mesoscale_fit - mesoscale surface fit function for the current griddat
%   xyh_points    - xyz point list with the z height being the residual
%                   height from the mesoscale surface

    %% Settings
    
    drawPlots = settings.drawPlots;

    rotateMat = settings.rotateMat;
    translateXYZ = settings.translateXYZ;
    meshAlign = settings.meshAlign;
    
    xrange = settings.xrange;
    yrange = settings.yrange;
    xyres  = settings.xyres;

    %% Surface Grid of all points: 
    
    halfXNum = ceil( (mean(xrange) - xrange(1)) / xyres );
    halfYNum = ceil( (mean(yrange) - yrange(1)) / xyres );
    
    xvals = (-halfXNum:halfXNum)*xyres + mean(xrange);
    yvals = (-halfYNum:halfYNum)*xyres + mean(yrange);
    
    [X, Y] = meshgrid(xvals, yvals);
    
    %% Combine Patches into lists of points
    % eliminate the overlapping points
    % Points on patch edges may be repeated, but are averaged by scattered
    % interpolation.
    
    xyz_orig = zeros(3, 0);

    for patch = 1:length(griddat.xyzPatches)

        % Eliminate Asperities that are not in griddat.vertsPatches{patch}
        in = inpolygon(griddat.xyzPatches{patch}(:, 1), griddat.xyzPatches{patch}(:, 2), ...
                        griddat.vertsPatches{patch}(:, 1), griddat.vertsPatches{patch}(:, 2));

        xyz_orig = [xyz_orig; griddat.xyzPatches{patch}(in, :)];

    end
    
    %% Calculation of Mesoscale

    % Eliminate NaNs in orig data
    mask = ~isnan(xyz_orig(:, 3));
    xyz_orig = xyz_orig(mask, :);

    % Flip Bottom / Second Surface
    if(settings.flip)
        xyz_orig = (rotateMat*(xyz_orig'))' + ones(size(xyz_orig(:, 1)))*translateXYZ;
    end

    xyz_orig = xyz_orig + ones(size(xyz_orig(:, 1)))*meshAlign;

    % Baseline Original Surface
    F = scatteredInterpolant(xyz_orig(:,1), xyz_orig(:,2), xyz_orig(:,3), 'natural', 'nearest');

    heights = F(X, Y);

    % Some background: https://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
    filtered_heights = imgaussfilt(heights, settings.mesoFilterSigma, ...
                                    'FilterSize', settings.filterSupportSigma*ceil(2*settings.mesoFilterSigma)+1);

    mesoscale_fit = @(Qps) interp2(X, Y, filtered_heights, Qps(:, 1), Qps(:, 2), 'linear');

    %% Residual Heights to Fit Asperities To
    
    xyh_points = xyz_orig;
    xyh_points(:, 3) = xyh_points(:, 3) - mesoscale_fit(xyh_points(:, 1:2));

    %% Plots for Verifying Results
    
    if(drawPlots)

        figure;
        subplot(1,2,1);
        scatter3(X(:), Y(:), filtered_heights(:), 9, filtered_heights(:), '.'); %9 is the point size
        hold on;
        title(sprintf('Fitted Surface'));

        subplot(1,2,2);
        scatter3(xyh_points(:, 1), xyh_points(:, 2), xyh_points(:, 3), 9, xyh_points(:, 3), '.'); %9 is the point size
        hold on;
        title(sprintf('Residual Heights: Surface'));

    end
    
end