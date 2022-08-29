function [curr_surf, Xmat, Ymat, Zmat, peakMat, Holes] = PROFILE_SURF(xyzpoints, settings, patchVerts)
% This function profiles an entire surface given the xyz points and
% settings.
%
% Inputs:
%   xyzpoints - matrix with rows=points and 3 columns for x,y,z. Units
%               don't matter so long as the appropriate setting is used to
%               convert to meters.
%   settings - structure with fields for settings. Fields are as follows.
%               Descriptions can be found in comments for appropriate
%               functions
%       ORGANIZE_POINTS:
%           (None)
%       FIND_PEAKS
%           settings.
%       FIT_ASPERITY
%           settings.
%       settings.partialArea - uses patchVerts to calculated a reduced area
%                               of the searched space
%   patchVerts - vertices of a patch. If settings.partialArea = true, then
%                   only the area of this patch is calculated
%
% Outputs:
%   curr_surf - struct with asperity statistics
%   Other outputs - various location matrices for looking at the process.
%

    % Clean up the data into matrices
    [Xmat, Ymat, Zmat, Holes] = ORGANIZE_POINTS(xyzpoints);

    % Identify asperity peaks
    peakMat = FIND_PEAKS_WATERSHED(Zmat, Holes, settings);

    % Processes Asperities into Radius Information

    num_asp = max(peakMat(:));

    valid = false(1, num_asp);
    radii = zeros(2, num_asp);
    xyzloc = zeros(num_asp, 3);
    alpha = zeros(1, num_asp);
    
    % reduce to only the points of interest -> Less communication and faster
    % May be mostly irrelevant with the watershed divisions now filling up
    % most of the area.
    Xmat_tmp = Xmat(peakMat > 0);
    Ymat_tmp = Ymat(peakMat > 0);
    Zmat_tmp = Zmat(peakMat > 0);
    peakMat_tmp = peakMat(peakMat > 0);

%     warning('Fitting asperities in serial');
    parfor ii = 1:num_asp

        [radii(:, ii), xyzloc(ii, :), alpha(ii), valid(ii), ...
            ~, ~, ~, ~, ~, ~, ~] = FIT_ELLIPSOID3D([Xmat_tmp(peakMat_tmp == ii), ...
                                                    Ymat_tmp(peakMat_tmp == ii), ...
                                                    Zmat_tmp(peakMat_tmp == ii)]);
    end
    
    
    if(settings.partialArea)
        % If the searched area exceeds the full area, it is desired to
        % reduce the area to only be the full area of the patch for later
        % calculations.
        
        in = inpolygon(Xmat, Ymat, patchVerts(:, 1), patchVerts(:, 2));
        
        area = trapz(Xmat(1, :), trapz(Ymat(:, 1), 1.0*(~Holes & in)));

    else
    
        % Calculate area, this is a bit less than the mesh area, but this is
        % more accurate since this is the area that was checked for asperities.
        area = trapz(Xmat(1, :), trapz(Ymat(:, 1), 1.0*(~Holes)));
    
    end
    
    %Create a data structure to output properties
    curr_surf.xyzloc = xyzloc;
    curr_surf.valid = valid;
    curr_surf.radii = radii;
    curr_surf.alpha = alpha;
    curr_surf.area = area;

    % Properties about processing
    curr_surf.settings = settings;

end