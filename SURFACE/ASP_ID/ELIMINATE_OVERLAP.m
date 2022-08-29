function [xyzPatches, vertsPatches, xyzPatches_ungridded, X, Y] ...
    = ELIMINATE_OVERLAP(RawDatabyElem, xrange, yrange, overlap, xyres, xyholes, hole_rad)
% Function to eliminate overlapped data between patches
% 
% Inputs:
%   RawDatabyElem - scan data divided between different scans and elements
%   xrange - [xmin, xmax]
%   yrange - [ymin, ymax]
%   overlap - the amount of overlap to leave between patches
%   xyres - resolution for the gridded data in the x and y directions.
%   xyholes - hole locations of [x, y] with rows = different holes
%   hole_rad - radius of holes
%
% Outputs:
%   xyzPatches - cell structure with a single matrix of xyzpoints for each
%                   scan with overlap
%
% NOTES:
%   1. Assumes two patches in the y-direction and that the mean(yrange) is
%   a good place to split all of them along the y-direction. 

    xyzPatches = cell(size(RawDatabyElem, 1), 1);
    vertsPatches = cell(size(RawDatabyElem, 1), 1);
    
    centerxy = [(1:length(xyzPatches))', zeros(length(xyzPatches), 2)]; % Patch Number, center x, center y

    %% Create Patch Info
    
    for patch = 1:length(xyzPatches)

       xyzPatches{patch} =  cell2mat(cellfun(@(c) c, RawDatabyElem(patch, :), 'UniformOutput', false)');    
       
       centerxy(patch, 2:end) = mean(xyzPatches{patch}(:, 1:2));

    end

    %% Search for Boundaries
    
    yrows = centerxy(:, 3) > mean(yrange);

    for patch = 1:length(xyzPatches)
        
        currentY = (yrows == yrows(patch));
        
        currCenters = centerxy(currentY, :);
        
        patchyrange = [yrange(1), mean(yrange)] .*(~yrows(patch)) ...
                        + [mean(yrange), yrange(2)] .*(yrows(patch));
        
        % Find the patch to larger x
        if(centerxy(patch, 2) == max(currCenters(:, 2)) )
            % Larger neighbor is the global upper limit
            smallerx = interp1(currCenters(:, 2), currCenters(:, 1), centerxy(patch, 2)-eps, 'previous');
            
            patchxrange = [(max(xyzPatches{smallerx}(:, 1)) + min(xyzPatches{patch}(:, 1)))/2, xrange(2)];
        elseif(centerxy(patch, 2) == min(currCenters(:, 2)) )
            % lower neightbor is the global lower limit
            largerx  = interp1(currCenters(:, 2), currCenters(:, 1), centerxy(patch, 2)+eps, 'next');
            
            patchxrange = [xrange(1), (min(xyzPatches{largerx}(:, 1)) + max(xyzPatches{patch}(:, 1)))/2];
        else
            % Have larger and smaller neighbor
            smallerx = interp1(currCenters(:, 2), currCenters(:, 1), centerxy(patch, 2)-eps, 'previous');
            largerx  = interp1(currCenters(:, 2), currCenters(:, 1), centerxy(patch, 2)+eps, 'next');
            
            patchxrange = [(max(xyzPatches{smallerx}(:, 1)) + min(xyzPatches{patch}(:, 1)))/2, ...
                           (min(xyzPatches{largerx}(:, 1)) + max(xyzPatches{patch}(:, 1)))/2];
        end
        
        vertsPatches{patch} = [patchxrange(1), patchyrange(1);
                               patchxrange(1), patchyrange(2);
                               patchxrange(2), patchyrange(2);
                               patchxrange(2), patchyrange(1)];
                             
    end
    
    %% Reduce the amount of points that are output + Grid the data up
    
    xyzPatches_ungridded = cell(size(xyzPatches));
    holerad_measured = zeros(size(xyzPatches_ungridded));
    
    halfXNum = ceil( (mean(xrange) - xrange(1)) / xyres );
    halfYNum = ceil( (mean(yrange) - yrange(1)) / xyres );
    
    xvals = (-halfXNum:halfXNum)*xyres + mean(xrange);
    yvals = (-halfYNum:halfYNum)*xyres + mean(yrange);
    
    [X, Y] = meshgrid(xvals, yvals);
    
    
    for patch = 1:length(xyzPatches)
        
        %%%%%%%%%% Reduce Data Points
        patchExpand = vertsPatches{patch}   + overlap*[-1, -1;
                                                       -1,  1;
                                                        1,  1;
                                                        1, -1];
        
        in = inpolygon(xyzPatches{patch}(:, 1), xyzPatches{patch}(:, 2), patchExpand(:, 1), patchExpand(:, 2));
        
        xyzPatches_ungridded{patch} = xyzPatches{patch}(in, :);
        
        %%%%%%%%%%%%% Put on Grid
        XYGrid = [X(:), Y(:)];
        
        GridIn = inpolygon(XYGrid(:, 1), XYGrid(:, 2), patchExpand(:, 1), patchExpand(:, 2));
        
        xyzPatches{patch} = [XYGrid(GridIn, :), zeros(sum(GridIn), 1)];
        
        F = scatteredInterpolant(xyzPatches_ungridded{patch}(:,1), ...
                                 xyzPatches_ungridded{patch}(:,2), ...
                                 xyzPatches_ungridded{patch}(:,3), ...
                                 'natural', 'none');
                             
        xyzPatches{patch}(:, 3) = F(xyzPatches{patch}(:, 1), xyzPatches{patch}(:, 2));
        
        %%%%%%%%%%%%%% Check the radius of the data that was given.
        
        dx = xyzPatches_ungridded{patch}(:, 1) - xyholes(:, 1)';
        dy = xyzPatches_ungridded{patch}(:, 2) - xyholes(:, 2)';
        
        holerad_measured(patch) = sqrt(min( dx.^2 + dy.^2, [], 'all'));
        
        %%%%%%%%%%%%%% Eliminate Holes from Gridded Data

        dx = xyzPatches{patch}(:, 1) - xyholes(:, 1)';
        dy = xyzPatches{patch}(:, 2) - xyholes(:, 2)';
        
        [distancesq, hole_ind] = min(dx.^2 + dy.^2, [], 2);
        
        xyzPatches{patch} = xyzPatches{patch}(sqrt(distancesq) > hole_rad(hole_ind), :);
        
    end
    
    fprintf('The minimum measured hole radius is %.3s. The average of the minimum 8 is %.3s. The min input radius is %.3s.\n', ...
        min(holerad_measured), mean(mink(holerad_measured, 6)), min(hole_rad));

    
end