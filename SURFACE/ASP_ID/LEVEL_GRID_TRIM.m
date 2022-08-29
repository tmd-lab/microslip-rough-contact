function [xyzPatches, vertsPatches, xyzPatches_ungridded, X, Y] ...
    = LEVEL_GRID_TRIM(RawDatabyElem, xrange, yrange, overlap, xyres, xyholes, hole_rad, levelData)
% Function levels the initial scan data, places each scan patch onto a grid
% (but keeps the patches separate), and then trims the holes and edges to
% minimize distortions where the surface drops off rapidly. 
%
% Inputs: 
%   RawDatabyElem - scan data divided between different scans and elements
%   xrange - [xmin, xmax]
%   yrange - [ymin, ymax]
%   overlap - the amount of overlap to leave between patches
%   xyres - resolution for the gridded data in the x and y directions.
%   xyholes - hole locations of [x, y] with rows = different holes
%   hole_rad - radius of holes
%   levelData - boolean flag, if false the data if not leveled.
%
% Outputs:
%   xyzPatches - cell structure with a single matrix of xyzpoints for each
%                   scan with overlap
%   vertsPatches - vertices of each patch so overlap can be easily
%                   eliminated after asperity ID.
%   xyzPatches_ungridded - Data for the patches without gridding
%   X - X grid coordinates in matrix
%   Y - Y grid coordinates in matrix
%   
%
% Notes:
%   1. Trims the holes so they have hole rad
%   2. Trims outside to match xrange and yrange. Change these values to
%   eliminate noisy points on edges of data




    %% Level the entire input data
    if(levelData)
        
        %%%% All data
        xyz_all =  cell2mat(cellfun(@(c) c, RawDatabyElem(:)', 'UniformOutput', false)');    


        %%%% Find Least Squares Plane
        %plane of coefs(1) + coefs(2)*x + coefs(3)*y = z
        coefs = [ones(length(xyz_all(:, 1)), 1), xyz_all(:, 1), xyz_all(:, 2)] \ xyz_all(:, 3);

        %%%% Apply rotation to level plane

        vecZ = [-coefs(2), -coefs(3), 1];
        vecZ = vecZ / norm(vecZ);

        vecY = cross(vecZ, [1, 0, 0]);

        vecX = cross(vecY, vecZ);

        % Transforms as [x;y;z]_current = QtoCurrent*[x; y; z]_leveled
        QtoCurrent = [vecX'/norm(vecX), vecY'/norm(vecY), vecZ'];

        LevelDatabyPatch = cell(size(RawDatabyElem, 1), 1);

        for patch = 1:length(LevelDatabyPatch)

            xyzPatch =  cell2mat(cellfun(@(c) c, RawDatabyElem(patch, :), 'UniformOutput', false)');

            % Apply: [x; y; z]_leveled = Qinv*[x;y;z]_curr
            %        [x; y; z]_leveled = QtoCurrent' *[x;y;z]_curr (inv = transpose)
            %        [x; y; z]_leveled' = (QtoCurrent' *[x;y;z]_curr)' (desired shape)
            %        [x; y; z]_leveled' = ([x;y;z]_curr')*QtoCurrent (transpose props)

            LevelDatabyPatch{patch} = xyzPatch*QtoCurrent;

        end
        
    else
        warning('Not leveling the data before placing on grid.');
    end
    
    %% Grid Data + Trim out Holes
    
    [xyzPatches, vertsPatches, xyzPatches_ungridded, X, Y] ...
        = ELIMINATE_OVERLAP(LevelDatabyPatch, xrange, yrange, overlap, xyres, xyholes, hole_rad);
    

end