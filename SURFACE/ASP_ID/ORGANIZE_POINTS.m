function [Xmat, Ymat, Zmat, Holes] = ORGANIZE_POINTS(xyzpoints)
% Function for organizing all of the data into a regular grid and applying
% smoothing.
%
% Inputs:
%   xyzpoints - (n x 3) matrix with columns xyz and row are points
%   smoothLengthX - size of smoothing filter for x direction, same units as 
%                   xyzpoints.
%   smoothLengthY - size of smoothing filter for y direction, same units as 
%                   xyzpoints.
%   
% Outputs:
%   Xmat - matrix of x values
%   Ymat - matrix of y values
%   Zmat - matrix of original Z values
%   Holes - matrix with trues corresponding to locations of holes. 
%
% NOTES:
%   1. This assumes that the x and y data is from a regular grid. If this
%       is not the case, the function may not work as appropriate. 

    
    %% Organize Data into grid.
    
    Zmat=Inf*ones(length(unique(xyzpoints(:,2))),length(unique(xyzpoints(:,1))));

    %changem uses mapping toolbox    
    Xinds=changem(xyzpoints(:,1),transpose(1:length(unique(xyzpoints(:,1)))),unique(xyzpoints(:,1)));
    Yinds=changem(xyzpoints(:,2),transpose(1:length(unique(xyzpoints(:,2)))),unique(xyzpoints(:,2)));
    
    ind=sub2ind(size(Zmat),Yinds,Xinds);
    Zmat(ind)=xyzpoints(:,3);

    [Xmat,Ymat]=meshgrid(unique(xyzpoints(:,1)),unique(xyzpoints(:,2)));
    
    Holes = ~isfinite(Zmat);
    
end