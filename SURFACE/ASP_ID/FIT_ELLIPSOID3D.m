function [radii, xyzloc, alpha, valid, center, axesq, csys, Serr, Bmat, Lvec, Csca, mxyz, Tm] = FIT_ELLIPSOID3D(xyz)
%FIT_ELLIPSOID3D fits 3D data to 3D ellipsoid (rotated) using TLS where only
%rotations about the z axis are permitted. The parametric form that is 
%fitted is:
%           (x-mxyz)^T*Bmat*(x-mxyz) + Lv^T*(x-mxyz) + Csca = 0
%     where x = [x;y;z] 3D cartesian coordinates vector that has been
%     scaled by Tm
%
%   USAGE:
%       [radii, xyzloc, alpha, valid, center, axesq, csys, Serr, Bmat, Lvec, Csca, mxyz, Tm] = FIT_ELLIPSOID3D(xyz)
%   INPUTS:
%       xyz     : (Npt, 3)
%   OUTPUTS:
%       radii   : (1, 2) Principal radii of curvature at peak, larger one 
%                   is the second entry
%       xyzloc  : (1, 3) Location of the peak of the asperity
%       alpha   : Rotation angle around z axis from the x+ axis to the
%                   short semi-axis range of (-pi/2, pi/2]
%       valid   : Indicates whether the outputs are valid (true) or not
%                   Checks include: Is the peak within the range of the xy
%                   coordinates? Is the fit an ellipsoid? Were sufficient
%                   points (>= 9) provided to have a unique solution?
%       center  : (3, 1) center
%       axesq    : (3, 1) semi-minor axes lengths squared 
%       csys    : (3, 3) coordinate transform matrix (affine)
%       Serr    : Scaled Least-squares error of fit 
%       Bmat    : (3, 3) B matrix 
%       Lvec    : (3, 1) L vector 
%       Csca    : C scaling 
%       mxyz    : mean location of points
%       Tm      : Scaling vector
%
% Notes:
%   1. Implementation provided by Nidish Balaji from his PhD Dissertation
%   with modifications to outputs
%   2. The Z direction is fixed vertical in the fit.
%   3. A minimum of 9 points is required so that it is 1 more than the
%   number of least squares coefficients

    if size(xyz,1)==3 && size(xyz,2)~=3
        xyz = xyz';
    end
    
    % Remove any potential infinite points e.g., from invalid boundary
    % points
    xyz = xyz(isfinite(sum(xyz, 2)), :);
    
    valid = size(xyz, 1) >= 9;
    
    if(valid)
        % Scaling
%     scalexyz = [1, 1, 1];
        scalexyz = range(xyz);
        
        valid = valid & scalexyz(1) > 0 & scalexyz(2) > 0;
    end
    
    % Recheck Valididity
    if(valid)  
        
        Xminmax = [min(xyz(:, 1)), max(xyz(:, 1))];
        Yminmax = [min(xyz(:, 2)), max(xyz(:, 2))];

        % Center + Scale Data
        mxyz = mean(xyz);
        
        Tm = diag(scalexyz);
        Tmi = diag(1./scalexyz);
        xyz = (xyz-mxyz)*Tmi;

        try
            [~, S, V] = svd([xyz.^2 xyz(:,1).*xyz(:,2) xyz ones(size(xyz,1),1)]);
        catch 
            keyboard
        end

        S = diag(S);
        Serr = S(end)/S(1); % Sum of Squares of errors is S(end)^2
        acofs = V(:, end); % The smallest output is the coefficients aligned with the smallest singular value

        % Assemble Matrices (rescaled)
        Bmat = Tmi'*[acofs(1) acofs(4)/2 0;
            acofs(4)/2 acofs(2) 0;
            0 0 acofs(3)]*Tmi;

        Lvec = Tmi*acofs(5:7);
        Csca = acofs(8);

        center = -Bmat\Lvec/2;

        [bU, bS, bV] = svd(Bmat);
        sc = center'*Bmat*center-Csca;
        eV = bU;
        eD = diag(bU'*Bmat*bU)/sc;

        si = [1 2 3]; 
        [~, si(1)] = max(abs(eV(1,:)));

        sim = setdiff([1 2 3], si(1));
        [~, si(2)] = max(abs(eV(2,sim))); si(2) = sim(si(2));

        si(3) = setdiff([1 2 3], si(1:2));

        eD = eD(si); eV = eV(:, si);


        eV = eV*diag(sign(diag(eV)));
        axesq = 1./eD;  % axes^2 of ellipse (or hyperboloid/paraboloid if negative)
        axesq(~isfinite(axesq)) = 0.0;
        csys = eV;

        if sum(axesq>0)~=3
            valid = false;
    %         error('Not an ellipse')
        end

        % Matrices for return
        Bmat = Bmat/sc;
        Lvec = Lvec/sc;
        Csca = Csca/sc;
        
        % Add back the mean
        center = center + mxyz';
        mxyz = mxyz';
        
        %% Additional outputs

        % Principal Radii
        if( axesq(1) <= axesq(2) )
            shortvec = csys(:, 1);
            radii = [axesq(1), axesq(2)] / sqrt(axesq(3));
        else
            shortvec = csys(:, 2);
            radii = [axesq(2), axesq(1)] / sqrt(axesq(3));
        end

        % Location
        xyzloc = center' + [0, 0, sqrt(axesq(3))];

        % Orientation
        alpha = atan(shortvec(2)/shortvec(1));
        alpha = alpha + pi*(alpha == -pi/2);

        % Check that the center lies within the bounding box
        valid = valid ...
                    && (xyzloc(1) >= Xminmax(1)) ...
                    && (xyzloc(1) <= Xminmax(2)) ...
                    && (xyzloc(2) >= Yminmax(1)) ...
                    && (xyzloc(2) <= Yminmax(2));

        % Verify Local Maximum
        valid = valid & (center(3) < mxyz(3));
        
        % Verify that there are at least 3 points in both directions
        valid = valid & length(unique(xyz(:, 1))) >=3 ...
            & length(unique(xyz(:, 2))) >=3;
        
    else
        
        radii = zeros(1, 2);
        xyzloc = zeros(1, 3);
        alpha = 0;
        center = zeros(3, 1);
        axesq = zeros(3, 1);
        csys = zeros(3, 3);
        Serr = 0;
        Bmat = zeros(3, 3);
        Lvec = zeros(3, 1);
        Csca = 0;
        mxyz = zeros(3, 1);
        
    end
end