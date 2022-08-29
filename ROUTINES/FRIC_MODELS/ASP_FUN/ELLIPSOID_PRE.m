function [fxyn, dfxynduxyn, rq0, tx0, ty0, varargout] = ELLIPSOID_PRE(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq)
% Function determines forces for the contact state of an asperity. Only
% considers normal contact for prestressed state. = Elliptical contact area
% 
% Inputs:
%   pars - structure of material and asperity parameters - intialize
%           elliptical contact parameters with: pars = INITIALIZE_ECCENTRIC(pars);
%   uxyn - displacements, un is the gap distance with positive being in
%           contact
%   uxyn0 - previous displacements, un0 is again gap distance
%   rq0   - quadrature radii at last instant
%   tx0   - tractions in x direction at rq0 at last instant
%   ty0   - tractions in y direction at rq0 at last instant
%   rq    - quadrature radii to be used normalized to a, so range [0, 1]
%   wq    - weights for the quadrature radii. 
%
% Outputs:
%   fxyn  - forces in xyn directions - zeros in xy directions
%   dfxynduxyn - Jacobian stiffness matrix of fyxn 
%   rq0 - quadrature radii at the end - equals rq
%   tx0 - tractions at each radii after instant - zeros
%   ty0 - tractions at each radii after this instant. - zeros
%   vargargout{1} - contact dimension b
%   vargargout{2} - derivative db/dun
%   vargargout{3} - contant dimension a
%   vargargout{4} - derivative da/dun
%   vargargout{5} - Maximum pressure p0
%   vargargout{6} - derivative dp0/dun
    
    %% Normal Contact Model 
    
    fxyn(1, 3) = 4*pars.Estar*sqrt(pars.Re)/3/(pars.F2)^(3/2).*uxyn(3).^(3/2).*(uxyn(3) > 0);
    
    dfxynduxyn(3,3) = 2*pars.Estar*sqrt(pars.Re)/(pars.F2)^(3/2).*uxyn(3)^(1/2).*(uxyn(3) > 0);

    
    %% Final Details for outputting previous information - ignored since doing zero tangential contact
    rq0 = rq;
    tx0 = zeros(size(rq));
    ty0 = zeros(size(rq));
    
    %% Contact parameters, not always needed
    if(nargout > 5)
        

%         tmp = (3*fxyn(3)*pars.Re/4/pars.Estar)^(1/3);
%         c = tmp*pars.F1;
%         dcdun = 1/3 /tmp^2*pars.F1*dfxynduxyn(3,3);
        
        %Intermediate parameter defined in Johnson
        c = (3*fxyn(3)*pars.Re/4/pars.Estar)^(1/3)*pars.F1;
        
        dcdun = 1/3*(3*pars.Re/4/pars.Estar)^(1/3)*fxyn(3)^(-2/3)*pars.F1*dfxynduxyn(3,3);
        
        if(c == 0)
            % Remove NaN, one sided derivative is correct
            dcdun = 0;
        end
        
        %Short axis length
        varargout{1} = c/sqrt(pars.aperb); %b
        varargout{2} = dcdun/sqrt(pars.aperb); %db/dun
        
        %Long axis length
        varargout{3} = varargout{1}*pars.aperb; %a
        varargout{4} = varargout{2}*pars.aperb; %da/dun
    end
    
    if(nargout > 9)
        %Maximum pressure
        varargout{5} = 3*fxyn(3)/2/pi/varargout{1}/varargout{3}; %p0
        
        varargout{6} = 3*dfxynduxyn(3,3)/2/pi/varargout{1}/varargout{3} ...
                            -3*fxyn(3)/pi/pars.aperb/varargout{1}^3*varargout{2}; %dp0/dun
        
    end
    
end