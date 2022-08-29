function [fxyn, dfxynduxyn, rq0, tx0, ty0] = ASP_ROTATION(ORTHO_ASP_FUN, pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq)
% Function determines forces for the contact state of an asperity as
% rotated with the short principal axis pars.alpha radians counter
% clockwise of the +x axis.
% 
% Inputs:
%   pars - structure of material and asperity parameters - intialize
%           elliptical contact parameters with: pars = INITIALIZE_ECCENTRIC(pars);
%   uxyn - displacements, un is the gap distance with positive being in
%           contact
%   uxyn0 - previous displacements, un0 is again gap distance
%   rq0   - quadrature radii at last instant - scaled to a circle of radius
%               a.
%   tx0   - tractions in x direction at rq0 at last instant
%   ty0   - tractions in y direction at rq0 at last instant
%   rq    - quadrature radii to be used normalized to a, so range [0, 1]
%   wq    - weights for the quadrature radii. 
%
% Outputs:
%   fxyn  - forces in xyn directions 
%   dfxynduxyn - Jacobian stiffness matrix of fyxn 
%   rq0 - quadrature radii at the end - equals rq*a
%   tx0 - tractions at each radii after instant 
%   ty0 - tractions at each radii after this instant. 
%
% Notes:
%   1. History variables that are output here are stored in the local 
%   coordinate system after rotation
    
    % Transform:
    %   u_{b,a,n} = Q u_{xyn} with column vectors
    Q = [cos(pars.alpha),  sin(pars.alpha), 0;
         -sin(pars.alpha), cos(pars.alpha), 0;
         0,                0,               1];
     
     
%      uban = (Q * uxyn')'; % (Derivation)
     uban = uxyn*Q';
     uban0 = uxyn0*Q';
     
     [fban, dfbanduban, rq0, tx0, ty0] = ORTHO_ASP_FUN(pars, uban, uban0, rq0, tx0, ty0, rq, wq);
     
%      Q' = inv(Q); % (Derivation)
%      fxyn = (Q'*fban')'; % (Derivation)
     fxyn = fban*Q;
     
     dfxynduxyn = Q' * dfbanduban * Q; 
    
end


