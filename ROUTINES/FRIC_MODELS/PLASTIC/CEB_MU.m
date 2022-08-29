function [mu, dmudun, dmudunmax] = CEB_MU(pars, un, unmax)
% Function for determining the slip limit of the CEB model from: 
% This comes from:
%       Eriten et al., 2010, Physics-based modeling for partial slip
%       behavior of spherical contacts
%
% Some assumptions on the application here differ in that the yield
% displacement is taken to be deltam since that is when the sphere
% would yield again according to the hardening model 
%
% Inputs:
%   pars - structure of material and asperity parameters
%   un   - positive being in contact. Difference between uxyn(3)-deltabar
%   unmax - previous maximum displacement minus deltabar :
%               unmax = deltam - deltabar
%
% Outputs:
%   mu      - slip limit from Eriten paper
%   dmudun  - derivative w.r.t. normal displacement
%   dmuddeltam - derivative w.r.t. normal history
%
% NOTES:
%   1. WARNING: At 0 normal interference, the returned friction coefficient
%   should be infinite. However it is set to zero to prevent NaN's since
%   Hertzian contact goes to zero faster than CEB goes to Inf. However, if
%   not using with a Hertzian contact this could be wrong. Though, this
%   should be the most physical solution.

    K = 0.454 + 0.41*pars.nu; %Top of page 3

    assert(pars.nu == 0.3, 'Value for zeta is only given for nu=0.3. Try a different slip type or update the implementation for your value.');
    zeta = 0.48;

    % Normalized normal displacement
    if(unmax < pars.delta_y)
        % Normalize to initial yielding
        % This can be used for plotting and testing, but should never be
        % reached by dynamic analysis. Dynamic analysis uses the strain
        % hardening model that defines yielding as starting at 1.9*delta_y,
        % so for consistency, 1.9*delta_y is the minimum input for unmax
        % from the dynamic friction functions.
        wstar           = un/pars.delta_y;
        dwstar_dun      = 1/pars.delta_y;
        dwstar_ddeltam  = 0;
    else
        % Normalize to reyielding if strain hardened
        wstar           = un/unmax;
        dwstar_dun      = 1/unmax;
        dwstar_ddeltam  = -un/unmax^2;
    end

    % CEB "Constants" - following eqn 16
    c1 = -1 + 1.5*zeta*atan(1/zeta) - zeta^2/(2 * (1 + zeta^2));
%     c2 = (1 + pars.nu)*(zeta*atan(1/zeta) - 1) + 3/(2*(1 + zeta^2)); % Not actually used
    c3 = 9*pi^2/16*(2-pars.nu/2 + 7/8*pars.nu^2);
    c4 = 9*pi/4*(1 - 2*pars.nu)*(1-pars.nu/2);

    % NOT A CONSTANT - derivative below (if needed)
    c5 = 1.5*(1 - 2*pars.nu)^2 - 0.56/K^2/wstar;

    % Two friction coefficients that the min is taken of
    mu1 = (0.2045/(K*abs(c1)))*sqrt(1/wstar - 1);
    mu2 = ( -c4 + sqrt(c4^2 - 4*c3*c5) )/(2*c3);

    if(mu1 <= mu2)
        % First CEB option is the minimum
        
        mu = mu1;
        dmudun     = (0.2045/(K*abs(c1)))/(2*sqrt(1/wstar - 1)) * ( -1/wstar^2 ) * dwstar_dun;
        dmudunmax = (0.2045/(K*abs(c1)))/(2*sqrt(1/wstar - 1)) * ( -1/wstar^2 ) * dwstar_ddeltam;

    else % mu2 < mu1 is known
        % Second CEB option is the minimum

        dc5dun      = (0.56/K^2/wstar^2)*dwstar_dun;
        dc5ddeltam  = (0.56/K^2/wstar^2)*dwstar_ddeltam;

        mu = mu2;
        dmudun     = (1/2/c3)*( -2*c3 / sqrt(c4^2 - 4*c3*c5) ) * dc5dun;
        dmudunmax = (1/2/c3)*( -2*c3 / sqrt(c4^2 - 4*c3*c5) ) * dc5ddeltam;

    end
    
    if(isinf(mu) || wstar >= 1) % need to include 1 otherwise get a Inf derivative
        
        % return zero since Hertzian contact goes to zero faster than CEB
        % goes to inf
        mu = 0;
        dmudun = 0;
        dmudunmax = 0;
    end

end