function [R, dRdaperb, K, E, e, dKdaperb, dEdaperb, dedaperb] = CONTACT_ECCENTRICITY_RES(aperb, Rprime, Rpprime)
% Function for calculating the residual to define the eccentricity of the
% the contact ellipse
% Inputs:
%   aperb - quantity (a/b) which describes the ratio of major and minor
%           axes of the elliptical contact area. This ratio is greater than
%           1 (will probably break if greater than 1). 
%   Rprime  - Radius R' as defined in Johnson. 
%   Rpprime - Radius R'' as defined in Johnson. R'>R''
% Outputs:
%   R  - Residual
%   dRdaperb - derivative of residual
%
% A good initial guess is R'/R''

    if(aperb >= 1)

        e = (1 - aperb^(-2))^(1/2);

        dedaperb = 1/(aperb^3)/sqrt(1 - aperb^(-2));

        %Calculate the elliptical values at e
        [K,E] = ellipke(e^2);

        % Derivatives of Elliptical contact, from:
        % https://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html
        % and
        % https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html
        dKdaperb = (E/e/(1- e^2) - K/e)*dedaperb;
        dEdaperb = ((E - K)/e)*dedaperb;

        %%%% Alternative approach
        %Residual
        R = (aperb^2*E - K)*Rpprime - (K - E)*Rprime;

        %Derivative
        dRdaperb = (2*aperb*E + aperb^2*dEdaperb - dKdaperb)*Rpprime - (dKdaperb - dEdaperb)*Rprime;
    
    else
        
        R = NaN;
        dRdaperb = NaN;
        
    end

end