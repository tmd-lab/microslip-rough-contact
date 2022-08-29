function [pars] = INITIALIZE_ECCENTRIC(pars)
% This function determines the costants related to eccentric elliptoid on
% elliptoid Hertzian contact. 
%
% Inputs:
%   pars - parameters, fields include Rprime and Rpprime
%       pars.Rprime  - Radius R' as defined in Johnson. 
%       pars.Rpprime - Radius R'' as defined in Johnson. R'>R''
%
%       pars may have other fields that are returned unmodified except those
%       under outputs.
% 
% Outputs:
%   pars - parameters with new fields
%       pars.aperb - values of (a/b)
%       pars.e - eccentricity of contact 
%       pars.Ee - First elliptic integral
%       pars.Ke - second elliptic integral
%       pars.F1 - as defined in Johnson
%       pars.F2 - as defined in Johnson
%       pars.Re - as defined in Johnson
%       pars.N  - as defined in Mindlin 1949 (removed)
%
% NOTES:
%   1. The MATLAB function ellipke takes input of e^2 (not e) to be
%   consistent with other references using the elliptic integral functions
%   2. Alternative forms of F1 and F2 could be obtained from Hale 1999 (PhD
%   Thesis), but are not clearly better that just using the equations from
%   Johnson with the approximated eccentricty of contact.
%   3. The tolerances may be too tight for fsolve to get a solution.
%   However, the solutions it returns are generally sufficiently converged.
%   The goal is just to get as good as possible within numerical tolerances
%   by setting the tolerances tighter.


    aperb0 = pars.Rprime/pars.Rpprime;

    % Try to tightnen fsolve tolerances, may not help much.
    fopts = optimoptions('fsolve','Display','iter',...
                                'FunctionTolerance', 1e-20, ...
                                'MaxIterations', 400, ...
                                'MaxFunctionEvaluations', 400, ...
                                'StepTolerance', 1e-20);

    aperb_sol = fsolve(@(x)CONTACT_ECCENTRICITY_RES(x, pars.Rprime, pars.Rpprime), aperb0, fopts);
    
    % Eccentricity
    e = (1 - aperb_sol^(-2))^(1/2);
    
    %Calculate the elliptical values at e
    [K,E] = ellipke(e^2);

    % Outputs
    pars.aperb = aperb_sol; % - values of (a/b)
    pars.e = e; % - eccentricity of contact 
    pars.Ee = E; % - First elliptic integral
    pars.Ke = K; % - second elliptic integral
    
    if(e == 0)
        pars.F1 = 1;
    else
%         pars.F1 = (4/pi/e^2*(aperb_sol^(-3/2)) ...
%                     *( (aperb_sol^2*E - K)*(K-E) ...            
%                       )^(1/2) )^(1/3); % - as defined in Johnson
                  
        % Multiply out exponents form in paper:
        pars.F1 = (4/pi/e^2)^(1/3)*(aperb_sol^(-1/2)) ...
                    *( (aperb_sol^2*E - K)*(K-E)  )^(1/6); % - as defined in Johnson
                  
%     pars.F1 = (1-e^2)^(1/4) * (pars.Rprime/pars.Rpprime)^(1/6) * (4 * (K - E)/(pi*e^2))^(1/3); % Hale 1999, pg430

    end

    pars.F2 = (2/pi)*(aperb_sol^(-1/2))*K/pars.F1; % - as defined in Johnson, pg96
    
    pars.Re = sqrt(pars.Rprime*pars.Rpprime); % - as defined in Johnson

    
    % This is given in Mindlin theory, but should not be used in the final
    % implementation.
%     pars.N = 4*pi*( (2/e- e)*K - 2*E/e);
    
    
    % Alernative Approximations, not used here.
%     pars.F1 = (1-e^2)^(1/4) * (pars.Rprime/pars.Rpprime)^(1/6) * (4 * (K - E)/(pi*e^2))^(1/3); % Hale 1999, pg430
%     pars.F2 = 2/pi*(1 - e^2)^(1/4) * K / pars.F1; % Hale 1999, pg430

end