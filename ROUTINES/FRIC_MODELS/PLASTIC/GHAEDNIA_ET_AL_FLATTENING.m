function [Fn, a, Fe, ae, Fp, ap, deltas_star, Upsilon] = GHAEDNIA_ET_AL_FLATTENING(Es, nus, Ef, nuf, Et, Sys, R, Delta)
% Implementation of the Asperity Flattening Model from:
%       H. Ghaednia, M.R.W. Brake, M. Berryhill, R.L. Jackson, 2019, Strain 
%       Hardening From Elastic-Perfectly Plastic to Perfectly Elastic 
%       Flattening Single Asperity Contact, Journal of Tribology.
%
% Inputs:
%   Es - elastic modulus sphere
%   nus - Poisson Ratio sphere
%   Ef - elastic modulus flat
%   nuf - Poisson Ratio flat
%   Et  - plastic tangent modulus of the sphere
%   Sys - yield strength of solfter material (I believe the model assumes
%           that the flat is fully elastic)
%   R  - effective radius of curvature
%   Delta - current displacement total between sphere and flat.
%
% Outputs:
%   Fn - normal force
%   a  - contact radius
%   Fe - force from the hetzian contact
%   ae - radius from the hetzian contact
%
% NOTES:
%   1. The plastic outputs for the elastic regime are just NaNs, these are
%   obviously not useful answers, but it is outside the relevant regime.

% Eqn 3
E = ( (1 - nus^2)/Es + (1 - nuf^2)/Ef )^(-1);

% Eqn 6
C = 1.295*exp(0.736*nus);

% Eqn 5
deltac = (pi * C * Sys / (2*E))^2 * R;

%%%% Elastic Equations

% Hertzian contact / Eqn 1
ae = sqrt(R*Delta);

% Hertzian Contact - Factor of pi is removed from Eqn 31
Fe = (4/3)*sqrt(R)*E*Delta^(1.5);


% Page 5, top col 2, (i)
deltase_star = E*(1-nus^2)/Es;

Upsilon = 0;

%%%% Elastic/Plastic Model choice:
if(Delta <= 1.9*deltac) %Start Sec 3.2
    %Elastic Regime
        
    a = ae;
    Fn = Fe;
    
    deltas_star = deltase_star;
    
    ap = NaN;
    Fp = NaN;
        
else %Plastic regime
    
    % Eqn 11
    Etstar = ( (1 - nus^2)/Et + (1 - nuf^2)/Ef )^(-1);
    deltasp_star = Etstar*(1-nus^2)/Et;
    
    % Eqn 21
    deltas_star = deltase_star + 0.8*(deltasp_star - deltase_star)*(1 - (Et/Es)^(pi/2))*(1 - exp(-Es/(2.7*Sys)*Delta/R));

    % Eqn 24
    B = 0.14*exp(23*Sys/E);

    % Eqn 23
    ap = R*(pi*C*Sys/(2*E))*sqrt( Delta/deltac * (Delta/deltac / 1.9)^B);

    % Eqn 26 - sqrt in exp() was verified against a draft where the
    % typesetting was not vague.
    % This equation is wrong in the published paper. The authors are
    % contacting the journal to add a note.
%     Upsilon = sqrt(Es/Sys)*(Et/Es)^(1 - Delta/R) / (0.5*( 4-3*exp(-2*sqrt(Es/Sys)*Et/Es ) )*(1-Et/Es) );

    % The correct equation from their implementation is:
    Upsilon = -2*sqrt(Es/Sys)*(Et/Es)^(1 - (Delta - deltac)/R) / (4-3*exp(-2*sqrt(Es/Sys)*Et/Es))/(1-(Et/Es));

    
    % Eqn 25
    a = ae + (ap - ae)*exp(Upsilon);
%     Upsilon
%     exp(Upsilon)
%     [ae, ap, a]
    
    % Eqn 28
    H_Sys = 2.84 - 0.92*(1-cos(pi*a/R));
    
    % Eqn 29
    Fc = (4/3)*(R/E)^2*(C*pi*Sys/2)^3;
    
    % Eqn 27
    Fp = Fc*( exp(-0.25*(Delta/deltac)^(5/12))*(Delta/deltac)^1.5 ...
                +(4*H_Sys/C)*(1-exp(-1/25*(Delta/deltac)^(5/9)))*Delta/deltac );
    
    % Eqn 30
    Fn = Fp + (Fe - Fp)*(1 - exp(-3.3*Et/Es));
    
    
    

end





end