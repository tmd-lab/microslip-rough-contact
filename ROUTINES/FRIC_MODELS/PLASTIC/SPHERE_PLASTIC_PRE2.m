function [fxyn, dfxynduxyn, rq0, tx0, ty0, deltam, Fm, am, ...
    dfduxyn0, dfdfxy0, dfddeltam, varargout] = SPHERE_PLASTIC_PRE2(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq, deltam, Fm, am, varargin)
% Function determines forces for the contact state of an asperity. Only
% considers normal contact for prestressed state. 
%
% Elastic Loading:
%       Hertzian Contact
% Elastic Plastic Loading:
%       H. Ghaednia, M.R.W. Brake, M. Berryhill, R.L. Jackson, 2019, Strain 
%       Hardening From Elastic-Perfectly Plastic to Perfectly Elastic 
%       Flattening Single Asperity Contact, Journal of Tribology.
% Elastic Unloading:
%       I. Etsion, Y. Kligerman and Y. Kadin, 2005, Unloading of an Elastic
%       -Plastic Loaded Spherical Contact, International Journal of Solids 
%       and Structures
% 
% Inputs:
%   pars - structure of material and asperity parameters
%   uxyn - displacements, un is the gap distance with positive being in
%           contact
%   uxyn0 - previous displacements, un0 is again gap distance
%   rq0   - quadrature radii at last instant
%   tx0   - tractions in x direction at rq0 at last instant
%   ty0   - tractions in y direction at rq0 at last instant
%   rq    - quadrature radii to be used normalized to a, so range [0, 1]
%   wq    - weights for the quadrature radii. 
%   deltam - previous maximum displacement
%   Fm - previous maximum normal force
%   am - previous maximum contact radius
%   dFmddeltam = varargin{1} - derivative of Fm w.r.t. deltam is
%                   dfxynduxyn(3,3) evaluated when deltam is initially
%                   reached
%   damddeltam = varargin{2}; - derivative of am w.r.t. deltam, is dadun
%                    when the force is evaluated for the first time at 
%                    deltam
%
% Outputs:
%   fxyn  - forces in xyn directions - zeros in xy directions
%   dfxynduxyn - Jacobian stiffness matrix of fyxn 
%   rq0 - quadrature radii at the end - equals rq
%   tx0 - tractions at each radii after instant - zeros
%   ty0 - tractions at each radii after this instant. - zeros
%   deltam - previous maximum displacement (updated)
%   Fm - previous maximum normal force (updated)
%   dfduxyn0 - 3x3 matrix of derivative of current force w.r.t. previous
%               displacement (zeros for this specific model)
%   dfdfxy0  - 3x2 matrix of derivatives - [dfxdfx0,0; 0,dfydfy0; 0,0]
%               (zeros for this specific model)
%   dfddeltam - derivative w.r.t. deltam with full dependencies on Fm,am
%               3x1 matrix. 
%   varargout{1} - contant dimension a
%   varargout{2} - derivative da/dun
%   varargout{3} - Yield stress (strain hardened)
%   varargout{4} - dYieldStress/dun (0 if not currently yielding).
%   varargout{5} - daddeltam - d(ContactRadius)/d(MaxNormalHistory)
%                   captures all dependencies of deltam,fm,am
%   varargout{6} - dSysddeltam - d(Yield)/d(MaxNormalHistory)
%                   captures all dependencies of deltam,fm,am
%   varargout{7} - deltabar - the point where the unloading separates
%   varargout{8} - ddeltabar_ddeltam - captures all dependencies of deltam,
%                   fm, am
%
% Additional Outputs to add for debugging:
%       , a, dadun, tx, dtxdun, kt, dktdun
%
% NOTES:
%   1. Version 2 of these plasticity functions use this normal force
%   function which differs in that it uses an unloading formulation which
%   has a continuous contact radius. - Etsion, Kligerman, Kadin "Unloading
%   of an elastic-plastic loaded spherical contact"
%
%   2. WARNING: the maximum normal history is taken solely from deltam and
%   completely ignores uxyn0(3). The user is responsible for ensuring that
%   deltam is always appropriately kept up to date, but this current
%   behavior allows for independent derivative tests of everything that is
%   passed out. uxyn0 is included as input here solely to match the call
%   usage of the other plasticity based models.

    
    %Initialize Outputs.
    fxyn = zeros(1,3);
    dfxynduxyn = zeros(3, 3);
    
    % Information for harmonic derivatives - Initialize
    dfduxyn0 = zeros(3, 3); % No dependencies
    dfdfxy0 = zeros(3, 2); % No dependencies
    dfddeltam = zeros(3,1);
    
    if nargin==11
        % These are generally wrong, but are fine for initialization or
        % when doing quasi-static calculations that do not require
        % derivatives w.r.t. perturbations in previous states.
        dFmddeltam = 0;
        damddeltam = 0;
    else
        % These can inputs can be easily found by calling this function
        % with un=deltam and un0 = 0. See initial comments
        dFmddeltam = varargin{1};
        damddeltam = varargin{2};
    end
    
    delta = uxyn(3);
        
    C = 1.295*exp(0.736*pars.nu);
            
    %Different cases for loading regimes
    if(delta < 1.9*pars.delta_y && deltam < 1.9*pars.delta_y) 
        %% Purely Elastic Regime
        
        % Contact Area Elastic
        a = sqrt( (2*pars.R) * (delta/2) );
        dadun = pars.R/(2*a);
        
        
        fxyn(1, 3) = 4*(2*pars.Estar)*sqrt(2*pars.R)/3.*(delta/2).^(3/2).*(delta > 0);
        dfxynduxyn(3,3) = 2*(2*pars.Estar)*sqrt(2*pars.R).*(delta/2).^(1/2)/2.*(delta > 0);
        
%         % Old
%         Sys = pars.Sys;
%         dSysdun = 0;
        
        % Yield stress is approximated based on the average pressure at the
        % transition to the other loading curves
        delta_trans = 1.9*pars.delta_y;
        
        
        a_trans = sqrt( (2*pars.R) * (delta_trans/2) );
        fn_trans = 4*(2*pars.Estar)*sqrt(2*pars.R)/3.*(delta_trans/2).^(3/2).*(delta_trans > 0);

%         Sys = fn_trans / (pi*a_trans^2); error('This is wrong');
%         dSysdun = 0;

        % Non-hysteretic portion of loading, no history dependencies
        daddeltam = 0;
        dSysddeltam = 0;
        dfnddeltam = 0;
        
        [Sys, dSysdun] = calcSys(a_trans, fn_trans, a, dadun, 0, 0, (delta > 0), dFmddeltam, damddeltam, daddeltam);

        % deltabar for the reloading solution offset
        deltabar = 0;
        ddeltabar_ddeltam = 0;

    elseif(delta > deltam) % Elastic-Plastic Loading
        
        %% Contact Area
        
        ae = sqrt( (2*pars.R)*(delta/2) );
        daedun = pars.R/(2*ae);
        
        B = 0.14*exp(23*pars.Sys/(2*pars.Estar));
        
        ap = (2*pars.R)*(pi*C*pars.Sys/(2*2*pars.Estar))*sqrt( delta/pars.delta_y * ( (delta/pars.delta_y)/1.9 )^B );
        dapdun = (2*pars.R)*(pi*C*pars.Sys/(2*2*pars.Estar))...
                    * ( 1/pars.delta_y * ( (delta/pars.delta_y)/1.9)^B ...
                    + 1/pars.delta_y * B*( (delta/pars.delta_y)/1.9)^B  )...
                    /(2*sqrt(delta/pars.delta_y * ( (delta/pars.delta_y)/1.9)^B));
        
        % Eqn 26 - not sure about the square root in the denominator.
        Upsilon = -2*sqrt(pars.E/pars.Sys)*(pars.Et/pars.E)^(1 - ((delta-pars.delta_y)/2)/(pars.R*2)) ...
                    ./ ((4 - 3*exp(-2*sqrt(pars.E/pars.Sys)*pars.Et/pars.E))*(1 - pars.Et/pars.E));
        
        % Derivative from first Upsilon
        dUpsilondun = -2*sqrt(pars.E/pars.Sys) ...
                    ./ ((4 - 3*exp(-2*sqrt(pars.E/pars.Sys)*pars.Et/pars.E))*(1 - pars.Et/pars.E))...
                    *(pars.Et/pars.E)^(1 - ((delta-pars.delta_y)/2)/(pars.R*2)) ...
                    *log(pars.Et/pars.E)*(-1/(4*pars.R));
                
        if(pars.Et == 0)
%             assert( (1 - ((delta-pars.delta_y)/2)/(pars.R*2)) > 0 , 'Did not check limit for other exponents');
            % This derivative is not correct if the exponent goes negative,
            % however that is expected to be rare and the residual function
            % based on the force is not affected.
            dUpsilondun = 0; % Limitting value as Et -> 0
        end
        
        % Eqn 25
        a = ae + (ap - ae)*exp(Upsilon);
        dadun = exp(Upsilon)*(ap-ae)*dUpsilondun + exp(Upsilon)*(dapdun - daedun) + daedun;
        
        %% Force Parameters
        
        H_Sys = 2.84 - 0.92*(1 - cos(pi*a/(2*pars.R)) );
        dH_Sys_dun = -0.92*pi*dadun*sin(pi*a/(2*pars.R))/(2*pars.R);

        Fc = 4/3*(2*pars.R/(2*pars.Estar))^2*(C/2*pi*pars.Sys)^3;
        
        %% Forces

        % Elastic Forces - No idea where the pi came from.
        Fe = 4*(2*pars.Estar)*sqrt(2*pars.R)/3.*(uxyn(3)/2).^(3/2).*(uxyn(3) > 0);
        dFedun = 2*(2*pars.Estar)*sqrt(2*pars.R).*(uxyn(3)/2).^(1/2)/2.*(uxyn(3) > 0);
        
%         % Remove an extra factor of pi, but no longer consistent with the
%         % paper.
%         warning('Removed factor of pi from Hamid''s paper');
%         Fe = 4*(2*pars.Estar)*sqrt(2*pars.R)/3.*(uxyn(3)/2).^(3/2).*(uxyn(3) > 0);
%         dFedun = 2*(2*pars.Estar)*sqrt(2*pars.R).*(uxyn(3)/2).^(1/2)/2.*(uxyn(3) > 0);

        % Plastic Forces
        Fp = Fc*(exp(-0.25*(delta/pars.delta_y)^(5/12))*(delta/pars.delta_y)^(3/2) ...
            + 4*H_Sys/C*(1-exp(-1/25*(delta/pars.delta_y)^(5/9)))*(delta/pars.delta_y));
        
        %Derivative of first term in parens of Fp
        dFp1dun = 3*sqrt(delta/pars.delta_y)*exp(-0.25*(delta/pars.delta_y)^(5/12)) / (2*pars.delta_y) ...
                    - 5*(delta/pars.delta_y)^(11/12)*0.25*exp(-0.25*(delta/pars.delta_y)^(5/12))/(12*pars.delta_y);
        %Second term in Fp
        dFp2dun = 20*delta^(5/9)*H_Sys/25*exp(-(delta/pars.delta_y)^(5/9)/25) / (9*C*pars.delta_y^(14/9))...
                    + 4*delta*dH_Sys_dun*(1-exp(-1/25*(delta/pars.delta_y)^(5/9)))/(C*pars.delta_y) ...
                    + 4*H_Sys*(1-exp(-1/25*(delta/pars.delta_y)^(5/9)))/(C*pars.delta_y);
        
        dFpdun = Fc*(dFp1dun + dFp2dun);
        
        fxyn(1, 3) = Fp + (Fe - Fp)*(1 - exp(-3.3*pars.Et/pars.E));

        dfxynduxyn(3,3) = dFpdun + (dFedun - dFpdun)*(1 - exp(-3.3*pars.Et/pars.E));
    
%         Sys = fxyn(1, 3)/(pi*a^2);
%         dSysdun = 0; %only case that matters is if unloading were to occur next, in which case this would stay contant, except the jump to the other curve.
        
        
        % Yield part of loading. All dependencies on deltam = un are
        % calculated as dependencies on un rather than deltam (no double
        % counting derivatives) thus, all derivatives are zero.
        daddeltam = 0;
        dSysddeltam = 0;
        dfnddeltam = 0;
        
        
        [Sys, dSysdun] = calcSys(a, fxyn(1, 3), a, dadun, dfxynduxyn(3,3), dadun, true, dFmddeltam, damddeltam, daddeltam);

        % deltabar for the reloading solution offset
        deltabar = 0;
        ddeltabar_ddeltam = 0;
        
    else
        %% Elastic Unloading (From Etsion Paper)
        
        
        deltabar = deltam*(1- 3*Fm/(4*pars.Estar*am*deltam));
        ddeltabar_ddeltamOnly = 1;
        ddeltabar_dFm = -3/(4*pars.Estar*am);
        ddeltabar_dam = 3*Fm/(4*pars.Estar*am^2);
        
        ddeltabar_ddeltamFull = ddeltabar_ddeltamOnly + ddeltabar_dFm*dFmddeltam + ddeltabar_dam*damddeltam;
        
        
        Rebar = 4*pars.Estar*(am^3)/(3*Fm);
%         dRebar_ddeltamOnly = 0; %direct dependency = None
        dRebar_dFm = -4*pars.Estar*(am^3)/(3*Fm^2);
        dRebar_dam = 4*pars.Estar*(am^2)/Fm;
        
        dRebar_ddeltamFull = dRebar_dFm*dFmddeltam + dRebar_dam*damddeltam;
        
        
        a = sqrt(Rebar*(delta - deltabar)).*(delta >= deltabar);
        dadun = 0;
        if(delta >= deltabar)
            dadun = Rebar/(2*a);
        end
        daddeltabar = -Rebar/(2*sqrt(Rebar*(delta - deltabar)) ).*(delta >= deltabar);
        dadRebar = (delta - deltabar)/( 2*sqrt(Rebar*(delta - deltabar)) ).*(delta >= deltabar);
        
        
        fxyn(1, 3) = 4/3*pars.Estar*sqrt(Rebar)*(delta - deltabar).^(3/2).*((delta - deltabar) >= 0);
        dfxynduxyn(3,3) = 2*pars.Estar*sqrt(Rebar)*(delta - deltabar).^(1/2).*((delta - deltabar) >= 0);
        dfnddeltabar = -2*pars.Estar*sqrt(Rebar)*(delta - deltabar).^(1/2).*((delta - deltabar) >= 0);
        dfndRebar = 2/3*pars.Estar/sqrt(Rebar)*(delta - deltabar).^(3/2).*((delta - deltabar) >= 0);
        
        
        % Derivatives w.r.t. History variables (deltam and all Fm(deltam)
        % and am(deltam)).
        daddeltam = daddeltabar*ddeltabar_ddeltamFull + dadRebar*dRebar_ddeltamFull;

        dfnddeltam = dfnddeltabar*ddeltabar_ddeltamFull + dfndRebar*dRebar_ddeltamFull;
        
        % Yield stress and appropriate derivatives.
        [Sys, dSysdun, dSysddeltam] = calcSys(am, Fm, a, dadun, 0, 0, (delta - deltabar) >= 0, dFmddeltam, damddeltam, daddeltam);
        
%         deltabar % Already defined
        ddeltabar_ddeltam = ddeltabar_ddeltamFull;
          
    end
    
    % Final derivative of force w.r.t. deltam
    dfddeltam(3) = dfnddeltam;
    
    %% Update History
    if(delta > deltam)
        
        deltam = delta;
        Fm = fxyn(1, 3);
        am = a;
        
    end
    
    %% Other outputs of interest
    if(nargout > 7)
        
        varargout{1} = a;
        varargout{2} = dadun;
        
%         warning('Wrong extra outputs'); %Useful for debugging
%         varargout{1} = Fe;
%         varargout{2} = dFedun;
%         
%         varargout{1} = Fp;
%         varargout{2} = dFpdun;
        
        %Things related to the normal pressure distribution
        
        % Strain hardened yield stress (if 1 quadrature radius output
        % average, if multiple output at quad radii).
        if(length(rq) == 1)
            
            varargout{3} = Sys;
            varargout{4} = dSysdun;
            
        else
            error('hardened yield stress not implemented at specific points yet');
        end
        
        varargout{5} = daddeltam;
        varargout{6} = dSysddeltam;
        
        varargout{7} = deltabar; % the point where the unloading separates
        varargout{8} = ddeltabar_ddeltam; % captures all dependencies of deltam,fm,am
        
    end
    
    %% Final Details for outputting previous information
    rq0 = rq;
    tx0 = zeros(size(rq));
    ty0 = zeros(size(rq));
    
%     tx = tx0;
    
end

function [Sys, dSysdun, dSysddeltam] = calcSys(am, Fm, a, dadun, dFmdun, damdun, incontact, dFmddeltam, damddeltam, daddeltam)
% Function for calculating the yield stress at the current point
% Inputs:
%   am - reversal contact radius
%   Fm - reversal force
%   a  - current contact radius
%   dadun - derivative of current radius w.r.t. un
%   incontact - boolean for in contact
%   dFmddeltam - see top of file
%   damddeltam - see top of file
%   daddeltam - derivative of current contact radius w.r.t. all dependency
%               on deltam (including through Fm/am)
%
% Outputs:
%   Sys - yield stress
%   dSysdun - derivative of yield stress w.r.t. normal displacement
%   dSysddeltam - history dependence derivative w.r.t. all dependencies on
%                   deltam including Fm(deltam) and am(deltam)

    
    if(a > am)
        if(abs(a - am)/a < 1e-10)
            %On the unloading curve in general with the recalculation of
            %contact radius can have some numerical error
            a = am;
        else
            error('In normal loading a > am, error of %f.', abs(a - am)/a);
        end
        
    end
    
    
    if(incontact)
        
        Sys = Fm/(pi*a^2)*(1 - (1 - (a/am)^2)^1.5);

        dSysda = -2*Fm/(pi*a^3)*(1 - (1 - (a/am)^2)^1.5) + 3*Fm/(pi*a*am^2)*( (1 - (a/am)^2)^0.5);
%         dSysddeltam = 0 %No Dependency
        dSysdFm = 1/(pi*a^2)*(1 - (1 - (a/am)^2)^1.5);
        dSysdam = -3*Fm/(pi*am^3)*(1 - (a/am)^2)^0.5;
        
        dSysdun = dSysda*dadun + dSysdFm*dFmdun + dSysdam*damdun;
        dSysddeltam = dSysdFm*dFmddeltam + dSysdam*damddeltam + dSysda*daddeltam;

    else
        
        Sys = 3*Fm/(2*pi*am^2); %Peak normal pressure at center
        dSysdun = 0;
        
%         dSysddeltam = 0 %No Dependency
        dSysdFm = 3/(2*pi*am^2);
        dSysdam = -3*Fm/(pi*am^3);
        
        dSysddeltam = dSysdFm*dFmddeltam + dSysdam*damddeltam;
        
    end
   
    
end
