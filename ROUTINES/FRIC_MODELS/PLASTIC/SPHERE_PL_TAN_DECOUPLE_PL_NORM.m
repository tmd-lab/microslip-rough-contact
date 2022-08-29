function [fxyn, dfxynduxyn, rq0, fx0, fy0, deltam, Fm, am, ...
    dfduxyn0, dfdfxy0, dfddeltam, a, dadun] = SPHERE_PL_TAN_DECOUPLE_PL_NORM(pars, uxyn, uxyn0, rq0, fx0, fy0, rq, wq, deltam, Fm, am, varargin)
% Function determines forces for the contact state of an asperity. 
% Plastic normal direction, elastic tangential directions
% 
% Inputs:
%   pars - structure of material and asperity parameters
%   uxyn - displacements, un is the gap distance with positive being in
%           contact
%   uxyn0 - previous displacements, un0 is again gap distance
%   rq0   - quadrature radii at last instant, ignored
%   tx0   - total force in x direction at last instant
%   ty0   - total force in y direction at last instant
%   rq    - quadrature radii to be used normalized to a, so range [0, 1]
%   wq    - weights for the quadrature radii. 
%   deltam - maximum displacement that has occured in the normal direction
%   Fm - maximum normal force that has occured
%   am - contact radius at the maximum past normal displacement/force
%   dFmddeltam = varargin{1} - derivative of Fm w.r.t. deltam is
%                   dfxynduxyn(3,3) evaluated when deltam is initially
%                   reached
%   damddeltam = varargin{2}; - derivative of am w.r.t. deltam, is dadun
%                    when the force is evaluated for the first time at 
%                    deltam
%
% Outputs:
%   fxyn  - forces in xyn directions
%   dfxynduxyn - Jacobian stiffness matrix of fyxn
%   rq0 - quadrature radii at the end
%   fx0 - total force in x direction
%   fy0 - total force in y direction
%   deltam - maximum displacement that has occured in the normal direction
%   Fm - maximum normal force that has occured
%   am - contact radius at the maximum past normal displacement/force
%   dfduxyn0 - 3x3 matrix of derivative of current force w.r.t. previous
%               displacement (zeros for this specific model)
%   dfdfxy0  - 3x2 matrix of derivatives - [dfxdfx0,0; 0,dfydfy0; 0,0]
%               (zeros for this specific model)
%   dfddeltam - derivative w.r.t. deltam with full dependencies on Fm,am
%               3x1 matrix. 
%   a - contact radius, passed out here so that output positions match
%       SPHERE_PLASTIC_PRE2 if/when dadun is needed.
%   dadun - derivative of contact radius, passed out here in case it is
%           needed to initialize damddeltam at some point.
%
% Additional Outputs to add for debugging:
%       , a, dadun, tx, dtxdun, kt, dktdun
%
% NOTES:
%   1. Since the derivatives returned from SPHERE_PLASTIC_PRE include all
%   dependencies on deltam (through deltam, fm, am), every derivative in
%   this function w.r.t. deltam also contains all of these dependencies in
%   one term.
%   2. The derivatives of tangential information w.r.t. normal information
%   are ill defined when ts = 0 due to plasticity. They appear to have
%   infinite derivatives in one direction and 0 in the other. These cases
%   just return 0's since they play nice with everything else, but they are
%   only correct in one direction. The 0 derivative direction is not
%   consistent with the derivative direction returned from
%   SPHERE_PLASTIC_PRE for dfndun, so it is very difficult to verify these
%   derivatives. 

    %% Handle varargin
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
    

    %% Re-establishment of contact conditions 
    
    %initialize derivatives of previous information w.r.t. new disps.
    duxyn0duxyn = zeros(3,3); %Independent unless new contact
    
    duxyn0duxyn0 = eye(3); % prev Displacements vary exactly the same when prev is in contact
    
    % This section should do an update as:
    if(uxyn0(1, 3) < 0)
        
        % Other asperity models use a re-establishment of contact method to
        % update uxyn0. This is omitted due to the complexity of harmonic
        % derivatives and the residual displacement of the sphere.

        fx0(1) = 0;
        fy0(1) = 0;
        
    end
    

    %% Normal force solution
    
    [fxyn, dfxynduxyn, ~, ~, ~, deltam, Fm, am, ...
    dfduxyn0, dfdfxy0, dfddeltam, a, dadun, Sys, ...
    dSysdun, daddeltam, dSysddeltam, deltabar, ddeltabar_ddeltam, ...
    Rebar, dRebar_ddeltam] ...
        = SPHERE_PLASTIC_PRE(pars, uxyn, uxyn0, rq0, fx0, fy0, rq, wq, deltam, Fm, am, dFmddeltam, damddeltam);

    
    % Basic Plasticity Model:
    % This is an ad hoc calculation. Using the CEB model provides a more
    % rigorous was of estimating the plasticity friction coefficient. 
    %
    % Consider average normal pressure and average shear stress in the
    % limit of yielding for f_s + strain hardening of surface
    tn = fxyn(3)/(pi*a^2);
    dtndun = -2*fxyn(3)/(pi*a^3)*dadun+dfxynduxyn(3,3)/(pi*a^2);
    
    dtnddeltam = -2*fxyn(3)/(pi*a^3)*daddeltam+dfddeltam(3,1)/(pi*a^2);
    
    % Strain hardened yield stress
    ts = sqrt((Sys^2 - tn^2)/3);
    dtsdun = 1/(2*ts*3)*(2*Sys*dSysdun-2*tn*dtndun);
    
    dtsddeltam = 1/(2*ts*3)*(2*Sys*dSysddeltam-2*tn*dtnddeltam);
    
    if(ts < pars.mu*tn && mod(pars.sliptype, 2) == 0 ) % NaN should correspond to setting Sys = inf to force a switch to mu. 
        %plastic limit
    
        if(uxyn(3) == deltam) % Better way to check that ts == 0 should be true (except this can actually happen)
            
            %This occurs because yielding has occured setting Sys == tn
            ts = 0;

            %Further increasing the normal load causes ts to track Sys and
            %therefore stay zero so the derivative in one direction is exactly
            %zero. If unloading, the derivative may be nonzero, but being
            %consistent in one direction is good enough. 

            dtsdun = 0;
            dtsddeltam = 0;
        end

        if(~isreal(ts))
            ts = 0;
            dtsdun = 0;
            dtsddeltam = 0;
        end

        fs = ts*(pi*a^2);
        dfsdun = ts*(2*pi*a*dadun) + dtsdun*(pi*a^2);
        
        dfsddeltam = ts*(2*pi*a*daddeltam) + dtsddeltam*(pi*a^2);
    else
        
        % Coulomb limit
        fs          = pars.mu*fxyn(1,3);
        dfsdun      = pars.mu*dfxynduxyn(3,3);
        dfsddeltam  = pars.mu*dfddeltam(3,1);
    end
    
    %% CEB Limit for force
    % This comes from:
    %       Eriten et al., 2010, Physics-based modeling for partial slip
    %       behavior of spherical contacts
    %
    % Some assumptions on the application here differ in that the yield
    % displacement is taken to be deltam since that is when the sphere
    % would yield again according to the hardening model 
    
    if(mod(pars.sliptype, 3) == 0)
        
        % Determine CEB Friction Coefficient + Derivative
        [mu_ceb, dmudun, dmudunmax] = CEB_MU(pars, uxyn(3)-deltabar, ...
                        (deltam-deltabar).*(deltam >= 1.9*pars.delta_y) + 1.9*pars.delta_y.*(deltam < 1.9*pars.delta_y));
        
        if(mu_ceb < pars.mu)
            % Force Limit
            fs          = mu_ceb*fxyn(1,3);
            dfsdun      = dmudun*fxyn(1,3) + mu_ceb*dfxynduxyn(3,3);
            
            % dmu(un(deltam), deltam-deltabar)/ddeltam * f + mu * dfddeltam + dmu(un(deltam))/ddeltam*f
            % Last one is because dun/ddeltam = -ddeltabar when CEB is evaluated.
            dfsddeltam  = dmudunmax*(1-ddeltabar_ddeltam)*fxyn(1,3).*(deltam >= 1.9*pars.delta_y)...
                            + mu_ceb*dfddeltam(3,1) - dmudun*ddeltabar_ddeltam*fxyn(1,3);
            
            % The present implementation of the previous two lines ignores
            % the case of the current normal interference setting deltam in
            % which case un=deltam. If that was the case, the derivatives
            % above would need to be combined or at least add:
            % dfsdun = dfsdun + dfsddeltam;
            %
            % However, I handle deltam initialization separately in the
            % Harmonic code, so this should be correct in that all parts of
            % the derivative get counted exactly once.
        end
    end
    
    %% Repeat CEB Limit, but with maximum of UTS
    
    if(mod(pars.sliptype, 5) == 0)
        
        [deltay_UTS, ddeltayddeltam_UTS] = UTS_CRIT(pars, Rebar, dRebar_ddeltam);
        
        
        % Determine CEB Friction Coefficient + Derivative
        [mu_ceb, dmudun, dmudunmax] = CEB_MU(pars, uxyn(3)-deltabar, deltay_UTS);
        
        if( (mu_ceb*fxyn(1,3)) < fs)
            
            % Force Limit
            fs          = mu_ceb*fxyn(1,3);
            dfsdun      = dmudun*fxyn(1,3) + mu_ceb*dfxynduxyn(3,3);
            
            dfsddeltam  = dmudunmax*ddeltayddeltam_UTS...
                            + mu_ceb*dfddeltam(3,1) - dmudun*ddeltabar_ddeltam*fxyn(1,3);
            
            % The present implementation of the previous two lines ignores
            % the case of the current normal interference setting deltam in
            % which case un=deltam. If that was the case, the derivatives
            % above would need to be combined or at least add:
            % dfsdun = dfsdun + dfsddeltam;
            %
            % However, I handle deltam initialization separately in the
            % Harmonic code, so this should be correct in that all parts of
            % the derivative get counted exactly once.
        end
    end
    
    %% Contact Stiffnesses
    
    %The assumption of spheres uses the same values of kt in both
    %directions.
    
    kt = 8.*pars.Gstar*a;
    dktdun = 8.*pars.Gstar.*dadun;

    dktddeltam = 8.*pars.Gstar.*daddeltam;

    %% Tangential Contact - X
    
    % %% Calculation of traction
    
    %stick and slip forces
    tx_stick = fx0(1) + kt.*(uxyn(1) - uxyn0(1));
    tx_slip  = fs.*(sign(tx_stick) + (sign(tx_stick) == 0)); %if stick is zero, just set sign to positive.

    fx0 = tx_slip.*(abs(tx_stick)>=abs(tx_slip)) + tx_stick.*(abs(tx_stick)<abs(tx_slip));
    
    %just the traction derivative.
    dtxdux = (kt - kt*duxyn0duxyn(1,1)).*(abs(tx_stick)<abs(tx_slip));

    %derivatives w.r.t. normal coordinates.
    dtxdun = dktdun.*(uxyn(1) - uxyn0(1)).*(abs(tx_stick)<abs(tx_slip)) ...
        - kt*duxyn0duxyn(1,3).*(abs(tx_stick)<abs(tx_slip)) ...
        + dfsdun.*(sign(tx_stick) + (sign(tx_stick) == 0)).*(abs(tx_stick)>=abs(tx_slip));
        
    % derivatives w.r.t. uxyn0
    dtxdux0 = (-kt*duxyn0duxyn0(1,1)).*(abs(tx_stick)<abs(tx_slip));
    dtxdun0 = -kt*duxyn0duxyn0(1,3).*(abs(tx_stick)<abs(tx_slip));
    dtxdfx0 = (abs(tx_stick)<abs(tx_slip));
    
    % %%integrate forces and derivatives of forces
    fxyn(1,1) = fx0;
    
    dfxynduxyn(1,1) = dtxdux;
    dfxynduxyn(1,3) = dtxdun;

    dfduxyn0(1,1) = dtxdux0;
    dfduxyn0(1,3) = dtxdun0;
    
    dfdfxy0(1,1) = dtxdfx0;
    
    dfddeltam(1,1) = dktddeltam.*(uxyn(1) - uxyn0(1)).*(abs(tx_stick)<abs(tx_slip)) ...
                     + dfsddeltam.*(sign(tx_stick) + (sign(tx_stick) == 0)).*(abs(tx_stick)>=abs(tx_slip));
    
    %% Tangential Contact - Y
    
    % %% Calculation of traction
    
    %stick and slip forces
    ty_stick = fy0(1) + kt.*(uxyn(2) - uxyn0(2));
    ty_slip  = fs.*(sign(ty_stick) + (sign(ty_stick) == 0)); %if stick is zero, just set sign to positive.

    fy0 = ty_slip.*(abs(ty_stick)>=abs(ty_slip)) + ty_stick.*(abs(ty_stick)<abs(ty_slip));
    
    %just the traction derivative.
    dtyduy = (kt - kt*duxyn0duxyn(2,2)).*(abs(ty_stick)<abs(ty_slip));

    %derivatives w.r.t. normal coordinates.
    dtydun = dktdun.*(uxyn(2) - uxyn0(2)).*(abs(ty_stick)<abs(ty_slip)) ...
        - kt*duxyn0duxyn(2,3).*(abs(ty_stick)<abs(ty_slip)) ...
        + dfsdun.*(sign(ty_stick) + (sign(ty_stick) == 0)).*(abs(ty_stick)>=abs(ty_slip));
    
    % derivatives w.r.t. uxyn0
    dtyduy0 = (-kt*duxyn0duxyn0(2,2)).*(abs(ty_stick)<abs(ty_slip));
    dtydun0 = -kt*duxyn0duxyn0(2,3).*(abs(ty_stick)<abs(ty_slip));
    dtydfy0 = (abs(ty_stick)<abs(ty_slip));
    
    % %%integrate forces and derivatives of forces
    fxyn(1,2) = fy0;
    dfxynduxyn(2,2) = dtyduy;

    dfxynduxyn(2, 3) = dtydun;


    dfduxyn0(2,2) = dtyduy0;
    dfduxyn0(2,3) = dtydun0;
    
    dfdfxy0(2,2) = dtydfy0;
    
    dfddeltam(2,1) = dktddeltam.*(uxyn(2) - uxyn0(2)).*(abs(ty_stick)<abs(ty_slip)) ...
                     + dfsddeltam.*(sign(ty_stick) + (sign(ty_stick) == 0)).*(abs(ty_stick)>=abs(ty_slip));

    
    %% Final Details for outputting previous information
    rq0 = 0;
    
    %% Check for separation
    if(fxyn(1,3) == 0)
        
        fxyn(1:2) = 0;
        dfxynduxyn = zeros(3,3);
        
        fx0 = 0;
        fy0 = 0;
        dfduxyn0 = zeros(3,3);
        dfdfxy0 = zeros(3,2);
        dfddeltam = zeros(3,1);
        
        a = 0;
        dadun = 0;
    end
    
end