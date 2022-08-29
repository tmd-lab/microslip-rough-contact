function [fxyn, dfxynduxyn, rq0, tx0, ty0] = ELLIPSOID_IWAN_FIT_COUPLE(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq)
% Function determines forces for the contact state of an asperity. Uses an
% Iwan model fit to the Mindlin solution. Couples x/y directions at each
% radius
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


    %% Normal Forces + derivatives
    
    [fxyn, dfxynduxyn, ~, ~, ~, b, dbdun, a, dadun] ...
        = ELLIPSOID_PRE(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq);
    
    %Slip force on reference domain of the unit circle
    ts = pars.mu*3*fxyn(1, 3)/(2*pi*a*b)*sqrt(1 - rq.^2);
    
    %consider replacing ab with b^2*pars.aperb in the denominator yielding
    dtsdun = pars.mu*3*dfxynduxyn(3, 3)/(2*pi*a*b)*sqrt(1 - rq.^2) ...
                -2*pars.mu*3*fxyn(1, 3)/(2*pi*b^3*pars.aperb)*sqrt(1 - rq.^2)*dbdun;
    
    %% Tangential Stiffnesses - Tangent Form
    
    %Model fit kt
    pars_opt = [0.8709    0.0629   -0.8915    0.6998];
    ktfun = @(pars, rq) pars(1) + (1- pars(1))./sqrt(1 - rq.^2) + pars(2)*rq + pars(3)*rq.^2 + pars(4)*rq.^3;

    %evaluate kt_tilde
    kttilde = ktfun(pars_opt, rq); %kt_tilde is fit on a unit circle, so return rq to the unit circle.
    
    kttilde(end) = 1; %At rq == a, the slip traction is always zero, so it doesn't matter what kt is just needs to not be NaN.
    
    %Assemble dimensional kt values
    % In paper, Phix = Phib (aligned with short semi-axis)
    %           Phiy = Phia (aligned with long semi-axis)
    Phix = 4/(pi*(2-pars.nu))*((1 - pars.nu/pars.e^2)*pars.Ke + pars.nu*pars.Ee/pars.e^2);
    Phiy = 4/(pi*(2-pars.nu))*( (1 - pars.nu + pars.nu/pars.e^2)*pars.Ke - pars.nu*pars.Ee/pars.e^2);
    
    if(pars.e == 0)
        Phix = 1;
        Phiy = 1;
    end
    
    ktx = 8*pars.Gstar/(Phix*pi*b)*kttilde;
    
    dktxdun = -8*pars.Gstar/(Phix*pi*b^2)*dbdun*kttilde;
    
    kty = 8*pars.Gstar/(Phiy*pi*b)*kttilde;

    dktydun = -8*pars.Gstar/(Phiy*pi*b^2)*dbdun*kttilde;
    
    if( b == 0 )
        
        % Remove NaN's at zero displacement
        ktx = zeros(size(ktx));
        dktxdun = zeros(size(dktxdun));
        kty = zeros(size(kty));
        dktydun = zeros(size(dktydun));
        
        ts = zeros(size(ts));
        dtsdun = zeros(size(dtsdun));
    end

    %% Re-establishment of contact conditions 
    % One could apply this to each radii instead of just the full asperity
    % having new contact, but that is a complicated nonlinear equation for
    % a(w). 
    
    %initialize derivatives of previous information w.r.t. new disps.
    duxyn0duxyn = zeros(3,3); %Independent unless new contact
    
    if(uxyn0(1, 3) < 0)
        
        %Derivatives of previous displacements - before overwriting the
        %uxyn0
        duxyn0duxyn(1, 1) = (0 - uxyn0(3))./(uxyn(3) - uxyn0(3));
        duxyn0duxyn(2, 2) = (0 - uxyn0(3))./(uxyn(3) - uxyn0(3));
        duxyn0duxyn(1:2, 3) = (uxyn0(3))./(uxyn(3) - uxyn0(3))^2.*(uxyn(1:2)' - uxyn0(1:2)');
        
        %now update uxyn0
        uxyn0(1, 1:2) = (0 - uxyn0(3))./(uxyn(3) - uxyn0(3)).*(uxyn(1:2) - uxyn0(1:2)) + uxyn0(1:2);
%         uxyn0(1, 3) = 0; %Should not matter.

        %new contact, so only rq0 is meaningless, fix so interpolation
        %doesn't error out. 
        rq0 = rq;
        
    end
    
    %% Interpolate Previous Information to Current 
    % (X & Y) should share the same boolean matrix, but differ in the slope
    % calculation.
    
    % %% Interpolation to new quadrature radii
    
    %set to zero if new contact location.
    old_contact = (rq*a+eps <= rq0(end)) & (a ~= 0);
    
    % Derivatives of interpolation calculated with table look up
    % Next and previous return the same value if rq*a == rq0 at a point, so
    % perturb by eps
    tx0_high = interp1(rq0, tx0, rq(old_contact)*a+eps, 'next');
    tx0_low  = interp1(rq0, tx0, rq(old_contact)*a-eps, 'previous');
    
    ty0_high = interp1(rq0, ty0, rq(old_contact)*a+eps, 'next');
    ty0_low  = interp1(rq0, ty0, rq(old_contact)*a-eps, 'previous');
    
    rq0_high = interp1(rq0, rq0, rq(old_contact)*a+eps, 'next');
    rq0_low  = interp1(rq0, rq0, rq(old_contact)*a-eps, 'previous');

    % The first point is always at 0 radius, so needs to be replaced to
    % eliminate NaN's for the case of -eps
    tx0_low(1) = tx0(1);
    ty0_low(1) = ty0(1);
    rq0_low(1) = rq0(1);
    
    % Initialize storage for derivatives, only nonzero when previous
    % contact
    dtx0da = zeros(size(tx0));
    dty0da = zeros(size(ty0));
    
    dtx0da(old_contact) = ((tx0_high - tx0_low)./(rq0_high-rq0_low)).*rq(old_contact);
    dty0da(old_contact) = ((ty0_high - ty0_low)./(rq0_high-rq0_low)).*rq(old_contact);
    
    % Update interpolation last overwriting data.
    %Interpolate data to new radii
    tx0(old_contact) = interp1(rq0, tx0, rq(old_contact)*a);
    ty0(old_contact) = interp1(rq0, ty0, rq(old_contact)*a);
    
    tx0(~old_contact) = 0;
    ty0(~old_contact) = 0;
    
    %% Stick Predictions
    
    tx_stick = tx0 + ktx.*(uxyn(1) - uxyn0(1));
    ty_stick = ty0 + kty.*(uxyn(2) - uxyn0(2));
    
    txynorm = sqrt((tx_stick).^2+(ty_stick).^2);
    
    % this is strictly greater to prevent NaNs, but this results in
    % including kt in the derivative at r=a in the derivative even though
    % it can never actually be included because ts is always exactly zero
    % at that point (normal derivative can be nonzero, tangential should be
    % zero).
    slipped = txynorm > ts;
    
%     % Some condition to prevent divide by zeros?
%     txynorm(abs(txynorm)<0.01*ts) = 1.0;
    
    tx_normed = tx_stick./txynorm; 
    ty_normed = ty_stick./txynorm;  
    
    %% Derivatives of stick/norm
    
    % Stick X
    dtx_stick_dux = (ktx - ktx.*duxyn0duxyn(1,1));
    
    dtx_stick_dun = dktxdun.*(uxyn(1) - uxyn0(1)) ...
                        + dadun.*dtx0da ...
                        - ktx.*duxyn0duxyn(1, 3);
    
    % Stick Y
    dty_stick_duy = (kty - kty.*duxyn0duxyn(2,2));
    
    dty_stick_dun = dktydun.*(uxyn(2) - uxyn0(2)) ...
                        + dadun.*dty0da ...
                        - kty.*duxyn0duxyn(2, 3);
                    
    % Norm Derivative - only of slipped points
    dtxynorm_dux = tx_stick(slipped).*dtx_stick_dux(slipped)./txynorm(slipped);
    dtxynorm_duy = ty_stick(slipped).*dty_stick_duy(slipped)./txynorm(slipped);
    
    
    dtxynorm_dun = (tx_stick(slipped).*dtx_stick_dun(slipped) ...
                    + ty_stick(slipped).*dty_stick_dun(slipped))./txynorm(slipped);
                
    % Derivatives of norm X
    dtx_normed_dux = dtx_stick_dux(slipped)./txynorm(slipped) - tx_stick(slipped)./txynorm(slipped).^2.*dtxynorm_dux;
    dtx_normed_duy = - tx_stick(slipped)./txynorm(slipped).^2.*dtxynorm_duy;
    dtx_normed_dun = dtx_stick_dun(slipped)./txynorm(slipped) - tx_stick(slipped)./txynorm(slipped).^2.*dtxynorm_dun;
    
    
    % Derivatives of norm Y
    dty_normed_duy = dty_stick_duy(slipped)./txynorm(slipped) - ty_stick(slipped)./txynorm(slipped).^2.*dtxynorm_duy;
    dty_normed_dux = - ty_stick(slipped)./txynorm(slipped).^2.*dtxynorm_dux;
    dty_normed_dun = dty_stick_dun(slipped)./txynorm(slipped) - ty_stick(slipped)./txynorm(slipped).^2.*dtxynorm_dun;

    %% Tractions in Each Direction
    
    tx0 = tx_stick.*(~slipped);
    ty0 = ty_stick.*(~slipped);
    
    tx0(slipped) = tx_normed(slipped).*ts(slipped);
    ty0(slipped) = ty_normed(slipped).*ts(slipped);
    
    %% Derivatives of Tractions
    
    % X Direction
    dtxdux = dtx_stick_dux.*(~slipped);
    dtxduy = zeros(size(slipped));
    dtxdun = dtx_stick_dun.*(~slipped);
    
    % Account for rq=a is always slipped, but strictly greater is used in slip
    % condition to ensure that NaN's are avoided.
    dtxdux(end) = 0;
    
    dtxdux(slipped) = dtx_normed_dux.*ts(slipped);
    dtxduy(slipped) = dtx_normed_duy.*ts(slipped);
    dtxdun(slipped) = dtx_normed_dun.*ts(slipped) + tx_normed(slipped).*dtsdun(slipped);
    
    % Y Direction
    dtyduy = dty_stick_duy.*(~slipped);
    dtydux = zeros(size(slipped));
    dtydun = dty_stick_dun.*(~slipped);
    
    % Account for rq=a is always slipped, but strictly greater is used in slip
    % condition to ensure that NaN's are avoided.
    dtyduy(end) = 0;
    
    dtyduy(slipped) = dty_normed_duy.*ts(slipped);
    dtydux(slipped) = dty_normed_dux.*ts(slipped);
    dtydun(slipped) = dty_normed_dun.*ts(slipped) + ty_normed(slipped).*dtsdun(slipped);
    
    
    
    %% Integrate X
    % %%integrate forces and derivatives of forces
    fxyn(1,1) = 2*pi*a*b*(tx0.*rq)*wq';
    
    dfxynduxyn(1,1) = 2*pi*a*b*(dtxdux.*rq)*wq';
        
    dfxynduxyn(1,2) = 2*pi*a*b*(dtxduy.*rq)*wq';

    dfxynduxyn(1, 3) = 2*pi*a*b*(dtxdun.*rq)*wq' ... 
                        + 4*pi*pars.aperb*b*(tx0.*rq)*wq'*dbdun; %Combine leading to b^2*pars.aperb + derivative

        
    %% Integrate Y
    % %%integrate forces and derivatives of forces
    fxyn(1,2) = 2*pi*a*b*(ty0.*rq)*wq';
    
    dfxynduxyn(2,2) = 2*pi*a*b*(dtyduy.*rq)*wq';
    
    dfxynduxyn(2,1) = 2*pi*a*b*(dtydux.*rq)*wq';

    dfxynduxyn(2,3) = 2*pi*a*b*(dtydun.*rq)*wq' ... 
                        + 4*pi*pars.aperb*b*(ty0.*rq)*wq'*dbdun; %Combine leading to b^2*pars.aperb + derivative


    %% Update outputs
    
    rq0 = rq*a;
    
end