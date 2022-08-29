function [fxyn, dfxynduxyn, rq0, tx0, ty0] = ELLIPSOID_IWAN_DECOUPLE(pars, uxyn, uxyn0, rq0, tx0, ty0, rq, wq)
% Function determines forces for the contact state of an asperity. Iwan
% with constant distributed stiffness to approximate Mindlin solution
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
    
    % In paper, Phix = Phib (aligned with short semi-axis)
    %           Phiy = Phia (aligned with long semi-axis)
    Phix = 4/(pi*(2-pars.nu))*((1 - pars.nu/pars.e^2)*pars.Ke + pars.nu*pars.Ee/pars.e^2);
    Phiy = 4/(pi*(2-pars.nu))*( (1 - pars.nu + pars.nu/pars.e^2)*pars.Ke - pars.nu*pars.Ee/pars.e^2);
    
    if(pars.e == 0)
        Phix = 1;
        Phiy = 1;
    end
    
    ktx = 8*pars.Gstar/(Phix*pi*b);
    
    dktxdun = -8*pars.Gstar/(Phix*pi*b^2)*dbdun;
    
    kty = 8*pars.Gstar/(Phiy*pi*b);

    dktydun = -8*pars.Gstar/(Phiy*pi*b^2)*dbdun;
    
    if( b == 0 )
        % Remove NaN's at zero displacement
        ktx = 0;
        dktxdun = 0;
        kty = 0;
        dktydun = 0;
        
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
    
    %% Tangential Contact - X

    % %% Calculation of traction
    
    %stick and slip forces
    tx_stick = tx0 + ktx.*(uxyn(1) - uxyn0(1));
    tx_slip  = ts.*(sign(tx_stick) + (sign(tx_stick) == 0)); %if stick is zero, just set sign to positive.

    tx0 = tx_slip.*(abs(tx_stick)>=abs(tx_slip)) + tx_stick.*(abs(tx_stick)<abs(tx_slip));
    
    %just the traction derivative.
    dtxdux = (ktx - ktx.*duxyn0duxyn(1,1)).*(abs(tx_stick)<abs(tx_slip));

    %derivatives w.r.t. normal coordinates.
    dtxdun = dktxdun.*(uxyn(1) - uxyn0(1)).*(abs(tx_stick)<abs(tx_slip)) ...
        + dadun.*dtx0da.*(abs(tx_stick)<abs(tx_slip)) ...
        - ktx.*duxyn0duxyn(1, 3).*(abs(tx_stick)<abs(tx_slip)) ...
        + dtsdun.*(sign(tx_stick) + (sign(tx_stick) == 0)).*(abs(tx_stick)>=abs(tx_slip));
        
    % %%integrate forces and derivatives of forces
    fxyn(1,1) = 2*pi*a*b*(tx0.*rq)*wq';
    dfxynduxyn(1,1) = 2*pi*a*b*(dtxdux.*rq)*wq';

    dfxynduxyn(1, 3) = 2*pi*a*b*(dtxdun.*rq)*wq' ... 
                        + 4*pi*pars.aperb*b*(tx0.*rq)*wq'*dbdun; %Combine leading to b^2*pars.aperb + derivative

    %% Tangential Contact - Y

    % %% Calculation of traction
    
    %stick and slip forces
    ty_stick = ty0 + kty.*(uxyn(2) - uxyn0(2));
    ty_slip  = ts.*(sign(ty_stick) + (sign(ty_stick) == 0)); %if stick is zero, just set sign to positive.

    ty0 = ty_slip.*(abs(ty_stick)>=abs(ty_slip)) + ty_stick.*(abs(ty_stick)<abs(ty_slip));
    
    %just the traction derivative.
    dtyduy = (kty - kty.*duxyn0duxyn(2,2)).*(abs(ty_stick)<abs(ty_slip));

    %derivatives w.r.t. normal coordinates.
    dtydun = dktydun.*(uxyn(2) - uxyn0(2)).*(abs(ty_stick)<abs(ty_slip)) ...
        + dadun.*dty0da.*(abs(ty_stick)<abs(ty_slip)) ...
        - kty.*duxyn0duxyn(2, 3).*(abs(ty_stick)<abs(ty_slip)) ...
        + dtsdun.*(sign(ty_stick) + (sign(ty_stick) == 0)).*(abs(ty_stick)>=abs(ty_slip));
        
    % %%integrate forces and derivatives of forces
    fxyn(1,2) = 2*pi*a*b*(ty0.*rq)*wq';
    dfxynduxyn(2,2) = 2*pi*a*b*(dtyduy.*rq)*wq';

    dfxynduxyn(2,3) = 2*pi*a*b*(dtydun.*rq)*wq' ... 
                        + 4*pi*pars.aperb*b*(ty0.*rq)*wq'*dbdun; %Combine leading to b^2*pars.aperb + derivative


    %% Update outputs
    
    rq0 = rq*a;
    
end