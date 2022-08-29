function [fxyn, dfxynduxyn, rq0, fx0, fy0] = ELLIPSOID_TAN_DECOUPLE_ADHESION(pars, uxyn, uxyn0, rq0, fx0, fy0, rq, wq)
% Function determines forces for the contact state of an asperity. Uses an
% elastic dry friction model with correct tangential stiffness.
% 
% Inputs:
%   pars - structure of material and asperity parameters - intialize
%           elliptical contact parameters with: pars = INITIALIZE_ECCENTRIC(pars);
%   uxyn - displacements, un is the gap distance with positive being in
%           contact
%   uxyn0 - previous displacements, un0 is again gap distance
%   rq0   - quadrature radii at last instant [0, a_prev]
%   tx0   - tractions in x direction at rq0 at last instant
%   ty0   - tractions in y direction at rq0 at last instant
%   rq    - quadrature radii to be used normalized to a, so range [0, 1]
%   wq    - weights for the quadrature radii. 
%
% Outputs:
%   fxyn  - forces in xyn directions 
%   dfxynduxyn - Jacobian stiffness matrix of fyxn 
%   rq0 - quadrature radii at the end - equals rq
%   tx0 - tractions at each radii after instant 
%   ty0 - tractions at each radii after this instant. 
%
% NOTES:
%   1. The quadrature integration scheme is ignored by this function since
%   the asperity is treated as one slider


    %% Re-establishment of contact conditions 
    % One could apply this to each radii instead of just the full asperity
    % having new contact, but that is a complicated nonlinear equation for
    % a(w). 
    
    %initialize derivatives of previous information w.r.t. new disps.
    duxyn0duxyn = zeros(3,3); %Independent unless new contact
    
    % This section should do an update as:
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
%         rq0 = rq; %Not relevant for single slider representation
        fx0(1) = 0;
        fy0(1) = 0;
        
    end
    
    %% Normal Forces + derivatives
    
    [fxyn, dfxynduxyn, ~, ~, ~, ~, ~, a, dadun] ...
        = ELLIPSOID_PRE(pars, uxyn, uxyn0, rq0, fx0, fy0, rq, wq);
    
%     fs = fxyn(3)*pars.mu;
%     dfsdun = dfxynduxyn(3,3)*pars.mu;
    
    % Full area created averaged over full distance
    fs = pi*a*pars.DeltaGamma;
    dfsdun = pi*dadun*pars.DeltaGamma;
    
%     % Linearized new area at initially overlapping contact
%     fs = pars.DeltaGamma*a;
%     dfsdun = dadun*pars.DeltaGamma;
    
    %% Tangential Stiffnesses 
    
    % In paper, Phix = Phib (aligned with short semi-axis)
    %           Phiy = Phia (aligned with long semi-axis)
    Phix = 4/(pi*(2-pars.nu))*((1 - pars.nu/pars.e^2)*pars.Ke + pars.nu*pars.Ee/pars.e^2);
    Phiy = 4/(pi*(2-pars.nu))*( (1 - pars.nu + pars.nu/pars.e^2)*pars.Ke - pars.nu*pars.Ee/pars.e^2);
    
    if(pars.e == 0)
        Phix = 1;
        Phiy = 1;
    end
    
    ktx = 4*pars.G*a/(2 - pars.nu)/Phix;
    
    dktxdun = 4*pars.G*dadun/(2 - pars.nu)/Phix;
    
    kty = 4*pars.G*a/(2 - pars.nu)/Phiy;

    dktydun = 4*pars.G*dadun/(2 - pars.nu)/Phiy;

    %% Tangential Contact - X
    
    % %% Calculation of traction
    
    %stick and slip forces
    tx_stick = fx0(1) + ktx.*(uxyn(1) - uxyn0(1));
    tx_slip  = fs.*(sign(tx_stick) + (sign(tx_stick) == 0)); %if stick is zero, just set sign to positive.

    fx0 = tx_slip.*(abs(tx_stick)>=abs(tx_slip)) + tx_stick.*(abs(tx_stick)<abs(tx_slip));
    
    %just the traction derivative.
    dtxdux = (ktx - ktx*duxyn0duxyn(1,1)).*(abs(tx_stick)<abs(tx_slip));

    %derivatives w.r.t. normal coordinates.
    dtxdun = dktxdun.*(uxyn(1) - uxyn0(1)).*(abs(tx_stick)<abs(tx_slip)) ...
        - ktx.*duxyn0duxyn(1,3).*(abs(tx_stick)<abs(tx_slip)) ...
        + dfsdun.*(sign(tx_stick) + (sign(tx_stick) == 0)).*(abs(tx_stick)>=abs(tx_slip));
    
    % %%integrate forces and derivatives of forces
    fxyn(1,1) = fx0;
    dfxynduxyn(1,1) = dtxdux;

    dfxynduxyn(1, 3) = dtxdun;

    %% Tangential Contact - Y
    
    % %% Calculation of traction
    
    %stick and slip forces
    ty_stick = fy0(1) + kty.*(uxyn(2) - uxyn0(2));
    ty_slip  = fs.*(sign(ty_stick) + (sign(ty_stick) == 0)); %if stick is zero, just set sign to positive.

    fy0 = ty_slip.*(abs(ty_stick)>=abs(ty_slip)) + ty_stick.*(abs(ty_stick)<abs(ty_slip));
    
    %just the traction derivative.
    dtyduy = (kty - kty*duxyn0duxyn(2,2)).*(abs(ty_stick)<abs(ty_slip));

    %derivatives w.r.t. normal coordinates.
    dtydun = dktydun.*(uxyn(2) - uxyn0(2)).*(abs(ty_stick)<abs(ty_slip)) ...
        - kty.*duxyn0duxyn(2,3).*(abs(ty_stick)<abs(ty_slip)) ...
        + dfsdun.*(sign(ty_stick) + (sign(ty_stick) == 0)).*(abs(ty_stick)>=abs(ty_slip));
    
    % %%integrate forces and derivatives of forces
    fxyn(1, 2) = fy0;
    dfxynduxyn(2, 2) = dtyduy;

    dfxynduxyn(2, 3) = dtydun;

end