function [fxyn, uxy_slip] = ELLIPSOID_MONOTONIC(uxyn, pars)
%This function returns the Ellipsoid on Ellipsoid force for a given
%displacement, only does monotonic loading from zero tangential
%displacement to inputted displacement at the fixed normal displacement
%uxyn(3)
%
%Input Parameters
%   uxyn    - Displacements with rows for different DOFs, columns for xyn
%   pars    - Parameters of the model 
%
%Output Parameters
%   fxyn        - force for each dof in directions xyn
%   uxy_slip    - displacement required to generate full slip
%
% NOTES:
%   Assumes that the short axis is aligned with x.
    
    fxyn = zeros(size(uxyn));

    %% Normal Contact, copied from ELLIPSOID PRE
    
    fxyn(:, 3) = 4*pars.Estar*sqrt(pars.Re)/3/(pars.F2)^(3/2).*uxyn(:, 3).^(3/2).*(uxyn(:, 3) > 0);
    
    %Intermediate parameter defined in Johnson
    c = (3.*fxyn(:, 3).*pars.Re./4./pars.Estar).^(1/3).*pars.F1;

    %Short axis length
    b = c/sqrt(pars.aperb); %b

    %Long axis length
    a = b*pars.aperb; %a

%     %Maximum pressure
%     p0 = 3*fxyn(3)/2/pi./a./b; %p0
    
    %% Tangential contact - X (Short Axis) - Parameter
    
    Phix = 4/(pi*(2-pars.nu))*( (1 - pars.nu/pars.e^2)*pars.Ke + pars.nu*pars.Ee/pars.e^2);
    
    %% Tangential contact - Y (Long Axis) - Parameter
    
    Phiy = 4/(pi*(2-pars.nu))*( (1 - pars.nu + pars.nu/pars.e^2)*pars.Ke - pars.nu*pars.Ee/pars.e^2);
    
    if(pars.e == 0)
        Phix = 1;
        Phiy = 1;
    end
    
    %% Calculate the Force vector
    
    %divide by the Phi's
    delta_Phi = uxyn(:, 1:2) ./ (ones(length(uxyn(:, 1)), 1)*[Phix, Phiy]);
    
    fxyn(:, 1:2) = pars.mu.*(fxyn(:, 3)*ones(1, 2)).*( 1 - (1 - 16.*pars.Gstar.*a./(3.*pars.mu.*(fxyn(:, 3).*ones(1, 2))).*delta_Phi).^(3/2));
    
    
    %%
    uxy_slip = 3*pars.mu.*fxyn(:, 3)*[Phix, Phiy]./(16*pars.Gstar*a);
    
    fxyn(:, 1:2) = fxyn(:, 1:2).*(uxyn(:, 1:2) < uxy_slip) + pars.mu*(fxyn(:, 3)*ones(1, 2)).*(uxyn(:, 1:2) >= uxy_slip);
    
end