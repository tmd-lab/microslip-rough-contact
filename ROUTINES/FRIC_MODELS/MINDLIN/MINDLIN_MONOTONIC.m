function [fxyn, uxy_slip] = MINDLIN_MONOTONIC(uxyn, pars)
%This function returns the Mindlin sphere on sphere force for a given
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
    
    
    %% Normal contact
    
    
    fxyn(:, 3) = 4.*pars.Estar.*sqrt(pars.R)/3.*uxyn(:, 3).^(3/2);

    %Radius of contact
    a = sqrt(uxyn(:, 3).*pars.R); 
    
    %displacement that causes macroslip
    uxy_slip = 3.*pars.mu.*fxyn(:, 3)./16./pars.Gstar./a;
    
    uxy = abs(uxyn(:, 1:2));
    suxy = sign(uxyn(:, 1:2));
    
    %% Tangental Contact X/Y
    
    %partial slip solution
    fxyn(:, 1:2) = pars.mu.*fxyn(:, 3).*(1 - (1 - 16.*pars.Gstar.*a.*abs(uxy)/3/pars.mu./fxyn(:, 3) ).^1.5);
    
    fxyn(:, 1:2) = fxyn(:, 1:2).*suxy.*(uxy < uxy_slip);
    
    %% Add fully slipped.
    
    fxyn(:, 1:2) = fxyn(:, 1:2) + suxy.*(fxyn(:, 3)*ones(1,2)).*pars.mu.*(uxy >= uxy_slip);
    
end