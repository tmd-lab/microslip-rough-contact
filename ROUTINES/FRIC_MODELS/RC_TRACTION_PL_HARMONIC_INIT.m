function [txyn, dtxynduxynh, prev] = RC_TRACTION_PL_HARMONIC_INIT(uxyn_init, uxyn0, h, t, ASP_FUN_PRE, ...
                    PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density, quadz)
% Function returns traction time history and derivatives w.r.t. harmonic
% coefficients for a rough contact model. 
%
% Inputs:
%   uxyn_init - Displacement with first tangential and maximum normal
%   uxyn0 - Full displacements at 0 time. - Used for the previous
%           information for subsequent calls to other things. For
%           consistency, current usuage is to have uxyn0 = uxyn_init with
%           tangential components equal to the zeroth harmonic.
%   h    - Harmonic coefficients that are used in the simulation
%   t    - time instants corresponding to rows of uxyn.
%   ASP_FUN_PRE - purely normal asperity model
%   PZFUN - Probability density function
%   pars - model parameters
%   prev - previous information for history
%   Nqp_heights - number of quadrature points for asperity heights
%   Nqp_radius - number of quadrature points for the radius
%   zmin - min asperity height
%   zmax - maximum asperity height
%   area_density - area density of asperities
%   quadz - mesoscale height of the quadrature point for the current
%           element.
%
% Outputs:
%   txyn - tractions at current time instant in the same shape as uxyn
%   dtxynduxynh - harmonic derivatives of tractions:
%                   dim 1 = time instant - current
%                   dim 2 = tractions or dt
%                   dim 3 = displacements of du
%                   dim 4 = harmonics w.r.t. which du is calculated
%   prev - updated history

    
    %% Asperity Height Quadrature
    zq = linspace(zmin, zmax, Nqp_heights);
    wq_z = ones(size(zq));
    wq_z(2:end-1) = 2;
    wq_z = wq_z/sum(wq_z);

    pz = PZFUN(zq);

    %% Asperity Radius Quadrature

    rq = linspace(0, 1, Nqp_radius);
    wq = ones(size(rq));
    wq(2:end-1) = 2;
    wq = wq/sum(wq);
    
    
    %% Change normal displacement to asperity interference displacement.
        
    uxyw = (uxyn_init - [0, 0, quadz]) - [zeros(length(zq), 2), zq'];
    uxyw0 = (uxyn0 - [0, 0, quadz]) - [zeros(length(zq), 2), zq'];
    

    %% Initialize Memory
    
    Nhc = sum((h==0)+2*(h~=0));
    
    txyn = zeros(size(uxyn_init));
    dtxynduxynh = zeros(length(t), size(txyn, 1), size(uxyn_init,1), Nhc);
    
    prev.dFmddeltam = zeros(Nqp_heights, 1);
    prev.damddeltam = zeros(Nqp_heights, 1);
    
    prev.dfxyn0duxynh = zeros(3, 3, Nhc, Nqp_heights);
    
    %% Collect data from each asperity 
    
    for ii = 1:Nqp_heights

        if(uxyw(ii, 3) > 0)
            
            
            [fxyn, dfxynduxyn, prev.rq0(ii, :), prev.tx0(ii, :), prev.ty0(ii, :), ...
                prev.deltam(ii), prev.Fm(ii), prev.am(ii), ...
                dfduxyn0, dfdfxy0, dfddeltam, a, dadun, Sys, dSysdun, daddeltam, dSysddeltam] ...
                    = ASP_FUN_PRE(pars, uxyw(ii, :), prev.uxyw0(ii, :), prev.rq0(ii, :), ...
                                prev.tx0(ii, :), prev.ty0(ii, :), rq, wq, prev.deltam(ii), prev.Fm(ii), prev.am(ii));
            
            txyn = txyn + fxyn*pz(ii)*wq_z(ii);
            
            prev.dFmddeltam(ii, 1) = dfxynduxyn(3,3);
            prev.damddeltam(ii, 1) = dadun;
            
        else
            %No contact, so update prev state appropriately.
            prev.rq0(ii, :) = rq;
            prev.tx0(ii, :) = zeros(size(rq));
            prev.ty0(ii, :) = zeros(size(rq));
        end

    end
    
    prev.uxyw0 = uxyw0;
        
    %% Integrate Forces into tractions
    
    txyn = area_density*txyn;    

end