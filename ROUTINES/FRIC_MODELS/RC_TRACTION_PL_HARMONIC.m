function [txyn, dtxynduxynh, prev] = RC_TRACTION_PL_HARMONIC(uxyn, h, t, ASP_FUN, ...
    PZFUN, pars, prev, Nqp_heights, Nqp_radius, zmin, zmax, area_density, quadz, cst)
% Function returns traction and derivative w.r.t. harmonic coefficients at 
% the current instant 
%
% Inputs:
%   uxyn - history of displacements to calculate forces on. Each row is a
%           different time instant.
%   h    - Harmonic coefficients that are used in the simulation
%   t    - time instants corresponding to rows of uxyn.
%   ASP_FUN - asperity model
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
%   cst - cosine, sine values at the current time for the appropriate
%           harmonics
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
    
    uxyw = (uxyn - [0, 0, quadz]) - [zeros(length(zq), 2), zq'] ;
    
    %% Initialize Memory
    
    
    Nhc = sum((h==0)+2*(h~=0));
    
    txyn = zeros(size(uxyn));
    dtxynduxynh = zeros(length(t), length(txyn), length(uxyn), Nhc);
    
        
    %% Collect data from each asperity 
    
    duxynduxynh = reshape(cst, 1, 1, []).*eye(3);
    
    for ii = 1:Nqp_heights

        if(uxyw(ii, 3) > 0)
            
            [fxyn, dfxynduxyn, prev.rq0(ii, :), prev.tx0(ii, :), prev.ty0(ii, :), ...
                prev.deltam(ii), prev.Fm(ii), prev.am(ii), dfduxyn0, dfdfxy0, dfddeltam] ...
                = ASP_FUN(pars, uxyw(ii, :), prev.uxyw0(ii, :), prev.rq0(ii, :), ...
                            prev.tx0(ii, :), prev.ty0(ii, :), rq, wq, prev.deltam(ii), ...
                            prev.Fm(ii), prev.am(ii), prev.dFmddeltam(ii, 1), prev.damddeltam(ii, 1));
            
            txyn = txyn + fxyn*pz(ii)*wq_z(ii);
            
            dfxynduxynh_curr = zeros(3,3,Nhc);
            
            for jj = 1:Nhc
                dfxynduxynh_curr(:, :, jj) = dfxynduxyn*duxynduxynh(:, :, jj) ...
                                            + dfduxyn0*prev.duxyn0duxynh(:, :, jj) ... 
                                            + (dfdfxy0*prev.dfxyn0duxynh(1:2, :, jj, ii)) ...
                                            + [zeros(3,2), dfddeltam].*prev.ddeltamduxynh(:, :, jj);
                
            end
            
            dtxynduxynh = dtxynduxynh + pz(ii)*wq_z(ii)*permute(dfxynduxynh_curr, [4, 1, 2, 3]);
                        
            prev.dfxyn0duxynh(:, :, :, ii) = dfxynduxynh_curr;
        else
            %No contact, so update prev state appropriately.
            prev.rq0(ii, :) = rq;
            prev.tx0(ii, :) = zeros(size(rq));
            prev.ty0(ii, :) = zeros(size(rq));
            
            prev.dfxyn0duxynh(:, :, :, ii) = 0.*prev.dfxyn0duxynh(:, :, :, ii);
        end

    end
    
    prev.uxyw0 = uxyw;
    prev.duxyn0duxynh = duxynduxynh;
        
    %% Integrate Forces into tractions
    
    txyn = area_density*txyn;
    dtxynduxynh = area_density*dtxynduxynh;
    


end