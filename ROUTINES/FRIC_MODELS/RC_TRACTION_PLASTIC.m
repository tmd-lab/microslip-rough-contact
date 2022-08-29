function [txyn, dtxynduxyn, prev] = RC_TRACTION_PLASTIC(pars, uxyn, prev, ASP_FUN, PZFUN, Nqp_heights, Nqp_radius, zmin, zmax, area_density)
% function determines an integrated traction over asperity heights for a
% contact function. 
%
% Inputs:
%   pars - parameters of Mindlin contact Gstar/mu/ect. 
%   uxyn - displacement [ux, uy, un] - single row. - normal displacement is
%           0 at initiation of contact. Positive into surface.
%   prev - previous state information. 
%   ASP_FUN - Function of forces in a direction (x or y) for an asperity
%   PZFUN - Probability Density function for asperity gaps. Minimimum gap
%           should be zero. 
%   Nqp_heights - number of quadrature points used for asperity heights,
%                   uses a trapezoidal quadrature weighting scheme.
%   Nqp_radius - number of quadrature points to use at each radius for the
%                Mindlin Iwan model.
%   zmin, zmax - min and max asperity gap distances
%   area_density - number of asperities per square meter. 
%
% Outputs:
%   txyn - tractions [tx, ty, tn] - single row as well. Integrated over
%       asperity heights, tractions using area density of asperities
%   dtxynduxyn - derivatives of tractions w.r.t. to displacements. Rows are
%        dtx/,dty/,dtn/. Columns are /dux,/duy,/dun.
%   
%   prev - updated states, same structure as input.


%     % Initialize Outputs:
%     txyn = zeros(1, 3);
%     dtxynduxyn = zeros(3, 3);

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
    
    uxyw = uxyn - [zeros(length(zq), 2), zq'];
    
    %% Collect data from each asperity 
    
    fxyn_all = zeros(3, length(zq));
    dfxynduxyn_all = zeros(9, length(zq));
    
    for ii = 1:Nqp_heights

        if(uxyw(ii, 3) > 0)
            
            [fxyn, dfxynduxyn, prev.rq0(ii, :), prev.tx0(ii, :), prev.ty0(ii, :), prev.deltam(ii), prev.Fm(ii), prev.am(ii)] ...
                    = ASP_FUN(pars, uxyw(ii, :), prev.uxyw0(ii, :), prev.rq0(ii, :), ...
                                prev.tx0(ii, :), prev.ty0(ii, :), rq, wq, prev.deltam(ii), prev.Fm(ii), prev.am(ii));
            
            fxyn_all(:, ii) = fxyn(:);
            dfxynduxyn_all(:, ii) = dfxynduxyn(:);
        
        else
            %No contact, so update prev state appropriately.
            prev.rq0(ii, :) = rq;
            prev.tx0(ii, :) = zeros(size(rq));
            prev.ty0(ii, :) = zeros(size(rq));
        end

    end
    
    prev.uxyw0 = uxyw;
    
    %% Integrate Forces into tractions
    
    txyn = ( area_density*(pz.*fxyn_all)*wq_z' )';
    
    dtxynduxyn = area_density*(pz.*dfxynduxyn_all)*wq_z';
    
    %Reformat into a matrix
    dtxynduxyn = reshape(dtxynduxyn(:), 3, 3);
    
end