function [tx,ty,tn,dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun,varargout] ...
        = RC_FULL_TR(uxyn, pn0, pars, ASP_FUN, PZFUN, ELEM_TRAC, Nqp_heights, Nqp_radius, zmin, zmax, area_density, prev, zte_gaps) 
% RC_FULL_TR - 3D integral rough contact (integrate over probability and
% Mindlin Iwan model
% USAGE:
%	[Ptx,Pty,Pn,Jttx,Jtnx,Jtty,Jtny,Jnn] = RC_FULL_TR(uxyn,pn0,pars,uxyntxynp);
% INPUTS:
%   uxyn	: (Npx3) 2 tangential and 1 normal displacement vectors 
%   pn0         : Ignored (support for legacy formatting)
%   pars        : structure with the parameters
%   ASP_FUN - Function of forces for all directions
%   PZFUN - Probability Density function for asperity heights (min height
%           should be 0 ,max equal to zmax.)
%   ELEM_TRAC - element force function @RC_TRACTION
%   Nqp_heights - number of quadrature points used for asperity heights,
%                   uses a trapezoidal quadrature weighting scheme.
%   Nqp_radius - number of quadrature points to use at each radius for the
%                Mindlin Iwan model.
%   zmin, zmax - min and max asperity heights
%   area_density - number of asperities per square meter. 
%   prev        : Prevous state information - cell array with information
%   zte_gaps - nominal gap in ZTE for each element. in order output from Qm*stuff
% OUTPUTS:
%   Ptx         : Npx1 x tangential traction
%   Pty         : Npx1 y tangential traction
%   Pn          : Npx1 normal pressure
%   dtxdux      : Npx1 dtx/dux
%   dtxduy      : Npx1 dtx/duy
%   dtxdun      : Npx1 dtx/dun    
%   dtydux      : Npx1 dty/dux
%   dtyduy      : Npx1 dty/duy
%   dtydun      : Npx1 dty/dun        
%   dtndun      : Npx1 dtn/dun
%   (Parameter Derivatives) : Not implemented, just support for previous
%                               arrangement of outputs. 
%   varargout{4} = prev : Updated history structure
%
% NOTES:
%   1. Parameter derivatives are not correct or supported since they were
%   useless for large models with cycled RQNMA in a previous study


    %% Initialization
    Np = size(uxyn,1);
    
    tx    	= zeros(Np,1);
    ty    	= zeros(Np,1);
    tn     	= zeros(Np,1);
    
    dtxdux 	= zeros(Np,1);	
    dtxduy 	= zeros(Np,1);
    dtxdun 	= zeros(Np,1);
    
    dtydux 	= zeros(Np,1);
    dtyduy 	= zeros(Np,1);	
    dtydun 	= zeros(Np,1);
    
    dtndun 	= zeros(Np,1);
	
    %% Loop over quadrature points and calculate the tractions
    
    parfor jj = 1:Np
%     for jj = 1:Np % For debugging

        uxyn_elem = uxyn(jj, :);
        uxyn_elem(3) = uxyn_elem(3) - zte_gaps(jj);
    
        [txyn_elem, dtxynduxyn_elem, prev_elem] = ELEM_TRAC(pars, uxyn_elem, prev{jj}, ASP_FUN, PZFUN, Nqp_heights, Nqp_radius, zmin, zmax, area_density);
        
        tx(jj) = txyn_elem(1, 1);
        ty(jj) = txyn_elem(1, 2);
        tn(jj) = txyn_elem(1, 3);
        
        dtxdux(jj) 	= dtxynduxyn_elem(1, 1);
        dtxduy(jj) 	= dtxynduxyn_elem(1, 2);
        dtxdun(jj) 	= dtxynduxyn_elem(1, 3);

        dtydux(jj) 	= dtxynduxyn_elem(2, 1);
        dtyduy(jj) 	= dtxynduxyn_elem(2, 2);
        dtydun(jj) 	= dtxynduxyn_elem(2, 3);

        dtndun(jj) 	= dtxynduxyn_elem(3, 3);
    
        prev{jj} = prev_elem;
    end
    
    %% Outputs for multiple varargs
    if nargout>10
        dtxdp = zeros(Np, 1);
        dtydp = zeros(Np, 1);
        dtndp = zeros(Np, 1);
        
        %% Putting it all together
        if size(pars,1)==1  % Same parameters for all contact elements
            varargout{1} = dtxdp;
            varargout{2} = dtydp;
            varargout{3} = dtndp;
        elseif size(pars,1)==Np  % Unique parameters for each element on interface
            varargout{1} = zeros(Np, 4*Np);
            varargout{2} = zeros(Np, 4*Np);
            varargout{3} = zeros(Np, 4*Np);
            for i=1:4
                varargout{1}(:,i:4:end) = diag(dtxdp(:,i));
                varargout{2}(:,i:4:end) = diag(dtydp(:,i));
            end
            varargout{3}(:,4:4:end) = diag(dtndp(:,4));
        end
    end
    
    %% Return updated history
    
    varargout{4} = prev;
end