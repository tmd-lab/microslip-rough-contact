function [Jd] = QUAD2D_JACDET(V,P)
%QUAD2D_JACDET Returns the Jacobian determinant evaluated at query
%points for a 2D quad element defined by given ordered node
%locations
% USAGE:
%	[Jd] = QUAD2D_JACDET(V,P);
% INPUTS:
%   V		: 4x2 order coordinates of vertices (global CS)
%   P		: Npx2 query points (natural CS)
% OUTPUTS:
%   Jd		: Npx1 stacked Jacobian determinants

    Np = size(P,1);
    J = zeros(2*Np,2);
    dN = QUAD2D_SD(P);
    
    J(1:2:end,1:2) = [dN(1:2:end,:)*V(:,1) dN(2:2:end,:)*V(:,1)];
    J(2:2:end,1:2) = [dN(1:2:end,:)*V(:,2) dN(2:2:end,:)*V(:,2)];
    Jd		   = zeros(Np,1);
    for i=1:Np
        Jd(i)	   = det(J((i-1)*2+(1:2),(1:2)));
    end
end

%% Can be combined with TRI2D_JACMAT %%