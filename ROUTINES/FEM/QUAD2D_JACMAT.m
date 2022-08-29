function [J] = QUAD2D_JACMAT(V,P)
%QUAD2D_JACMAT Returns the Jacobian matrix evaluated at query
%points for a 2D quad element defined by given ordered node
%locations
% USAGE:
%	[J] = QUAD2D_JACMAT(V,P);
% INPUTS:
%   V		: 4x2 order coordinates of vertices (global CS)
%   P		: Npx2 query points (natural CS)
% OUTPUTS:
%   J		: 2Npx2 stacked Jacobian matrices

    Np = size(P,1);
    J = zeros(2*Np,2);
    dN = QUAD2D_SD(P);
    
    J(1:2:end,1:2) = [dN(1:2:end,:)*V(:,1) dN(2:2:end,:)*V(:,1)];
    J(2:2:end,1:2) = [dN(1:2:end,:)*V(:,2) dN(2:2:end,:)*V(:,2)];
end

%% Can be combined with TRI2D_JACMAT %%