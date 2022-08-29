function [J] = TRI2D_JACMAT(V,P)
%TRI2D_JACMAT Returns the Jacobian matrix evaluated at query points
%for a 2D Tri element defined by given ordered node locations
% USAGE:
%	[J] = TRI2D_JACMAT(V,P);
% INPUTS:
%   V		: 3x2 order coordinates of vertices (global CS)
%   P		: Npx2 query points (natural CS)
% OUTPUTS:
%   J		: 2Npx2 stacked Jacobian matrices

    Np = size(P,1);
    J = zeros(2*Np,2);
    dN = TRI2D_SD(P);
    
    J(1:2:end,1:2) = [dN(1:2:end,:)*V(:,1) dN(2:2:end,:)*V(:,1)];
    J(2:2:end,1:2) = [dN(1:2:end,:)*V(:,2) dN(2:2:end,:)*V(:,2)];
end

%% Can be combined with QUAD2D_JACMAT %%