function [N] = TRI2D_SF(P)
%TRI2D_SF Returns the shape functions of a 2D Tri element evaluated
%at query points
% USAGE:
%	[N] = TRI2D_SF(P);
% INPUTS:
%   P 		: Npx2 query points (natural CS)
% OUTPUTS:
%   N		: Npx3 shape function evaluates
    
    Np = size(P,1);
    N = zeros(Np,3);
    
    N(:,1) = (1.0-P(:,1)-P(:,2));
    N(:,2) = P(:,1);
    N(:,3) = P(:,2);
end