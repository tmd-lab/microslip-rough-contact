function [N] = QUAD2D_SF(P)
%QUAD2D_SF Returns the shape functions of a 2D Quad element
%evaluated at query points
% USAGE:
%	[N] = QUAD2D_SF(P);
% INPUTS:
%   P		: Npx2 query points (natural CS) 
% OUTPUTS:
%   N		: Npx4 shape function evaluates
    
    Np = size(P,1);
    N = zeros(Np,4);
    
    N(:,1) = (1.0-P(:,1)).*(1.0-P(:,2))/4.0;
    N(:,2) = (1.0+P(:,1)).*(1.0-P(:,2))/4.0;
    N(:,3) = (1.0+P(:,1)).*(1.0+P(:,2))/4.0;
    N(:,4) = (1.0-P(:,1)).*(1.0+P(:,2))/4.0;
end