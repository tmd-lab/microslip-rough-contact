function [dN] = QUAD2D_SD(P)
%QUAD2D_SD Returns the first derivatives of the shape functions
%evaluated at query points
% USAGE:
%	[dN] = QUAD2D_SD(P);
% INPUTS:
%   P		: Npx2 query points (natural CS)
% OUTPUTS:
%   dN		: 2Npx4 derivatives in the format,
%	  [N1,x N2,x N3,x N4,x; _xi derivatives at point p1_
%	   N1,y N2,y N3,y N4,y; _eta derivatives at point p1_
%	   N1,x N2,x N3,x N4,x; _xi derivatives at point p2_
%	   N1,y N2,y N3,y N4,y; _eta derivatives at point p2_
%          ...];

    Np = size(P,1);
    dN = zeros(2*Np,4);
    
    dN(1:2:end,1) = -(1.0-P(:,2))/4.0;
    dN(1:2:end,2) =  (1.0-P(:,2))/4.0;
    dN(1:2:end,3) =  (1.0+P(:,2))/4.0;
    dN(1:2:end,4) = -(1.0+P(:,2))/4.0;

    dN(2:2:end,1) = -(1.0-P(:,1))/4.0;
    dN(2:2:end,2) = -(1.0+P(:,1))/4.0;
    dN(2:2:end,3) =  (1.0+P(:,1))/4.0;
    dN(2:2:end,4) =  (1.0-P(:,1))/4.0;
end