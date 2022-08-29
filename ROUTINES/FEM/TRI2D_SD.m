function [dN] = TRI2D_SD(P)
%TRI2D_SD Returns the first derivatives of the shape functions at
%provided points
% USAGE:
%	[dN] = TRI2D_SD(P);
% INPUTS:
%   P 		: Npx2 query points (natural CS)
% OUTPUTS:
%   dN		: 2Npx3 derivatives in the format,
%	  [N1,x N2,x N3,x; _xi derivatives at point p1_
%	   N1,y N2,y N3,y; _eta derivatives at point p1_
%	   N1,x N2,x N3,x; _xi derivatives at point p2_
%	   N1,y N2,y N3,y; _eta derivatives at point p2_
%          ...];

    Np = size(P,1);
    dN = kron(ones(Np,1),[-1 1 0;-1 0 1]);
end