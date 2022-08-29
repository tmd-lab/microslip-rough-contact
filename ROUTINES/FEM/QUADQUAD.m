function [X,Y,W,Jd] = QUADQUAD(No,V)
%QUADQUAD Returns GL quadrature points and weights for a 2D quadrilateral
% USAGE:
%	[X,Y,W,Jd] = QUADQUAD(No,V);
% INPUTS:
%   No		: Number of points in each dimension 
%   V		: 4x2 vector of vertex node coordinates
% OUTPUTS:
%   X		: NxN matrix of x values 
%   Y		: NxN matrix of y values 
%   W    	: Nx1 weight vector
%   Jd      : NxN matrix of Jacobian determinants
% EXAMPLE USAGE:
% 	func = @(X) X(:,1)+X(:,2);
%	No = 10; % Order=10    
%	[X,Y,W,Jd] = QUADQUAD(No,V);
%	X = reshape(X,No^2,1);
%	Y = reshape(Y,No^2,1);
%	F = reshape(func([X Y]),No,No);
%	Integral = W'*(F.*Jd)*W;
    
    [x,W] = LGWT(No,-1,1);
    [X,Y] = meshgrid(x,x);

    Xr = reshape(X,No^2,1);
    Yr = reshape(Y,No^2,1);
    
    Jr = QUAD2D_JACMAT(V,[Xr Yr]);
    Jd = zeros(No^2,1);
    for k=1:(No^2)
        Jd(k) = det(Jr((k-1)*2+(1:2),:));
    end
    Jd = reshape(Jd, No,No);
end