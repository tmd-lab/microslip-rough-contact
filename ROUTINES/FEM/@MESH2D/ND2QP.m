function [Q,T] = ND2QP(MESH)
%ZTE_ND2QP Returns Matrix transforming nodal degrees of freedom to
%their INTERPOLATED correspondents at the required number of GLQ
%points and the matrix transforming quadrature point estimates to
%their INTEGRATED correspondents at the nodal locations. Matrix Q
%can also be used for super-convergent pressure exraction.
% USAGE:
%	[Q,T] = ZTE_ND2QP(MESH,No);
% INPUTS: 
%   MESH            : Structure containing interface mesh information
%       Nn      : Number of interface nodes
%       Nds     : (Nnx2) [x y] coordinates of interface nodes
%       Ne      : Number of interface elements
%       Ne_Tri	: Number of interface triangular elements
%       Ne_Quad	: Number of interface quadrilateral elements
%       Tri     : (Ne_Trix4) [eid n1 n2 n3]
%       Quad	: (Ne_Quadx4) [eid n1 n2 n3 n4]
%   No 		   : Order of GLQ in each direction
% OUTPUTS:
%   Q		   : (Ne(No^2)xNn) Transformation matrix from nodal
%   			to quadrature points
%   T		   : (NnxNe(No^2)) Transformation matrix from
%   			quadrature points to integrated nodal

  No = MESH.Nq;
  Nq = No^2;
  
  Q = sparse(MESH.Ne*Nq,MESH.Nn);
  T = sparse(MESH.Nn,MESH.Ne*Nq);

					       % Triangular Elements
  [X,Y,Wx,Wy]   = TRIQUAD(No,[0 0;1 0;0 1]); % Quadrature Locations
  X   	= reshape(X,No^2,1);
  Y   	= reshape(Y,No^2,1);
  Ws 		= reshape(Wx*Wy',Nq,1);
  Nv		= TRI2D_SF([X Y]);
  for e=1:MESH.Ne_Tri
    nds = MESH.Tri(e,2:4);
    V   = MESH.Nds(nds,:);
    
				% Nodal to QP
    Q((e-1)*Nq+(1:Nq),nds) = Nv;
    
				      % QP to integrated Nodal
    Jd   = TRI2D_JACDET(V,[X Y]); % Jacobian
    T(nds,(e-1)*Nq+(1:Nq)) = T(nds,(e-1)*Nq+(1:Nq)) + ...
			     repmat(Jd',3,1).*Q((e-1)*Nq+(1:Nq),nds)'.*repmat(Ws',3,1);
  end

				% Quadrilateral Elements
  [X,Y,W,~]= QUADQUAD(No,[-1 -1;1 -1;1 1;-1 1]);
  X = reshape(X,No^2, 1);
  Y = reshape(Y,No^2, 1);
  Ws= reshape(W*W', No^2, 1);
  Nv= QUAD2D_SF([X Y]);
  for e=1:MESH.Ne_Quad
    nds  	= MESH.Quad(e,2:5);
    V    	= MESH.Nds(nds,:);
    
				% Nodal to QP
    Q((MESH.Ne_Tri+e-1)*Nq+(1:Nq),nds) = Nv;

				% QP to integrated nodal
    Jd	= QUAD2D_JACDET(V,[X Y]);
    T(nds,(MESH.Ne_Tri+e-1)*Nq+(1:Nq)) = T(nds,(MESH.Ne_Tri+e-1)*Nq+(1:Nq)) + ...
					 repmat(Jd',4,1).*Q((MESH.Ne_Tri+e-1)*Nq+(1:Nq),nds)'.*repmat(Ws', 4,1);
  end
end
