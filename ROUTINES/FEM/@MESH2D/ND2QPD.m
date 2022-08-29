function [Qx,Qy,Tx,Ty] = ND2QPD(MESH,No)
%ZTE_ND2QPD Returns Matrix transforming nodal degrees of freedom to
%their INTERPOLATED derivatives at the required number of GLQ
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
    
    Nq = No^2;
    Qx = sparse(MESH.Ne*Nq,MESH.Nn);
    Qy = sparse(MESH.Ne*Nq,MESH.Nn);
    Tx = sparse(MESH.Nn,MESH.Ne*Nq);
    Ty = sparse(MESH.Nn,MESH.Ne*Nq);
    % Triangular Elements
    [X,Y,Wx,Wy]   = TRIQUAD(No,[0 0;1 0;0 1]); % Quadrature Locations
    X   	= reshape(X,No^2,1);
    Y   	= reshape(Y,No^2,1);
    Ws 		= reshape(Wx*Wy',Nq,1);
    Nv		= TRI2D_SF([X Y]);
    Nd      = TRI2D_SD([X Y]);
    Ndx     = Nd(1:2:end,:);
    Ndy     = Nd(2:2:end,:);
    for e=1:MESH.Ne_Tri
        nds = MESH.Tri(e,2:4);
        V   = MESH.Nds(nds,:);
        [J,Ji]   = TRI2D_JACMAT(V, [X Y]);
        
        % Nodal to QP
        Qx((e-1)*Nq+(1:Nq),nds) = Ndx.*Ji(1:2:end,1)+Ndy.*Ji(2:2:end,1);
        Qy((e-1)*Nq+(1:Nq),nds) = Ndx.*Ji(1:2:end,2)+Ndy.*Ji(2:2:end,2);
        
        % QP to integrated Nodal
        Jd   = TRI2D_JACDET(V,[X Y]); % Jacobian
        Tx(nds,(e-1)*Nq+(1:Nq)) = Tx(nds,(e-1)*Nq+(1:Nq))+...
				  repmat(Jd',3,1).*(Qx((e-1)*Nq+(1:Nq),nds))'.*repmat(Ws',3,1);
        Ty(nds,(e-1)*Nq+(1:Nq)) = Ty(nds,(e-1)*Nq+(1:Nq))+...
				  repmat(Jd',3,1).*(Qy((e-1)*Nq+(1:Nq),nds))'.*repmat(Ws',3,1);
    end

    % Quadrilateral Elements
    [X,Y,W,~]   = QUADQUAD(No,[-1 -1;1 -1;1 1;-1 1]);
    X    	= reshape(X,No^2, 1);
    Y    	= reshape(Y,No^2, 1);
    Ws 		= reshape(W*W', No^2, 1);
    Nv 		= QUAD2D_SF([X Y]);
    Nd      = QUAD2D_SD([X Y]);
    Ndx     = Nd(1:2:end,:);
    Ndy     = Nd(2:2:end,:);
    for e=1:MESH.Ne_Quad
        nds  	= MESH.Quad(e,2:5);
        V    	= MESH.Nds(nds,:);
        [J,Ji]  = QUAD2D_JACMAT(V, [X Y]);
        
        % Nodal to QP
        Qx((MESH.Ne_Tri+e-1)*Nq+(1:Nq),nds) = Ndx.*Ji(1:2:end,1)+Ndy.*Ji(2:2:end,1);
        Qy((MESH.Ne_Tri+e-1)*Nq+(1:Nq),nds) = Ndx.*Ji(1:2:end,2)+Ndy.*Ji(2:2:end,2);
        
        % QP to integrated nodal
        Jd	= QUAD2D_JACDET(V,[X Y]);
        Tx(nds,(MESH.Ne_Tri+e-1)*Nq+(1:Nq)) = Tx(nds,(MESH.Ne_Tri+e-1)*Nq+(1:Nq))+ ...
					      repmat(Jd',4,1).*(Qx((MESH.Ne_Tri+e-1)*Nq+(1:Nq),nds))'.*repmat(Ws',4,1);
        Ty(nds,(MESH.Ne_Tri+e-1)*Nq+(1:Nq)) = Ty(nds,(MESH.Ne_Tri+e-1)*Nq+(1:Nq))+ ...
					      repmat(Jd',4,1).*(Qy((MESH.Ne_Tri+e-1)*Nq+(1:Nq),nds))'.*repmat(Ws',4,1);
    end
end
