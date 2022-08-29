function [R,dRdU,dRda,dRdX,varargout] = NLRES_TMP(contactfunc, U, K, L, Fs, Fv, prev, pars, p0, QuadMats, MESH, varargin)
%NLRES Returns the nonlinear residue and its Jacobians for a
%quasi-static simulation. The parameter by default is the
%scaling applied to the given forcing vector.
% USAGE:
%	[R,dRdU,dRda] = NLRES(contactfunc, U, K, L, Fl, prev, Cmat); 
% INPUTS:
%   contactfunc	:
%   U		: ((Ndof+1)x1) displacement vector [U; alpha]
%   K		: (NdofxNdof) stiffness matrix
%   L		: (NphxNdof) Null transformation matrix
%   Fs		: (Ndofx1) force vector (to be held constant)
%   Fv		: (Ndofx1) force vector (to be scaled by alpha)
%   prev        : Previous state structure
%   	ux      : (Nnx1) vectors of interface nodal rel. x displacements
%       uy      : (Nnx1) vectors of interface nodal rel. y displacements  
%       un      : (Nnx1) vectors of interface nodal rel. normal displacements
%   	tx	: (Nnx1) vectors of interface nodal rel. x tractions
%       ty	: (Nnx1) vectors of interface nodal rel. y tractions
%       tn      : (Nnx1) vectors of interface nodal rel. normal tractions
% OUTPUTS:
%   R		: (Ndofx1) vector of residue
%   dRdX 	: (Ndofx(Ndof+1)) matrix of combined Jacobians    
%   dRdU	: (NdofxNdof) matrix of derivatives of R w.r.t
%   			displacements
%   dRda	: (Ndofx1) vector of derivative of R w.r.t alpha
%   dRdp    : (NdofxNpars) Jacobian wrt parameters in pars
%
% NOTES:
%   1. Normal traction dependency on tangential displacements can be
%   implemented in the contact function. However, the Jacobian does not
%   currently consider these effects.
%   2. Parameter derivatives are not currently being used, so should be
%   verified before using. 

    if isempty(varargin)
        Cvec = ones(length(U)-1,1);
    else
        Cvec = varargin{1};
    end

    Uph = L*(Cvec.*U(1:end-1));
    % Non-Linear Force Evaluation
    Ndof    = length(Uph);
    Nq      = size(QuadMats.Q,1);
    uxyn_qp = QuadMats.Q*reshape(Uph(1:(MESH.Nn*MESH.dpn)), MESH.dpn, MESH.Nn)';
    
    if(~iscell(prev))
        if size(prev.uxyntxyn,1)==MESH.Nn
            prev.uxyntxyn	= QuadMats.Q*prev.uxyntxyn;
        end
    end
    
    if length(p0)==MESH.Nn
        p0 	= QuadMats.Q*p0;
    end
    if length(p0)==1
        p0 = ones(size(QuadMats.Q, 1),1)*p0;
    end
    if nargout<=4
        [tx,ty,tn,dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun] = ...
            contactfunc(uxyn_qp, p0, pars, prev);
    elseif nargout==5
        [tx,ty,tn,dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun,dtxdp,dtydp,dtndp] = ...
            contactfunc(uxyn_qp, p0, pars, prev);
    end
        
    
    % Integrating & Assembling - Forces
    Fnl = zeros(Ndof,1);
    Fnl(1:(3*MESH.Nn)) = reshape((QuadMats.T*[tx ty tn])',3*MESH.Nn,1);
    
    % Integrating & Assembling - Jacobians
    Jnl = sparse(Ndof,Ndof);
    
    Jnl(1:3:(MESH.Nn*3), 1:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtxdux,1,size(QuadMats.Q,2)).*QuadMats.Q);
    Jnl(1:3:(MESH.Nn*3), 2:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtxduy,1,size(QuadMats.Q,2)).*QuadMats.Q);
    Jnl(1:3:(MESH.Nn*3), 3:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtxdun,1,size(QuadMats.Q,2)).*QuadMats.Q);
    
    Jnl(2:3:(MESH.Nn*3), 1:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtydux,1,size(QuadMats.Q,2)).*QuadMats.Q);
    Jnl(2:3:(MESH.Nn*3), 2:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtyduy,1,size(QuadMats.Q,2)).*QuadMats.Q);
    Jnl(2:3:(MESH.Nn*3), 3:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtydun,1,size(QuadMats.Q,2)).*QuadMats.Q);
    
    Jnl(3:3:(MESH.Nn*3), 3:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtndun,1,size(QuadMats.Q,2)).*QuadMats.Q);
	if length(varargin)==2
        Fnl = varargin{2}*Fnl;
        Jnl = varargin{2}*Jnl;
    end    
    % --------------------
    R 		= K*(Cvec.*U(1:end-1)) + L'*Fnl - U(end)*Fv - Fs;
    dRdU 	= K.*Cvec + (L'*Jnl*L).*Cvec;
    dRda 	= -Fv;
    dRdX 	= [dRdU dRda];
    % --------------------
    if nargout==5
        dFnldpars = zeros(Ndof, numel(pars));
        dFnldpars(1:3:MESH.Nn*3,:) = QuadMats.T*dtxdp;
        dFnldpars(2:3:MESH.Nn*3,:) = QuadMats.T*dtydp;
        dFnldpars(3:3:MESH.Nn*3,:) = QuadMats.T*dtndp;
        
        varargout{1} = L'*dFnldpars;
    end
end