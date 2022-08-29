function [R,dRdUl,dRdq,dRdX,varargout] = NLRES_RQNM_TMP(contactfunc, Ulq, Us, dUsdp, M, K, L, Fs, prev, pars, p0, QuadMats, MESH, varargin)
%NLRES Returns the nonlinear residue and its Jacobians for a
%quasi-static simulation. The parameter by default is the
%scaling applied to the given forcing vector.
% USAGE:
%
% INPUTS:
%   contactfunc	: Friction model contact function
%   Ulq     : Displacements for the current level q excitation [U; lambda; log10(q)]
%   Us		: Static Displacement Vector
%   dUsdp   : Parameter derivatives of Us
%   M       : (NdofxNdof) mass matrix
%   K		: (NdofxNdof) stiffness matrix
%   L		: (NphxNdof) Null transformation matrix
%   Fs		: (Ndofx1) force vector (to be held constant)
%   prev        : Previous state structure - whatever form contactfunc
%                   takes
% OUTPUTS:
%   R		: (Ndofx1) vector of residue
%   dRdUl	: (NdofxNdof) matrix of derivatives of R w.r.t
%   			displacements and lambda
%   dRdq	: (NdofxNdof) matrix of derivatives of R w.r.t
%   			amplitude
%   dRdX	: (Ndofx1) vector of derivative of R w.r.t full unknown vector
%   varargout{1}    : (NdofxNpars) Jacobian wrt parameters in pars (not
%                     fully supported)
%
% NOTES:
% 1. Parameter Derivatives are not verified because they are currently
% unused.


    if isempty(varargin)
        Cvec = ones(length(Ulq)-2,1);
    else
        Cvec = varargin{1};
    end

    Uph = L*(Cvec.*Ulq(1:end-2));
    lam = Ulq(end-1);
    q = 10^Ulq(end);
    dqdu = q*log(10);
    % Non-Linear Force Evaluation
    Ndof    = length(Uph);
    uxyn_qp = QuadMats.Q*reshape(Uph(1:(MESH.Nn*MESH.dpn)), MESH.dpn, MESH.Nn)';
    
    if ~iscell(prev) && size(prev.uxyntxyn,1)==MESH.Nn
        prev.uxyntxyn	= QuadMats.Q*prev.uxyntxyn;
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
    % Residue
    % --------------------
    R 		= [K*(Cvec.*Ulq(1:end-2))+L'*Fnl-lam*M*(Ulq(1:end-2)-Us)-Fs;
        (Ulq(1:end-2)-Us)'*M*(Ulq(1:end-2)-Us)-q^2];
    dRdUl 	= [K.*Cvec+L'*Jnl*L.*Cvec-lam*M, -M*(Ulq(1:end-2)-Us);
        2*(Ulq(1:end-2)-Us)'*M, 0];
    dRdq 	= [zeros(size(Ulq(1:end-2))); -2*dqdu*q];
    dRdX 	= [dRdUl dRdq];
    
    if nargout==5
        dFnldpars = zeros(Ndof, numel(pars));
        dFnldpars(1:3:MESH.Nn*3,:) = QuadMats.T*dtxdp;
        dFnldpars(2:3:MESH.Nn*3,:) = QuadMats.T*dtydp;
        dFnldpars(3:3:MESH.Nn*3,:) = QuadMats.T*dtndp;
        
        varargout{1} = [L'*dFnldpars + lam*M*dUsdp; 
            -2*(Ulq(1:end-2)-Us)'*M*dUsdp];
    end
end