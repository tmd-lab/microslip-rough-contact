function [Fn, Z, dFndUn, dFndUnd, dFndPars, Fxyn_qp] = CONTACTEVAL(m, Un, Z, Und, Pars, varargin)
%CONTACTEVAL evaluates the given hysteretic contact forces at quadrature points and integrates them to provide nodal forces
%
% USAGE:
% ------
%   MESH.SETCFUN(fcont);  % Intialize
%   [Fn, Z, dFndUn, dFndUnd, dFndPars] = MESH.CONTACTEVAL(Un, Und, Z, Pars)
% INPUTS:
% -------
%   (Initialization)  
%   fcont	: Contact function handle of the form,
%		 	[fxyn, z, DfxynDuxyn, DfxynDuxynd, DfxynDpars] = fcont(uxyn, z, uxynd, Pars);
%   (Inputs)  
%   Un		: (MESH.Nn,1) Nodal displacements
%   Z		: (Nz,MESH2D.Ne*MESH2D.Nq^2)
%   Und		: (MESH.Nn,1) Nodal velocities  
%   Pars	: (Npars,1)  Parameters array
%   pA(optional): (xxx,Npars) Parameter access array. fcont will be called using
%			  reshape(pA*Pars, [], MESH.Ne*MESH.Nq^2) as Pars
%   pU(optional): (xxx,length(Un)) Projection/Mask for Un such that the Un that will be
%			  used for the calculations is mU*Un.
% OUTPUTS:
% --------
%   Fn		:
%   Z		:
%   dFndUn	:
%   dFndUnd	:
%   dFndPars	:

  if isempty(m.fcont)
    error('Contact force function not set');
  end

  if length(varargin)>=1
    pA = varargin{1};
    Parsc = reshape(pA*Pars, [], m.Ne*m.Nq^2);
  else
    pA = speye(length(Pars));
    Parsc = reshape(Pars, [], m.Ne*m.Nq^2);
  end
  Npars = size(Parsc, 1);

  if length(varargin)>=2
    pU = varargin{2};
    Un = pU*Un;
    Und = pU*Und;
  end

  Nu = length(Un);     % Number of dofs in U vector
  Ncdof = m.Nn*m.dpn;  % Number of contact DoFs
  
  Uxyn_qp = m.INTERP_QP(reshape(Un(1:Ncdof), m.dpn,m.Nn)');  % Uxyn_qp=Qm*reshape(Un, dpn, [])';
  Uxynd_qp = m.INTERP_QP(reshape(Und(1:Ncdof), m.dpn,m.Nn)');% Uxynd_qp=Qm*reshape(Und, dpn, [])';

  [Fxyn_qp, Z, DfxynDuxyn_qp, DfxynDuxynd_qp, DfxynDpars_qp] = m.fcont(Uxyn_qp', Z, Uxynd_qp', Parsc);

  Fn = reshape( m.INTEG_QP( Fxyn_qp' )', Ncdof, 1);  % Fn = reshape((Tm*Fxyn_qp')', [], 1);
  Fn = [Fn; zeros(Nu-Ncdof,1)];  % Pad rest with zeros

  dFndUn = zeros(Nu);
  dFndUnd = zeros(Nu);
  dFndPars = zeros(Nu, length(Pars(:)));

  % CAN BE PARALLELIZED SINCE ACCESS IS EXCLUSIVE
  for i=1:m.dpn
    for j=1:m.dpn
      % Tm*der*Qm
      dFndUn(i:m.dpn:Ncdof, j:m.dpn:Ncdof) = m.MATTFM_QP2ND(diag(squeeze(DfxynDuxyn_qp(i, j, :))));
      dFndUnd(i:m.dpn:Ncdof, j:m.dpn:Ncdof) = m.MATTFM_QP2ND(diag(squeeze(DfxynDuxynd_qp(i, j, :))));
    end

    for j=1:Npars
      dFndPars(i:m.dpn:Ncdof, :) = dFndPars(i:m.dpn:Ncdof, :) + ...
				   m.INTEG_QP(diag(squeeze(DfxynDpars_qp(i, j, :))))*pA(j:Npars:end, :);
    end
  end

  if length(varargin)>=2  % Project/Mask everything back
    Fn = pU'*Fn;
    dFndUn = pU'*dFndUn*pU;
    dFndUnd = pU'*dFndUnd*pU;
    dFndPars = pU'*dFndPars;
  end
end
