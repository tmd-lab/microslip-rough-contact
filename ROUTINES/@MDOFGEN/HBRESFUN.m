function [R, dRdU, dRdw, FNL] = HBRESFUN(m, Uw, Fl, h, Nt, tol, varargin)
%HBRESFUN 
%
%   USAGE: 
%       [R, dRdU, dRdw, FNL] = HBRESFUN(m, Uw, Fl, h, Nt, tol)
%   INPUTS:
%       MDOFGEN class
%       Uw  
%       Fl
%       h
%       Nt
%       tol 
%   OUTPUTS:
%       R
%       dRdU
%       dRdw
%       FNL

  Nhc = sum((h==0)+2*(h~=0));
  
  w = Uw(end);
  
  [E, dEdw] = HARMONICSTIFFNESS(m.M, m.C, m.K, w, h);
  
  t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
  cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);  
  sct = w*TIMESERIES_DERIV(Nt, h, eye(Nhc), 1);
  
  FNL = zeros(m.Ndofs*Nhc, 1);
  dFNL = zeros(m.Ndofs*Nhc);
  for ni=1:length(m.NLTs)
    Unl = (m.NLTs(ni).L*reshape(Uw(1:end-1), m.Ndofs, Nhc))';  % Nhc x Ndnl
    Ndnl = size(m.NLTs(ni).L, 1);
    
    unlt = TIMESERIES_DERIV(Nt, h, Unl, 0);  % Nt x Ndnl
    unldot = w*TIMESERIES_DERIV(Nt, h, Unl, 1);  % Nt x Ndnl
    
    if mod(m.NLTs(ni).type, 2)==0  % Instantaneous force
      [ft, dfdu, dfdud] = m.NLTs(ni).func(t, unlt, unldot);
			% (Nt,Ndnl); (Nt,Ndnl); (Nt,Ndnl) (point-wise)
      F = GETFOURIERCOEFF(h, ft);
      J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
      dFdU = reshape(GETFOURIERCOEFF(h, reshape(dfdu.*permute(cst, [1, 3, 2]), Nt, Ndnl*Nhc)),...
                     Nhc, Ndnl, Nhc);
      for di=1:Ndnl
        J(di:Ndnl:end, di:Ndnl:end) = dFdU(:, di, :);
      end
    else  % Hysteretic force
      ft = zeros(Nt, Ndnl);
      dfdu = zeros(Nt, Ndnl, Ndnl, Nhc);
      
      its = 0;
      while its==0 || max(abs(fprev-ft(end, :)))>tol
        fprev = ft(end, :);
        for ti=1:Nt
            tm1 = mod(ti-2, Nt)+1;

            [ft(ti,:), dfdu(ti,:,:,:)] = ...
                m.NLTs(ni).func(t(ti), unlt(ti,:), h, t(tm1), ...
                unlt(tm1,:), ft(tm1,:), dfdu(tm1,:,:,:));
        end
        its = its+1;
      end
      F = GETFOURIERCOEFF(h, ft);
      J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
      for di=1:Ndnl
          for dj=1:Ndnl
              tmp = squeeze(dfdu(:, di, dj, :));
              if ~isempty(find(tmp~=0, 1))
                  J(di:Ndnl:end, dj:Ndnl:end) = ...
                      GETFOURIERCOEFF(h, tmp);
              end
          end
      end
    end
    
    if m.NLTs(ni).type<=5  % Self adjoint forcing
      FNL = FNL + reshape(m.NLTs(ni).L'*F', Nhc*m.Ndofs, 1);
      dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).L')*J*kron(eye(Nhc), m.NLTs(ni).L);
    else  % Non-self adjoint forcing
      FNL = FNL + reshape(m.NLTs(ni).Lf*F', Nhc*m.Ndofs, 1);
      dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).Lf)*J*kron(eye(Nhc), m.NLTs(ni).L);                
    end
  end
  
  % Residue
  if length(m.Rsc)~=length(Fl)
      m.Rsc = (1/max(abs(Fl)))*ones(length(Fl),1);
  end
  R = [E*Uw(1:end-1) + FNL - Fl].*m.Rsc;
  dRdU = (E+dFNL).*m.Rsc;
  dRdw = (dEdw*Uw(1:end-1)).*m.Rsc;
  
  % All Gradient terms in one matrix
  if ~isempty(varargin)
      dRdU = [dRdU dRdw];
  end
end
