function [R, dRdUl, dRdq] = RQMRESFUN(m, Ulq, lsc, varargin)
%RQMRESFUN returns the residue function for RQNM
%
%   In the presence of static states (Fstat, Ustat), the following residual
%   is returned:
%       [Ku+fnl - \lambda M(u-Ustat) - Fstat;
%        (u-Ustat)^T(Ku+fnl - Fstat) - \lambda q^2]
%  USAGE:
%   [R, dRdUl, dRdq] = m.RQMRESFUN(Ulq, lsc);
%       (or)
%   [R, dRdUl, dRdq] = m.RQMRESFUN(Ulq, lsc, Fstat, Ustat);
%  INPUTS:
%   Ulq     : (Nd+2 x 1)
%   lsc     : 1 or 0 - Log scale amplitude or not
%   Fstat   : (Nd x 1) Static force
%   Ustat   : (Nd x 1) Static Solution
%  OUTPUTS:
%   R       : (Nd+1 x 1)
%   dRdUl   : (Nd+1 x Nd+1)
%   dRdq    : (Nd+1 x 1)
        
    if lsc==1
        Q = 10^Ulq(end);
        dQdlq = Q*log(10);
    else
        Q = Ulq(end);
        dQdlq = 1.0;
    end
    
    FNL = zeros(m.Ndofs,1);
    dFNL = zeros(m.Ndofs);
    for ni=1:length(m.NLTs)
        if mod(m.NLTs(ni).type,2)==0
            [Fnl, dFnl, ~] = m.NLTs(ni).func(0,  m.NLTs(ni).L*Ulq(1:end-2));  % Additional arguments ignored, implying zeros
        else
            [Fnl, dFnl] = m.NLTs(ni).func(0,  m.NLTs(ni).L*Ulq(1:end-2));  % Additional arguments ignored, implying zeros
        end
        
        if m.NLTs(ni).type<=5  % Self-adjoint forcing
            FNL = FNL + m.NLTs(ni).L'*Fnl;
            dFNL = dFNL + m.NLTs(ni).L'*dFnl*m.NLTs(ni).L;
        else
            FNL = FNL + m.NLTs(ni).Lf*Fnl;
            dFNL = dFNL + m.NLTs(ni).Lf*dFnl*m.NLTs(ni).L;
        end
    end
    
%     % Residue - "Vanilla" version
%     R = [m.K*Ulq(1:end-2)+FNL-Ulq(end-1)*m.M*Ulq(1:end-2);
%         0.5*Ulq(1:end-2)'*m.M*Ulq(1:end-2)-0.5*Q^2];
%     dRdUl = [m.K+dFNL-Ulq(end-1)*m.M, -m.M*Ulq(1:end-2);
%         Ulq(1:end-2)'*m.M, 0];
%     dRdq = [zeros(m.Ndofs,1);-Q*dQdlq]; 
    
    % Residue - Better conditioned version (Jacobian not nearly singular at
    % solution)    
    if length(varargin)==2
        Ustat = varargin{2};
        Fstat = varargin{1};
        
        R = [m.K*Ulq(1:end-2)+FNL-Ulq(end-1)*m.M*(Ulq(1:end-2)-Ustat) - Fstat;
            (Ulq(1:end-2)-Ustat)'*(m.K*Ulq(1:end-2)+FNL-Fstat)-Ulq(end-1)*Q^2];

        dRdUl = [m.K+dFNL-Ulq(end-1)*m.M, -m.M*Ulq(1:end-2);
            (Ulq(1:end-2)-Ustat)'*(K+dFNL)+Ulq(1:end-2)'*K+FNL'-Fstat', -Q^2];

        dRdq = [zeros(m.Ndofs,1);-2*Ulq(end-1)*Q*dQdlq]; 
    else
        R = [m.K*Ulq(1:end-2)+FNL-Ulq(end-1)*m.M*Ulq(1:end-2);
            Ulq(1:end-2)'*(m.K*Ulq(1:end-2)+FNL)-Ulq(end-1)*Q^2];
    
        dRdUl = [m.K+dFNL-Ulq(end-1)*m.M, -m.M*Ulq(1:end-2);
            Ulq(1:end-2)'*(2*m.K+dFNL)+FNL', -Q^2];
        dRdq = [zeros(m.Ndofs,1);-2*Ulq(end-1)*Q*dQdlq]; 
    end
end