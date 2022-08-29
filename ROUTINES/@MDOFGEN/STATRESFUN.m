function [R, dRdU, m] = STATRESFUN(m, U, Fstat)

%     FNL = zeros(m.Ndofs,1);
%     dFNL = zeros(m.Ndofs);
%     for ni=1:length(m.NLTs)
%         if mod(m.NLTs(ni).type,2)==0
%             [Fnl, dFnl, ~] = m.NLTs(ni).func(0,  m.NLTs(ni).L*U);  % Additional arguments ignored, implying zeros
%         else
%             [Fnl, dFnl] = m.NLTs(ni).func(0,  m.NLTs(ni).L*U);  % Additional arguments ignored, implying zeros
%         end
%         
%         if m.NLTs(ni).type<=5  % Self-adjoint forcing
%             FNL = FNL + m.NLTs(ni).L'*Fnl;
%             dFNL = dFNL + m.NLTs(ni).L'*dFnl*m.NLTs(ni).L;
%         else
%             FNL = FNL + m.NLTs(ni).Lf*Fnl;
%             dFNL = dFNL + m.NLTs(ni).Lf*dFnl*m.NLTs(ni).L;
%         end
%     end

%     [FNL, dFNL, ~, m] = m.NLFORCE(0, U, U*0, -1, 1);
    [FNL, dFNL, ~, m] = m.NLFORCE(0, U, U*0, 0, 1);
    
    % Residue
    R = m.K*U+FNL-Fstat;
    dRdU = m.K+dFNL;
end