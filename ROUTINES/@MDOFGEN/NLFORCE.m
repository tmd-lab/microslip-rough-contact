function [F, dFdU, dFdUd, m] = NLFORCE(m, t, U, Ud, tp, varargin)
%NLFORCE evaluates force and jacobians. Member of MDOFGEN class.
%     
%  USAGE:
%    [F, dFdU, dFdUd, m] = m.NLFORCE(t, U, Ud, tp, init);
%  INPUTS:
%    t		: scalar time
%    U, Ud	: Ndx1 disp & vel vectors
%    tp		: scalar previous time instant
%    init	: [optional] 0 if initializing
%    
    if length(varargin)==1
        init = varargin{1};  % 1 if this is the initalization run, 0 if not
    else
      init = 0;
    end
    
    if length(varargin)==2
        cjac = varargin{2};
    else
        cjac = 1;
    end
        

    F = zeros(m.Ndofs, 1);
    dFdU = zeros(m.Ndofs);
    dFdUd = zeros(m.Ndofs);
    
    for ni=1:length(m.NLTs)
        if mod(m.NLTs(ni).type,2)==0  % Inst. force
            [f, dfdu, dfdud] = m.NLTs(ni).func(t, m.NLTs(ni).L*U, m.NLTs(ni).L*Ud);
        else
            if ~isfield(m.NLTs(ni), 'up') || ~isfield(m.NLTs(ni), 'fp') || ~isfield(m.NLTs(ni), 'dfdup') || init 
                [f, dfdu] = m.NLTs(ni).func(t, m.NLTs(ni).L*U);
            else
                [f, dfdu] = m.NLTs(ni).func(t, m.NLTs(ni).L*U, 0, tp, m.NLTs(ni).up, ...
                    m.NLTs(ni).fp, m.NLTs(ni).dfdup);
            end

            dfdu = sparse(dfdu);
            
            m.NLTs(ni).up = m.NLTs(ni).L*U;
            m.NLTs(ni).fp = f;
            m.NLTs(ni).dfdup = full(dfdu);
        end

        if m.NLTs(ni).type<=5
            F = F + m.NLTs(ni).L'*f;
            
            if cjac
                dFdU = dFdU + m.NLTs(ni).L'*dfdu*m.NLTs(ni).L;
                if mod(m.NLTs(ni).type,2)==0
                    dFdUd = dFdUd + m.NLTs(ni).L'*dfdud*m.NLTs(ni).L;
                end
            end
        else
            F = F + m.NLTs(ni).Lf*f;
            
            if cjac
                dFdU = dFdU + m.NLTs(ni).Lf*dfdu*m.NLTs(ni).L;
                if mod(m.NLTs(ni).type,2)==0
                    dFdUd = dFdUd + m.NLTs(ni).Lf*dfdud*m.NLTs(ni).L;
                end
            end
        end
    end
end
