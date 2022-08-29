function [Fnl] = NLEVAL(m, t, U, Udot,tol)
%NLEVAL evaluates the nonlinearities in the time domain for given set of
%points

    Nt = length(t);
    if size(U,1)==Nt
        U = U';
    end  % (Nd, Nt)
    if size(Udot,1)==Nt
        Udot = Udot';
    end
    
    Fnl = zeros(m.Ndofs, Nt);
    for ni=1:length(m.NLTs)        
        unlt = (m.NLTs(ni).L*U)';
        unldot = (m.NLTs(ni).L*U)';
        
        Ndnl = size(m.NLTs(ni).L, 1);
        
        if mod(m.NLTs(ni).type,2)==0  % Inst. force
            [ft, dfdu, dfdud] = m.NLTs(ni).func(t(:), unlt, unldot);
            % (Nt,Ndnl); (Nt,Ndnl); (Nt,Ndnl) (point-wise)
            
        else  % Hysteretic force
            ft = zeros(Nt, Ndnl);
            dfdu = zeros(Nt, Ndnl, Ndnl, 1);
            
            its  = 0;
            while its==0 || abs(fprev-ft(end, :))>tol
                fprev = ft(end, :);
                for ti=1:Nt
                    tm1 = mod(ti-2, Nt)+1;
                    [ft(ti,:), dfdu(ti,:,:,:)] = ...
                        m.NLTs(ni).func(t(ti), unlt(ti,:), 0, t(tm1), ...
                        unlt(tm1,:), ft(tm1,:), dfdu(tm1,:,:,:));
                end
                its = its+1;
            end
        end
        
        if m.NLTs(ni).type<=5
            Fnl = Fnl + m.NLTs(ni).L'*ft';
        else
            Fnl = Fnl + m.NLTs(ni).Lf*ft';
        end
    end    
    
    Fnl = Fnl';  % Ndofs x Nt
end