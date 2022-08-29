function [unlt_quad, fnlt_quad] = EPMC_HYST(m, Uwxa, Fl, h, Nt, tol, varargin)

    Nhc = sum((h==0)+2*(h~=0));

    la = Uwxa(end);
    A = 10^la;
    dAdla = A*log(10);
    
    h0 = double(h(1)==0);
    Asc = kron([ones(h0,1); A*ones(Nhc-h0,1)], ones(m.Ndofs,1));
    dAscdA = kron([zeros(h0,1); ones(Nhc-h0,1)], ones(m.Ndofs,1));
    
    xi = Uwxa(end-1);
    w = Uwxa(end-2);
    
    [E, dEdw] = HARMONICSTIFFNESS(m.M, m.C-xi*m.M, m.K, w, h);
    dEdxi = HARMONICSTIFFNESS(m.M*0, -m.M, m.M*0, w, h);
    
    t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);  
    sct = w*TIMESERIES_DERIV(Nt, h, eye(Nhc), 1);
    
    unlt_quad = cell(length(m.NLTs), 1);
    fnlt_quad = cell(length(m.NLTs), 1);
    
%     warning('Set to serial version');
%     for ni=1:length(m.NLTs)
    parfor ni=1:length(m.NLTs)
        Unl = (m.NLTs(ni).L*reshape(Asc.*Uwxa(1:end-3), m.Ndofs, Nhc))';  % Nhc x Ndnl
        Ndnl = size(m.NLTs(ni).L, 1);
        
        unlt = TIMESERIES_DERIV(Nt, h, Unl, 0);  % Nt x Ndnl
        unldot = w*TIMESERIES_DERIV(Nt, h, Unl, 1);  % Nt x Ndnl
        
        unlt_quad{ni} = unlt;
        
        if mod(m.NLTs(ni).type, 2)==0  % Instantaneous force
            [ft, dfdu, dfdud] = m.NLTs(ni).func(t, unlt, unldot);
            % (Nt,Ndnl); (Nt,Ndnl); (Nt,Ndnl) (point-wise)
            
            fnlt_quad{ni} = ft;
        else  % Hysteretic force
            ft = zeros(Nt, Ndnl);
            dfdu = zeros(Nt, Ndnl, Ndnl, Nhc);
            
            its = 0;
            
            usePrev = false;
            if(isfield(m.NLTs(ni), 'func_init'))
                
                usePrev = true;
                prev = m.NLTs(ni).prev;
                
%                 % Find maximum displacement
%                 uxyn_init = unlt(1, :);
%                 [uxyn_init(3), unmax_ind] = max(unlt(:, 3));
%                 uxyn0 = unlt(1, :);
% 
%                 % derivative of maximum normal w.r.t. the harmonic coefficients
%                 ddeltamduxynh = reshape(cst(unmax_ind, :), 1, 1, []);
% 
%                 prev.ddeltamduxynh = ddeltamduxynh;
%                 prev.duxyn0duxynh = reshape(cst(1, :), 1, 1, []).*eye(3);


                % Start Tangential Models with Zeroth Harmonic 
                % displacements (phase invariance if does not slip).
                uxyn_init = Unl(1, :); 
                [uxyn_init(3), unmax_ind] = max(unlt(:, 3));
                uxyn0 = uxyn_init;

                % derivative of maximum normal w.r.t. the harmonic coefficients
                ddeltamduxynh = reshape(cst(unmax_ind, :), 1, 1, []);

                prev.ddeltamduxynh = ddeltamduxynh;
                prev.duxyn0duxynh = zeros(3,3,size(cst, 2));
                prev.duxyn0duxynh(:, :, 1) = eye(3); % Zeroth harmonic sets initial tangential position
                prev.duxyn0duxynh(3, 3, :) = prev.ddeltamduxynh; % Normal is the derivative at the max
    

                [~, ~, prev] = m.NLTs(ni).func_init(uxyn_init, uxyn0, h, t(unmax_ind), prev);
                
            end
            
            
            
            while its==0 || max(abs(fprev-ft(end, :)))>tol && its < 10
                fprev = ft(end, :);
                for ti=1:Nt
                    if(usePrev)
                        [ft(ti,:), dfdu(ti,:,:,:), prev] = m.NLTs(ni).func(unlt(ti,:), h, t(ti), cst(ti, :), prev);
                    else
                        tm1 = mod(ti-2, Nt)+1;
                        [ft(ti,:), dfdu(ti,:,:,:)] = ...
                            m.NLTs(ni).func(t(ti), unlt(ti,:), h, t(tm1), ...
                            unlt(tm1,:), ft(tm1,:), dfdu(tm1,:,:,:));
                    end
                end
                its = its+1;
%                keyboard 
            end
            
            fnlt_quad{ni} = ft;

        end
        
    end
    

end