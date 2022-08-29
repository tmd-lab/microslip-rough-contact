function prev = COPY_HISTORY(m, Uwxa, h, Nt)
% Function takes the final results from an EPMC Run at high amplitude and
% recreates the prev structure based on that instant. 
%
% Inputs:
%   m - same as for EPMCRESFUN
%   Uwxa - same as for EPMCRESFUN
%
% Output: 
%   prev - history structure for friction model.


    Nhc = sum((h==0)+2*(h~=0));

    la = Uwxa(end);
    A = 10^la;
    
    h0 = double(h(1)==0);
    Asc = kron([ones(h0,1); A*ones(Nhc-h0,1)], ones(m.Ndofs,1));
    
    t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);  
    
    prev = cell(length(m.NLTs), 1);
    
    warning('Set to serial version');
    for ni=1:length(m.NLTs)
%     parfor ni=1:length(m.NLTs)
        Unl = (m.NLTs(ni).L*reshape(Asc.*Uwxa(1:end-3), m.Ndofs, Nhc))';  % Nhc x Ndnl
        
        unlt = TIMESERIES_DERIV(Nt, h, Unl, 0);  % Nt x Ndnl
        
        if mod(m.NLTs(ni).type, 2)==0  % Instantaneous force
            error('Why are you updating history for an instantaneous force?');
        else  % Hysteretic force
            
            usePrev = false;
            
            if(isfield(m.NLTs(ni), 'func_init'))
                
                usePrev = true;
                prev{ni} = m.NLTs(ni).prev;
                
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

                prev{ni}.ddeltamduxynh = ddeltamduxynh;
                prev{ni}.duxyn0duxynh = zeros(3,3,size(cst, 2));
                prev{ni}.duxyn0duxynh(:, :, 1) = eye(3); % Zeroth harmonic sets initial tangential position
                prev{ni}.duxyn0duxynh(3, 3, :) = prev{ni}.ddeltamduxynh; % Normal is the derivative at the max
    

                [~, ~, prev{ni}] = m.NLTs(ni).func_init(uxyn_init, uxyn0, h, t(unmax_ind), prev{ni});
                
            end
            
            assert(usePrev, 'Only implemented for using a prev structure.');
            

        end
    end
    

end