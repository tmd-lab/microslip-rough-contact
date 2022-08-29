function [R, dRdUwx, dRda] = EPMCRESFUN(m, Uwxa, Fl, h, Nt, optEPMC, varargin)
% Residual Function for the Extended Periodic Motion Concept (EPMC)
%
% Function call:
%   m.EPMCRESFUN(Uwxa, Fl, h, Nt, tol, varargin)
%
% Inputs:
%   m    : model object that it is called on
%   Uwxa : Unknowns vector [U0; U1c/A; U1s/A; ... ; 
%                           omega (freq); xi (damping); log10(A=amplitude)
%   Fl   : Forcing vector [Fstatic; Fdyn] - Fstatic has size(U0) and is an
%           applied prestress force - Fdyn has size(U)-size(U0) and is used
%           for the phase constraint
%   h    : List of harmonics to use (e.g., h = 0:3)
%   Nt   : Number of time steps to use for AFT procedure
%   optEPMC : Settings for AFT convergence:
%               .tol - Absolute tolerance of AFT friction traction / force
%               .maxAFT - maximum number of allowed cycles for AFT
%                           procedure to converge. Simple slider models
%                           expect that 2 should be sufficient.
%
% Outputs:
%   R      : Residual Vector
%   dRdUwx : Jacobian of Residual w.r.t. unknowns
%   dRda   : Jacobian w.r.t. amplitude for arc length continuation
%
% NOTES:
%   1. dRdUwx does not support velocity dependent nonlinear forces where
%       dFnl/dw =/= 0.

    Nhc = sum((h==0)+2*(h~=0));

    la = Uwxa(end);
    A = 10^la;
    dAdla = A*log(10);
    
    % Amplitude so that amplitude normalization can be selectively removed:
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
    
    FNL  = zeros(m.Ndofs*Nhc, 1);
    dFNL = zeros(m.Ndofs*Nhc);
    
    % Remove Broadcast of full m object
    NLTs = m.NLTs;
    Ndofs = m.Ndofs;
    
%     warning('Set to serial version');
%     for ni=1:length(NLTs)
    parfor ni=1:length(NLTs)
        Unl = (NLTs(ni).L*reshape(Asc.*Uwxa(1:end-3), Ndofs, Nhc))';  % Nhc x Ndnl
        Ndnl = size(NLTs(ni).L, 1);
        
        unlt = TIMESERIES_DERIV(Nt, h, Unl, 0);  % Nt x Ndnl
        unldot = w*TIMESERIES_DERIV(Nt, h, Unl, 1);  % Nt x Ndnl
        
        if mod(NLTs(ni).type, 2)==0  % Instantaneous force
            [ft, dfdu, dfdud] = NLTs(ni).func(t, unlt, unldot);
            % (Nt,Ndnl); (Nt,Ndnl); (Nt,Ndnl) (point-wise)
            F = GETFOURIERCOEFF(h, ft);
            J = zeros(size(NLTs(ni).L,1)*Nhc, size(NLTs(ni).L,1)*Nhc);
            dFdU = reshape(GETFOURIERCOEFF(h, reshape(dfdu.*permute(cst, [1, 3, 2]), Nt, Ndnl*Nhc)),...
                Nhc, Ndnl, Nhc);
            for di=1:Ndnl
                J(di:Ndnl:end, di:Ndnl:end) = dFdU(:, di, :);
            end
        else  % Hysteretic force
            ft = zeros(Nt, Ndnl);
            dfdu = zeros(Nt, Ndnl, Ndnl, Nhc);
            
            its = 0;
            
            usePrev = false;
            if(isfield(NLTs(ni), 'func_init'))
                
                usePrev = true;
                prev = NLTs(ni).prev;
                
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
    
                [~, ~, prev] = NLTs(ni).func_init(uxyn_init, uxyn0, h, t(unmax_ind), prev);
                
            end
            
            
            while its==0 || max(abs(fprev-ft(end, :)))> optEPMC.tol && its <= optEPMC.maxAFT
                fprev = ft(end, :);
                for ti=1:Nt
                    if(usePrev)
                        [ft(ti,:), dfdu(ti,:,:,:), prev] = NLTs(ni).func(unlt(ti,:), h, t(ti), cst(ti, :), prev);
                    else
                        tm1 = mod(ti-2, Nt)+1;
                        [ft(ti,:), dfdu(ti,:,:,:)] = ...
                            NLTs(ni).func(t(ti), unlt(ti,:), h, t(tm1), ...
                            unlt(tm1,:), ft(tm1,:), dfdu(tm1,:,:,:));
                    end
                end
                its = its+1;
%                keyboard 
            end
%             fprintf('%d\n', its);
            F = GETFOURIERCOEFF(h, ft);
            J = zeros(size(NLTs(ni).L,1)*Nhc, size(NLTs(ni).L,1)*Nhc);
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
        
        if NLTs(ni).type<=5  % Self adjoint forcing
            FNL = FNL + reshape(NLTs(ni).L'*F', Nhc*Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), NLTs(ni).L')*J*kron(eye(Nhc), NLTs(ni).L);
	    else  % Non-self adjoint forcing
            FNL = FNL + reshape(NLTs(ni).Lf*F', Nhc*Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), NLTs(ni).Lf)*J*kron(eye(Nhc), NLTs(ni).L);                
        end
    end
    
    % Residue
    h0 = double(h(1)==0);
    Fstat = kron([ones(h0,1); zeros(Nhc-h0,1)], ones(Ndofs,1)).*Fl;  % Only static loads if it's there
    Fdyn = kron([zeros(h0,1); ones(Nhc-h0, 1)], ones(Ndofs,1)).*Fl;  % Only dynamic loads (for phase constraint)
    
%   Use all harmonic amplitudes as mode shape and normalize by sum of all harmonics
%     R = [E*(Asc.*Uwxa(1:end-3))+FNL-Fstat;
%         Uwxa(1:end-3)'*kron(diag([0 ones(1,Nhc-1)]), m.M)*Uwxa(1:end-3)-1;
%         Fdyn'*Uwxa(1:end-3)];
%     dRdUwx = [(E+dFNL)*diag(Asc), dEdw*(Asc.*Uwxa(1:end-3)), dEdxi*(Asc.*Uwxa(1:end-3));
%         2*Uwxa(1:end-3)'*kron(diag([0 ones(1,Nhc-1)]), m.M), 0, 0;
%         Fdyn', 0, 0];
%     dRda = [(E+dFNL)*(dAscdA.*Uwxa(1:end-3))*dAdla;
%         0; 
%         0];

%   Use first harmonic amplitudes as mode shape
    R = [E*(Asc.*Uwxa(1:end-3))+FNL - Fstat;
        Uwxa(Ndofs*h0+(1:2*Ndofs))'*blkdiag(m.M,m.M)*Uwxa(Ndofs*h0+(1:2*Ndofs))-1;
        Fdyn'*Uwxa(1:end-3)];
    
    dRdUwx = [(E+dFNL)*diag(Asc), dEdw*(Asc.*Uwxa(1:end-3)), dEdxi*(Asc.*Uwxa(1:end-3));
        zeros(1, h0*Ndofs), 2*Uwxa(Ndofs*h0+(1:2*Ndofs))'*kron(eye(2), m.M), zeros(1,Ndofs*(Nhc-3)), 0, 0;
        Fdyn', 0, 0];
    
    dRda = [(E+dFNL)*(dAscdA.*Uwxa(1:end-3))*dAdla;
        0; 
        0];
    
    if ~isempty(varargin)
        dRdUwx = [dRdUwx dRda];
    end
end