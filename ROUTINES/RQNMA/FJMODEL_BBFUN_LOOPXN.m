function BB = FJMODEL_BBFUN_LOOPXN(pars, mdi, Qamps, K, M, X0, Fv, L, QuadMats, MESH, CFUN, Nqp, opt, optFJ, varargin)
%FJMODEL_BBFUN_LOOPXN returns the backbone after cycling hysteresis loops
%to achieve convergence
%
% Inputs: 
%   pars - model parameters
%   mdi - mode of interest (total mode number not just counting bending
%           modes)
%   Qamps - modal amplitudes to create the backbone at
%   K - stiffness matrix
%   M - mass matrix
%   X0 - initial guess for the prestress solution
%   Fv - prestress force Vector (scaled to appropriate magnitude)
%   L - Null space projection
%   QuadMats - has Matrices Q and T to transform to quad points and
%           integrate tractions to nodal forces
%   MESH - FE Mesh details
%   CFUN - contact function to be evaluated at quad points
%   Nqp - number of quadrature points per loading/unloading portion of a
%       cycle
%   opt - fsolve options
%   optFJ - settings for this function (see below)
%   varargin - settings for conditioning improvements, implementation is
%           untested/not currentely used.
%
% Results:
%   BB - backbone structure with history of the hysteresis loops
% 
% Initializing Settings Example / Details
%
% optFJ.viscousDamping    = 0; %magnitude of viscous damping to add (e.g., 1e-3 or 0)
% optFJ.repeatLoop        = 3; %Max number of loop iterations including final one
% optFJ.zRelTol           = -1; %damping relative tolerance, negative means
%                                   it will be forced to do the max number of iterations.
% optFJ.calcQuads         = true; %Calculate all quad points on initial iterations, 
%                                   needed if using a relative tolerance or a penalty type model
% optFJ.RC_prev           = 2; % Rough contact history variables are used
% optFJ.Nqp_heights       = 100; % Number of asperity heights to integrate 
%                                   over, used if optFJ.RC_prev > 0
% optFJ.Nqp_radius        = 100; % Number of asperity contact radii to integrate 
%                                   over, used if optFJ.RC_prev > 0
% optFJ.frictionlessPrestress = true; % Do not consider tangential friction
%                                           in prestress
% optFJ.CFUN_PRE          = CFUN_PRE; % Function for forces with 0
%                                       tangential forces to be used in prestress
%
% NOTES:
%   1. Parameter derivatives are not implemented. Some places pass around
%   variables for these derivatives to support expected inputs/outputs of
%   older verious of the code. The multiple steps of fsolve tended to 
%   2. Normal Plasticity Models are not supported due to a lack of 
%   convergence for RQNMA solution method. The history variables are not 
%   initialized here. 


    %% Prestress Analysis
    
    if(optFJ.RC_prev == 1 )
        
        error('This Rough Contact History Type is no longer supported');
        
    elseif(optFJ.RC_prev == 2 )
        
        prevS1.tx0 = zeros(optFJ.Nqp_heights, optFJ.Nqp_radius); %traction for each slider of each asperity height - X
        prevS1.ty0 = zeros(optFJ.Nqp_heights, optFJ.Nqp_radius); %traction for each slider of each asperity height - X
        prevS1.rq0 = ones(optFJ.Nqp_heights, 1)*linspace(0, 1, optFJ.Nqp_radius); %radii for each slip radius used
        prevS1.uxyw0 = ones(optFJ.Nqp_heights, 1)*[0, 0, 0]; %displacement at previous step. 
        
        prevS = repmat({prevS1}, size(QuadMats.Q, 1), 1);

    else
        prevS.uxyntxyn = zeros(size(QuadMats.Q,1),6);
        
        % Parameter derivatives are not supported, but inputs are passed
        % around to prevent older versions from breaking
        prevS.duxdp = zeros(size(QuadMats.Q,1),numel(pars));
        prevS.duydp = zeros(size(QuadMats.Q,1),numel(pars));
        prevS.dtxdp = zeros(size(QuadMats.Q,1),numel(pars));
        prevS.dtydp = zeros(size(QuadMats.Q,1),numel(pars));
        prevS.uxykp0 = zeros(size(prevS.uxyntxyn, 1), 2); %Tracks reference point of kp spring
    end
    
    if(optFJ.frictionlessPrestress)
        CFUNS = optFJ.CFUN_PRE;
    else
        CFUNS = CFUN;
    end
    
    if isempty(varargin)
        [Xstat, eflg] = fsolve(@(X) NLRES_TMP(CFUNS, [X;0], K, L, Fv, Fv*0, prevS, pars, 0, QuadMats, MESH), X0, opt);
    else
        [Xstat, eflg] = fsolve(@(X) NLRES_TMP(CFUNS, [X;0], K, L, Fv, Fv*0, prevS, pars, 0, QuadMats, MESH, ones(size(X0)), varargin{1}), X0, opt);
    end
    if eflg<0
        disp(pars)
        error('Non-convergent prestress!');
    end
    
    %These use CFUN rather than CFUNS to give an accurate stiffness matrix
    %for the following eigen analysis, therefore these do not really
    %provide a good way to get derivatives w.r.t. parameters of prestress.
    if isempty(varargin)
        [~,dRstat,~,~,dRdpstat] = NLRES_TMP(CFUNS, [Xstat;0], K, L, Fv, Fv*0, prevS, pars, 0, QuadMats, MESH);
    else
        [~,dRstat,~,~,dRdpstat] = NLRES_TMP(CFUNS, [Xstat;0], K, L, Fv, Fv*0, prevS, pars, 0, QuadMats, MESH, ones(size(X0)), varargin{1});
    end
    dXdpstat = -dRstat\dRdpstat; 
    
    % Update History Variables
    if(optFJ.RC_prev > 0)
        [tx, ty, tn, ~, ~, ~, ~, ~, ~, ~, dtxdp, dtydp, dtndp, prevS] = ...
            CFUNS(QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xstat, 3, [])', 0, pars, prevS);
    else
        [tx, ty, tn, ~, ~, ~, ~, ~, ~, ~, dtxdp, dtydp, dtndp] = ...
            CFUNS(QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xstat, 3, [])', 0, pars, prevS);

        prevS.uxyntxyn = [QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xstat, 3, [])', tx, ty, tn];
        prevS.uxykp0 = zeros(size(prevS.uxyntxyn, 1), 2); %Tracks reference point of kp spring
    end
    
    % recalculate the residual stiffness with the prev set slipped
    % prestress. This gives a Jacobian that can be used in initial 
    % eigenanalysis
    if isempty(varargin)
        [~,dRstat,~,~,~] = NLRES_TMP(CFUN, [Xstat;0], K, L, Fv, Fv*0, prevS, pars, 0, QuadMats, MESH);
    else
        [~,dRstat,~,~,~] = NLRES_TMP(CFUN, [Xstat;0], K, L, Fv, Fv*0, prevS, pars, 0, QuadMats, MESH, ones(size(X0)), varargin{1});
    end
    
    %% Modal Analysis
    [Vst, Wst] = eigs(dRstat, M, 10, 'SM');
    [Wst,si] = sort(sqrt(diag(Wst)));
    Vst = Vst(:,si)./sqrt(diag(Vst(:,si)'*M*Vst(:,si))');
    
    Vmi = Vst(:, mdi);
    Wmi = Wst(mdi);
    %% Initialization 
    BB.Q = Qamps;
    BB.W = zeros(size(BB.Q)); % Averaged Frequency
    BB.W1 = zeros(size(BB.Q)); % Frequency for convergence checking
    BB.W2 = zeros(size(BB.Q));
    BB.W3 = zeros(size(BB.Q));
    BB.Z = zeros(size(BB.Q)); % Damping Factor (Fraction)
    BB.D = zeros(size(BB.Q)); % Dissipation
    BB.Xavg = zeros(length(Xstat)+1, length(BB.Q)); % Mode Shape
    BB.R1 = zeros(length(Xstat)+1, length(BB.Q)); % Finaly Residuals
    BB.R2 = zeros(length(Xstat)+1, length(BB.Q));
    BB.R3 = zeros(length(Xstat)+1, length(BB.Q));
    
    % Full Displacement History
    BB.XS = Xstat;
    BB.X1 = zeros(length(Xstat)+1, optFJ.repeatLoop, length(BB.Q));
    BB.X2 = zeros(length(Xstat)+1, optFJ.repeatLoop, length(BB.Q));
    BB.X3 = zeros(length(Xstat)+1, optFJ.repeatLoop, length(BB.Q));
    
    % Local Hysteresis loop history
    BB.uxynHist = zeros(size(tx, 1), 3, (optFJ.repeatLoop*2+1)*Nqp, length(BB.Q)); %store the displacements for each RQNMA evaluation
    BB.txynHist = zeros(size(tx, 1), 3, (optFJ.repeatLoop*2+1)*Nqp, length(BB.Q)); %store the tractions at each RQNMA evaluation as well. 
    BB.qHist = zeros((optFJ.repeatLoop*2+1)*Nqp, length(BB.Q)); %store q at each of the history points for other vars.
    BB.RHist = zeros((optFJ.repeatLoop*2+1)*Nqp, length(BB.Q)); %store the residual norms
%     uxynHist_ind = 1;
    
    % Repeated Loop Convergence Details
    BB.Wconverge = zeros(length(BB.Q), optFJ.repeatLoop, 3); %rows = amplitude levels, columns iterations, pages W1;W2;W3
    BB.Dconverge = zeros(length(BB.Q), optFJ.repeatLoop); %rows = amplitude levels, columns iterations
    BB.Zconverge = zeros(length(BB.Q), optFJ.repeatLoop); %rows = amplitude levels, columns iterations
    
    
    %% LGL Quadrature Points
    
    [xi, wi] = LGLWT(Nqp, 0, 1);
    qbb = zeros(size(xi));  abb = zeros(size(xi));
    qul = zeros(size(xi));  aul = zeros(size(xi));
    qrl = zeros(size(xi));  arl = zeros(size(xi));
    
    
    BB.q_hyst = zeros(length(xi)*3, length(BB.Q));
    BB.a_hyst = zeros(length(xi)*3, length(BB.Q));
    
%     opt.Display = 'off';
    
    %% RQNMA
    
    fprintf('Starting RQNM....................................\n');
    for i=1:length(BB.Q)
        
        uxynHist_ind = 1;
        prev = prevS;
        
        if(optFJ.calcQuads) %if calculating quadrature do BB since it is probably a penalty model.
            xi_curr = xi;
            wi_curr = wi;
        else
            xi_curr = [0, 1];
            wi_curr = [0, 0];
        end
        
        % Backbone Loading
        for k=2:length(xi_curr)
            qbb(k) = BB.Q(i)*xi_curr(k);
            X0 = [Xstat+Vmi*qbb(k); Wmi^2];
            if isempty(varargin)
                Xb = fsolve(@(X) NLRES_RQNM_TMP(CFUN, [X; log10(qbb(k))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH), X0, opt);
                [R3, ~, ~, ~, ~] = NLRES_RQNM_TMP(CFUN, [Xb; log10(qbb(k))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH);
            else
                Xb = fsolve(@(X) NLRES_RQNM_TMP(CFUN, [X; log10(qbb(k))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH, ones(size(X0(1:end-1))), varargin{1}), X0, opt);
                [R3, ~, ~, ~, ~] = NLRES_RQNM_TMP(CFUN, [Xb; log10(qbb(k))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH, ones(size(X0(1:end-1))), varargin{1});
            end
            abb(k) = qbb(k)*Xb(end);
            
            %Update prev
            uxyn = QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xb(1:end-1), 3, [])';
            [tx, ty, tn, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, prev] = ...
                CFUN(uxyn, 0, pars, prev);
            
            %store history of loading
            BB.uxynHist(:, :, uxynHist_ind, i) = uxyn; %store the displacements for each RQNMA evaluation
            BB.txynHist(:, :, uxynHist_ind, i) = [tx, ty, tn]; %store the tractions at each RQNMA evaluation as well. 
            BB.qHist(uxynHist_ind, i) = qbb(k);
            BB.RHist(uxynHist_ind, i) = norm(R3);
            uxynHist_ind = uxynHist_ind + 1;
            
        end
        Xd = Xb;
        
        %%%%%% START of ITERATION LOOP
        for j = 1:optFJ.repeatLoop
        Xb = Xd; %update backbone to the previous forward iteration
        R1 = R3; %initial residuals stored from previous loop.
        BB.D(i) = 0; %reset dissipation to zero to calculate on every loop iteration
        
        BB.X1(:, j, i) = Xb;
        
        if(optFJ.calcQuads || j == optFJ.repeatLoop) % if calculating quadrature or if at the last iteration use actual quadrature
            xi_curr = xi;
            wi_curr = wi;
        else
            xi_curr = [0, 1];
            wi_curr = [0, 0];
        end
        
        
        if(optFJ.RC_prev > 0)
            % Rough contact history information is already updated
        else
            uxyn = QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xb(1:end-1), 3, [])';

            [tx, ty, tn, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                CFUN(uxyn, 0, pars, prev);
        
            prev.uxyntxyn = [uxyn, tx, ty, tn];
        end
        
        % for when k = 1, don't need to calc RQNM since already have the
        % solution 
        qul(1) = BB.Q(i) - 2*BB.Q(i)*xi_curr(1);
        aul(1) = qul(1)*Xb(end);
        BB.D(i) = BB.D(i) - qul(1)*wi_curr(1)*(2*BB.Q(i))*Xb(end);
        
        
        for k=2:length(xi_curr)
            qul(k) = BB.Q(i) - 2*BB.Q(i)*xi_curr(k);
            X0 = [Xstat+Vmi*qul(k); Wmi^2];
            if isempty(varargin)
                Xd = fsolve(@(X) NLRES_RQNM_TMP(CFUN, [X; log10(abs(qul(k)))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH), X0, opt);
                [R2, dRdX, ~, ~, dRdp] = NLRES_RQNM_TMP(CFUN, [Xd; log10(abs(qul(k)))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH);
            else
                Xd = fsolve(@(X) NLRES_RQNM_TMP(CFUN, [X; log10(abs(qul(k)))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH, ones(size(X0(1:end-1))), varargin{1}), X0, opt);
                [R2, dRdX, ~, ~, dRdp] = NLRES_RQNM_TMP(CFUN, [Xd; log10(abs(qul(k)))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH, ones(size(X0(1:end-1))), varargin{1});
            end
            
            %Update prev
            uxyn = QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xd(1:end-1), 3, [])';
            [tx, ty, tn, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, prev] = ...
                CFUN(uxyn, 0, pars, prev);
            
            aul(k) = qul(k)*Xd(end);
            
            BB.D(i) = BB.D(i) - qul(k)*wi_curr(k)*(2*BB.Q(i))*Xd(end);            
            
            %store history of loading
            BB.uxynHist(:, :, uxynHist_ind, i) = uxyn; %store the displacements for each RQNMA evaluation
            BB.txynHist(:, :, uxynHist_ind, i) = [tx, ty, tn]; %store the tractions at each RQNMA evaluation as well. 
            BB.qHist(uxynHist_ind, i) = qul(k);
            BB.RHist(uxynHist_ind, i) = norm(R2);
            uxynHist_ind = uxynHist_ind + 1;
            
        end
        
        
        BB.X2(:, j, i) = Xd;
        BB.W2(i) = sqrt(Xd(end));
        
        % Hysteretic Reloading
        
        if(optFJ.RC_prev > 0)
            % History already updated
        else
            [tx, ty, tn, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                CFUN(QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xd(1:end-1), 3, [])', 0, pars, prev);
            prev.uxyntxyn = [QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xd(1:end-1), 3, [])', tx, ty, tn];
            % prev.uxykp0 continues without modification
        end
        
        %for when k = 1, don't need to calc RQNM since already have the
        %solution
        qrl(1) = -BB.Q(i) + 2*BB.Q(i)*xi_curr(1);
        arl(1) = qrl(1)*Xd(end);
        BB.D(i) = BB.D(i) + Xd(end)*qrl(1)*wi_curr(1)*(2*BB.Q(i));
        
        for k=2:length(xi_curr)
            qrl(k) = -BB.Q(i) + 2*BB.Q(i)*xi_curr(k);
            X0 = [Xstat+Vmi*qrl(k); Wmi^2];
            if isempty(varargin)
                Xd = fsolve(@(X) NLRES_RQNM_TMP(CFUN, [X; log10(abs(qrl(k)))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH), X0, opt);
                [R3, ~, ~, ~, ~] = NLRES_RQNM_TMP(CFUN, [Xd; log10(abs(qrl(k)))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH);
            else
                Xd = fsolve(@(X) NLRES_RQNM_TMP(CFUN, [X; log10(abs(qrl(k)))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH, ones(size(X0(1:end-1))), varargin{1}), X0, opt);
                [R3, ~, ~, ~, ~] = NLRES_RQNM_TMP(CFUN, [Xd; log10(abs(qrl(k)))], Xstat, dXdpstat, M, K, L, Fv, prev, pars, 0, QuadMats, MESH, ones(size(X0(1:end-1))), varargin{1});
            end
            
            
            %Update prev
            uxyn = QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xd(1:end-1), 3, [])';
            [tx, ty, tn, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, prev] = ...
                CFUN(uxyn, 0, pars, prev);
            
            arl(k) = qrl(k)*Xd(end);
            
            BB.D(i) = BB.D(i) + Xd(end)*qrl(k)*wi_curr(k)*(2*BB.Q(i));
            
            %store history of loading
            BB.uxynHist(:, :, uxynHist_ind, i) = uxyn; %store the displacements for each RQNMA evaluation
            BB.txynHist(:, :, uxynHist_ind, i) = [tx, ty, tn]; %store the tractions at each RQNMA evaluation as well. 
            BB.qHist(uxynHist_ind, i) = qrl(k);
            BB.RHist(uxynHist_ind, i) = norm(R3);
            uxynHist_ind = uxynHist_ind + 1;
            
        end
        
        % Processing / Final Storage
        BB.W1(i) = sqrt(Xb(end));
        BB.W3(i) = sqrt(Xd(end));
        
        BB.X3(:, j, i) = Xd;
        
        BB.W(i) = (BB.W3(i) + BB.W2(i))/2;
        
        BB.Z(i) = BB.D(i)/(2*pi*(BB.Q(i)*BB.W(i))^2) + optFJ.viscousDamping;
        
        BB.Xavg(:, i) = (BB.X3(:, j, i) - BB.X2(:, j, i))/2;
        
        BB.Wconverge(i, j, 1) = BB.W1(i); %rows = amplitude levels, columns iterations, pages W1;W2;W3
        BB.Wconverge(i, j, 2) = BB.W2(i); %rows = amplitude levels, columns iterations, pages W1;W2;W3
        BB.Wconverge(i, j, 3) = BB.W3(i); %rows = amplitude levels, columns iterations, pages W1;W2;W3
        BB.Dconverge(i, j) = BB.D(i); %rows = amplitude levels, columns iterations
        BB.Zconverge(i, j) = BB.Z(i); %rows = amplitude levels, columns iterations
    
        
        if(j > 1 && optFJ.calcQuads) %can start checking convergence criteria if quadrature points are used
            if( (abs(BB.Zconverge(i, j-1) - BB.Zconverge(i, j)) / abs(BB.Zconverge(i, j))) < optFJ.zRelTol)
                %simulation has converged at this amplitude level
                break; %break the j loop, continue with the i loop.
            end
        end
        
        %%%%%%%%%%%%%%%%%%End of Repeat Area Calc Iteration loop
        end
        
        BB.R1(:, i) = R1;
        BB.R2(:, i) = R2;
        BB.R3(:, i) = R3;
        fprintf('%d...', i);
        
        % These outputs are of limitted use, mainly just for debugging:
        BB.q_hyst(:, i) = [qbb; qul; qrl];
        BB.a_hyst(:, i) = [abb; aul; arl];
        
    end
    fprintf('\nCompleted RQNM....................................');
    
end