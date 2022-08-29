classdef MDOFGEN
    %MDOFGEN creates a generic mdof model class with given Mass, Stiffness
    %and Damping matrices and specified nonlinearities
    properties
        Ndofs  % Number of DoFs
        
        M   % Mass
        K   % Stiffness
        C   % Damping
        
        L  % displacement transformation matrix
        
        NLTs % Nonlinear functions
        Rsc
    end
    
    methods
        function m = MDOFGEN(M, K, C, L)
        %MDOFGEN initializes the MDOFGEN object
        % 
        % USAGE:
        %   m = MDOFGEN(M, K, C, L);
            m.M = M;
            m.K = K;
            m.C = C;
            m.L = L;
            
            m.Ndofs = size(m.M,1);
        end
        
        function m = SETNLFUN(m, type, Ls, func, Lf, prev, func_init)
        %SETNLFUN sets nonlinear functions to be applied to the model. Two
        %types are allowed: 'inst' and 'hyst' (see below)
        %
        % USAGE:
        %   type    : a+b (possible values: {4, 6, 5, 7}),
        %               with "a" being,
        %     [1] for instantaneous nonlinearity (in time domain) or
        %     [2] for hysteretic nonlinearity (in time domain)
        %               AND "b" being,
        %     [3] for self adjoint force application (F = Ls'*func(Ls*U))
        %     [5] for non-self adjoint force application (F = Lf*func(Ls*u)
        %   Ls      : (Nldofs, Ndofs) selection matrix
        %   func    : if 'inst', [ft, dfdut, dfdudt] = @(t, u, ud) func;
        %                           returning (Nldofs-qtties); assumed to
        %                           be vectorized in time and u (arranged
        %                           properly)
        %             if 'hyst', [ft, dfdut, dfdudt] = @(t, u, tp, up, fp, h)
        %                   func; returning (Nldofs-qtties); assumed to be
        %                   vectorized in u only.
        %   Lf      : [reqd only for type = 7,8] (Ndofs,Nldofs)
        %               "integration" matrix
            nlfun.L    = Ls;
            nlfun.func = func;
            nlfun.type = type;
            if nlfun.type > 5  % Non-self adjoint forcing
                nlfun.Lf = Lf;
            end
            
            if nargin>5  % Adding "generalized history structure" (for hysteretic nonlinearities)
                nlfun.prev = prev;
                nlfun.func_init = func_init; %Function for initializing history before the time marching.
            end
            
            m.NLTs = [m.NLTs; nlfun];
        end
        
        function [m, Fshape] = ATTACHSHAKER(m, E, A, B, K, Kd, Nshape)
        %ATTACHSHAKER stores the information to attach (one or more) shaker
        %model(s) to the model. Only applicable to explicit (RK) time
        %stepping solvers.
        %  Shaker-attached Model assumed to be written in state-space form as,
        %       E Xd = A X + B U - K*(X-Nshape Y) - Kd*(X - Nshape Yd)
        %       M Ydd + C Yd + K Y + Fnl + Nshape'*K*(Nshape' Y - X) + Nshape'*Kd*(Nshape' Yd - X)
        %   
        %   USAGE:
        %       E, A    : (nX, nX) A matrix 
        %       B       : (nX, nU) B matrix 
        %       K, Kd   : (nX, nX) K matrix 
        %       Nshape  : (nX, nY) Nshape matrix
            
            n = size(A,1);
            
            m.M = blkdiag(m.M, zeros(n));
            m.C = [m.C+Nshape'*Kd*Nshape, -Nshape'*Kd;
                -Kd*Nshape, E];
            m.K = [m.K+Nshape'*K*Nshape, -Nshape'*K;
                -K*Nshape, -A+K+Kd];
            Fshape = [zeros(m.Ndofs,size(B,2)); B];
            
            if ~isempty(m.NLTs)
                for i=1:length(m.NLTs)
                    m.NLTs(i).L = [m.NLTs(i).L, zeros(size(m.NLTs(i).L,1), n)];
                    if m.NLTs(i).type > 5
                        m.NLTs(i).Lf = [m.NLTs(i).Lf; zeros(n, size(m.NLTs(i).Lf,2))];
                    end
                end
            end
            
            m.Ndofs = size(m.M,1);
        end
        
        function m = ATTACHSHAKER2(m, Ms, Cs, Ks, K, Kd, Nshape)
        %ATTACHSHAKER stores the information to attach (one or more) shaker
        %model(s) to the model. Only applicable to explicit (RK) time
        %stepping solvers.
        %  Shaker-attached Model assumed to be written in second order form as
        %       Ms Xdd + Cs Xd + Ks X + K*(X-Nshape*Y) + Kd*(Xd-Nshape*Yd) = 0
        %       M Ydd + C Yd + K Y + Fnl + Nshape'*K*(Nshape' Y - X) + Nshape'*Kd*(Nshape' Yd - X) = 0
        %   
        %   USAGE:
        %       Ms,Cs,Ks: (nX, nX) A matrix 
        %       K, Kd   : (nX, nX) K matrix 
        %       Nshape  : (nX, nY) Nshape matrix
        
            n = size(Ms,1);
            
            m.M = blkdiag(m.M, Ms);
            m.C = [m.C+Nshape'*Kd*Nshape, -Nshape'*Kd;
                -Kd*Nshape, Cs+Kd];
            m.K = [m.K+Nshape'*K*Nshape, -Nshape'*K;
                -K*Nshape, Ks+K];
            
            if ~isempty(m.NLTs)
                for i=1:length(m.NLTs)
                    m.NLTs(i).L = [m.NLTs(i).L, zeros(size(m.NLTs(i).L,1), n)];
                    if m.NLTs(i).type > 5
                        m.NLTs(i).Lf = [m.NLTs(i).Lf; zeros(n, size(m.NLTs(i).Lf,2))];
                    end
                end
            end
            
            m.Ndofs = size(m.M,1);
        end        
    end
end

