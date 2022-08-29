function [E,dEdw] = HARMONICSTIFFNESS(M,C,K,w,h)
%HARMONICSTIFFNESS Returns the harmonic stiffness and its
%derivative w.r.t the frequency w
% USAGE:
%	[E,dEdw] = HARMONICSTIFFNESS(M,C,K,w,h);
% INPUTS:
%   M,C,K	: Mass, Damping, & Stiffness matrices
%   w		: Fundamental Frequency
%   h		: List of harmonics
% OUTPUTS:
%   E		: (nd(2Nh+1)) square stiffness matrix
%   dEdw	: (nd(2Nh+1)) square derivative matrix
    
    zi = find(h==0);
    h = reshape(h,length(h),1);
    h([1 zi]) = h([zi 1]);
    nd = size(M,1);
    
    if ~isempty(zi)
        n = length(h)-1;
        E = sparse((2*n+1)*nd,(2*n+1)*nd);
        dEdw = sparse((2*n+1)*nd,(2*n+1)*nd);
        E(1:nd, 1:nd) = K;
        zi = 1;
    else
        n = length(h);
        E = sparse((2*n)*nd,(2*n)*nd);
        dEdw = sparse((2*n)*nd,(2*n)*nd);        
        zi = 0;
    end

    % Without loops (Bad on memory)
    if strcmp(class(M), 'double')
        E((nd*zi+1):end, (nd*zi+1):end) = kron(eye(2*n),K) - ...
            kron(kron(diag(h((zi+1):end)*w).^2,eye(2)),M) +...
            kron(diag(h(zi+1:end)*w),kron([0,1;-1,0],C));
        dEdw(nd*zi+1:end, nd*zi+1:end) = ...
            -kron(kron(2*w*diag(h(zi+1:end)).^2,eye(2)),M) +...
            kron(diag(h(zi+1:end)),kron([0,1;-1,0],C));
    elseif strcmp(class(M), 'single')
        h = double(h);
        E((nd*zi+1):end, (nd*zi+1):end) = double(kron(eye(2*n),K) - ...
            kron(kron(diag(h((zi+1):end)*w).^2,eye(2)),M) +...
            kron(diag(h(zi+1:end)*w),kron([0,1;-1,0],C)));
        dEdw(nd*zi+1:end, nd*zi+1:end) = ...
            double(-kron(kron(2*w*diag(h(zi+1:end)).^2,eye(2)),M) +...
                   kron(diag(h(zi+1:end)),kron([0,1;-1,0],C)));
    end
end
