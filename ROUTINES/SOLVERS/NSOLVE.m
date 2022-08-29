function [U, R, eflag, it, dRc, reu] = NSOLVE(func, U0, varargin)
%NSOLVE Uses Newton iterations to solve
%
% USAGE:
%  [U, R, eflag, it, jc] = NSOLVE(func, U0, opts);
% INPUTS:
%  func		: Function handle [F, dFdU] = func(u);
%  U0		: Initial Guess
%  opts		: Options structure with,
% 	reletol (float)
% 	etol    (float)
% 	utol    (float)
% 	rtol    (float)
% 	Display (boolean)
% OUTPUTS:
%  U		:
%  R		:
%  eflag	:
%  it		:
%  Jc		:
% error('This is obsolete now - please use GSOLVE.m')
  %% Default Parameters
  opts = struct('reletol', 1e-6, 'rtol', 1e-6, 'etol', 1e-6, 'utol', ...
                1e-6, 'Display', false, 'Dscale', ones(size(U0)), ...
	       'ITMAX', 10, 'lsrch', 0, 'crit', 12);
  if nargin==3
      nflds = fieldnames(varargin{1});
      for i = 1:length(nflds)
	opts.(nflds{i}) = varargin{1}.(nflds{i});
      end
  end  

  %% Starting Point Residual, Jacobian, and error metrics
  Nu = length(U0);
  U0 = U0./opts.Dscale(1:Nu);
  
  [R0, dR0] = func(opts.Dscale(1:Nu).*U0);
  dU0 = -(dR0*diag(opts.Dscale(1:Nu)))\R0;
  e0  = abs(R0'*dU0);
  
  if (e0 < eps)
    e0 = 1.0;
    R0 = ones(Nu, 1);
  end
  r0 = R0'*R0;
  u0 = dU0'*dU0;
  
  %% Termination Flag
    flagfun = @(e, e0, r, u) 4*(e<opts.etol) + 8*(r<opts.rtol) +...
        2*(e/e0<opts.reletol) + 1*(u<opts.utol) + 16*(e/eps<1e3);
    
  %% Initializing for iterations  
  r   = r0;
  e   = e0;
  u   = u0;

  reu = [r e u];
  
  dRc = dR0*diag(opts.Dscale(1:Nu));
  R   = R0;
  U   = U0;
  dU  = dU0;
  it  = 0;

  eflag = flagfun(e, e0, r, u); 
  if opts.Display
    fprintf('ITN, E, E/E0, r, du\n%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, eflag);
  end
  eflagp = 0;

  %% Start Iterations  
  while (eflag<opts.crit || it<1) && it<=opts.ITMAX
      % Line search code 
      if opts.lsrch~=0          
          ls_e0 = e;
          
%           for il=1:opts.lsrch
%             [R0, J0] = func(opts.Dscale(1:Nu).*(U+dU));
%             ls_e1 = R0'*(-(J0*diag(opts.Dscale(1:Nu)))\R0);
%             if ls_e0.*ls_e1<=0
%                 ls_s = ls_e0/(ls_e0-ls_e1);
%                 dU = ls_s*dU;
%             else
%                 break;
%             end
%           end
          
          % Dog-Leg Algorithm
          dUgn = dU;  % Gauss-Newton step
          al   = (R'*(dRc*dRc')*R)/(R'*(dRc*dRc')*(dRc*dRc')*R);
          dUc  =  -al*dRc'*R;  % Cauchy step
          
          [R0, dR0] = func(opts.Dscale(1:Nu).*(U+dUc));
          ls_ec = R0'*(-(dR0*diag(opts.Dscale(1:Nu)))\R0);
              
          [R0, dR0] = func(opts.Dscale(1:Nu).*(U+dUgn));
          ls_egn = R0'*(-(dR0*diag(opts.Dscale(1:Nu)))\R0);
              
          ls_l = ls_ec/(ls_ec-ls_egn);  % minimizer
          
%           ls_l = min(ls_l, 1);
%           ls_l = max(ls_l, 0);
          
%           fprintf('%e\t', ls_l);

          dU = dUc + ls_l*(dUgn-dUc);
      end
      % Line search code 
      
    U  = U + dU;
    it = it+1;
    
    [R, dRc] = func(opts.Dscale(1:Nu).*U);
%     dU = (-Jc\R)./opts.Dscale(1:Nu);
    dRc = dRc*diag(opts.Dscale(1:Nu));
    dU = -dRc\R;
    
    e = abs(R'*dU);
    r = R'*R;
    u = dU'*dU;
    
    reu = [reu; [r e u]];

    eflag = flagfun(e, e0, r, u); 
    if opts.Display
      fprintf('%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, eflag);
    end
    
    if it>opts.ITMAX
%       eflag = 0; % Change to allow the true flag to be passed out even if
%       it does not meet the criteria for convergence that was given to end
%       this loop early. 
      break;
    end
  end
  
  % Rescale Solution
  U = opts.Dscale(1:Nu).*U;
  dRc = dRc*diag(1./opts.Dscale(1:Nu));
  
  if eflag == 0
    disp('No Convergence : Returning')
    return;
  end

  if opts.Display
    disp('Successful Campaign!')
  end
end
