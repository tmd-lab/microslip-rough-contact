function [Tpres, U, Ud, Udd, m] = RKGENMARCH(m, T0, T1, dt, U0, Ud0, Fex, varargin)
%RKGENMARCH marches in time using generalized Runge Kutta scheme.
%ERROR NORM UNRELIABLE FOR SMALL DISPLACEMENT APPLICATIONS.
%
%  USAGE:
%    [X, T] = m.RKGENMARCH(T0, T1, dt, U0, Ud0, Fex, opts);
%  INPUTS:
%    T0, T1, dt	: Starting, ending, and interval times
%    U0, Ud0 	: Ndx1 displacement and velocity initial conditions
%    Fex	: (1x1)->(Ndx1) forcing function
%    opts	: [optional] Structure with parameters (defaults to rk45)
%  	a, b, bs, c, stepadapt [1], abstol [1e-6], pow [1/4],
%       maxstep [1e-3], minstep [1e-6], Display ['on'] 
%  OUTPUTS:
%    T		: 1xNt time vector
%    U, Ud, Udd	: NdxNt displacement, velocity, and acceleration 
%    m		: MDOFGEN struct (with updated hysteretic states if
%                 relevant)
  opts = struct('a', [0 0 0 0 0 0;
		      1/4 0 0 0 0 0;
		      3/32 9/32 0 0 0 0;
		      1932/2197 -7200/2197 7296/2197 0 0 0;
		      439/216 -8 3680/513 -845/4104 0 0;
		      -8/27 2 -3544/2565 1859/4104 -11/40 0], ...
		'b', [16/135 0 6656/12825 28561/56430 -9/50 2/55], ...
		'bs', [25/216 0 1408/2565 2197/4104 -1/5 0], ...
		'c', [0 1/4 3/8 12/13 1 1/2], ...
		'stepadapt', 1, ...
		'abstol', 1e-6, 'pow', 1/4, 'maxstep', 1e-3, 'minstep', 1e-6, ...
		'Display', 'on');
  if length(varargin)==1
    nflds = fieldnames(varargin{1});
      for i = 1:length(nflds)
        opts.(nflds{i}) = varargin{1}.(nflds{i});
      end
  end
  %% Solver options
  a  = opts.a;
  b  = opts.b;
  bs = opts.bs;
  c  = opts.c;
  abstol  = opts.abstol;
  pow     = opts.pow;
  maxstep = opts.maxstep;
  minstep = opts.minstep;
  ord     = length(opts.c);
  
  %% State-Space: Y0 = [U0; Ud0];

  %% Allocate vectors
  Treq = T0:dt:T1;
  Tpres = Treq;
  Nt = length(Treq);

  U   = zeros(m.Ndofs, Nt);
  Ud  = zeros(m.Ndofs, Nt);
  Udd = zeros(m.Ndofs, Nt);

  U(:, 1)  = U0;
  Ud(:, 1) = U0;
  [~, Udd(:, 1)] = ROCFUN(m, Tpres(1), U(:, 1), Ud(:, 1), Fex(Tpres(1)), Tpres(1)-dt);

  kU  = zeros(ord, m.Ndofs);
  kUd = zeros(ord, m.Ndofs);
  n  = 1;
  cflag = 0;

  if strcmp(opts.Display, 'waitbar')
    wb = waitbar(Tpres(1)/T1, sprintf('Progress: %.e/%.e dt: %.e', Tpres(1), T1, dt), ...
		 'createcancelbtn', "setappdata(gcbf, 'interrupt', true)");
  end  
  while Tpres(n)<=Treq(end) && n<Nt
    kU  = zeros(ord, m.Ndofs);
    kUd = zeros(ord, m.Ndofs);
    Tp = Tpres(n);
    for i=1:ord
%       Tp = Tpres(n)-dt;
      % Option 1: Do not update the "previous states"
      [kU(i, :), kUd(i, :), ~] = ROCFUN(m, Tpres(n)+c(i)*dt, ...
					  U(:, n)+dt*(a(i,:)*kU)', ...
					  Ud(:, n)+dt*(a(i,:)*kUd)', ...
					  Fex(Tpres(n)+c(i)*dt), Tp);
%       % Option 2: Keep updating the "previous states"
%       [kU(i, :), kUd(i, :), m] = ROCFUN(m, Tpres(n)+c(i)*dt, ...
% 					  U(:, n)+dt*(a(i,:)*kU)', ...
% 					  Ud(:, n)+dt*(a(i,:)*kUd)', ...
% 					  Fex(Tpres(n)+c(i)*dt), Tp);
      Tp = Tpres(n) + c(i)*dt;
    end
    en = dt*( norm((b-bs)*kU)+norm((b-bs)*kUd) );
    er = abstol/en;
    if en==0
      en  = dt*eps;
      er  = abstol/en;
      chi = 50;
    else
      chi = er^pow;
    end

    if chi>0.8
      if Treq(n)-Tpres(n)<=dt
        U(:, n+1)  = U(:, n)  + dt*(b*kU)';
        Ud(:, n+1) = Ud(:, n) + dt*(b*kUd)';

        Tpres(n+1) = Tpres(n) + dt;
        [~, Udd(:, n+1), m] = ROCFUN(m, Tpres(n+1), U(:, n+1), Ud(:, n+1), Fex(Tpres(n+1)), Tpres(n));
        
        n = n+1;
      else
        U(:, n)  = U(:, n)  + dt*(b*kU)';
        Ud(:, n) = Ud(:, n) + dt*(b*kUd)';

        Tpres(n) = Tpres(n) + dt;	
        [~, Udd(:, n), m] = ROCFUN(m, Tpres(n), U(:, n), Ud(:, n), Fex(Tpres(n)), Tpres(n)-dt);
      end
    end
    if strcmp(opts.Display, 'progress') || strcmp(opts.Display, 'both')
      fprintf('%d: %.4e/%.4e %.4e rer: %e\n', n, Tpres(n), T1, dt, 1./er);
    end

    if strcmp(opts.Display, 'waitbar')
      waitbar(Tpres(n)/T1, wb, sprintf('Progress: %.e/%.e dt: %.e', Tpres(n), T1, dt));

      if (~ishandle(wb))
        break;
      elseif getappdata(wb, 'interrupt')
        delete(wb);

        U     = U(:, 1:n);
        Ud    = Ud(:, 1:n);
        Udd   = Udd(:, 1:n);
        Tpres = Tpres(1:n);
        return;
      end
    end
    
    if opts.stepadapt
        dt = min(maxstep, chi*dt);
    end
  end
  
  if strcmp(opts.Display, 'waitbar')
    waitbar(1.0, wb, 'COMPLETED!');
    delete(wb);
  end
  
  U     = U(:, 1:n);
  Ud    = Ud(:, 1:n);
  Udd   = Udd(:, 1:n);
  Tpres = Tpres(1:n);
% Command to delete stray waitbars: delete(findall(0,'type','figure','tag','TMWWaitbar'))    
end
function [Ud, Udd, m] = ROCFUN(m, T, U, Ud, Fvec, Tp)
%ROCFUN is the rate of change function
  [Fnl, ~, ~, m] = m.NLFORCE(T, U, Ud, Tp, 0, 0);
  Udd = m.M\(-m.C*Ud-m.K*U-Fnl+Fvec);
end
