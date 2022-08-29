function Uqp = INTERP_QP(m, U, varargin)
%INTERP_QP interpolates the given nodal data to quadrature points
  if nargin==2
    Uqp = m.Qm*double(U);
  else
    Uqp = m.Qm*double(varargin{1}*U);
  end
end
