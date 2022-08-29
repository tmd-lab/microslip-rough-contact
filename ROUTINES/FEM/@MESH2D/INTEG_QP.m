function Un = INTEG_QP(m, U, varargin)
%INTEG_QP integrates given qp data to nodal points
  if nargin==2
    Un = m.Tm*U;
  else
    Un = varargin{1}*m.Tm*U;
  end
end
