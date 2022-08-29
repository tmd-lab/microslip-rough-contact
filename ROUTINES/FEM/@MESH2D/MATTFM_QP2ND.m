function Mx = MATTFM_QP2ND(m, Mx)
%MATTFM_QP2ND transforms matrix given in qp-coordinates to nodal coordinates.
  Mx = m.Tm*Mx*m.Qm;
end
