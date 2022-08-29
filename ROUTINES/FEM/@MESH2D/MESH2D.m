classdef MESH2D
  %MESH2D creates a MESH class.

  properties
    Nds		% Nodal locations [x y]
    Tri		% Triangular elements [eid n1 n2 n3]
    Quad	% Quadrilateral elements [eid n1 n2 n3 n4]
    Nn		% Number of nodes
    Ne_Tri	% Number triangular elements
    Ne_Quad	% Number of quadrilateral elements
    Ne		% Number of elements

    dpn		% Number of dofs per node

    Nq=2	% Quadrature order
    Qm 		% Quadrature point interpolation matrix
    Tm		% Quadrature point integration matrix

    fcont 	% Contact Function
    z       % Slider States for fcont
  end

  methods
    function m = MESH2D(nds, dpn, tri, quad, nq)
    %MESH2D initializes the 2D Mesh object.
    %
    % USAGE:
    %  m = MESH2D(Nds, Tri, Quad);
      % Nodes & Dofs per node
      m.Nds = nds;
      m.Nn  = size(nds, 1);
      m.dpn = dpn;
      
      % Triangular Elements
      m.Tri = tri;
      if ~isempty(tri)  % Non-empty triangular elements: conduct checks
	if size(tri,2) ~= 4 || max(max(tri(:, 2:end))) > m.Nn
	  error('Triangular elements error');
	end
      end
      m.Ne_Tri = size(tri, 1);

      % Quadrilateral Elements
      m.Quad = quad;
      if ~isempty(quad)  % Non-empty quad elements: conduct checks
	if size(quad,2) ~= 5 || max(max(quad(:, 2:end))) > m.Nn
	  error('Quadrilateral elements error');
	end
      end
      m.Ne_Quad = size(quad, 1);
      
      m.Ne = m.Ne_Tri+m.Ne_Quad;

      m = m.SETQUAD(nq);
    end

    function m = SETQUAD(m, nq)
      m.Nq = nq;

      % Generate Quadrature points and weights for the reference Quad
      [m.Qm, m.Tm] = m.ND2QP();
    end

    function m = SETCFUN(m, fcont, varargin)
      m.fcont = fcont;
      if length(varargin)==1
          m.z = varargin{1};
      end
    end

    function m = SINGLEPREC(m)
      m.Nds    = single(m.Nds);
      m.Tri    = uint16(m.Tri);
      m.Quad   = uint16(m.Quad);
      m.Nn     = uint16(m.Nn);
      m.Ne_Tri = uint16(m.Ne_Tri);
      m.Ne_Quad= uint16(m.Ne_Quad);
      m.Ne     = uint16(m.Ne);

      m.dpn = uint16(m.dpn);
      
      m.Nq = uint16(m.Nq);
      m.Qm = m.Qm;
      m.Tm = m.Tm;
      
      z = single(m.z);
    end
  end
end
