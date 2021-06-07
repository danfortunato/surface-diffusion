function [gc, mc] = curvature(varargin)
%CURVATURE   Gaussian and mean curvature of a surface.
%   [GC, MC] = CURVATURE(X, Y, Z) computes the Gaussian curvature GC and
%   mean curvature MC of the surface defined by the Cartesian coordinates
%   (X, Y, Z). All inputs and outputs are given as chebfun2 objects over
%   the domain [-pi pi] x [0 pi].
%
%   [GC, MC] = CURVATURE(RHO, THETA) computes the Gaussian curvature GC and
%   mean curvature MC of the surface of revolution whose generating curve
%   is defined by the polar coordinates (RHO, THETA). RHO and THETA are
%   given as chebfun objects over the domain [0 pi] and GC and MC are
%   returned as chebfun2 objects over the domain [-pi pi] x [0 pi].
%
%   NOTE: If the surface contains singularities at the poles (e.g., if it
%   was defined as a surface of revolution in spherical coordinates), the
%   latter form should be used, as this form takes care to avoid sampling
%   at the poles.

narginchk(2, 3);

if ( nargin == 2 )

    % Polar coordinates
    rho   = varargin{1};
    theta = varargin{2};

    % The curvature of a surface of revolution has a one-dimensional
    % formula
    dom = [-pi pi 0 pi];
    phi = chebfun(@(t) rho(t).*sin(theta(t)), dom(3:4));
    psi = chebfun(@(t) rho(t).*cos(theta(t)), dom(3:4));
    dphi = diff(phi); dphi2 = diff(dphi);
    dpsi = diff(psi); dpsi2 = diff(dpsi);
    g = dphi.^2 + dpsi.^2;
    gc = (-dpsi.^2.*dphi2 + dphi.*dpsi.*dpsi2) ./ (phi.*g.^2);
    mc = (phi.*(dphi2.*dpsi - dphi.*dpsi2) - dpsi.*g) ./ (2.*abs(phi).*g.^(3/2));

    % Evaluate at first-kind Chebyshev points to avoid sampling the poles
    n = max(length(gc), length(mc));
    [~, vv] = chebpts2(n, n, dom, 1);
    pref = chebfunpref(); pref.tech = @chebtech1;
    gc = chebfun2(gc(vv), dom, pref);
    mc = chebfun2(mc(vv), dom, pref);

else

    % Cartesian coordinates
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    r = [x; y; z];

    % First fundamental form
    ru = diff(r, 1, 1);
    rv = diff(r, 1, 2);
    E = dot(ru, ru);
    F = dot(ru, rv);
    G = dot(rv, rv);
    D = E.*G - F.^2;

    % Second fundamental form
    ruu = diff(ru, 1, 1);
    ruv = diff(ru, 1, 2);
    rvv = diff(rv, 1, 2);
    m = -cross(ru, rv) ./ sqrt(D);
    L = dot(ruu, m);
    M = dot(ruv, m);
    N = dot(rvv, m);

    % Curvature
    gc = (L.*N - M.^2) ./ D;
    mc = (E.*N + G.*L - 2*F.*M) ./ (2*D);
end

end
