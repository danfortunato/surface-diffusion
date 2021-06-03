function [grad_s, div_s, curl_s, lap_s, D] = diffs(x, y, z)
%DIFFS   Differential operators on surfaces.
%   [GRAD, DIV, CURL, LAP] = DIFFS(X, Y, Z) returns the surface gradient,
%   divergence, curl, and Laplacian operators over the parametric surface
%   defined by the chebfun2 objects (X, Y, Z). GRAD, DIV, CURL, and LAP are
%   function handles which accept as input chebfun2 and chebfun2v objects.

r = [x; y; z];
ru = diff(r, 1, 1);
rv = diff(r, 1, 2);
E = ru'*ru;
F = ru'*rv;
G = rv'*rv;
D = E.*G - F.^2;
if ( min2(D) < 0 )
    D = D + abs(min2(D)) + 1e-14;
end
r1 = (G*ru-F*rv)/D;
r2 = (E*rv-F*ru)/D;

    function out = GRAD(f)
        if ( ~isa(f, 'chebfun2') )
            error('f must be a chebfun2.')
        end
        out = r1.*diff(f,1,1) + r2.*diff(f,1,2);
    end

    function out = DIV(f)
        if ( ~isa(f, 'chebfun2v') )
            error('f must be a chebfun2v.')
        end
        out = dot(diff(f,1,1), r1) + dot(diff(f,1,2), r2);
    end

    function out = CURL(f)
        if ( ~isa(f, 'chebfun2v') )
            error('f must be a chebfun2v.')
        end
        gf1 = GRAD(f(1)); gf2 = GRAD(f(2)); gf3 = GRAD(f(3));
        out = [gf3(2)-gf2(3); gf1(3)-gf3(1); gf2(1)-gf1(2)];
    end

    function out = LAP(f)
        if ( ~isa(f, 'chebfun2') )
            error('f must be a chebfun2.')
        end
        out = DIV(GRAD(f));
    end

grad_s = @GRAD;
div_s  = @DIV;
curl_s = @CURL;
lap_s  = @LAP;

end
