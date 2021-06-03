classdef LaplaceBeltrami %#ok<*PROPLC>

    properties

        A     % Cell array of discretized ODOs for each longitudinal mode
        scl   % Scale factor for all ODEs
        ns    % Number of longitudinal Fourier modes
        nt    % Number of latitudinal modes
        c     % Helmholtz constant

        x
        y
        z
        normal % Normal vector to surface (a chebfun2v)
        det    % Determinant of the metric tensor (a chebfun2)

        rf     % Radial parametrization in cylindrical coordinates
        zf     % Vertical parametrization in cylindrical coordinates
        dom = [0 2*pi 0 2*pi]

    end

    properties ( Hidden )

        ss
        tt
        ops
        bc

    end

    methods

        function L = LaplaceBeltrami(r, z, ns, nt, c)

            % The operator is lap(u) + c*u
            if ( nargin < 5 )
                L.c = 0;
            else
                L.c = c;
            end

            L.ns = ns;
            L.nt = nt;

            ss = trigpts(ns, L.dom(1:2));
            tt = trigpts(nt, L.dom(3:4));
            [L.ss, L.tt] = meshgrid(ss, tt);

            % Indexing origins
            os = floor(ns/2) + 1;
            ot = floor(nt/2) + 1;

            L.rf = r;
            L.zf = z;

            if ( isa(r, 'chebfun') )
                if ( ~isPeriodicTech(z) )
                    r = chebfun(@(t) r(t), domain(r), 'trig');
                end
            elseif ( isa(r, 'function_handle') )
                r = chebfun(r, L.dom(3:4), 'trig');
            else
                error('Expected r(t) as a function handle or chebfun.');
            end

            if ( isa(z, 'chebfun') )
                if ( ~isPeriodicTech(z) )
                    z = chebfun(@(t) z(t), domain(z), 'trig');
                end
            elseif ( isa(z, 'function_handle') )
                z = chebfun(z, L.dom(3:4), 'trig');
            else
                error('Expected z(t) as a function handle or chebfun.');
            end

            if ( ~all(domain(r) == domain(z)) )
                error('Inconsistent domains.');
            end
            L.dom(3:4) = domain(r);

            dr = diff(r); drr = diff(dr);
            dz = diff(z); dzz = diff(dz);

            I   = speye(nt);
            D2  = trigspec.diffmat(nt, 2);
            D1  = trigspec.diffmat(nt, 1, 1);
            M2  = trigspec.multmat(nt, r.^2./(dr.^2+dz.^2));
            M1  = trigspec.multmat(nt, (r.*dr-r.^2.*(dr.*drr+dz.*dzz))./(dr.^2+dz.^2));
            Mr2 = trigspec.multmat(nt, r.^2);

            L.A = cell(ns, 1);
            negk = -floor(L.ns/2):-1;
            posk = 1:(floor(L.ns/2)-(mod(L.ns,2)==0));
            for k = [negk posk]
                L.A{os+k} = M2*D2 + M1*D1 - k^2*I + L.c*Mr2;
                %L.A{os+k} = decomposition(L.A{os+k});
            end

            A = M2*D2 + M1*D1 + L.c*Mr2;
            L.bc = diff(L.dom(1:2)) * diff(L.dom(3:4)) * trigcoeffs(r, nt).';
            if ( L.c == 0 )
                % The k=0 mode has an integral constraint
                % Replace the zero-mode row with this constraint
                A(ot,:) = L.bc;
            end
            L.A{os} = decomposition(A);

            % We scale the equations by this factor
            L.scl = Mr2;

            L.x = chebfun2(@(s,t) L.rf(t).*cos(s), L.dom, 'trig');
            L.y = chebfun2(@(s,t) L.rf(t).*sin(s), L.dom, 'trig');
            L.z = chebfun2(@(s,t) L.zf(t) + 0*s,   L.dom, 'trig');
            [L.ops.grad, L.ops.div, L.ops.curl, L.ops.lap, L.det] = util.diffs(L.x, L.y, L.z);
            L.normal = normal(chebfun2v(L.x, L.y, L.z));

        end

        function u = solve(L, f, const)

            if ( nargin < 3 )
                const = 0;
            end

            % Indexing origins
            os = floor(L.ns/2) + 1;
            ot = floor(L.nt/2) + 1;

            % Get the bivariate Fourier coefficients of the righthand side
            if ( isa(f, 'chebfun2') || isa(f, 'function_handle') )
                F = f(L.ss, L.tt);
                F = trigtech.vals2coeffs( trigtech.vals2coeffs( F ).' ).';
            elseif ( isnumeric(f) && isscalar(f) )
                F = zeros(L.nt, L.ns);
                F(ot,os) = f;
            else
                F = f;
            end

            U = zeros(L.nt, L.ns);
            negk = -floor(L.ns/2):-1;
            posk = 1:(floor(L.ns/2)-(mod(L.ns,2)==0));
            for k = [negk posk]
                U(:,os+k) = L.A{os+k} \ (L.scl * F(:,os+k));
            end

            ff = L.scl * F(:,os);
            if ( L.c == 0 )
                % The k=0 mode has an integral constraint
                ff(ot) = 0;
            end
            U(:,os) = L.A{os} \ ff;
            U(ot,os) = U(ot,os) + const;

            if ( isa(f, 'chebfun2') )
                % Convert to chebfun2 if we were given a chebfun2
                u = chebfun2(U, L.dom, 'trig', 'coeffs');
                u = real(u);
            else
                u = U;
            end

        end

        function u = mldivide(L, f)
            u = solve(L, f);
        end

        function f = mtimes(L, u)
            f = lap(L, u) + L.c * u;
        end

        function I = integral2(L, u)
            if ( isa(u, 'chebfun2') )
                I = integral2(u .* sqrt(L.det));
            else
                I = L.bc * u(:,floor(L.ns/2)+1);
                I = real(I);
            end
        end

        function I = integral(L, u)
            I = integral2(L, u);
        end

        function I = sum2(L, u)
            I = integral2(L, u);
        end

        function I = mean2(L, u)
            I = integral2(L, u) / surfacearea(L);
        end

        function A = surfacearea(L)
            A = integral2(sqrt(L.det));
        end

        function V = volume(L)
            F = chebfun2v(L.x, L.y, L.z) / 3;
            V = integral2( F' * L.normal );
        end

        function u = grad(L, u)
            u = L.ops.grad(u);
        end

        function u = div(L, u)
            u = L.ops.div(u);
        end

        function u = curl(L, u)
            u = L.ops.curl(u);
        end

        function u = lap(L, u)
            u = L.ops.lap(u);
        end

        function dom = domain(L)
            dom = L.dom;
        end

    end

end
