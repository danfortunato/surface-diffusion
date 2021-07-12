classdef LaplaceBeltramiDFS %#ok<*PROPLC>

    properties

        scl   % Scale factor for all ODEs
        ns    % Number of longitudinal Fourier modes
        nt    % Number of latitudinal modes (doubled up)
        nt1   % Number of latitudinal modes
        c     % Helmholtz constant

        x
        y
        z
        normal % Normal vector to surface (a chebfun2v)
        det    % Determinant of the metric tensor (a chebfun2)
        sdet   % Sqrt determinant of the metric tensor (a chebfun2)

        rho    % Radial parametrization in spherical coordinates
        theta  % Polar angle parametrization in spherical coordinates
        dom  = [-pi pi -pi pi]
        dom1 = [-pi pi  0  pi]

    end

    properties ( Hidden )

        ss
        tt
        ss1
        tt1
        bc
        A  % Banded LU decompositions batched over longitudinal modes
        A0 % LU Decomposition for the zero mode (used only when c=0)

    end

    methods

        function L = LaplaceBeltramiDFS(rho, theta, ns, nt, c, nthreads)

            % The operator is lap(u) + c*u
            if ( nargin < 5 )
                L.c = 0;
            else
                L.c = c;
            end

            if ( nargin < 6 )
                nthreads = maxNumCompThreads;
            end

            L.ns  = ns;
            L.nt  = 2*(nt-1);
            L.nt1 = nt;

            L.rho = rho;
            L.theta = theta;
            if ( isa(rho, 'function_handle') )
                rho = chebfun(rho, L.dom1(3:4));
            end
            if ( isa(theta, 'function_handle') )
                theta = chebfun(theta, L.dom1(3:4));
            end

            % Doubled-up grid
            ss = trigpts(L.ns, L.dom(1:2));
            tt = trigpts(L.nt, L.dom(3:4));
            [L.ss, L.tt] = meshgrid(ss, tt);

            % Non-doubled-up grid
            ss1 = trigpts(L.ns, L.dom1(1:2));
            tt1 = linspace(L.dom1(3), L.dom1(4), L.nt1);
            [L.ss1, L.tt1] = meshgrid(ss1, tt1);

            % Indexing origins
            os = floor(L.ns/2) + 1;
            ot = floor(L.nt/2) + 1;

            dr = diff(rho);   drr = diff(dr);
            dt = diff(theta); dtt = diff(dt);
            a00 = (rho.^2.*dt.^2 + dr.^2).^2;
            a01 = sin(theta).*rho.*( cos(theta).*rho.*dt.*(rho.^2.*dt.^2 + dr.^2) + ...
                                     sin(theta).*(dr.^3 - rho.^3.*dt.*dtt - rho.*dr.*drr) );
            a02 = rho.^2.*sin(theta).^2.*(rho.^2.*dt.^2 + dr.^2);
            scl = rho.^2.*sin(theta).^2.*(rho.^2.*dt.^2 + dr.^2).^2;

            % Convert to trig
            a00 = periodize(a00);
            a01 = periodize(a01, 'flip');
            a02 = periodize(a02);
            scl = periodize(scl);

            % Construct the differential operators
            D1   = trigspec.diffmat(L.nt, 1, 1);
            D2   = trigspec.diffmat(L.nt, 2);
            M0   = trigspec.multmat(L.nt, a00);
            M1   = trigspec.multmat(L.nt, a01);
            M2   = trigspec.multmat(L.nt, a02);
            Mscl = trigspec.multmat(L.nt, scl);

            A = cell(L.ns, 1);
            negk = -floor(L.ns/2):-1;
            posk = 1:(floor(L.ns/2)-(mod(L.ns,2)==0));
            prealloc = M2*D2 + M1*D1 - M0 + L.c*Mscl;
            for k = [negk 0 posk]
                A{os+k} = prealloc + (1-k^2)*M0;
            end

            fac = chebfun(@(t) sqrt(rho(t).^2.*dt(t).^2 + dr(t).^2).*rho(t).*sin(theta(t)), [0 pi]);
            mm = -floor(L.nt/2):ceil(L.nt/2)-1;
            trigs = chebfun(@(t) exp(1i*mm.*t), [0 pi]);
            sarea = 2*pi*integral(fac);
            L.bc  = 2*pi*integral(trigs .* fac) ./ sarea;
            if ( L.c == 0 )
                % The k=0 mode has an integral constraint
                % Replace the zero-mode row with this constraint
                A{os}(ot,:) = L.bc;
                L.A0 = decomposition(A{os});
                A{os} = A{1}; % We will eat the cost of one extra inversion
            end
            L.A = BandedBatch(A, nthreads);
            toc

            % We scale the equations by this factor
            L.scl = full(Mscl); % Dense matrix multiplication is faster

            L.x = chebfun2(@(s,t) L.rho(t).*sin(L.theta(t)).*cos(s), L.dom1);
            L.y = chebfun2(@(s,t) L.rho(t).*sin(L.theta(t)).*sin(s), L.dom1);
            L.z = chebfun2(@(s,t) L.rho(t).*cos(L.theta(t)) + 0*s,   L.dom1);
            L.normal = -normal(chebfun2v(L.x, L.y, L.z));

            r = chebfun2v(L.x, L.y, L.z);
            ru = diff(r, 1, 1);
            rv = diff(r, 1, 2);
            E = ru'*ru;
            F = ru'*rv;
            G = rv'*rv;
            D = E.*G - F.^2;
            L.det = D;
            if ( min2(L.det) < 0 )
                L.det = L.det + abs(min2(L.det));
            end
            L.sdet = chebfun2(@(s,t) fac(t)+0*s, L.dom1);

        end

        function u = solve(L, f, const)

            if ( nargin < 3 )
                const = 0;
            end

            % Indexing origins
            os = floor(L.ns/2) + 1;
            ot = floor(L.nt/2) + 1;

            % Get the bivariate Fourier coefficients of the righthand side
            if ( isa(f, 'chebfun2') || isa(f, 'spherefun') || isa(f, 'function_handle') )
                % Underlying discretization grid:
                F = f(L.ss, L.tt);
                F = trigtech.vals2coeffs( trigtech.vals2coeffs( F ).' ).';
            elseif ( isnumeric(f) && isscalar(f) )
                F = zeros(L.nt, L.ns);
                F(ot,os) = f;
            else
                F = f;
            end

            % First, let's project the rhs to have mean zero:
            % meanF = L.bc * F(:,os);
            % F(ot,os) = F(ot,os) - meanF;

            % Scale the RHS by scl
            F = L.scl * F;

            % Do a batch solve
            U = L.A \ F;

            if ( L.c == 0 )
                % The k=0 mode has an integral constraint
                ff = F(:,os);
                ff(ot) = 0;
                U(:,os) = L.A0 \ ff;
            end

            U(ot,os) = U(ot,os) + const;
            u = U;

        end

        function u = mldivide(L, f)
            u = solve(L, f);
        end

        function I = integral2(L, u)
            if ( isa(u, 'chebfun2') )
                I = integral2(u .* L.sdet);
            else
                I = L.bc * u(:,floor(L.ns/2)+1) * surfacearea(L);
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
            A = integral2(L.sdet);
        end

        function V = volume(L)
            F = chebfun2v(L.x, L.y, L.z) / 3;
            V = integral2( F' * L.normal );
        end

        function dom = domain(L)
            dom = L.dom1;
        end

    end

end

function f = periodize(f, flipFlag)

if ( nargin < 2 )
    flipFlag = '';
end

if ( strcmpi(flipFlag, 'flip') )
    f = chebfun({@(t) -f(-t), @(t) f(t)}, [-pi 0 pi]);
else
    f = chebfun({@(t)  f(-t), @(t) f(t)}, [-pi 0 pi]);
end

f = chebfun(@(t) f(t), [-pi pi], 'trig');

end
