function [U, dt] = SphericalDiffusion(U, params)
%SPHERICALDIFFUSION   Solve the spherical diffusion problem.
%   U = SPHERICALDIFFUSION(U, PARAMS) solves the PDE
%
%      u_t = delta^2*lap_s(u) - u
%            + (beta + u^nu/(u^nu + gamma^nu))*(1 - 2*pi*alpha*mean2(u))
%
%   on the surface of the sphere, S^2. Here, lap_s is the Laplace-Beltrami
%   operator on S^2 and mean2(u) = 1/(4*pi) int_{S^2} u dS. U should be
%   given as a vector of spherical harmonic coefficients.
%
%   The parameters are:
%
%      Biological parameters
%      ----------------------
%      PARAMS.DELTA:        Diffusivity
%      PARAMS.ALPHA:        Coupling parameter
%      PARAMS.BETA:         Constitutive rate
%      PARAMS.GAMMA:        Membrane recruitment threshold
%      PARAMS.NU:           Cooperativity parameter
%
%      Discretization parameters
%      -------------------------
%      PARAMS.LMAX:         Maximum spherical harmonic degree
%      PARAMS.MMAX:         Maximum spherical harmonic order
%      PARAMS.NLAT:         Number of latitudinal grid points
%      PARAMS.NLON:         Number of longitudinal grid points
%      PARAMS.PLAN:         Spherical harmonic transform plan
%      PARAMS.DT:           Time step
%      PARAMS.TEND:         End time
%      PARAMS.USEHEAVISIDE: Replace the nonlinearity with a Heaviside
%
%      Visualization parameters
%      ------------------------
%      PARAMS.QUIET:        Flag to keep quiet
%      PARAMS.MOVIE:        Flag to plot a movie during simulation
%      PARAMS.COLORMAP:     Colormap for plotting
%      PARAMS.MOVFILE:      Filename of movie to output

% Biological parameters:
delta = params.delta;
alpha = params.alpha;
beta  = params.beta;
nu    = params.nu;
gamma = params.gamma;

% Discretization parameters:
tend = params.tend;
dt   = params.dt;
lmax = params.lmax;
mmax = params.mmax;
nlat = params.nlat;
nlon = params.nlon;
plan = params.plan;
assert(nlat >= lmax+1, 'Bad grid size.');

% Make sure dt divides into tend nicely
nsteps = ceil(tend/dt);
dt = tend/nsteps;

if ( ~params.quiet )
    fprintf('\n')
    fprintf('   Spectral size: %d x %d\n', lmax+1, 2*mmax+1)
    fprintf('   Grid size:     %d x %d\n', nlat, nlon)
    fprintf('   Time step:     %g\n\n', dt)
end

outputMovie = isfield(params, 'movfile') && ~isempty(params.movfile);
if ( outputMovie )
    open(params.movfile);
end

if ( params.movie )
    % Wrap around longitude for plotting:
    [lon, colat] = meshgrid(plan.grid.lon, plan.grid.lat);
    lat = pi/2 - colat;
    lon = wrap(lon);
    lat = wrap(lat);
    [x, y, z] = sph2cart(lon, lat, 1);

    % Plot the initial condition:
    V = plan.coeffs2vals(U);
    V = wrap(V);
    surf(x, y, z, V)
    shading interp, axis equal
    colormap(params.colormap)
    title('$t = 0$', 'interpreter', 'latex')
    colorbar
    shg

    if ( outputMovie )
        frame = getframe(gcf);
        writeVideo(params.movfile, frame);
    end

    if ( ~params.quiet )
        fprintf('   Press any key when ready.\n\n')
        pause
    end
end

    function u = N_direct(u)
        v = plan.coeffs2vals(u);
        v = beta + (v.^nu)./(v.^nu + gamma^nu);
        U2 = plan.vals2coeffs(v);
        scl = 1 - alpha*u(1)/sqrt(pi);
        u = scl*U2 - u;
    end

    function u = N_heaviside(u)
        v = plan.coeffs2vals(u);
        mu = params.nlat/(5*log(params.nlat));
        v = beta + 1./(1+exp(-2*mu*(v-gamma)));
        U2 = plan.vals2coeffs(v);
        scl = 1 - alpha*u(1)/sqrt(pi);
        u = scl*U2 - u;
    end

if ( params.useHeaviside )
    N = @N_heaviside;
else
    N = @N_direct;
end

l = util.fromPyramid(repmat((0:lmax).', 1, 2*lmax+1));
t = 0;
for k = 1:nsteps

    rhs = U + dt*N(U);
    U = rhs ./ (1 + dt*delta^2*l.*(l+1));
    t = t + dt;

    if ( params.movie )
        V = plan.coeffs2vals(U);
        V = wrap(V);
        surf(x, y, z, V)
        shading interp, axis equal
        title(sprintf('$t = %g$', t), 'interpreter', 'latex')
        colorbar
        drawnow

        if ( outputMovie )
            frame = getframe(gcf);
            writeVideo(params.movfile, frame);
        end
    end
end

if ( outputMovie )
    close(params.movfile);
end

end

function X = wrap(X)

X = [X X(:,1)];

end
