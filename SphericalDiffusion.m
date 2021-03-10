function [U, dt] = SphericalDiffusion(U, params)
%SPHERICALDIFFUSION   Solve the spherical diffusion problem.
%   U = SPHERICALDIFFUSION(U, PARAMS) solves the PDE
%
%      u_t = delta^2*lap_s(u) - u
%            + (beta + u^nu/(u^nu + gamma^nu))*(1 - 2*pi*alpha*mean2(u))
%
%   on the surface of the sphere, S^2. Here, lap_s is the Laplace-Beltrami
%   operator on S^2 and mean2(u) = 1/(4*pi) int_{S^2} u dS.

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

fprintf('\n')
fprintf('   Spectral size: %d x %d\n', lmax+1, 2*mmax+1)
fprintf('   Grid size:     %d x %d\n', nlat, nlon)
fprintf('   Time step:     %g\n\n', dt)

if ( params.movie )
    % Wrap around longitude for plotting:
    lon = plan.grid.lon([1:end, 1]);
    lat = plan.grid.lat;
    [lon, colat] = meshgrid(lon, lat);
    lat = pi/2 - colat;
    [x, y, z] = sph2cart(lon, lat, 1);

    % Plot the initial condition:
    V = plan.coeffs2vals(U);
    surf(x, y, z, [V V(:,1)])
    shading interp, axis equal
    colormap(params.colormap)
    title(sprintf('t = 0'))
    shg
    
    fprintf('   Press any key when ready.\n\n')
    pause
end

function u = N(u)
    v = plan.coeffs2vals(u);
    v = beta + (v.^nu)./(v.^nu + gamma^nu);
    U2 = plan.vals2coeffs(v);
    scl = 1 - alpha*u(1)/sqrt(pi);
    u = scl*U2 - u;
end

l = util.fromPyramid(repmat((0:lmax).', 1, 2*lmax+1));
t = 0;
for k = 1:nsteps

    rhs = U + dt*N(U);
    U = rhs ./ (1 + dt*delta^2*l.*(l+1));
    t = t + dt;

    if ( params.movie )
        V = plan.coeffs2vals(U);
        surf(x, y, z, [V V(:,1)])
        shading interp, axis equal
        title(sprintf('t = %g', t))
        drawnow
    end
end

end
