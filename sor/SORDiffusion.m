function [U, dt] = SORDiffusion(U, params)
%SORDIFFUSION   Solve the diffusion problem on a surface of revolution.
%   U = SORDIFFUSION(U, PARAMS) solves the PDE
%
%      u_t = delta^2*lap_s(u) - u
%            + (beta + u^nu/(u^nu + gamma^nu))*(1 - 2*pi*alpha*mean2(u))
%
%   on a surface of revolution. Here, lap_s is the Laplace-Beltrami
%   operator on the surface and mean2(u) = (int_{surface} u dS) /
%   (int_{surface} dS). U should be given as a matrix of double Fourier
%   coefficients.
%
%   The parameters are:
%
%      Geometry parameters
%      -------------------
%      PARAMS.RHO:              Radial spherical coordinate of generating curve
%      PARAMS.THETA:            Angular spherical coordinate of generating curve
%
%      Biological parameters
%      ----------------------
%      PARAMS.DELTA:            Diffusivity
%      PARAMS.ALPHA:            Coupling parameter
%      PARAMS.BETA:             Constitutive rate
%      PARAMS.GAMMA:            Membrane recruitment threshold
%      PARAMS.NU:               Cooperativity parameter
%
%      Discretization parameters
%      -------------------------
%      PARAMS.NLAT:             Number of latitudinal coefficients
%      PARAMS.NLON:             Number of longitudinal coefficients
%      PARAMS.DT:               Time step
%      PARAMS.TEND:             End time
%      PARAMS.USEHEAVISIDE:     Replace the nonlinearity with a Heaviside
%      PARAMS.TRUECONSERVATION: Enforce conservation according to the true volume
%
%      Visualization parameters
%      ------------------------
%      PARAMS.QUIET:            Flag to keep quiet
%      PARAMS.MOVIE:            Flag to plot a movie during simulation
%      PARAMS.COLORMAP:         Colormap for plotting
%      PARAMS.MOVFILE:          Filename of movie to output

% Geometry parameters:
rho   = params.rho;
theta = params.theta;

% Biological parameters:
delta = params.delta;
alpha = params.alpha;
beta  = params.beta;
nu    = params.nu;
gamma = params.gamma;

% Discretization parameters:
nlat  = params.nlat;
nlon  = params.nlon;
tend  = params.tend;
dt    = params.dt;

% Make sure dt divides into tend nicely
nsteps = ceil(tend/dt);
dt = tend/nsteps;

% Build the Laplace-Beltrami solver
cscl = -1./(dt*delta^2);
L = LaplaceBeltramiDFS(rho, theta, nlon, nlat, cscl);

if ( ~params.quiet )
    fprintf('\n')
    fprintf('   Spectral size: %d x %d\n', nlat, nlon)
    fprintf('   Grid size:     %d x %d\n', nlat, nlon)
    fprintf('   Time step:     %g\n\n', dt)
end

outputMovie = isfield(params, 'movfile') && ~isempty(params.movfile);
if ( outputMovie )
    open(params.movfile);
end

if ( params.movie )
    % Wrap around longitude for plotting:
    ss = wrap(L.ss1);
    tt = wrap(L.tt1);
    xx = L.x(ss,tt);
    yy = L.y(ss,tt);
    zz = L.z(ss,tt);

    % Plot the initial condition:
    V = util.coeffs2valsDbl(U);
    V = real(V);
    V = wrap(V);
    surf(xx, yy, zz, V)
    colormap(params.colormap)
    setupfigure(0)
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

% The true conservation law changes depending on the volume of the
% geometry, which has the effect of rescaling alpha.
if ( params.trueConservation )
    % Conservation based on the true volume
    fac = 3/2 * volume(L);
else
    % Conservation based on the volume of a sphere
    fac = 2*pi;
end

    function u = N_direct(u)
        v = util.coeffs2valsDbl(u);
        v = beta + (v.^nu)./(v.^nu + gamma^nu);
        U2 = util.vals2coeffsDbl(v);
        scl = 1 - alpha/fac * L.integral2(u);
        u = scl*U2 - u;
    end

    function u = N_heaviside(u)
        v = util.coeffs2valsDbl(u);
        mu = params.nlat/(5*log(params.nlat));
        v = beta + 1./(1+exp(-2*mu*(v-gamma)));
        U2 = util.vals2coeffsDbl(v);
        scl = 1 - alpha/fac * L.integral2(u);
        u = scl*U2 - u;
    end

if ( params.useHeaviside )
    N = @N_heaviside;
else
    N = @N_direct;
end

t = 0;
for k = 1:nsteps

    rhs = U + dt*N(U);
    U = L \ (cscl * rhs);
    t = t + dt;

    if ( params.movie )
        V = util.coeffs2valsDbl(U);
        V = real(V);
        V = wrap(V);
        surf(xx, yy, zz, V)
        setupfigure(t)
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

function setupfigure(t)
    shading interp
    axis equal off
    title(sprintf('$t = %.2f$', t), 'Interpreter', 'Latex', 'FontSize', 18)
    colorbar('FontSize', 14)
end
