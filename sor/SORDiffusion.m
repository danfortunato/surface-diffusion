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
%      PARAMS.SCHEME:           Time-stepping scheme ('bdf1', 'bdf2', or 'bdf4')
%      PARAMS.USEHEAVISIDE:     Replace the nonlinearity with a Heaviside
%      PARAMS.TRUECONSERVATION: Enforce conservation according to the true volume
%
%      Visualization parameters
%      ------------------------
%      PARAMS.QUIET:            Flag to keep quiet
%      PARAMS.MOVIE:            Flag to plot a movie during simulation
%      PARAMS.COLORMAP:         Colormap for plotting
%      PARAMS.MOVFILE:          Filename of movie to output
%      PARAMS.KEEPALL:          Flag to return an array of solutions at all time points.
%                               If true, U is returned as a tensor of DFS coefficients.
%                               U(:,:,k) is the solution at time step k-1.

% Geometry parameters:
rho   = params.rho;
theta = params.theta;

% Biological parameters:
delta = params.delta;
alpha = params.alpha;
beta  = params.beta;
nu    = params.nu;
gamma = params.gamma;
gamma_nu = gamma^nu;

% Discretization parameters:
nlat  = params.nlat;
nlon  = params.nlon;
tend  = params.tend;
dt    = params.dt;

% Make sure dt divides into tend nicely
nsteps = ceil(tend/dt);
dt = tend/nsteps;

% Build the Laplace-Beltrami solver
switch ( lower(params.scheme) )
    case 'bdf1'
        cscl = -1/(dt*delta^2);
        step = @bdf1;
        multistep = 1;
    case 'bdf2'
        cscl = -3/(2*dt*delta^2);
        step = @bdf2;
        multistep = 2;
    case 'bdf3'
        cscl = -11/(6*dt*delta^2);
        step = @bdf3;
        multistep = 3;
    case 'bdf4'
        cscl = -25/(12*dt*delta^2);
        step = @bdf4;
        multistep = 4;
    otherwise
        error('Unknown time-stepping scheme.');
end
L = LaplaceBeltramiDFS(rho, theta, nlon, nlat, cscl, params.nthreads);

keepAll = false;
if ( isfield(params, 'keepAll') )
    keepAll = params.keepAll;
end

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
        v_nu = v.^nu;
        v = beta + v_nu./(v_nu + gamma_nu);
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

if ( keepAll )
    UU = zeros([size(U) nsteps]);
    UU(:,:,1) = U;
end

t = 0;
Uold = cell(multistep-1, 1);
if ( multistep > 1 )
    % Do a few steps of LIRK4 to start the multistep method
    cscl1 = -4/(dt*delta^2);
    L1 = LaplaceBeltramiDFS(rho, theta, nlon, nlat, cscl1, params.nthreads);
    for k = 1:multistep-1
        Uold{multistep-k} = U;
        rhs = U + dt*1/4*N(U);
        a = L1 \ (cscl1 * rhs);
        rhs = U + dt*(delta^2*L1.lap(1/2*a) - 1/4*N(U) + N(a));
        b = L1 \ (cscl1 * rhs);
        rhs = U + dt*(delta^2*L1.lap(17/50*a - 1/25*b) - 13/100*N(U) ...
            + 43/75*N(a) + 8/75*N(b));
        c = L1 \ (cscl1 * rhs);
        rhs = U + dt*(delta^2*L1.lap(371/1360*a - 137/2720*b + 15/544*c) ...
            - 6/85*N(U) + 42/85*N(a) + 179/1360*N(b) - 15/272*N(c));
        d = L1 \ (cscl1 * rhs);
        rhs = U + dt*(delta^2*L1.lap(25/24*a - 49/48*b + 125/16*c - 85/12*d) ...
            + 79/24*N(a) - 5/8*N(b) + 25/2*N(c) - 85/6*N(d));
        e = L1 \ (cscl1 * rhs);
        U = U + dt*(delta^2*L1.lap(25/24*a - 49/48*b + 125/16*c - 85/12*d + 1/4*e) ...
            + 25/24*N(a) - 49/48*N(b) + 125/16*N(c) - 85/12*N(d) + 1/4*N(e));
        t = t + dt;
        nsteps = nsteps - 1;

        if ( keepAll )
            UU(:,:,k+1) = U;
        end

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
end

for k = 1:nsteps

    [rhs, Uold] = step(N, dt, U, Uold);
    U = L \ (cscl * rhs);
    t = t + dt;

    if ( keepAll )
        UU(:,:,k+multistep) = U;
    end

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

if ( keepAll )
    U = UU;
end

end

%% Time-stepping schemes

function [rhs, Uold] = bdf1(N, dt, U, Uold)
    rhs = U + dt*N(U);
end

function [rhs, Uold] = bdf2(N, dt, U, Uold)
    rhs = 4/3*U - 1/3*Uold{1} + dt*(4/3*N(U) - 2/3*N(Uold{1}));
    Uold{1} = U;
end

function [rhs, Uold] = bdf3(N, dt, U, Uold)
    rhs = 18/11*U - 9/11*Uold{1} + 2/11*Uold{2} + ...
        dt*(18/11*N(U) - 18/11*N(Uold{1}) + 6/11*N(Uold{2}));
    Uold(2) = Uold(1);
    Uold{1} = U;
end

function [rhs, Uold] = bdf4(N, dt, U, Uold)
    rhs = 48/25*U - 36/25*Uold{1} + 16/25*Uold{2} - 3/25*Uold{3} + ...
          dt*(48/25*N(U) - 72/25*N(Uold{1}) + 48/25*N(Uold{2}) - 12/25*N(Uold{3}));
    Uold(2:3) = Uold(1:2);
    Uold{1} = U;
end

%% Plotting utilities

function X = wrap(X)
    X = [X X(:,1)];
end

function setupfigure(t)
    shading interp
    axis equal off
    title(sprintf('$t = %.2f$', t), 'Interpreter', 'Latex', 'FontSize', 18)
    colorbar('FontSize', 14)
end
