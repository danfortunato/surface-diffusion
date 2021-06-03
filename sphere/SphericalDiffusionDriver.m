%% Parameters
params = struct();

% Biological parameters:
params.delta = sqrt(0.05);
params.alpha = 1;
params.beta  = 0.1;
params.nu    = 20;
params.gamma = 0.3;

% Discretization parameters:
params.tend = 50;
params.dt   = 0.1;
params.lmax = 128-1;
params.mmax = params.lmax;
params.nlat = params.lmax+1;
params.nlon = 2*params.mmax+1;
params.plan = sht_plan([params.lmax params.mmax], [params.nlat params.nlon], 'shtns');
params.useHeaviside = false;

% Visualization parameters:
params.quiet    = false;
params.movie    = true;
params.colormap = 'jet';
% params.movfile = VideoWriter('data/sphere.mp4', 'MPEG-4');

%% Initial conditions
init = 'gaussian';

switch lower(init)
    case 'random'
        rng(0)
        U = randnfunsphere(0.1); % Random spherical harmonic expansion
        U = U - min2(U) + 1e-8;  % Restrict to positive values
        U = U ./ max2(U);        % Normalize
    case 'gaussian'
        % Gaussian bump from "Spherical caps in cell polarization"
        U = spherefun(@(x,y,z) exp(-2*(x.^2+y.^2+(z-1).^2)));
    otherwise
        error('Unknown initial condition.');
end

if ( isa(U, 'spherefun') )
    % Get spherical harmonic coefficients from a spherefun
    V = feval(U, params.plan.grid.lon, params.plan.grid.lat);
    U = params.plan.vals2coeffs(V);
end

%% Simulation
tic
U = SphericalDiffusion(U, params);
toc
