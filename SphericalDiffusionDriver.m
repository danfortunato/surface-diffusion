% Biological parameters:
params.delta = sqrt(0.05);
params.alpha = 1;
params.beta  = 0.1;
params.nu    = 10;
params.gamma = 0.3;

% Discretization parameters:
params.tend = 10;
params.dt   = 0.1;
params.lmax = 128-1;
params.mmax = params.lmax;
params.nlat = params.lmax+1;
params.nlon = 2*params.mmax+1;
params.plan = sht_plan([params.lmax params.mmax], [params.nlat params.nlon], 'shtns');

% Visualization parameters:
params.movie    = true;
params.colormap = 'jet';

% Initial conditions:
% Random spherical harmonic expansion
U = randnfunsphere(0.1);
% Gaussian bump from "Spherical Caps in Cell Polarization"
% U = spherefun(@(x,y,z) exp(-2*(x.^2+y.^2+(z-1).^2)));

if ( isa(U, 'spherefun') )
    % Get spherical harmonic coefficients from a spherefun
    V = feval(U, params.plan.grid.lon, params.plan.grid.lat);
    U = params.plan.vals2coeffs(V);
end

tic
U = SphericalDiffusion(U, params);
toc
