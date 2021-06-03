%% Geometry
shape = 'torus';

switch lower(shape)
    case 'torus'
        a = 3; % Major radius
        b = 1; % Minor radius
        r = @(t) a + b*cos(t);
        z = @(t) b*sin(t);
    case 'star-torus'
        a = 3; % Major radius
        b = 1; % Minor radius
        rho = @(t) 0.5*(1 + 0.3*cos(4*t));
        r = @(t) a + b*rho(t).*cos(t);
        z = @(t) b*rho(t).*sin(t);
    otherwise
        error('Unknown shape.');
end

%% Parameters
params = struct();

% Geometry parameters:
params.r = r;
params.z = z;

% Biological parameters:
params.delta = sqrt(0.05);
params.alpha = 0.7;
params.beta  = 0.1;
params.nu    = 100;
params.gamma = 0.17;

% Discretization parameters:
params.nlat = 128;
params.nlon = 256;
params.tend = 800;
params.dt   = 0.05;
params.useHeaviside     = false;
params.trueConservation = true;

% Visualization parameters:
params.quiet    = false;
params.movie    = true;
params.colormap = 'jet';
% params.movfile = VideoWriter('torus.mp4', 'MPEG-4');

%% Initial conditions
init = 'gaussian';

switch lower(init)
    case 'random'
        rng(0)
        U = randnfun3(1, [-a-b a+b -a-b a+b -b b]); % Random smooth function
        U = U - min2(U) + 1e-8;                     % Restrict to positive values
        U = U ./ max2(U);                           % Normalize
    case 'gaussian'
        U = @(x,y,z) exp(-4*((x-a).^2+y.^2+(z-b).^2));
    case 'step'
        U = @(x,y,z) 0.4*(x>4);
    otherwise
        error('Unknown initial condition.');
end

% Evaluate on the surface
ss = trigpts(params.nlon, [0 2*pi]);
tt = trigpts(params.nlat, [0 2*pi]);
[ss, tt] = meshgrid(ss, tt);
xx = params.r(tt).*cos(ss);
yy = params.r(tt).*sin(ss);
zz = params.z(tt) + 0*ss;
U = U(xx, yy, zz);

% Convert to double Fourier coefficients
U = trigtech.vals2coeffs( trigtech.vals2coeffs(U).' ).';

%% Simulation
tic
U = TorusDiffusion(U, params);
toc
