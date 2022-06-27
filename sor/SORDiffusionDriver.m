%% Geometry
shape = 'prolate';

switch lower(shape)
    case 'sphere'
        a = 1;
        b = 1;
        rho = @(t) 1+0*t;
        theta = @(t) t;
    case 'prolate'
        a = 1.44;
        b = 0.77;
        rho = @(t) a*b./sqrt((a*sin(t)).^2 + (b*cos(t)).^2);
        theta = @(t) t;
    case 'oblate'
        a = 1;
        b = 2;
        rho = @(t) a*b./sqrt((a*sin(t)).^2 + (b*cos(t)).^2);
        theta = @(t) t;
    case 'ufo'
        a = 1.3;
        b = 1.3;
        rho = @(t) 1 + 0.3*cos(4*t);
        theta = @(t) t;
    case 'star'
        a = 1.3;
        b = 1.3;
        rho = @(t) 1 + 0.3*cos(8*t);
        theta = @(t) t;
    case 'bowl'
        a = 0.2;
        b = pi/6;
        sabs = @(x) sqrt(x.^2 + a^2);
        rho = @(t) 1 - a*erf((t-pi/2)/a);
        theta = @(t) b - a + 2*(1-(b-a)/pi)*sabs(t-pi/2);
    case 'egg'
        b = 1;
        d = 0.1;
        a = 2;
        x = @(t) (sqrt(a^2 - d^2.*sin(t).^2)  + d.*cos(t)).*cos(t);
        y = @(t) b.*sin(t);
        rho = @(t) sqrt(x(t).^2 + y(t).^2);
        theta = @(t) acos(x(t)./rho(t));
    case 'hourglass'
        a = 0.97;
        b = 1;
        rho = @(t) 2.*a.^2 .*(cos(2.*t) + sqrt((b/a).^4 - sin(2.*t).^2));
        theta = @(t) t;
    otherwise
        error('Unknown shape.');
end

%% Parameters
params = struct();

% Geometry parameters:
params.rho   = rho;
params.theta = theta;

% Biological parameters:
params.delta = sqrt(0.03);
params.alpha = 3;
params.beta  = 0.1;
params.nu    = 100;
params.gamma = 0.23;

% Discretization parameters:
params.tend = 300;
params.dt   = 0.1;
params.nlat = 128;
params.nlon = 128;
params.nthreads = 8;
params.scheme = 'bdf4'; % or 'bdf1' or 'bdf2'
params.useHeaviside     = false;
params.trueConservation = false;

% Visualization parameters:
params.quiet    = false;
params.movie    = true;
params.colormap = 'jet';
% params.movfile  = VideoWriter('prolate_polarization.mp4', 'MPEG-4');
params.keepAll = false;

%% Initial conditions:
init = 'gaussian';

switch lower(init)
    case 'random'
        rng(0)
        U = randnfun3(0.6, [-b b -b b -a a]);
        U = U - min3(U) + 1e-8;  % Restrict to positive values
        U = U ./ max3(U);        % Normalize
    case 'gaussian'
        U = @(x,y,z) exp(-4*(x.^2+(y+0.7).^2+(z-0.6).^2));
    otherwise
        error('Unknown initial condition.');
end

% Evaluate on the surface
ss = trigpts(params.nlon, [-pi pi]);
tt = linspace(0, pi, params.nlat);
[ss, tt] = meshgrid(ss, tt);
xx = rho(tt).*sin(theta(tt)).*cos(ss);
yy = rho(tt).*sin(theta(tt)).*sin(ss);
zz = rho(tt).*cos(theta(tt));
U = U(xx, yy, zz);

% Convert to double Fourier coefficients
U = util.vals2coeffsDbl(U);

%% Simulation
tic
U = SORDiffusion(U, params);
toc
