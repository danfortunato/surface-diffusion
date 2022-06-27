function Gamma = interface_length(U, xx, yy, zz, gamma)

% Remap U from DFS coordinates to spherical coordinate domain
vals = real( trigtech.coeffs2vals( trigtech.coeffs2vals(U).' ).' );
u_dbl = chebfun2(vals, [-pi pi -pi pi], 'trig');
u = chebfun2(@(x,y) u_dbl(x,y), [-pi pi 0 pi], 'trigx');

% Same for x, y, z coordinates
vals = real( trigtech.coeffs2vals( trigtech.coeffs2vals(util.vals2coeffsDbl(xx)).' ).' );
x_dbl = chebfun2(vals, [-pi pi -pi pi], 'trig');
x_fun = chebfun2(@(x,y) x_dbl(x,y), [-pi pi 0 pi], 'trigx');

vals = real( trigtech.coeffs2vals( trigtech.coeffs2vals(util.vals2coeffsDbl(yy)).' ).' );
y_dbl = chebfun2(vals, [-pi pi -pi pi], 'trig');
y_fun = chebfun2(@(x,y) y_dbl(x,y), [-pi pi 0 pi], 'trigx');

vals = real( trigtech.coeffs2vals( trigtech.coeffs2vals(util.vals2coeffsDbl(zz)).' ).' );
z_dbl = chebfun2(vals, [-pi pi -pi pi], 'trig');
z_fun = chebfun2(@(x,y) z_dbl(x,y), [-pi pi 0 pi], 'trigz');

% Use chebfun to calculate interfaces where u = gamma
rr = roots(u - gamma);

% Check how many disjoint curves this returned
[~, N_curve] = size(rr);

% Set length to zero
Gamma = 0;

% For loop to find length of each curve
for ci = 1:N_curve

    % Real and imaginary coordinates of curve are spherical coordinates
    rs = real(rr(:, ci));
    rt = imag(rr(:, ci));

    % Smooth out resultant curve to avoid numerical issues
    rs = chebfun(rs, [-1 1], 10);
    rt = chebfun(rt, [-1 1], 10);

    % Extract xyz coordinates from above
    x_curve = x_fun(rs, rt); % Parameterized on [-1 1]
    y_curve = y_fun(rs, rt);
    z_curve = z_fun(rs, rt);

    % Apply another smoothing
    x_curve = chebfun(x_curve,[-1 1], 20);
    y_curve = chebfun(y_curve,[-1 1], 20);
    z_curve = chebfun(z_curve,[-1 1], 20);

    % Calculate tangent vector of each point on curve in xyz space
    dxdt = diff(x_curve);
    dydt = diff(y_curve);
    dzdt = diff(z_curve);

    % Norm of tangent vector, and sum norms to get total curve length
    ds2 = dxdt^2 +dydt^2 + dzdt^2;
    ds_curve = sqrt(ds2);
    ds_curve = chebfun(ds_curve,[-1 1], 40);
    Gamma_loc = sum(ds_curve);
    Gamma = Gamma + Gamma_loc;

end

end
