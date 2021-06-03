function [f, gamma, eta_c] = analytical_solution(params)
%ANALYTICAL_SOLUTION   Analytical solution to cell polarization on a sphere.
%   F = ANALYTICAL_SOLUTION(PARAMS) computes the analytical solution from
%   [1] for the cell polarization problem on a sphere. PARAMS is a struct
%   containing ALPHA, BETA, DELTA, and either GAMMA or ETA_C. If GAMMA is
%   given, ETA_C is computed by rootfinding on equation (11). If ETA_C is
%   given, GAMMA is computed directly from equation (11). The function
%   handle F represents one longitidunal slice of the axisymmetric
%   solution.
%
%   References:
%
%      [1] R. Diegmiller et al., "Spherical caps in cell polarization",
%          Biophys. J., 115 (2018), pp. 26-30.

alpha = params.alpha;
beta  = params.beta;
delta = params.delta;
gamma = params.gamma;
eta_c = params.eta_c;

mu = (sqrt(1-4/delta^2)-1)/2;
P = @(x) hypergeom([-mu mu+1], 1, (1-x)/2);
dP = matlabFunction(diff(P(sym('x'))));

    function [gamma, k0] = gammaFun(eta_c)
        k0 = 1 + dP(eta_c)*P(-eta_c)/(dP(-eta_c)*P(eta_c));
        gamma = (1 + beta*k0) / (k0 * (1 + alpha*(1+eta_c) + 2*alpha*beta));
    end

if ( ~isempty(eta_c) && isempty(gamma) )

    % Compute gamma from eta_c
    [gamma, k0] = gammaFun(eta_c);

elseif ( isempty(eta_c) && ~isempty(gamma) )

    % Compute eta_c from gamma by root-finding
    % Initial guess:
    eta0 = (1+2*beta-2*gamma-2*alpha*gamma-4*alpha*beta*gamma) / (2*alpha*gamma);
    eta_c = fzero(@(eta) gammaFun(eta)-gamma, eta0);
    [~, k0] = gammaFun(eta_c);

else
    error('Either gamma or eta_c must be specified.');
end

k = gamma * k0 / (1 + beta*k0);

    function out = solutionFun(th)
        eta = -cos(th);
        out = zeros(size(th));
        out(eta<eta_c) = (gamma-(1+beta)*k)*P(-eta(eta<eta_c))/P(-eta_c) + (1+beta)*k;
        out(eta>=eta_c) = (gamma-beta*k)*P(eta(eta>=eta_c))/P(eta_c) + beta*k;
    end

f = @solutionFun;

end
