function values = coeffs2valsDbl(coeffs)

% Get the length of the input:
[n,m] = size(coeffs);

% The coefficients are for interpolation defined on [-pi,pi), but the FFT
% works for values on [0,2*pi). To fix the coefficients for this we just need to
% assign c_k = (-1)^k c_k, with k=-(N-1)/2:(N-1)/2 for N odd, and
% k = -N/2:N/2-1 for N even.
if ( mod(n, 2) )
    even_odd_fix = (-1).^(-(n-1)/2:(n-1)/2).';
else
    even_odd_fix = (-1).^((-n/2):(n/2-1)).';
end
coeffs = bsxfun(@times, coeffs, even_odd_fix);

% Shift the coefficients properly.
coeffs = ifft(ifftshift(n*coeffs, 1), [], 1);

if ( mod(m, 2) )
    even_odd_fix = (-1).^(-(m-1)/2:(m-1)/2);
else
    even_odd_fix = (-1).^((-m/2):(m/2-1));
end
coeffs = bsxfun(@times, coeffs, even_odd_fix);

% Shift the coefficients properly.
values = ifft(ifftshift(m*coeffs, 2), [], 2);

values = values([n/2+1:n 1],:);

end
