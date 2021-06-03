function coeffs = vals2coeffsDbl(F)

% Get the length of the input:
[n,m] = size(F);
F = [[flipud(F(2:n,m/2+1:m)) flipud(F(2:n,1:m/2))];F(1:n-1,:)];
n = 2*(n-1);

coeffs = (1/n)*fftshift(fft(F, [], 1), 1);

% These coefficients are for interpolation defined on [0,2*pi), but we want
% to work on [-pi,pi). To fix the coefficients for this we just need to
% assign c_k = (-1)^k c_k, for k=-(N-1)/2:(N-1)/2 for N odd and
% k = -N/2:N/2-1 for N even.
if ( mod(n, 2) )
    even_odd_fix = (-1).^(-(n-1)/2:(n-1)/2).';
else
    even_odd_fix = (-1).^((-n/2):(n/2-1)).';
end

coeffs = bsxfun(@times, coeffs, even_odd_fix);

coeffs = (1/m)*fftshift(fft(coeffs, [], 2), 2);

% These coefficients are for interpolation defined on [0,2*pi), but we want
% to work on [-pi,pi). To fix the coefficients for this we just need to
% assign c_k = (-1)^k c_k, for k=-(N-1)/2:(N-1)/2 for N odd and
% k = -N/2:N/2-1 for N even.
if ( mod(m, 2) )
    even_odd_fix = (-1).^(-(m-1)/2:(m-1)/2);
else
    even_odd_fix = (-1).^((-m/2):(m/2-1));
end

coeffs = bsxfun(@times, coeffs, even_odd_fix);

end
