function F = vals2dbl(F)

% Get the length of the input:
[n,m] = size(F);
F = [[flipud(F(2:n,m/2+1:m)) flipud(F(2:n,1:m/2))];F(1:n,:)];

end
