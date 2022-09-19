function [H] = build_hankel_OL(markov_open_loop,q,k, n_cols, nz, nu)

% n_cols is the number of times to take the column. Doesn't mean the exact
% number of columns. 

H = zeros(q*nz, n_cols*nu); %hankel

%{
if (k-n_cols < 0)
    disp('request invalid');
end
%}

for i = 0:q %rows of hankel
    
    rows = (k+i)*nz + 1: (k+i+1)*nz;
    cols = nu + i*nu + 1: nu + (i+n_cols)*nu;
    H(i*nz+1:(i+1)*nz,:) = markov_open_loop(rows,cols);
end

