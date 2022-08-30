function [H] = build_hankel(alpha_beta,q,k,n_cols, nz, nu)

% n_cols is the number of times to take the column. Doesn't mean the exact
% number of columns. 

H = zeros(q*nz, n_cols*(nz + nu)); %hankel

if (k-n_cols <= 0)
    disp('request invalid');
end

for j = 1:n_cols
    for i = 0:q-1
        i
        if ((k+i) <= q)

            rows = (k + i - 1)*nz + 1: (k + i)*nz;

            a_cols = nu + ((k-2) + i + j-1)*nz + 1:nu + ((k-1) + i + j-1)*nz;

            alpha = alpha_beta(rows, a_cols);

            b_cols = nu + (k + i -1)*nz + ((k-2) + i + (j-1))*nu + 1:...
                    nu + (k + i -1)*nz + ((k-1) + i + (j-1))*nu;

            beta = alpha_beta(rows, b_cols);

        else
            
            rows = (k + i - 1)*nz + 1: (k + i)*nz;

            a_cols = nu + (i + j-1)*nz + 1:nu + ((i + 1) + j-1)*nz;

            alpha = alpha_beta(rows, a_cols);

            b_cols = nu + q*nz + (i + (j-1))*nu + 1:...
                    nu + q*nz + ((i+1) + (j-1))*nu;

            beta = alpha_beta(rows, b_cols);
        end
        H(i*nz + 1: (i+1)*nz,(j-1)*(nu+nz) + 1:j*(nu+nz)) = [alpha, beta];
    end
    
end
end

