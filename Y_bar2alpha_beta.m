function [alpha, beta] = Y_bar2alpha_beta(Y_bar, nu, nz, q)

alpha = zeros(nz, q*nz);
beta = zeros(nz, q*nu);

for i=1:q
    
    beta(:,(i-1)*nu + 1: i*nu) = Y_bar(:,nu + (i-1)*(nu+nz) + 1: nu + (i-1)*(nu+nz) + nu);
    alpha(:,(i-1)*nz + 1: i*nz) = Y_bar(:,nu + (i-1)*(nu+nz) + nu + 1: nu + i*(nu+nz));
    
end


end

