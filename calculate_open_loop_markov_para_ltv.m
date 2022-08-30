function [markov_open_loop] = calculate_open_loop_markov_para_ltv(nu, nz,...
                                 num_mp, alpha_beta, ID_time_idxs, t_steps, q)

D_est = alpha_beta((ID_time_idxs(1)-1)*nz + 1: ID_time_idxs(1)*nz,1:nu);

markov_open_loop = zeros(nz*t_steps, num_mp*nu);

markov_open_loop(:,1:nu) = alpha_beta(:,1:nu);

for k = ID_time_idxs % k is the time step at which the ARMA is identified.
    k
    for gamma = 1:num_mp-1
        gamma
        
        if gamma <= q
            beta_k_gamma = alpha_beta((k-1)*nz + 1: k*nz, nu + q*nz + (gamma-1)*nu + 1:nu + q*nz + gamma*nu);
            alpha_k_gamma = alpha_beta((k-1)*nz + 1: k*nz, nu + (gamma-1)*nz + 1:nu + gamma*nz);
        else
            beta_k_gamma = zeros(nz,nu);
            alpha_k_gamma = zeros(nz,nz);
        end
        
        temp = beta_k_gamma - alpha_k_gamma*D_est;
        
        if (gamma ~= 1)
            if gamma<= q
                sum_limit = gamma-1;
            else
                sum_limit = q;
            end
            
            for j = 1:sum_limit
                j
                
                alpha_k_j = alpha_beta((k-1)*nz + 1: k*nz, nu + (j-1)*nz + 1:nu + j*nz);
                temp = temp -  alpha_k_j*markov_open_loop((k-1-j)*nz + 1: (k-j)*nz, gamma*nu + 1: (gamma + 1)*nu);
                
            end
        end
        
        markov_open_loop((k-1)*nz + 1: k*nz, gamma*nu + 1: (gamma + 1)*nu) = temp;
                
    end
end


end
