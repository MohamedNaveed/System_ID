function [markov_open_loop] = calculate_open_loop_markov_para_ltv(nu, nz,...
                                 num_mp, alpha_beta, ID_time_idxs, t_steps, q)

D_est = alpha_beta((ID_time_idxs(1)-1)*nz + 1: ID_time_idxs(1)*nz,1:nu);

markov_open_loop = zeros(nz*t_steps, num_mp*nu);

markov_open_loop(:,1:nu) = alpha_beta(:,1:nu); %% assigning D

for k = 1:q-1
   k
   for gamma = 1:k
       gamma
       beta_k_gamma = alpha_beta((k)*nz + 1: (k+1)*nz, nu + (k)*nz + (gamma-1)*nu + 1:nu + (k)*nz + gamma*nu);
       alpha_k_gamma = alpha_beta((k)*nz + 1: (k+1)*nz, nu + (gamma-1)*nz + 1:nu + gamma*nz);
       
       temp = beta_k_gamma - alpha_k_gamma*D_est;
       
       if (gamma ~= 1)
           sum_limit = gamma - 1;
           
           for j = 1:sum_limit
                j
                
                alpha_k_j = alpha_beta((k)*nz + 1: (k+1)*nz, nu + (j-1)*nz + 1:nu + j*nz);
                
                mp_cols = (num_mp - (k-gamma) - 1)*nu + 1: (num_mp - (k-gamma))*nu;
                temp = temp -  alpha_k_j*markov_open_loop((k-j)*nz + 1: (k+1-j)*nz, mp_cols);
                
            end
       end
       
       mp_cols = (num_mp - (k-gamma) - 1)*nu + 1: (num_mp - (k-gamma))*nu; %columns for markov parameters. 
       markov_open_loop((k)*nz + 1: (k+1)*nz, mp_cols) = temp;
   end
    
end

for k = q:t_steps-1  % k is the time step at which the ARMA is identified.
    k
    
    if k>num_mp
        gamma_idx = k-num_mp+1:k;
    else
        gamma_idx = 1:k;
    end
    for gamma = gamma_idx
        gamma
        
        if gamma <= q
            beta_k_gamma = alpha_beta((k)*nz + 1: (k+1)*nz, nu + q*nz + (gamma-1)*nu + 1:nu + q*nz + gamma*nu);
            alpha_k_gamma = alpha_beta((k)*nz + 1: (k+1)*nz, nu + (gamma-1)*nz + 1:nu + gamma*nz);
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
  
                alpha_k_j = alpha_beta((k)*nz + 1: (k+1)*nz, nu + (j-1)*nz + 1:nu + j*nz);
               
                mp_cols = (num_mp - (k-gamma) - 1)*nu + 1: (num_mp - (k-gamma))*nu;
                temp = temp -  alpha_k_j*markov_open_loop((k-j)*nz + 1: (k-j + 1)*nz,  mp_cols);
        
                
            end
        end
        mp_cols = (num_mp - (k-gamma) - 1)*nu + 1: (num_mp - (k-gamma))*nu;
        markov_open_loop((k)*nz + 1: (k+1)*nz, mp_cols) = temp;
                
    end
end


end

