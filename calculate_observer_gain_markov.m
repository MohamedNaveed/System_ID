function [h_o] = calculate_observer_gain_markov(nu, nz, num_mp, alpha_beta, ...
                            t_steps, q)


h_o = zeros(nz*t_steps, num_mp*nz);

for k = 1:q-1 %time step
    
    for gamma = 1:k
        alpha_k_gamma = -alpha_beta((k)*nz + 1: (k+1)*nz, nu + (gamma-1)*nz + 1:nu + gamma*nz);
        
        temp = alpha_k_gamma;
        
        if (gamma ~=1)
            sum_limit = gamma -1;
            
            for j = 1:sum_limit
               
                alpha_k_j = -alpha_beta((k)*nz + 1: (k+1)*nz, nu + (j-1)*nz + 1:nu + j*nz);
                
                cols = (gamma -j - 1)*nz + 1:(gamma - j)*nz; 
                temp = temp - alpha_k_j*h_o((k-j)*nz + 1: (k+1-j)*nz, cols);
                
            end
            
        end
        cols = (gamma - 1)*nz + 1:(gamma)*nz;
        h_o((k)*nz + 1: (k+1)*nz, cols) = temp;
    end
end

for k = q:t_steps-1
   
    if k>num_mp
        gamma_idx = k-num_mp+1:k;
    else
        gamma_idx = 1:k;
    end
    
    for gamma = gamma_idx
       
        if gamma <=q
            alpha_k_gamma = -alpha_beta((k)*nz + 1: (k+1)*nz, nu + (gamma-1)*nz + 1:nu + gamma*nz);
        else
            alpha_k_gamma = zeros(nz,nz);
        end
        
        temp = alpha_k_gamma;
        
        if (gamma ~= 1)
            if gamma<= q
                sum_limit = gamma-1;
            else
                sum_limit = q;
            end
            
            for j = 1:sum_limit
                
                alpha_k_j = -alpha_beta((k)*nz + 1: (k+1)*nz, nu + (j-1)*nz + 1:nu + j*nz);
                
                cols = (gamma -j - 1)*nz + 1:(gamma - j)*nz; 
                
                temp = temp - alpha_k_j*h_o((k-j)*nz + 1: (k+1-j)*nz, cols);
                
            end
        end
        cols = (gamma - 1)*nz + 1:(gamma)*nz;
        h_o((k)*nz + 1: (k+1)*nz, cols) = temp;
                  
    end
end

end

