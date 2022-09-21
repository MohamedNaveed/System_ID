function [alpha_beta] = calculate_markov_from_ABCG(...
                A_hat, B_hat, C_hat, D_hat, G_hat, q,t_steps,nz,nu, true_arma)

alpha_beta = zeros(nz*t_steps, q*(nz + nu) +  nu);

n = size(A_hat,1);

for t = 0:t_steps-1
    
    Phi = eye(n);
    
    C_k = C_hat(:,:,t+1);
    
    if t+1 <= q
        
        i_uplimit = t;
        i_lowlimit = 0;
    else
        i_uplimit = t;
        i_lowlimit = t - q;
    end
    
    for i = i_uplimit : -1 : i_lowlimit
       
        if i == i_uplimit
            
            alpha_beta(t*nz + 1:(t+1)*nz, 1:nu) = D_hat;
        
        elseif i == i_uplimit - 1
            
            B_bar_k_1 = [B_hat(:,:,i+1) + G_hat(:,:,i+1)*D_hat, -G_hat(:,:,i+1)];
            
            h_bar = C_k*B_bar_k_1;
         
            alpha = h_bar(:,nu+1:end);
            beta = h_bar(:,1:nu);
            
            alpha_beta(t*nz + 1:(t+1)*nz, nu + 1:nu+nz) = alpha;
            
            if t<q
                alpha_beta(t*nz + 1:(t+1)*nz, nu + t*nz + 1:nu + t*nz + nu) = beta;
            else
                alpha_beta(t*nz + 1:(t+1)*nz, nu + q*nz + 1:nu + q*nz + nu) = beta;
            end
        elseif i < i_uplimit - 1
            
            A_bar = A_hat(:,:,i+2) + G_hat(:,:,i+2)*C_hat(:,:,i+2);
            
            B_bar = [B_hat(:,:,i+1) + G_hat(:,:,i+1)*D_hat, -G_hat(:,:,i+1)];
            
            Phi = Phi*A_bar;
            
            h_bar = C_k*Phi*B_bar;
            
            alpha = h_bar(:,nu+1:end);
            beta = h_bar(:,1:nu);
            
            alpha_cols = nu + (i_uplimit-i-1)*nz + 1 : nu + (i_uplimit-i)*nz;
           
            alpha_beta(t*nz + 1:(t+1)*nz, alpha_cols) = alpha;
            
            if t<q
                beta_cols = nu + t*nz + (i_uplimit-i-1)*nu + 1 : nu + t*nz + (i_uplimit-i)*nu; 
            else
                beta_cols = nu + q*nz + (i_uplimit-i-1)*nu + 1 : nu + q*nz + (i_uplimit-i)*nu; 
            end
            
            alpha_beta(t*nz + 1:(t+1)*nz, beta_cols) = beta;
        end
            
    end
        
end
end

