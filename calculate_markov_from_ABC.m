function [markov_parameters_ABC] = calculate_markov_from_ABC(A_hat, B_hat, C_hat,q,t_steps,num_mp)

nu = size(B_hat,2);
nz = size(C_hat,1);
n = size(A_hat,1);

markov_parameters_ABC = zeros(t_steps*nz,num_mp*nu);

for t = q:t_steps-q-1
    t
    Phi = eye(n);
    
    C_k = C_hat(:,:,t+1);
    
    if t+1 < num_mp
        
        i_uplimit = t;
        i_lowlimit = 0;
    else
        i_uplimit = t;
        i_lowlimit = t - num_mp+1;
    end
    
    for i = i_uplimit:-1:i_lowlimit
        i
        if i == i_uplimit

            markov_parameters_ABC(t*nz + 1:(t+1)*nz, 1:nu) = zeros(nz,nu);

        elseif i == i_uplimit - 1

            B_k_1 = B_hat(:,:,i+1); %B(k-1)
            cols = nu + 1 : nu + nu;
            markov_parameters_ABC(t*nz + 1:(t+1)*nz, cols) = C_k*B_k_1;

        elseif i< i_uplimit - 1

            A = A_hat(:,:,i+2);
            Phi = Phi*A; %output is recorded from t=1
            
            B = B_hat(:,:,i+1);
            cols = nu + (i_uplimit-i-1)*nu + 1 : nu + (i_uplimit-i)*nu;
            markov_parameters_ABC(t*nz + 1:(t+1)*nz,cols) = C_k*Phi*B;
        end
    end
end
end

