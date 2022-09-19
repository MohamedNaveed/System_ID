function [A_hat, B_hat, C_hat, D_hat, G_hat] = TVERA(system, markov_open_loop, q, nu, nz, t_steps, h_o)

n_cols = q*nz;

D_hat = markov_open_loop(1:nz,1:nu);

%% free response experiment to calculate the markov parameters for first q steps. 

[A_hat, B_hat, C_hat, G_hat, O_hat_ref] = free_response_exp(system, q, markov_open_loop, nu, nz, h_o);

%% for k>=q
for k = q:t_steps-q-1
    
    Hankel = build_hankel_OL(markov_open_loop,q,k,n_cols, nz, nu);
    
    if k == q
        order = rank(Hankel);
    end
    
    
    
    [U,Sig,S] = svd(Hankel);
    
    root_Sig = Sig(1:order,1:order)^(1/2);
    O_k = U(:,1:order)*root_Sig;
    R_k_1 = root_Sig*S(:,1:order)';
    
    T_tilde = pinv(O_hat_ref(:,1:order))*O_k(1:q*nz,1:order);
    
    if k > q
       
        A_temp = pinv(O_k(1:q*nz,1:order))*O_prev(nz+1:nz+q*nz,1:order);
        A_hat(:,:,k) = T_tilde*A_temp*inv(T_tilde_prev);
        
        B_temp = R_k_1(1:order,1:nu);
        B_hat(:,:,k) = T_tilde*B_temp;
        
        P_k = calculate_P(h_o, k, q, nz);
        G_hat(:,:,k) = T_tilde*pinv(O_k(1:q*nz,1:order))*P_k;
    end
    
    O_prev = O_k;
    
    C_hat(:,:,k+1) = O_k(1:nz,1:order)*inv(T_tilde);
    
    
    T_tilde_prev = T_tilde;
end

% assuming time invarient in the last q steps.
for k = t_steps-q : t_steps
   
    A_hat(:,:,k) = A_hat(:,:,k-1);
    B_hat(:,:,k) = B_hat(:,:,k-1);
    G_hat(:,:,k) = G_hat(:,:,k-1);
    if k~=t_steps
        C_hat(:,:,k+1) = C_hat(:,:,k);
    end
end

end

