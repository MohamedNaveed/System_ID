function [A_hat, B_hat, C_hat] = TVERA(markov_open_loop, q, nu, nz, t_steps)

n_cols = q;


for k = q:t_steps-q-1
    
    Hankel = build_hankel_OL(markov_open_loop,q,k,n_cols, nz, nu);
    
    if k == q
        order = rank(Hankel);
    end
    
    [U,Sig,S] = svd(Hankel);
    O_k = U(:,1:order)*(Sig(1:order,1:order)^(1/2));
    R_k_1 = (Sig(1:order,1:order)^(1/2))*S(:,1:order)';
    
    if k > q
       
        A_hat(:,:,k) = pinv(O_k(1:q*nz,1:order))*O_prev(nz+1:nz+q*nz,1:order);
    
    end
    
    O_prev = O_k;
    
    C_hat(:,:,k+1) = O_k(1:nz,1:order);
    B_hat(:,:,k) = R_k_1(1:order,1:nu);
end

end

