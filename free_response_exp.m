function [A_hat_tilde, B_hat_tilde, C_hat_tilde, G_hat_tilde, O_hat_ref] = ...
    free_response_exp(system, q, markov_open_loop, nu, nz, h_o)

if strcmp(system,'oscillator')
    sysd = oscillator(0);
    
    rng(0);
    
    n = size(sysd.A,1); % order of the system
    
    N = q+10; %number of experiments

    X0 = normrnd(0, 20, n, N);

    t_steps = q+10;
    y_matrix = zeros(nz*t_steps, N);

    n_cols = 1; %number of time steps in Hankel columns
    u_vec = zeros(nu, t_steps);

    for i=1:N

            x0 = X0(:,i); 
            if strcmp(system,'oscillator')
                [x, y] = generate_response_oscillator(x0, u_vec, n, nz, sysd.Ts);%output
            end

            y_matrix(:,i) = reshape(y, nz*t_steps,1); 

    end


elseif strcmp(system, 'cartpole')
   
    load('data/cartpole_data_free_response_sept21.mat');
    n = Model.nsys; %states in the system
    U = flipud(delta_u);
    
    n_cols = 1; %number of time steps in Hankel columns
    t_steps = Task.horizon+1;
 
    y_matrix = flipud(delta_y(1:nz*t_steps,:));

end

%calculate A,B,C, G
for k= 0:q
    
    
    [U,Sig,V] = svd(y_matrix((k)*nz + 1: (k+q)*nz,:));

    rank_Sig = rank(Sig,1e-5);

    root_Sig = Sig(1:rank_Sig,1:rank_Sig)^(1/2);

    O_hat = U(:,1:rank_Sig)*root_Sig;
        
    X_hat = root_Sig*V(:,1:rank_Sig)';
    
    
    if k~=0
        A_hat(:,:,k) = X_hat*pinv(X_hat_prev); %A calculation
        
        Hankel = build_hankel_OL(markov_open_loop,q, k, n_cols, nz, nu);
        
        B_hat(:,:,k) = pinv(O_hat)*Hankel(1:q*nz,:);
        
        P_k = calculate_P(h_o, k, q, nz);
        
        G_hat(:,:,k) = pinv(O_hat)*P_k;
        
    end

    if k<=q
        C_hat(:,:,k+1) = O_hat(1:nz,:);
        D_hat_tilde(:,:,k+1) = markov_open_loop((k)*nz + 1 : (k+1)*nz, 1:nu);
    end
    X_hat_prev = X_hat;
    
    % transforming to reference coordinates. 
    
    if k==0
        
        C_hat_tilde(:,:,k+1) = C_hat(:,:,k+1);
        O_hat_ref = O_hat;
       
    elseif k==1
        
        T_tilde = pinv(O_hat_ref)*O_hat;
        C_hat_tilde(:,:,k+1) = C_hat(:,:,k+1)*inv(T_tilde);
        
        
        B_hat_tilde(:,:,k) = T_tilde*B_hat(:,:,k);
        
        
        A_hat_tilde(:,:,k) = T_tilde*A_hat(:,:,k);
        
        G_hat_tilde(:,:,k) = T_tilde*G_hat(:,:,k);
        
        T_tilde_prev = T_tilde;
        
    elseif k>=2
        
        T_tilde = pinv(O_hat_ref)*O_hat;
        
        if k<=q-1
            C_hat_tilde(:,:,k+1) = C_hat(:,:,k+1)*inv(T_tilde);
            
        end
        B_hat_tilde(:,:,k) = T_tilde*B_hat(:,:,k);
        G_hat_tilde(:,:,k) = T_tilde*G_hat(:,:,k);
        
        A_hat_tilde(:,:,k) = T_tilde*A_hat(:,:,k)*inv(T_tilde_prev);
        
        T_tilde_prev = T_tilde;
    end
    
end

end

