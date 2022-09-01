function [A_hat_tilde, B_hat_tilde, C_hat_tilde] = free_response_exp(system, q, alpha_beta)

if strcmp(system,'oscillator')
    sysd = oscillator(0);
end
rng(0);
n = size(sysd.A,1); % order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

N = q+10; %number of experiments

X0 = normrnd(0, 20, n, N);

t_steps = q+10;
y_matrix = zeros(nz*t_steps, N);

n_cols = 1; %number of time steps in Hankel columns
u_vec = zeros(nu, t_steps);

%calculate A,B,C
for k= 1:q+1
    k
    for i=1:N

        x0 = X0(:,i); 
        if strcmp(system,'oscillator')
            [x, y] = generate_response_oscillator(x0, u_vec, n, nz, sysd.Ts);%output
        end

        y_matrix(:,i) = reshape(y, nz*t_steps,1); 

    end

    [U,Sig,V] = svd(y_matrix((k-1)*nz + 1: (k+q-1)*nz,:));

    rank_Sig = rank(Sig);

    root_Sig = Sig(1:rank_Sig,1:rank_Sig)^(1/2);

    O_hat = U(:,1:rank_Sig)*root_Sig;
        
    X_hat = root_Sig*V(:,1:rank_Sig)';
    
    Hankel = zeros(nz*q, nz + nu);
    
    if k~=1 
        A_hat(:,:,k-1) = X_hat*pinv(X_hat_prev); %A calculation
        
        Hankel = build_hankel(alpha_beta,q, k, n_cols, nz, nu);
        
        B_hat(:,:,k-1) = pinv(O_hat)*Hankel;
    end

    if k<=q
        C_hat(:,:,k) = O_hat(1:nz,:);
    end
    X_hat_prev = X_hat;
    
    % transforming to reference coordinates. 
    
    if k==1
        C_hat_tilde(:,:,k) = C_hat(:,:,k);
        O_hat_k_1 = O_hat;
        
    elseif k==2
        
        T_tilde = pinv(O_hat_k_1)*O_hat;
        C_hat_tilde(:,:,k) = C_hat(:,:,k)*inv(T_tilde);
        
        B_hat_tilde(:,:,k-1) = B_hat(:,:,k-1);
        A_hat_tilde(:,:,k-1) = A_hat(:,:,k-1);
        T_tilde_prev = T_tilde;
        
    elseif k>=3
        
        T_tilde = pinv(O_hat_k_1)*O_hat;
        
        if k<=q
            C_hat_tilde(:,:,k) = C_hat(:,:,k)*inv(T_tilde);
        end
        B_hat_tilde(:,:,k-1) = T_tilde*B_hat(:,:,k-1);
        A_hat_tilde(:,:,k-1) = T_tilde*A_hat(:,:,k-1)*inv(T_tilde_prev);
        T_tilde_prev = T_tilde;
    end
    
end
end

