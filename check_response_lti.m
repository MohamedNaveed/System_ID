function [err_y_arma,err_y_OKID,err_y_ABC] = check_response_lti(sysd, q, Y_bar,...
    markov_open_loop, A_est, B_est, C_est, D_est, M)

N_mc = 100;

n = size(sysd.A,1); % order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs
t_steps = 30;

for n_mc = 1:N_mc
    
    x0 = rand(n,1);

    
    u_vec = normrnd(0, 20, 1, t_steps); %perturbation

    ADD_PROC_NOISE = false;
    ADD_MEAS_NOISE = false;

    [x, y_true] = generate_response(x0, u_vec, sysd, ADD_PROC_NOISE, ADD_MEAS_NOISE);

    y_predicted_arma = zeros(nz, t_steps);
    y_predicted_ABC  = zeros(nz, t_steps);
    y_predicted_OKID = zeros(nz, t_steps);
    
    % ARMA
    for k =1:t_steps
        
        
        if k>q
            info_state = zeros((nz + nu)*q + nu, 1);
            info_state(1:nu,:) = u_vec(:,k);
            
            for i = 1:q
                info_state(nu + (i-1)*(nu + nz) + 1: nu + i*(nu + nz)) = [u_vec(:,k-i);y_true(:,k-i)];
            end
            
            y_predicted_arma(:,k) = Y_bar(:,1:nu+q*(nu+nz))*info_state;
        else
            y_predicted_arma(:,k) = y_true(:,k);
        end
        
    end
    
    % ABC, ABCMD
    y_stack = zeros(q*nz,1);
    Gq_matrix = zeros(q*nz,q*nu);
    u_stack = zeros(q*nu,1);
    Oq = zeros(q*nz,n);
    
    for k =1:q
        y_stack((q-k)*nz + 1:(q-k + 1)*nz,:) = y_true(:,k);
        
        Gq_matrix((q-k)*nz +1:(q-k+1)*nz,(q-k)*nu + 1:end) = ...
                markov_open_loop(:,1:k*nu);
            
        u_stack((q-k)*nu + 1:(q-k+1)*nu,:) = u_vec(:,k);
    end
    
    

    for i = 1:q
    
        Oq((q - i)*nz+1:(q-i+1)*nz,:) = C_est*A_est^(i-1);
    end
    x0_hat = pinv(Oq)*(y_stack - Gq_matrix*u_stack);
    
    x_hat = x0_hat; 
    
    for k = 1:t_steps

        y_predicted_OKID(:,k) = C_est*x_hat + D_est*u_vec(:,k);

        x_hat = (A_est + M*C_est)*x_hat + (B_est + M*D_est)*u_vec(:,k) -...
                    M*y_true(:,k);
    end
    
    x_hat = x0_hat;
    for k = 1:t_steps

        y_predicted_ABC(:,k) = C_est*x_hat + D_est*u_vec(:,k);

        x_hat = A_est*x_hat + B_est*u_vec(:,k);

    end
    
    
    err_y_arma(:,:,n_mc) = (y_true - y_predicted_arma)./abs(y_true);
    err_y_ABC(:,:,n_mc) = (y_true - y_predicted_ABC)./abs(y_true);
    err_y_OKID(:,:,n_mc) = (y_true - y_predicted_OKID)./abs(y_true);
end
end

