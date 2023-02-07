function [err_y_arma, err_y_OKID, err_y_open_loop] = check_response(system, alpha_beta, markov_open_loop,...
    markov_parameters_ABC, t_steps, q, nu, nz, n, Ts, num_mp, U, y_matrix,...
    A_hat, B_hat, C_hat, D_hat, G_hat, ZERO_INIT)


N_mc = 100; %monte carlo

for n_mc = 1:N_mc
    
    
    y_predicted_arma = zeros(nz, t_steps);
    y_predicted_open_loop = zeros(nz, t_steps);
    y_predicted_ABC  = zeros(nz, t_steps);
    y_predicted_OKID = zeros(nz, t_steps);
    
    if strcmp(system,'oscillator')
            
            if ZERO_INIT
                x0 = zeros(n,1);
            else
                x0 = randn(n,1);
            end
            u_vec = normrnd(0, 20, nu, t_steps);
            [x, y_true] = generate_response_oscillator(x0, u_vec, n, nz, Ts);%output

    elseif strcmp(system,'cartpole')
        
        y_true = reshape(y_matrix(:,n_mc),[nz,t_steps]);
        u_vec = reshape(U(:,n_mc),[nu,t_steps]);
    
    elseif strcmp(system,'fish')
        
        y_true = reshape(y_matrix(:,n_mc),[nz,t_steps]);
        u_vec = reshape(U(:,n_mc),[nu,t_steps]);
    end
    

    %prediction using markov.
    for k = 1:t_steps 

        if k > q

            info_state = zeros((nz + nu)*q + nu, 1);
            info_state(1:nu,:) = u_vec(:,k);

            control_state = zeros(nu*num_mp,1); % for moving average model
            control_state(1:nu,:) = u_vec(:,k);

            for i = 1:q

                info_state(nu + (i-1)*nz + 1: nu + i*nz) = y_true(:,k-i);

                info_state(nu + q*nz + (i-1)*nu + 1:nu + q*nz + i*nu) = u_vec(:,k-i);
            end


            y_predicted_arma(:,k) = alpha_beta((k-1)*nz + 1: k*nz,:)*info_state;

            if k>num_mp
                i_idx = k-num_mp:k-1;
            else
                i_idx = 1:k-1;
            end
            for i = i_idx

                    %mp_rows = (num_mp - (k-i-1) - 1)*nu + 1: (num_mp - (k-i-1))*nu;
                    mp_rows = nu + (i - 1)*nu + 1:nu + (i)*nu;
                    control_state(mp_rows,:) = u_vec(:,k-i);
            end

            y_predicted_open_loop(:,k) = markov_open_loop((k-1)*nz + 1: (k)*nz, :)*control_state;
            %y_predicted_ABC(:,k) = markov_parameters_ABC((k-1)*nz + 1: (k)*nz, :)*control_state;

        else

            %info_state = zeros((nz + nu)*(k-1) + nu, 1); %for arma model
            %info_state(1:nu,:) = u_vec(:,k);

            control_state = zeros(nu*num_mp,1); % for moving average model
            control_state(1:nu,:) = u_vec(:,k);

            if k~=1
                for i = 1:k-1

                    %info_state(nu + (i-1)*nz + 1: nu + i*nz) = y_predicted_arma(:,k-i);
                    %info_state(nu + (k-1)*nz + (i-1)*nu + 1:nu + (k-1)*nz + i*nu) = u_vec(:,k-i);

                    %mp_rows = (num_mp - (k-i-1) - 1)*nu + 1: (num_mp - (k-i-1))*nu;
                    mp_rows = nu + (i - 1)*nu + 1:nu + (i)*nu;
                    control_state(mp_rows,:) = u_vec(:,k-i);
                end
            end

            %y_predicted_arma(:,k) = alpha_beta((k-1)*nz + 1: k*nz,1:(k-1)*(nu+nz)+nu)*info_state;

            y_predicted_open_loop(:,k) = markov_open_loop((k-1)*nz + 1: (k)*nz, :)*control_state;
            %y_predicted_ABC(:,k) = markov_parameters_ABC((k-1)*nz + 1: (k)*nz, :)*control_state;
            
            y_predicted_arma(:,k) = y_true(:,k);
            
        end
    end

    % prediction using A,B,C,D,G
    
    y_stack = zeros(q*nz,1);
    Gq_matrix = zeros(q*nz,q*nu);
    u_stack = zeros(q*nu,1);
    Oq = zeros(q*nz,n);
    
    % building data to reconstruct initial condition
    for k = 1:q
       
        y_stack((q-k)*nz + 1:(q-k + 1)*nz,:) = y_true(:,k);
            
        Gq_matrix((q-k)*nz +1:(q-k+1)*nz,(q-k)*nu + 1:end) = ...
                markov_open_loop((k-1)*nz + 1:(k)*nz,1:k*nu);
        
        
        u_stack((q-k)*nu + 1:(q-k+1)*nu,:) = u_vec(:,k);
  
        Phi = eye(n,n);
        if k>1
            for j = k-1:-1:1

               Phi = Phi*A_hat(:,:,j);
 
            end
        end
        Oq((q-k)*nz + 1:(q-k+1)*nz,:) = C_hat(:,:,k)*Phi;
    end
    x0_hat = pinv(Oq)*(y_stack - Gq_matrix*u_stack);%reconstructing initial condition
    
    x_hat = x0_hat; 
    
    for k = 1:t_steps

        y_predicted_OKID(:,k) = C_hat(:,:,k)*x_hat + D_hat*u_vec(:,k);

        x_hat = (A_hat(:,:,k) + G_hat(:,:,k)*C_hat(:,:,k))*x_hat + ...
            (B_hat(:,:,k) + G_hat(:,:,k)*D_hat)*u_vec(:,k) - ...
            G_hat(:,:,k)*y_true(:,k);
    end
    
    x_hat = x0_hat;
    for k = 1:t_steps

        y_predicted_ABC(:,k) = C_hat(:,:,k)*x_hat + D_hat*u_vec(:,k);

        x_hat = A_hat(:,:,k)*x_hat + B_hat(:,:,k)*u_vec(:,k);

    end
    err_y_arma(:,:,n_mc) = (y_true - y_predicted_arma)./abs(y_true);
    err_y_open_loop(:,:,n_mc) = (y_true - y_predicted_open_loop)./abs(y_true);
    err_y_ABC(:,:,n_mc) = (y_true - y_predicted_ABC)./abs(y_true);
    err_y_OKID(:,:,n_mc) = (y_true - y_predicted_OKID)./abs(y_true);

end
end

