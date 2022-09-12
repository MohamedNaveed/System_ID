function [] = check_response(system, alpha_beta, markov_open_loop,markov_parameters_ABC, t_steps, q, nu, nz, n, Ts, num_mp)


u_vec = normrnd(0, 20, nu, t_steps);
y_predicted_arma = zeros(nz, t_steps);
y_predicted_open_loop = zeros(nz, t_steps);
y_predicted_ABC  = zeros(nz, t_steps);
x0 = zeros(n,1);

if strcmp(system,'oscillator')
        [x, y_true] = generate_response_oscillator(x0, u_vec, n, nz, Ts);%output
end

for k = 1:t_steps 
    
    if k > q
        
        info_state = zeros((nz + nu)*q + nu, 1);
        info_state(1:nu,:) = u_vec(:,k);
        
        control_state = zeros(nu*num_mp,1); % for moving average model
        control_state(1:nu,:) = u_vec(:,k);
        
        for i = 1:q

            info_state(nu + (i-1)*nz + 1: nu + i*nz) = y_predicted_arma(:,k-i);

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
        y_predicted_ABC(:,k) = markov_parameters_ABC((k-1)*nz + 1: (k)*nz, :)*control_state;
        
    else
        
        info_state = zeros((nz + nu)*(k-1) + nu, 1); %for arma model
        info_state(1:nu,:) = u_vec(:,k);
        
        control_state = zeros(nu*num_mp,1); % for moving average model
        control_state(1:nu,:) = u_vec(:,k);
        
        if k~=1
            for i = 1:k-1
                
                info_state(nu + (i-1)*nz + 1: nu + i*nz) = y_predicted_arma(:,k-i);
                info_state(nu + (k-1)*nz + (i-1)*nu + 1:nu + (k-1)*nz + i*nu) = u_vec(:,k-i);
                
                %mp_rows = (num_mp - (k-i-1) - 1)*nu + 1: (num_mp - (k-i-1))*nu;
                mp_rows = nu + (i - 1)*nu + 1:nu + (i)*nu;
                control_state(mp_rows,:) = u_vec(:,k-i);
            end
        end
        
        y_predicted_arma(:,k) = alpha_beta((k-1)*nz + 1: k*nz,1:(k-1)*(nu+nz)+nu)*info_state;
        
        y_predicted_open_loop(:,k) = markov_open_loop((k-1)*nz + 1: (k)*nz, :)*control_state;
        y_predicted_ABC(:,k) = markov_parameters_ABC((k-1)*nz + 1: (k)*nz, :)*control_state;
    end
end

err_y_arma = y_true - y_predicted_arma;
err_y_open_loop = y_true - y_predicted_open_loop;
err_y_ABC = y_true - y_predicted_ABC;

figure;
subplot(2,1,1)
plot(1:t_steps, err_y_arma(1,:),'b', 'Linewidth',3, 'DisplayName', 'ARMA output 1');
hold on;
plot(1:t_steps, err_y_open_loop(1,:),'--r', 'Linewidth',3, 'DisplayName', 'Open-loop out 1');
plot(1:(t_steps-q-1), err_y_ABC(1,1:(t_steps-q-1)),':','Color',	'[0.4660 0.6740 0.1880]', 'Linewidth',3, 'DisplayName', 'ABC model out 1');
ylabel('Error in prediction');
legend();

subplot(2,1,2)
plot(1:t_steps, err_y_arma(2,:),'b', 'Linewidth',3, 'DisplayName', 'ARMA output 2');
hold on;
plot(1:t_steps, err_y_open_loop(2,:),'--r', 'Linewidth',3, 'DisplayName', 'Open-loop out 2');
plot(1:(t_steps-q-1), err_y_ABC(2,1:(t_steps-q-1)),':','Color',	'[0.4660 0.6740 0.1880]', 'Linewidth',3, 'DisplayName', 'ABC model out 2');

xlabel('time steps');
ylabel('Error in prediction');

legend();

end

