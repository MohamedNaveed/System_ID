function [] = check_response(system, alpha_beta, markov_open_loop,...
    markov_parameters_ABC, t_steps, q, nu, nz, n, Ts, num_mp, U, y_matrix)



y_predicted_arma = zeros(nz, t_steps);
y_predicted_open_loop = zeros(nz, t_steps);
y_predicted_ABC  = zeros(nz, t_steps);

rng(0);
if strcmp(system,'oscillator')
        x0 = zeros(n,1);
        u_vec = normrnd(0, 20, nu, t_steps);
        [x, y_true] = generate_response_oscillator(x0, u_vec, n, nz, Ts);%output
        
elseif strcmp(system,'cartpole')
    roll_out_id = 50;
    y_true = reshape(y_matrix(:,roll_out_id),[nz,t_steps]);
    u_vec = reshape(U(:,roll_out_id),[nu,t_steps]);
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

err_y_arma = (y_true - y_predicted_arma)./abs(y_true);
err_y_open_loop = (y_true - y_predicted_open_loop)./abs(y_true);
err_y_ABC = (y_true - y_predicted_ABC)./abs(y_true);

figure;
subplot(2,1,1)
plot(1:t_steps, err_y_arma(1,:),'b', 'Linewidth',3, 'DisplayName', 'ARMA');
hold on;
plot(1:t_steps, err_y_open_loop(1,:),'--r', 'Linewidth',3, 'DisplayName', 'Open-loop markov');
plot(1:(t_steps-q-1), err_y_ABC(1,1:(t_steps-q-1)),':','Color',	'[0.4660 0.6740 0.1880]', 'Linewidth',3, 'DisplayName', 'TV-OKID');
ylabel('Error in output-1');
%ylim([-0.1, 0.4]);
legend('Location', 'NorthWest');
%legend();
subplot(2,1,2)
plot(1:t_steps, err_y_arma(2,:),'b', 'Linewidth',3);
hold on;
plot(1:t_steps, err_y_open_loop(2,:),'--r', 'Linewidth',3);
plot(1:(t_steps-q-1), err_y_ABC(2,1:(t_steps-q-1)),':','Color',	'[0.4660 0.6740 0.1880]', 'Linewidth',3);
%ylim([-5e-12, 5e-12]);
xlabel('time steps');
ylabel('Error in output-2');

figure;
subplot(2,1,1)
plot(1:t_steps, err_y_arma(1,:),'b', 'Linewidth',3, 'DisplayName', 'ARMA');
hold on;
plot(1:t_steps, err_y_open_loop(1,:),'--r', 'Linewidth',3, 'DisplayName', 'Open-loop markov');
%plot(1:(t_steps-q-1), err_y_ABC(1,1:(t_steps-q-1)),':','Color',	'[0.4660 0.6740 0.1880]', 'Linewidth',3, 'DisplayName', 'TV-OKID');
ylabel('Error in output-1');
%ylim([-0.1, 0.4]);
legend('Location', 'NorthWest');
%legend();
subplot(2,1,2)
plot(1:t_steps, err_y_arma(2,:),'b', 'Linewidth',3);
hold on;
plot(1:t_steps, err_y_open_loop(2,:),'--r', 'Linewidth',3);
%plot(1:(t_steps-q-1), err_y_ABC(2,1:(t_steps-q-1)),':','Color',	'[0.4660 0.6740 0.1880]', 'Linewidth',3);
%ylim([-5e-12, 5e-12]);
xlabel('time steps');
ylabel('Error in output-2');
end

