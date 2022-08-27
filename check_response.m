function [] = check_response(system, alpha_beta, t_steps, q, nu, nz, n, Ts)


u_vec = normrnd(0, 20, nu, t_steps);
y_predicted = zeros(nz, t_steps);
x0 = zeros(n,1);

if strcmp(system,'oscillator')
        [x, y_true] = generate_response_oscillator(x0, u_vec, n, nz, Ts);%output
end

for k = 1:t_steps 
    
    if k > q
        
        info_state = zeros((nz + nu)*q + nu, 1);
        info_state(1:nu,:) = u_vec(:,k);

        for i = 1:q

            info_state(nu + (i-1)*nz + 1: nu + i*nz) = y_true(:,k-i);

            info_state(nu + q*nz + (i-1)*nu + 1:nu + q*nz + i*nu) = u_vec(:,k-i);
        end

        y_predicted(:,k) = alpha_beta((k-1)*nz + 1: k*nz,:)*info_state;
    else
        
        info_state = zeros((nz + nu)*(k-1) + nu, 1);
        info_state(1:nu,:) = u_vec(:,k);
        
        if k~=1
            for i = 1:k-1
                
                info_state(nu + (i-1)*nz + 1: nu + i*nz) = y_true(:,k-i);
                info_state(nu + (k-1)*nz + (i-1)*nu + 1:nu + (k-1)*nz + i*nu) = u_vec(:,k-i);
            end
        end
        
        y_predicted(:,k) = alpha_beta((k-1)*nz + 1: k*nz,1:(k-1)*(nu+nz)+nu)*info_state;
    end
end

err_y = y_true - y_predicted;

figure;
plot(1:t_steps, err_y, 'Linewidth',2);

xlabel('time steps');
ylabel('Error in prediction');

end

