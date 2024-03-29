function [err_y_occ, err_y_arma] = response_observer_form(sysd, q, Y_bar)


n = size(sysd.A,1); % order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs
t_steps = 30;

x0 = rand(n,1);
u_vec = normrnd(0, 20, 1, t_steps); %perturbation

ADD_PROC_NOISE = false;
ADD_MEAS_NOISE = false;

[x, y_true] = generate_response(x0, u_vec, sysd, ADD_PROC_NOISE, ADD_MEAS_NOISE);

y_predicted_arma = zeros(nz, t_steps);
y_predicted_occ = zeros(nz, t_steps);

A_occ = zeros(q*nz,q*nz);
B_occ = zeros(q*nz,nu);
C_occ = zeros(nz,q*nz);

for i = 1:q
    A_occ((i-1)*nz + 1: i*nz,1:nz) = Y_bar(:,nu + (i-1)*(nu+nz)+ nu + 1:nu + (i)*(nu+nz));
    if i~=q
        A_occ((i-1)*nz + 1: i*nz,i*nz + 1: (i+1)*nz) = eye(nz,nz);
    end
    B_occ((i-1)*nz + 1: i*nz,:) = Y_bar(:,nu + (i-1)*(nu+nz)+ 1:nu + (i-1)*(nu+nz) + nu);
end
C_occ(:,1:nz) = eye(nz,nz);

X0 = zeros(q*nz,1);
% ARMA
for k =1:t_steps


    if k>q
        info_state = zeros((nz + nu)*q + nu, 1);
        info_state(1:nu,:) = u_vec(:,k);

        for i = 1:q
            info_state(nu + (i-1)*(nu + nz) + 1: nu + i*(nu + nz)) = [u_vec(:,k-i);y_true(:,k-i)];
        end

        y_predicted_arma(:,k) = Y_bar(:,1:nu+q*(nu+nz))*info_state;
        
        X_hat = A_occ*X_hat + B_occ*u_vec(:,k-1);
        y_predicted_occ(:,k) = C_occ*X_hat;
    else
        y_predicted_arma(:,k) = y_true(:,k);
        
        if k==q
            info_state = zeros((nz + nu)*q + nu, 1);
            %info_state(1:nu,:) = u_vec(:,k);

            for i = 0:q-1
                info_state(nu + (i)*(nu + nz) + 1: nu + (i+1)*(nu + nz)) = [u_vec(:,k-i);y_true(:,k-i)];
            end
            
            for i = 1:q
                
                if i == 1
                    X0((i-1)*nz +1:i*nz) = y_true(:,q);%Y_bar(:,nu + (i-1)*(nu+nz) + nu+ 1:nu + (i)*(nu+nz))*
                else
                    X0((i-1)*nz +1:i*nz) = Y_bar(:,nu + (i-1)*(nu+nz)+1:nu + (q)*(nu+nz))...
                        *info_state(nu + (nu+nz) +1:nu + (q - i +2)*(nu+nz));
                end
            end
            X_hat = X0;
        end
    end

end

err_y_occ = y_predicted_occ - y_true;
err_y_arma = y_predicted_arma - y_true;
end

