function [Y_true] = calculate_true_markov_parameters_ltv(system,num_mp, t_steps)

if strcmp(system,'oscillator')
    sysd = oscillator(0);
end
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs
n = size(sysd.A,1);

Y_true = zeros(t_steps*nz, num_mp*nu);
 
del_t = sysd.Ts;


for t = 0:t_steps-1
    t
    Phi = eye(n);
    
    temp = oscillator(t*del_t);
    C_k = temp.C;
    
    if t+1 < num_mp
        
        i_uplimit = t;
        i_lowlimit = 0;
    else
        i_uplimit = t;
        i_lowlimit = t - num_mp+1;
    end
    
    for i = i_uplimit:-1:i_lowlimit
        i
        if i == i_uplimit

            Y_true(t*nz + 1:(t+1)*nz, 1:nu) = temp.D;

        elseif i == i_uplimit - 1

            temp = oscillator(i*del_t);

            B_k_1 = temp.B; %B(k-1)
            cols = nu + 1 : nu + nu;
            Y_true(t*nz + 1:(t+1)*nz, cols) = C_k*B_k_1;

        elseif i< i_uplimit - 1

            temp = oscillator((i+1)*del_t);
            A = temp.A;
            Phi = Phi*A; %output is recorded from t=1

            temp = oscillator((i)*del_t);
            B = temp.B;
            cols = nu + (i_uplimit-i-1)*nu + 1 : nu + (i_uplimit-i)*nu;
            Y_true(t*nz + 1:(t+1)*nz,cols) = C_k*Phi*B;
        end
    end
end
end

