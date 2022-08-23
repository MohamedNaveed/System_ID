function [Y_true] = calculate_true_markov_parameters_ltv(system,num_mp)

if strcmp(system,'oscillator')
    sysd = oscillator(0);
end
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs
n = size(sysd.A,1);

Y_true = zeros(nz, num_mp*nu);

Y_true(:,1:nu) = sysd.D;
del_t = sysd.Ts;
Phi = eye(n);


B = sysd.B;

for i = 1:num_mp-1
    
    if strcmp(system,'oscillator')
        
        temp = oscillator(i*del_t);
        C = temp.C;
        
        if i>1
            temp = oscillator((i-1)*del_t);
            A = temp.A;
            Phi = A*Phi;
        end
        
    end
    
    
    Y_true(:,nu + 1 + (i-1)*nu:nu + i*nu) = C*Phi*B;
    
end

end

