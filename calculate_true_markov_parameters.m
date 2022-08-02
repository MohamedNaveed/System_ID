function [Y_true] = calculate_true_markov_parameters(sysd,num_mp)

nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

Y_true = zeros(nz, num_mp*nu);

Y_true(:,1:nu) = sysd.D;

for i = 1:num_mp-1
   
    Y_true(:,nu + 1 + (i-1)*nu:nu + i*nu) = sysd.C*(sysd.A^(i-1))*sysd.B;
    
end
end
