function [Y] = calculate_open_loop_markov_para(sysd, num_mp, Y_bar)

nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

Y = zeros(nz, num_mp*nu);

D_est = Y_bar(:,1);
Y(:,1:nu) = D_est;

for k = 1:num_mp-1
    
    temp = Y_bar(:, nu + 1 + (k - 1)*(nu + nz):nu + (k - 1)*(nu + nz) + nu) + ...
            Y_bar(:, nu + (k - 1)*(nu + nz)+ nu + 1:nu + (k)*(nu + nz))*D_est;
    
    if (k ~= 1)
        for i=1:k-1
            temp = temp + Y_bar(:, nu + (i - 1)*(nu + nz)+ nu + 1:...
                       nu + (i)*(nu + nz))*Y(:,nu + (k-i-1)*nu + 1:nu + (k-i)*nu);
        end
    end
    
    Y(:,nu + (k-1)*nu + 1:nu + k*nu) = temp;
    
end
end

