function [M, LS_error] = calculate_M(sysd, num_mp, Y_bar, A_est, C_est)

n = size(sysd.A,1);
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

Y_o = zeros(num_mp*nz, nz);

for k = 1:num_mp
    
    temp = -Y_bar(:, nu + (k - 1)*(nu + nz)+ nu + 1:nu + k*(nu + nz));
    
    if (k ~= 1)
        for i=1:k-1
            temp = temp + Y_bar(:, nu + (i - 1)*(nu + nz)+ nu + 1:...
                nu + (i)*(nu + nz))*Y_o((k-i -1)*nz + 1:(k-i)*nz,:);
        end
    end
    
    Y_o((k-1)*nz + 1:k*nz,:) = temp;

end

O_mat = zeros(nz*(num_mp), n);

for i = 1:num_mp
    
    O_mat((i-1)*nz+1:i*nz,:) = C_est*A_est^(i-1);
end

M = inv(O_mat'*O_mat)*O_mat'*Y_o

LS_error = O_mat*M - Y_o;


end
