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

%{ 
%batch implementation of the previous loop. Exactly the same answer.
P_mat = zeros((num_mp-1)*nz, (num_mp-1)*nu);
Q_mat = zeros((num_mp-1)*nz, nu);

for k = 1:num_mp -1
   
    P_mat((k-1)*nz + 1: k*nz,(k-1)*nz + 1: k*nz) = eye(nz,nz);
    
    if k >1
        for j = 1:k-1

            P_mat((k-1)*nz + 1: k*nz,(j-1)*nz + 1:j*nz) = -Y_bar(:, nu + (k-j-1)*(nu + nz)+ nu + 1:...
                       nu + (k-j)*(nu + nz));
        end
    end
    
    Q_mat((k-1)*nz + 1: k*nz,:) = Y_bar(:, nu + 1 + (k - 1)*(nu + nz):nu + (k - 1)*(nu + nz) + nu) + ...
            Y_bar(:, nu + (k - 1)*(nu + nz)+ nu + 1:nu + (k)*(nu + nz))*D_est;
end

temp = P_mat\Q_mat;

Y_hat = [D_est, reshape(temp,[nz,num_mp-1])];
    
%}
end













