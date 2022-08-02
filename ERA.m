function [A_est, B_est, C_est] = ERA(sysd, Y)

n = size(sysd.A,1);
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

alpha = n + 5;
beta = n + 5;

H_0 = zeros(alpha*nz, beta*nu);

for i=1:alpha
    
    H_0((i-1)*nz + 1:i*nz,:)  = Y(:,nu + i*nu: nu + (beta + i-1)*nu);

end

[R,Sig,S] = svd(H_0);

order = n;
P_b = R(:,1:order)*(Sig(1:order,1:order)^(1/2));
Q_b = S(:,1:order)*(Sig(1:order,1:order)^(1/2));

%% Estimation of A,B,C

H_1 = zeros(alpha*nz, beta*nu);

for i=1:alpha
    
    H_1((i-1)*nz + 1:i*nz,:)  = Y(:,nu + (i+1)*nu: nu + (beta + i)*nu);

end
    
inv_Sig = inv(Sig(1:order,1:order)^(1/2));
A_est = inv_Sig*R(:,1:order)'*H_1*S(:,1:order)*inv_Sig;
C_est = P_b(1:nz,:);
B_est = Q_b(1:nu,:)';


end
