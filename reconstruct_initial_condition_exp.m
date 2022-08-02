function [Y_bar_from_xtq] = reconstruct_initial_condition_exp(sysd, q, x, y, u_vec, Y_bar)

n = size(sysd.A,1); % order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

Oq = zeros(nz*q,n);

for i=1:q
    Oq((i-1)*nz + 1: i*nz,:) = sysd.C*sysd.A^(q-i);
end

Oq_dagger = inv(Oq'*Oq)*Oq';

t = 100; %time instant

samples_taken = t-1:-1:t-q;

Zq = zeros(nz*q,1);
Uq = zeros(nu*q,1);

for i = 1:length(samples_taken)
    
    Zq((i-1)*nz +1: i*nz) = y(:,samples_taken(i));
    Uq((i-1)*nu +1: i*nu) = u_vec(:,samples_taken(i));
end

Gq = zeros(nz*q, nu*q);

for i=1:q %rows
    for j=1:q %cols
        
        if (j - i -1 < 0)
            Gq((i-1)*nz +1:i*nz,(j-1)*nu +1: j*nu) = zeros(nz,nu);
        else
            Gq((i-1)*nz +1:i*nz,(j-1)*nu +1: j*nu) = sysd.C*sysd.A^(j-i-1)*sysd.B;
        end
        
    end
end
  
estimated_init_cond = Oq_dagger*(Zq - Gq*Uq)

disp('true initial condition')
x(:,t-q)

Y_true = calculate_true_markov_parameters(sysd,q + 1); %true markov parameters (+1 to account for D)

alpha_bar = sysd.C*(sysd.A^q)*Oq_dagger;
beta_bar = Y_true(:,nu + 1:nu + q*nu) - alpha_bar*Gq;

Y_bar_from_xtq = zeros(nz,nu + q*(nu+nz));
Y_bar_from_xtq(:,1:nu) = sysd.D;

for i=1:q
    Y_bar_from_xtq(:,nu + (i-1)*(nu+nz) + 1: nu + i*(nu+nz)) = ...
        [beta_bar(:,(i-1)*(nu) + 1: (i)*(nu)),alpha_bar(:,(i-1)*nz + 1: i*nz)];
    
end

figure;
plot(1:length(Y_bar_from_xtq),Y_bar(:,1:nu + q*(nu+nz)) - Y_bar_from_xtq,'Linewidth',2);
ylabel('error');

end

