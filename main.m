% System ID main. 
clc;clear;
% get system 

%sysd = spring_mass_damper();
sysd = OKID_paper_sys();

n = size(sysd.A,1); %order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs


x0 = zeros(n,1);

t_steps = 1000;
u_vec = normrnd(0, 1, 1, t_steps); %perturbation

ADD_NOISE = false;
[x, y] = generate_response(x0, u_vec, sysd, ADD_NOISE);%output


q = 40; % number of markov parameters to estimate

%% true open loop markov parameters

num_mp = 100; % number of markov parameters
Y_true = zeros(size(sysd.C,1), num_mp*nu);

Y_true(:,1:nu) = sysd.D;

for i = 1:num_mp-1
   
    Y_true(:,nu + 1 + (i-1)*nu:nu + i*nu) = sysd.C*(sysd.A^(i-1))*sysd.B;
    
end

%% build V matrix

V = build_V(u_vec, y, q);

%% calculate markov parameters

% leave out first p parameters for initial condition to die out. 
%sample_idxs = [50:100,200:250,350:400]; %choose samples time steps
sample_idxs = [q+1:t_steps];

y_bar = y(:,sample_idxs); %choosing samples
V_bar = V(:,sample_idxs);

disp('Number of rows in V_bar');
size(V_bar,1)

disp('rank of matrix V_bar');
rank_V_bar = rank(V_bar)

%idxs = [1,2,4,5,6,7];
%V_bar_rows = V_bar(idxs,:);
%alpha = V_bar(3,:)*V_bar_rows'*inv(V_bar_rows*V_bar_rows')

%Y_bar = y_bar*V_bar'*inv(V_bar*V_bar');
%Y_bar = y_bar/V_bar;

Y_bar = y_bar*pinv(V_bar);

%% Find D

D_est = Y_bar(:,1);

%% Find open loop markov parameters

Y = zeros(size(sysd.C,1), num_mp*nu);

Y(:,1:nu) = D_est;
Y_bar = [Y_bar zeros(size(sysd.C,1), (nu+nz)*(num_mp - q))]; %appending 0's to find the other markov parameters

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

%% error in markov parameters

err_Y = (Y_true - Y);

%% Hankel 
alpha = n + 20;
beta = n + 20;

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

Y_est = zeros(size(sysd.C,1), (num_mp -1)*nu);

for i=1:num_mp-1
   
    Y_est(:,(i-1)*nu+1:i*nu) = C_est*(A_est^(i-1))*B_est;
    
end
Y_est = [D_est Y_est];

err_Y_est = (Y_true - Y_est);

%% estimation of observer gain M

Y_o = zeros((num_mp)*nz, nz);

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

O_mat = zeros(nz*(num_mp), order);

for i = 1:num_mp
    
    O_mat((i-1)*nz+1:i*nz,:) = C_est*A_est^(i-1);
end

M = inv(O_mat'*O_mat)*O_mat'*Y_o

%% observer system 

%A_bar = sysd.A + M*sysd.C;
A_bar = A_est + M*C_est;
%B_bar = [sysd.B + M*sysd.D, -M];
B_bar = [B_est + M*D_est, -M];

disp('Eigenvalues of A_bar = ');
eig(A_bar)
disp('A_bar^q =');
A_bar^q


%%
Y_bar_est = D_est;

for i = 1:num_mp
    
    Y_bar_est = [Y_bar_est, C_est*(A_bar)^(i-1)*B_bar];
end

err_Y_bar_est = (Y_bar_est - Y_bar);
%%
%{
figure;
plot(1:t_steps, y, 'Linewidth',2);
xlabel('time steps');
ylabel('output');
%}

figure(1);
subplot(3,1,1);
hold on;
for j=1:size(err_Y,1)
    plot(1:length(err_Y),err_Y(j,:),'Linewidth',2); 
end
ylabel('error $Y$','Interpreter','latex','fontsize',16);

subplot(3,1,2);
hold on;
for j=1:size(err_Y_est,1)
    plot(1:length(err_Y_est),err_Y_est(j,:),'Linewidth',2);
end
ylabel('error $Y_{est}$','Interpreter','latex','fontsize',16);

subplot(3,1,3);
hold on;
for j=1:size(err_Y_bar_est,1)
    plot(1:length(err_Y_bar_est),err_Y_bar_est(j,:),'Linewidth',2);
end
xlabel('Markov parameters');
ylabel('error $\bar{Y}_{est}$','Interpreter','latex','fontsize',16);

%% reconstruct initial condition
%{
Oq = zeros(nz*q,n);

for i=1:q
    Oq((i-1)*nz + 1: i*nz,:) = sysd.C*sysd.A^(q-i);
end

Oq_dagger = inv(Oq'*Oq)*Oq';

samples_taken = q+5:-1:1+5;
Zq = y(samples_taken);
Uq = u_vec(samples_taken);

Gq = zeros(nz*q,m*q);

Gq = [0, sysd.C*sysd.B, sysd.C*sysd.A*sysd.B, sysd.C*sysd.A^2*sysd.B, sysd.C*sysd.A^3*sysd.B; ...
      0, 0, sysd.C*sysd.B, sysd.C*sysd.A*sysd.B, sysd.C*sysd.A^2*sysd.B;...
      0, 0, 0, sysd.C*sysd.B, sysd.C*sysd.A*sysd.B; ...
      0, 0, 0, 0, sysd.C*sysd.B;
      0, 0, 0, 0, 0];
  
init_cond = Oq_dagger*(Zq' - Gq*Uq')

alphas = sysd.C*(sysd.A^q)*Oq_dagger;
%}