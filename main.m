% System ID main. 
clc;clear;

% get system 
%sysd = spring_mass_damper();
%sysd = OKID_paper_sys();
sysd = one_dim_sys();

rng(0);
n = size(sysd.A,1); % order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

x0 = zeros(n,1);

t_steps = 5;
u_vec = normrnd(0, 20, 1, t_steps); %perturbation

ADD_PROC_NOISE = false;
ADD_MEAS_NOISE = false;

[x, y] = generate_response(x0, u_vec, sysd, ADD_PROC_NOISE, ADD_MEAS_NOISE);%output

q = 2; % number of markov parameters to estimate

%% true open loop markov parameters

num_mp = 40; % number of markov parameters

Y_true = calculate_true_markov_parameters(sysd,num_mp);

%% build data (V) matrix

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

Y_bar = [Y_bar zeros(nz, (nu+nz)*(num_mp - q))]; %appending 0's to find the other markov parameters

%% Find open loop markov parameters

Y = calculate_open_loop_markov_para(sysd,num_mp, Y_bar);

%% error in markov parameters

err_Y = (Y_true - Y);

%% Eigen realization algorithm ERA

[A_est, B_est, C_est] = ERA(sysd, Y);

D_est = Y_bar(:,1);

Y_est = zeros(size(sysd.C,1), (num_mp -1)*nu);

for i=1:num_mp-1
   
    Y_est(:,(i-1)*nu+1:i*nu) = C_est*(A_est^(i-1))*B_est;
    
end
Y_est = [D_est Y_est];

err_Y_est = (Y_true - Y_est);

%% estimation of observer gain M

[M, LS_error_M] = calculate_M(sysd, num_mp, Y_bar, A_est, C_est);

%% observer system 

A_bar = A_est + M*C_est;
B_bar = [B_est + M*D_est, -M];

disp('Eigenvalues of A_bar = ');
eig(A_bar)
disp('A_bar^q =');
A_bar^q


%% Y_bar_est

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
hold off;

%% reconstruct initial condition experiment

Y_bar_from_xtq = reconstruct_initial_condition_exp(sysd, q, x, y, u_vec, Y_bar);


%% control using the ARMA model.
% same q but using different arma parameters

%[Z, UU, K] = apply_control_ARMA(Y_bar, sysd, q, y, u_vec);
[Z_1, UU_1, K_1] = apply_control_ARMA(Y_bar_from_xtq, sysd, q, y, u_vec);

%%
%{
figure(2);

subplot(3,1,1);
plot(1:length(K),K-K_1,'LineWidth',3);
ylabel('Gain K');
title('Difference between the ARMA systems');
subplot(3,1,2);
plot(1:size(Z,2),Z(1,:)-Z_1(1,:),'LineWidth',3);
ylabel('Response');

subplot(3,1,3);
plot(1:size(UU,2),UU(1,:)-UU_1(1,:),'LineWidth',3);
ylabel('Control');
xlabel('time steps');


figure(2);

subplot(2,1,1);
plot(1:size(Z_1,2),Z_1_q(1,:)-Z_1(1,:),'LineWidth',3);
title('Difference between the ARMA systems');
ylabel('Response');

subplot(2,1,2);
plot(1:size(UU_1,2),UU_1_q(1,:)-UU_1(1,:),'LineWidth',3);
ylabel('Control');
xlabel('time steps');
%}

figure(2)
plot(1:size(Z_1,2),Z_1(1,:),'LineWidth',3)
ylabel('Response');