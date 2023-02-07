% System ID main. LTI case
clc;clear;

% get system 
%sysd = spring_mass_damper();
sysd = OKID_paper_sys();
%sysd = one_dim_sys();

rng(0);
n = size(sysd.A,1); % order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

x0 = zeros(n,1);

t_steps = 40;
u_vec = normrnd(0, 20, 1, t_steps); %perturbation

ADD_PROC_NOISE = false;
ADD_MEAS_NOISE = false;

[x, y] = generate_response(x0, u_vec, sysd, ADD_PROC_NOISE, ADD_MEAS_NOISE);%output

disp('ARMA parameter:');
q = 10 % number of markov parameters to estimate


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

Y_bar = y_bar*pinv(V_bar); % markov parameters 

%% reconstruct initial condition experiment

Y_bar_from_xtq = reconstruct_initial_condition_exp(sysd, q, x, y, u_vec, Y_bar);

error_Y_bar = Y_bar - Y_bar_from_xtq;
%% control using the ARMA model.
q_bar = ceil(n/nz); %minimal q

control_steps = 50;
[Z, U_Z, K_Z, z_Z, cost_Z] = apply_control_ARMA(Y_bar, sysd, q,q_bar, y, u_vec, control_steps);

[Z_1, U_Z_1, K_Z_1, z_Z_1, cost_Z_1] = apply_control_ARMA(Y_bar_from_xtq, sysd, q,q_bar, y, u_vec, control_steps);

disp('difference in cost armas(normalised):');
(cost_Z - cost_Z_1)/cost_Z_1

%% control on the true model

[X, U_X, K_X, z_X, cost_X] = apply_control_true_sys(sysd, q_bar, x, y, u_vec, control_steps);

disp('difference in cost (normalised):');
(cost_Z - cost_X)/cost_X

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
%}

%% plotting error
font_size = 14;
error_z = (z_Z - z_X)./abs(z_X);
error_u = (U_Z - U_X)./abs(U_X);
figure(2);

subplot(2,1,1);
plot(1:size(z_Z,2),error_z(1,:),'LineWidth',3);
title('Error in output(normalized)','FontSize',font_size);
ylabel('Error output 1','FontSize',font_size);

subplot(2,1,2);
plot(1:size(z_Z,2),error_z(2,:),'LineWidth',3);
ylabel('Error output 2','FontSize',font_size);
xlabel('time steps','FontSize',font_size);

figure(3);

plot(1:size(error_u,2),error_u,'LineWidth',3);
title('Error in control (normalized)','FontSize',font_size);
ylabel('Error','FontSize',font_size);
xlabel('time steps','FontSize',font_size);


%%
error_z_ARMAs =  (z_Z - z_Z_1)./abs(z_Z_1);
error_u_ARMAs = (U_Z - U_Z_1)./abs(U_Z_1);

figure(4);

subplot(2,1,1);
plot(1:size(error_z_ARMAs,2),error_z_ARMAs(1,:),'LineWidth',3);
title('Error in output ARMA(normalized)','FontSize',font_size);
ylabel('Error output 1','FontSize',font_size);

subplot(2,1,2);
plot(1:size(error_z_ARMAs,2),error_z_ARMAs(2,:),'LineWidth',3);
ylabel('Error output 2','FontSize',font_size);
xlabel('time steps','FontSize',font_size);

figure(5);

plot(1:size(error_u_ARMAs,2),error_u_ARMAs,'LineWidth',3);
title('Error in control ARMA(normalized)','FontSize',font_size);
ylabel('Error','FontSize',font_size);
xlabel('time steps','FontSize',font_size);