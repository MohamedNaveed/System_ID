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

q = 3; % number of markov parameters to estimate

%% true open loop markov parameters

num_mp = t_steps; % number of markov parameters

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


%% Eigen realization algorithm ERA

[A_est, B_est, C_est] = ERA(sysd, Y);

D_est = Y_bar(:,1);

Y_est = zeros(size(sysd.C,1), (num_mp -1)*nu);

for i=1:num_mp-1
   
    Y_est(:,(i-1)*nu+1:i*nu) = C_est*(A_est^(i-1))*B_est;
    
end
Y_est = [D_est Y_est];


%% estimation of observer gain M

[M, LS_error_M] = calculate_M(sysd, num_mp, Y_bar, A_est, C_est);

%% observer system 

A_bar = A_est + M*C_est;
B_bar = [B_est + M*D_est, -M];

disp('Eigenvalues of A_bar = ');
eig(A_bar)
disp('C_est A_bar^q B_bar =');
C_est*(A_bar^q)*B_bar



%% Y_bar_est

Y_bar_est = D_est;

for i = 1:q
    
    Y_bar_est = [Y_bar_est, C_est*(A_bar)^(i-1)*B_bar];
end


%%
%{
figure;
plot(1:t_steps, y, 'Linewidth',2);
xlabel('time steps');
ylabel('output');
%}

%%data refining
thres = 1e-12;
Y_true(abs(Y_true) < thres) = 0;
Y_bar(abs(Y_bar) < thres ) = 0;
Y_bar_est(abs(Y_bar_est) < thres) = 0;
Y_est(abs(Y_est) < thres) = 0;
Y(abs(Y) < thres) = 0;

err_Y = (Y_true - Y)./abs(Y_true);%% error in markov parameters
err_Y(isnan(err_Y)) = 0;

err_Y_bar_est = (Y_bar_est - Y_bar(:,1:nu+q*(nu+nz)))./abs(Y_bar(:,1:nu+q*(nu+nz)));
err_Y_bar_est(isnan(err_Y_bar_est)) = 0;

err_Y_est = (Y_true - Y_est)./abs(Y_true);

norm_err_Y = vecnorm(err_Y,1);
norm_err_Y_bar = vecnorm(err_Y_bar_est,1);
%%plotting

SAVE_PLOT = true;


fig = figure;
subplot(2,1,1);
hold on;
plot(1:length(norm_err_Y_bar),norm_err_Y_bar,'Linewidth',2);
ylabel('$||error\ \bar{Y}||_1$','Interpreter','latex','fontsize',16);
ylim([-0.5,3])
subplot(2,1,2);
plot(1:length(norm_err_Y),norm_err_Y,'Linewidth',2); 
ylim([-1e-11,4e-11])
xlabel('Index');
ylabel('$||error\ Y||_1$','Interpreter','latex','fontsize',16);

%{
subplot(3,1,2);
hold on;
for j=1:size(err_Y_est,1)
    plot(1:length(err_Y_est),err_Y_est(j,:),'Linewidth',2);
end
ylabel('error $Y_{est}$','Interpreter','latex','fontsize',16);
%}

hold off;
if SAVE_PLOT
    set(fig,'Units','inches');
    fig.Position = [100,100,4.5,3.5];
    screenposition = get(fig,'Position');
    set(fig,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters '/home/naveed/Dropbox/Research/Manuscripts/ACC23/plots/okid_sys_markov_error_q=4.pdf'
end
%% 

figure;
for j=1:size(err_Y_bar_est,1)
    plot(1:length(err_Y_bar_est(:,1:nu+q*(nu+nz))),err_Y_bar_est(j,1:nu+q*(nu+nz)),'Linewidth',2);
end
xlabel('Markov parameters');
ylabel('error $\bar{Y}_{est}$','Interpreter','latex','fontsize',16);


%% check response

[err_y_arma,err_y_OKID,err_y_ABC] = check_response_lti(sysd, q, Y_bar, Y,...
                    A_est, B_est, C_est, D_est, M);

%% observer canonical form
[err_y_occ, err_y_arma] = response_observer_form(sysd, q, Y_bar);

%% error analysis

%remove outliers
err_y_arma = rmoutliers(err_y_arma,3);
err_y_OKID = rmoutliers(err_y_OKID,3);
err_y_ABC = rmoutliers(err_y_ABC,3);

mean_err_y_arma = mean(err_y_arma,3,'omitnan');
mean_err_y_OKID = mean(err_y_OKID,3,'omitnan');
mean_err_y_ABC = mean(err_y_ABC,3,'omitnan');

std_err_y_arma = std(err_y_arma,0,3,'omitnan');
std_err_y_OKID = std(err_y_OKID,0,3,'omitnan');
std_err_y_ABC = std(err_y_ABC,0,3,'omitnan');

SAVE_PLOT = false;
plot_error_stats(mean_err_y_arma, std_err_y_arma, mean_err_y_OKID, std_err_y_OKID,...
    mean_err_y_ABC, std_err_y_ABC, q, SAVE_PLOT);

%% reconstruct initial condition experiment

Y_bar_from_xtq = reconstruct_initial_condition_exp(sysd, q, x, y, u_vec, Y_bar);


%% control using the ARMA model.
% same q but using different arma parameters

%[Z, UU, K] = apply_control_ARMA(Y_bar, sysd, q, y, u_vec);
[Z_1, UU_1, K_1] = apply_control_ARMA(sysd, q, Y_bar_from_xtq, sysd, q, y, u_vec);

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