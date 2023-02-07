clc; clear;

%load('data/fish_outputs11.mat');
load('data/cartpole_data_sept21.mat');
Model.name
nu = Model.nu; % number of control inputs
nz = Task.nm; % number of outputs
n = Model.nsys; %states in the system
q = 4; % steps to loop back

t_steps = Task.horizon+1; % horizon

no_rollouts = 300; 

U = flipud(delta_u(:,1:no_rollouts));
y_matrix = flipud(delta_z(1:nz*t_steps,1:no_rollouts));

%% build data (V) matrix and calculate the ARMA parameters

alpha_beta = zeros(nz*t_steps, q*(nz + nu) +  nu);

rank_V = [];
ID_time_idxs = 1:t_steps;

for k = ID_time_idxs
    
    V = build_data_mat_ltv(U, y_matrix, q, nu, nz, k, no_rollouts);
    
    if k<= q
        alpha_beta((k-1)*nz + 1: k*nz,1:(k-1)*(nu+nz)+nu) = y_matrix((k - 1)*nz + 1: (k)*nz, :)*pinv(V);%moore penros e inverse.
    else
        alpha_beta((k-1)*nz + 1: k*nz,:) = y_matrix((k - 1)*nz + 1: (k)*nz, :)*pinv(V);  
    end
    rank_V = [rank_V rank(V)];
    
end

fprintf('Estimated ARMA parameters\n\n');

%% Find open loop markov parameters
num_mp = t_steps;
markov_open_loop = calculate_open_loop_markov_para_ltv(nu, nz,...
                             num_mp, alpha_beta, ID_time_idxs, t_steps, q);

fprintf('Calculated open-loop markov parameters\n\n');

%% calculate observer gain matrix

h_o = calculate_observer_gain_markov(nu, nz, num_mp, alpha_beta,...
        t_steps, q);

fprintf('Calculated observer markov parameters \n\n');

%% build hankel to estimate A,B,C

[A_hat, B_hat, C_hat, D_hat, G_hat] = TVERA(Model.name, markov_open_loop, q, nu, nz, t_steps, h_o);

fprintf('Calculated A, B, C, G \n\n');

%% calculate open-loop markov parameters from A,B,C

markov_parameters_ABC = calculate_markov_from_ABC(A_hat, B_hat, C_hat, D_hat, q,t_steps,num_mp);

fprintf('Calculated open-loop markov parameters from A, B, C\n\n');

%% calculate observer in the loop markov parameters

closed_loop_markov_parameters = calculate_markov_from_ABCG(A_hat, B_hat,...
    C_hat, D_hat, G_hat, q, t_steps,nz,nu, alpha_beta);

fprintf('Calculated closed-loop markov parameters from A, B, C, G\n\n');

%% checking response for ARMA model
sysd.Ts = 0.1; %dummy
ZERO_INIT = 0;
n = size(A_hat,1);
y_matrix_test = flipud(delta_z(1:nz*t_steps,1+no_rollouts:end));
U_test = flipud(delta_u(:,1+no_rollouts:end));

[err_y_arma, err_y_OKID, err_y_ABC] = check_response(Model.name, alpha_beta, markov_open_loop,markov_parameters_ABC,...
                t_steps, q, nu, nz, n, sysd.Ts, num_mp, U_test, y_matrix_test, ...
                A_hat, B_hat, C_hat, D_hat, G_hat, ZERO_INIT);

fprintf(' Experiment done \n\n');
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

norm_err_arma = vecnorm(mean_err_y_arma,1);
norm_err_OKID = vecnorm(mean_err_y_OKID,1);
norm_err_ABC = vecnorm(mean_err_y_ABC,1);

%% plot
sample_id = 1;
plot_response(err_y_arma(:,:,1), err_y_OKID(:,:,1), t_steps, q);

%%
SAVE_PLOT = true;

plot_error_norm(norm_err_arma, norm_err_OKID, norm_err_ABC,q, SAVE_PLOT);

%%
%plot_error_stats(mean_err_y_arma, std_err_y_arma, mean_err_y_OKID,...
%    std_err_y_OKID, mean_err_y_ABC, std_err_y_ABC,q, SAVE_PLOT);
            