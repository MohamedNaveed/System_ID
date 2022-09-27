clc; clear;

load('data/cartpole_data_sept21.mat');
Model.name
nu = Model.nu; % number of control inputs
nz = Task.nm; % number of outputs
n = Model.nsys; %states in the system
q = 4; % steps to loop back

t_steps = size(delta_u,1); % horizon

no_rollouts = size(delta_u,2); 

U = flipud(delta_u);
y_matrix = flipud(delta_z(1:nz*t_steps,:));

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
[err_y_arma, err_y_OKID] = check_response(Model.name, alpha_beta, markov_open_loop,markov_parameters_ABC,...
                t_steps, q, nu, nz, n, sysd.Ts, num_mp, U, y_matrix, ...
                A_hat, B_hat, C_hat, D_hat, G_hat, ZERO_INIT);

fprintf(' Experiment done \n\n');
%% error analysis

%remove outliers
err_y_arma = rmoutliers(err_y_arma,3);
err_y_OKID = rmoutliers(err_y_OKID,3);

mean_err_y_arma = mean(err_y_arma,3,'omitnan');
mean_err_y_OKID = mean(err_y_OKID,3,'omitnan');

std_err_y_arma = std(err_y_arma,0,3,'omitnan');
std_err_y_OKID = std(err_y_OKID,0,3,'omitnan');

%% plot
sample_id = 1;
plot_response(err_y_arma(:,:,1), err_y_OKID(:,:,1), t_steps, q);

%%
SAVE_PLOT = false;
plot_error_stats(mean_err_y_arma, std_err_y_arma, mean_err_y_OKID, std_err_y_OKID,q, SAVE_PLOT);
            