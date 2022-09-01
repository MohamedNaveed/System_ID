% System ID main. 
clc;clear;

% get system 
system = 'oscillator'

if strcmp(system,'oscillator')
    sysd = oscillator(0);
end
rng(0);
n = size(sysd.A,1); % order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

x0 = zeros(n,1);
%x0 = ones(n,1);

t_steps = 20;

no_rollouts = 50; 

U = zeros(nu*t_steps,no_rollouts);
y_matrix = zeros(nz*t_steps, no_rollouts);

ADD_PROC_NOISE = false;
ADD_MEAS_NOISE = false;

for i=1:no_rollouts
    
    u_vec = normrnd(0, 20, nu, t_steps); %perturbation
    %u_vec = zeros(nu, t_steps); u_vec(:,1) = ones(nu,1); x0 = zeros(n,1); %impulse
    
    if strcmp(system,'oscillator')
        [x, y] = generate_response_oscillator(x0, u_vec, n, nz, sysd.Ts);%output
    end
    
    U(:,i) = reshape(u_vec,nu*t_steps,1); 
    y_matrix(:,i) = reshape(y, nz*t_steps,1); 

end
%% plot response
%{
figure;
plot(1:t_steps, y, 'Linewidth',2);
xlabel('time steps');
ylabel('output');
%}
%% Choose q value.

q = 3; % number of markov parameters to estimate


%% true open loop markov parameters

num_mp = 10; % number of markov parameters

Y_true = calculate_true_markov_parameters_ltv(system,num_mp);

%% build data (V) matrix and calculate the ARMA parameters
alpha_beta = zeros(nz*t_steps, q*(nz + nu) +  nu);

rank_V = [];
ID_time_idxs = 1:t_steps;

for k = ID_time_idxs
    
    V = build_data_mat_ltv(U, y_matrix, q, nu, nz, k, no_rollouts);
    
    if k<= q
        alpha_beta((k-1)*nz + 1: k*nz,1:(k-1)*(nu+nz)+nu) = y_matrix((k - 1)*nz + 1: (k)*nz, :)*pinv(V);%moore penrose inverse.
    else
        alpha_beta((k-1)*nz + 1: k*nz,:) = y_matrix((k - 1)*nz + 1: (k)*nz, :)*pinv(V);  
    end
    rank_V = [rank_V rank(V)];
    
end

%% checking response for ARMA model

check_response(system, alpha_beta, t_steps, q, nu, nz, n, sysd.Ts);

%% free response experiment to calculate the markov parameters for first q steps. 

[A_hat, B_hat, C_hat] = free_response_exp(system, q, alpha_beta);

%% Find open loop markov parameters

markov_open_loop = calculate_open_loop_markov_para_ltv(nu, nz,...
                             num_mp, alpha_beta, ID_time_idxs, t_steps, q);