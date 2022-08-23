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

x0 = ones(n,1);

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

num_mp = 40; % number of markov parameters

Y_true = calculate_true_markov_parameters_ltv(system,num_mp);

%% build data (V) matrix and calculate the ARMA parameters

alpha_beta = zeros(nz*(t_steps -q), q*(nz + nu) +nu);

rank_V = [];
ID_time_idxs = q+1:t_steps;

for k = ID_time_idxs % k is the time step at which the ARMA is identified.

    V = build_data_mat_ltv(U, y_matrix, q, nu, nz, k, no_rollouts);

    alpha_beta((k-1)*nz + 1: k*nz,:) = y_matrix((k - 1)*nz + 1: (k)*nz, :)*pinv(V);
    rank_V = [rank_V rank(V)];
end

%% checking response. 

rollout_id = 1; % which rollout are you checking the data for

info_state = zeros((nz + nu)*q + nu, 1);

y_predicted = zeros(nz, t_steps);


for k = q+1:t_steps 
    
    info_state(1:nu,:) = U((k-1)*nu +1: k*nu,rollout_id)
    
    for i = 1:q

        info_state(nu + (i-1)*nz + 1: nu + i*nz) = y_matrix((k - 1 - i)*nz + 1: (k-i)*nz, rollout_id)

        info_state(nu + q*nz + (i-1)*nu + 1:nu + q*nz + i*nu) = U((k - 1 - i)*nu +1: (k-i)*nu, rollout_id)
    end

    y_predicted(:,k) = alpha_beta((k-1)*nz + 1: k*nz,:)*info_state

end

y_true = reshape(y_matrix(:,rollout_id), nz, t_steps);
err_y = y_true - y_predicted;

figure;
plot(1:t_steps, err_y, 'Linewidth',2);

xlabel('time steps');
ylabel('Error in prediction');