% testing control optimality of ARMA models with different q values. 

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

t_steps = 30;

no_rollouts = 50; 

U = zeros(nu*t_steps,no_rollouts);
y_matrix = zeros(nz*t_steps, no_rollouts);

ADD_PROC_NOISE = false;
ADD_MEAS_NOISE = false;

for i=1:no_rollouts
    %x0 = randn(n,1);
    
    x0 = zeros(n,1);
    
    %x0 = ones(n,1);
    
    u_vec = normrnd(0, 20, nu, t_steps); %perturbation
    %u_vec = zeros(nu, t_steps); u_vec(:,1) = ones(nu,1); x0 = zeros(n,1); %impulse
    
    if strcmp(system,'oscillator')
        [x, y] = generate_response_oscillator(x0, u_vec, n, nz, sysd.Ts);%output
    end
    
    U(:,i) = reshape(u_vec,nu*t_steps,1); 
    y_matrix(:,i) = reshape(y, nz*t_steps,1); 

end
fprintf('Response calculated for %i rollouts and %i time steps\n\n', no_rollouts, t_steps);

%% Choose q value.

q = 4; % number of markov parameters to estimate

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

