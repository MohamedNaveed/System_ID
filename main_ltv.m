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

t_steps = 100;
u_vec = normrnd(0, 20, nu, t_steps); %perturbation
%u_vec = zeros(nu, t_steps); u_vec(:,1) = ones(nu,1); x0 = zeros(n,1); %impulse


ADD_PROC_NOISE = false;
ADD_MEAS_NOISE = false;

if strcmp(system,'oscillator')
    [x, y] = generate_response_oscillator(x0, u_vec, n, nz, sysd.Ts);%output
end
%% plot response

figure;
plot(1:t_steps, y, 'Linewidth',2);
xlabel('time steps');
ylabel('output');

%% Choose q value.

q = 2; % number of markov parameters to estimate


%% true open loop markov parameters

num_mp = 40; % number of markov parameters

Y_true = calculate_true_markov_parameters_ltv(system,num_mp);




