function [outputArg1,outputArg2] = free_response_exp(system, q)

if strcmp(system,'oscillator')
    sysd = oscillator(0);
end
rng(0);
n = size(sysd.A,1); % order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

N = 10; %number of experiments

X0 = normrnd(0, 20, n, N);

y_matrix = zeros(nz*t_steps, N);
t_steps = 10;

u_vec = zeros(nu, t_steps);

for i=1:N
    
    x0 = X0(:,i); 
    if strcmp(system,'oscillator')
        [x, y] = generate_response_oscillator(x0, u_vec, n, nz, sysd.Ts);%output
    end
   
    y_matrix(:,i) = reshape(y, nz*t_steps,1); 

end

end

