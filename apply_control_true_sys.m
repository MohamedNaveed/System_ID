function [X, UU, K, z, cost] = apply_control_true_sys(sysd, q_bar, x, y, u_vec, control_t_steps)
% lqr control on the true system 

n = size(sysd.A,1); % order of the system
nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

Oq = zeros(nz*q_bar,n);

for i=1:q_bar
    Oq((i-1)*nz + 1: i*nz,:) = sysd.C*sysd.A^(q_bar-i);
end

Oq_dagger = inv(Oq'*Oq)*Oq';


t = q_bar+1; %time instant

samples_taken = t-1:-1:t-q_bar;

Zq = zeros(nz*q_bar,1);
Uq = zeros(nu*q_bar,1);

for i = 1:length(samples_taken)
    
    Zq((i-1)*nz +1: i*nz) = y(:,samples_taken(i));
    Uq((i-1)*nu +1: i*nu) = u_vec(:,samples_taken(i));
end

Gq = zeros(nz*q_bar, nu*q_bar);

for i=1:q_bar %rows
    for j=1:q_bar %cols
        
        if (j - i -1 < 0)
            Gq((i-1)*nz +1:i*nz,(j-1)*nu +1: j*nu) = zeros(nz,nu);
        else
            Gq((i-1)*nz +1:i*nz,(j-1)*nu +1: j*nu) = sysd.C*sysd.A^(j-i-1)*sysd.B;
        end
        
    end
end
  
estimated_init_cond = Oq_dagger*(Zq - Gq*Uq)

disp('true initial condition')
x(:,t-q_bar)

%propagate initial condition till q
X = estimated_init_cond;

for i = 1:q_bar-1
    X = sysd.A*X + sysd.B*u_vec(:,i);
end

X_q = X; % initial condition for optimal control problem from q

% design lqr
Q = 10*eye(nz);
Q_N = 10*eye(nz); %terminal cost.
R = eye(nu);

Q_X = sysd.C'*Q*sysd.C;
Q_X_N = sysd.C'*Q_N*sysd.C;

K = finite_horizon_lqr(sysd.A, sysd.B, Q_X, R, Q_X_N, control_t_steps); 

%applying control

X = zeros(n, control_t_steps+1);
X(:,1) = X_q;

z = zeros(nz,control_t_steps+1);
z(:,1) = sysd.C*X_q;

UU = zeros(nu, control_t_steps);
cost = 0;

for i = 1:control_t_steps
    
    UU(:,i) = K(:,:,i)*X(:,i);
    cost = cost + 0.5*(z(:,i)'*Q*z(:,i) + UU(:,i)'*R*UU(:,i));
    
    X(:,i+1) = sysd.A*X(:,i) + sysd.B*UU(:,i); 
    z(:,i+1) = sysd.C*X(:,i+1);
end
cost = cost + 0.5*(z(:,control_t_steps+1)'*Q_N*z(:,control_t_steps+1));

disp('Done');
end

