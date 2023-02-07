function [Z, UU, K, z, cost] = apply_control_ARMA(Y_bar, sysd, q, q_bar, y, u_vec, control_t_steps)

nu = size(sysd.B,2); % number of control inputs
nz = size(sysd.C,1); % number of outputs

[alpha, beta] = Y_bar2alpha_beta(Y_bar, nu, nz, q); 

% POD2C realization
 
%building AA eq. 22
AA = zeros(q*nz + (q-1)*nu, q*nz + (q-1)*nu);

AA(1:nz,:) = [alpha, beta(:, nu+1:end)]; 

AA(nz + 1: q*nz,1:(q-1)*nz) = eye((q-1)*nz);

AA (q*nz + nu + 1:end, q*nz + 1:q*nz + (q-2)*nu) = eye((q-2)*nu);

%building BB eq. 22

BB = zeros(q*nz + (q-1)*nu,nu);

BB(1:nz,:) = beta(:,1:nu);

BB(q*nz + 1:q*nz + nu) = eye(nu);

CC = [eye(nz), zeros(nz,(q-1)*nz + (q-1)*nu)];
DD = zeros(q*nz + (q-1)*nu , nu);

% design lqr
Q = 10*eye(nz);
Q_N = 10*eye(nz); %terminal cost.
R = eye(nu);

Q_Z = CC'*Q*CC;
Q_Z_N = CC'*Q_N*CC;

K = finite_horizon_lqr(AA,BB,Q_Z, R,Q_Z_N,control_t_steps-(q-q_bar)); 


% choose initial condition
t = q;
Z0 = zeros(q*nz + (q-1)*nu,1);

for i=1:q
    
    Z0((i-1)*nz + 1: i*nz) = y(:,t-i+1);
end

for i=1:q-1
    
    Z0(q*nz + (i-1)*nu + 1:q*nz + i*nu) = u_vec(:,t-i);
end

% apply control

cost = 0;

if q > q_bar
   for j= q_bar:q-1
       
        cost = cost + 0.5*(y(:,j)'*Q*y(:,j) + u_vec(:,j)'*R*u_vec(:,j));
   end
end

Z = zeros(q*nz + (q-1)*nu, control_t_steps-(q-q_bar)+1);
z = zeros(nz, control_t_steps-(q-q_bar)+1); 

Z(:,1) = Z0;
z(:,1) = CC*Z0;

UU = zeros(nu, control_t_steps-(q-q_bar));


for i = 1:control_t_steps-(q-q_bar)

    UU(:,i) = K(:,:,i)*Z(:,i);
    cost = cost + 0.5*(z(:,i)'*Q*z(:,i) + UU(:,i)'*R*UU(:,i));

    Z(:,i+1) = AA*Z(:,i) + BB*UU(:,i); 
    z(:,i+1) = CC*Z(:,i+1);
  
end
cost = cost + 0.5*(z(:,i+1)'*Q_N*z(:,i+1));

%concatenating initial controls and outputs

if q > q_bar
    UU = [u_vec(:,q_bar:q-1), UU];
    z = [y(:,q_bar:q-1),z];
    
end