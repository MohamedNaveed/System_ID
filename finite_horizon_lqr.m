function [K] = finite_horizon_lqr(A,B,Q,R,Q_N,N)

% implementation of discrete-time finite horizon lqr (time-invarient)
% A, B system matrices. Q and R cost matrices. 
% Q_N terminal cost. N - horizon
nx = size(A,1);
nu = size(B,2);

K = zeros(nu,nx,N);

P = Q_N; %terminal condition.

for k = N:-1:1
    
    K(:,:,k) = -inv(R + B'*P*B)*B'*P*A;
    
    P = Q + A'*P*A + A'*P*B*K(:,:,k);
    
end

end


