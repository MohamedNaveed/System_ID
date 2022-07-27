function [x,y] = generate_response(x0, u_vec, sysd, ADD_PROC_NOISE,ADD_MEAS_NOISE)

t_steps = length(u_vec);

y = zeros(size(sysd.C,1),t_steps);
x = zeros(size(sysd.A,1),t_steps+1);
x(:,1) = x0;

for i = 1:t_steps
    B_u = sysd.B*u_vec(i);
    x(:,i+1) = sysd.A*x(:,i) + B_u;
    
    if ADD_PROC_NOISE
        %Q = 0.05*diag(abs(B_u));
        Q = diag([0.0242, 3.592, 0.0534, 1.034, 0.0226, 0.2279])*10^-4;
        x(:,i+1) = x(:,i+1) + mvnrnd(zeros(size(B_u)), Q, 1)';
    end
    
    y(:,i) = sysd.C*x(:,i) + sysd.D*u_vec(i);
    
    if ADD_MEAS_NOISE
        
        %R = 0.05*diag(abs(y(:,i)));
        R = diag([2.785, 2.785])*(10^-2);
        y(:,i) = y(:,i) + mvnrnd(zeros(size(y(:,i))), R, 1)';
    end
end
    
end

