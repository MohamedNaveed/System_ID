function [x,y] = generate_response(x0, u_vec, sysd, ADD_NOISE)

t_steps = length(u_vec);

y = zeros(size(sysd.C,1),t_steps);
x = zeros(size(sysd.A,1),t_steps+1);
x(:,1) = x0;

for i = 1:t_steps
    
    x(:,i+1) = sysd.A*x(:,i) + sysd.B*u_vec(i);
    y(:,i) = sysd.C*x(:,i) + sysd.D*u_vec(i);
    
    if ADD_NOISE
        y(:,i) = y(:,i) + normrnd(0, 0.001,size(sysd.C,1),1);
    
end
    
end

