function [x,y] = generate_response_oscillator(x0, u_vec, n, nz, del_t)


t_steps = size(u_vec,2);

y = zeros(nz,t_steps); %output is recorded from t=1
x = zeros(n,t_steps+1);
x(:,1) = x0;

for i = 1:t_steps
    
    sysd = oscillator((i-1)*del_t);
    B_u = sysd.B*u_vec(:,i);
    
    x(:,i+1) = sysd.A*x(:,i) + B_u;
    
    y(:,i) = sysd.C*x(:,i) + sysd.D*u_vec(:,i);
    
   
end
    
end

