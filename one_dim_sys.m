function [sysd] = one_dim_sys()

% simple spring mass damper system with position output
% discrete time system. 



%continuous time system
A = 0.1;
B = 1;
C = 1;
D = 0;

Ts = 1/10; %sampling time. (s)
sysd = ss(A,B,C,D,Ts);

end

