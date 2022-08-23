function [sysd] = oscillator(tk)

% simple spring mass damper system with position output
% discrete time system. 


Ts = 1.0/10.0; %sampling time. (s)

%discrete time system
tau = sin(10*tk); 
tau_dash = cos(10*tk);
Kt = [4 + 3*tau, 1; 1, 7 + 3*tau_dash];
Ac = [zeros(2,2) eye(2); -Kt zeros(2,2)];

A = expm(Ac*Ts);
B = [1, 0;1, -1;0, 1;-1, 0];
C = [1 0 1 0.2;1 -1 0 -0.5];
D = 0.1*eye(2);

sysd.A = A;
sysd.B = B;
sysd.C = C;
sysd.D = D;
sysd.Ts = Ts;

end

