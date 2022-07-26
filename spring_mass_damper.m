function [sysd] = spring_mass_damper()

% simple spring mass damper system with position output
% discrete time system. 

k = 1;
c = 0.01;
M = 1;
Ts = 0.2; %sampling time. (s)

%continuous time system
A = [0, 1; -k/M -c/M];
B = [0; 1/M];
%C = [1 0;-k/M -c/M];
C = [1 0];
%D = [0;1/M];
D = 0;

sysc = ss(A,B,C,D);
%sysc = rss(2); %random system. 

sysd = c2d(sysc,Ts,'zoh');

end

