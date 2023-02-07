% initial condition decay experiment. 

clc;clear;

load('data/cartpole_initial_condition_decay_100.mat');

Z_actual = Task.Ck*x_actual;
error_z = (Z_actual - Z_NORM);
error_z_norm = abs(error_z./Z_NORM);

time_idxs = 0:horizon;
%%
fig = figure(1); 
outputs = 1:2;
semilogy(time_idxs, error_z_norm(1,:),'b', 'LineWidth',3,'DisplayName','Cart pos');
hold on;
semilogy(time_idxs, error_z_norm(2,:),'r', 'LineWidth',3,'DisplayName','Pole pos');
title('Error in response');
xlabel('time steps');
ylabel('Normalized Error');
legend();
grid on;
%%
fig = figure(2);
plot(time_idxs, Z_actual(:,:),'r', 'LineWidth',2,'DisplayName', 'Actual');
hold on; 
plot(time_idxs, Z_NORM(:,:),'--b', 'LineWidth',2,'DisplayName', 'True');
xlabel('time steps');
ylabel('state');
title('Response');