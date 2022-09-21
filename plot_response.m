function [] = plot_response(err_y_arma, err_y_OKID, t_steps, q)

fig=figure;
subplot(2,1,1)
plot(1:t_steps, err_y_arma(1,:),'Marker','.','MarkerSize',15,'Color','b', 'Linewidth',1,...
                    'DisplayName', 'ARMA');
hold on;
%plot(1:t_steps, err_y_open_loop(1,:),'--r', 'Linewidth',3, 'DisplayName', 'Open-loop markov');
plot(1:(t_steps-q), err_y_OKID(1,1:(t_steps-q)),'o','Color',...
          '[0.4660 0.6740 0.1880]', 'Linewidth',2, 'DisplayName', 'TV-OKID','MarkerSize',10);
ylabel('Error in output-1');
%ylim([-0.1, 0.4]);
legend('Location', 'SouthEast');
%legend();
subplot(2,1,2)
plot(1:t_steps, err_y_arma(2,:),'Marker','.','MarkerSize',15,'Color','b', 'Linewidth',1);
hold on;
%plot(1:t_steps, err_y_open_loop(2,:),'--r', 'Linewidth',3);
plot(1:(t_steps-q), err_y_OKID(2,1:(t_steps-q)),'o','Color',...
            '[0.4660 0.6740 0.1880]', 'Linewidth',2,'MarkerSize',10);
%ylim([-5e-12, 5e-12]);
xlabel('time steps');
ylabel('Error in output-2');

%{
set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters '/home/naveed/Dropbox/Research/Manuscripts/ACC23/plots/oscillator_q_4.pdf'
%}

%{
figure;
subplot(2,1,1)
plot(1:t_steps, err_y_arma(1,:),'b', 'Linewidth',3, 'DisplayName', 'ARMA');
hold on;
plot(1:t_steps, err_y_open_loop(1,:),'--r', 'Linewidth',3, 'DisplayName', 'Open-loop markov');
%plot(1:(t_steps-q-1), err_y_ABC(1,1:(t_steps-q-1)),':','Color',	'[0.4660 0.6740 0.1880]', 'Linewidth',3, 'DisplayName', 'TV-OKID');
ylabel('Error in output-1');
%ylim([-0.1, 0.4]);
legend('Location', 'NorthWest');
%legend();
subplot(2,1,2)
plot(1:t_steps, err_y_arma(2,:),'b', 'Linewidth',3);
hold on;
plot(1:t_steps, err_y_open_loop(2,:),'--r', 'Linewidth',3);
%plot(1:(t_steps-q-1), err_y_ABC(2,1:(t_steps-q-1)),':','Color',	'[0.4660 0.6740 0.1880]', 'Linewidth',3);
%ylim([-5e-12, 5e-12]);
xlabel('time steps');
ylabel('Error in output-2');
%}

end

