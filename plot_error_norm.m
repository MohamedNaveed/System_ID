function [] = plot_error_norm(norm_err_arma, norm_err_OKID, norm_err_ABC,q, SAVE_PLOT)

t_steps = size(norm_err_arma,2);
font_size = 10;

plot_idx_info = q+1:t_steps;
plot_idx_okid = q+1:(t_steps-q);
clf;
fig=figure(1);
hold off;


semilogy(plot_idx_info, norm_err_arma(:,plot_idx_info),'Marker','.','MarkerSize',15,'Color','[0.4660 0.6740 0.1880]', 'Linewidth',1,...
                    'DisplayName', 'Information state');
hold on;


semilogy(plot_idx_okid, norm_err_OKID(:,plot_idx_okid),'o','Color',...
          'r', 'Linewidth',1, 'DisplayName', 'TV-OKID','MarkerSize',6);

semilogy(plot_idx_okid, norm_err_ABC(:,plot_idx_okid),'o','MarkerSize',6,'Color','b', 'Linewidth',1,...
            'DisplayName', 'TV-OKID (without obs.)');

hold on; 
ylabel('$||Error||_1 $','Interpreter','latex','FontSize',font_size);
xlabel('time steps','FontSize',font_size);
xlim([q,t_steps])
ylim([1e-18,1e-12]);
yticks([1e-18, 1e-16, 1e-14, 1e-12]);
%h = legend('Location', 'SouthEast');
h = legend();
set(h,'FontSize',font_size,'Position',[0.6 0.2410 0.1 0.1]);


if SAVE_PLOT
    set(fig,'Units','inches');
    fig.Position = [80,80,4.5,2.5];
    screenposition = get(fig,'Position');
    set(fig,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters '/home/naveed/Dropbox/ACC23/plots/oscillator_norm_zeroini_q=4_with_openloop.pdf'
end
end


