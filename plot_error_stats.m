function [] = plot_error_stats(mean_err_y_arma, std_err_y_arma, ...
    mean_err_y_OKID, std_err_y_OKID, mean_err_y_ABC, std_err_y_ABC, q, SAVE_PLOT)

t_steps = size(mean_err_y_arma,2);
font_size = 14;

plot_idx_info = q+1:t_steps;
plot_idx_okid = q+1:(t_steps);

fig=figure;
subplot(2,1,1)


plot(plot_idx_info, mean_err_y_arma(1,plot_idx_info),'Marker','.','MarkerSize',15,'Color','b', 'Linewidth',1,...
                    'DisplayName', 'Information state');
hold on;
x = q+1:numel(mean_err_y_arma(1,plot_idx_info))+q;
curve_1 = mean_err_y_arma(1,plot_idx_info) + 3*std_err_y_arma(1,plot_idx_info);
curve_2 = mean_err_y_arma(1,plot_idx_info) - 3*std_err_y_arma(1,plot_idx_info);
inbetween = [curve_1, fliplr(curve_2)];
x2 = [x, fliplr(x)];
%fill(x2,inbetween,'b','FaceAlpha',0.3,'HandleVisibility','off');
    

plot(plot_idx_okid, mean_err_y_OKID(1,plot_idx_okid),'o','Color',...
          '[0.4660 0.6740 0.1880]', 'Linewidth',2, 'DisplayName', 'OKID','MarkerSize',8);

x = q+1:numel(mean_err_y_OKID(1,plot_idx_okid))+q;
curve_1 = mean_err_y_OKID(1,plot_idx_okid) + 3*std_err_y_OKID(1,plot_idx_okid);
curve_2 = mean_err_y_OKID(1,plot_idx_okid) - 3*std_err_y_OKID(1,plot_idx_okid);
inbetween = [curve_1, fliplr(curve_2)];
x2 = [x, fliplr(x)];
%fill(x2,inbetween,[0.4660 0.6740 0.1880],'FaceAlpha',0.3,'HandleVisibility','off');

plot(plot_idx_okid, mean_err_y_ABC(1,plot_idx_okid),'o','MarkerSize',10,'Color','r', 'Linewidth',2,...
                    'DisplayName', 'A,B,C,D');
hold on;
ylabel('Error in output-1','FontSize',font_size);
xlim([q,t_steps])
%ylim([-1e-3,4e-3]);
%h = legend('Location', 'SouthWest');
%h =legend();
%set(h,'FontSize',font_size);


subplot(2,1,2)
plot(plot_idx_info, mean_err_y_arma(2,plot_idx_info),'Marker','.','MarkerSize',15,'Color','b', 'Linewidth',1);
hold on;

x = q + 1 : numel(mean_err_y_arma(2,plot_idx_info)) + q;
curve_1 = mean_err_y_arma(2,plot_idx_info) + 3*std_err_y_arma(2,plot_idx_info);
curve_2 = mean_err_y_arma(2,plot_idx_info) - 3*std_err_y_arma(2,plot_idx_info);
inbetween = [curve_1, fliplr(curve_2)];
x2 = [x, fliplr(x)];
%fill(x2,inbetween,'b','FaceAlpha',0.3,'HandleVisibility','off');

plot(plot_idx_okid, mean_err_y_OKID(2,plot_idx_okid),'o','Color',...
            '[0.4660 0.6740 0.1880]', 'Linewidth',2,'MarkerSize',10);
        
x = q + 1 : numel(mean_err_y_OKID(2,plot_idx_okid)) + q;
curve_1 = mean_err_y_OKID(2,plot_idx_okid) + 3*std_err_y_OKID(2,plot_idx_okid);
curve_2 = mean_err_y_OKID(2,plot_idx_okid) - 3*std_err_y_OKID(2,plot_idx_okid);
inbetween = [curve_1, fliplr(curve_2)];
x2 = [x, fliplr(x)];
%fill(x2,inbetween,[0.4660 0.6740 0.1880],'FaceAlpha',0.3,'HandleVisibility','off');
%ylim([-5e-12, 5e-12]);

plot(plot_idx_okid, mean_err_y_ABC(2,plot_idx_okid),'o','MarkerSize',10,'Color','r', 'Linewidth',2);
                
xlim([q,t_steps])
xlabel('time steps','FontSize',font_size);
ylabel('Error in output-2','FontSize',font_size);

if SAVE_PLOT
    set(fig,'Units','inches');
    screenposition = get(fig,'Position');
    set(fig,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters '/home/naveed/Dropbox/Research/Manuscripts/ACC23/plots/oscillator_q_4_mean_error_nonzero_ini.pdf'
end
end


