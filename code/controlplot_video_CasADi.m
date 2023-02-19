function controlplot_video_CasADi(t,u_cl,uu1,N,Var)


set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)
gray = [128 128 128]/255;
graydark = [90 90 90]/255;

line_width = 2;
fontsize_labels = 14;

% initialize control history
t1_1 = [];
t2_1 = [];
t3_1 = [];
fx_1 = [];
fy_1 = [];
fz_1 = [];

figure()

% Animate the motion
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55*0.7 1*0.7]);

fileName = append('control_vid-',datestr(datetime('now')),'.gif');
fileName = strrep(fileName, ' ', '.');
fileName = strrep(fileName, ':', '.');

for k = 1:size(u_cl,1)

    x1 = u_cl(k,4);
    y1 = u_cl(k,4);
    z1 = u_cl(k,4);
    t1 = u_cl(k,1);
    t2 = u_cl(k,2);
    t3 = u_cl(k,3);

    fx_1 = [fx_1 x1];
    fy_1 = [fy_1 y1];
    fz_1 = [fz_1 z1];
    t1_1 = [t1_1 t1];
    t2_1 = [t2_1 t2];
    t3_1 = [t3_1 t3];
    
    subplot(611)
    
    plot(t(1:k), fx_1,'b'); grid on; hold on
    if k < size(u_cl,1) 
        plot(t(k:k+N-1), uu1(:,4,k), 'g-', 'linewidth', line_width)
    end
    ylim([-1e-3 1e-3])
    xlim([0 7])
    ylabel('$f_x$','interpreter','latex','FontSize',fontsize_labels)
    %legend('Control History', 'Control Prediction')
    hold off

    subplot(612)
    
    plot(t(1:k), fy_1,'b'); grid on; hold on
    if k < size(u_cl,1) 
        plot(t(k:k+N-1), uu1(:,5,k), 'g-', 'linewidth', line_width)
    end
    ylim([-1e-3 1e-3])
    xlim([0 7])
    ylabel('$f_y$','interpreter','latex','FontSize',fontsize_labels)
    hold off

    subplot(613)
    
    plot(t(1:k), fz_1, 'b'); grid on; hold on
    if k < size(u_cl,1) 
        plot(t(k:k+N-1), uu1(:,6,k), 'g-', 'linewidth', line_width)
    end
    ylim([-1e-3 1e-3])
    xlim([0 7])
    ylabel('$f_z$','interpreter','latex','FontSize',fontsize_labels)
    hold off

    subplot(614)
    
    plot(t(1:k), t1_1, 'b'); grid on; hold on
    if k < size(u_cl,1) 
        plot(t(k:k+N-1), uu1(:,1,k), 'g-', 'linewidth', line_width)
    end
    ylim([-1e-4 1e-4])
    xlim([0 7])
    ylabel('$\tau_1$','interpreter','latex','FontSize',fontsize_labels)
    hold off

    subplot(615)
     
    plot(t(1:k), t2_1, 'b'); grid on; hold on
    if k < size(u_cl,1) 
        plot(t(k:k+N-1), uu1(:,2,k), 'g-', 'linewidth', line_width)
    end
    ylim([-1e-4 1e-4])
    xlim([0 7])
    ylabel('$\tau_2$','interpreter','latex','FontSize',fontsize_labels)
    hold off

    subplot(616)
    
    plot(t(1:k), t3_1, 'b'); grid on; hold on
    if k < size(u_cl,1) 
        plot(t(k:k+N-1), uu1(:,3,k), 'g-', 'linewidth', line_width)
    end
    ylim([-1e-4 1e-4])
    xlim([0 7])
    ylabel('$\tau_3$','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$t$','interpreter','latex','FontSize',fontsize_labels)

    hold off

    drawnow
    % for video generation

    %lgd.FontSize = 7;
    %lgd.Location = 'eastoutside';
    exportgraphics(gcf,fileName,'Append',true); % save animation as GIF
    F(k) = getframe(gcf); % to get the current frame

end
