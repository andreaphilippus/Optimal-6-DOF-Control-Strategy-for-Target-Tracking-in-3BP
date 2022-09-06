function draw_video_CR3BP(t,xx,xx1,u_cl,xs,N,xd,plotMPC)
%% Draw Satellite Trajectory in Lyapunov Orbit around Lagrange Point (2D)                                                                                                                                                                                 
% Mark Batistich, Andrew Kim, Joseph Le, Philip Ra                                            
% AAE 568 - Applied Optimal Control and Estimation - Professor Inseok Hwang
%---------------------------------------------------------------------------------------------
%
% [ This program animates the satellite trajectory using data from casadi_CR3BP_RK4 in 2D. ]
%
%---------------------------------------------------------------------------------------------

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 3;
fontsize_labels = 14;

% initialize reference trajectory
x_r_1 = [];
y_r_1 = [];

m1 = 5.9722e24; % kg, mass of earth
m2 = 1.98847e30; % kg, mass of sun
mu = m1 / (m1 + m2); % mass ratio
L1 = [1 - nthroot(mu/3, 3), 0, 0]';
figure(500)
% Animate the motion
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55*0.7 1*0.7]);

fileName = append('animation-',datestr(datetime('now')),'.gif');
fileName = strrep(fileName, ' ', '.');
fileName = strrep(fileName, ':', '.');

for k = 1:size(xx,2)
    
    plot(L1(1),L1(2),'r*');
    hold on;
    plot(xd(1,:), xd(2,:), 'k-');
    x1 = xx(1,k,1); y1 = xx(2,k,1); z1 = xx(3,k,1);
    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];

    plot(x_r_1,y_r_1,'-b','linewidth',line_width);hold on % plot exhibited trajectory
    if plotMPC
        if k < size(xx,2) % plot prediction
            plot(xx1(1:N,1,k),xx1(1:N,2,k),'g-','linewidth',line_width)
        end
    end
    
    plot(x1, y1, 'k*', 'MarkerSize', 10); % plot position
   
    hold off
    ylabel('$y$-position','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$-position','interpreter','latex','FontSize',fontsize_labels)
    axis([0.985, 1.015 -0.010, 0.010]);
    axis equal
    box on;
    grid on
    xlim([0.985, 0.995]);
    ylim([-0.005; 0.005]);
    drawnow
    % for video generation
    if plotMPC
        lgd = legend("L1 point", "Reference Trajectory", "Spacecraft Trajectory", "MPC predicted trajectory", "Spacecraft Position");
    else
        lgd = legend("L1 point", "Reference Trajectory", "Spacecraft Trajectory", "Spacecraft Position");
    end
    %lgd.FontSize = 7;
    %lgd.Location = 'eastoutside';
    exportgraphics(gcf,fileName,'Append',true); % save animation as GIF
    F(k) = getframe(gcf); % to get the current frame
    
end
close(gcf)
