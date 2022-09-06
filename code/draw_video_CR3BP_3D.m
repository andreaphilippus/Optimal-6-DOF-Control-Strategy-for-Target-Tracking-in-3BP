function draw_video_CR3BP_3D(t,xx,xx1,u_cl,xs,N,xd,plotMPC)
%% Draw Satellite Trajectory in Lyapunov Orbit around Lagrange Point (3D)                                                                                              
% Mark Batistich, Andrew Kim, Joseph Le, Philip Ra                                            
% AAE 568 - Applied Optimal Control and Estimation - Professor Inseok Hwang
%---------------------------------------------------------------------------------------------
%
% [ This program animates the satellite trajectory using data from casadi_CR3BP_RK4 in 3D. ]
%
%---------------------------------------------------------------------------------------------
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 2;
fontsize_labels = 14;

% initialize reference trajectory
x_r_1 = [];
y_r_1 = [];
z_r_1 = [];

m1 = 4902.800582147765; % kg, mass of Earth
m2 = 398600.4328969393; % kg, mass of Sun
mu = m1 / (m1 + m2); % mass ratio
L1 = [1 + nthroot(mu/3, 3), 0, 0]';
figure(500)
% Animate the motion
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55*0.7 1*0.7]);

fileName = append('animation-',datestr(datetime('now')),'.gif');
fileName = strrep(fileName, ' ', '.');
fileName = strrep(fileName, ':', '.');

for k = 1:size(xx,2)
    
    plot3(L1(1),L1(2),0,'r*');
    hold on;
    plot3(xd(1,:), xd(2,:),xd(3,:), 'k-');
    x1 = xx(1,k,1); y1 = xx(2,k,1); z1 = xx(3,k,1);
    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];
    z_r_1 = [z_r_1 z1];

    plot3(x_r_1,y_r_1,z_r_1,'-b','linewidth',line_width);hold on % plot exhibited trajectory
    if plotMPC == 1
        if k < size(xx,2) % plot prediction
            plot3(xx1(1:N,1,k),xx1(1:N,2,k),xx1(1:N,3,k),'g-','linewidth',line_width)
        end
    end
    
    plot3(x1, y1, z1, 'k*', 'MarkerSize', 10); % plot position
   
    hold off
    ylabel('$y$-position','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$-position','interpreter','latex','FontSize',fontsize_labels)
    zlabel('$z$-position','interpreter','latex','FontSize',fontsize_labels)

    %axis([0.985, 1.015 -0.010, 0.010]);
    axis equal
    box on;
    grid on
    
    xlim([0.8 1.2]); ylim([-0.2 0.2]); zlim([-0.2 0.2]);
    view([45 15])
    drawnow
    % for video generation
    if plotMPC == 1
         lgd = legend("L2 point", "Reference Trajectory", "Spacecraft Trajectory", "MPC predicted trajectory", "Spacecraft Position");
    else
         lgd = legend("L2 point", "Reference Trajectory", "Spacecraft Trajectory", "Spacecraft Position");
    end
    %lgd.FontSize = 7;
    %lgd.Location = 'eastoutside';
    set(gca, 'CameraPosition', [1.5 -0.4 0.4]);
    exportgraphics(gcf,fileName,'Append',true); % save animation as GIF
    F(k) = getframe(gcf); % to get the current frame
    
end
close(gcf)
