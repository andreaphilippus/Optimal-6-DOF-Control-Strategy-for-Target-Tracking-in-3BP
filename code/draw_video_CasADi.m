function draw_video_CasADi(xx,xx1,N,xd,Var,figNum)

% Draws in dimensional units

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 2;
fontsize_labels = 14;

% initialize reference trajectory
x_r_1 = [];
y_r_1 = [];
z_r_1 = [];

L2 = [Var.L2x * Var.lstar, 0, 0]';

figure(figNum)

% Animate the motion
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55*0.7 1*0.7]);

fileName = append('animation-',datestr(datetime('now')),'.gif');
fileName = strrep(fileName, ' ', '.');
fileName = strrep(fileName, ':', '.');

% Extract Data

for i = 1:size(xd,3)
    [~,Rd(:,i)] = invSE3(xd(:,:,i));
    %[~,R(:,i)] = invSE3(xx(:,:,i));
end

Temp.lstar = Var.lstar;
Temp.vecSize = Var.vecSize * 3;

for k = 1:size(xx,3)
    
    % Plot L2 Point
    plot3(L2(1),L2(2),L2(3),'r*');  hold on;

    % Plot Reference Attitude
    DrawAttitude(k, xd, Temp, 'bk')

    % Plot Reference Trajectory
    plot3(Rd(1,:)*Var.lstar, Rd(2,:)*Var.lstar,Rd(3,:)*Var.lstar, 'k-')
    
    % Plot Exhibited Attitude
    DrawAttitude(k, xx, Temp, 'cl')

    % Plot Exhibited Trajectory
    x1 = xx(1,4,k) * Var.lstar;
    y1 = xx(2,4,k) * Var.lstar; 
    z1 = xx(3,4,k) * Var.lstar;

    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];
    z_r_1 = [z_r_1 z1];

    plot3(x_r_1,y_r_1,z_r_1,'-b','linewidth',line_width)

    % Plot Trajectory Prediction
    if k < size(xx,3) 
        plot3(xx1(1:N,10,k)*Var.lstar,xx1(1:N,11,k)*Var.lstar,...
            xx1(1:N,12,k)*Var.lstar, 'g-', 'linewidth', line_width)
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
    
    xlim([0.8 1.2]*Var.lstar);
    ylim([-0.2 0.2]*Var.lstar);
    zlim([-0.2 0.2]*Var.lstar);
    view([45 15])
    drawnow
    % for video generation
<<<<<<< Updated upstream
    lgd = legend("L2 point", "Reference Trajectory", ...
        "Spacecraft Trajectory", "MPC predicted trajectory", ...
        "Spacecraft Position");
=======
    %leg = legend([l2 ref_traj sc_traj mpc_traj sc_pos], "L2 point", ...
    %    "Reference Trajectory", "Spacecraft Trajectory", ...
    %    "MPC predicted trajectory", "Spacecraft Position");
>>>>>>> Stashed changes
    %lgd.FontSize = 7;
    %lgd.Location = 'eastoutside';
    set(gca, 'CameraPosition', [1.5 -0.4 0.4]*Var.lstar);
    exportgraphics(gcf,fileName,'Append',true); % save animation as GIF
    F(k) = getframe(gcf); % to get the current frame
    
end
close(gcf)
