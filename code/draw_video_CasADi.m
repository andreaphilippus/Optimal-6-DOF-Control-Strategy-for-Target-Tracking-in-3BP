function draw_video_CasADi(xx,xx1,xxun1,N,xd,Var)

% Draws in dimensional units

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)
gray = [128 128 128]/255;
graydark = [90 90 90]/255;

line_width = 2;
fontsize_labels = 14;

% initialize exhibited trajectory
x_r_1 = [];
y_r_1 = [];
z_r_1 = [];

% initialize uncontrolled trajectory
x_r_u_1 = [];
y_r_u_1 = [];
z_r_u_1 = [];

L2 = [Var.L2x, 0, 0]';

figure()

% Animate the motion
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55*0.7 1*0.7]);

fileName = append('animation-',datestr(datetime('now')),'.gif');
fileName = strrep(fileName, ' ', '.');
fileName = strrep(fileName, ':', '.');

% Extract Data

for i = 1:size(xd,3)
    [~, Rd(:,i)] = invSE3(xd(:,:,i));
    %[NB, R] = invSE3(xx(:,:,i));
end

Temp.lstar = Var.lstar;
Temp.vecSize = Var.vecSize;

for k = 1:size(xx,3)
    
    
    % Plot L2 Point
    %plot3(L2(1),L2(2),L2(3),'ro');  

    % Plot Reference Attitude
    %DrawAttitude(k, xd, Temp, 'bk')

    % Plot Reference Trajectory
    plot3(Rd(1,:), Rd(2,:),Rd(3,:), 'k-'); hold on

    % Plot Uncontrolled Attitude
    DrawAttitude(k, xxun1, Temp, 'bk')

    % Plot Uncontrolled Trajectory
    xu1 = xxun1(1,4,k);
    yu1 = xxun1(2,4,k);
    zu1 = xxun1(3,4,k);

    x_r_u_1 = [x_r_u_1 xu1];
    y_r_u_1 = [y_r_u_1 yu1];
    z_r_u_1 = [z_r_u_1 zu1];

    plot3(x_r_u_1, y_r_u_1, z_r_u_1,'--', 'Color', gray, 'LineWidth', line_width)

    % Uncontrolled S/C Position
    plot3(xu1, yu1, zu1, 'k*', 'MarkerSize',10)
    
    % Plot Exhibited Attitude
    DrawAttitude(k, xx, Temp, 'cl')

    % Plot Exhibited Trajectory
    x1 = xx(1,4,k);
    y1 = xx(2,4,k); 
    z1 = xx(3,4,k);

    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];
    z_r_1 = [z_r_1 z1];

    plot3(x_r_1,y_r_1,z_r_1,'-b','linewidth',line_width)

    % Plot Trajectory Prediction
    if k < size(xx,3) 
        plot3(xx1(1:N,10,k),xx1(1:N,11,k),...
            xx1(1:N,12,k), 'g-', 'linewidth', line_width)
    end
    
    % Controlled S/C Position
    plot3(x1, y1, z1, 'r*', 'MarkerSize', 10);
    
    
    %Plot Moon
    plot3(1-Var.mu, 0, 0, '.', 'MarkerSize', 70, 'Color',graydark)

    hold off
    ylabel('$y$','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$','interpreter','latex','FontSize',fontsize_labels)
    zlabel('$z$','interpreter','latex','FontSize',fontsize_labels)

    %axis([0.985, 1.015 -0.010, 0.010]);
    axis equal
    box on;
    grid on

    %
    xlim([0.8 1.2]);
    ylim([-0.1 0.1]);
    zlim([-0.1 0.1]);
    view([45 15])
    %}
    campos([1.2,-0.1,0.1])
    drawnow
    % for video generation
    lgd = legend("Reference Trajectory", "","","","Uncontrolled Trajectory",...
        "Uncontrolled S/C Position", "","","",...
        "Controlled Trajectory", "MPC Prediction Window", ...
        "Controlled S/C Position", "Moon");
    %lgd.FontSize = 7;
    %lgd.Location = 'eastoutside';
    set(gca, 'CameraPosition', [1.2 -0.1 0.1]);
    exportgraphics(gcf,fileName,'Append',true); % save animation as GIF
    F(k) = getframe(gcf); % to get the current frame
    

end

