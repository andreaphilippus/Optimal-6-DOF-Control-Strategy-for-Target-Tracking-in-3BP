function DrawAttitude(index,g,Var,bkcl)

    [NB, R] = invSE3(g(:,:,index)); 

    R = Var.lstar * R;
    NB = Var.vecSize * NB;
    xhat = NB(:,1); yhat = NB(:,2); zhat = NB(:,3);

    if bkcl == 'cl'
        quiver3(R(1), R(2), R(3), xhat(1), xhat(2), xhat(3),...
            'c', 'LineWidth', 1)
        quiver3(R(1), R(2), R(3), yhat(1), yhat(2), yhat(3),...
            'g', 'LineWidth', 1)
        quiver3(R(1), R(2), R(3), zhat(1), zhat(2), zhat(3),...
            'm', 'LineWidth', 1)
    else
        quiver3(R(1), R(2), R(3), xhat(1), xhat(2), xhat(3),...
            'k', 'LineWidth', 1)
        quiver3(R(1), R(2), R(3), yhat(1), yhat(2), yhat(3),...
            'k', 'LineWidth', 1)
        quiver3(R(1), R(2), R(3), zhat(1), zhat(2), zhat(3),...
            'k', 'LineWidth', 1)
    end
    
    
end