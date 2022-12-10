function out = vee(wedge)

    eps = invSO3(wedge(1:3,1:3));
    delta = wedge(1:3,4);

    out = [eps; delta];
end