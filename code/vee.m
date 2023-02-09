function out = vee(in)

    out = [CrossProd(in(1:3)) in(4:6); zeros(1,4)];
end