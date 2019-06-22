function [newx, newy] = transform(x, y)
    k = [1.2 0.3; 0.2 1.3];
    new = k*[x, y]';
    newx = new(1);
    newy = new(2);
end

