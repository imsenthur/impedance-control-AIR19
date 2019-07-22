function [newx, newy] = transform(x, y)
    k = [4 3; 3 0];
    new = k*[x, y]';
    newx = new(1);
    newy = new(2);
end

