function out = coprod(vec, Coordinate)
n = size(vec, 3);
out = zeros(1,3,n);
if size(Coordinate,3) == 1
    for it = 1:n
        Coord = Coordinate(1,:,1);
        out(1,:,it) = vec(1,1,it)*Coord;
    end
else
    if size(Coordinate,3)~= n
        error('Size of vectors do not match');
    end
    for it = 1:n
        out(1,:,it) = vec(1,1,it)*Coordinate(1,:,it);
    end
end

end

