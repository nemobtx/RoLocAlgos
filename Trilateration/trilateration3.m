function x = trilateration3(ranges, stds, positions)
% ranges = vector of the range measurements
% stds = vector of the standard deviations of the range measurements
% positions = matrix of the beacon positions (xi,yi,zi) (each row correspond 
% to a beacon position: [x1,y1,z1; x2,y2,z2; ...; xn,yn,zn])

nbRanges = length(ranges);
nbStds = length(stds);
nbPositions = size(positions,1);
dim = size(positions,2);
x = [];

if (nbRanges ~= nbStds) || (nbRanges ~= nbPositions)
    disp('trilateration3: dim error');
    return;
end

if nbRanges < dim+1
    disp('trilateration3: not enough ranges');
    return;
end

if dim == 2
    A = zeros(nbRanges,4);
    for i=1:nbRanges
        A(i,1) = -2*positions(i,1);
        A(i,2) = -2*positions(i,2);
        A(i,3) = 1;
        A(i,4) = positions(i,1)^2 + positions(i,2)^2 - ranges(i)^2;
    end
    
    [U,D,V] = svd(A);
    
    x = V(:,end);
    x = x/x(end);
    x = x(1:2);
elseif dim == 3
    disp('trilateration3: dim==3 not implemented yet!');
else
    disp('trilateration3: error, dim not 2 or 3!');
end


end