function x = trilateration1(ranges, stds, positions)
% ranges = vector of the range measurements
% stds = vector of the standard deviations of the range measurements
% positions = matrix of the beacon positions (xi,yi,zi) (each row correspond 
% to a beacon position: [x1,y1,z1; x2,y2,z2; ...; xn,yn,zn])

x = [];
nbRanges = length(ranges);
nbStds = length(stds);
nbPositions = size(positions,1);
dim = size(positions,2);

if (nbRanges ~= nbStds) || (nbRanges ~= nbPositions)
    disp('trilateration1: dim error');
    return;
end

if nbRanges < dim+1
    disp('trilateration1: not enough range measurements');
    return;
end

A = zeros(nbRanges-1,dim);
b = zeros(nbRanges-1,1);

for i=1:nbRanges-1
    for j=1:dim
        A(i,j) = positions(i+1,j) - positions(1,j);
        b(i) = b(i) + positions(i+1,j)^2 - positions(1,j)^2;
    end
    b(i) = b(i) + ranges(1)^2 - ranges(i+1)^2;
    b(i) = .5*b(i);
end

AtA = A'*A;
x = inv(AtA)*A'*b;

end