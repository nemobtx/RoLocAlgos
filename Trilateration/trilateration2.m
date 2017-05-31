function [x, C] = trilateration2(ranges, stds, positions)
% ranges = vector of the range measurements
% stds = vector of the standard deviations of the range measurements
% positions = matrix of the beacon positions (xi,yi,zi) (each row correspond 
% to a beacon position: [x1,y1,z1; x2,y2,z2; ...; xn,yn,zn])


nbRanges = length(ranges);
nbStds = length(stds);
nbPositions = size(positions,1);
dim = size(positions,2);

x = [];
C = [];

if (nbRanges ~= nbStds) || (nbRanges ~= nbPositions)
    disp('trilateration2: dim error');
    return;
end

if nbRanges < dim+1
    disp('trilateration2: not enough ranges');
    return;
end

A = zeros(nbRanges,dim+1);
W = eye(nbRanges, nbRanges);
invW = eye(nbRanges, nbRanges);
R = eye(nbRanges, nbRanges);
b = zeros(nbRanges,1);

for i=1:nbRanges
    si = 0;
    for j=1:dim
        A(i,j) = 2*positions(i,j);
        si = si + positions(i,j)^2;
    end
    W(i,i) = stds(i)^2;
    invW(i,i) = 1/W(i,i);
    R(i,i) = -2*ranges(i);
    A(i,dim+1) = -1;
    b(i) = si - ranges(i)^2;
end

G = A'*invW*A;
d = A'*invW*b;
AAt = A*A';
invA = A'*inv(AAt);
invG = inv(G);

x = invG*d;
x = x(1:dim);

T = invA*R;
C = T*W*T';
%disp(C);
C = C(1:dim,1:dim);
for i=1:dim
    for j=1:dim
        if isnan(C(i,j)) || ~isfinite(C(i,j)) || C(i,j)<0
            C(i,j) = 100;
        end
    end
end

end





