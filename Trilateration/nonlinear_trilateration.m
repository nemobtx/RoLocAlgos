function x = nonlinear_trilateration(x0, rs, bs)
% x0=(x0,y0) or (x0,y0,z0): initial position estimate
% rs=[r1,r2,...,rn]: vector of the range measurements
% bs=[x1,y1,z1; x2,y2;z2; ...; xn,yn,zn]: matrix of the beacon positions 
% (one beacon per row)

min_e = 1e-6;
max_iter = 400;
max_iter_follow = 10;
diff_e = 1e-6;

dim = length(x0);
nbBeacons = size(bs,1);
nbRanges = length(rs);

if nbBeacons ~= nbRanges
    disp('nonlinear_trilateration: nbBeacons ~= nbRanges');
    return;
end

e = zeros(nbBeacons,1);
J = zeros(nbBeacons,dim);
x = x0;

iter = 1;
iter_follow = 0;
sum_e = 1e6;
sum_e_old = 1e6;
while abs(sum_e)>min_e && iter<max_iter && iter_follow<max_iter_follow
    for i=1:nbBeacons
        pb = bs(i,1:dim);
        dx = pb-x;
        d = sqrt(sum(dx.^2));
        e(i) = d - rs(i);
        J(i,:) = -dx/d;
    end
    gp = J'*e;
    gpp = J'*J;
    delta = -inv(gpp)*gp;
    x = x + delta';
    sum_e = sum(e);
    iter = iter + 1;
    diff = sum_e - sum_e_old;
    if abs(diff) < diff_e
        iter_follow = iter_follow + 1;
    else
        iter_follow = 0;
    end
    sum_e_old = sum_e;
    %disp(['iter=',num2str(iter),',sum_e=',num2str(sum_e),...
    %    ',iter_follow=',num2str(iter_follow),',x=',num2str(x)]);
end


end
