function [c, r] = minimalCircle(points)
%MINIMALCIRCLE Get minimal area enclosing circle
%   [c,r] = MINIMALCIRCLE(points) Radius 'r' and center 'c' of the minimal
%   area enclosing circle around the given points are returned. The input
%   parameter 'points' is expected to be a complex vector.
%
%   Algorithm by Emo Welzl to find the minimal covering circle for given 
%   points in expected linear time.
%   Welzl E. (1991) Smallest enclosing disks (balls and ellipsoids). 
%   In: Maurer H. (eds) New Results and New Trends in Computer Science. 
%   Lecture Notes in Computer Science, vol 555. Springer, Berlin, 
%   Heidelberg

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

assert(isvector(points),'The input must be a complex vector');
points = points(randperm(length(points)));
[c, r] = mCircle(length(points), []);
    function [c, r] = mCircle(n, R)
        if n==0 || length(R)==3
            [c, r] = trivial(R);
        else
            p = points(n);
            [c, r] = mCircle(n-1, R);
            if abs(p-c)>r
                [c, r] = mCircle(n-1, [R p]);
                points(1:n) = points([n,1:n-1]);
            end
        end
    end
end

function [c, r] = trivial(points)
nPoints = length(points);
switch nPoints
    case 0
        c = 0;
        r = -1;
    case 1
        c = points;
        r = 0;
    case 2
        d = points(1)-points(2);
        c = points(1)-d/2;
        r = abs(d)/2;
    case 3
        x = real(points); y = imag(points);
        A = 2*[x(2)-x(1) y(2)-y(1); x(3)-x(1) y(3)-y(1)];
        b = [x(2)^2-x(1)^2+y(2)^2-y(1)^2; x(3)^2-x(1)^2+y(3)^2-y(1)^2];
        c = A\b;
        c = c(1) + 1i*c(2);
        r = abs(c-points(1));
end    
end
