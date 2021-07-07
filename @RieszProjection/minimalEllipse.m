function E = minimalEllipse(points)
%MINIMALELLIPSE get minimal area enclosing ellipse
%   E = MINIMALILLIPSE(points) returns the parameters [c1 c2 c3 c4 c5] for 
%   the parameterization of the minimal area ellipse enclosing the given
%   points of the form: 
%   E(phi) = c1 + c2*1i + c3*exp(1i*(phi+c5)) + c4*exp(-1i*(phi-c5))
%   The input argument points is expected to be either a complex vector or
%   a real matrix of shape (n,2).
%
%   Algorithm by Emo Welzl to find the minimal covering ellipse for given 
%   points in expected linear time. 
%   Welzl E. (1991) in 
%   Maurer H. (eds) New Results and New Trends in Computer Science. 
%   Lecture Notes in Computer Science, vol 555. Springer
%   with primitives from G??rtner et al. in 
%   Information Processing Letters 68 (1998) 33-38 
%   and Rublev et al. in
%   Ukrainian Mathematical Journal, Vol. 50, No. 7 (1998).

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

assert(size(points,2)==2 || isvector(points),...
    'The input must either be of size (n,2) or a complex vector')
if ismember(1,size(points)), points = [real(points(:)) imag(points(:))]; end
M = max(points(:)); points = points/M;
points = points(randperm(size(points,1)),:);
E = me(mEllipse(size(points,1), []));
E(1:4) = E(1:4)*M;
    % If 3 or 5 points are at the boundary E has the form (r,s,t,u,v,w)
    % defining the coefficients of the ellipse equation. Otherwise only
    % the points at the boundary are returned. 
    function E = mEllipse(n,R)
        if n==0 || size(R,2)==4
            E = me_(R);
        else
            q = points(n,:);
            E = mEllipse(n-1,R);
            [outside,E_] = test(E,q,n,R);
            if outside
                E = E_;
                points(1:n,:) = points([n,1:n-1],:);
            end
        end
    end
    % Test if q is outside the ellipse defined by either by points on its
    % boundary or the coefficients of the ellipse equation.
    function [outside,C] = test(E,q,n,R)
        k = size(E,1); C = []; outside = true;
        if size(E,2)==6, k = 5; end
        if k==1
            outside = any(E~=q);
        elseif k==2
            t = (q-E(1,:))./(E(2,:)-E(1,:));
            outside = t(1)~=t(2) || abs(t(1)-1/2)>1/2;
        elseif k==3
            outside = (q-E(3,:))*E(1:2,:)*(q-E(3,:)).'-1>0;
        elseif k==4
            qxy = [q.^2 2*q(1)*q(2) 2*q 1];
            ch = convhull(E,'simplify',true); % polygon
            p1 = E(ch(1),:); p2 = E(ch(2),:); 
            p3 = E(ch(3),:); p4 = E(ch(4),:);
            % the minimal ellipse is a linear combination of C1 and C2
            v = p1-p2; w = p3-p4; c = [1/2 -1/2];
            a = p2(1)*p1(2)-p1(1)*p2(2); b = p4(1)*p3(2)-p3(1)*p4(2);
            C1 = [a*b c.*(v*b+w*a) -1/2*(w(1)*v(2)+v(1)*w(2)) v.*w];
            v = p2-p3; w = p4-p1;
            a = p3(1)*p2(2)-p2(1)*p3(2); b = p1(1)*p4(2)-p4(1)*p1(2);
            C2 = [a*b c.*(v*b+w*a) -1/2*(w(1)*v(2)+v(1)*w(2)) v.*w];
            C1 = C1(end:-1:1); C2 = C2(end:-1:1);
            lambda_0 =  dot(C2,qxy); mu_0 = -dot(C1,qxy);
            C = lambda_0*C1+mu_0*C2; if C(1)<0, C = -1*C; end % normalize
            d = (C(1)*C(2)-C(3)^2);
            if d<=0
                alpha = C1(1)*C1(2)-C1(3)^2; gamma = C2(1)*C2(2)-C2(3)^2;
                beta = C1(1)*C2(2)+C2(1)*C1(2)-2*C1(3)*C2(3);
                outside = dot((2*gamma-beta)*C1+(2*alpha-beta)*C2,qxy)>0;
            else
                rho = lambda_0*C1(1)+mu_0*C2(1);
                dC0 = [0 C1(1)*C2(2:end)-C2(1)*C1(2:end)];
                dd = dC0(1)*C(2)+C(1)*dC0(2)-2*C(3)*dC0(3);
                vu = C(5:-1:4); sr = C(2:-1:1);
                dvu = dC0(5:-1:4); dsr = dC0(2:-1:1);
                Z_ = C(4:5).*sr-vu*C(3); Z = dot(C(4:5),Z_);
                dZ_ = dC0(4:5).*sr+C(4:5).*dsr-vu*dC0(3)-dvu.*C(3);
                dZ = dC0(4)*Z_(1)+C(4)*dZ_(1)+dC0(5)*Z_(2)+C(5)*dZ_(2);
                outside = rho*(3*dd*Z+d*(2*d*dC0(end)-dd*C(end)-2*dZ))>0;
            end
        elseif k==5
            outside = dot(E,[q.^2 2*q(1)*q(2) 2*q 1])>0;
        end
        if outside&&(size(R,1)~=4), C = mEllipse(n-1,[R;q]); end
    end
end

% return parameters [c1 c2 c3 c4 c5] for the parameterization:
% E(phi) = c1 + c2*1i + c3*exp(1i*(phi+c5)) + c4*exp(-1i*(phi-c5))
function E = me(p)
nP = size(p,1); if size(p,2)==6, nP = 5; end
switch nP
    case 0
        E = zeros(1,5);
    case 1
        E = [p 0 0 0];
    case 2
        v = p(1,:)-p(2,:); phi = atan(v(2)/v(1));
        E = [p(2,:)+v/2 norm(v)/4 norm(v)/4 phi];
    case 3 % the ellipse has center form
        E = {p(1:2,:) p(3,:)};
    case 4 % get the ellipse in center form
        ch = convhull(p,'simplify',true); % define quadrangle
        a = p(ch(1),:); b = p(ch(2),:); c = p(ch(3),:); d = p(ch(4),:);
        % evaluate the intersection point between diagonals
        c1 = a(1)*c(2)-a(2)*c(1); c2 = b(1)*d(2)-b(2)*d(1);
        c3 = (a(1)-c(1))*(b(2)-d(2))-(a(2)-c(2))*(b(1)-d(1));
        o = (c1*(b-d)-(a-c)*c2)/c3;
        % define and solve cubic polynomial
        k = sqrt(norm(c-o)/norm(o-a)); n = sqrt(norm(o-b)/norm(d-o));
        I = (n^2*k^2+1)*(n^2+k^2); J = 2*n*k*(1-n^2)*(1-k^2);
        L = 4*n^2*k^2; g = [L 2*J (2*L-3*I) J];
        t0 = roots(g); t0 = t0(min(abs(t0))==abs(t0)); phi = acos(t0);
        % construct affine transfomation and transform circle
        v = (c.' + k^2*a.')/(k^2+1);
        U = [(a.'-v) ((v-b.')/(k*n)+(v-a.')*t0)/sin(phi)]; U_ = inv(U);
        center = [(1-k^2)/2;(k*(1-n^2)-t0*(1-k^2)*n)/(2*n*sin(phi))];
        R2 = (I-J*t0-L*t0^2)/(4*n^2*sin(phi)^2);
        center = U*center+v; M = U_.'*U_; E = {M/R2,center.',1};
    case 5 % convert the linear form to center form
        M = [p(1) p(3); p(3) p(2)];
        c = -inv(M)*[p(4);p(5)];
        z = c.'*M*c-p(6);
        E = {M/z,c.',1};
end
if nP>2 % convert center form to parameterization
    [~,D,V] = svd(E{1});
    ab = 1./sqrt(diag(D)); phi = sign(V(end))*acos(V(2));
    E = [E{2} sum(ab)/2 diff(ab)/2 phi];
end
end

function E = me_(points)
nP = size(points,1); 
switch nP
    case {1 2 4}
        E = points;
    case 3
        p1 = points(1,:); p2 = points(2,:); p3 = points(3,:);
        c = 1/3*(p1+p2+p3);
        M = 2/3*((p1-c).'*(p1-c)+(p2-c).'*(p2-c)+(p3-c).'*(p3-c));
        E = [inv(M);c];
    otherwise
        E = [];
end
end
