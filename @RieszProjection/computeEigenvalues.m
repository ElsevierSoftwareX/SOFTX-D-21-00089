function res = computeEigenvalues(RP,values,vars)
%COMPUTEEIGENVALUES solve the nonlinear eigenvalue problem using the
% Rieszprojection based eigensolver which has been proposed in 2020 by
% Binkowski et al.:
% https://www.sciencedirect.com/science/article/pii/S0021999120304526.
% This private method is called by integrate.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

rEvs = vars.recalculateEvs; m = [RP.selectedPoles{:}];
if ~RP.useEigenvalues, m = []; end
center = sum(RP.expansionPoints([1 end]))/2;
radius = diff(RP.expansionPoints([1 end]))/2+RP.minD(1);
nEvs = max(RP.nEigenvaluesExpected,length(m)); nEqs = 2*nEvs;
mu = zeros(nEqs,1); vars.bg = true;
for it1 = 1:nEqs
    vars.g = @(w) 1/(2i*pi) * ((w-center)/radius).^(it1-1);
    mu(it1) = RP.quad(values,vars,false);
%    if ~RP.useEigenvalues
%        vars.bg = false;
%        for it2 = 1:RP.nModes
%            values1 = RP.derivedQuantities.(vars.quantity){it2};
%            mu(it1) = mu(it1)-RP.quad(values1,vars,false,it2);
%        end
%        vars.bg = true;
%    end
end
% guess for the starting point
normx = radius; normy = mu(1);
radius = radius/normx; center = center/normx; mu = mu/normy;
x0 = rand(1,nEqs) - rand(1,nEqs) - 1i*rand(1,nEqs);
x0 = center + radius*x0;
if rEvs
    x0(1:length(m)) = RP.poles(m)/normx;
end

options = optimoptions('fsolve',...
    'Display', 'off', ...
    'MaxFunctionEvaluations', 10000000, ...
    'MaxIterations', 10000000, ... ...
    'OptimalityTolerance', 1e-16, ...
    'StepTolerance', 1e-16, ...
    'FunctionTolerance', 1e-16, ...
    'ScaleProblem', 'jacobian', ...
    'SpecifyObjectiveGradient',true);

% solving the NLSE
x1 = fsolve(@F,x0,options);

% results
ps = x1(1:nEvs)*normx; rs = x1(nEvs+1:end)*normy;
sel = abs(rs/max(rs))>min(1e-3,RP.precision_(end));
ps = ps(sel); rs = rs(sel); 
[ps,idx] = sort(ps,'ComparisonMethod','real'); rs = rs(idx);
inC = RP.inside(ps); if rEvs, RP.poles_(m) = []; ps = ps(inC); end
RP.poles_ = sort([ps RP.poles_],'ComparisonMethod','real');
rs = rs(inC); ps = ps(inC); n = length(rs);
RP.newBackground = false; convD = cell(1,n+1);
if ~rEvs, RP.useEigenvalues_ = true; end
RP.selectedPoles = num2cell(find(ismember(RP.poles_,ps)));
RP.nModes = n; m = length(RP.expansionPoints); vars.g = RP.g;
RP.convergenceData = rmfield(RP.convergenceData,'nPoints');
[bg,convD{end}] = RP.quad(values,vars,true);
RP.convergenceData.(vars.quantity) = convD; 
RP.up2date(3:4) = [true false];
res = mat2cell([rs./(ps-RP.expansionPoints) bg],m,ones(1,n+1));
RP.residues.(vars.quantity) = rs; 

% update plot
try 
    if ishandle(RP.plotSettings.fignumber.contours)
        RP.plot('contours')
    end
catch
end

    function [Y,J] = F(x)
        k = (-1:nEqs-1).'; fk = ((x(1:nEvs)-center)/radius).^k;
        Y = sum(x(nEvs+1:end).*fk(2:end,:),2)-mu;
        J = [k(2:end)/radius.*x(nEvs+1:end).*fk(1:end-1,:) fk(2:end,:)];
    end
end
