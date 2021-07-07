function [res, convData] = quad(RP, values, vars, conv, it)
%QUAD apply quadrature rule. Either the trapezoidal rule or Gauss-Kronrod.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

% Dimensions of the argument 'values':
% axis 1: w0, axis 2: data, axis 3: nodes, axis 4: subintervals
w0 = vars.w0;
if nargin<5, it = RP.nModes+1; end
if RP.useEigenvalues || vars.bg
    w = RP.contours_{end}; weights = RP.weights{end};
else
    w = RP.contours_{it}; weights = RP.weights{it};
end
vars.quadgk = size(w,2)>1;
w = reshape(w,[1 1 size(w)]);
if iscell(weights)
    weights = cellfun(@(x){reshape(x,size(w))},weights);
else
    weights = reshape(weights,size(w));
end
nSolutions = length(w0);
% get the weights
if nargin(vars.g)==2, gw = vars.g(w,w0); else, gw = vars.g(w); end

% if the residues are known use them
if isfield(vars,'residues') && ~vars.bg && RP.useEigenvalues
    res = squeeze(gw.*vars.residues(it)); convData = {}; return;
end

if iscell(weights)
    weights = cellfun(@(x){gw.*reshape(x,size(w))},weights);
else
    weights = reshape(weights,size(w)).*gw;
end

% integrate
if iscell(values) && isa(RP.f,'Scattering')
    % here a special case for the JCMsuite specific method vewFields is
    % handled. In the future this should be replaced by a more general
    % method that allows to perform the summation with an external program
    res = cell(1,nSolutions);
    if iscell(weights), weights = weights{end}; end
    for it1 = 1:nSolutions
        args = {weights(it1,:) w0(it1) ''};
        if ~isempty(vars.outNames), args{3} = vars.outNames{it1}; end
        res{it1} = squeeze(RP.f.sum(values,args{:}));
        convData = {};
    end
    res = [res{:}];
elseif vars.quadgk
    if iscell(weights)
        errw = weights{2};  weights = weights{1};
    else
        errw = [];
    end
    res = sum(weights.*values,[3 4]);
    if ~isempty(errw) && conv
        convData = squeeze(sum(errw.*values,3));
    else
        convData = {};
    end
elseif conv
    nv = size(values,3);
    while mod(nv,2)==0 && nv>2, nv = nv/2; end
    k = log2((size(values,3))/nv);
    skip = 2^k; T = zeros(nSolutions,size(values,2),k);
    for it2 = 1:k+1
        vs = values(:,:,1:skip:end);
        weights_ = weights(:,:,1:skip:end)*skip;
        T(:,:,it2) = sum(weights_.*vs,3);
        skip = skip/2;
    end
    res = T(:,:,end);
    convData = T;
elseif isfloat(values)
    res = sum(weights.*values,3);
    convData = {};
end

% Add field with the number of integration points to convergence data.
% This has to be done here, as this function is called by integrate and
% computeEigenvalues.
if conv && ~isfield(RP.convergenceData,'nPoints')
    nPoints = cell(1,RP.nModes+1);
    for it = [1:RP.nModes length(RP.contours_)]
        if RP.shape{3}==1
            nv = size(RP.contours_{it},1);
            while mod(nv,2)==0 && nv>2, nv = nv/2; end
            k = log2((size(RP.contours_{it},1))/nv);
            nPoints{min(end,it)} = nv*2.^(0:k);
        elseif iscell(RP.weights)
            m = size(RP.contours_{it},1)*size(RP.contours_{it},2);
            n = [m/2-size(RP.contours_{it},2)/2 m];
            nPoints{min(end,it)} = n;
        end
        if RP.useEigenvalues
            nPoints(:) = repmat(nPoints(1),1,RP.nModes+1); break;
        end
    end
    RP.convergenceData.nPoints = nPoints;
end
end
