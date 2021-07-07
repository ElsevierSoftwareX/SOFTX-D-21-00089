function out = integrate(RP, varargin)
%INTEGRATE Integrate along the contours using the trapezoidal rule
%   You will usually not need to use this method directly but rather use
%   expand of viewFields
%
%   see also expand, viewFields

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

persistent parser;
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'RieszProjection/integrate';
    addOptional(parser,'quantity','DipoleEmission',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addOptional(parser,'w0',0,...
        @(x)validateattributes(x,{'numeric'},{'real'}));
    addParameter(parser,'precision',0,...
        @(x)validateattributes(x,{'numeric'},{'nonnegative'}));  
    addParameter(parser,'upperbound',128,...
        @(x)validateattributes(x,{'numeric'},{'nonnegative','scalar'})); 
    addParameter(parser,'index',0,...
        @(x)validateattributes(x,{'numeric'},{'nonnegative','scalar'})); 
    addParameter(parser,'outNames',{},...
        @(x)validateattributes(x,{'cell'},{}));
    addParameter(parser,'recalculateEvs',false,...
        @(x)validateattributes(x,{'logical'},{'scalar'}));
        
end
parse(parser,varargin{:});
vars = parser.Results;
if ~vars.precision, vars.precision = RP.precision; end % for adaptivity
if ~vars.w0, vars.w0 = RP.expansionPoints; end % default
if isempty(vars.w0), out = {}; return; % check if validity of range
elseif min(vars.w0)<min(RP.expansionPoints) || ...
        max(vars.w0)>max(RP.expansionPoints)
    error(['You must select expansion points within the range '...
        'of the property ''expansionPoints'' [%e,%e]'],...
        min(RP.expansionPoints),max(RP.expansionPoints));
end

% if necessary update scattering solutions and derived quantities
if ~RP.up2date(1), RP.evaluate; end
vars.quadgk = RP.shape{3}>1; % number of gk nodes

if strcmpi(vars.quantity, 'Field')
    if RP.useEigenvalues
        values = RP.fields(end);
    else
        values = [RP.fields(1:RP.nModes) RP.fields(end)];
    end
else
    if isempty(RP.derivedQuantities.(vars.quantity))
        RP.updateDerivedQuantities;
    end
    values = RP.derivedQuantities.(vars.quantity);
end

lQs = RP.linearQ; % quantities linear in the field
q2 = ~any(strcmp(vars.quantity,lQs)); % quantity quadratic in the field

% if the eigenvalues are not known get them
update = ~(RP.up2date(3)||RP.up2date(4)) && RP.useEigenvalues;
if isempty(RP.poles) || update || ...
        (vars.recalculateEvs&&~isempty(RP.selectedPoles))
    RP.precision = Inf; RP.up2date(2) = false;
    if ~(RP.up2date(3)||isempty(RP.poles))
        vars.recalculateEvs = true; 
    end
    if strcmp(vars.quantity,lQs{1})
        out = RP.computeEigenvalues(values{end},vars); return
    else
        RP.derivedQuantities.(RP.linearQ{1}) = {};
        RP.qKeys.(RP.linearQ{1}) = struct;
        RP.updateDerivedQuantities;
        valuesL = RP.derivedQuantities.(RP.linearQ{1});
        varsL = vars; varsL.quantity = RP.linearQ{1};
        RP.computeEigenvalues(valuesL{end},varsL);
    end
end

nContours = 1 + RP.nModes + RP.nModes*q2;
if vars.index % the contours to be integrated
    inds = vars.index; % user defined (must be scalar)
else
    inds = 1:nContours;
end
% collect convergence data if the values are floats
conv = isfloat(values{1});
out = cell(1,min(RP.nModes+1,length(inds)));
if conv, convData = cell(1,RP.nModes+1); end

for it1 = inds
    vars.bg = it1==nContours;
    vars.g = RP.g;
    if vars.bg
        vs = values{end};
    elseif isfield(RP.residues,vars.quantity) && RP.useEigenvalues
        vars.g = @(w,w0) 1./(RP.poles(RP.selectedPoles{it1})-w0);
        vars.residues = RP.residues.(vars.quantity); vs = [];
    elseif RP.useEigenvalues && RP.selectedPoles{1}(1)
        % use the eigenvalues instead of contours
        idx_wn = zeros(size(RP.selectedPoles),'logical');
        idx_wn(mod(it1-1,RP.nModes)+1) = true;
        wn = RP.poles(RP.selectedPoles{idx_wn});
        wns = RP.poles([RP.selectedPoles{~idx_wn}]);
        if it1>RP.nModes, wn = conj(wn); end
        if q2, wns = [wns conj(wns) conj(wn)]; end %#ok<AGROW>
        if ~isempty(wns)
            a = @(w) prod((w-wns)./(wn-wns),2); 
        else
            a = @(w) 1;
        end
        vars.g = @(w,w0) a(w).*RP.g(wn,w0);
        vs = values{end};
    else
        vs = values{it1};
    end
    [res,convD] = RP.quad(vs,vars,conv,it1); % perform integration
    if vars.index
        out{1} = res; 
    elseif it1<=RP.nModes
        out{it1} = res; if conv, convData{it1} = convD; end
    elseif ~vars.bg % add the contributions from mirrored contours
        it = mod(it1-1,RP.nModes)+1; out{it} = out{it} + res;
        if conv&&~isempty(convD), convData{it} = convData{it}+convD; end
    else % the contribution of the background contour
        out{end} = res; if conv, convData{end} = convD; end 
    end
end
if conv
    RP.convergenceData.(vars.quantity) = convData;
end
if ~conv || ismember(Inf,vars.precision), return; end

% adaptive refinement
% check if the target precision is reached 
tol = vars.precision; % tol(1) absolute and tol(end) relative tolerance
RP.tmp_ = out; % RP.expansion is not yet set
finished = RP.error(vars.quantity,[0 tol],'m');
if isfield(RP.reference,vars.quantity)
    if finished(end), return; end % comparison with physical quanity
    finished = finished(1:end-1); gl = false;
end
if isempty(finished) || (all(finished) && length(finished)==RP.nModes)
    % add Kronrod nodes to background, which was integrated with Gauss
    finished = [finished false]; gl = true;
elseif length(finished)==RP.nModes && ~RP.useEigenvalues
    finished = [finished true]; % refine modal contours first
elseif all(finished), return;
end
if ~RP.useEigenvalues % modal contributions
    nPnts = [RP.nPoints_ RP.nPointsB];
    finished(~finished) = nPnts(~finished)>=vars.upperbound;
    if all(finished), return; end
    RP.nPoints_(~finished(1:end-1)) = RP.nPoints_(~finished(1:end-1))*2;
elseif numel(RP.contours{1})>vars.upperbound, return;
end
if length(finished)>RP.nModes || RP.useEigenvalues
    if ~finished(end) && vars.quadgk
        if ~gl
            err_abs = abs(real(RP.convergenceData.(vars.quantity){end}));
            err_rel = err_abs./abs(real(out{end}));
            errBnds = min(cat(3,err_abs/tol(1),err_rel/tol(end)),[],3);
            errBnds = squeeze(max(errBnds,[],1));
            % the weightet errorbounds are below tolerance if < 1
            v = RP.shape{1}; cr = []; sp = RP.shape{2};
            if ismember(sp,'ec'), cr = v(1:2); v = v(3:end); end
            h = abs(diff(v)); pathlen = sum(h);
            ndx = real(errBnds) > h/pathlen; 
            if ~any(ndx), return; end
            k = 1; c = RP.includeConjugatedPoles;
            if ismember(sp,'ec') && c, ndx = ndx|ndx(end:-1:1); end
            c = c && sp=='p';
            v_new = zeros(1,length(v)+sum(ndx)+c*sum(ndx)-1);
            for it = 1:length(ndx)
                v_new(k) = v(it); k = k+1;
                if ndx(it)
                    v_new(k) = (v(it+1)+v(it))/2; k = k+1;
                    if c, v_new(k) = conj(v_new(k-1)); k = k+1; end
                end
            end
            if c % sort the vertices to from a polygon
                phi = mod(2*pi+eps(10)+angle(v_new-mean(v_new)),2*pi);
                [~,idx] = sort(phi);
                d = abs(diff(v_new(idx)));
                idx = idx([d>eps(max(d))*10 true]);
                v_new = v_new(idx);
            end
            if c, endPoint = v_new(1); else, endPoint = v(end); end
            RP.shape{1} = [cr v_new endPoint];
        end
        RP.getContours('gk');
        RP.nPointsB_ = numel(RP.contours{end});
        if ~RP.up2date(4), RP.up2date(3) = false; end
    elseif ~finished(end)
        RP.nPointsB=RP.nPointsB*2;
    end
end
out = RP.integrate(vars.quantity, vars.w0, ...
    'precision', vars.precision, 'upperbound', vars.upperbound);
end
