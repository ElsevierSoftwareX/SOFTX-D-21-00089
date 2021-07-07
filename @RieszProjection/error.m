function [err,n] = error(RP, quantity, index, kind)
%ERROR compute an estimate of the integration error
%   [err,x] = ERROR(RP, quantity, index, kind) The errors of the integrals
%   along contours corresponding to the specified indices are returned. The
%   kind of error can be 'absolute' ('abs'), 'relative' ('rel') or 
%   'minimum' ('min') which is the minimum of absolute and relative error
%   weighted with absolute and relative error tolerances defined in the
%   property 'precision'.
%   The returned values are cell arrays of shape (1,m) if the argument
%   'index' has m elements. If a Gauss-Kronrod quadrature rule is used for
%   integration, the error of the background contour is given by the
%   difference between Gauss and Gauss-Kronrod quadrature.
%   The indices can take values between 1 and n+2 where n is the number of
%   modes within the background contour which has the index n+1. The error
%   of the physical quantity is err{k} if index(k) = n+2. The output x{k}
%   contains in this case the number of points of the background contour.
%   The output err{k} is a vector and err{k}(m) is the error of the
%   integration with x{k}(m) abscissae.
%   The default kind of error is minimum and the default index is n+2.
%
%   err = ERROR(...) only the error bound(s) of the best available
%   approximation(s) is returned. A scalar or a vector (if length(index>1)
%   is returned.
%
%   The returned error is always the maximal error along the expansion
%   points or the reference points. 

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

% check which data is available
dat = isfield(RP.expansions,quantity);
if ~isempty(RP.tmp_)
    expansion = RP.tmp_; dat = 2; RP.tmp_ = [];
elseif dat 
    expansion = RP.expansions.(quantity);
end
refD = isfield(RP.reference,quantity);
convD = isfield(RP.convergenceData,quantity);
if convD, conv = RP.convergenceData.(quantity); end
useEvs = RP.useEigenvalues && length(RP.fields)==1;

% if the first element of index is 0 the other values are expected to
% define a tolerance and logicals are returned
if ~index(1) && length(index)<4
    tol = index(2:end); index = 0; convD = 2;
    if isscalar(tol), tol = [tol tol]; end
else
    tol = RP.precision_;
end

if ~tol(1), kind = 'relative'; elseif ~tol(2), kind = 'absolute'; end

% set defaults end check if there is enough data to compute an error
if nargin<4, kind = 'rel'; if nargin<3, index = RP.nModes+2; end; end
if any(isinf(RP.precision)) && kind(1)=='m', kind = 'relative'; end
if ~((dat&&all(index)&&refD)||convD)
    fprintf(['There is not enough data for an error estimation of ' ...
        'the expansion of the %s.\n'],quantity);
    err = {}; n = {}; return
elseif ismember(RP.nModes+2,index) && ~refD
    fprintf(['The implemented estimation of the error bound of the '...
        'physical quanity ''%s'' requires a reference solution. \n'...
        'Yet, estimations of error bounds of individual '...
        'contributions are available.\n'],quantity)
    index(index==RP.nModes+2) = []; index(index==RP.nModes+2) = [];
    if isempty(index), err = {}; n = {}; return; end
elseif convD==2 % check if error is below tol
    index = 1:RP.nModes+2;
    if useEvs, index(index<RP.nModes+1) = []; end
    if ~refD, index(index==RP.nModes+2) = []; end
    if isempty(conv{end}), index(index==RP.nModes+1) = []; end
    if isempty(index), err = true; return; end
end

pointEval = strcmpi(quantity, 'DipoleEmission');

compSum = ismember(RP.nModes+2,index) && refD; gk = RP.shape{3}>1; % quadgk
if compSum % prepare what is needed for a comparison with the reference
    inds = 1:RP.nModes+1; tot = 0;
    w0ref = RP.w0.(quantity);
    w0 = RP.expansionPoints;
    if strcmp(quantity,'RadiationPattern'), w0 = w0ref; end
    [wr,pos] = ismembertol(w0,w0ref);
    if isempty(pos(wr)), error('No matching reference solution'); end
    ref_ = real(RP.reference.(quantity)(pos(wr),:));
    if pointEval % bulk emission is not included
        try
            bulk = RP.dipoleBulkEmission(RP.expansionPoints(wr));
        catch
            warning('No bulk emission available'); bulk = 0;
        end
    end
else
    inds = index;
end

% compute the error 
err = cell(1,length(index)); it = 1;
for it1 = 1:length(inds)
    ndx = inds(it1);
    if convD && ~isempty(conv{ndx})
        data = real(conv{ndx});
        if compSum && ndx<RP.nModes+1
            if gk && useEvs
                tot = tot - expansion{ndx}(wr,:);
            elseif useEvs
                tot = tot - data(wr,:,:);
            else
                tot = tot - data(wr,:,end);
            end
        elseif compSum
            it_ = it+ismember(ndx,index); ref = ref_;
            if gk
                tot = real(expansion{end}(wr,:))+tot;
                if dat==2 && pointEval, tot = tot+bulk; end
                err{it_} = maxError(tot,ref);
                if convD==1
                    err2 = maxError(data(wr,:,:),ref,1); err1 = err{it_};
                    err{it_} = [err1+err2 err1];
                end
            else
                tot = data(wr,:,:)+tot;
                if pointEval, tot = tot+bulk; end
                err{it_} = maxError(tot,ref);
            end
        end
        if ismember(ndx,index)
            if ndx==RP.nModes+1 && gk
                ref = real(expansion{end});
                err{it} = maxError(data,ref,1); 
            elseif size(data,3)>1
                ref = real(data(:,:,end));
                err{it} = maxError(data(:,:,1:end-1),ref);
            else
                fprintf(['There is not enough data to compute ' ...
                    'the error of Mode %d (%s).\n'],ndx,quantity)
                if convD==2
                    if all([err{1:it-1}])
                        err{it} = false;
                    else
                        err{it} = true;
                    end
                else
                    err{it} = NaN; inds(it1) = 0;
                end
            end
        end
    else
        if ismember(ndx,index)
                fprintf(['There is not enough data to compute ' ...
                    'the error of Mode %d (%s).\n'],ndx,quantity)
            err{it} = NaN; inds(it1) = 0;
        end
        if compSum && ndx<RP.nModes+1
            tot = tot - real(expansion{ndx}(wr,:));
        elseif compSum
            ref = ref_; tot = tot+real(expansion{ndx}(wr,:));
            if dat==2 && pointEval, tot = tot+bulk; end
            err{it} = maxError(tot,ref);
        end
    end
    if ismember(ndx,index), it = it+1; end
end

% prepare output
if ~(any(inds)||compSum), err = {}; n = {}; return; end
if nargout>1
    index(index==RP.nModes+2) = index(index==RP.nModes+2)-1;
    n = RP.convergenceData.nPoints(index);
    if length(n)==1 && length(err{1})==1
        if isempty(n{1})
            n = {length(RP.contours{end}(:))};
        else
            n = {n{1}(end)};
        end
    else
        n = cellfun(@(y,z){y(1:length(z))},n,err);
    end
else
    err = cellfun(@(x)x(end),err);
end

    function err = maxError(data,ref,gk)
        if nargin<3, gk = false; end
        if gk
            err_abs = sum(abs(data),2);
        else
            err_abs = abs(data-ref);
        end
        err_rel = err_abs./abs(ref);
        err_rel = reshape(err_rel,[],size(err_rel,3));
        err_abs = reshape(err_abs,size(err_rel));
        if kind(1)=='m'
            [a,b] = min(cat(3,err_abs/tol(1),err_rel/tol(end)),[],3);
            [c,d] = max(a,[],1); d = d+(0:size(b,2)-1)*size(b,1);
            err = c.*tol(b(d));
            if convD==2, err = err(end)<tol(b(d(end))); end
        elseif kind(1)=='r'
            err = max(err_rel,[],1);
            if convD==2, err = err(end)<tol(end); end
        elseif kind(1)=='a'
            err = max(err_abs,[],1);
            if convD==2, err = err(end)<tol(1); end
        end
    end
end
