function updateDerivedQuantities(RP, d)
%UPDATEDERIVEDQUANTITIES Update all quantities
%    UPDATEDERIVEDQUANTITIES(RP) All quantities saved as fields of the
%        property 'derivedQuantities' are updated.
%
%    UPDATEDERIVEDQUANTITIES(RP,d) If the number of contours has changed
%        the cell array is shrunk such that its length matches the number 
%        of modes.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

% The derived quantities may either be scalars or vectors of length m. They
% are avaluated at p*q integration points, where, in case of a higer order
% quadrature, q>1 is the number of subintervals. The dimension of the 
% derived quantity of a given contour must be (1,m,p,q). The first
% dimension is reserved for the number of expansion points.

if nargin==1, d=0; end % d elements will be removed to shrink the size
derivedquantities = fieldnames(RP.derivedQuantities); gk = RP.shape{3}>1;
nContours = length(RP.fields); nMC = min(nContours,RP.nModes+1);

for it1 = 1:length(derivedquantities)
    q = derivedquantities{it1};
    isLinear = ismember(q,RP.linearQ);
    getAll = isempty(RP.derivedQuantities.(q));
    if ~getAll
        if isempty(RP.undefinedValues), continue; end
        quantity = RP.derivedQuantities.(q);
        quantity = cellfun(@(x){reshape(x,size(x,2),[])},quantity);
        if isLinear, shift = sum(cellfun(@length,quantity(1:end-1))); end
        quantity = cat(2,quantity{:});
        if d, RP.derivedQuantities.(q)(1:d/(1+isLinear)) = []; end
        undefined = RP.undefinedValues([1:nMC-1 end]);
    end
    
    % prepare the keys and iterate over the contours
    qkeys = RP.qKeys.(q);
    fs = cell(nContours,2-isLinear);
    for it2 = 1:nMC
        if it2==nMC, it = nContours; else, it = it2; end
        vsp = RP.fields{it};
        if getAll
            ndef = [];
        elseif any(undefined{it2})
            ndef = undefined{it2};
            vsp = vsp(ndef);
        else
            continue
        end
        if ~isLinear
            if it2<nMC, vsn = reverseV(vsp,it,ndef); end
            vsp_ = RP.fields{mod(it2-1,nMC)+nMC};
            if ~isempty(ndef), vsp_ = vsp_(ndef); end
            vsn_ = reverseV(vsp_,it,ndef);
        end
        % linear case
        if isLinear, fs{it2,1} = vsp; continue; end
        % quadratic case
        fs(it,:) = {vsp,vsn_};
        if it2<nMC, fs(it2+nMC-1,:) = {vsp_,vsn}; end
    end
    
    notEmpty = ~cellfun('isempty',fs(:,1)); v = cell(size(fs,1));
    v(notEmpty) = RP.f(fs(notEmpty,:),q,qkeys);
    
    % distribute the results
    for it2 = 1:nContours
        if it2==nMC&&isLinear, it = nContours; else, it = it2; end
        if isLinear && it2>nMC, break; end
        sz = size(RP.contours_{it});
        if ~getAll
            ndef = RP.undefinedValues{it};
            locV = RP.definedValues{it};
            if (it2==nMC&&isLinear) && ~isempty(locV) % linear case
                locV = locV-shift;
            end
            RP.derivedQuantities.(q){it2} = ...
                zeros([1,size(quantity,1),sz]);
            RP.derivedQuantities.(q){it2}(1,:,~ndef) = ...
                quantity(:,locV);
            if isempty(v{it2}), continue; end
        end
        if getAll 
            % axis 1: w0, axis 2: data, axis 3: nodes, axis 4: subs
            sz_ = [1 size(v{it2},1) sz];
            RP.derivedQuantities.(q){it2} = reshape(v{it2},sz_);
        else
            RP.derivedQuantities.(q){it2}(1,:,ndef) = v{it2};
        end
    end
end

    function vr = reverseV(values,ndx,ndef)
        c1 = RP.contours{ndx}(:);
        M = max(abs(c1));
        if ~isempty(ndef), c1 = c1(ndef); end
        c1 = c1([1 end]);
        vr = values(end:-1:1);
        if gk && ndx==length(RP.contours)
            ngk = RP.shape{3};
            if imag(c1(2))>0 && isempty(ndef)
                vr = circshift(vr,-ngk);
            elseif imag(c1(2))>0
                n = sum(mod(find(~ndef)+ngk/2,length(ndef))<ngk);
                vr = circshift(vr,-ngk+n);
            end
        else
            if ndx<nMC
                c2 = RP.contours{ndx+nMC-1}; 
                if ~isempty(ndef), c2 = c2(ndef); end
                c2 = c2(1);
            else
                c2 = c1(1);
            end
            if imag(c2)+imag(c1(1))<2*eps(M)
                vr = circshift(vr,1);
            end
        end
    end
end

