function getReferenceSolution(RP, quantity, w0, keys)
%GETREFERENCESOLUTION Compute the reference solutions on the real line
%   GETREFERENCESOLUTION(RP, quantity) Computes the reference solution of
%   the given quantity.
%   The results are saved in the property 'reference', which is a struct
%   with field names corresponding to the quantities.
%
%   GETREFERENCESOLUTION(RP, quantity, w0) The reference solutions are
%   computed at the given real frequencies.
%
%   GETREFERENCESOLUTION(RP, quantity, keys) The struct keys is passed to
%   the interface may contain fields for parameter substitution. Be aware
%   that the keys are not passed to the call of the interface used to solve
%   the linear systems but only for the derivation of the target quantity.
%
%   GETREFERENCESOLUTION(RP, quantity, w0, keys) As above.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

% check the input arguments and set defaults if necessary
keys_ = struct;
if nargin < 3 || isstruct(w0)
    if nargin==3 && isstruct(w0), keys_ = w0; end
    if isfield(RP.w0,quantity)
        w0 = RP.w0.(quantity);
    elseif ~isempty(RP.referencePoints)
        w0 = RP.referencePoints;
    elseif ~isempty(RP.expansionPoints)
        w0 = RP.expansionPoints.'; 
    else
        warning('No reference points available.'); return
    end
elseif nargin==3
    if ~isreal(w0), error('The reference frequencies must be real'); end
elseif nargin == 4
    keys_ = keys;
end  

% if necessary evaluate the scattering problems at the expansion points
% query only those which are not already present
if ~isfield(RP.reference, 'fields')
    fs = RP.f({w0}); vs = fs{1};
    RP.reference.fields = vs;
    RP.reference.frequencies = w0;
elseif length(RP.reference.fields)~=length(w0) || ...
        any(RP.reference.frequencies~=w0)
    w0_ = uniquetol([w0;RP.reference.frequencies]);
    [ix,loc] = ismembertol(w0_,RP.reference.frequencies);
    if any(~ix)
        fs = eval([class(RP.reference.fields) '.empty']);
        fs(ix) = RP.reference.fields(loc(ix));
        fs_ = RP.f({w0_(~ix)}); fs(~ix) = fs_{1};
        RP.reference.fields = fs;
    end  
    [~,loc] = ismembertol(w0,w0_);
    vs = RP.reference.fields(loc);
    RP.reference.frequencies = sort(w0_);
else
    [~,loc] = ismembertol(w0,RP.reference.frequencies.');
    vs = RP.reference.fields(loc);
end

% check if the existing results are up to date and which can be reused
if isempty(fieldnames(keys_))
    if isfield(RP.rKeys,quantity), keys_ = RP.rKeys.(quantity);
    elseif isfield(RP.qKeys,quantity), keys_ = RP.qKeys.(quantity);
    end
end
if isfield(RP.w0,quantity) && isequal(keys_,RP.rKeys.(quantity))
    [ix,loc] = ismembertol(w0,RP.w0.(quantity));
    if all(ix)
        if length(RP.reference.(quantity))>length(w0)
            ref = RP.reference.(quantity)(loc(ix),:);
            RP.reference.(quantity) = ref;
            RP.w0.(quantity) = w0;
        end
        return;
    end
    vs = vs(~ix);
    v_ = RP.reference.(quantity)(loc(ix),:);
else
    ix = zeros(size(w0),'logical');
end

% save keys
RP.rKeys.(quantity) = keys_;
q = ['reference' quantity];

% get the referece solutions
v = RP.f({vs},q,keys_);

ref = zeros(length(w0),size(v{1},1));
if any(ix), ref(ix,:) = v_; end 
ref(~ix,:) = v{1}.';

% save the results

% If two NA are defined for the post process 'Radiation.jcmpt', the second
% one is used to get the photon collection efficiency which is the
% normalized radiated power. If the dipole emission is available the
% quantity dipole power collection efficiency is evaluated.
if strcmpi(quantity,'Radiation')
    if isfield(RP.reference,'DipoleEmission')
        w0ref = RP.w0.DipoleEmission;
        [ndxref,loc] = ismembertol(w0,w0ref);
        dpe = RP.reference.DipoleEmission(loc(ndxref));
        RP.reference.DipolePowerCollectionEfficiency = ref(ndxref,1)./dpe;
        RP.w0.DipolePowerCollectionEfficiency = w0(ndxref);
    end
% The normalized decay rate is not stored sperately but the dipole bulk
% emission returned by the post process is passed to the class 
% 'Scattering' from which it can be obtained calling the method 
% 'dipoleBulkEmission'.
elseif strcmpi(quantity,'DipoleEmission')
    if isfield(RP.reference,'Radiation')
        w0ref = RP.w0.Radiation;
        [ndxref,loc] = ismembertol(w0,w0ref);
        rad = RP.reference.Radiation(loc(ndxref),1);
        RP.reference.DipolePowerCollectionEfficiency = rad./ref(ndxref);
        RP.w0.DipolePowerCollectionEfficiency = w0(ndxref);
    end
end
RP.reference.(q(10:end)) = ref; RP.w0.(quantity) = w0;
end
