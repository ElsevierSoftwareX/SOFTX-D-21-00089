function evaluate(RP)
%EVALUATE Evaluate scattering problems and update derived quantities
%   EVALUATE(RP) The scattering poblems along the contours are solved. If
%   the property 'derivedQuantities' has fields, these are updated.
%
%   see also RieszProjection/updateDerivedQuantities

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

if isempty(RP.contours)
    error('No contours: Please call the method ''getContours''.')
end
d = length(RP.fields)-length(RP.contours);
if isempty(RP.undefinedValues)
    if ~isempty(RP.fields), return; end
    RP.fields = RP.f(RP.contours);
else
    locs = find(cellfun(@any, RP.undefinedValues));
    frequencies = cell(1,length(locs));
    for it = 1:length(locs)
        l = locs(it);
        frequencies{it} = RP.contours{l}(RP.undefinedValues{l});
    end
    newV = RP.f(frequencies);
    fV = [RP.fields{:}];
    if length(RP.fields)>length(RP.selectedPoles)
        RP.fields(1:d) = []; 
    end
    for it = 1:length(RP.contours)
        locV = RP.definedValues{it};
        ndef = RP.undefinedValues{it};
        getEmpty = [class(RP.fields{1}) '.empty(0,length(ndef))'];
        RP.fields{it} = eval(getEmpty);
        locnV = locs==it;
        if any(ndef)
            RP.fields{it}(ndef) = newV{locnV};
        end
        if ~isempty(locV)
            RP.fields{it}(~ndef) = fV(locV);
        end
    end
end

% update and reset quantities
RP.nModes = 0;
if ~isempty(RP.selectedPoles) && RP.selectedPoles{1}(1)
    RP.nModes = length(RP.selectedPoles);
end
RP.updateDerivedQuantities(d);
RP.contours_f = RP.contours;
RP.undefinedValues = {};
RP.definedValues = {};
RP.up2date(1:2) = [true false];
RP.convergenceData = struct;
% The property 'fieldExpansion' is reset to an empty cell as for given 
% expansionPoints entries are not reevaluated if present.
RP.fieldExpansion = {};
end

