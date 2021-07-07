function varargout = getDir(RP,index,omega)
%GETDIR get directory of scattering problem
%   dir = GETDIR(RP,index,omega) The directory of the scattering problem is
%   returned whose frequency is the closest to omega. The index selects
%   the contour. If there are n contributions the index n+1 refers to the
%   background contour and the index n+2 to the reference solutions.
%   [dir,frequency] = GETDIR(RP,index,omega) Additionally the frequency of
%   the corresponding scattering solution is returned.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021
if index<RP.nModes+1
    if length(RP.fields)>RP.nModes+1 && imag(omega)>0
        index = index + RP.nModes;
    end
    tags = RP.fields{index};
    frequencies = RP.contours_f{index};
elseif index==RP.nModes+1
    tags = RP.fields{end};
    frequencies = RP.contours_f{end};
elseif index==RP.nModes+2 
    if isfield(RP.reference,'fields')
        tags = RP.reference.fields;
        frequencies = RP.reference.frequencies;
    else
        error('No reference solution available.')
    end
else
    error(['The index %d is invalid. It must be larger 1 and smaller '...
        'or equal %d.'],index,length(RP.contours_f)+1);
end
rb = RP.f.resultbags{1+(index==RP.nModes+2)};
[~,index] = min(abs(frequencies));
directory = rb.results_.(tags{index}).result{1}.file;
a = strfind(directory,filesep);
directory = directory(1:a(end-1)-1);
varargout{1} = directory;
if nargout == 2
    varargout{2} = tags.frequencies(index);
end
end

