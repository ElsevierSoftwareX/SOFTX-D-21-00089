function varargout = export(sc,tags,pointListOrKeys,quantity)
%EXPORT Export field values at selected points or on a cartesian grid
%   fieldvalues = EXPORT(fb, tags, pointList, quantity, resultbag) returns
%   the field values at the points in pointList whose size must be Nx3.
%   The quantity can be chosen from those listed in the <a href="matlab:
%   web(['https://docs.jcmwave.com/JCMsuite/html/ParameterReference' ...
%   '/a062a9a76a59a014fec31f494b7447a1.html'],'-browser')" 
%   >parameter reference</a>
%   of JCMsuite. The corresponding resultbag is sc.pprbs.pointList
%
%   cartesianFields = EXPORT(fb, tags, keys, quantity) returns the 
%   fieldvalues on a cartesian grid defined by the fields of the struct
%   keys. The fields can be (N)GridPoints* where the asterix represents
%   X, Y or Z for further details please refer to the <a href="matlab:
%   web(['https://docs.jcmwave.com/JCMsuite/html/ParameterReference/'...
%   '0354993cda1716b2afb7a6af10ec11c7.html'],'-browser')"
%   >parameter reference</a>
%   of JCMsuite. The returned value is a cell array containing structs
%   with the fields: 'field', 'header', 'X', 'Y', 'Z' and 'domainIds'. The
%   correspnding resultbag is sc.pprbs.cartesian
%
%   [ids,tags] = EXPORT(fb, tags, keys, quantity) instead of the
%   results the ids are returned which are assigned to the jobs by the 
%   daemon at a call of jcmwave_solve. To get the results you must call
%   jcmwave_daemon_wait(ids,resultbag). The result k is than obtained
%   calling cartesianField = resultbag.results_(tags(k)).result.
%
%   see also RieszProjection, jcmwave_resultbag, jcmwave_solve

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

if isempty(tags), varargout(1) = {[]}; return; end
if nargin<4, quantity = sc.quantity; end
if isnumeric(pointListOrKeys) && size(pointListOrKeys,2)==3
    out = exportPoints(sc,tags,pointListOrKeys,quantity);
elseif isstruct(pointListOrKeys)
    out = exportCartesian(sc,tags,pointListOrKeys,quantity,nargout);
else
    error('You must pass a list of points (Nx3) or a struct.');
end
if nargout > 0, varargout(1:nargout) = out(1:nargout); end
end

% export points from point list
function out = exportPoints(sc,tags,pointList,q)

% create jcmpt file
exportList = sc.exportPoints;
jcmpt = [sc.postproDir filesep 'points.jcmpt'];
id = fopen(jcmpt, 'w');
fprintf(id, exportList);
fclose(id);
jcmpt = jcmpt(1:end-1);

% get resultbag with scattering results
rbsc = sc.resultbags{1+isfield(sc.resultbags{2}.results_,tags{1})};

% prepare the keys
pp = zeros(1,length(tags),'int32');
keys.px = pointList(:,1);
keys.py = pointList(:,2);
keys.pz = pointList(:,3);
keys.quantity = q;
keys.addSingularFields = 'yes';
keys.oPath = 'out.jcm';

% get resultbags and export
if ~isfield(sc.pprbs,'pointList')
    filename = [sc.resultbagDir filesep 'pointList' '.mat'];
    relevantFields = sort([fieldnames(keys);'fPath']);
    sc.pprbs.('pointList') = jcmwave_resultbag(filename,relevantFields);
end
rb = sc.pprbs.pointList;
for it1 = 1:length(tags)
    keys.fPath = rbsc.results_.(tags{it1}).result{1}.file;
    tags{it1} = rb.get_tag(keys);
    pp(it1) = jcmwave_solve(jcmpt,keys,rb,'temporary','yes');
end

% collect results
jcmwave_daemon_wait(pp,rb);
points = cell(size(tags)); 
for it1 = 1:length(tags)
    try
        points{it1} = rb.results_.(tags{it1}).result{1}.q;
        points{it1} = points{it1}{1}.ElectricFieldStrength{1};
    catch
        log = resultbag.get_log(keys);
        error(log.Log.Error)
    end
end
out = cat(3,points{:});
end

% export on a cartesian grid
function out = exportCartesian(sc,tags,keys,q,n)

fieldNames = fields(keys);
nameX = 'NGridPointsX'; nameY = 'NGridPointsY'; nameZ = 'GridPointsZ';

if any(contains(fieldNames, 'GridPointsX'))
    nameX = fieldNames{contains(fieldNames, 'GridPointsX')};
else, keys.(nameX) = 100;
end
if any(contains(fieldNames, 'GridPointsY'))
    nameY = fieldNames{contains(fieldNames, 'GridPointsY')};
else, keys.(nameY) = 100;
end
if any(contains(fieldNames, 'GridPointsZ'))
    nameZ = fieldNames{contains(fieldNames, 'GridPointsZ')};
else, keys.(nameZ) = 0;
end
Z = sprintf('%1$s=%%(%1$s)e', nameZ);
Y = sprintf('%1$s=%%(%1$s)e', nameY);
X = sprintf('%1$s=%%(%1$s)e', nameX);

if ~isfield(sc.pprbs,'cartesian')
    filename = [sc.resultbagDir filesep 'cartesian.mat'];
    relevantFields = sort([fieldnames(keys);'fPath';'oPath';'quantity']);
    sc.pprbs.('cartesian') = jcmwave_resultbag(filename,relevantFields);
end
rb = sc.pprbs.('cartesian');

pos = strfind(rb.filepath_, filesep); 
pDir = rb.filepath_(1:pos(end-1)-1);
outDir = [pDir filesep 'superpositions'];
if ~exist(outDir, 'dir'), mkdir(outDir); end

% create jcmpt file
jcmpt = [outDir filesep 'cartesian.jcmpt'];
if ~exist(jcmpt,'file')
    exportCartesian = sc.exportCartesian;
    id = fopen(jcmpt, 'w');
    fprintf(id, exportCartesian, X, Y, Z);
    fclose(id);
end
jcmpt = jcmpt(1:end-1);

% get resultbag with scattering results
if isfield(sc.resultbags{2}.results_,tags{1})
    rbsc = sc.resultbags{2};
    files = cellfun(@(x){rbsc.results_.(x).result{1}.file},tags);
else
    rbsc = sc.pprbs.superposition;
    files = cellfun(@(x){rbsc.results_.(x).keys.oPath},tags);
end

pp = zeros(1,length(tags),'int32');
keys.oPath = 'out.jcm';
keys.quantity = q;
for it = 1:length(tags)
    keys.fPath = files{it};
    pp(it) = jcmwave_solve(jcmpt,keys,rb,'temporary','yes');
    tags{it} = rb.get_tag(keys);
end

if n==2, out = {pp tags}; return; end
% collect results
jcmwave_daemon_wait(pp,rb);
out = cell(size(tags));
for it1 = 1:length(tags)
    out(it1) = rb.results_.(tags{it1}).result(1);
end
out = {out};
end
