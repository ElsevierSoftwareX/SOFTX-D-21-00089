function out = f(sc, varargin)
%F Solves scattering problems at complex frequencies
%   F(sc, contours) This function is called when an instance sc of class
%   Scattering is called with (). The input argument contours is 
%   expected to be a cell array of size=(1,:) containing a vector of
%   complex frequencies. It solves the scattering problems on the given
%   contours. If the input is an empty cell 0 is returned. Otherwise a
%   cell array of the same shape as contours containing tags corresponding
%   to solutions of the scattering problems.
%
%   F(sc, tags, quantity, qkeys) Derive the specified quantity from the
%   solution of the linear systems represented by tags. The qkeys my
%   contain parameters for the corresponding post process.
%
%   See also SCATTERING, RIESZPROJECTION

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

if length(varargin)==1
    out = solveScatteringProblems(sc,varargin{1});
else
    [tags,q,qkeys] = varargin{:};
    out = cell(1,size(tags,1));
    for it = 1:size(tags,1)
        [out{it},tags{it}] = sc.submit(tags{it,:},q,qkeys);
    end
    for it = 1:size(tags,1)
        out{it} = sc.collect(out{it},tags{it},q,it==size(tags,1));
    end
end
end

function tags = solveScatteringProblems(sc,contours)
% check if a resource has been set up for computation
if isempty(jcmwave_daemon_resource_info)
    error('Please start a daemon');
end

% select or create resultbag
if isempty(contours), tags=0; return; end
isRef = length(contours)==1 && isreal(contours{1});
if isRef, scDir = sc.referenceDir; else, scDir = sc.scatteringDir; end
index = isRef+1;
if isempty(sc.resultbags{index})
    bag = sc.bags{index};
    filename = [sc.resultbagDir filesep bag];
    relevantfields = [fieldnames(sc.keys).' sc.fields];
    resultbag = jcmwave_resultbag(filename,relevantfields);
    sc.resultbags{index} = resultbag;
else
    resultbag = sc.resultbags{index};
end

% define the points which must be exported
keys = sc.keys;
if isfield(keys, 'position') && ~isRef
    pointList = keys.position;
    keys.addSingularFields = 'no';
else
    pointList = sc.pointList;
    keys.addSingularFields = 'yes';
end

nContours = length(contours);
nPoints = sum(cellfun(@length, contours));
jcmpSc = sc.projectFile;

% check if the postprocess has been removed properly
exportPointList = sc.exportPoints;
id = fopen(jcmpSc, 'r');
jcmpt = fread(id,'*char')';
fclose(id);
s = regexp(jcmpt,exportPointList(1:29),'once');
if ~isempty(s), jcmpt(s:end) = []; end
if isempty(pointList) && ~isempty(s)
    id = fopen(jcmpSc,'w');
    fwrite(id,jcmpt,'char');
    fclose(id);
% append a post process to the project file if needed
elseif ~isempty(pointList)
    if isempty(s)
        id = fopen(jcmpSc,'w');
        fwrite(id,[jcmpt sprintf(exportPointList)],'char');
        fclose(id);
    end
    resdir = [sc.projectFile(length(sc.projectDir)+2:end-6) '_results'];
    keys.fPath = [resdir filesep 'fieldbag.jcm'];
    keys.oPath = [resdir filesep 'points.jcm'];
    keys.px = pointList(:,1);
    keys.py = pointList(:,2);
    keys.pz = pointList(:,3);
    keys.quantity = sc.quantity;
end

jcmpSc = jcmpSc(1:end-1);
job_ids = zeros(1,nPoints,'int32');

% evaluate first point to get common PMLs if there is no existing log file
logFile = sc.pml; idx = 0; waitForLogfile = true;
if isempty(logFile)
    logFile = [sc.projectDir filesep 'pml.log'];
end
if ~exist(logFile, 'file')
    keys.pml = sprintf('LogFile = "%s"', logFile);
    keys.omega = complex(contours{1}(1));
    if isfield(keys, 'relPerm') && isa(keys.relPerm, 'function_handle')
        keys.relPermittivity = keys.relPerm(keys.omega);
    end
    wdir = [scDir filesep resultbag.get_tag(keys)];
    job_ids(1) = jcmwave_solve(jcmpSc,keys,resultbag,'workingdir',wdir);
    if ~job_ids(1), waitForLogfile = false; end
    idx = idx+1;
end

% evaluate the remaining points
keys.pml = sprintf('InputFile = "%s"', logFile);
while ~exist(logFile, 'file') && waitForLogfile, pause(1); end
tags = cell(1,nContours);
for it1 = 1:nContours
    tags_it1 = cell(1,numel(contours{it1}));
    for it2 = 1:numel(contours{it1})
        keys.omega = complex(contours{it1}(it2));
        if isfield(keys, 'relPerm')&&isa(keys.relPerm, 'function_handle')
            keys.relPermittivity = keys.relPerm(keys.omega);
        end
        tags_it1{it2} = resultbag.get_tag(keys);
        wdir = [scDir filesep tags_it1{it2}]; idx = idx+1;
        job_ids(idx) = jcmwave_solve(jcmpSc,keys,resultbag, ...
            'workingdir',wdir);
    end
    tags{it1} = tags_it1;
end

% delete temporary jcmpt file
if ~isempty(pointList)
    id = fopen([jcmpSc 't'],'w');
    fwrite(id,jcmpt,'char');
    fclose(id);
end
jcmwave_daemon_wait(job_ids, resultbag);
end
