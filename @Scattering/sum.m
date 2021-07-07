function out = sum(sc,tags,weights,omega0,outName)
%SUM Compute a weighted sum of field bags or field values
%   out = SUM(sc, tags, weights, omega0, outName) returns the weighted sum 
%   of the fields represented by the tags. The post process
%   <a href="matlab:web(...
%   ['https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/' ...
%   '721244647add8d817b99d2e4661c6ca8.html'], '-browser')"
%   >Superposition</a> of JCMsuite is used. If specified the outNames
%   determine the output filename. Omega0 is the real expansion frequency.
%   The job ids and the corresponding tags are returned.
%
%   see also RieszProjection, jcmwave_resultbag,
%            jcmwave_solve, jcmwave_daemon_wait

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

% check if a resource has been set up for computation
if isempty(jcmwave_daemon_resource_info)
    error('Please start a daemon');
end

% resultbag for superposition post process
if ~isfield(sc.pprbs,'superposition')
    filename = [sc.resultbagDir filesep 'superposition.mat'];
    sc.pprbs.('superposition') = jcmwave_resultbag(filename); 
end
resultbag = sc.pprbs.('superposition');
if nargin<4, outName = ''; end

% get resultbag with scattering results
rbsc = sc.resultbags{1+isfield(sc.resultbags{2}.results_,tags{1})};

% create jcmpt file
wDir = [sc.projectDir filesep 'superpositions'];
jcmpt = [wDir filesep 'superposition.jcmpt'];
if ~exist(jcmpt,'file')
    superposition = sc.superposition;
    id = fopen(jcmpt, 'w');
    fprintf(id, superposition);
    fclose(id);
end
jcmpt = jcmpt(1:end-1);

fieldBags = cell(size(tags));
md = java.security.MessageDigest.getInstance('MD5');
hash = md.digest(double([tags{:}]));
bi = java.math.BigInteger(1, hash);
md5 = char(bi.toString(16));
for it = 1:length(tags)
    field = rbsc.results_.(tags{it}).result{1}.file;
    fieldBags{it} = sprintf(sc.fieldBag, field, ...
        real(weights(it)), imag(weights(it)));
end
fieldBags = [fieldBags{:}];
if isempty(outName)
    outDir = [wDir filesep md5];
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    outName = [outDir filesep ['s' num2hex(omega0)]];
    wDir = outName; outName = [outName 'fieldbag.jcm'];
else
    a = strfind(outName,filesep); wDir = outName(1:a(end-1)-1);
end
keys.oPath = outName; keys.fieldbags = fieldBags;
keys.omega = omega0; keys.md5 = md5; tag = resultbag.get_tag(keys);
out = {jcmwave_solve(jcmpt,keys,resultbag,'workingdir',wDir);tag};
end
