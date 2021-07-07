function [ids,tags] = submit(sc,tags,varargin)
%SUBMIT Apply post processes to the fields
%   [ids,tags] = SUBMIT(sc,tags,quantity) A post process of JCMsuite is 
%   applied to each field asigned to a tag in tags. The input parameter 
%   'quantity' must match the name of the output quantity and the filename
%   of the post process where it is defined in.
%   If e.g. you have quantity = 'ElectricFieldEnergy', there must be a
%   template file 'ElectricFieldEnergy.jcmpt' which, making use of
%   parameter substitution, contains the parameter definitions 
%   FieldBagPath = "%(field)s" and OutputFileName = "%(oPath)s". In the
%   definition of the <a href="matlab:
%   web(['https://docs.jcmwave.com/JCMsuite/html/ParameterReference' ...
%   '/d904047e6b00d14da719fb37a669e7eb.html'],'-browser')" 
%   >post process</a>. Its results are summed over all domain
%   ids and returned as a doulbe array.
%   As the reference solution can usually be obtained with a native
%   JCMsuite post process whereas for the expansion a python expression
%   must be evaluated, input parameter and filename for reference
%   solutions may differ. In the above example the input would change to
%   quantity = 'referenceElectricFieldEnergy' and the file name to 
%   'referenceElectricFieldEnergy.jcmpt'.
%   The ouput is an instance of this class which must be passed to the
%   method 'collect' in order to get the numeric results. You can call
%   'submit' many times before you call 'collect' for the first time. This
%   way parallelization is improved. The returned tags refer to the
%   resultbags wich hold the results (sc.pprbs.(quantity)) and the ids to
%   jobids.
%
%   [ods,tags] = SUBMIT(fb1,fb2,tags,quantity) If there are two field bags
%   the corresponding post process must additionally contain the template
%   parameter circfield i.e. FieldBagPath = "%(circfield)s".
%
%   [ids,tags] = SUBMIT(...,keys) Any of the above parameter combinations
%   can be combined with a parameter 'keys', which is passed to 
%   jcmwave_solve. 
%
%   If a python expression is used for integration, the template file of
%   the post process must, within the section Python {...}, include the
%   definition IntegrationOrder = %(integrationOrder)e. This parameter
%   will be set automatically to twice the finite element degree of the
%   FEM solutions. In this section you can also define the name of the
%   output quantity, e.g., IntegralName = "ElectricFieldEnergy". A python
%   expression is currently always needed if you want to expand a quadratic 
%   functional like the energy of the electric field. As complex 
%   conjugation, being not holomorphic, does not allow for a straight
%   forward application of this method, <a href="matlab: web(...
%   'https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.035432',...
%   '-browser')">circ fields</a> are used. Instead of
%   dot(E,E^*), dot(E,E^o) is computed, with E^o(w) = E(-w). To envoke the
%   field bag of a circ field within a python expression you can define
%   FieldBagPath = "%(circfield)e". 
%   Please also refer to the <a href="matlab:
%   web(['https://docs.jcmwave.com/JCMsuite/html/ParameterReference' ...
%   '/d904047e6b00d14da719fb37a669e7eb.html'],'-browser')" 
%   >parameter reference</a> of JCMsuite and the 
%   documentation of the <a href="matlab:web(['https://docs.jcmwave.' ...
%   'com/JCMsuite/html/MatlabInterface/index.html'],'-browser')" 
%   >matlab interface</a> if you are interested in the 
%   definitions of the post processes in JCMsuite. The corresponding
%   project files must be stored in a directory whose path is passed to
%   the constructor of this class. Its default name is 'postprocesses' and
%   it is expected to be located in the same folder as @Scattering. Please
%   also refer to the existing post processes.
%
%   see also RieszProjection, Scattering,
%            jcmwave_resultbag, jcmwave_solve

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

% check if a resource has been set up for computation
if isempty(jcmwave_daemon_resource_info)
    error('Please start a daemon');
end

% define defaults
tags_ = []; q = ''; keys = struct;
for it = 1:length(varargin)
    switch class(varargin{it})
        case 'cell'
            tags_ = varargin{it};
        case {'char','string'}
            q = varargin{it};
        case 'struct'
            keys = varargin{it};
    end
end

% make sure that the input parameters are as expected
if ~isa(tags, 'cell')
    error('Cell expected but got %s', class(tags));
elseif ~isempty(tags_) && length(tags)~=length(tags_)
    error('Dimension missmatch (%d,%d)',length(tags),length(tags_))
end
if isempty(q)
    error('You must provide the quantity to be evaluated')
end

% the dipole emission does not require submitting jobs to jcmwave_solve
% unless the reference is required
if strcmp(q,'DipoleEmission'), ids = []; return; end

% some quantities need special default keys
if ismember(q,{'Radiation','referenceRadiation'})
    keys = sc.radiationKeys(keys); 
elseif ismember(q,{'RadiationPattern' 'referenceRadiationPattern'})
    keys = sc.radiationPatternKeys(keys);
end

% get resultbags
if ~isfield(sc.pprbs,q)
    filename = [sc.resultbagDir filesep [q '.mat']];
    relevantFields = [fieldnames(keys);'field';'omega'];
    if ~isempty(tags_), relevantFields = [relevantFields;'circfield']; end
    sc.pprbs.(q) = jcmwave_resultbag(filename,sort(relevantFields));
end
rb = sc.pprbs.(q);

% get resultbag with scattering results
rbsc = sc.resultbags{1+isfield(sc.resultbags{2}.results_,tags{1})};

% check if the corresponding post process is available
pDir = sc.postproDir; jcmpt = [q '.jcmpt'];
if ~exist(pDir, 'dir')
    error('The directory %s does not exist',pDir);
end
d = dir(pDir); d = d(~[d.isdir]);
indx = strcmpi({d.name},jcmpt);
if startsWith(jcmpt,'reference')
    if ~any(indx)
        jcmpt = jcmpt(10:end); indx = strcmpi({d.name},jcmpt);
    end
end
if ~any(indx)
    error('Scattering:evaluate:missingPostProcess',...
        'Please provide a file %s located in:\n%s',jcmpt,pDir)
end
jcmpt = [pDir filesep d(indx).name(1:end-1)];

% submit jobs to jcmwave_solve
ids = zeros(size(tags), 'int32');
for it1 = 1:length(tags)
    keys.integrationOrder = 2*sc.keys.finiteElementDegree+1;
    keys.oPath = 'out.jcm';
    keys.omega = rbsc.results_.(tags{it1}).keys.omega;
    keys.field = rbsc.results_.(tags{it1}).result{1}.file;
    if ~isempty(tags_)
        keys.circfield = rbsc.results_.(tags_{it1}).result{1}.file;
    end
    ids(it1) = jcmwave_solve(jcmpt,keys,rb,'temporary','yes');
    tags{it1} = rb.get_tag(keys);
end
tags = tags(1:it1);
end

