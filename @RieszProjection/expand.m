function expand(RP, varargin)
%EXPAND Perform a modal expansion
%   EXPAND(RP) All quantities available for an expansion are listed and can
%   be selected with the mouse. A left click will call this method with
%   default parameters. This method is optimized to work with the interface
%   to JCMsuite, yet can be used independently. 
%
%   EXPAND(RP, quantity) An expansion of the given quantity is evaluated.
%   The interface which has been passed to the constructor, must be able to
%   derive this quantity from the solutions of the linear systems.
%   The argument 'quantity' can either be a string or a cell array 
%   containing strings. If the property 'precision' is set, the number of
%   points is refined adaptively until the maximal error is sufficiently
%   small or the number of points on the corresponding contour is larger
%   than the value of the keyword argument 'upperbound', whose default is
%   1024. The property 'precision' can be either a scalar or a vector with
%   two elements. In the latter case, the first value is used as the bound 
%   of the absolute error and the second one as the bound of the relative
%   error. For absolute and relative error control you can set the first or 
%   second value to zero, respectively. If you pass the keyword argument 
%   precision, it is used instead of the property. If the precision is Inf
%   no refinements will be made.
%   The result is saved as a field of the property 'expansions' with the
%   same name as the quantity. If the property 'referencePoints' is not
%   empty, the scattering problems at these points will be evaluated and
%   used as reference solutions for the evaluation of the error.
%   
%   EXPAND(...,PARAM1,VAL1,PARAM2,VAL2,...) The above arguments can be
%   followd by key value pairs, which include:
%      precision (numeric): Determines when to stop the adaptive
%         refinement. The default is the value of the property of the same
%         name, which is Inf by default.
%      upperbound (numeric): Fixes the maximal number of integration points
%         a single contour can have. The default is 1024.
%      keys (struct): Keys which will be passed to corresponding post
%         processes, e.g., for parameter substitutions. In case of the
%         radiation post process this can, e.g., be the NA. Be aware that
%         if the input parameter 'quantity' is a cell array, you must
%         provide a cell array of keys. The keys corresponding to a given
%         quantity are stored and will be reused if expand is called for
%         the same quantity again.
%      reference (logical): If true scattering problems at reference
%         or expansion points are evaluated. Make sure that the
%         corresponding post processes exist. The default is false unless
%         you have set the property 'referencePoints'.
%
%   see also RieszProjection

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

persistent parser;
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'RieszProjection/expand';
    addOptional(parser,'quantity','',@validatequantity);
    addParameter(parser,'precision',0,...
        @(x)validateattributes(x,{'numeric'},{'nonnegative','scalar'}));  
    addParameter(parser,'upperbound',2^10,...
        @(x)validateattributes(x,{'numeric'},{'nonnegative','scalar'}));
    addParameter(parser,'keys',struct,...
        @(x)validateattributes(x,{'struct'},{'scalar'})); 
    addParameter(parser,'reference',false,...
        @(x)validateattributes(x,{'logical'},{'scalar'})); 
    addParameter(parser,'recalculateEvs',false,...
        @(x)validateattributes(x,{'logical'},{'scalar'}));
end

parse(parser,varargin{:});
vars = parser.Results;

% Check if there are any expansion points defined
if ~isempty(RP.expansionPoints)
    w0 = RP.expansionPoints;
elseif ~isempty(RP.referencePoints)
    w0 = RP.referencePoints;
else
    fprintf(1,'No expansion points defined\n')
    return;
end
if ~vars.precision, vars.precision = RP.precision; end

% if RP.referencePoints is not empty, reference solutions will be evaluated
if ~isempty(RP.referencePoints), vars.reference = true; end
% validate the quantities provided by the user
allQs = RP.getQuantities; 
if isempty(vars.quantity), printPossibilities(inputname(1)); return; end
expand_(vars.quantity,vars.precision,vars.recalculateEvs)

plothandles = struct2cell(RP.plotSettings.fignumber);
for h = [plothandles{:}]
    if ~ishghandle(h), continue; end
    qi = get(h,'UserData');
    if iscell(qi) && ~strcmp(qi{1},'contours')
        RP.plot(qi{1:2})
    end
end

    function expand_(quantity,precision,recalculateEvs)
        [quantity,useReferencePoints] = validateQuantities(quantity);
        parameters = {w0,'precision',precision,'upperbound',...
            vars.upperbound,'recalculateEvs',recalculateEvs};
        for it1 = 1:length(quantity)
            q = quantity{it1};
            if it1==useReferencePoints
                parameters{1} = RP.referencePoints; 
            end
            % evaluate reference solution
            if vars.reference
                try RP.getReferenceSolution(q,RP.referencePoints);
                catch me
                    fprintf(1,['\nCould not evaluate reference '...
                        'solution for %s\n%s\nIn %s line: %d\n'],q,...
                        me.message,me.stack(1).name,me.stack(1).line)
                end
            end
            RP.expansions.(q) = RP.integrate(q,parameters{:});
            % update the existing quantities (without refinement)
            if ~RP.up2date(2) 
                q_ = fieldnames(RP.expansions);
                q_ = q_(~strcmp(q_,q)); RP.up2date(2) = true;
                if ~isempty(q_),expand_(q_,Inf,false); end
            end
        end
        
        % This part of the nested function is JCMsuite specific
        if any(ismember({'Radiation','DipoleEmission'},quantity))
            if contains('DipoleEmission',quantity)
                try bulk = RP.dipoleBulkEmission(w0);
                catch, warning('No bulk emission available'); bulk = 0; end
                dpBg = RP.expansions.DipoleEmission{end};
                RP.expansions.DipoleEmission{end} = dpBg+bulk;
            end
            expansions = fieldnames(RP.expansions);
            if all(ismember({'Radiation' 'DipoleEmission'},expansions))
                % combine existing quantities
                dp = RP.expansions.DipoleEmission;
                dp = real(dp{end}-sum(cat(2,dp{1:end-1}),2));
                rad = RP.expansions.Radiation;
                dce = cellfun(@(x){real(x(:,1))./dp},rad);
                RP.expansions.DipolePowerCollectionEfficiency = dce;
            end
        end
    end


    
    % Print all quantities that can be expanded to the command window
    function printPossibilities(instancename)
        href = '<a href="matlab: %s.expand(''%%s'')">%%s</a>\n';
        href = sprintf(href,instancename);
        str = ['You can expand the following '...
            'quantities:\n' repmat(href,1,length(allQs))];
        qs = cell(1,2*length(allQs));
        qs(1:2:end-1) = allQs; qs(2:2:end) = allQs;
        str = sprintf(str,qs{:});
        str = str(1:end-1);
        disp(str);
    end

    % If you have set the properties 'linearQ' and 'quadraticQ' this
    % function will work, however it "knows" about the interface to
    % JCMsuite
    function [qs,useReferencePoints] = validateQuantities(qs)
        if ~iscell(qs), qs = {qs}; end
        if ~iscell(vars.keys), vars.keys = {vars.keys}; end
        uRefP = zeros(size(qs),'logical'); dP = false;
        for it = 1:length(qs)
            qs{it} = validatestring(qs{it},allQs,mfilename,'quantity');
            uRefP(it) = strcmp(qs{it},'RadiationPattern');
            if strcmp(qs,'NormalizedDecayRate')
                qs{it} = 'DipoleEmission'; 
            elseif contains(qs{it},'CollectionEfficiency')
                dP = strcmp('DipoleCollectionEfficiency',qs{it});
                qs{it} = 'Radiation';
            elseif strcmp(qs,'RadiationPattern'), uRefP(it) = true;
            end
            % set or remember keys of given quantity and initialize derived
            % quantities if necessary
            try keys = RP.qKeys.(qs{it}); catch, keys = struct; end
            if ismember('keys',parser.UsingDefaults)
                vars.keys{it} = keys; 
            end
            if ~(isfield(RP.qKeys,qs{it})&&isequal(keys,vars.keys{it}))
                RP.derivedQuantities.(qs{it}) = {};
                if vars.reference && isfield(RP.rKeys,qs{it})
                    RP.rKeys = rmfield(RP.rKeys,qs{it}); 
                    RP.w0 = rmfield(RP.w0,qs{it});
                end
            end
            % also accessed by getReferenceSolution
            RP.qKeys.(qs{it}) = vars.keys{it};
        end
        if dP, qs = [qs 'DipoleEmission']; end
        qs = unique(qs); useReferencePoints = length(qs)+1;
        if any(uRefP) 
            if isempty(RP.referencePoints)
                fprintf(1,['The expansion of the quantity %s is only '...
                    'supported at the reference points.\n'], qs{uRefP}); 
                qs(uRefP) = [];
            else
                ndx = find(uRefP); n = sum(uRefP)-1;
                qs([ndx end-n+1:end]) = qs([end-n+1:end ndx]);
                useReferencePoints = length(qs)-n;
            end
        end
    end
end

function validatequantity(quantity)
if ischar(quantity) || (iscell(quantity)&&all(cellfun(@ischar,quantity)))
    return;
end
error(['The input parameter quantity is expected to be a string' ...
    'or a cell array containing strings.'])
end
