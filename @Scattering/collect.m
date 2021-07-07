function v = collect(sc,ids,tags,q,last)
%COLLECT Collect results of jobs started by a call of the method 'submit'
%   v = COLLECT(tags) The input argument 'tags' must be the results of a
%   call to the method submit. The returned array v is of shape (n,m) with
%   n being the data dimension and m being the length of the input
%   argument.
%
%   see also submit, RieszProjection, Scattering,
%            jcmwave_resultbag, jcmwave_solve

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

q_ = q; if startsWith(q,'reference'), q = q(10:end); end
dE = strcmp(q,'DipoleEmission');
if isempty(ids) && dE
    rbsc = sc.resultbags{1+isfield(sc.resultbags{2}.results_,tags{1})};
    v = zeros(size(tags));
    for it = 1:length(tags)
        fvs = rbsc.results_.(tags{it}).result{2}.ElectricFieldStrength{1};
        v(it) = -0.5*sum(fvs.*sc.keys.strength(:).');
    end
    return;
end
rb = sc.pprbs.(q_); 
jcmwave_daemon_wait(ids,rb);
v = cell(size(ids)); bE = []; if dE, bE = zeros(length(v),2); end
for it = 1:length(v)
    tag = tags{it};
    try
        % Extract numeric results
        % v{it3} must be of shape (:,1)
        if strcmp(q,'Radiation')
            v{it} = integrateFarfield(tag,last&&it==length(v));
        elseif strcmp(q,'RadiationPattern')
            v{it} = pointingVector(tag);
        else
            v1 = rb.results_.(tag).result{end};
            if dE
                v{it} = v1.DipoleEmission;
                bE(it,2) = v1.DipoleBulkEmission; %#ok<AGROW>
                bE(it,1) = rb.results_.(tag).keys.omega; %#ok<AGROW>
            elseif any(startsWith(fieldnames(v1),'DomainId'))
                v{it} = sum(v1.(q){1}); % sum over all domains
            else
                v{it} = v1.(q){1}(:); % flatten the result
            end
        end
    catch ME
        err = rb.results_.(tag).log.Log.Error;
        if ~isempty(err), error(err); end
        rethrow(ME);
    end
end
v = cell2mat(v);

if ~isempty(bE)
    sc.bulkEmission = bE;
end

% Evaluate surface integral of results with given tag in the
% corresponding resultbag whose memory usage is reduced replacing
% the field values with their integral. The post process exports
% the far-field on a surface defined by the a NA. The points are
% distributed according to the abcissae of a Gauss-Kronrod
% quadrature rule.
    function out = integrateFarfield(tag,last)
        persistent results resize;
        if isempty(resize), resize = false; results = rb.results_; end
        farfield = rb.results_.(tag).result;
        % Check if the integration has already been done
        if isnumeric(farfield)
            out = farfield;
        else
            keys_ = rb.results_.(tag).keys;
            pwt = RieszProjection.weights15;
            WT = [pwt(end:-1:2), pwt];
            r = keys_.radius; ndx = keys_.ndx;
            theta = reshape(keys_.points_the,15,[]);
            sz = [size(theta) length(keys_.points_phi) 3];
            field1 = reshape(farfield{1}.ElectricFieldStrength{1},sz);
            try field2 = reshape(farfield{2}.ElectricFieldStrength{1},sz);
            catch, field2 = field1; end
            hlfhp = deg2rad(keys_.halfhp(1)); 
            hlfht = deg2rad(keys_.halfht); c1 = r^2*hlfhp;
            eps = farfield{1}.header.RelPermittivity*sc.eps0;
            mu = farfield{1}.header.RelPermeability*sc.mu0;
            S = 0.5*sqrt(eps/mu)*sum(field1.*conj(field2),4).*sind(theta);
            P = c1.*sum(S.*repmat(reshape(WT,1,1,[]),1,1,size(S,3)/15),3);
            out = sum(P(:,1:ndx).*WT.'.*hlfht(1:ndx),[1 2]);
            if size(P,2)>ndx, out(2,1) = sum(P.*WT.'.*hlfht,[1 2]); end
            resize = true; 
        end
        results.(tag).result = out; % replace fields whith integrals
        if last % reload resultbag with smaller size
            if resize
                fprintf('\nPrepare resized resultbag.')
                results.fieldnames = rb.fieldnames_;
                results.keys = struct;
                results.source_files = rb.source_files_;
                delete(rb.filepath_)
                save(rb.filepath_,'-struct','results')
                [~,quantity,~] = fileparts(rb.filepath_);
                sc.pprbs.(quantity) = jcmwave_resultbag( ...
                    rb.filepath_, rb.fieldnames_); fprintf(1,'\n');
            end
            results = []; resize = [];
        end
    end

% evaluate pointing vector from far-field post process
    function out = pointingVector(tag)
        farfield = rb.results_.(tag).result;
        field1 = farfield{1}.ElectricFieldStrength{1};
        try
            field2 = farfield{2}.ElectricFieldStrength{1};
        catch
            field2 = field1;
        end
        eps = farfield{1}.header.RelPermittivity*sc.eps0;
        mu = farfield{1}.header.RelPermeability*sc.mu0;
        out = 0.5*sqrt(eps/mu)*sum(field1.*conj(field2),2);
        out = out(:);
    end

end
