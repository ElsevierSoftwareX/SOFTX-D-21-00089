% plot an expansion

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021
function fig = plotQuantity(RP, vars)
q = {}; q_ref = {}; fig = []; bulk = [];
requiredq = [true true]; M = RP.nModes+3;
if any(vars.index)
    requiredq(1) = any(vars.index<M);
    requiredq(2) = any(vars.index==M);
end
% get data
if requiredq(1)
    if isfield(RP.expansions, vars.quantity)
        q = RP.expansions.(vars.quantity);
        if size(q{1},2)>1, q = cellfun(@(x){real(x(:,1))},q); end
    elseif strcmpi(vars.quantity,'NormalizedDecayRate')
        if vars.pointEvaluation && isfield(RP.expansions,'DipoleEmission')
            q = RP.expansions.DipoleEmission;
        elseif isfield(RP.expansions,'ElectromagneticFieldEnergyFlux')
            q = RP.expansions.ElectromagneticFieldEnergyFlux;
        end
        if ~isempty(q)
            try
                bulk = RP.dipoleBulkEmission(RP.expansionPoints);
            catch
                error('No normalization available')
            end
            q = cellfun(@(x)x./bulk,q,'UniformOutput',false);
        end 
    elseif strcmpi(vars.quantity,'PhotonCollectionEfficiency')
        rad = RP.expansions.Radiation;
        if size(rad{1},2)==2
            nm = cellfun(@(x){x(:,2)},rad); 
            nm = real(nm{end}-sum(cat(2,nm{1:end-1}),2));
            q = cellfun(@(x){real(x(:,1))./nm},rad);
        end
    elseif vars.warning
        fprintf('No expansion available.\n');
    end
end
if requiredq(2)
    refq = vars.quantity;
    dp = {'NormalizedDecayRate' 'ElectromagneticFieldEnergyFlux'};
    rad = {'DipolePowerCollectionEfficiency' 'PhotonCollectionEfficiency'};
    refdp = isfield(RP.reference,'DipoleEmission');
    if any(strcmp(refq,dp))&&refdp, refq = 'DipoleEmission'; end
    if isfield(RP.w0,vars.quantity)
        w0_ref = RP.w0.(vars.quantity);
    elseif any(strcmpi(vars.quantity,dp))&&isfield(RP.w0,'DipoleEmission')
        w0_ref = RP.w0.DipoleEmission;
    elseif any(strcmpi(vars.quantity,rad))&&isfield(RP.w0,'Radiation')
        w0_ref = RP.w0.Radiation;
    end
    if isfield(RP.reference, refq)
        q_ref = RP.reference.(refq);
        if strcmp(vars.quantity,'NormalizedDecayRate') 
            if isempty(bulk) || length(q_ref)~=length(q{end})
                try
                    bulk = RP.dipoleBulkEmission(w0_ref);
                catch
                    error('No normalization available')
                end
            end
            q_ref = q_ref./bulk;
        elseif size(q_ref,2)>1
            q_ref = q_ref(:,1);
        end
    elseif strcmp(refq,rad{2}) && size(RP.reference.Radiation,2)==2
        q_ref = RP.reference.Radiation(:,1)./RP.reference.Radiation(:,2);
    elseif vars.warning
        fprintf('No reference solution available.\n');
    end
    if ~isempty(q_ref)
        xref = RieszProjection.unit(w0_ref,vars.xunit); 
    end
end
% set indices
emptyq = [isempty(q) isempty(q_ref)]; if all(emptyq), return; end
emptyq = emptyq & requiredq;
if vars.index, inds = vars.index; else, inds = 1:M-(~requiredq(2)); end
if emptyq(1),inds=inds(inds==M); elseif emptyq(2), inds=inds(inds<M); end
if ~vars.fignumber
    vars.fignumber=vars.quantityindex;
end
% get ylabel from quantity name
a = find(ismember(vars.quantity,'A':'Z'));
a = [1 a(a>1) length(vars.quantity)+1];
qname = join(arrayfun(@(x,y)vars.quantity(x:y),a(1:end-1),a(2:end)-1,...
    'UniformOutput',false));
% unitcomversions
[x,xu] = RieszProjection.unit(RP.expansionPoints,vars.xunit);
if strcmp(xu,'Hz')
    wl = '$\omega_0$ '; 
else
    wl = '$\lambda_0$ ';
end
% plot
leg = cell(1,length(inds)-length(vars.add2bg(vars.add2bg~=RP.nModes+1)));
it2 = 1; fig = figure(vars.fignumber);
clf(fig)
ax1 = axes(fig); 
ax1.TickLabelInterpreter = vars.opts{2}; ax1.FontSize = vars.opts{4};
hold on;
if any(inds<M)
    final = zeros(size(q{1}));
    for it1 = 1:RP.nModes
        if ismember(it1,vars.add2bg)
            q{end} = q{end}-q{it1}; continue;
        elseif ~ismember(it1,inds)
            continue; 
        end
        final = final - real(q{it1});
        plot(x, -real(q{it1}), '.-')
        leg{it2} = sprintf('Mode %d',it2); it2=it2+1;
    end
end
if ismember(RP.nModes+1,inds)
    final = final + real(q{end});
    plot(x, real(q{end}), '--'); leg{it2} = 'Background'; it2=it2+1;
end
if ismember(RP.nModes+2,inds)
    plot(x, final, '.-'); leg{it2} = 'Sum';
end
if ismember(M,inds)
    plot(xref, real(q_ref), 'o');
    leg{end} = 'Reference';
end
hold off
if any(vars.ylim), ylim(vars.ylim); end
if any(vars.xlim), xlim(vars.xlim); end
if strcmpi(vars.interpreter,'latex')
    xlabel(ax1,[wl '[$\mathrm{' xu '}$]'],vars.opts{:})
    if ~isempty(vars.yunit)
        vars.yunit = ['$\mathrm{' vars.yunit '}$'];
    end
else
    wl = split(wl,{'$'});
    xlabel(ax1,[wl{2} ' [' xu ']'],vars.opts{:})
end
if ~isempty(vars.yunit), qname = [qname{1} ' [' vars.yunit ']']; end
ylabel(qname,vars.opts{:})
if length(leg)<2, return; end
legend(leg,'location',vars.legendLocation,vars.opts{:}); legend('boxoff')
end
