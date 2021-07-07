% plot the error of the contour integrations

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021
function fig = plotError(RP, vars)
exclude = {'DipolePowerCollectionEfficiency',...
    'PhotonCollectionEfficiency'};
if ~any(vars.index), vars.index = RP.nModes+2; end
errortype = lower(vars.errortype(1));
if errortype == 'r', label_y = 'Maximal Relative Error';
elseif errortype == 'a', label_y = 'Maximal Absolute Error';
else, label_y = 'Error'; 
end
label_x = 'Number of Integration Points';
if all(vars.index==RP.nModes+2)
    label_x = [label_x ' (Background Contour)'];
end
if isempty(vars.quantities{1})
    quantities = fieldnames(RP.expansions); 
    quantities = quantities(~ismember(quantities,exclude));
else
    quantities = vars.quantities;
end
if isempty(quantities), fig = []; return; 
elseif all(vars.index==RP.nModes+2) && isempty(fieldnames(RP.reference))
    vars.index = RP.nModes+1; label_y = [label_y ' (Background)'];
end
if ~vars.fignumber, vars.fignumber=2; end
fig = figure(vars.fignumber);
clf(vars.fignumber)
ax = axes(fig); ax.NextPlot = 'add';
ax.TickLabelInterpreter = vars.opts{2}; ax.FontSize = vars.opts{4};
set(ax, 'YScale', 'log');
xlabel(label_x,vars.opts{:})
ylabel(label_y,vars.opts{:})
colors = ax.ColorOrder;
markers = ['o';'s';'d';'^';'v';'>';'<'];
MC = {'Marker', 'Color'}; mc = {markers colors};
if length(quantities)==1, MC = MC(end:-1:1); mc = mc(end:-1:1); end
lines = cell(1,length(quantities)+length(vars.index)); idx = 1;
lines = lines(1:end-length(quantities)==1-length(vars.index)==1);
qnames = cell(1,length(quantities)); conv_ = NaN(1,0);
for it1 = 1:length(quantities)
    n1 = MC{2}; v1 = mc{2}(idx,:);
    q = quantities{it1};
    [conv,x] = error(RP, q, vars.index, errortype); 
    if isempty(conv), continue; end
    useBarPlot = isscalar(vars.index) && length(x{1})<3;
    [~, ind] = max(cellfun(@length, x));
    a = find(ismember(q,'A':'Z')); a = [1 a(a>1) length(q)+1];
    qnames(it1) = join(arrayfun(@(x,y){q(x:y)},a(1:end-1),a(2:end)-1));
    if length(quantities)>1 && ~useBarPlot
        lines{idx} = plot(ax,NaN,NaN,'o-',n1,v1,'DisplayName',qnames{it1});
        idx = idx+1;
    end
    if it1==1
        step = max(x{ind}(end)/2^3,2);
        xticks(step:step:x{ind}(end));
    end
    for it2 = 1:length(vars.index)
        if it2>length(conv), break; end
        n2 = MC{1}; v2 = mc{1}(mod(it2-1,7)+1,:);
        if isempty(conv_) && useBarPlot
            conv_ = NaN(2,length(quantities));
            names = cell(1,length(quantities));
            x_ = x; if isscalar(x{1}), x_ = {[NaN x{1}]}; end
            conv_(end-1+isscalar(conv{1}):end,1) = conv{1};
            names(1) = qnames(it1); ndx = 2;
        elseif useBarPlot
            conv_(end-1+isscalar(conv{1}):end,ndx) = conv{1};
            names(ndx) = qnames(it1); ndx = ndx+1;
        else
            plot(ax,x{it2},conv{it2},'-',n1,v1,n2,v2);
        end
        if it1==length(quantities) && length(vars.index)>1
            if it2<RP.nModes+1
                nm = sprintf('Mode %d',it2);
            else
                nm = 'Background';
            end
            if it1>1, v3 = 'k'; a = 'o'; else, v3 = v1; a = 'o-'; end
            lines{idx} = plot(ax,NaN,NaN,a,n1,v3,n2,v2,'DisplayName',nm);
            idx = idx+1;
        end
    end
end
if ~isempty(conv_)
    n = vars.fignumber;
    if ~isempty(ax.Children), n = n + 110*n; end
    fig_ = figure(n); clf(n);
    ax_ = axes(fig_);
    conv_(:,ndx:end) = []; names(ndx:end) = []; % remove empty content
    bar(ax_,x_{1},conv_)
    if any(isnan(x_{1})), ax_.XTick = x_{1}(2); end
    set(ax_, 'YScale', 'log');
    ylabel(label_y,vars.opts{:})
    ylim([min(conv_,[],'all')*2e-1 max(conv_,[],'all')*11])
    legend(ax_,names,vars.opts{:}); legend('boxoff')
    xlabel(ax_,label_x,vars.opts{:})
    ax_.TickLabelInterpreter = vars.opts{2};
    ax_.FontSize = vars.opts{4};
    figure(vars.fignumber); fig = {fig fig_};
end
if ~isempty(lines)
    legend([lines{:}],'location',vars.legendLocation,vars.opts{:})
    legend('boxoff')
end
if any(vars.ylim), ylim(vars.ylim); end
if any(vars.xlim), xlim(vars.xlim); end
hold off
end
