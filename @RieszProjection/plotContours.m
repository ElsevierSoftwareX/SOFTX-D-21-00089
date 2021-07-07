% plot contours, expansion points and eigenvalues

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021
function fig = plotContours(RP, vars)
if vars.maxIndex > 1
    nModes=length(RP.selectedPoles);
else
    nModes = 0; % if there is not any contour yet or no poles are selected
end
% plot
if ~vars.fignumber, RP.plotSettings.fignumber.contours = 1; end
fig = figure(RP.plotSettings.fignumber.contours);
clf(vars.fignumber);
ax1 = axes(fig); ax1.NextPlot = 'add'; ax1.Box = 'on'; axis(ax1,'equal')
ax1.TickLabelInterpreter = vars.opts{2}; ax1.FontSize = vars.opts{4};
if strcmpi(vars.interpreter,'latex')
    xlabel(ax1,'Re($\omega$) [Hz]',vars.opts{:})
    ylabel(ax1,'Im($\omega$) [Hz]',vars.opts{:})
else
    xlabel(ax1,'Re(\omega) [Hz]',vars.opts{:})
    ylabel(ax1,'Im(\omega) [Hz]',vars.opts{:})
end
if vars.index
    inds = vars.index;
    if all(inds==nModes+1), inds = 1:nModes+1; end
else
    inds = 1:nModes+1;
end
if RP.includeConjugatedPoles
    inds = [inds(inds<nModes+1) nModes+inds];
end
if ~isempty(RP.poles)
    if RP.includeConjugatedPoles
        conjPoles = conj(RP.poles);
        if ~isempty(RP.selectedPoles) && RP.selectedPoles{1}(1)
            conjPoles = conjPoles([RP.selectedPoles{:}]);
            plot(ax1,conjPoles,'b+','markersize',4,'linewidth',1,...
                'displayname','$\omega_n^\ast$')
        end
    end
    l1 = plot(ax1,RP.poles, 'b.','DisplayName','$\omega_n$');
    l1.MarkerSize = 12;
end
if ~isempty(RP.expansionPoints)
    l2 = plot(ax1,RP.expansionPoints, zeros(size(RP.expansionPoints)),...
        'k.','DisplayName','$\omega_0$');
    l2.MarkerSize = 12;
end
for it = inds
    try contour = RP.contours{it}; catch, break; end
    if isempty(contour), continue; end
    if it==inds(end)
        if isreal(contour)
            plot(ax1,contour,0,'ro','DisplayName','Integration Points');
        else
            plot(ax1,contour(:),'r.','DisplayName','Integration Points');
        end
    else
        plot(ax1,contour(:),'r.','DisplayName','');
    end
end
selection = ax1.Children;
if length(selection)>(3+RP.includeConjugatedPoles)
    selection = selection([1 end-1-RP.includeConjugatedPoles:end]);
end
legend(ax1,selection,'location',vars.legendLocation,vars.opts{:})
legend(ax1,'boxoff')
cindex = max(vars.index); 
if cindex>length(RP.selectedPoles), cindex = length(RP.contours); end
if any(vars.ylim)
    ylim(vars.ylim)
elseif vars.index
    M = max(imag(RP.contours{cindex})); 
    m = min(imag(RP.contours{cindex}));
    shift = (M-m)/10; ylim([m-shift M+shift]);
end
if any(vars.xlim)
    xlim(vars.xlim)
elseif vars.index
    shift = (M-m)/2;
    M = max(real(RP.contours{cindex}));
    m = min(real(RP.contours{cindex}));
    shift = max(shift-(M-m)/2,(M-m)/10); xlim([m-shift M+shift]);
end
if isempty(RP.contours_), return; end

pos = get(ax1, 'Position');
if vars.contourNumbers && ~isempty(RP.selectedPoles)
    ax2 = axes('Position',pos,'XAxisLocation','top',...
        'Color','none','XColor','k','YColor','none');
    ax2.TickLabelInterpreter = vars.opts{2}; ax2.FontSize = vars.opts{4};
    axis(ax2,'equal')
    xlim(get(ax1,'xlim')); ylim(get(ax1,'ylim'));
    inds = 1:length(RP.selectedPoles);
    if all(RP.selectedPoles{1})
        xticklabels(arrayfun(@num2str,inds,'UniformOutput',false))
        ps = cellfun(@(x) x(1), RP.selectedPoles);
        xticks(ax2,real(RP.poles(sort(ps(inds)))));
        yticks(ax2,[]);
    end
    linkaxes([ax1 ax2],'xy'); 
    fig.CurrentAxes = ax1; ax1.Box = 'off';
end
hold off
end
