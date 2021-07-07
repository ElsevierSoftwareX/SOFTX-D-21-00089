% plot the radiation pattern either 2d or 3d

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021
function fig = plotRadiationPattern(RP,vars)

% get the data and the keys to reconstruct the mesh in theta and phi 
rad = []; ref = []; keys = struct; keys_ref = [];
hz = strcmpi(vars.xunit,'hz'); % if x-unit is not Hz revert order
if isfield(RP.expansions,'RadiationPattern')
    if isfield(RP.qKeys,'RadiationPattern')
        keys = RP.qKeys.RadiationPattern;
    end
    rad = RP.expansions.RadiationPattern;
    if ~hz, rad = cellfun(@(x){x(end:-1:1,:)},rad); end
end
if isfield(RP.reference,'RadiationPattern')
    if isfield(RP.rKeys,'RadiationPattern')
        keys_ref = RP.rKeys.RadiationPattern;
    end
    if ~isempty(rad) && ~isempty(keys_ref)&&~isequal(keys,keys_ref)
        error('Keys used for expansion and reference values differ.')
    end
    ref = RP.reference.RadiationPattern;
    if ~hz, ref = ref(end:-1:1,:); end
    if ~isempty(rad) && size(rad{1},1)~=size(ref,1)
        error('Expansion points of reference solution must coincide.')
    end
    if ~isempty(keys_ref), keys = keys_ref; end
end
keys = Scattering.radiationPatternKeys(keys);
theta = keys.gridPointsTheta; phi = keys.gridPointsPhi;

% set defaults
if ~vars.index, index = 1:RP.nModes+3; else, index = vars.index; end
index(ismember(index,vars.add2bg(vars.add2bg~=RP.nModes+1))) = []; 

if length(phi)==2
    [theta,idx] = sort([-theta theta]);
    rad = cellfun(@(x){x(:,idx)},rad); ref = ref(:,idx);
    [theta,idx] = unique(theta); 
    rad = cellfun(@(x){x(:,idx)},rad); ref = ref(:,idx);
end
rad(1:end-1) = cellfun(@(x){-x},rad(1:end-1));

% get normalization
tot = sum(cat(3,rad{:}),3); if isempty(ref), ref = tot; end
m = size(tot,1); idx = zeros(1,RP.nModes+1,'logical');
idx(end) = true; idx(vars.add2bg) = true; bg = sum(cat(3,rad{idx}),3);
nm_ = max(ref(:)); rad = [rad(1:end-1) {bg tot ref}]; 

% plot
if vars.fignumber, figN = vars.fignumber; else, figN = 20; end
fig = figure(figN); clf; xL = [Inf -Inf]; yL = [Inf -Inf];
is3D = all([length(theta) length(phi)]>2); ax = cell(1,m); 
if is3D
    [Phi,Theta] = meshgrid(deg2rad(phi),pi/2-deg2rad(theta));
    sz = size(Theta); zL = [Inf -Inf];
else
    theta = deg2rad(theta)+pi/2; % 2d
end
ax_ = subplot(1,m+1,1); disableDefaultInteractivity(ax_) % for legend
ax_.CLim = [1 RP.nModes+3];
for it1 = 1:m
    ax{it1} = subplot(1,m+1,it1+1); ax{it1}.DataAspectRatio = [1 1 1];
    ax{it1}.NextPlot = 'add'; ax{it1}.Visible = 'off';
    % in case that the axis is made visible
    ax{it1}.TickLabelInterpreter = vars.interpreter;
    ax{it1}.FontSize = vars.fontSize; ax{it1}.CLim = [1 RP.nModes+3];
    ax{it1}.Clipping = 'off'; ax{it1}.Projection = 'perspective';
    % plot
    
    for idx = index
        if is3D
            R = real(reshape(rad{idx}(it1,:)./nm_,sz));
            [X,Y,Z] = sph2cart(Phi,Theta,R);
            s = mesh(ax{it1},X,Y,Z); s.CData = idx*ones(sz);
            s.FaceAlpha = 0.1; s.FaceColor = 'flat';
        else
            [X,Y] = pol2cart(theta,real(rad{idx}(it1,:))./nm_);
            l = plot(ax{it1},X,Y);
        end
        if idx==RP.nModes+2 && is3D
            s.FaceColor = 'none';
        elseif idx==RP.nModes+3
            if is3D
                s.EdgeColor = 'none'; s.FaceAlpha = 0.5;
            else
                l.LineStyle = ':'; l.LineWidth = 2;
            end
        end
        xL = [min(xL(1),min(X(:))) max(xL(2),max(X(:)))];
        yL = [min(yL(1),min(Y(:))) max(yL(2),max(Y(:)))];
        if is3D, zL = [min(zL(1),min(Z(:))),max(zL(2),max(Z(:)))]; end
    end
end

% Link the aspect ratios, reduce distances and put labels

ax = [ax{:}]; 
if is3D
    link = linkprop(ax,{'CameraPosition','CameraUpVector',...
        'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(fig, 'StoreTheLink', link);
    ax(end).ZLim = zL; ax(end).YLim = yL; ax(end).XLim = xL;
    ax(end).CameraPosition = [2.5 -2.5 2.5];
    ax(end).CameraTarget = [0 0 zL(2)/2];
    ax(end).CameraUpVector = [0 0 1];
else
    linkaxes(ax);
    ax(end).YLim = yL; ax(end).XLim = xL;
end
opts = ['Units','normalized','Visible','on',vars.opts];
[x,xu] = RieszProjection.unit(RP.referencePoints,vars.xunit);
x = sort(x); if hz, x = x*1e-15; end
if strcmp(xu,'nm'), pt = '%.0f'; else, pt = '%.3g'; end
sh = 1/(10*m); % move to correct positions
for it1 = 1:m
    pos = it1/m-1/m/2-0.03*is3D;
    text(ax_,pos,0.1,sprintf(pt,x(it1)),opts{:});
    ax(it1).Position(1:3) = [sh+(it1-1)*(1/m-sh/m) 0.15 1/m-2*sh/m];
end
if hz, xu = ['P' xu]; wl = '$\omega_0$ '; else, wl = '$\lambda_0$ '; end
if strcmpi(vars.interpreter,'latex')
    t = [wl '[$\mathrm{' xu '}$]'];
else
    wl = split(wl,{'$'}); t = [wl{2} ' [' xu ']'];
end
text(ax_,0.47,0.05,t,opts{:});

ax_.Visible = 'off'; ax_.NextPlot = 'add';
ax_.Position = [sh-0.02 0 1-2*(sh-0.02) 0.95];
ax_.XTick = 0.5:1:m; ax_.XLim = [0 m]; ax_.XTickLabels = string(x);

% legend
leg = cell(1,length(index));
for it1 = 1:length(index)
    idx = index(it1);
    if is3D
        s = mesh(ax_,zeros(2),zeros(2),zeros(2)); 
        s.CData = ones(2)*idx; s.FaceAlpha = 0.1;
        s.FaceColor = 'flat';
    else
        l = plot(ax_,NaN);
    end
    if idx<RP.nModes+1
        leg{it1} = sprintf('Mode %d',idx);
    elseif idx==RP.nModes+1
        leg{it1} = 'Background';
    elseif idx==RP.nModes+2
        leg{it1} = 'Sum';
        if is3D, s.FaceColor = 'none'; end
    elseif idx==RP.nModes+3
        leg{it1} = 'Reference';
        if is3D
            s.EdgeColor = 'none'; s.FaceColor = 'flat';
            s.FaceAlpha = 0.5; s.CData = s.CData+1;
        else
            l.LineWidth = 2; l.LineStyle = ':';
        end
    end
end
legend(ax_,leg,vars.opts{:},'location',vars.legendLocation);
