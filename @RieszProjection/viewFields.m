function viewFields(RP, varargin)
%VIEWFIELDS Expansion and visualisation of the total field
%   VIEWFIELDS(RP) Visualizes the reference field at the
%   frequencies specified in the property 'referencePoints'. If it is a 
%   scalar the field will be opened in JCMview. Otherwise a cartesian 
%   export of the x-y plane at z=0 with 100x100 equidistantly distributed
%   points is performed and the field intensities are plotted in a row in 
%   ascending order with the smallest expansion frequency to the left.
%
%   VIEWFIELDS(RP,w0) Visualizes the reference field at a selection
%   of real expansion frequencies w0. 
%
%   VIEWFIELDS(RP,w0,index) Expands and visualizes the modal 
%   contributions corresponding to the contours of the specified index.
%   Zero refers to the reference solution. If index is an array the
%   visualizations of the fields will be displayed as a matrix where each
%   row corresponds to the mode with the same index. The last row gives the
%   reference solution. If there are n modal contributions, the index n+1
%   gives the contribution of the background, n+2 the sum of the displayed
%   contributions and n+3 the reference solution.
%
%   VIEWFIELDS(___,PARAM1,VAL1,PARAM2,VAL2,___) Specifies
%   parameters with name value pairs which include:
%      quantity (string): The specified quantity will be exported in case
%         of an cartesian export. Please refere to the JCMsuite
%         <a href="matlab: web(['https://www.docs.jcmwave.com/' ...
%         'JCMsuite/html/ParameterReference/a062a9a76a59a014fec3' ...
%         '1f494b7447a1.html'],'-browser')"
%         >parameter reference</a> for a list of quantities which can be
%         exported. The default is ElectricFieldStrength.
%      viewer (string): jcm or matlab, if jcm is chosen all fields are 
%         opened individually with JCMview. The default is jcm for a
%         scalar w0 and matlab otherwise.
%      keys (struct): Struct whose fields may be GridPoints* or 
%         NGridPoints* where the asterix replaces X, Y or Z. They
%         correspond to the input parameters of the same name used for a
%         <a href="matlab: web(['https://www.docs.jcmwave.com/' ...
%         'JCMsuite/html/ParameterReference/0354993cda1716b2af' ...
%         'b7a6af10ec11c7.html'],'-browser')"
%         >cartesian export</a> with JCMsuite. The asterix stands for X, 
%         Y or Z. The default is:
%         struct('NGridPointsX',100,'NGridPointsY',100,GridPointsZ,0).
%         This defines a regular grid with 100x100 points in the
%         x-y-plane at z=0.
%      selection (integer): Vector of integers being a subset of the
%         argument 'index'. The selected modes are summed up.
%      cmap (numeric): Colormap as returned by colormap(...). The colormap
%         is used for the visualization with matlab (see IMSHOW). The
%         default is a grayscale image.
%      component (string): Defines the component to be visualized. The
%         possibilities are real*, imag*, abs and log. The asterix
%         stands for x, y and z. Be aware that the intensities will not
%         sum up to the reference solution as taking the absolute value 
%         is not a holomorphic operation. The default is realx
%      normalize (logical): If true, the values will be normalized to
%         range between -1 and 1. The default is true if the component is
%         neither abs nor log.
%      displayRange (numeric): Two-element vector [LOW HIGH] that controls
%         the display range. The default is [-1 1] if normalize = true
%         and [min max] otherwise.
%      d (numeric): A scalar that gives the width of a colored frame which
%         helps to identify the rows.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

if isempty(RP.poles)
    error('Poles are required for an expansion of the total field');
elseif isempty(RP.fields) % RP.nModes must be known beforehand
    RP.evaluate;
end

% construct the parser if not yet there and parse the input arguments
persistent parser
if isempty(parser)
    red = [zeros(4,1);linspace(0,1,28).';ones(27,1);linspace(1,0.8,5).'];
    blue = [linspace(0.8,1,5).';ones(27,1);linspace(1,0,28).';zeros(4,1)];
    green = [zeros(4,1);linspace(0,1,28).';linspace(1,0,28).';zeros(4,1)];
    cm = [red green blue];
    parser = inputParser;
    parser.FunctionName = 'RieszProjection/plot';
    addOptional(parser,'w0',RP.referencePoints,...
        @(x)validateattributes(x,{'numeric'},...
        {'vector' 'nonempty' 'positive' 'increasing'}));
    addOptional(parser,'index',0,...
        @(x)validateattributes(x,{'numeric'},...
        {'integer' 'vector' 'nonnegative' 'increasing'}));
    addParameter(parser,'quantity','ElectricFieldStrength',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addParameter(parser,'viewer','?',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addParameter(parser,'keys',...
        struct('NGridPointsX',100,'NGridPointsY',100),...
        @(x)validateattributes(x,{'struct'},{'scalar'}));
    addParameter(parser,'cmap',cm,...
        @(x)validateattributes(x,{'numeric' 'char'},{})); 
    addParameter(parser,'component','realx',...
        @(x)validateattributes(x,{'char'},{'nonempty'}))
    addParameter(parser,'normalize',true,...
        @(x)validateattributes(x,{'logical'},{'scalar'})); 
    addParameter(parser,'displayRange',[0 0],...
        @(x)validateattributes(x,{'numeric'},...
        {'size', [1 2], 'increasing'}));
    addParameter(parser,'selection',[],...
        @(x)validateattributes(x,{'numeric'},...
        {'vector' 'integer' 'positive' 'increasing'}));
    addParameter(parser,'d',3,...
        @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    addParameter(parser,'interpreter','latex',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addParameter(parser,'fontSize',12,...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'})); 
    addParameter(parser,'xunit','Hz',...
        @(x)validateattributes(x,{'char'},{})); 
    addParameter(parser,'fignumber',40,...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'}));
end
parse(parser,varargin{:}); vars = parser.Results;
% check wether the index is valid
index = vars.index; index(~index) = RP.nModes+3;
if any(index) && ~all(index)
    error('All indices must be >0.')
elseif any(index) && ~all(index<=RP.nModes+3)
        error(['The index must be >0 and <=%d, which is '...
            'the index of the reference solution'], RP.nModes+3)
end
% check whether the selection is valid
if ~all(ismember(vars.selection,index)) ||...
        ~all(vars.selection<=RP.nModes+1)
    error(['The selection must contain a subset of the argument '...
        '''index'' and must be smaller than %d, which is the index '...
        'of the background contour.'], RP.nModes+1)
end
getFromPlotSettings = {'xunit','interpreter','fontSize'};
for option = getFromPlotSettings
    s = option{1};
    if ismember(s,parser.UsingDefaults) && isfield(RP.plotSettings,s)
        vars.(s) = RP.plotSettings.(s);
    end
end
if ismember('fignumber',parser.UsingDefaults) 
    try vars.fignumber = RP.plotSettings.fignumber.viewFields; catch; end
end
RP.plotSettings.fignumber.viewFields = vars.fignumber;
if isempty(vars.w0), fprintf('No real frequencies specified\n'); end
displayRange = vars.displayRange;
if strcmp(vars.viewer,'?') 
    if isscalar(vars.w0) && isscalar(index)
        vars.viewer = 'jcm';
    else
        vars.viewer = 'matlab';
    end
end
viewer = vars.viewer; vars.w0 = vars.w0(:);
w0hex = cellfun(@(x){['s' x]},cellstr(num2hex(vars.w0)));
vars.w0hex = w0hex;
if strcmp(viewer,'matlab')
    cartesianFields = cell(1,RP.nModes+3);
    if ~isequal(vars.keys,RP.exportkeys)
        for it1 = 1:length(RP.fieldExpansion)
            for it2 = 1:length(w0hex)
                RP.fieldExpansion{it1}.(w0hex{it2}){2} = [];
            end
        end
        RP.exportkeys = vars.keys;
    end
end
vars.w0 = round(vars.w0);

% collect data

fieldExpansions = getData(RP,vars); 

% visualize

if strcmp(vars.viewer,'jcm')
    index(index==RP.nModes+3) = RP.nModes+2;
    for it1 = 1:size(fieldExpansions(index,:),1)
        tag = fieldExpansions{1,it1};
        if isfield(RP.f.resultbags{2}.results_,tag)
            rb = RP.f.resultbags{2};
        else
            rb = RP.f.pprbs.superposition;
        end
        jcmwave_view(rb.results_.(tag).keys.oPath)
    end
    return
end
cs = {'realx' 'realy' 'realz' 'imagx' 'imagy' 'imagz' 'abs' 'log'};
c = validatestring(vars.component,cs); 
c = [int16(c(1)) int16(c(end))-119];
sum = 0; selection = 0;
for idx = 1:RP.nModes+3
    if ismember(RP.nModes+2,index) && idx==RP.nModes+2
        cartesianFields{idx} = sum; continue; 
    end
    fieldExpansion = squeeze(fieldExpansions(min(idx,end),:,:));
    if size(fieldExpansions,3)==1, fieldExpansion = fieldExpansion.'; end
    if strcmp(viewer,'jcm')
        for it2 = 1:size(fieldExpansion,2)
            jcmwave_view(fieldExpansion{1,it2})
        end
    elseif strcmp(viewer,'matlab')
        field = cellfun(@(x)permute(x.field{1},[2 1 3]), ...
            fieldExpansion(2,:),'UniformOutput',false);
        if strcmpi(vars.xunit,'hz')
            field = [field{:}];
        else
            field = [field{end:-1:1}];
        end
        if idx<RP.nModes+1, field = -field; end
        if c(1) == 'r'
            field = real(field(:,:,c(2)));
        elseif c(1) == 'i'
            field = imag(field(:,:,c(2)));
        elseif c(1) == 'l'
            field = log(vecnorm(field,2,3));
        else
            field = vecnorm(field,2,3);
        end
        if ismember(RP.nModes+2,index), sum = field + sum; end
        if ismember(idx,vars.selection)
            selection = field + selection;
            if idx==vars.selection(end)
                cartesianFields{idx} = selection;
            end
        elseif ismember(idx,index)
            cartesianFields{idx} = field;
        end
    end
end
fig = figure(vars.fignumber);
clf(fig); ax = axes(fig,'Visible','off');
field = cat(1,cartesianFields{end:-1:1});
if ~any(displayRange), displayRange = []; end
if vars.normalize && all(c(1)~=['a' 'l'])
    m = min(field,[],'all'); M = max(field,[],'all');
    mM = max(-m,M);
    field = field/mM;
    if isempty(displayRange), displayRange = [-1 1]; end
end
args = {'Colormap',vars.cmap,'DisplayRange',displayRange,'Parent',ax};
imshow(field, args{:}); ax.YDir = 'normal';
ax.NextPlot = 'add'; ax.DataAspectRatio = [1 1 1];
sel = ismember(index,vars.selection(1:end-1)); index(sel) = [];
nFields = length(index); nFrequencies = length(vars.w0);
height = size(field,1); width = size(field,2); d = vars.d;
x = {[0 width width 0] [1*d width-d width-1*d d]};
y = [[height height]-height/nFields [height height]];
y = {y [y(1:2)+d y(3:end)-d]};
name = cell(1,nFields); ndx = 1; line = cell(1,nFields);
args = {'FaceColor',vars.cmap(floor(end/2),:),'LineWidth',1.5,...
    'FaceAlpha',1}; c = ax.ColorOrder;
for it = index
    if ~vars.d, break; end
    if ~isempty(vars.selection) && it == vars.selection(end)
        name{ndx} = 'Selection';
    elseif it<RP.nModes+1
        name{ndx} = sprintf('Mode %d',it);
    elseif it==RP.nModes+1
        name{ndx} = 'Background';
    elseif it==RP.nModes+2
        name{ndx} = 'Sum';
    elseif it==RP.nModes+3
        name{ndx} = 'Reference Solution';
    end
    p = polyshape(x,y); y = cellfun(@(x){x-height/nFields},y);
    ps = plot(ax,p,'FaceColor',c(it,:),'EdgeColor','none');
    line{ndx} = plot(ax,polyshape,args{:}); % for legend
    line{ndx}.EdgeColor = ps.FaceColor; ps.FaceAlpha = 1;
    ndx = ndx+1;
end
ax.XLim = [0,width]; ax.YLim = [0,height];
args = {'FontSize',vars.fontSize,'Interpreter',vars.interpreter};
d = width/nFrequencies; xy = [d/2 0]; wl = '\lambda_0';
if strcmp(vars.interpreter,'latex'), wl = ['$' wl '$']; end
[x,xu] = RieszProjection.unit(vars.w0,vars.xunit); x = sort(x);
if strcmp(xu,'Hz'), xu = ['P' xu]; x = x*1e-15; wl = '$\omega_0$'; end
if strcmp(xu,'nm'), pt = '%.0f'; else, pt = '%.2g'; end
for it = 1:nFrequencies
    t = text(xy(1),xy(2),sprintf(pt,x(it)),args{:});
    t.Position(2) = t.Position(2)-t.Extent(end); xy(1) = xy(1)+d;
    t.Position(1) = t.Position(1)-t.Extent(3)/2;
end
t_label = text(width/2,-t.Extent(end),[wl ' [' xu ']'],args{:});
t_label.Position(2) = t_label.Position(2)-t_label.Extent(end);
t_label.Position(1) = t_label.Position(1)-t_label.Extent(3)/2;
if nFields>=nFrequencies
    args_ = {'NumColumns',1,'Location','WestOutside'};
else
    args_ = {'NumColumns',2,'Location','NorthOutside'};
end
legend([line{:}],name,args{:},args_{:});
if nFields>=nFrequencies, ax.OuterPosition = [0 0 0.9 0.9]; end
hold off
end

function fieldExpansions = getData(RP,vars)
index = 1:RP.nModes+2;
w0 = vars.w0; quantity = vars.quantity; w0hex = vars.w0hex;
ndef = ones([length(index),size(w0)],'logical');
fieldExpansions = cell(length(index),2,length(w0));
fE = cell(length(index),2);
for it1 = 1:length(index)
    w0 = vars.w0; w0hex = vars.w0hex;
    if length(RP.fieldExpansion)>=it1 && ~isempty(RP.fieldExpansion{it1})
        def = isfield(RP.fieldExpansion{it1},w0hex);
        fields = cellfun(@(x)RP.fieldExpansion{it1}.(x), ...
            w0hex(def),'UniformOutput',false);
        fieldExpansions(it1,:,def) = vertcat(fields{:})';
        w0(def) = []; w0hex(def) = []; ndef(it1,:) = ~def;
    end
    if any(ndef(it1,:))
        % prepare working directory
        pDir = [RP.f.projectDir filesep 'superpositions'];
        if ~exist(pDir,'dir'), mkdir(pDir); end
        if it1<=RP.nModes+1
            pDir = RP.f.projectDir;
            material = [pDir filesep 'materials.jcmt'];
            source = [pDir filesep 'sources.jcmt'];
            if it1==RP.nModes+1
                name = 'background'; 
            else
                name = sprintf('mode_%d',it1);
            end
            outNames = cellfun(@(x)[pDir filesep 'superpositions' ...
                filesep name filesep x],w0hex,'UniformOutput',false);
            for it2 = 1:length(w0)
                if ~exist(outNames{it2},'dir'), mkdir(outNames{it2}); end
                keys = RP.f.keys; keys.omega = w0(it2);
                if isfield(keys, 'relPerm') && ...
                        isa(keys.relPerm, 'function_handle')
                    keys.relPermittivity = keys.relPerm(keys.omega);
                end
                outfile = [outNames{it2} filesep 'materials.jcm'];
                if exist(material,'file')
                    jcmwave_jcmt2jcm(material,keys,'outputfile',outfile);
                else
                    copyfile(material(1:end-1),outfile);
                end
                if it1==RP.nModes+1
                    outfile = [outNames{it2} filesep 'sources.jcm'];
                    jcmwave_jcmt2jcm(source,keys,'outputfile',outfile);
                end
                outNames{it2} = [outNames{it2} filesep 'project_results'];
                if ~exist(outNames{it2},'dir'), mkdir(outNames{it2}); end
                outNames{it2} = [outNames{it2} filesep 'fieldbag.jcm'];
            end
            opts = {'index',it1,'outNames',outNames};
            fE(it1,1)= RP.integrate('Field',w0,opts{:}); 
            fE(it1,:) = {[fE{it1,1}{1,:}] fE{it1,1}(2,:)};
        elseif it1==RP.nModes+2
            tags = cell(1,length(w0));
            idx = zeros(size(w0));
            if isfield(RP.reference,'fields')
                frequencies = RP.reference.frequencies;
                [idx,loc] = ismembertol(w0,frequencies);
                if any(idx)
                    loc = loc(idx);
                    tags(idx) = RP.reference.fields(loc);
                end
            end
            if ~all(idx)
                fs = RP.f({w0(~idx)});
                tags(~idx) = fs{1}.fields;
            end
            fieldExpansions(it1,1,ndef(it1,:)) = tags;
        end   
    end
end

% collect superimposed fields and export cartesian fields if needed
for it1=1:length(index)
    if ~isempty(fE{it1,1})
        jcmwave_daemon_wait(fE{it1,1},RP.f.pprbs.superposition);
        fieldExpansions(it1,1,ndef(it1,:)) = fE{it1,2}; fE{it1,1} = [];
    end
    if strcmp(vars.viewer,'matlab')
        needed = cellfun(@isempty,fieldExpansions(it1,2,:));
        ws = vars.w0(needed); if isempty(ws), continue; end
        tags = squeeze(fieldExpansions(it1,1,needed));
        [fE{it1,:}] = export(RP.f,tags,vars.keys,quantity);
        ndef(it1,:) = needed;
    end
end
% collect the exports and save results
for it1 = 1:length(index)
    if ~isempty(fE{it1,1})
        rbe = RP.f.pprbs.cartesian;
        jcmwave_daemon_wait(fE{it1,1},rbe);
        fs = arrayfun(@(x)rbe.results_.(x{1}).result,fE{it1,2});
        fieldExpansions(it1,2,ndef(it1,:)) = fs;
    end
    positions = find(ndef(it1,:));
    if isempty(positions), continue; end
    for it2 = 1:length(vars.w0(ndef(it1,:)))
        idx = index(it1); if idx==RP.nModes+3, idx = idx-1; end
        pos = positions(it2); w0hex = vars.w0hex{it2};
        RP.fieldExpansion{idx}.(w0hex) = fieldExpansions(it1,:,pos);
    end
end
end
