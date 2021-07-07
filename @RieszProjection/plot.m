function varargout = plot(RP, varargin)
%PLOT plot results
%   PLOT(RP) A list of quanities is displayed which are available for
%   plotting. Left click with the mouse will plot the quantity under the
%   cursor.
%
%   PLOT(RP,quantity) Plots the provided quantity, which must be contained
%   in the list that is displayed in the case above.
%
%   PLOT(RP,quantity,index) Plots the data corresponding to index, which
%   can be an integer or an array of integers. The indices correspond to
%   contours and are displayed on the upper axis of the plot showing the
%   integration contours. If there are n modal contributions, the index n+1
%   gives the contribution of the background, n+2 the sum of the plotted
%   contributions and n+3 the reference solution.
%
%   PLOT(...,PARAM1,VAL1,PARAM2,VAL2,...) Specifies
%   parameters with name value pairs which include:
%      fignumber (integer): Sets the number of the figure.
%      ylim (numeric): Sets the ylim-property of the axis. If [0 0]
%         (default) the limits will be determined automatically.
%      xlim (numeric): Sets the xlim-property of the axis. If [0 0]
%         (default) the limits will be determined automatically.
%      pointEvaluation (logical): If true (default) the quantity
%         NormalizedDecayRate is derived based on a point evaluation
%         if data is available otherwise it is obtained from the quantity
%         ElectromagneticFieldEnergyFlux.
%      interpreter (string): Defines the interpreter for all the text
%         displayed in the figures. The default is 'latex'.
%      fontSize (numeric): Defines the fontSize for all text displayed
%         in the figures. The default is 12.
%      legendLocation (string): Sets the location of the legend. The
%         default is 'northeast'.
%      xunit (string): When plotting an expansion the unit on the x-axis
%         can be changed from frequency (Hz, default) to a wavelength
%         possible inputs are 'nm', 'um', 'mm' and 'm'.
%      yunit (string): When plotting an expansion a unit can be added to
%         the axis label. E.g. J in the case of ElectricFieldEnergy. It
%         will appear automatically in square brackets.
%      contourNumbers (logical): If true (default) on the upper axis of the
%         plot showing the integration contours, the indices corresponding
%         to the contours are displayed.
%      quantities (cell): Select which quantities you want to include to
%         the error plot.
%      add2bg (int): Integer accounting for modes which are not of
%         interest and can be added to the background.
%      errortype (string): Defines the error type. Possible values are
%         'relative' (default), 'absolute' and 'minimum'. The last option
%         displays the minimum ob absolute and relative errror. 
%
%   fig = PLOT(...) Returns the figure handle.
%
%   The parameters you provide are saved such that you do not have to
%   provide them again if you replot the same quantity. Some are changed
%   for all quantities and others only for those they are relavant for. If
%   you do not provide a quantity and/or an index, these parameters will
%   be taken from the current figure. E.g. if you have plotted the
%   normalized decay rate with an instance 'rp' of the class 
%   'RieszProjection', you can change the unit of the x-axis from Hz to nm
%   by simply typing: rp.plot('xunit','nm').
%
%   see also figure, legend, ylim, xlim

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

% construct the parser if not yet there and parse the input arguments
persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'RieszProjection/plot';
    addOptional(parser,'quantity','?',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addOptional(parser,'index',0,...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'}));
    addParameter(parser,'fignumber',0,...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'}));
    addParameter(parser,'ylim',[0 0],...
        @(x)validateattributes(x,{'numeric'},...
        {'2d' 'numel' 2 'nondecreasing' 'real'}));
    addParameter(parser,'xlim',[0 0],...
        @(x)validateattributes(x,{'numeric'},...
        {'2d' 'numel' 2 'nondecreasing' 'real'}));
    addParameter(parser,'pointEvaluation',true,...
        @(x)validateattributes(x,{'logical'},{'scalar'})); 
    addParameter(parser,'interpreter','latex',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addParameter(parser,'fontSize',12,...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'})); 
    addParameter(parser,'legendLocation','northeast',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addParameter(parser,'xunit','Hz',...
        @(x)validateattributes(x,{'char'},{}));
    addParameter(parser,'yunit','',...
        @(x)validateattributes(x,{'char'},{}));
    addParameter(parser,'contourNumbers',true,...
        @(x)validateattributes(x,{'logical'},{'scalar'})); 
    addParameter(parser,'quantities',{[]},...
        @(x)validateattributes(x,{'cell'},{'vector'})); 
    addParameter(parser,'add2bg',[],...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'}));
    addParameter(parser,'errortype','relative',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
end

% check which quantities can be plotted
quantities = fieldnames(RP.expansions)';
q_ref = fieldnames(RP.reference)';
if isempty(q_ref), q_ref = {}; end
idx = ismember(q_ref,[quantities 'fields' 'frequencies']);
quantities = [q_ref(~idx) quantities];
normalizedQs = {'NormalizedDecayRate' 'PhotonCollectionEfficiency'};
idx = ismember({'DipoleEmission' 'Radiation'},quantities);
quantities = [normalizedQs(idx) quantities];
allQuantities = quantities;
if ~isempty(quantities) && ~isempty(fieldnames(RP.convergenceData))
    allQuantities = ['error' quantities];
end
allQuantities = ['contours' allQuantities];

% check which optional arguments are given
if nargin>1
    try quantity = validatestring(varargin{1},allQuantities); 
    catch, quantity = []; end
    try fn = RP.plotSettings.fignumber.(quantity); catch, fn = 0; end
    fig = get(groot,'CurrentFigure'); % get current figure if exists
    if isempty(quantity) % add to current figure
        try qidx = fig.UserData; catch, qidx = {'?' 0}; end
        if isnumeric(varargin{1})
            varargin = [qidx(1) varargin];
        else
            varargin = [qidx varargin];
        end
    elseif ~isempty(fig) && fn~=fig.Number
    elseif nargin>2 && ~isnumeric(varargin{2})
        try idx = fig.UserData{2}; catch, idx = 0; end
        varargin = [quantity idx varargin(2:end)];
    else
        varargin{1} = quantity;
    end
end

parse(parser,varargin{:});
vars = parser.Results; 
vars.quantityindex = find(strcmp(vars.quantity,quantities))+2;

% check wether the index is valid
if any(vars.index) && ~all(vars.index)
    error('All indices must be >0.')
elseif strcmp(vars.quantity,'contours')
    vars.maxIndex = length(RP.contours); % passed to plotContours
    if RP.includeConjugatedPoles % RP.nModes is probably not yet defined
        vars.maxIndex = ceil(vars.maxIndex/2);
    end
    if any(vars.index) && ~all(vars.index<=vars.maxIndex)
        error(['The index must be >0 and <=%d, which is '...
            'the index of the background contour'], vars.maxIndex)
    end
elseif strcmp(vars.quantity,'error') && any(vars.index>RP.nModes+2)
    error(['The index must be <=%d, which is the index '...
        'of the physical quantity'],RP.nModes+2)
elseif any(vars.index) && ~all(vars.index<=RP.nModes+3)
        error(['The index must be >0 and <=%d, which is '...
            'the index of the reference solution'], RP.nModes+3)
end
  
% load and save settings for the plot
usingDefaults = parser.UsingDefaults;
parameters = parser.Parameters;
pS = RP.plotSettings;
q = vars.quantity; 
locals = {'fignumber','yunit','ylim','index','xlim','legendLocation'};
ce = any(strcmp(q,{'contours' 'error'}));
if ~any(strcmp(usingDefaults,'xunit'))
    try pS.xlim=rmfield(pS.xlim,'expansion'); catch, end
end
for it = 1:length(parameters)
    parameter = parameters{it}; islocal = strcmp(parameter,locals);
    if ismember(parameter,usingDefaults)
        if ~isfield(pS,parameter) % check if contained in pS
            if strcmp(parameter,'legendLocation')
                if strcmp(q,'contours') % other default value for contours
                    vars.legendLocation = 'northwest';
                end
            end
            continue; 
        end
        if ~any(islocal)
            vars.(parameter) = pS.(parameter);
        elseif isfield(pS.(parameter),q)
            vars.(parameter) = pS.(parameter).(q);
        elseif isfield(pS.(parameter),'expansion') && ~ce
            vars.(parameter) = pS.(parameter).expansion;
        end
    else
        if strcmp(parameter,'quantity'), continue; end
        if ~any(islocal)
            pS.(parameter) = vars.(parameter);
        elseif any(islocal(1:3)) || (any(islocal(4:6)) && ce)
            pS.(parameter).(q) = vars.(parameter);
        elseif any(islocal(4:6))
            pS.(parameter).expansion = vars.(parameter);
        end
    end
end
RP.plotSettings = pS;
vars.opts = {'Interpreter' vars.interpreter 'FontSize' vars.fontSize};
vars.warning = ~any(strcmp(usingDefaults,'index'));
quantities = quantities(~strcmp(quantities,'RadiationPattern'));

% plot
switch lower(q)
    case 'contours'
        fig = RP.plotContours(vars);
    case lower(quantities)
        fig = RP.plotQuantity(vars);
    case 'radiationpattern'
        fig = RP.plotRadiationPattern(vars);
    case 'error'
        fig = RP.plotError(vars);
    otherwise
        href = '<a href="matlab: %s.plot(''%%s'')">%%s</a>\n';
        href = sprintf(href,inputname(1));
        str = ['You can plot the following '...
            'quantities:\n' repmat(href,1,length(allQuantities))];
        qs = cell(1,2*length(allQuantities));
        qs(1:2:end-1) = allQuantities; qs(2:2:end) = allQuantities;
        str = sprintf(str,qs{:});
        if ~strcmp(vars.quantity,'?')
            str = sprintf(['Unknown quantity %s. ' str], vars.quantity); 
        end
        str = str(1:end-1);
        disp(str);
        return;
end
if ishandle(fig), fig = {fig}; end
for it = 1:length(fig)
    fig{it}.UserData = {vars.quantity vars.index};
end
if ~isempty(fig), RP.plotSettings.fignumber.(q) = fig{1}.Number; end 
if nargout>0, varargout = fig; end
end
