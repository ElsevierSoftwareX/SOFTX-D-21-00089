function getContours(RP, varargin) 
%GETCONTOURS Get the contours for integration
%   GETCONTOURS(RP) If this method has been called before with additional
%   arguments, these are used for all subsequent calls. Otherwise elliptic
%   contours are constructed with equally distributed points.
%
%   GETCONTOURS(RP,shape) The parameter shape defines the background
%   contour. Possible values are 'ellipse' ('e'), 'circle' ('c') and 
%   'polygon' ('p'). The default is 'ellipse'.
%
%   GETCONTOURS(RP,shape,quadgk) The parameter quadgk determines if a
%   Gauss-Kronrod quadrature rule is used instead of the trapezoidal rule
%   for the integration of the background contour. 
%   If true, the integration path is devided in subsections which are
%   integrated with a 7/15-Point rule and can be refined adaptively.
%   Alternatively it can be set to 'local' ('l') or 'global' ('g'). The
%   former is equivalent to true. In the latter case a high order Gauss 
%   quadrature rule is selected based on the given number of integration
%   points. If the desired precision compared to the direct scattering
%   solutions at the reference points is not reached, inserting the Kronrod
%   nodes allows for a global refinement step. Further refinements will be
%   based on the local errors of the subintervals which than are available.
%
%   If the first Element of the property 'selectedPoles' is zero the
%   expansion points will be enclosed with a contour. If it is empty, you
%   will be asked to select poles unless the poles are not yet given.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

if any(contains(string(varargin),'updatenobg'))
elseif ~RP.newBackground
    answer = questdlg('Do you want to alter the outer contour?');
    if strcmp(answer,'Yes')
        RP.newBackground = true;
        RP.getContours(varargin{:})
    end
    return
end
if nargin<3
    if isfield(RP.contourSettings,'quadgk') % reuse previous parameters
        varargin{2} = RP.contourSettings.quadgk;
        if nargin==1 || contains(varargin{1},'update') % update contours
            varargin{1} = RP.contourSettings.currentShape;
        elseif nargin==2 && strcmp(varargin{1},'gk') % gk update
            RP.contourSettings.quadgk = 'l'; gl = false; gk = true;
            [RP.contours{end},RP.weights{end}] = distribute; return;
        end
    elseif nargin==2 && contains(varargin{1},'update')
        return
    end      
end

persistent parser;
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'RieszProjection/getContours';
    addOptional(parser,'shape','ellipse',@validateshape);
    addOptional(parser,'quadgk',false,@(x)validatequadgk(x));
end

parse(parser,varargin{:});
sp = parser.Results.shape(1); % shape of the contour (circle or polygon)
gk = parser.Results.quadgk(1); % use Gauss-Kronrod nodes and weights
gl = ischar(gk) && gk=='g'; % use only Gauss nodes and weights
% the trapezoidal rule is not supported along the sides of a polygon
if sp == 'p' && ~gk
    gk = true; gl = false; 
    fprintf(['Using Gauss-Kronrod quadrature rule, as the trapezoidal '...
        'rule is not implemented for polygons\n'])
end

% remember optional parameters
RP.contourSettings.currentShape = sp;
RP.contourSettings.quadgk = gk;
gk = logical(gk);

% make sure that there are selected poles
if isempty(RP.poles)||(~isempty(RP.selectedPoles)&&~RP.selectedPoles{1}(1))
    % if no poles are available, expansion points define the contour
    if isempty(RP.expansionPoints)
        fprintf('No expansion points defined.\n'); return; 
    end
    points = RP.expansionPoints([1 end]).'; minD = RP.minD;
    if sp == 'p'
        shift = diff(minD); 
        if shift==0, shift = minD(1); minD(2) = 0; end
        points = points+shift*[-1i;1i]; 
    end
    RP.shape = getShape(points(:),sp,minD);
    if gk, RP.shape = splitPath; end % intervals for quadgk
    [newContours,newWeights] = distribute;
    RP.contours = {newContours}; RP.weights = {newWeights};
    return
elseif isempty(RP.selectedPoles)
    % if poles are known, a selection is used to define contours
    radius = max(min(diff(RP.poles))/2,real(RP.poles(1))*1e-3);
    href = '<a href="matlab: %s.selectPoles(%e)">%s</a>';
    href = sprintf(href,inputname(1),radius,'interactively');
    fprintf('You must select Poles. This can be done %s.\n',href);
    return
end

radius = [num2cell(RP.radius_) RP.minD];
nPoints = [RP.nPoints_ RP.nPointsB];

% extract the poles
selectedPoles = cellfun(@(x){RP.poles(x)}, RP.selectedPoles);
% the background must enclose all selected poles and expansion points
background = RP.defineBg;
if isempty(RP.defineBg), background = [selectedPoles{:}]; end
if RP.includeConjugatedPoles
    background = [background conj(background(imag(background)~=0))];
end
if isscalar(RP.expansionPoints)
    background = [background RP.expansionPoints];
elseif ~isempty(RP.expansionPoints)
    background = [background RP.expansionPoints([1 end]).'];
end
% if the eigenvalues are used a single contour suffices
if RP.useEigenvalues
    selectedPoles = {background}; nContours = 1; 
    nPoints = RP.nPointsB; radius = {RP.minD};
else
    selectedPoles = [selectedPoles background]; 
    nContours = length(selectedPoles);
end
RP.weights = cell(1,nContours);
newContours = cell(1,nContours); % set dependent property last

for it = 1:nContours
    bg = it==nContours;
    if bg && ~RP.newBackground % if no new outer contour is required,
        shape = RP.shape; sh_ = shape{2}; % stick to the previous shape
    else
        % use always circles for modal contributions of a single pole
        if bg, sh_ = sp; else, sh_ = 'c'; end % allow for other shapes in future
        if ~bg && length(selectedPoles{it})==1, sh_ = 'c'; end
        shape = getShape(selectedPoles{it},sh_,radius{it});
    end
    % for quadgk the integration path must be split in intervals
    if gk && bg || sh_=='p', shape = splitPath(shape,nPoints(it)); end
    % distribute integration points along the path
    [newContours{it},RP.weights{it}] = distribute(nPoints(it),shape);
end
RP.shape = shape;
% mirror the modal contours at the real axis
if RP.includeConjugatedPoles && ~RP.useEigenvalues
    RP.weights = horzcat(repmat(RP.weights(1:end-1),1,2),RP.weights(end));
    shift = @(x) x-2i*imag(x(1));
    nC_conj = cellfun(shift,newContours(1:end-1),'UniformOutput',false);
    newContours = [newContours(1:end-1) nC_conj newContours(end)];
end
RP.contours = newContours;

    function shape = getShape(points,sp,r)
        points = reshape(points,1,[]);
        if length(points)>3 % reduce number of points
            ch = convhull(real(points),imag(points),'simplify',true);
            points = points(ch(1:end-1));
        end
        if sp == 'c' % circle defined by center an radius
            [c, R] = RieszProjection.minimalCircle(points);
            shape = {[c R+r(1)] sp 1};
        elseif sp == 'p'
            points = rmpointstol(points,1e-1*min(RP.minD));
            if length(points) == 2
                d = diff(r); if isempty(d)||d==0, d = r(1); r(2)=0; end
                points = points+d*[-1i;1i]; points = points(:).';
            end
            phi = mod(eps(10)+angle(points-mean(points)),2*pi);
            [~,ndx] = sort(phi);
            points = points([ndx(end) ndx]);
            s = diff(points); % sides
            n = s./abs(s)*exp(-1i*pi/2); % normals
            a = abs(angle([n(2:end) n(1)].*exp(-1i*angle(n)))); % angles
            v = points(2:end); % vertices
            a_ = angle(n)+a/2;
            v = v + 1./cos(a/2).*(r(1)*cos(a_)+1i*r(2)*sin(a_));
            shape = {[v v(1)] sp 1 [v v(1)]};
        elseif sp == 'e'
            E = RieszProjection.minimalEllipse(points);
            if isscalar(r)
                E(3) = E(3)+r(1);
            else
                E(3) = E(3) + sum(r)/2; E(4) = E(4)-diff(r)/2;
            end
            shape = {[0 1] sp 1 E};
        end
    end

    function shape = splitPath(shape,nP)
        if nargin == 0
            shape = RP.shape; nP = RP.nPointsB;
        end
        ngk = 15; % default for Gauss-Kronrod quadrature
        if shape{2} == 'c' || shape{2} == 'e'
            % choose the order of the Gauss quadrature if needed
            if gl, ngk = RP.gn(min(abs(RP.gn-nP))==abs(RP.gn-nP)); end
            h = 2*pi/round(nP/ngk); v = 0:h:2*pi; 
            shape{1} = [shape{1} v];
        elseif shape{2} == 'p'
            v = shape{4};
            ivslen = abs(diff(v)); pathlen = sum(ivslen);
            if gl
                m = min(ivslen)/pathlen*nP;
                ngk = RP.gn(min(abs(RP.gn-m))==abs(RP.gn-m));
            end
            nSplits = round(nP/ngk)-1;
            minlen = pathlen/nSplits;
            nNew = round(ivslen/minlen)-1; 
            idxNew = find(nNew>0);
            for I = idxNew(end:-1:1)
                newPts = v(I) + (1:nNew(I))./(nNew(I)+1)*(v(I+1)-v(I));
                v = [v(1:I) newPts v(I+1:end)];
            end
            if RP.includeConjugatedPoles && imag(v(end-1))>0
                phi = mod(eps(10)+angle(v(2:end)-mean(v(2:end))),2*pi);
                [~,ndx] = sort(phi); v = v(ndx([1:end 1])+1);
            end
            shape{1} = v;
        end
        shape{3} = ngk+(ngk+1)*gl; % here ngk is always given for gk nodes
    end

    function [contour,weights] = distribute(nP,shape)
        if nargin==0
            shape = RP.shape; nP = RP.nPointsB; bg = true;
        end
        [v,s,ngk] = shape{1:3};
        gk_ = (gk && bg) || s=='p'; % get gk weights
        if s == 'c' || s == 'e' % circular shape
            c = v(1); r = v(2); 
            if RP.includeConjugatedPoles && bg, c = real(c); end
            if gk_
                v = v(3:end);
            else 
                v = [0 2*pi]; weights = 2/nP; 
                nodes = (-1:2/nP:1-2/nP).';
            end
        end
        if gk_ % get Gauss-Kronrod nodes and weights
            nodes = RP.(sprintf('nodes%d',ngk));
            weights = RP.(sprintf('weights%d',ngk));
            w_err = RP.(sprintf('weights%dg',floor(ngk/2)));
            if gl
                nodes = nodes(2:2:end); 
                weights = w_err(1:2:end); 
            else
                w_err = weights-w_err;
                w_err = [w_err(end:-1:2) w_err].';
            end
            nodes = [-nodes(end:-1:1);0;nodes];
            weights = [weights(end:-1:2) weights].';
        end
        ivs = [v(1:end-1);v(2:end)];
        midpts = sum(ivs)/2; h = diff(ivs);
        contour = midpts + nodes*h/2;
        if s == 'c'
            alpha = contour;
            contour = c + r*exp(1i*alpha);
            weights = 1i*r*exp(1i*alpha).*weights;
            if gk_ && ~gl, w_err = 1i*exp(1i*contour/r).*w_err; end
        elseif s == 'e'
            alpha = contour; E = shape{4}; c = E(1)+1i*E(2); 
            if RP.includeConjugatedPoles, alpha = alpha-E(5); end
            A = E(3); B = E(4); a = alpha+E(5); b = alpha-E(5);
            contour = c + A*exp(1i*a)+B*exp(-1i*b);
            weights = 1i*(A*exp(1i*a)-B*exp(-1i*b)).*weights;
            if gk_&&~gl, w_err = 1i*(A*exp(1i*a)-B*exp(-1i*b)).*w_err; end
        end
        weights = h/2.*weights;
        if gk_ && ~gl, weights = {weights h/2.*w_err}; end
        
    end
end

function points = rmpointstol(points,tol)
if length(points)<3, return; end
s1 = [points(2:end) points(1)]-points;
s2 = [points(3:end) points(1:2)]-points;
s1 = s1.*exp(-1i*angle(s2));
d = abs(sin(angle(s1)).*abs(s1));
selection = circshift(d>tol,1);
if sum(selection)<2
    re = real(points); im = imag(points);
    [~,ire] = sort(re); [~,iim] = sort(im);
    if abs(diff(re(ire([1 end]))))>abs(diff(im(iim([1 end]))))
        points = points(ire([end 1]));
    else
        points = points(iim([end 1]));
    end
else
    points = points(circshift(d>tol,1));
end
end

function validateshape(x)
if ischar(x) && any(startsWith(x,{'c','p','e'})), return; end
error(['The value of this parameter (default: ''cirlce'') ' ...
    'can be set to ''circle'', ''ellipse'' or ''polygon''.'])
end

function validatequadgk(x)
if islogical(x) && isscalar(x), return; end
if ischar(x) && any(strcmp(x,{'l','g'})), return; end
error(['The value of this parameter (default: false) '...
    'can be a logical scalar or ''l'' or ''g'', '...
    'where ''l'' is equivalent to true.'])
end
