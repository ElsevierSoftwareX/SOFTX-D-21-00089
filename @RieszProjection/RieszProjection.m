classdef RieszProjection < handle
    %RIESZPROJECTION Class aiming to organize Riesz projection projects
    %   On the one hand this class is designed to allow for an easy 
    %   construction of well suited contours in the complex plain and for
    %   the application of different quadrature rules with adaptive
    %   refinement. On the other hand it provides together with the class
    %   '<a href="matlab: help Scattering"
    %   >Scattering</a>' used to solve scattering problems and apply post
    %   processes an interface to the FEM solver JCMsuite.
    %   The aim is to make the Riesz projection method for modal expansions
    %   easily available and thereby make use of symmetries and
    %   parallelization for the sake of efficiency.
    %   If the number of relevant eigenvalues is small (i.e., <= 5) the
    %   evaluation of the eigenvalues and the expansion can be done with a
    %   single contour enclosing the relevant region of the complex plane.
    %   If this is not the case, you currently have to solve the eigenvalue
    %   problem externally and provide the complex eigenfrequencies as a
    %   vector. This can be done with the <a href="matlab: web(['https' ...
    %   '://jeos.springeropen.com/articles/10.1186/s41476-019-0098-z'], ...
    %   '-browser')">solver of JCMsuite</a>. Please refer
    %   to its <a href="matlab: web(['https://www.docs.jcmwave.com/' ...
    %   'JCMsuite/html/ParameterReference/9184954f6f210148ef1d27247' ...
    %   'b8d7d47.html'],'-browser')">parameter reference</a> for details.
    %   The idea to use an approach based on Riesz projections (RPs) has
    %   been presented 2018 by Zschiedrich et al. <a
    %   href="matlab: web(['https://journals.aps.org/pra/abstract/' ...
    %   '10.1103/PhysRevA.98.043806'],'-browser')">PhysRevA</a>. 
    %   For further details on the expansion of quadratic functionals using
    %   RPs including such related to the far field of the scatterer,  
    %   please refer to an article by Binkowski et al. (2020) <a
    %   href="matlab: web(['https://journals.aps.org/prb/abstract' ...
    %   '/10.1103/PhysRevB.102.035432'],'-browser')">PhysRevB</a>. 
    %   The implemented eigensolver is based on an algorithm published 2020
    %   by Binkowski et al. <a
    %   href="matlab: web(['https://www.sciencedirect.com/science' ...
    %   '/article/pii/S0021999120304526?via%3Dihub'],'-browser')"
    %   >jcp</a>. It is used by default if no eigenvalues
    %   are given. Please refer to the properties <a href="matlab: 
    %   help('RieszProjection.poles')">poles</a> and <a href="matlab: 
    %   help('RieszProjection.nEigenvaluesExpected')"
    %   >nEigenvaluesExpected</a>.
    %   
    %RIESZPROJECTION Properties:
    %   poles - Poles of f (see below) if they are known 
    %   selectedPoles - Grouped poles defining the modal contributions
    %   expansionPoints - Expansion points on the real axis
    %   nPoints - Number of points of the contours except the background
    %   nPointsB - Number of points of the background contour
    %   radius - Radius of circles enclosing the selected poles
    %   minD - Minimal distance of the background to the enclosed poles
    %   contours - Complex frequencies defining the contours
    %   precision - Precision demanded for adaptive refinement
    %   referencePoints - Real frequencies of reference solutions
    %   includeConjugatedPoles - Include the complex conjugate poles
    %   useEigenvalues - Expansion based on a single contour
    %   nEigenvaluesExpected - Estimated number of eigenvalues
    %   f - Expensive part of the function to be integrated
    %   weights - Integration weights
    %   fields - Cell array as returned by f
    %   derivedQuantities - Values at integration points, e.g., energies
    %   expansions - Struct holding data from expansions
    %   reference - Struct containing reference solutions
    %   fieldExpansion - Cell array holding data for visualization
    %
    %RIESZPROJECTION Methods:
    %   RIESZPROJECTION - Construct an instance with given properties
    %   selectPoles - Select poles interactively
    %   getContours - Get the contours for integration
    %   expand - Expand a selected quantity
    %   getReferenceSolution - Evaluate reference solutions at real omega
    %   plot - Plot results and errors
    %   viewFields - Expand the total field and visualize
    %   getError - Get the error of the integrals
    %
    %RIESZPROJECTION static Methods:
    %   minimalCircle - compute minimal area enclosing circle
    %   minimalEllipse - compute minimal area enclosing ellipse
    %
    %DEPENDENCIES:
    %   Though the construction of the contours and the integration are
    %   largely independent of JCMsuite, it requires some changes of the
    %   code to use it independently. For a use with 'Scattering' as an
    %   interface to JCMsuite, the third party support has to be set up.
    %   Please refer to the documentation of the <a href="matlab:
    %   web(['https://docs.jcmwave.'com/JCMsuite/html/MatlabInterface' ...
    %   '/index.html'],'-browser')">matlab interface</a> to set up
    %   the hird party support and to the homepage of <a href="matlab:
    %   web('https://jcmwave.com/jcmsuite/jcmsuite-documentation', ...
    %   '-browser')">JCMwave</a> for
    %   installation instructions.
    %
    %   see also Scattering

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021
      
    properties (Dependent)
        % poles - Poles of the system if known. Especially if more than
        % five eigenvalues are located in the outer contour, you currently
        % need to solve the eigenvalue problem externally and provide the
        % complex eigenfrequencies as an array, which will be sorted 
        % automatically when the property is set. If no poles are given,
        % the eigenvalue problem is solved using a contour integral
        % algorithm. This requires an estimated number of eigenvalues
        % expected to be found inside the contour. This can be set in the
        % property <a href="matlab: help( ...
        % 'RieszProjection/nEigenvaluesExpected')">nEigenvaluesExpected</a>.
        poles (1,:) double;
        
        % selectedPoles - Grouped poles defining the modes to project on.
        % This cell array contains vectors of integers which are indices to
        % the property 'poles'. If you do not want to  select any, set this 
        % property to {0}.
        selectedPoles (1,:) cell;
        
        % expansionsPoints - Expansion points on the real axis.
        expansionPoints (:,1) double {mustBeReal};
        
        % nPoints - Number of points on the contours around the selected 
        % poles. It can be a scalar or an array of the same length as the
        % property 'selectedPoles'.
        nPoints (1,:) double ...
            {mustBePositive,mustBeInteger};
       
        % Number of points defining the background contour.
        nPointsB (1,1) double ...
            {mustBePositive, mustBeInteger};
        
        % radius - Radius of the circles enclosing the selected poles. If
        % more than one pole must be enclosed it defines the minimal
        % distance of the contour to any of the poles
        radius (1,:) double {mustBeNonnegative};
        
        % minD - Approximate minimal distance of the background contour to 
        % enclosed poles. If the contour is an ellipse it can be a
        % two-element vector. The first value is than added to the
        % semi-major axis and the second to the semi-minor axis. In the
        % case of a polygon the first entry results in a shift in x
        % direction and the second entry in a shift in y direction.
        minD (1,:) double;
        
        % contours - Contours associated with selected poles. The
        % background contour RP.contours{end} is the contour enclosing all
        % poles of interest and all expansion points if given.
        contours (1,:) cell;
        
        % precision - This property can be a scalar or a vector with two
        % elements. In the latter case the first value is used as upper
        % bound of the absolute error and the second one as upper bound of
        % the relative error. In the former case both bounds have the same
        % value. This property is used for an adaptive refinement of the
        % integration contours. For absolute and relative error control
        % you can set the first or the second value to zero respectively.
        % If the precision is Inf (default), no refinements are made. 
        precision (1,2) double;
        
        % useEigenvalues - If this property is set true only the background
        % contour is used to obtain the expansion. Setting this property
        % will cause the creation or deletion of the other contours.
        useEigenvalues (1,1) logical;
        
        % includeConjugatedPoles - If set true, the poles of f(-w) are
        % taken into account which is necessary for the expansion of
        % quadratic forms such as the electric field energy.
        includeConjugatedPoles (1,1) logical;
        
        % nEigenvaluesExpected - If the eigenvalues are not given, a guess
        % of the expected number of eigenvalues is needed. It should better
        % be too large, than too small. Eigenvalues whose normalized
        % residues are smaller than max(1e-2,<a href="matlab: help(...
        % 'RieszProjection/precision')">precision</a>(end)) are ignored. If
        % this property is set even though poles are already given, these
        % will be reset to an empty list. 
        nEigenvaluesExpected (1,1) double {mustBePositive,mustBeInteger};
    end
    
    properties
        % referencePoints - If not empty the scattering problems at these
        % points are solved and used to estimate the error of the
        % expansion.
        referencePoints (:,1) double {mustBeReal};
        
        % defineBg - If not empty, this points are used instead of the
        % selected poles in Order to define the background contour
        defineBg (:,1) double;
    end
    
    properties (SetAccess=private)
        % f - Expensive part of the function to be integrated given either
        % as function handle or callable class which takes a cell array
        % of doubles with shape=(1,:) containing complex doubles.
        % Working with JCMsuite you should use an instance of the class
        % Scattering which returns a cell array containing cells of
        % character arrays.
        % f({}) must return 0 to pass the simple test mustTakeCell.
        f {mustTakeCell};
        
        % weights - Defined by the contours and needed for the integration
        % with the trapezoidal rule. For a contour gamma those are the
        % derivative of gamma at the corresponding point times h being the 
        % distance to the next point.
        weights (1,:) cell;
        
        % fields - Cell array as returned by f
        fields (1,:) cell;
        
        % derivedQuantities - Quantities derived from the fields e.g. 
        % energies.
        derivedQuantities (1,1) struct;
        
        % expansions - The struct expansions holds data collected with the 
        % method expand.
        expansions (1,1) struct;
        
        % fieldExpansion - The cell array fieldExpansion contains structs 
        % whose fields are hexadecimal representations of correspoinding 
        % expansion points with a leading s to make it evaluate to a valid
        % field name. The index of the cell array corresponds to a contour 
        % index.
        fieldExpansion (1,:) cell;
        
        % reference - Struct containing reference solutions
        reference (1,1) struct;
    end
    
    properties (Hidden)
        % g - Cheap part of the function to be integrated, if expansion
        % points are defined g must take two arguments: first frequency and 
        % second expansion points. For Riesz projections the default value
        % does not have to be changed. If you want to use the contour
        % integrations for an eigensolver you might wish to evaluate
        % integrals for various of these functions.
        g function_handle = @(w, w0) 1/(2i*pi).*1./(w-w0);
        
        % dipoleBulkEmission - Used to normalize the decay rate of a dipole
        % source placed in a resonator. Working with JCMsuite and the class
        % Scattering this property is set automatically. 
        dipoleBulkEmission function_handle;
        
        % quadraticQ - Quantities which are available for expansion and
        % quadratic in the electric Field.
        quadraticQ (1,:) cell = {};
        
        % linearQ - Quantities which are linear in the electric field are
        % given here.
        linearQ (1,:) cell = {'DipoleEmission' 'Field'};
        
        % newBackground - If this value is set to false, the
        % parameterization of the background contour will not be changed
        % when the method 'getContours' is called. 
        % This can be used if small changes do not require an 
        % adaption of the integration path, e.g., the recalculation of the
        % eigenvalues will automatically set this property to false.
        newBackground (1,1) logical = true;
    end
    
    properties (Constant, Hidden)
        c0 = 299792458; % speed of light
        gn = [7 15 25];
        % list of properties
        % nodes and weights for common Gauss-Kronrod quadrature rules
        % https://www.advanpix.com/2011/11/07/
        % gauss-kronrod-quadrature-nodes-weights/
        nodes15 = [ ...
            0.2077849550078985; 0.4058451513773972; ...
            0.5860872354676911; 0.7415311855993944; ...
            0.8648644233597691; 0.9491079123427585; ...
            0.9914553711208126];
        weights15 = [ ...
            0.2094821410847278, 0.2044329400752989, ...
            0.1903505780647854, 0.1690047266392679, ...
            0.1406532597155259, 0.1047900103222502, ...
            0.06309209262997855, 0.02293532201052922];
        weights7g = [ ...
            0.4179591836734694, 0, 0.3818300505051189, 0, ...
            0.2797053914892767, 0, 0.1294849661688697, 0];
        nodes31 = [ ...
            0.1011420669187175; 0.2011940939974345; ...
            0.2991800071531688; 0.3941513470775634; ...
            0.4850818636402397; 0.5709721726085388; ...
            0.6509967412974170; 0.7244177313601700; ...
            0.7904185014424659; 0.8482065834104272; ...
            0.8972645323440819; 0.9372733924007059; ...
            0.9677390756791391; 0.9879925180204854; ...
            0.9980022986933971];
        weights31 = [ ...
            0.1013300070147915, 0.1076984552387560, ...
            0.09917359872179196, 0.09664272698362368, ...
            0.09312659817082532, 0.08856444305621177, ...
            0.08308050282313302, 0.07684968075772038, ...
            0.06985412131872826, 0.06200956780067064, ...
            0.05348152469092809, 0.04458975132476488, ...
            0.03534636079137585, 0.02546084732671532, ...
            0.01500794732931612, 0.005377479872923349];
        weights15g = [ ...
            0.2025782419255613, 0, 0.1984314853271116, 0, ...
            0.1861610000155622, 0, 0.1662692058169939, 0, ...
            0.1395706779261543, 0, 0.1071592204671719, 0, ...
            0.07036604748810812, 0, 0.03075324199611727, 0];
        nodes51 = [ ...
            0.06154448300568508; 0.1228646926107104; ...
            0.1837189394210489; 0.2438668837209884; ...
            0.3030895389311078; 0.3611723058093878; ...
            0.4178853821930377; 0.4730027314457150; ...
            0.5263252843347192; 0.5776629302412230; ...
            0.6268100990103174; 0.6735663684734684; ...
            0.7177664068130844; 0.7592592630373576; ...
            0.7978737979985001; 0.8334426287608340; ...
            0.8658470652932756; 0.8949919978782754; ...
            0.9207471152817016; 0.9429745712289743; ...
            0.9616149864258425; 0.9766639214595175; ...
            0.9880357945340772; 0.9955569697904981; ...
            0.9992621049926098];
        weights51 = [ ...
            0.06158081806783294, 0.06147118987142532, ...
            0.06112850971705305, 0.06053945537604586, ...
            0.05972034032417406, 0.05868968002239421, ...
            0.05743711636156783, 0.05595081122041232, ...
            0.05425112988854549, 0.05236288580640748, ...
            0.05027767908071567, 0.04798253713883671, ...
            0.04550291304992179, 0.04287284502017005, ...
            0.04008382550403238, 0.03711627148341554, ...
            0.03400213027432934, 0.03079230016738749, ...
            0.02747531758785174, 0.02400994560695321, ...
            0.02043537114588284, 0.01684781770912830, ...
            0.01323622919557167, 0.009473973386174152, ...
            0.005561932135356714, 0.001987383892330316];
        weights25g = [ ...
            0.1231760537267155, 0, 0.1222424429903100, 0, ...
            0.1194557635357848, 0, 0.1148582591457116, 0, ...
            0.1085196244742637, 0, 0.1005359490670506, 0, ...
            0.09102826198296365, 0, 0.08014070033500102, 0, ...
            0.06803833381235692, 0, 0.05490469597583519, 0, ...
            0.04093915670130631, 0, 0.02635498661503214, 0, ...
            0.01139379850102629, 0];
    end
    
    % Get and set methods for the dependent properties.
    % If a dependent property is set to a new value, contours will be
    % reevaluated and function values which need to be
    % recalculated are marked in the cell array undefinedValues.
    methods
        function set.nPoints(RP, nPoints)
            if isequal(nPoints, RP.nPoints), return; end
            if ~isempty(RP.selectedPoles)
                if isscalar(nPoints)
                    nPoints = repmat(nPoints,1,length(RP.selectedPoles));
                elseif length(nPoints)~=length(RP.selectedPoles)
                    error(['To set the property ''nPoints'' provide a '...
                        'scalar or an array which matches the length '...
                        'of the property ''selectedPoles'' (%d)'], ...
                       	length(RP.selectedPoles));
                end
            end
            RP.nPoints_ = nPoints; 
            if ~RP.useEigenvalues, RP.getContours('updatenobg'); end
        end
        function value = get.nPoints(RP)
            if all(RP.nPoints_==RP.nPoints_(1))
                value = RP.nPoints_(1);
            else
                value = RP.nPoints_;
            end
        end
        
        function set.radius(RP, radius)
            if isequal(radius,RP.radius), return; end
            if ~isempty(RP.selectedPoles)
                if isscalar(radius)
                    radius = repmat(radius,1,length(RP.selectedPoles));
                elseif length(radius)~=length(RP.selectedPoles)
                    error(['To set the property ''radius'' provide a '...
                        'scalar or an array which matches the length '...
                        'of the property ''selectedPoles'' (%d)'], ...
                        length(RP.selectedPoles));
                end
            end
            RP.radius_ = radius; 
            if ~RP.useEigenvalues, RP.getContours('updatenobg'); end
        end
        function value = get.radius(RP)
            if all(RP.radius_==RP.radius_(1))
                value = RP.radius_(1);
            else
                value = RP.radius_;
            end
        end
        
        function set.minD(RP, minD)
            if all(minD==RP.minD), return; end
            if isscalar(minD), minD = [minD minD]; end
            RP.minD_ = minD;
            RP.getContours('update');
        end
        function value = get.minD(RP)
            value = RP.minD_;
        end
        
        function set.nPointsB(RP, nPointsB)
            if ~isempty(RP.contours)
                nP = numel(RP.contours{end});
            else
                nP = RP.nPointsB;
            end
            if nPointsB==nP, return; end
            if nPointsB>RP.nPointsB_ && ~RP.up2date(4)
                RP.up2date(3) = false; 
            end
            RP.nPointsB_ = nPointsB;
            % only the number of nodes will change but not the shape
            RP.getContours('updatenobg'); 
        end
        function value = get.nPointsB(RP)
            value = RP.nPointsB_;
        end
        
        function set.expansionPoints(RP, expansionPoints)
            if isequal(expansionPoints,RP.expansionPoints_), return; end
            expansionPoints = sort(expansionPoints);
            try ps1 = RP.expansionPoints_([1 end]); catch, ps1=[]; end
            try ps2 = expansionPoints([1 end]); catch, ps2 = []; end
            tolerance = -RP.minD + min(1e-2,RP.precision(1))*RP.minD;
            noBg = ~RP.newBackground && all(RP.inside(ps2,tolerance));
            RP.expansionPoints_ = expansionPoints(:);
            if noBg || isempty(RP.contours) || isequal(ps1,ps2)
                try fignumber = RP.plotSettings.fignumber.contours;
                    if ishandle(fignumber), RP.plot('contours'); end
                catch
                end
            else
                RP.getContours('update');
                inC = find(RP.inside(RP.poles));
                sP = inC(~ismember(inC,[RP.selectedPoles{:}]));
                if ~isempty(RP.selectedPoles) && RP.selectedPoles{1}(1)
                    RP.selectedPoles = [RP.selectedPoles num2cell(sP)];
                elseif ~isempty(RP.poles)
                    RP.selectedPoles = num2cell(sP);
                end
            end
            RP.up2date(2) = false;
        end
        function value = get.expansionPoints(RP)
            value = RP.expansionPoints_;
        end
        
        function set.poles(RP, poles)
            poles = sort(poles,'ComparisonMethod','real');
            if isequal(RP.poles, poles), return; end
            RP.poles_ = poles; 
            RP.up2date(3:4) = true;
            if isempty(RP.contours), return; end
            RP.newBackground = false;
            sP = find(RP.inside(poles));
            RP.selectedPoles = num2cell(sP);
        end
        function value = get.poles(RP)
            value = RP.poles_;
        end
        
        function set.selectedPoles(RP, selectedPoles)
            if isequal(selectedPoles, RP.selectedPoles), return; end
            selectedPoles = cellfun(@(x){x(:).'},selectedPoles);
            RP.selectedPoles_ = selectedPoles;
            n = max(length(selectedPoles),1);
            if length(selectedPoles)~=length(RP.nPoints)
                if any(RP.nPoints_~=RP.nPoints_(1))
                    warning(['The value of property '...
                        '''nPoints'' has been changed.'])
                end
                RP.nPoints_ = repmat(min(RP.nPoints),1,n);
            end
            if length(selectedPoles)~=length(RP.radius)
                if any(RP.radius_~=RP.radius_(1))
                    warning(['The value of property '...
                        '''radius'' has been changed.'])
                end
                RP.radius_ = repmat(min(RP.radius),1,n);
            end
            if isfield(RP.plotSettings,'fignumber') && ...
                    n~=length(RP.selectedPoles)
                RP.plotSettings = rmfield(RP.plotSettings,'index');
            end
            if ~isempty(selectedPoles)
                if ~selectedPoles{1}(1), selectedPoles = {}; end
                selectedPoles = [selectedPoles{:}];
                tolerance = min(1e-2,RP.precision(1))*RP.minD;
                inC = RP.inside(RP.poles(selectedPoles),tolerance);
                if all(inC) && sum(inC)==sum(RP.inside(RP.poles))
                    RP.getContours('updatenobg');
                else
                    RP.getContours('update');
                end
                if ~RP.up2date(4)
                    RP.up2date(3) = false; RP.residues = struct;
                end
                if isempty(selectedPoles)
                    RP.up2date(3) = true;
                elseif length(selectedPoles)<sum(RP.inside(RP.poles))
                    warning('Unselected poles inside the contour')
                end
            end
        end
        function value = get.selectedPoles(RP)
            value = RP.selectedPoles_;
        end
        
        function set.contours(RP, contours)
            try plt = ishandle(RP.plotSettings.fignumber.contours); 
            catch, plt = false; end
            if isempty(RP.contours_) || isequal(contours, RP.contours_)
                RP.contours_ = contours;
                if plt, RP.plot('contours'); end
                return;
            end
            RP.contours_ = contours; 
            if plt, RP.plot('contours'); end
            if isempty(RP.contours_f), return; end
            RP.undefinedValues = cell(1,length(contours));
            RP.definedValues = cell(1,length(contours));
            c = cellfun(@(x){x(:)},RP.contours_f); c = cat(1,c{:});
            c = cat(2,real(c),imag(c));
            for it = 1:length(contours)
                c_it = contours{it}(:);
                if isempty(c_it), continue; end
                c_it = cat(2,real(c_it),imag(c_it));
                [def,loc] = ismembertol(c_it, c, 'ByRows', true);
                RP.undefinedValues{it}=~def.';
                RP.definedValues{it} = loc(def).';
            end
            RP.up2date(1) = false;
        end
        function value = get.contours(RP)
            value = RP.contours_;
        end
        
        function set.precision(RP, p)
            this = 'RieszProjection.precision';
            if isscalar(p)
                if ~isfloat(p) && isreal(p) && (p > 0)
                    error('Invalid value for %s', this);
                end
                RP.precision_ = [p p];
            elseif numel(p)==2
                pabs = p(1); prel = p(2);
                if ~isfloat(pabs) && isreal(pabs) && (pabs > 0)
                    error('Invalid absolute tolerance for %s', this);
                elseif ~isfloat(prel) && isreal(prel) && (prel > 0)
                    error('Invalid relative tolerance for %s', this);
                end
                RP.precision_ = [pabs prel];
            end
        end
        function value = get.precision(RP)
            if RP.precision_(1)==RP.precision_(2)
                value = RP.precision_(1);
            else
                value = RP.precision_;
            end
        end
        function set.includeConjugatedPoles(RP,iCP)
            if iCP==RP.includeConjugatedPoles, return; end
            RP.includeConjugatedPoles_ = iCP;
            RP.getContours('update');
        end
        function value = get.includeConjugatedPoles(RP)
            value = RP.includeConjugatedPoles_;
        end
        function set.useEigenvalues(RP,uEV)
            if uEV==RP.useEigenvalues, return; end
            RP.useEigenvalues_ = uEV;
            if isempty(RP.contours), return; end
            if uEV
                RP.contours = RP.contours(end);
                RP.weights = RP.weights(end);
            else
                RP.getContours('updatenobg');
            end
        end
        function value = get.useEigenvalues(RP)
            value = RP.useEigenvalues_;
        end
        function set.nEigenvaluesExpected(RP,nEv)
            if nEv==RP.nEigenvaluesExpected, return; end
            RP.nEigenvaluesExpected_ = nEv;
            if ~isempty(RP.poles)
                RP.up2date(3:4) = false;
            end
        end
        function value = get.nEigenvaluesExpected(RP)
            value = RP.nEigenvaluesExpected_; 
        end
    end
    
    properties (Access=private)
        % See dependent properties
        nPoints_ = 4;
        nPointsB_ = 64;
        radius_ = 1e13;
        minD_ = 2e13;
        expansionPoints_;
        poles_;
        selectedPoles_ = {};
        contours_ = {};
        precision_ = Inf(1,2);
        useEigenvalues_ = false;
        includeConjugatedPoles_ = false;
        nEigenvaluesExpected_ = 5;
        tmp_; % for temporary storage
        
        % contours_f - The contours corresponding to the current fields
        % which are updated at each call of evaluate.
        contours_f;
        
        % nModes - The number of modes contributing to the modal picture
        % This property is updated at each call of evaluate.
        nModes = 0;
        
        contourSettings (1,1) struct = struct;
        plotSettings (1,1) struct = struct('fignumber',struct);
        qKeys (1,1) struct = struct;
        rKeys (1,1) struct = struct;
        exportkeys (1,1) struct = struct;
        w0 (1,1) struct = struct;
        % [integr. points, expansion, eigenvalues, from external source]
        up2date (1,4) logical = [false false false true];

        % Values which are not defined in fields are marked and
        % locations of already existing values are saved in definedValues
        undefinedValues (1,:) cell;
        definedValues (1,:) cell;
        
        % convergenceData - Convergence data is stored in this struct
        % The fieldnames refer to the quantity
        convergenceData (1,1) struct;
        
        % shape - saves information about the shape of the background
        % shape{1} - circle: [c r], ellipse: [0 1], polygon: vertices
        % shape{2} - 'c', 'e' or 'p'
        % shape{3} - order of Gauss-Kronrod quadrature rule
        % shape{4} - ellipse parameters or vertices without splits (pol.)
        shape (1,:) cell;
        
        % residues - If the eigenvalues are calculated the residues are
        % determined as well and can be used for expansions.
        residues (1,1) struct;
    end
    
    methods
        function RP=RieszProjection(f, varargin)
            % Construct an instance of the class
            %    rp = RIESZPROJECTION(f) Constructs an instance with
            %       default properties based on evaluations of f at
            %       complex frequencies.
            %    rp = RIESZPROJECTION(f,parameters) Constructs an instance
            %       whose properties are set according to the fields of
            %       the struct parameters. 
            %    rp = RIESZPROJECTION(f,PARAM1,VAL1,PARAM2,VAL2,...)
            %       Instead of a struct you can set the properties of this
            %       class using name value pairs where the names must be
            %       unique caseinsensitive partial matches of the property
            %       names
            %        
            % The function f defines the computationally expensive part of 
            % the function to be integrated. The analytic part is given by
            % the function handle g whose value is fixed for Riesz
            % projections but can be chosen freely if required for special
            % cases.
            
            RP.f = f;
            if isa(f,'Scattering')
                RP.dipoleBulkEmission = @f.dipoleBulkEmission;
            end
            if nargin==1, return; end
            ps = ['linearQ';'quadraticQ';properties(RP)];
            ps_ = ps(1:find(strcmp(ps,'f'))-1);
            if isstruct(varargin{1})
                fs = fieldnames(varargin{1});
                vs = struct2cell(varargin{1});
            elseif length(varargin)>1
                fs = varargin{1:2:end};
                vs = varargin{2:2:end};
            end
            for it = 1:length(fs)
                p = validatestring(fs{it},ps_);
                RP.(p) = vs{it};
            end
        end
        
        % selectPoles - Select poles interactively
        selectPoles(RP, radius)
        
        % getContours - Get the contours according to the selected poles    
        getContours(RP, varargin);
        
        % Evaluate contour integrals for given expansion points omega,
        % quantity can be 'field' or 'energy'.
        out = integrate(RP, varargin)
        
        % Plot contour and poles
        varargout = plot(RP, varargin)
            
        % Expand the given quantity (field or energy).
        expand(RP, varargin)
        
        % Get the modal field corresponding to the contour with the given
        % index and at the specified expansion point(s) w0.
        viewFields(RP, varargin)
        
        % Plot modal radiation diagrams at specified expansion point(s) w0
        viewRadiation(RP, index, w0, theta, phi, selection)
            
        % Calculates the reference solution on the real axis.
        getReferenceSolution(RP, quantity, w0, keys)
        
        % Calculate the error of the expansion. The maximum error along all
        % expansion points is returned. 
        [err,x] = error(RP, quantity, index, kind)
        
        % Get directory of integration point or reference solution by 
        % contour index and the frequency of the point. The first argument 
        % is the index. If there are n modal contributions the index n+1 
        % refers to the background contour and the index n+2 to the
        % reference solution. The path to the scattering problem with the
        % frequency closest to the given one is returned.
        varargout = getDir(RP,index,omega)
    end
    
    methods (Static)
        % Get the circle with minimal area enclosing all points.
        [c, r] = minimalCircle(points);
        % Get the circle with minimal area enclosing all points.
        E = minimalEllipse(points);
    end
    
    methods (Access=private)
        % Evaluate all function values on the contours.
        evaluate(RP)
        
        % Evaluate derived quantities.
        updateDerivedQuantities(RP, d)
        
        % compute the eigenvalues
        res = computeEigenvalues(RP,values,vars)
        
        % apply quadrature rule
        [res, convData] = quad(RP, values, vars, conv, it)
        
        % plot contours
        fig = plotContours(RP, vars)
        
        % plot error
        fig = plotError(RP, vars)
        
        % plot expansion
        fig = plotQuantity(RP, vars)
        
        % plot radiation pattern
        fig = plotRadiationPattern(RP, vars)
        
        % get a list of all quantities that can be expanded
        function qs = getQuantities(RP)
            if ~isa(RP.f,'Scattering')
                qs = unique([RP.linearQ RP.quadraticQ...
                    fieldnames(RP.expansions).']); 
                return;
            end
            addQs = {'DipoleEmission' 'NormalizedDecayRate'};
            if RP.includeConjugatedPoles
                qs = dir(RP.f.postproDir); qs = {qs.name};
                s = endsWith(qs,'.jcmpt')&~startsWith(qs,'ref');
                qs = qs(s);
                qs = cellfun(@(x)x(1:end-6),qs,'UniformOutput',0);
                qs = [addQs qs];
                if any(strcmpi(qs,'Radiation'))
                    nq = {'DipolePowerCollectionEfficiency' ...
                        'PhotonCollectionEfficiency'};
                    qs = [qs nq];
                end
            else
                qs = addQs;
            end
        end
        
        % check wether points are inside the contour and keep the distance
        % returns a vector of booleans
        function in = inside(RP,ps,tol)
            if isempty(ps) || isempty(RP.contours), in = true; return; end
            sh = RP.shape{2}; if nargin==2, tol = [0 0]; end
            if sh=='c'
                cr = RP.shape{1};
                in = abs(ps-cr(1))/(cr(2)+tol(1))<=1;
            elseif sh=='p'
                ps1 = RP.shape{4}; s = diff(ps1);
                n = s./abs(s)*exp(-1i*pi/2);
                a = abs(angle([n(2:end) n(1)].*exp(-1i*angle(n))));
                v = ps1(2:end); a_ = angle(n)+a/2;
                ps1 = v+1./cos(a/2).*(tol(1)*cos(a_)+1i*tol(2)*sin(a_));
                in = inpolygon(real(ps),imag(ps),real(ps1),imag(ps1));
            elseif sh=='e'
                A = RP.shape{4}; ps = (ps-A(1)-1i*A(2))*exp(-1i*A(5));
                a = real(ps)/(A(3)+A(4)+tol(1)); 
                b = imag(ps)/(A(3)-A(4)+tol(1));
                in = a.^2+b.^2<=1;
            end
        end
    end
    
    methods (Static,Access=private)
        % unit conversion
        function [lambda, xunit] = unit(w,xunit)
            c0 = RieszProjection.c0;
            switch xunit
                case 'nm'
                    lambda = 2e9*pi*c0./w;
                case 'um'
                    lambda = 2e6*pi*c0./w;
                    xunit = '\mu m';
                case 'mm'
                    lambda = 2e3*pi*c0./w;
                case 'm'
                    lambda = 2*pi*co./w;
                otherwise
                    lambda = w;
                    xunit = 'Hz';
            end
        end
    end
end

% Minimal test for f
function mustTakeCell(f)
if isempty(f), return; end
try a=f({}); catch, a=1; end
if ~(isequal(a,0) || isempty(a))
    error(['f should be callable taking cell arrays of doubles with '...
        'shape=(1,:) containing complex frequencies and return a '...
        ' cell array of the same size containing the results of the'...
        ' linear systems. At a second call with the previous output as'...
        ' input plus a character array defining a quantity, this'...
        ' quantity should be derived from the solutions of the'
        ' linear systems. f({}) must return 0']);
end
end
