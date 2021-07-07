classdef Scattering < handle
    %SCATTERING Evaluation of scattering problems with JCMsuite
    %   This class parallelizes the evaluation of scattering problems with
    %   all parameters fixed except the frequency. An instance sc of this
    %   class can be used as a callable which takes a cell array of 
    %   doubles and returns tags representing results in the property
    %   'resultbags'. sc({}) returns 0. sc(contours) calls the hidden
    %   method 'f' which passes the scattering problems to JCMsolve. A
    %   second call with tags and a character array specifiying a quantity
    %   as input arguments will apply post processes to the scattering
    %   solutions.
    %   If the property 'keys' has a field 'position', this is assumed to
    %   be the position of a dipole located within the computational
    %   domain and the field value of the scattered field at its position
    %   is exported. Please refer to the <a href="matlab:web( ...
    %   ['https://docs.jcmwave.com/JCMsuite/html/ParameterReference' ...
    %   '/d77be858da38e480d80e9f59c1a3f25f.html'],'-browser')"
    %   >parameter reference</a> 
    %   of JCMsuite if you want to learn more about the definition of  
    %   scattering problems. The project file, e.g., project.jcmpt, must
    %   contain the section PML { %(pml)s } and the primitive 
    %   FiniteElementDegree = %(finiteElementDegree)e for parameter
    %   substitution. Please refer to the documentation of the
    %   <a href="matlab:web(['https://docs.jcmwave.com/JCMsuite/html/' ...
    %   'MatlabInterface/index.html'],'-browser')" >matlab interface</a>.
    %   Your custom parameters for keyword substitution can be passed via
    %   the property 'keys' of this class which must at least have the
    %   field 'finiteElementDegree'. Apart from the project file, you need 
    %   the files <a href="matlab:web(['https://docs.jcmwave.com/' ...
    %   'JCMsuite/html/ParameterReference/b61236968b3822be5ffbfee' ...
    %   '6564f23da.html'],'-browser')"
    %   >layout.jcm</a> or layout.jcmt, <a href="matlab:web(['https:' ...
    %   '//docs.jcmwave.com/JCMsuite/html/ParameterReference/' ...
    %   '3df274a2924c89630ff2393cc22b686e.html'],'-browser')"
    %   >materials.jcm</a> or materials.jcmt
    %   and <a href="matlab:web(...
    %   ['https://docs.jcmwave.com/JCMsuite/html/ParameterReference/' ...
    %   'f3e666a5067147d3cd45b67773bb77ae.html'], '-browser')"
    %   >sources.jcmt</a>. The latter has to contain the 
    %   primitive Omega = %(omega)e and in case of a <a href="matlab:
    %   web(['https://docs.jcmwave.com/JCMsuite/html/Parameter' ...
    %   'Reference/67f12453f439cffbc214bb0a2362a853.html'],'-browser')"
    %   >point source</a>
    %   additionally Position = %(position)e and Strength = %(strength)e
    %   are required which then must also be fields of the property 'keys'.
    %
    %SCATTERING Properties:
    %   keys - Struct passed to JCMsolve for parameter substitution
    %   projectFile - Path to a JCMsuite template file
    %   The file must end with .jcmpt and its location must be a project
    %   directory containing all files needed for the scattering problem.
    %   pointList - Points to be exported
    %   If a list of points is given, a call to this class will add a post
    %   process to the project file which exports the corresponding points.
    %   If the property 'keys' contains the fields position and strength 
    %   for a dipole source, this property is ignored and the position of
    %   the dipole is exported instead.
    %   quantity - quantity to be exported
    %   You can specify the quantity that you want to have exported. Be
    %   aware that it must be a holomorphic quantity to allow for Riesz
    %   projections.
    %   pml - Path to PML logfile
    %
    %JCMSCATTERING Methods: 
    %   Scattering - Constructor
    %   clean - Delete Files that have been created using the class
    %   'RieszProjection'
    %
    %DEPENDENCIES:
    %   The JCMwave third party support has to be setup for this class.
    %   And you have to run a daemon in order to apply post processes.
    %   Please refer to the documenation of the <a href="matlab:
    %   web(['https://docs.jcmwave.'com/JCMsuite/html/MatlabInterface' ...
    %   '/index.html'],'-browser')">matlab interface</a> to set up
    %   the third party support and to the homepage of <a href="matlab:
    %   web('https://jcmwave.com/jcmsuite/jcmsuite-documentation', ...
    %   '-browser')">JCMwave</a> for
    %   installation instructions.
    %
    %   See also RIESZPROJECTION

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021
    
    
    properties (SetAccess=private)
        keys (1,1) struct {mustHaveField} = ...
            struct('finiteElementDegree',0);
        projectFile (1,:) char {mustBeValidFile};
    end
    
    properties
        pointList (:,3) double;
        quantity char = 'ElectricFieldStrength';
        pml (1,:) char {mustBeValidFile};
    end
    
    properties (Hidden)
        % projectDir - The Directory where the project files are located
        projectDir (1,:) char {mustBeValidDirectory};
        
        % postproDir - The Directory where post processes are located
        postproDir (1,:) char {mustBeValidDirectory};
        
        % scatteringDir - Working directory for the scattering problems
        scatteringDir (1,:) char {mustBeValidDirectory}; 
        
        % referenceDir - Working directory for the reference solutions
        referenceDir (1,:) char {mustBeValidDirectory};
        
        % resultbagDir - The directory containing all resultbags
        resultbagDir (1,:) char {mustBeValidDirectory};
        
        % resultbags - Resultbags of the scattering problems
        resultbags (1,2) cell;
        
        % pprbs - Resultbags for the postprocesses. 
        pprbs (1,1) struct;
        
        % bulkEmission - Dipole bulk emission of the given setting
        bulkEmission (:,2) double;
    end
    
    methods
        function sc = Scattering(projectFile, keys)
            %SCATTERING Construct an instance of this class
            %   SCATTERING(projectFile, keys) constructs an object with
            %   working directories located in the parent directory of the
            %   .jcmpt file 'projectFile', which must define a scattering
            %   problem. The struct keys contains fields which are
            %   used for the parameter substitution of JCMsuite. It must a
            %   least contain the field 'finiteElementDegree'.
            
            % make shure all keys exist, that could be relevant
            fields = Scattering.fields;
            for it = 1:length(fields)
                if isfield(keys,fields{it})
                    error(['The fieldname %s is reserved for ', ...
                        'internal use.'], fields{it})
                else
                    keys.(fields{it}) = [];
                end
            end
            sc.keys = keys;
            sc.setPaths(projectFile);
        end
        
        function clean(sc,deleteResultbags)
            %CLEAN Delete files
            %   CLEAN(sc) Removes all scattering problems in the working
            %   directories integration_points and reference_points as well
            %   as the directory 'superpositions' and 'bulk_emission' if
            %   they exist. The file 'pml.log' is deleted too. If you want
            %   to use a custom file which you do not want to be deleted,
            %   you should place it outside the project directory and set
            %   the property 'pml' to the file path.
            %
            %   CLEAN(sc, deleteResultbags) If deleteResultbags is false
            %   the directory containing the resultbags is kept. The
            %   default is true.
            
            if nargin<2, deleteResultbags = true; end
            sc.resultbags = cell(1,2);
            if exist(sc.scatteringDir,'dir')
                rmdir(sc.scatteringDir,'s');
                mkdir(sc.scatteringDir);
            end
            if exist(sc.referenceDir,'dir')
                rmdir(sc.referenceDir,'s');
                mkdir(sc.referenceDir);
            end
            if deleteResultbags && exist(sc.resultbagDir,'dir')
                rmdir(sc.resultbagDir,'s');
                mkdir(sc.resultbagDir);
                sc.resultbags = cell(1,2);
            end
            sp = [sc.projectDir filesep 'superpositions'];
            if exist(sp,'dir'), rmdir(sp,'s'); end
            bk = [sc.projectDir filesep 'bulk_emission'];
            if exist(bk,'dir'), rmdir(bk,'s'); end
            if exist(sc.pml,'file'), delete(sc.pml); end
        end
        
        emission = dipoleBulkEmission(sc,omega);
        
        function set.bulkEmission(sc,bE)
            if size(sc.bulkEmission,1)>size(bE,1)
                idx = ismembertol(bE(:,1),sc.bulkEmission(:,1));
                bE = sortrows([bE(~idx,:); sc.bulkEmission]);
            end
            sc.bulkEmission = bE;
        end
    end
    
    methods (Hidden)
        function varargout = subsref(sc, S)
            %SUBSREF Subscripted reference
            %   To allow for calls sc(contours) of an istance sc of this 
            %   class the subsref method has to be redefined. sc(contours)
            %   calls the method f which can be redefined and is supposed
            %   to solve scattering problems on given contours in the
            %   frequency plane.
            %
            %   See also Scattering.F, SUBSREF
            if strcmp(S(1).type, '()') && iscell(S(1).subs{1})
                varargout{1} = sc.f(S(1).subs{:});
            else
                if nargout>0
                    varargout{1:nargout} = builtin('subsref', sc, S);
                else
                    builtin('subsref', sc, S);
                end
            end
        end
        
        % This function is called by subsref and can be redefined for 
        % different scattering problems. The input parameter contours is 
        % expected to be a cell of size=(1,:) containing doubles of 
        % shape=(1,:).
        out = f(sc, varargin)
        
        out = sum(sc, tags, weights, omega0, outName)
        
        [ids,tags] = submit(sc, tags, varargin)
        
        v = collect(sc, ids, tags, q, last)
        
        varargout = export(sc, tags, pointListOrKeys, quantity)
    end
    
    methods (Static, Hidden)
        % get keys for radiation post process with quadgk
        function keys = radiationKeys(keys)
            % set defaults
            if ~isfield(keys,'radius'), keys.radius = 1e-5; end
            if ~isfield(keys,'NA'), keys.NA = [0.8 1]; end
            if ~isfield(keys,'step'), keys.step = 20; end
            if ~isfield(keys,'direction'), keys.direction = 'up'; end
            % get nodes for Gauss quadrature
            gknodes = RieszProjection.nodes15;
            nodes = [-gknodes(end:-1:1); 0; gknodes];
            % subdevide intervalls in phi
            subs = linspace(0,360,ceil(360/keys.step)+1);
            subs = [subs(1:end-1); subs(2:end)];
            midpt = sum(subs)/2; keys.halfhp = diff(subs)/2;
            keys.points_phi = reshape(nodes*keys.halfhp + midpt,1,[]);
            % subdevide intervalls in theta
            NA = asind(keys.NA); s = 0; step = keys.step;
            if strcmp(keys.direction,'down'), s = 180; NA = 180-NA; end
            s1 = linspace(s,NA(1),ceil(abs(s-NA(1))/step)+1); 
            s2 = s1(end); keys.ndx = length(s1); % end of the first NA
            if diff(NA)
                s2 = linspace(NA(1),NA(2),ceil(abs(diff(NA))/step)+1);
            end
            subs = [s1(1:end-1) s2]; subs = [subs(1:end-1); subs(2:end)];
            midpt = sum(subs)/2; keys.halfht = abs(diff(subs)/2);
            keys.points_the = reshape(nodes*keys.halfht + midpt,1,[]);
        end
        
        function keys = radiationPatternKeys(keys)
            % set defaults
            if ~isfield(keys,'radius'), keys.radius = 1e-5; end
            if ~isfield(keys,'gridPointsTheta')
                keys.gridPointsTheta = 1:90;
            end
            if ~isfield(keys,'gridPointsPhi')
                keys.gridPointsPhi = [-90 90];
            end
        end
    end
    
    methods (Access=private)
        function sc = setPaths(sc, projectFile)
            %SETPATHS Create directories to work in
            %   SETPATHS(sc, projectFile) Creates directories located in
            %   the parent directory of projectFile.
            
            pF = java.io.File(projectFile);
            if ~pF.isAbsolute
                pF = java.io.File([pwd filesep projectFile]);
            end
            pF = pF.getCanonicalFile;
            sc.projectFile = char(pF.getPath);
            pD = char(pF.getParent);
            sc.projectDir = pD;
            wDir = [pD filesep 'integration_points'];
            if ~exist(wDir,'dir'), mkdir(wDir); end
            refDir = [pD filesep 'reference_points'];
            if ~exist(refDir,'dir'), mkdir(refDir); end
            resDir = [pD filesep 'resultbags'];
            if ~exist(resDir,'dir'), mkdir(resDir); end
            sc.scatteringDir = wDir;
            sc.referenceDir = refDir;
            sc.resultbagDir = resDir;
            pD = mfilename('fullpath'); a = strfind(pD, filesep);
            pD = pD(1:a(end-1)-1);
            sc.postproDir = [pD filesep 'postprocesses'];
            if ~exist(sc.postproDir,'dir')
                warning(['The directory containing the post processes '...
                    'is expected to be found in the directory '...
                    'containing the folder with the definition of this '...
                    'class. You must either copy it there or set the '...
                    'hidden property postproDir of this class.'])
            end
        end
    end
    
    properties (Constant, Access=private)
        eps0 = 8.85418781762039e-12; % permittivity in vacuum
        mu0 = 4e-7*pi; % permeability in vacuum
        bags = {'integration_points.mat' 'reference_points.mat'};
        fields = {'addSingularFields' 'px' 'py' 'pz' 'quantity' 'omega'};
        
        pP = ['PostProcess {\n' ...
            '  BulkEmission {\n' ...
            '    Omegas=%%(omega)e\n' ...
            '    ProjectPath="%s"\n' ...
            '    OutputFileName="bulk_emissions.jcm"\n' ...
            '  } \n' ...
            '}'...
            ];
        
        superposition = [ ... superposition post rprocess
            'PostProcess { \n' ...
            '  Superposition { \n' ...
            '    OutputFileName = "%%(oPath)s"\n' ...
            '%%(fieldbags)s\n' ...
            '    Omega = %%(omega)e\n' ...
            '  } \n' ...
            '}\n'];

        fieldBag = [ ... field bags for superposition post process
            '    FieldBag { \n' ...
            '      FieldBagPath = "%s"\n' ...
            '      Weights = %.8e%+.8ei \n' ...
            '    }\n'];
        
        exportPoints = [ ... point list
            'PostProcess {\n' ...
            '  ExportFields {\n' ...
            '    FieldBagPath = "%%(fPath)s"\n' ...
            '    OutputFileName = "%%(oPath)s"\n' ...
            '    OutputQuantity = "%%(quantity)s"\n' ...
            '    AddSingularFields = %%(addSingularFields)s\n' ...
            '    PointList {\n' ...
            '      PointsX = %%(px)e\n' ...
            '      PointsY = %%(py)e\n' ...
            '      PointsZ = %%(pz)e\n' ...
            '    }\n' ...
            '  }\n' ...
            '}\n'];
                    
        exportCartesian = [ ... cartesian
            'PostProcess {\n' ...
            '  ExportFields {\n' ...
            '    FieldBagPath = "%%(fPath)s"\n' ...
            '    OutputFileName = "%%(oPath)s"\n' ...
            '    OutputQuantity = "%%(quantity)s"\n' ...
            '    Cartesian {\n' ...
            '      %s\n' ...
            '      %s\n' ...
            '      %s\n' ...
            '    }\n' ...
            '  }\n' ...
            '}\n'];
    end
    
end

function mustBeValidFile(filepath)
if ~exist(filepath, 'file') && ~isempty(filepath)
    error('The file \"%s\" does not exist', filepath);
end
end

function mustBeValidDirectory(filepath)
if ~exist(filepath, 'dir') && ~isempty(filepath)
    error('The directory \"%s\" does not exist', filepath);
end
end


function mustHaveField(keys)
if ~isfield(keys,'finiteElementDegree')
    error('The struct keys must contain the field finiteElementDegree');
end
end
