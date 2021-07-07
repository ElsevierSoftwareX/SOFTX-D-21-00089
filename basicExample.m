% Minimal example to get familiar with the method of Riesz projection
% expansions.
% Please have in mind that you can always use the help function of matlab
% to get detailed discriptions of the properties and methods you are
% interested in. E.g., execute help RieszProjection to get a list of its
% methods. You can follow the corresponding link to get more details about
% a particular method.

% clear the workspace
clearvars
clc

%% Set up the matlab interface to JCMsuite
% Activate the third party support and start a daemon. Please refer to the
% README for details.

% Provide the path to the installation directory of JCMsuite and add it to
% the matlab path. The version shoud be 4.4.0 or higher.
jcm_root = '/path/to/local/installation/JCMsuite.x.x.x';
addpath(fullfile(jcm_root, 'ThirdPartySupport', 'Matlab'));

% start a daemon which handles the jobs submitted to jcmwave_solve
options = struct('Hostname', 'localhost', ...
                'Multiplicity',1, ...
                'NThreads',2 ... 
                );
% Shutdown a possibly running daemon and register a new computer ressource.
jcmwave_daemon_shutdown;
jcmwave_daemon_add_workstation(options);

%% Solve the eigenvalue problem
% There is a contour integral solver implemented in the RieszProjection 
% class that allows for obtaining the eigenvalues from a contour enclosing 
% the region of interest in the frequency plane, which can be recycled for 
% the expansion. Yet, it does not work properly if more than five 
% resonances within the contour couple to the target quantity. Even though
% it perfectly works for the given example, for the sake of simplicity, the
% built-in eigensolver of JCMsuite is used here. For the use of this
% code without a-priori knowledge of the eigenvalues, please refer to the
% example singleContour.m.

% Set the finite element degree for the resonance problem defined in
% 'resonance/project.jcmpt'
keys_res.finiteElementDegree = 4;

% In this example, we are interested in frequencies in the range from 
% 3.5e15 to 4.5e15 Hz. JCMsuite will use the following guess as a starting 
% point.
keys_res.guess = 4.0e15;

% The number of resonance frequencies provided with this field of the struct
% keys determines when the calculation is finished.
keys_res.numberEigenvalues = 10;

% The results returned by jcmwave_daemon_wait are not used in this script
% but loaded directly from file. Once it exists, you can comment out this
% block.
res = jcmwave_solve(['resonance' filesep 'project.jcmp'],keys_res);
[res,log] = jcmwave_daemon_wait(res);

% Wait until the result is written to a file.
resDir = ['resonance' filesep 'project_results'];
while ~exist([resDir filesep 'eigenvalues.jcm'],'file'), pause(1); end

%% Construct instance of Scattering
% This object takes care of solving the scattering problems for different
% frequencies while all other parameters remain constant. You can add all
% parameters you might need for parameter substitution to the structure
% 'keys' which is passed to the constructor together with the path to the
% project file.

% The keys must at least have the field 'finiteElementDegree'. If you want
% to simulate a point source and want to expand its dipole emission Gamma
% you must define position and strength in order to use the expression
% Gamma(w) = -real(E(w,position).'*strength) which assumes the strength to
% be real. Alternatively, you can expand  Gamma using the electromagnetic
% field energy flux, which is quadratic in the electric field strength E.

% Set the finite element degree for the scattering problem defined in
% 'example/project.jcmpt'.
keys.finiteElementDegree = 4;

% In this example, the dipole emission Gamma will be expanded. As stated
% above, in order to do this based on a point evaluation of the electric
% field, you have to define a strength vector, i.e., amplitude and
% polarization of the point source and its position within the 
% computational domain. Both parameters are used in the source file 
% 'example/sources.jcmt' and the latter additionally defines the position 
% where E must be evaluated.
keys.position = [0 145e-9 0];
keys.strength = [1 0 0];

% Path to the project file. The line of code that follos this comment
% assumes, that it is located in the directory 'example' whose parent is
% the current working directory. If this is not the case, please enter
% either the relative or the absolute path to the project file.
projectFile = ['example' filesep 'project.jcmpt'];

% Call the contructor of Scattering
sc = Scattering(projectFile, keys);

%% Construct an instance of RieszProjection
% You can pass a struct with options to the constructor which will be
% assigned to the corresponding properties. Alternatively, you can set them
% after construction. For a complete list of public properties, type
% help RieszProjection in the command window.

% The eigenmodes computed previously have been saved in the file
% 'resonance/project_results/eigenvalues.jcm'. They are loaded and assigned
% to the field 'poles' of the parameter struct. If no poles are given, the
% eigenvalue problem is solved with a contour integral solver. Please refer
% to the example 'singleContour.m'.
evDir = ['resonance' filesep 'project_results'];
ev = jcmwave_load([evDir filesep 'eigenvalues.jcm']);
parameters.poles = ev.eigenmode; 
% You can view the eigenmodes executing in the command window: 
% jcmwave_view([evDir filesep 'fieldbag.jcm'])

% In this example, we are only interest in the expansion of the dipole
% emission given as an expression which is linear in the electric field.
% Symmetry with respect to the real axis is therefore not required and the
% property 'includeConjugatedPoles' can be set to false. If you want to 
% expand quantities which are quadratic in E, such as the electric field 
% energy, this property must be true. For further details refer to the 
% example quadraticForms.m.
parameters.includeConjugatedPoles = false;

% The radius of the circles around poles which are expected to contribute
% to the target quantity. It has to be small enough to avoid the enclosure
% of additional poles and large enough to ensure that, given its error, the 
% corresponding eigenvalue lies inside the contour.
parameters.radius = 1e13; 

% The value will be added to the the axes of the ellipse or to the radius
% of the circle. In the latter case, it defines the minimal distance of the
% background contour to the enclosed poles. In the former case, also a two
% element vector can be given. The first value is then added to the
% semi-major axis and the second to the semi-minor axis. The convergence is
% optimal if the distance of the background contour to any pole is maximal.
% Please be aware that you have to control the distance to poles outside it
% manually, i.e, when setting the property 'minD' you must check if the
% distances to poles outside the contour are large enough using the plot
% function.
parameters.minD = [2e14 3e14];

% Number of integration points for the small contours.
parameters.nPoints = 4;

% Number of integration points for the large contour (background contour).
parameters.nPointsB = 32;

% You can change this value from its default (Inf) to allow for an adaptive
% refinement of the integration points on the contours until the target
% precision is reached. The convergence of the individual contours is
% estimated comparing the best numerical solution with less
% integration points and the sum of all modal components is compared to
% selected direct solutions on the real line if refernce solutions are
% demanded. This can be done either by setting the property 
% 'referencePoints' (see below) or by calling the method 
% 'getReferenceSolution' directly. 
parameters.precision = Inf;

% If the property 'referencePoints' is not empty, the scattering problems
% at the specified real frequencies are solved and compared to the
% expansion. 
parameters.referencePoints = [3.8e15 4e15 4.15e15 4.3e15 4.5e15];

% The expansion points define the positions where the expansion is
% evaluated. As additional expansion points do not require solutions of
% additional scattering problems, the number of expansion points does not
% significantly affect the computational costs. 
parameters.expansionPoints = 3.5e15:1e13:4.5e15;

% Construct an instance of the class RieszProjection.
rp = RieszProjection(sc, parameters);

%% Define the contours
% The contours are defined based on a selection of poles which are expected
% to contribute to the target quantity and the expansion points. The
% individual poles are enclosed in circular contours, the background must 
% enclose all selected poles as well as the expansion points. It can have 
% elliptic, circular or polygonial shape.

% The selected poles are given as a cell array containing indices to the
% property 'poles'. These can be scalars or vectors if the contributions of
% several poles shall be summed up. You can interactively select poles with
% the mouse executing: rp.selectPoles in the command window. 
rp.selectedPoles = {3 5 6};

% Get the contours. By default, the background contour is an ellipse. For
% more options: execute help rp.getContours in the command window. 
rp.getContours

% Plot contours, expansion points, and poles. If you omit the argument, a
% list of quantities that can be plotted will be printed to the command
% window and you can select with the mouse. 
rp.plot('contours')

%% Expand
% Expand the dipole emission. In order to see a list of quantities which 
% can be expanded, omit the argument. Be aware that if you have set the 
% property includeConjugatedPoles to false only quantities are listed which 
% are linear in E. 
rp.expand('DipoleEmission')

%% View results
% You can view the results using the plot function. For a list of
% quantities that can be plotted, execute rp.plot in the command window. 
rp.plot('NormalizedDecayRate','legendLocation','northwest')
rp.plot('error')

%% Visualize the total field
% You can visualize the whole field as a function of the real frequency. By
% default, a section through the center of the x-component of the electric 
% field strength is exported on a cartesian grid and then visualized. The 
% rows correspond to contours. If there are n modal contributions, n+1 is 
% the index of the background contour, n+2 the index of the summed
% contributions, and n+3 is the index of the reference solutions. The 
% colums correspond to the real frequencies. 

% As the expansion of the dipole emission shows that the first mode merely
% contributes to the dipole emission, we define a selection, containing
% indices to modes which we want to be visualized together, i.e., summed
% up. In this case first mode and background. Their sum is displayed in the 
% row with index n+2-length(selection). Consequently the fields in the
% first row correspond to the second mode which dominates in the given
% range. The second row shows the third mode which gains relevance at
% higher frequencies. The third row shows the sum of first mode and
% background. Fourth and fith row show the sum of the first three rows and
% the direct scattering solutions at the real line, respectively. 
% Negative values are displayed in blue and cyan, positive values in yellow
% and red. 
selection = [1 3 4];
rp.viewFields(rp.referencePoints,1:6,'selection',selection)

% In order to remove all files that have been added to project directory
% you can call sc.clean.
