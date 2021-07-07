% This example demonstrate how in case of few resonances within the region
% of interest all the required information can be extracted from the 
% solution of scattering problems on a single contour. 
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

% Give the path to the installation directory of JCMsuite and add it to
% the matlab path. The version should be 4.4.0 or higher.
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

%% Construct an instance of Scattering
% Define keys which are passed to jcmwave_solve and call the constructor.
% As in the basic example, the light source is a point source.
keys.finiteElementDegree = 4;
keys.position = [0 145e-9 0];
keys.strength = [1 0 0];

% Path to the project file. The line of code that follos this comment
% assumes, that it is located in the directory 'example' whose parent is
% the current working directory. If this is not the case, please enter
% either the relative or the absolute path to the project file.
projectDir = 'example';
projectFile = [projectDir filesep 'project.jcmpt'];
sc = Scattering(projectFile, keys);
sc.resultbagDir = [projectDir filesep 'resultbags_sC'];

%% Construct an instance of RieszProjection
% As in the example for quadratic forms we aim to expand quantities
% quadratic in the electricfield and therefore need to consider the complex
% conjugates of the (unknown) poles. 
parameters.includeConjugatedPoles = true;

% As we aim to work with a single contour we can focus on the parameters 
% for the background contour. Once more we are interested in a range from 
% 3.5e15 to 4.5e15 Hz and we want to focus on resonances with Q-factors 
% larger than ten. The range is defined by the expansion points and we add 
% 2e14 to the semi-major axis to get some distance between these arteficial
% poles and the contour. As there are no poles to span the ellipse the 
% second value of the property 'minD' defines the semi-minor axis which we 
% set to 5e14. Consequently the smallest possible Q-Factor of a pole inside
% the contour is larger than eight. 
parameters.minD = [5e14 5e14];
parameters.nPointsB = 32;
parameters.referencePoints = [3.8e15 4e15 4.15e15 4.3e15 4.5e15];
parameters.expansionPoints = 3.5e15:1e13:4.5e15;

% The contour integral algorithm used to solve the eigenvalue problem
% requires an esitmate of the number of eigenvalues expected to be found
% inside the contour. Currently the number of modes which couple to the
% source should not exceed five for relyable results. However improvements
% are planned for future versions. Either with a variation of the current
% algorithm or implementing another contour integral method (Beyn, FEAST).
parameters.nEigenvaluesExpected = 4;

% Construct an instance of the class RieszProjection.
rp = RieszProjection(sc, parameters);

%% Define the contours
% The contour is contructed based on the expansion points and the values
% defined in the property 'minD'.
rp.getContours

% Plot contours and expansion points. 
rp.plot('contours')

%% Expand and solve eigenvalue problem
% The eigenvalues are obtained solving a nonlinear system of equations with
% 2*n unknowns (where n is the number of expected eigenvalues). In
% principle n unknowns would suffice but using the residues as additional
% unknwons allows for an easy access to the derivatives which speed up the
% computation considerably. Additionally the residues anable to directly
% evalutate the expansion without further integrations. For further
% quantities the exact positions of the eigenvalues are used to obtain the
% expansion from the existing contour. 
rp.expand({'DipoleEmission' 'Radiation'})

% Replotting the results shows the eigenvalues, which have been found. The
% fourth solution has been removed as the normalized residue is smaller
% than a threshold which is either given by the property precision or at
% most 1e-2.
rp.plot('contours')

%% View expansions
% You can view the results using the plot function. For a list of
% quantities that can be plotted, execute rp.plot in the command window. 
rp.plot('NormalizedDecayRate','legendLocation','northwest')
rp.plot('DipolePowerCollectionEfficiency')
rp.plot('error')

% Comparison to the example 'quadraticForms.m' shows that the error differs
% by an order of magnitude. Yet, one has to bear in mind that for each of
% the smaller contours eight scattering problems had to be solved. The
% advantage of knowing the eigenvalues beforehand is the possibility to
% place the contour in an optimal way, such that its distance to the poles
% is maximized. 
