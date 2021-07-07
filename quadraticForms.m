% In this example, expansions of quadratic forms and quantities related 
% to the far field are demonstrated. The former are well suited for this 
% method as no crossterms appear. An expansion of the latter leads to
% finite results as the contour integration allows to circumvent the
% exponential divergence of quasi normal modes outside the resonator.

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

%% Solve the eigenvalue problem
% If you have executed the basic example, the required results already
% exists. Otherwise, you must uncomment and execute the following lines,
% which are a short version of the code in the corresponding section of the
% basic example.

% k = struct('finiteElementDegree',4,'guess',4e15,'numberEigenvalues',10);
% res = jcmwave_solve(['resonance' filesep 'project.jcmp'],k);
% [res,log] = jcmwave_daemon_wait(res);
% resDir = ['resonance' filesep 'project_results'];
% while ~exist([resDir filesep 'eigenvalues.jcm'],'file'), pause(1); end

%% Construct instance of Scattering and RieszProjection
% Define keys which are passed to jcmwave_solve and call the constructor.
% As in the basic example, the light source is a point source.
keys.finiteElementDegree = 4;
keys.position = [0 145e-9 0];
keys.strength = [1 0 0];

% Path to the project file. The line of code that follos this comment
% assumes, that it is located in the directory 'example' whose parent is
% the current working directory. If this is not the case, please enter
% either the relative or the absolute path to the project file.
projectFile =['example' filesep 'project.jcmpt'];
sc = Scattering(projectFile, keys);
sc.resultbagDir = ['example' filesep 'resultbags_qF'];

%% Construct an instance of RieszProjection
% You can pass a struct with options to the constructor, which will be
% assigned to the corresponding properties. Alternatively, you can set them
% after construction. For a complete list of public properties, type
% help RieszProjection in the command window.

% Set the eigenmodes which correspond to poles of the electric field.
evDir = ['resonance' filesep 'project_results'];
ev = jcmwave_load([evDir filesep 'eigenvalues.jcm']);
parameters.poles = ev.eigenmode;

% In this example, we will expand different quantities including the
% electric field energy and the radiation to the far field. Both are based 
% on functionals quadratic in the electric field. In order to avoid complex 
% conjugation, which does not meet the requirement of holomorphicity for 
% the application of the residue theorem, the complex conjugate of a field 
% E is replaced by the field evaluated at the negative frequency, i.e., 
% E^*(w) -> E(-w), which at real frequencies is the same. The physical 
% quantity is real in the time domain. Hence, E(-w) = E^*(w^*) := E^o(w)
% since E^o(w) has the same poles as E(w) but mirrored at the real line we
% have to consider these 'conjugated poles' as well. With contours being
% symmetric with respect to the real line, no extra scattering problems have
% to be solved to obtain E^o(w).
parameters.includeConjugatedPoles = true;

% All other parameters are the same as in the basic example
parameters.radius = 1e13; 
parameters.minD = [2e14 3e14];
parameters.nPoints = 4;
parameters.nPointsB = 16;
parameters.referencePoints = [3.8e15 4e15 4.15e15 4.3e15 4.5e15];
parameters.expansionPoints = 3.5e15:1e13:4.5e15;

% Construct an instance of the class RieszProjection.
rp = RieszProjection(sc, parameters);

%% Define the contours
% The same poles as in the basic example are selected. As the property
% 'includeConjugatedPoles' is true, the background contour will also
% enclose the mirrored poles above the real line.
rp.selectedPoles = {3 5 6};

% Get and plot the contours
rp.getContours
rp.plot('contours')

%% Expand
% If you execute rp.expand in the command window, a list of quantities is
% printed to the command window which can be expanded with the provided
% post processes. Additional to dipole emission and normalized decay rate,
% which have already been available in the basic example, quantities based 
% on quadratic forms are now included. Since the post process for the
% radiation to the far field includes the option to simultaneously evaluate
% a normalization and since the dipole emission is available as well,
% photon collection efficieny (i.e., the power radiated in a predefined
% numerical aperture (NA) devided by the total power radiated to the upper 
% hemisphere) and the dipole power collection efficiency (i.e., the power 
% radiated in the given NA devided by the total emitted power of the
% source) are appended to the list.

% You can pass all the quantities to expand at once. Here, dipole
% emission and electromagnetic field energy flux are essentially the same
% quantities. The first of them is based on the expression
% Gamma(w) = -real(E(w,position).'*strength), where a real-valued strength is
% defined as the source and therefore fixed even though it is related to E
% via the Maxwell's equations. The second is based on an expression
% quadratic in E namely the surface integral of the complex pointing vector 
% S (proportional to cross(E,rot(E)) over a closed surface containing the 
% dipole.
rp.expand({'DipoleEmission' 'ElectromagneticFieldEnergyFlux'})

% As mentioned above, the post process 'Radiation' allows for normalization.
% That means you can pass a keys structure to the method with a field 'NA'
% where you can define the NA you are interested in and the NA for
% normalization which for the photon collection efficiency would simply be
% one. The possibility to pass a keys structure to the expand method which
% is used for parameter substitution in the post processes can for instance
% also be used if you want to apply the post process to different domains 
% of the computational domain separately.
rp.expand('Radiation','keys',struct('NA',[0.8 1]))

%% Plot and control the error
% If you plot the error without specifying the index of a contour, the
% error of the physical quantity is evalutated. As the property 
% 'referencePoints' is not empty, the method 'expand' has called
% 'getReferencesolution' at these points. The reference solutions are used
% to estimate the error. The relative error with respect to them is
% evaluated and the maximum is plotted as a function of the number of
% integration points on the background contour. 
rp.plot('error')

%%

% If you want the error to be below a certain threshold, you can set the
% property 'precision' to this value. This property can be a scalar or a 
% vector with two elements. In the latter case, the first value is used as 
% upper bound of the absolute error and the second one as upper bound of
% the relative error. In the former case, both bounds have the same value.
% If we want to drop the relative error of the radiation below 1e-4 we
% should also set a value for the absolute error as it is of the order of
% 1e15 and values smaller than 1e11 can be considered close to zero and
% should be controlled by the absolute error. Otherwise modal contributions 
% close to zero will spoil the convergence. A call of the property expand
% will now refine the contours until the precision is reached or the number
% of integration points on any contourexceeds 1024. All other quantities
% which have already been expanded will be updated to take advantage of the
% new integration points but without initializing additional refinements if
% they do not reach the target precision. 
rp.expand('Radiation','precision',1e-3)
rp.plot('error','errortype','relative')

%% View results
% You can visualize results using the plot function. For a list of
% quantities that can be plotted, execute rp.plot in the command window and
% select the quantity you are interested in with the mouse.
% You can customize the plots passing name value pairs. For a complete list
% execute help RieszProjection/plot in the command window.
rp.plot('NormalizedDecayRate','legendLocation','northwest')

%% Change units
% If you want to change the current plot you can omit the quantity and just
% pass the corresponding keyword argument, e.g.:
rp.plot('xunit','nm','legendLocation','northeast','xlim',[416 540])

% The unit used for the x-axis is saved and used for other quantities as
% well. If you want to add the unit of a specific quantity, this can be
% done with the keyword argument 'yunit'.
rp.plot('Radiation','yunit','W')
% The default interpreter is LaTeX:
rp.plot('ElectromagneticFieldEnergyFlux','yunit','\frac{J}{s}')

%% Visualize radiation pattern
% To get an impression to what directions energy is radiated by a point
% source located in a nanoresonator, you can give discretizations of theta
% and phi in spherical coordinates which are used to create a cartesian 
% grid. The expansion of the radiation pattern is mainly used for 
% visualization purposes, therefore, it is expanded at the reference points 
% instead of the expansion points.
rad_keys.gridPointsTheta = 0:5:90;
rad_keys.gridPointsPhi = 0:5:180;
rp.expand('RadiationPattern','keys',rad_keys)

% The plot shows the energy density radiated in radial direction at the 
% reference frequencies. The plots must be read in spherical coordinates. 
% A large radius means that much energy is radiated in the corresponding 
% direction. The modal contributions can have negative radii. This can be
% understood as the suppression of radiation into the opposite direction. 
% The reference solution is displayed with faces instead of edges. 
rp.plot('RadiationPattern','fignumber',20,'legendLocation','northeast')

%% Radiation pattern in 2D 
% If you provide only two values of the azimuthal anlge, which define a 
% cross section, the plot function will display 2d plots.
rad_keys.gridPointsTheta = 0:1:90;
rad_keys.gridPointsPhi = [0 180];
rp.expand('RadiationPattern','keys',rad_keys)
rp.plot('RadiationPattern')
