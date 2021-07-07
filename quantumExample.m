%% problem description
% We want to solve the Schroedinger equation in 1D for a quantum particle 
% of effective mass m* = 1 given the potential V(x) = -10 for -L<=x<=L and 
% 0 otherwise. Furthermore, h-bar is set to sqrt(2) and L = pi/sqrt(2). A
% thourough description of the physical problem and its discretization
% can be found in the paper 'An eigenvalue method for open-boudary quantum 
% transmission problems' of Shao et al. (1995, doi: 10.1063/1.360132). The
% numerical constants are taken from Gavin et al.: 'FEAST eigensolver for
% nonlinear eigenvalue problems' (2018, doi: 10.1016/j.jocs.2018.05.006).
% 
% From the physical point of view, we are interested in the question, how
% the modes of the open system contribute to the transition amplitude at a
% potential well.
%
% This simple problem runs very fast and is well suited to get familiar
% with the capabilities of the class RieszProjection.
%
% This example covers two typical usecases. First, we assume that we are
% interested in a small energy interval, where we expect less than ten
% modes to dominate the spectrum. This allows for solving the eigenvalue 
% problem and performing the expansion with solutions of linear systems
% along a single contour. In the second case, the eigenvalues are given
% from an external source and the remaining task is the evaluation of the
% expansion. Here, the number of modes which are taken into account can be 
% much larger than in the first case. 

%% prepare the workspace

% In most cases, it is nice to start with a clean workspace.
% close all; clearvars; clc

% In order to get reproducible results, the random number generation must
% be reset before running this script;
rng(15)

% If your current folder is not the directory containing the class
% definition '@RieszProjection', you should add its parent directory to the
% matlab path.
% addpath('/path/to/your/copy')

%% construct an instance

% As we want to test the eigensolver, which is expected to work well for a
% small number of eigenvalues (about 5), we do not need to pass eigenvalues
% from an external source to the constructor. 
% The physical quantity we are interested in is the transmission which is 
% quadratic in the solutions of the Schreodinger equation. Therefore, the 
% conjugated poles must be included.
parameters.includeConjugatedPoles = true;
% Furthermore, define a range of frequencies (the expansion points), which
% represent the energy interval we ar interested in, 
parameters.expansionPoints = 1:0.1:6;
% the number of integration points on the background contour 
parameters.nPointsB = 64;
% and their minimal distance to the expansion points. In the case of an
% ellipse, the first value will be added to the mayor and the second value
% to the minor semiaxis. If the shape is selected to be a polygon, first 
% and second value correspond to shifts in x and y direction, respectively.
% If you want the outer contour to be a circle, you can provide a scalar
% value which is added to the radius.
parameters.minD = [0.4 2]; 
% As we expect to find poles below the real axis, we should add enough 
% space between expansion points and contour in horizontal direction and we
% must guess how many poles the contour will enclose.
parameters.nEigenvaluesExpected = 10;
% Eventually, we define some reference points which allow us to check if
% the result converges to the expected solution. 
parameters.referencePoints = parameters.expansionPoints(1:2:end);
% As we provide a custom interface, we must tell which quantities can be
% expanded and which of them are linear. A linear quantity is required to
% solve the eigenvalue problem and if RieszProjection knows the available
% quantities, it will suggest them if you call the method 'expand' without 
% arguments.
parameters.linearQ = {'RightBoundary'};
parameters.quadraticQ = {'Reflection', 'Transmission'};
% Now, an instance of the class RieszProjection can be constructed. Apart
% from the parameters we also hand a function handle to the constructor,
% which will be used to solve the linear systems, i.e., the Schroedinger
% equations, at the integration points and to derive the target quantities.
% It is defined at the bottom of this script.
rp = RieszProjection(@f, parameters);

%% get the contours and expand

% First, we must construct the contours. We want an elliptic shape and the
% trapezoidal rule, which is the default.
rp.getContours
% Now, we can plot the contour and, subsequently, expand the target quantity.
% As our target quanity is a quadratic form, but the property 'poles' is
% empty, the linear quantity 'RightBoundary' is expanded automatically in 
% order to solve the eigenvalue problem. The solutions of the linear
% systems at the integration points are subsequently used to evaluate the
% expansions.
f1 = rp.plot('contours');
rp.expand('Transmission');
% In order to see if the error decreases, we first plot the error. In
% contrast to other methods, we do not need the eigenvectors to estimate
% the error. As we are particularly interested in the expansion, we compare
% the result of the expansion to solutions of the linear system at the 
% reference points, i.e., solutions of the Schroedinger equation at the 
% real line.
f2 = rp.plot('error');

% Apparently, the error is indeed decreasing and the eigenvalues, displayed 
% in figure 1, which has been updated automatically, does not seem to be
% distributed randomly. There might be one pole that does not fit into the
% pattern but if it is a mere artefact, we can expect it to vanish once
% we further decrease the error. There are two ways to improve the
% accuracy: On the one hand side, we can simply increase the number of 
% integration points. On the other hand, since it is known that
% the convergence rate depends on the distance of the contour to the
% closest pole, we could first define a new shape. As the contour passes 
% very closely one of the poles, this might indeed be a good idea. We will 
% do both and compare the results.

%% refine the contour

% As we wish to reuse the existing solutions of the linear system, we
% double the number of integration points twice. Subsequently, we just
% expand the quantity we are interested in again. As the existing poles are
% used as initial guess for the eigenvalues, it might matter if the number
% of integration points is increased at once or in two subsequent steps.
% For simplicity, we will just do the first. As the eigenvalues have not
% been set by the user, they will be recalculated automatically if expand
% is called after the property 'nPointsB' has been increased. If you wand
% to improve the precision of the eigenvalues from an external source, you
% can call expand with the additional argument 'recalculateEvs'. E.g.,
% rp.expand('Transmission','recalculateEvs',true)
rp.nPointsB = 512;
% The figure containing the contours is updated automatically, yet, as the
% solutions of the linear systems might be computationally expensive, you
% must call the method 'expand' explicitly in order to update the
% expansions. Please note: If you call rp.expand in the command window,
% a list of all available quantities is displayed and you can select the
% one you are interested in with the mouse. The same holds for the plot
% function. 
rp.expand('Transmission')
% The plot of the error is updated automatically. If you want to visualize
% the result of the expansion, execute rp.plot in the command window and
% select the quanity you are interested in. 

%%
% At this point, the precision is limited by the precision of the
% eigenvalues which are used to evaluate the modal contributions. We can
% further improve the precision by adding small contours. From the
% convergence plot we can estimate the precision of the eigenvalues to be
% below 1e-4. Thus, circles with radii 1e-2 and centers at the estimated
% eigenvalues should enclose the true poles nicely. As a circular contour
% around a single pole converges very fast if its radius is small and if 
% the error of the eigenvalue is small compared to its radius, 8
% integration points will suffice.
rp.radius = 1e-2;
rp.nPoints = 8;
rp.useEigenvalues = false;
rp.expand('Transmission')
% Now the error of the expansion is limited by the error of the numerical
% integration.

%%
% You can address the errors of the individual contributions by adding the 
% index as a positional argument. If n modes contribute, the index of the 
% background contour is n+1 and the error of the total expansion can be 
% obtained with the index 0 or n+2. To convince ourselfs that the small
% circular contours converge well, we cann have, e.g. a look at the first
% two of them.
rp.plot('error',1:2)
% As no reference is available for the individual contributions, the
% solution with the most integration points is considered to be the best
% approximation. Consequently, the error of the soltion with 8 integration
% points can expected to be below 1e-8;

%% 
% For a later comparison to the deformed contour, we go back to the
% expansion error without small contours before we proceed.
rp.useEigenvalues = true;
rp.plot('error',0)
rp.expand('Transmission')

%% adjust and refine the contour

% In order to make a fair comparison, we reset the random number generator,
% construct a second instance, and proceed to the point where we made the
% refinement. 
rng(15);
rpc = RieszProjection(@f, parameters);
rpc.getContours;
fc1 = rpc.plot('contours','fignumber',11);
rpc.expand('Transmission');
fc2 = rpc.plot('error','fignumber',22);
% As the poles are positioned on a tilted line and one can expect more
% poles to the right, the minimal area enclosing ellipse cannot be easily
% defined in a way, that the distance to a pole outside the contour is
% maximized. Therefore, we select a circular shape and include the pole
% outside the contour to the selected poles. To avoid a change of the
% contour by accidance, we must confirm that we want to change the shape
% of the outer contour. If you do not set the hidden property
% 'newBackground', a dialogbox will pop up and ask for confirmation.
rpc.newBackground = true;
rpc.getContours('c');
rpc.selectedPoles = [rpc.selectedPoles length(rpc.poles)];
% As the new contour already required 128 new integration points we will
% set the property 'nPointsB' to 512 - 128 = 384
rpc.nPointsB = 384;
rpc.expand('Transmission');
% Here, the rate of convergence is indeed much faster and using the same
% number of integrations points results in a lower error.

%%
% As before, we can further increase the precision adding small contours.
% Here, we could decrease the radius even further. Yet, we need to have in
% mind that the conditioning of the linear problem becomes worse the
% closer we get to a pole.
rpc.radius = 1e-2;
rpc.nPoints = 8;
rpc.useEigenvalues = false;
rpc.expand('Transmission')


%% use eigenvalues from external source

% This time we assume that the eigenvalues are given and want to
% investigate the expansion over a larger energy interval.
% Fist, the quasi exact reference solution is obtained from the linearized
% problem. We request 30 eigenvalues centered at 6. 
ref = f([30,6],'Eigensystem');
evs = diag(ref.eigenvalues);
% Now we can equip our instances with the quasi exact eigenvalues. The
% outer contour is kept and all eigenvalues inside the contour are
% selected automatically.
rp.poles = evs;
rpc.poles = evs;
% We will start with a small number of integration points and refine the
% contours adaptively setting the property 'precision'. This property can
% be given as a scalar as below or as a vector with two entries. In the
% latter case, the first entry is used as an upper bound of the relative
% error and the second as an upper bound of the absolute value. For
% adaptive refinement, the minimum of absolute and relative error is used.
rpc.nPointsB = 60;
rp.nPointsB = 60;
rp.precision = 1e-4;
rpc.precision = 1e-4;
% In this usecase, we will compare the trapezoidal rule to the
% Gauss-Kronrod quadrature rule on a polygon.
rp.newBackground = true;
rp.getContours('p',true) 
% Here, the second argument refers to the use of the Gauss-Kronrod
% quadrature rule.

%% 
% Now, we are interested in the frequency range from 2 to 13.
rpc.newBackground = true;
rpc.expansionPoints = 2:0.1:15;
rpc.referencePoints = rp.expansionPoints(1:4:end);
% We want to use the same poles and the same expansion points for the
% polygon and we must add the small contours again.
rp.selectedPoles = rpc.selectedPoles;
rp.expansionPoints = 2:0.1:15;
rp.useEigenvalues = false;
rp.referencePoints = rp.expansionPoints(1:4:end);

%% expand

rp.expand('Transmission');
rpc.expand('Transmission');
% If we have not closed the error plots, they will have been updated
% (otherwise we can just replot them) and we see that in this particular
% example the target precision has been reached after the same number of
% linear systems has been solved as one has to take into account the 90
% solutions in the polygon case, that could not be reused after refinement.
% Nevertheless, the exponential convergence of the trapezoidal rule, given
% circular or elliptic contours, will in practice often beat higher order
% quadrature rules.

%%
% The minimum of relative and absolute error, which has been used for the
% adaptive refinement can be plotted as well.
rp.plot('error','errortype','m')

%% function used to evaluate the scattering problem
function out = f(varargin)
% For custom solvers you need a similar function that, depending on the
% input arguments, either solves a linear system (nargin==1) or applies a
% post process. In this particular case there is third option that returns
% the eigenspectrum of the linearized system. 
% The function is expected to return fields {{f11 f12 ...} {f21 ...} ...}
% if the input have been frequencies {[w11 ...] [w21 ...] ...} and if the
% input are fields {{f11 f21 ...} {f21 ...} ...} it is expected to return
% the corresponding target quantities {[q11 q12 ...] [q21 ...] ...}. If the
% quantity is quadratic in the fields the input is of the form:
% {{f111 f121 ...} {f112 f122 ...};{f211 f221 ...} {f212 f222 ...};...}. So
% the first index indicates the contour, the second the integration point
% on this contour, and the third distinguishes the normal field from the
% circfield. If the first input contains fields a second argument is
% expected that specifies the post process. 

persistent Ah Bh Ch jh s a;
if isempty(Ah)
    % define sparse matrices used as descrete representation of the problem
    n = 300; L = pi*sqrt(2)/2; h=2*L/(n+1); V0 = 10;
    ndx1 = [2:n+2 1:n+1 1:n+2]; ndx2 = [1:n+1 2:n+2 1:n+2]; a = 1;
    Ah = sparse(ndx1,ndx2,h/6*[ones(1,2*n+2) 2 4*ones(1,n) 2]);
    Bh = sparse([1 n+2],[1 n+2],[1 1]);
    Ch = sparse(ndx1,ndx2,1/h*[-1*ones(1,2*n+2) 1 2*ones(1,n) 1]) - V0*Ah;
    % define physical source as a function of the frequency
    jh = @(w) [2i*w*exp(-1i*w*pi/sqrt(2)); zeros(n+1,1)]; 
    s = @(w,x) a*exp(1i*w*x);
end
% Get quasi exact solution of the linearized problem (nc: [# center])
if nargin==2 && strcmpi('eigensystem',varargin{2})
    out = struct; nc = varargin{1}; if length(nc)==1, nc = [nc 0]; end
    E = speye(size(Ah)); O = sparse(size(Ah,1),size(Ah,1));
    A = [-1i*Bh, Ch; E,  O];  B = [Ah,  O;  O, E];
    [out.eigenvectors,out.eigenvalues]=eigs(A,B, nc(1), nc(2));
elseif iscell(varargin{1})
    out = cell(1,size(varargin{1},1+(nargin==1)));
    for it = 1:length(out)
        if nargin==1
            out{it} = solveScatteringProblems(varargin{1}{it});
        else
            out{it} = postProcess(varargin{1}(it,:));
        end
    end
end

    function out = solveScatteringProblems(ws)
        % Solve the inhomegeneous equation
        if isempty(ws), return; end
        out = cell(1,numel(ws));
        for it2=1:numel(ws)
            w = ws(it2);
            out{it2} = {(w^2*Ah+1i*w*Bh-Ch)\(a*jh(w)) w};
        end
    end


    function out = postProcess(fs)
        if contains(lower(varargin{end-1}),'rightboundary')
            out = cellfun(@(x)x{1}(end),fs{1});
        elseif contains(lower(varargin{end-1}),'transmission')
            % The probability that the particle has passed the quantum well
            f1 = cellfun(@(x)x{1}(end),fs{1});
            if length(fs)==2
                % Because of the relation f(-w) = conj(f(conj(w)) the 
                % problem is not solved at negative but at complex
                % conjugate frequencies and therefore f(-w) we obtained
                % taking the complex conjugate.
                f2 = conj(cellfun(@(x)x{1}(end),fs{2}));
            else
                % In the case of the reference solution i.e. if the 
                % frequencies are real.
                f2 = conj(f1);
            end
            out = f1.*f2;
        elseif contains(lower(varargin{end-1}),'reflection')
            % The probability that the particle is scattered back
            f1 = cellfun(@(x)x{1}(1),fs{1});
            if length(fs)==2
                % As in the case of transmission
                f2 = conj(cellfun(@(x)x{1}(1),fs{2}));
            else
                f2 = conj(f1);
            end
            w = cellfun(@(x)x{2},fs{1});
            out = (f1-s(w,-pi/sqrt(2))).*(f2-s(-w,-pi/sqrt(2)));
        end
    end
end
