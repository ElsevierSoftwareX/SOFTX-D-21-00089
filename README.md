# RPExpand

The software RPExpand is designed to give an easy and efficient access to 
modal expansions of electromagnetic fields in optical resonators based on 
Riesz projections [1]. Given the solutions of scattering problems,
individual contributions of eigenmodes to the target quantity are made
accessible. These measurables do not need to be linear in the electric or
magnetic field but can be based on quadratic forms [2] or related to the 
far field [3], which enables for direct comparisons to experimental data.
Results can be directly visualized with a built-in plot function and a 
comprehensive control of the error is available. The expansion can be
compared to solutions of the scattering problem at selected real-valued 
frequencies and to results with reduced numbers of integration points.
If only few modes are present in the spectral region of interest, the
implementation of a contour integral method [4] for nonlinear eigenvalue 
problems makes it possible to extract all modal informations from one
single contour without the need to solve the eigenvalue problem beforehand.
An advantage of contour integral methods is parallelizability. This has
been exploited providing an interface to the finite element method (FEM) 
solver [JCMsuite](https://jcmwave.com/). The interface handles the parallel 
submission of jobs keeping the interaction with the solver as easy as 
possible while allowing advanced users to exploit its full capabilities. 

## Installation

Get a copy of the two class folders and the directory containing post
processes of JCMsuite used to derive the target quantities from the 
electric field:

- @Scattering (interface to JCMsuite)
- @RieszProjection 
- postprocesses

and add their parent directory to your matlab path or copy them to the
directory from which your scripts are executed.

## Documentation

Comments explaining the usage of classes and functions can be displayed
using the help function. E.g., `help RieszProjection` will print a short
description and a list of properties and methods, which are linked to more
detailed instructions on the corresponding items. Detailed examples are
provided to make you familiar with the functionalities and to serve as
a guideline for setting up your own projects. Additionally, the following 
sections of this README summarize the prerequisites needed in order to use
this code with the FEM solver JCMsuite. 

The class 'RieszProjection' can be used independently of JCMsuite. An 
example is given in the script 'quantumExample.m'. This example runs very
fast, does not require an installation of JCMsuite and covers the most 
important features of the code. Therefore it is recommanded to run this 
script to get familiar with RPExpand.

## Requirements

RPExpand has been tested with MATLAB R2018b under Linux but should be
independent of the platform. The class 'RieszProjection' is independent of 
JCMsuite. You can provide an interface to any software capable to solve the 
scattering problem at complex frequencies. This interface must have the form 
of a callable object. A rather complex example is the class 'Scattering'
which is the interface to JCMsuite whose version should be 4.4.0 or higher.
A much simpler example is provided in the file 'quantumExample.m'. There, 
the interface is a function handle. Input and output of such a custom interface 
are described in more detail in the section 'Notes on a usage independent
of JCMsuite' at the end of this README.

### Install JCMsuite

The current release can be downloaded following this 
[link](https://installation.jcmwave.com/).
Please follow the installation and activation instructions provided there.
A 14 day trial licence is available free of charge. It does 
not allow for the parallel computation of different scattering problems, 
however, multi threading is supported. Nevertheless, you need to start a 
daemon.

### Start a daemon

As a prerequisite for the usage of RPExpand the matlab interface to
JCMsuite must be set up and and a daemon has to be started. A minimal 
example given a personal laptop or desktop machine is:

```matlab
% Set up thirdparty support
jcm_root = '/path/to/local/installation/JCMsuite.x.x.x'
addpath(fullfile(jcm_root,'ThirdPartySupport','Matlab'))

% shutdown a possibly running daemon
jcmwave_daemon_shutdown();

% register a new resource
options = struct(...
                 'Hostname', 'localhost', ...
                 'Multiplicity', 1, ...
                 'NThreads', 2, ...
                 );
jcmwave_daemon_add_workstation(options);
```

If you have access to a remote machine via ssh, you can provide the 
corresponding hostname, user name, etc. For information about additional 
parameters or the usage of a queue or a cluster, please read the
[daemon command reference](https://www.docs.jcmwave.com/JCMsuite/html/MatlabInterface/6b7e17ce176f4ffde673cced4f12eeee.html).

If you have access to a license which allows for the parallel submission of
different jobs, you can increase the multiplicity to a value larger than one
depending on the available kernels of the computer. The parameter 'NThreads'
refers to the number of threads per job.

### Custom project requirements

This section gives a brief overview on the creation of a directory hosting 
a scattering project of JCMsuite compatible with the present code. It must
contain the following files: 
- [project.jcmpt](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/db9062933554f66e7fb21c46d53fcca2.html)
- [sources.jcmt](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/f3e666a5067147d3cd45b67773bb77ae.html)
- [materials.jcm(t)](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/3df274a2924c89630ff2393cc22b686e.html)
- [layout.jcm(t)](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/b61236968b3822be5ffbfee6564f23da.html)

The trailing 't' at the file extentions marks template files, which are used
for parameter substitution. Whereas in the last two cases it is optional, the
first to files must contain predefined parameters. A minimal project file is: 

```
Project {
    Electromagnetics {
        TimeHarmonic {
            Scattering {
                PML {
                    %(pml)s
                }
                FieldComponents = Electric
                Accuracy {
                    FiniteElementDegree = %(finiteElementDegree)e
                }
            }
        }
    }
}
```
The parameters 'pml' and 'finiteElementDegree' will be set by the program. The 
former as the perfectly matched layers (PML) must be the same for each 
integration point. The latter, as this way it can be conveniently set in the 
script and as some post processes require this information in order to 
determine the integration order. 
If there exists a file 'pml.log' in the project directory, the PML will be 
based on the parameters in this file, otherwise it will be created with the 
first scattering problem.
As the scattering problem has to be solved for many different frequencies,
the source file must contain the definition `Omega = %(omega)e`.

Custom parameters can be passed to the FEM solver via the property 'keys'
of Scattering. Please refer to the 
[parameter reference](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/index.html)
and the documentation of the
[matlab interface](https://www.docs.jcmwave.com/JCMsuite/html/MatlabInterface/index.html)
of JCMsuite. Further details are provided in the examples.

## Expand a custom quantity

Currently the following quantities are available for expansion: 

- the total field
- dipole emission and normalized decay rate of a point source 
- electromagnetic field energy flux
- electric field energy 
- power radiated to the far field and photon collection efficiency
- radiation pattern

The electromagnetic field energy flux density is currently integrated over
the boundaries of the computational domain, i.e., the result gives the dipole
emission. Yet, it can be easily adapted to expand, e.g., the energy flux between 
two neighboring domains. The radiated power together with the dipole emission 
is used to get the dipole power collection efficiency. The imaginary part of 
the expansion of the electric field energy can be used to investigate the 
absorption of a given domain.

If you want to expand a quantity which is currently not available, you can 
either contact us or add a new post process to the existing ones. Doing the 
latter, the following paragraph will give you some guidance. 

For an efficient integration, not the total field is integrated along the
contours but the target quantities, which often are scalars. An array must be 
flattened to the size nx1 for integration and can be reshaped to its
original size afterwards. This has been done for the expansion of the
radiation pattern. The name of the quantity must match the name of the 
jcmpt-file located in the directrory 'postprocesses'. It is defined in JCMsuite
with the primitive 'OutputQuantity' or, in cases where the integrand is a 
python expression, 'IntegralName'. If a scalar is returned, the creation of the 
new jcmpt-file is all you need to do. If the result needs further post 
processing this can be done in the methods 'updateDerivedQuantities' and 
'getReferenceSolution', below the comment `% Extract numeric results`. Please be 
aware that you can only expand holomorphic expressions. For the visualization of
your result you can add a new private method which is called by 'plot'. 

## Notes on a usage independent of JCMsuite

If you want to use your own software to solve the linear systems at the abszissae, 
you can still use the class 'RieszProjection'. The integration is based on 
quadrature rules of the form: 

![\begin{equation*}
I(\omega_0) = \sum_i c_i(\omega_0) f(\omega_i)
\end{equation*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Bequation%2A%7D%0AI%28%5Comega_0%29+%3D+%5Csum_i+c_i%28%5Comega_0%29+f%28%5Comega_i%29%0A%5Cend%7Bequation%2A%7D%0A)

Using the method 'getContours' you generate the weights ![c_i(\omega_0)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+c_i%28%5Comega_0%29)
and the frequencies ![omega](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Comega_i%0A). 
All you have to provide are the function values 
![f(\omega_i)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+f%28%5Comega_i%29), 
e.g., solutions of the second order Maxwell's equation:

![\begin{equation*}
\nabla \times \mu^{-1}\nabla\times\bf{E}(\bf{r},\omega)-\omega^2\epsilon(\omega)\bf{E}(\bf{r},\omega)=i\omega\bf{J}(\bf{r}).
\end{equation*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Bequation%2A%7D%0A%5Cnabla+%5Ctimes+%5Cmu%5E%7B-1%7D%5Cnabla%5Ctimes%5Cbf%7BE%7D%28%5Cbf%7Br%7D%2C%5Comega%29-%5Comega%5E2%5Cepsilon%28%5Comega%29%5Cbf%7BE%7D%28%5Cbf%7Br%7D%2C%5Comega%29%3Di%5Comega%5Cbf%7BJ%7D%28%5Cbf%7Br%7D%29.%0A%5Cend%7Bequation%2A%7D%0A)

Usually, it is not necessary to integrate the total field. RPExpand is
designed to evaluate the target quantities at the integration points and, 
subsequently, to integrate them seperately. As quantities based on quadratic 
forms require solutions of the two independent function values 
![f(\omega_i)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+f%28%5Comega_i%29)
and ![f(-\omega_i) = f^*(\omega_i^*)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+f%28-%5Comega_i%29+%3D+f%5E%2A%28%5Comega_i%5E%2A%29), the solution of the 
partial differential equation (the scattering problem) and the postprocess to 
extract the quantity of interest must be done in two independent steps. 

The interface to your favourite solver must be defined as a callable meeting
the following criteria:

The frequencies which discretize the contours are saved in a cell array 
'contours' with size (1, n) where n is the number of contours. Each element is a 
complex double array of shape (k,m), which corresponds to a contour. Here, k is the 
number of nodes per subinterval. Unless using higher order quadrature methods, there 
is only one interval, i.e., m = 1 and k is the total number of integration points. 
In a first step, this cell array is passed to your function, which is expected to
return a cell array of size (1, n) or (n, 1) whose elements are objects of size (1, k\*m), 
e.g., a cell array containing numeric arrays, each of them representing the 
solution of a partial differential equation. In a second step, each element of 
the returned cell vector is passed to your custom function a second time with
some additional arguments: `v = f(sc_results,quantity,keys)`. The additional
arguments 'quantity' and 'keys' define the quantity which has to be exctracted from 
the results returned previously and additional meta data (e.g., about the shape), 
respectively. The quantity is defined as a 'char' vector (e.g., 'ElectricFieldEnergy') 
and the last argument is a scalar struct. The returned value 'v' is expected to be a 
cell containing numeric arrays of size (d, k\*m), where d is a custom number, which in many 
cases will be one, but can be larger as the examples for the expansion of the 
radiation pattern shows. In the case of a quantity quadratic in the solutions of the
scattering problems, the first input is of size (n, 2) and the second column
contains the solutions of the conjugated scattering problems. Otherwise it is of 
shape (n, 1).
An example for a function meeting the described criteria is given in the file 
'quantumExample.m' which implements the quantum transition problem in 1D.

In a last step, you have to edit the constant properties 'linearQ' and 'quadraticQ'.
Insert all quantities available for expansion. This way, depending on whether
they are linear or quadratic in ![f(\omega_i)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+f%28%5Comega_i%29), either one or two function values are 
passed and a call of the method 'expand' without argument will suggest the 
available expansions. 

## References

[[1]](https://doi.org/10.1103/PhysRevA.98.043806) L. Zschiedrich et al.,
*Phys. Rev. A*, vol. 98, p. 043806 (2018)

[[2]](https://doi.org/10.1103/PhysRevB.100.155406) F. Binkowski et al., 
*Phys. Rev. B*, vol. 100, p. 155406 (2019)

[[3]](https://doi.org/10.1103/PhysRevB.102.035432) F. Binkowski et al., 
*Phys. Rev. B*, vol. 102, p. 035432 (2020)

[[4]](https://doi.org/10.1016/j.jcp.2020.109678) F. Binkowski et al., 
*J. Comput. Phys.*, vol. 419, p. 109678 (2020)
