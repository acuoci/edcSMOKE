edcSMOKE
========

EDC (Eddy Dissipation Concept) for OpenFOAM based on the OpenSMOKE++ Library

Compulsory libraries
--------------------
Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
RapidXML (http://rapidxml.sourceforge.net/)
Boost C++ (http://www.boost.org/)

Optional libraries
------------------
Intel MKL (https://software.intel.com/en-us/intel-mkl)
ODEPACK (http://computation.llnl.gov/casc/odepack/odepack_home.html)
DVODE (http://computation.llnl.gov/casc/odepack/odepack_home.html)
DASPK (http://www.engineering.ucsb.edu/~cse/software.html)
Sundials (http://computation.llnl.gov/casc/sundials/main.html)
MEBDF (http://wwwf.imperial.ac.uk/~jcash/IVP_software/readme.html)
RADAU (http://www.unige.ch/~hairer/software.html)

Compilation
-----------
Three different options are available to compile the code, according to the level of support for the solution of ODE systems
1. Minimalist: no external, optional libraries are required. Only the native OpenSMOKE++ ODE solver can be used.
2. Minimalist + Intel MKL: only the native OpenSMOKE++ ODE solver can be used, but linear algebra operations are managed by the Intel MKL libraries
3. Complete: all the optional libraries are linked to the code, in order to have the possibility to work with different ODE solvers

1. Instructions to compile the Minimalist version
----------------------------------------------
a. Open the `mybashrc.minimalist` and adjust the paths to the compulsory external libraries
b. Type: `source mybashrc.minimalist`
c. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder type `wmake`
d. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`

2. Instructions to compile the Minimalist+MKL version
-----------------------------------------------------
a. Open the `mybashrc.minimalist.mkl` and adjust the paths to the compulsory external libraries and the paths to the Intel MKL library
b. Type: `source mybashrc.minimalist.mkl`
c. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder type `wmake`
d. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`

3. Instructions to compile the Complete version
-----------------------------------------------------
a. Open the `mybashrc.complete` and adjust the paths to the compulsory external libraries and the Intel MKL library. You can choose the additional external libraries you want to add to edcSMOKE, by modifying the `EXTERNAL_ODE_SOLVERS` variable: in particular `1` means that the support is requested, while `0` means that no support is requested. Obviously, for each requested library, you need to provide the correct path.
b. Type: `source mybashrc.complete`
c. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder type `wmake`
d. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`


