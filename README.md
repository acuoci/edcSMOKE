edcSMOKE
========

EDC (Eddy Dissipation Concept) for OpenFOAM based on the [OpenSMOKE++ framework][1]
The edcSMOKE requires one of the following OpenFOAM versions:
- OpenFOAM-4.x
- OpenFOAM-2.4
- OpenFOAM-2.3
- OpenFOAM-2.2

If you use edcSMOKE for your publications, we kindly ask you to cite the following two papers:

> Parente, A., Malik, R.M., Contino, F., Cuoci, A., Dally, B., 
> Extension of the Eddy Dissipation Concept for turbulence/chemistry interactions to MILD combustion
> (2016) Fuel, 163, pp. 98-111, DOI: 10.1016/j.fuel.2015.09.020

> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms 
> (2015) Computer Physics Communications, 192, pp. 237-264, DOI: 10.1016/j.cpc.2015.02.014

Compulsory libraries
--------------------
- OpenSMOKE++ (already included in edcSMOKE)
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- RapidXML (http://rapidxml.sourceforge.net/)
- Boost C++ (http://www.boost.org/)

Optional libraries
------------------
- Intel MKL (https://software.intel.com/en-us/intel-mkl)
- ODEPACK (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DVODE (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DASPK (http://www.engineering.ucsb.edu/~cse/software.html)
- Sundials (http://computation.llnl.gov/casc/sundials/main.html)
- MEBDF (http://wwwf.imperial.ac.uk/~jcash/IVP_software/readme.html)
- RADAU (http://www.unige.ch/~hairer/software.html)

Compilation
-----------
Three different options are available to compile the code, according to the level of support for the solution of ODE systems
1. Minimalist: no external, optional libraries are required. Only the native OpenSMOKE++ ODE solver can be used.
2. Minimalist + Intel MKL: only the native OpenSMOKE++ ODE solver can be used, but linear algebra operations are managed by the Intel MKL libraries
3. Complete: all the optional libraries are linked to the code, in order to have the possibility to work with different ODE solvers

1. Instructions to compile the Minimalist version
-------------------------------------------------
1. Open the `mybashrc.minimalist` and adjust the paths to the compulsory external libraries (in particular choose the OpenFOAM version you are working with)
2. Type: `source mybashrc.minimalist`
3. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder type `wmake`
4. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`

2. Instructions to compile the Minimalist+MKL version
-----------------------------------------------------
1. Open the `mybashrc.minimalist.mkl` and adjust the paths to the compulsory external libraries and the paths to the Intel MKL library (in particular choose the OpenFOAM version you are working with)
2. Type: `source mybashrc.minimalist.mkl`
3. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder type `wmake`
4. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`

3. Instructions to compile the Complete version
-----------------------------------------------------
1. Open the `mybashrc.complete` and adjust the paths to the compulsory external libraries and the Intel MKL library (in particular choose the OpenFOAM version you are working with). You can choose the additional external libraries you want to add to edcSMOKE, by modifying the `EXTERNAL_ODE_SOLVERS` variable: in particular `1` means that the support is requested, while `0` means that no support is requested. Obviously, for each requested library, you need to provide the correct path.
2. Type: `source mybashrc.complete`
3. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder type `wmake`
4. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`

[1]: https://www.opensmokepp.polimi.it/
