edcSMOKE++
==========

EDC (Eddy Dissipation Concept) for OpenFOAM based on the [OpenSMOKE++ framework][1].
The edcSMOKE++ framework is available for the following OpenFOAM versions:
- OpenFOAM-10

If you use edcSMOKE++ for your publications, we kindly ask you to cite the following two papers:

> Parente, A., Malik, R.M., Contino, F., Cuoci, A., Dally, B., 
> Extension of the Eddy Dissipation Concept for turbulence/chemistry interactions to MILD combustion
> (2016) Fuel, 163, pp. 98-111, DOI: 10.1016/j.fuel.2015.09.020

> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms 
> (2015) Computer Physics Communications, 192, pp. 237-264, DOI: 10.1016/j.cpc.2015.02.014

Compulsory libraries
--------------------
- OpenSMOKE++ (already included in edcSMOKE++)
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- Boost C++ (http://www.boost.org/)

Optional libraries
------------------
- Intel OneAPI MKL (https://software.intel.com/en-us/intel-mkl)
- ODEPACK (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DVODE (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DASPK (http://www.engineering.ucsb.edu/~cse/software.html)
- Sundials (http://computation.llnl.gov/casc/sundials/main.html)
- MEBDF (http://wwwf.imperial.ac.uk/~jcash/IVP_software/readme.html)
- RADAU (http://www.unige.ch/~hairer/software.html)

Compilation
-----------
Two main different options are available to compile the code, according to the level of support for the solution of ODE systems:
1. Minimalist: no external, optional libraries are required. Only the native OpenSMOKE++ ODE solver can be used.
2. Minimalist + Intel OneAPI MKL: only the native OpenSMOKE++ ODE solver can be used, but linear algebra operations are managed by the Intel OneAPI MKL libraries
3. Complete: all the optional libraries are linked to the code, in order to have the possibility to work with different ODE solvers

### 1. Instructions to compile the Minimalist version
-------------------------------------------------
1. Copy the `mybashrc.minimalist` file as `mybashrc.local` by typing: `cp mybashrc.minimalist source mybashrc.local`
2. Adjust the paths to the compulsory libraries in the `mybashrc.local` file
3. Type: `source mybashrc.local`
4. Compile the kinetic pre-processor: from the `solver/chemkin2OpenSMOKE++PreProcessor` folder, type `wmake`
5. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder, type `wmake`
6. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`

### 2. Instructions to compile the Minimalist+MKL version
-----------------------------------------------------
1. Copy the `mybashrc.minimalist.mkl` file as `mybashrc.local` by typing: `cp mybashrc.minimalist.mkl source mybashrc.local`
2. Adjust the paths to the compulsory libraries in the `mybashrc.local` file
3. Adjust the paths to the Intel OneAPI MKL libraries in the `mybashrc.local` file
4. Type: `source mybashrc.local`
5. Compile the kinetic pre-processor: from the `solver/chemkin2OpenSMOKE++PreProcessor` folder, type `wmake`
6. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder, type `wmake`
7. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`

Preprocessing of CHEMKIN files
-----------------------------------------------------
In order to run a simulation with edcSMOKE++, a CHEMKIN mechanism (kinetics, thermodynamic and transport properties) has to be pre-processed using the `chemkin2OpenSMOKE++PreProcessor` utility (see above). 
Examples of mechanisms ready to be pre-processed are available in the `run/kinetic-mechanisms` folder. In particular, for in each mechanism folder you can find the three files corresponding to the CHEMKIN input (kinetics, thermodynamics, and transport properties) and an additional `input.dic` file, containing the instructions for the `chemkin2OpenSMOKE++PreProcessor`.
In order to pre-process a kinetic mechanisms, the operations to carry out are very simple. As an example, for `POLIMI_H2CO_1412` mechanism:
1. Go to the `run/kinetic-mechanisms/POLIMI_H2CO_1412`
2. Type `chemkin2OpenSMOKE++PreProcessor`
3. If everything works correctly, a `kinetics-POLIMI_H2CO_1412` folder will be created, including the preprocessed CHEMKIN files (in XML folder). This is the folder which has to be supplied to the `edcSMOKE++` solvers.

[1]: https://www.opensmokepp.polimi.it/
