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

# Tutorials
## 1. Turbulent jet flame fed with H2 (unsteady simulation)
In the first tutorial we are going to simulate a turbulent jet flame fed with pure H2, burning in regular air. The flame is unconfined. Since the fuel and oxidizer nozzles have a circular section, a 2D axysymmetric simulation is carried out.

The tutorial is available in the [run/tutorials/JetFlameH2/unsteady/](run/tutorials/JetFlameH2/unsteady/) folder.

### Flame description

* Fuel composition (molar basis): 100% H2

* Nozzle inner diameter: Di = 3.75 mm

* Nozzle thickness: s = 0.55 mm

* Outer diameter: De =  150 mm

* Fuel inlet velocity: Ujet = 60 m/s

* Coflow velocity: Ucof = 2 m/s 

### Detailed instructions
#### 1. ED (Eddy Dissipation) simulation
The first step consists in running the unsteady simulation using the simplest turbulent combustion model available in edcSMOKE++, which is the Eddy Dissipation (ED) model. A 1-step global kinetic mechanism is adopted. It is always a good idea to start a new simulation from scratch using the ED model to speedup the simulation to reach a steady-state solution. This solution will become the starting point for the application of a more refined turbulent combustion model, such as the Eddy Dissipation Concept (EDC).
1. Go to the `run/tutorials/JetFlameH2/unsteady/01-ED-1step` folder.
2. Build the mesh using the `blockMesh` utility.
3. Prepare the kinetic mechanism in the OpenSMOKE++ format. The kinetic mechanism is provided in the `kinetics` folder, which includes the usual CHEMKIN files for chemical reactions (`kinetics.cki`), thermodynamics (`therm.dat`), and transport properties (`tran.trc`). From the kinetics folder, run the `chemkin2OpenSMOKE++PreProcessor` utility to transform the kinetic mechanism from the CHEMKIN format to the OpenSMOKE++ format. If everything works well, a new subfolder with name `kinetics` will be created in the main `kinetics` folder.
4. Prepare the kinetic mechanism in the OpenFOAM format. Unfortunately, the current version of edcSMOKE++ requires a double version of the kinetic mechanism, which must be included in the `chemkin` folder. This folder includes three file: the thermodynamic file which is exactly the same available in the `kinetics` folder (`therm.dat`); the kinetic mechanism file already provided in the `kinetics` folder, but without chemical reactions (i.e., it must include the list of species only, exactly in the same order); a file called `transportProperties` which is always the same, independently of the kinetic mechanism adopted. In order to process these three files and transform them in the OpenFOAM format, the `chemkinToFoam` utility must be used. From the main folder of the case, run the following command: `chemkinToFoam chemkin/kinetics.cki chemkin/therm.dat chemkin/transportProperties constant/reactions constant/speciesThermo`. If everything works properly, a `speciesThermo` and a `reactions` file will be created in the `constant` folder.
5. We are now ready to run the simulation by means of the `edcPimpleSMOKE` solver. From the main folder, simply type: `edcPimpleSMOKE`. In this specific example, about 0.6 s of simulation time are needed to reach a steady-state solution.
   

[1]: https://www.opensmokepp.polimi.it/
