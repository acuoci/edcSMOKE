edcSMOKE++
==========

EDC (Eddy Dissipation Concept) for OpenFOAM based on the [OpenSMOKE++ framework][1].
The edcSMOKE++ framework is available for the following OpenFOAM versions:
- OpenFOAM 10

If you use edcSMOKE++ for your publications, we kindly ask you to cite the following two papers:

> Li, Z., Malik, M.R., Cuoci, A., Parente, A.,
> edcSMOKE: A new combustion solver for stiff chemistry based on OpenFOAM,
> (2017) AIP Conf. Proc., 1863(1), p. 180004, DOI: 10.1063/1.4992364

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

Recommended libraries
---------------------
- Intel OneAPI MKL (https://software.intel.com/en-us/intel-mkl)

Optional libraries
------------------
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
2. Minimalist + Intel OneAPI MKL (recommended): only the native OpenSMOKE++ ODE solver can be used, but linear algebra operations are managed by the Intel OneAPI MKL libraries
3. Complete: all the optional libraries are linked to the code, in order to have the possibility to work with different ODE solvers

### 1. Instructions to compile the Minimalist version
-----------------------------------------------------
1. Copy the `mybashrc.minimalist` file as `mybashrc.local` by typing: `cp mybashrc.minimalist mybashrc.local`
2. Adjust the paths to the compulsory libraries in the `mybashrc.local` file
3. Type: `source mybashrc.local`
4. Compile the kinetic pre-processor: from the `solver/chemkin2OpenSMOKE++PreProcessor` folder, type `wmake`
5. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder, type `wmake`
6. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`

### 2. Instructions to compile the Minimalist+MKL version
---------------------------------------------------------
1. Copy the `mybashrc.minimalist.mkl` file as `mybashrc.local` by typing: `cp mybashrc.minimalist mybashrc.local`
2. Adjust the paths to the compulsory libraries in the `mybashrc.local` file
3. Adjust the paths to the Intel OneAPI MKL libraries in the `mybashrc.local` file
4. Type: `source mybashrc.local`
5. Compile the kinetic pre-processor: from the `solver/chemkin2OpenSMOKE++PreProcessor` folder, type `wmake`
6. Compile the steady-state solver: from the `solver/edcSimpleSMOKE` folder, type `wmake`
7. Compile the unsteady solver: from the `solver/edcPimpleSMOKE` folder type `wmake`

Preprocessing of CHEMKIN files
------------------------------
In order to run a simulation with edcSMOKE++, a CHEMKIN mechanism (kinetics, thermodynamic and transport properties) has to be pre-processed using the `chemkin2OpenSMOKE++PreProcessor` utility (which is part of edcSMOKE++). 
Examples of mechanisms ready to be pre-processed are available in the `run/kineticMechanisms` folder. In particular, in each mechanism folder three files are present, corresponding to the conventional CHEMKIN input files (kinetics, thermodynamics, and transport properties). An additional `input.dic` file, containing the instructions for the `chemkin2OpenSMOKE++PreProcessor`, is also included.
In order to pre-process a kinetic mechanisms, the operations to carry out are very simple. As an example, for the `CRECK_2003_SYNGAS` mechanism:
1. Go to the `run/kineticMechanisms/CRECK_2003_SYNGAS`
2. Type `chemkin2OpenSMOKE++PreProcessor`
3. If everything works correctly, a `kinetics` folder will be created, including the preprocessed CHEMKIN files (in XML format). This is the folder which has to be supplied to the edcSMOKE++ solvers.

# Tutorials
## 1. Turbulent jet flame fed with H2
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
The first step consists in running the unsteady simulation using the simplest turbulent combustion model available in edcSMOKE++, which is the Eddy Dissipation (ED) model. A 1-step global kinetic mechanism is adopted. It is always a good idea to start a new simulation from scratch using the ED model to speedup the simulation to reach a steady-state solution. Then, this solution will become the starting point for the application of a more refined turbulent combustion model, such as the Eddy Dissipation Concept (EDC).
1. Go to the `run/tutorials/JetFlameH2/unsteady/01-ED-1step` folder.
2. Build the mesh using the `blockMesh` utility.
3. Prepare the kinetic mechanism in the OpenSMOKE++ format. The kinetic mechanism is provided in the `kinetics` folder, which includes the usual CHEMKIN files for chemical reactions (`kinetics.cki`), thermodynamics (`therm.dat`), and transport properties (`tran.trc`). From the `kinetics` folder, run the `chemkin2OpenSMOKE++PreProcessor` utility to transform the kinetic mechanism from the CHEMKIN format to the OpenSMOKE++ format. If everything works well, a new subfolder with name `kinetics` will be created in the main `kinetics` folder.
4. Prepare the kinetic mechanism in the OpenFOAM format. Unfortunately, the current version of edcSMOKE++ requires a double version of the kinetic mechanism, which must be included in the `chemkin` folder. This folder includes three files: the thermodynamic file (`therm.dat`) which is exactly the same available in the `kinetics` folder; the kinetic mechanism file (`kinetics.cki`) already provided in the `kinetics` folder, but without chemical reactions (i.e., it must include the list of species only, exactly in the same order); a file called `transportProperties` which is always the same, independently of the kinetic mechanism adopted. In order to process these three files and transform them into the OpenFOAM format, the `chemkinToFoam` utility must be used. From the main folder, run the following command: `chemkinToFoam chemkin/kinetics.cki chemkin/therm.dat chemkin/transportProperties constant/reactions constant/speciesThermo`. If everything works properly, two files with names `speciesThermo` and `reactions` will be created in the `constant` folder.
5. We are now ready to run the simulation by means of the `edcPimpleSMOKE` solver. From the main folder, simply type: `edcPimpleSMOKE`. In this specific example, about 0.6 s of simulation time are needed to reach a steady-state solution.

#### 2. EDC (Eddy Dissipation Concept) simulation
We are now ready to run the EDC simulation with a detailed kinetic mechanism (32 species and 173 rections) including also NOx chemistry. The starting point is obviously the steady-state solution obtained using the ED model (see step 1 above).    
1. Go to the `run/tutorials/JetFlameH2/unsteady/02-EDC-Polimi-H2NOX` folder.
2. Build the mesh using the `blockMesh` utility.
3. Prepare the kinetic mechanism in the OpenSMOKE++ format (please follow exactly the same procedure described for the ED model).
4. Prepare the kinetic mechanism in the OpenFOAM format (please follow exactly the same procedure described for the ED model).
5. Prepare the initial solution. Copy the steady-state solution from the ED simulation (which we are assuming is at time 0.6 s): `cp -r ../01-ED-1step/0.6 0`. Remove the `uniform` folder: `rm -r 0/uniform`. Copy the `Ydefault` file: `cp ../01-ED-1step/0/Ydefault 0/`.
6. We are now ready to run the simulation by means of the `edcPimpleSMOKE` solver. From the main folder, simply type: `edcPimpleSMOKE`. Since we start already from a steady-state simulation, the time needed to reach the new steady state conditions is smaller. In this specific example, about 0.2 s of simulation time are sufficient.

#### 3. EDC (Eddy Dissipation Concept) steady-state simulation
The EDC model is computationally very expensive. Thus, if we are interested in the steady-state solution only, it is more convenient to solve the trasnport equations of mass, momentum, species and energy directly in their steady-state formulation, using the `edcSimpleSMOKE` steady-state solver. However, it is important that the `edcSimpleSMOKE` solver is applied only if a have a reasonable first-guess solution. 
1. Go to the `run/tutorials/JetFlameH2/steady/02-EDC-Polimi-H2NOX` folder.
2. Build the mesh using the `blockMesh` utility.
3. Prepare the kinetic mechanism in the OpenSMOKE++ format (please follow exactly the same procedure described for the ED model).
4. Prepare the kinetic mechanism in the OpenFOAM format (please follow exactly the same procedure described for the ED model).
5. Prepare the initial solution. Copy the latest solution from the unsteady EDC simulation (which we are assuming is at time 0.2 s): `cp -r ../../unsteady/02-EDC-Polimi-H2NOX/0.2/ 0`. Remove the `uniform` folder: `rm -r 0/uniform`.
6. We are now ready to run the simulation by means of the `edcSimpleSMOKE` solver. From the main folder, simply type: `edcSimpleSMOKE`.
   
   
[1]: https://www.opensmokepp.polimi.it/
