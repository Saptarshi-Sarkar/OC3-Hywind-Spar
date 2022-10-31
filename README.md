# 1. Introduction
This package provides a MATLAB-based simulation tool to analyse the dynamics of a spar-type floating offshore wind turbine using Kane's method. The model is set up with the 5 MW OC3 Hywind Spar type FOWT. However, the codes are generic and will work for wind turbines with different ratings. For detail on the theory and derivation refer to [1].
The inputs are formatted as in [OpenFAST](https://github.com/OpenFAST/openfast). All structural and aerodynamic input data is read from OpenFAST input files. However, inputs related to the wind field, wave field and control systems are initialized separately as described in section 3.  

The mooring cables can be modelled using the [MoorDyn](https://www.nrel.gov/wind/nwtc/moordyn.html) or the [OpenMoor](https://github.com/chen-lin/openmoor) program. When using [OpenMoor](https://github.com/chen-lin/openmoor) the users have the additional capability of modelling wave-current interaction. For details on the theory of wave-current interaction refer to [2].
# 2. How to use
Three steps to run the code. 
## Step 1 - Create wind field using `TurbSim`
Use NREL's [TurbSim](https://www.nrel.gov/wind/nwtc/turbsim.html) tool to generate the desired 3D wind field. Please ensure that the rotor-swept area is entirely enclosed. Once generated, place the `*.sum` and `*.wnd` files in the `Wind_Dataset` folder. 

## Step 2 - Create a wave field using the `CreateWaveDataset.m` script
Use the MATLAB script `CreateWaveDataset.m` to create the wave velocity field. The created wave field in placed in the `Wave_Dataset` folder by the script. Use the script to create **regular wave** (sinusoidal wave) field or use the **PiersonMoskowitz** or **JONSWAP** spectrum. You can also include wave-current interaction in the generated wave field. For details on the wave-current interaction model refer to [2].

The `Wave` class has been created by [Lin Chen](https://github.com/chen-lin). You can download [OpenMoor](https://github.com/chen-lin/openmoor) from his page as well.

## Step 3 - Run the main simulation using the `OffshoreWindTurbine.m` script
The MATLAB script `OffshoreWindTurbine.m` is used to run the simulations. The script is described briefly below.
The code starts with some initial housekeeping.
```
clc; clear -globals;
clear SystemMatrices BaselineControllers;
addpath('./Wind_Dataset','./Wave_Dataset','./src') 
pause(.5);
global GRAVACC FLUIDDENSITY;      
GRAVACC = 9.80655; FLUIDDENSITY = 1025;
```
``MoorDyn`` of ``OpenMoor`` needs to be unloaded if it's loaded.
```
if libisloaded('MoorDyn') || libisloaded('MoorApiwin64')
    if libisloaded('MoorApiwin64')
    calllib('MoorApiwin64','finish'); unloadlibrary('MoorApiwin64');
    else
    calllib('MoorDyn','LinesClose');      unloadlibrary MoorDyn;  
    end                             % unload library (never forget to do this!)end
end
```
Read FAST input files
```
% Read FAST input data files
[Airfoils, Geometry] = ReadWindTurbineAeroDataInterp('rad');
[Blade, Tower]       = ReadWindTurbineStructuralData();
[ElastoDyn]          = ReadElastoDyn();
[Servo]              = ServoDyn();
```
Create the modal properties required for the dynamic analysis. Also, select the mooring cable model and select if wave loads are to be calculated or not.
```
% Create tower and blade structures
Twr = CreateTwr(Tower,ElastoDyn); Bld = CreateBld(ElastoDyn,Geometry,Blade); Platform = CreatePlatform();
Platform.Mooring     = 1;    % 1 for Moordyn, 2 for OpenMoor
Platform.WaveLoads   = 1;    % 1 to calculate wave load on spar using Morisson's equation, 0 to not.
```
If wave loads are included, place the correct file name (name of the ``.mat`` file created using the `CreateWaveDataset.m` script in the below line
```
    waveopt.wave_file = 'WaveHs3_Tp11_Dir90_NoCur0_Tf50';
```

String array that contains the list of all available DOFs
```
% DOFs available
DOFsStr = {'Sg','Sw','Hv','R','P','Y','TFA1','TSS1','TFA2','TSS2','NacYaw','GeAz','DrTr','B1F1','B1E1','B1F2','B2F1','B2E1','B2F2','B3F1','B3E1','B3F2'};
```
If you want to turn off certain DOFs, include them in the string array. An example showing how to turn off all blade second flapwise modes ``B1F2``, ``B2F2`` and ``B3F2``
```
TurnedOffDOFsStr = {'B1F2','B2F2','B2F2'};
```
Next include the correct filename name of the TurbSim generated wind file in the argument of the `readBLgrid.m` function
```
% Get the TurbSim-generated wind field and grid
[velocity, Wind.y, Wind.z, Wind.nz, Wind.ny, Wind.dz,...
          Wind.dy, Wind.dt, Wind.zHub, Wind.z1, Wind.SummVars] = readBLgrid('1ETMC');
```
For steady wind field, uncomment the next line
```
% [Wind.t_TurbSim, velocity] = SteadyWind(Wind.y, Wind.z, 11.4, 0.2);
```
Turn on/off **Pitt and peters** correction for the skewed wake, and **Blade pitch control**. The default settings are given below. In the default setting the **blade pitch control** must all be turned on.
```
WindNom.PittandPeters = true(0);  % Place 0 in the argument for FALSE or 1 for TRUE
WindNom.PitchControl  = true(1);  % Place 0 in the argument for FALSE or 1 for TRUE
WindNom.AeroElastic   = true(1);  % 1 to include aeroelastic effect in inplane direction, 0 to ignore 
```
Next, select start time ``t0``, final time ``tf`` and time step ``deltat``
```
t0  = 0;
tf  = 1;
deltat = 0.0125;
```
This concludes the set-up of the `OffshoreWindTurbine.m` script. The script runs using ``ODE 4`` RK $4^{th}$ order method and plots the results.
## Step 4* - Build  `./src/NominalSystemMatrix_mex.mex64` 
*This is an optional step. 
This package uses MATLAB mex files to accelerate the simulation time. The `NominalSystemMatrix.m` function is identified as the bottleneck function. Therefore, it has been built into a mex (MATLAB executable) function. The mex function is already provided in the `src` folder.  It should work normally without the need for rebuilding it. However, in case your MATLAB version does not support the pre-build mex file follow the following steps
1. In the `OffshoreWindTurbine.m` script change the final time `tf` to 0.1 sec.
2. In the `SystemMatrices.m` function change the function call in line 109 from 
``` 
[IM_nom, f_nom, Controls] = NominalSystemMatrix_mex(q_Nom, Controls, ElastoDyn, Airfoils, Twr, Bld, Platform, WindNom, mooring_load, f_Morison);
 ```
  to 
  ```
  [IM_nom, f_nom, Controls] = NominalSystemMatrix(q_Nom, Controls, ElastoDyn, Airfoils, Twr, Bld, Platform, WindNom, mooring_load, f_Morison);
```
3. Open the app **MATLAB Coder**.
4. In the **MATLAB Coder** app select the `NominalSystemMatrix.m` function.
5. In the next page, to autodefine the input arguments select the `OffshoreWindTurbine.m` script. Let the app run the `OffshoreWindTurbine.m` script and autodefine the inputs. 
6. In the next page, generate the code and verify it.
7. Finally, generate a `mex` file. The coder app with generate a `mex` file that's ready to be used. The file is called `NominalSystemMatrix_mex`. 
8. Return to the `SystemMatrices.m` function and change line 109 back to 
```
[IM_nom, f_nom, Controls] = NominalSystemMatrix_mex(q_Nom, Controls, ElastoDyn, Airfoils, Twr, Bld, Platform, WindNom, mooring_load, f_Morison);
```
9. Note: `MATLAB Support for MinGW-w64 C/C++ Compiler` must be installed in your MATLAB.

# Get in touch
- 👋 Hi, I’m [Saptarshi Sarkar](https://www.chalmers.se/en/Staff/Pages/Saptarshi-Sarkar.aspx), a postdoctoral researcher at the Chalmers University of Technology, Göteborg, Sweden.
- 👀 I’m interested in renewable energy systems and football :soccer:.
- :hammer: I’m currently working on the dynamics and control of wind turbines. This includes dynamics of onshore and offshore wind turbines and, structural vibration control. 
- :raised_hands: :muscle: I’m looking to collaborate on topics related to the dynamics, control and reliability of renewable energy systems, mainly wind turbines and wave energy converters.
-  :computer: Check out my  [Google Scholar](https://scholar.google.com/citations?user=3lJndhcAAAAJ&hl=en) profile or my [ORCiD](https://orcid.org/0000-0002-2111-2154) or my [ResearchGate](https://www.researchgate.net/profile/Saptarshi-Sarkar-5) profile. 
- 📫 If my work interests you and you want to collaborate, reach me at [ssarkar@chalmers.se](mailto:ssarkar@chalmers.se).

# 3. Brief description of the functions in `src`
A brief description of the functions is provided here. The description is brief and often more detail is available as comments in the functions themselves.
## `AddedMass.m`
Create the added mass matrix for the floating platform using the dimensions of the floating platform, later added to the inertia (mass) matrix at the platform DOFs. 
## `BaselineControllers.m`
This is the MATLAB representation of the DISCON controller used in FAST. This function uses the same input parameters as the DISCON controller in FAST.
## `BEMTMex.m`
Evaluate the aerodynamic loads on the blades using the Blade Element Momentum theory using [Ning's](https://doi.org/10.1002/we.1636) method [5].
## `BladeModeShapes.m`
Function to create the twisted mode shapes of the blades using the structural twist angle as the parameter. For a mathematical definition of the twisted mode shape refer to [1] and [3]. 
## `CaseOC3.xml`
Input file for the [OpenMoor](https://github.com/chen-lin/openmoor) program used to model the mooring cables. 
## `CheckInterpPoints.m`
The wind speeds are interpolated on the blades. Therefore, the grid of the wind field created using `TurbSim` must encompass the point. This function checks whether the point on the blade lies within the grid.
## `Coordinate_systems.m`
Create the inertial coordinate system and all other local coordinate systems required for the simulation.
## `coprod.m`
Multiply the position vector of a point of a turbine with a coordinate system.
## `CreateBld.m`
Create blade structural properties for simulation.
## `CreatePlatform.m`
Initiate hydrodynamic load parameters and the dimensions of the platform.
## `CreateTwr.m`
Create tower structural properties for simulation.
## `current.dat`
The x, y, z coordinates below MSL and corresponding current velocities in the x, y and z coordinates.
## `FASTTransMat.m`
Equation (2-2) in reference [4].
## `LiftDragCoeffInterp.m`
Interpolate the lift and drag coefficient of the airfoils for any give angle of attack using linear interpolation.
## `InitializeInflowAngle.m`
Initialize the inflow angles at all blade nodes. The step is required to help with the convergence at the first evaluation of the BEMT function.
## `moorapi.h`
Headfile for the OpenMoor dll. Can also be downloaded from [OpenMoor](https://github.com/chen-lin/openmoor).
## `MoorApiwin64.dll`
OpenMoor dll is loaded and called by MATLAB to estimate the mooring loads. This dll can also be downloaded from [OpenMoor](https://github.com/chen-lin/openmoor).
## `MoorDyn.dll`
MoorDyn dll that is loaded and called by MATLAB to estimate the mooring loads. This dll can also be downloaded from [MoorDyn](https://www.nrel.gov/wind/nwtc/moordyn.html). 
## `MoorDyn.h`
Headfile for the MoorDyn dll. Can also be downloaded from [MoorDyn](https://www.nrel.gov/wind/nwtc/moordyn.html). 
## `Morisons.m`
Function to evaluate the wave loads on the floating platform using Morison's equations.
## `NominalSystemMatrix_mex.mexw64`
Nominal System Matrices mex function. Mex function is used to speed up execution.
## `ode4.m`
Runge-Kutta fourth-order method for numerical integration. This is the main function that performs the time integration.
## `readBLgrid.m`
MATLAB function distributed by [NREL](https://github.com/OpenFAST/matlab-toolbox) used to read `TurbSim` generated 3D wind field in MATLAB.
## `ReadElastoDyn.m`
Function to read the `./5MW_Baseline/NRELOffshrBsline5MW_OC3Hywind_ElastoDyn.dat` file.
## `ReadWindTurbineAeroData.m`
Function to read the blade aerodynamic properties from the `./5MW_Baseline/AeroData` folder.
## `ReadWindTurbineStructuralData.m`
Read the blade's and tower's structural data from the `./5MW_Baseline/NRELOffshrBsline5MW_Blade.dat` and the `./5MW_Baseline/NRELOffshrBsline5MW_OC3Hywind_ElastoDyn_Tower.dat` respectively.
## `Results.m`
Function to evaluate the response time histories from the integrated state vector.
## `RHS.m`
Stores the right hand side of the equation $\dot{q} = f(t, q, \dot{q}, u)$
## `ServoDyn.m`
Function to initialize all control parameters. This function is used instead for the `...ServoDyn.dat` file used in FAST\OpenFAST.
## `SteadyWind.m`
If you want to simulate the turbine for steady wind, this function creates a steady wind field. 
## `SystemMatrices.m`
By adding together the Nominal System Matrix, the loads and the control inputs, this function evaluates the final system matrices i.e., the $\mathbf{M}(\mathbf{q}, t)$ and $\mathbf{f}(\mathbf{\dot{q}}, \mathbf{q}, t)$ matrices in equation 11 of reference [1].
## `Transform1.m`, `Transform2.m` and `Transform2.m`
These three functions provide three different options for applying coordinate transformations based on the dimensions of the matrices. The comments provided in the matrices will make it clear. 

_Like any set of codes, even after thorough verification, bugs may exist. If not bugs, may be possible room for improvement. In case you spot any, please don't hesitate to reach out_.
# References
[1] Sarkar, S., & Fitzgerald, B. (2021). [Use of kane’s method for multi-body dynamic modelling and control of spar-type floating offshore wind turbines](https://www.mdpi.com/1996-1073/14/20/6635). _Energies_, _14_(20), 6635.  

[2] Sarkar, S., Chen, L., Fitzgerald, B., & Basu, B. (2020). [Multi-resolution wavelet pitch controller for spar-type floating offshore wind turbines including wave-current interactions](https://doi.org/10.1016/j.jsv.2020.115170). _Journal of Sound and Vibration_, _470_, 115170.  

[3] Sarkar, S. (2020). [Individual blade pitch control strategies for spar-type floating offshore wind turbines](http://hdl.handle.net/2262/92495). _Trinity College Dublin_.  

[4] Jonkman, J. M. (2007). [_Dynamics modeling and loads analysis of an offshore floating wind turbine_](https://www.nrel.gov/docs/fy08osti/41958.pdf). University of Colorado at Boulder.

[5] Ning, S. A. (2014). [A simple solution method for the blade element momentum equations with guaranteed convergence](https://doi.org/10.1002/we.1636). _Wind Energy_, _17_(9), 1327-1345.
