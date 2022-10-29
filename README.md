# 1. Indroduction
This package provides a MATLAB based simulation tool to analyse the dynamics of a spar-type floating offshore wind turbine using Kane's method. The model is set-up with the 5 MW OC3 Hywind Spar type FOWT. However, the codes are generic and will work for wind turbines with a different rating. For detail on the theory and derivation refer [1].
The inputs are formatted as in [OpenFAST](https://github.com/OpenFAST/openfast). All structural and aerodynamic input data is read from OpenFAST input files. However, inputs related to the wind field, wave field and control systems are initialized separety as described in section 3.  

The mooring cables can be modelled using the [MoorDyn](https://www.nrel.gov/wind/nwtc/moordyn.html) or the [OpenMoor](https://github.com/chen-lin/openmoor) program. When using [OpenMoor](https://github.com/chen-lin/openmoor) the users have the additional capability of modelling wave-current interaction. For details on the theory of wave-current interaction refer [2].
# 2. How to use
Three steps to run the code. 
## Step 1 - Create wind field using `TurbSim`
Use NREL's [TurbSim](https://www.nrel.gov/wind/nwtc/turbsim.html) tool to generate the desired 3D wind field. Please ensure that the rotor swept area is entirely enclosed. Once generated, place the `*.sum` and `*.wnd` files in the `Wind_Dataset` folder. 

## Step 2 - Create wave field using `CreateWaveDataset.m` script
In the folder `Wave_Dataset` find the MATLAB script `CreateWaveDataset.m`.  Use the script to create **regular wave** (sinusoidal wave) field or use the **PiersonMoskowitz** or **JONSWAP** spectrum. You can also include wave-current interaction is the generate wave field. For details on the wave current interaction model refer [2].

## Step 3 - Run the main simulation using `OffshoreWindTurbine.m` script
The MATLAB script `OffshoreWindTurbine.m` is used to run the simulations. The script is described breifly below.
The code starts with some initial house keeping.
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
Create the modal propertise required for the dynamic analysis. Also, select the mooring cable model and select if wave loads are to be calculated or not.
```
% Create tower and blade structures
Twr = CreateTwr(Tower,ElastoDyn); Bld = CreateBld(ElastoDyn,Geometry,Blade); Platform = CreatePlatform();
Platform.Mooring     = 1;    % 1 for Moordyn, 2 for OpenMoor
Platform.WaveLoads   = 1;    % 1 to calculate wave load on spar using Morisson's equation, 0 to not.
```
If wave loads are include, place the correct file name (name of the ``.mat`` file created using the `CreateWaveDataset.m` script in the below line
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
% Get the TurbSim generated wind field and grid
[velocity, Wind.y, Wind.z, Wind.nz, Wind.ny, Wind.dz,...
          Wind.dy, Wind.dt, Wind.zHub, Wind.z1, Wind.SummVars] = readBLgrid('1ETMC');
```
For steady wind field, uncomment the next line
```
% [Wind.t_TurbSim, velocity] = SteadyWind(Wind.y, Wind.z, 11.4, 0.2);
```
Turn on/off **Pitt and peters** correction for skewed wake, and **Blade pitch control**. The default settings are given below. In the default setting the **blade pitch control** must all be turned on.
```
WindNom.PittandPeters = true(0);
WindNom.PitchControl  = true(1);
```
Next, select start time ``t0``, final time ``tf`` and time step ``deltat``
```
t0  = 0;
tf  = 1;
deltat = 0.0125;
```
This concludes the set-up of the `OffshoreWindTurbine.m` script. The scripts runs using ``ODE 4`` RK $4^{th}$ order method and plots the results.
## Step 4* - Build  `./src/NominalSystemMatrix_mex.mex64` 
*This is an optional step. 
This package uses MATLAB mex files to accelerate the simulation time. The `NominalSystemMatrix.m` function is identified as the bottleneck function. Therefore, it has been build into a mex (MATLAB executable) function. The mex function is already provided in the `src` folder.  It should should work normally without the need for re-building it. However, in case your MATLAB version does not support the pre-build mex file follow the following steps
1.  Move into the `src` folder.
2. Type `mex NominalSystemMatrix.m` on the MATLAB command window.
3. Note: `MATLAB Support for MinGW-w64 C/C++ Compiler` must be installed in your MATLAB.

# Get in touch
- 👋 Hi, I’m [Saptarshi Sarkar](https://www.chalmers.se/en/Staff/Pages/Saptarshi-Sarkar.aspx), postdoctoral researcher at Chalmers University of Technology, Göteborg, Sweden.
- 👀 I’m interested in renewable energy systems and football :soccer:.
- :hammer: I’m currently working on dynamics and control of wind turbines. This includes dynamics onshore and offshore wind turbines and, structural vibration control. 
- :raised_hands: :muscle: I’m looking to collaborate on topics related to dynamics, control and reliability of renewable energy systems, mainly wind turbines and wave energy converters.
-  :computer: Check out my  [Google Scholar](https://scholar.google.com/citations?user=3lJndhcAAAAJ&hl=en) profile or my [ORCiD](https://orcid.org/0000-0002-2111-2154) or my [ResearchGate](https://www.researchgate.net/profile/Saptarshi-Sarkar-5) profile. 
- 📫 If my work interests you and you want to collaborate, reach me at [ssarkar@chalmers.se](mailto:ssarkar@chalmers.se).

# 3. Breif description of the functions in `src` folder
A breif description of the functions are provided here. The description is breif and often more detail is available as comments in the functions iteself.
## `AddedMass.m`

## `BaselineControllers.m`

## `BladeModeShapes.m`

## `CaseOC3.xml`

## `Coordinate_systems.m`

## `CreateBld.m`

## `CreatePlatform.m`

## `CreateTwr.m`

## `current.dat`

## `FASTTransMat.m`

## `InitializeInflowAngle.m`

## `moorapi.h`

## `MoorApiwin64.dll`

## `MoorDyn.dll`

## `MoorDyn.h`

## `Morisons.m`

## `NominalSystemMatrix_mex.mexw64`

## `ode4.m`

## `ReadElastoDyn.m`

## `ReadWindTurbineAeroData.m`

## `ReadWindTurbineStructuralData.m`

## `Results.m`

## `RHS.m`

## `ServoDyn.m`

## `SteadyWind.m`

## `SystemMatrices.m`

## `Transform1.m`, `Transform2.m` and `Transform2.m`

_Please don't hesitate to reach out with bugs or typos in the code_.
# References
[1] Sarkar, S., & Fitzgerald, B. (2021). Use of kane’s method for multi-body dynamic modelling and control of spar-type floating offshore wind turbines. _Energies_, _14_(20), 6635.
[2] Sarkar, S., Chen, L., Fitzgerald, B., & Basu, B. (2020). Multi-resolution wavelet pitch controller for spar-type floating offshore wind turbines including wave-current interactions. _Journal of Sound and Vibration_, _470_, 115170.
