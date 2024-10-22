This solver is based on the OpenFOAM solver multiphaseCompressibleInterFoam, as modified by Hendrik Reese.

To run a simulation, one must first properly install OpenFOAMv2006:
https://www.openfoam.com/download/release-history#v2006

In the following instructions it is assumed that OpenFOAM was installed on an Ubuntu operating system.
To execute OpenFOAM commands, one must first enter the OpenFOAM environment by 
executing the alias command for the installed OpenFOAM version (usually 'of2006'),
which has to have been sourced after OpenFOAM installation.

1. Extract the folder "MultiphaseCavBubbleFoam" in the OpenFOAM installation folder under "applications/solvers/multiphase".
2. Open a terminal in "MultiphaseCavBubbleFoam" and compile it via './Allwmake'.
3. Extract the folder "exampleCases" in any location.
4. Run the simulation by opening a terminal in any simulation folder and executing './Allrun' for serial computation or './AllrunPar' for parallel computation.
5. If executed in parallel computation, reconstruct the simulation data using 'reconstructPar'.
6. To view the simulation results, open the file fluid.foam using ParaView.

The progress of the simulation may be monitored by running 'tail -f log.MultiphaseCavBubbleFoam' or 'tail -f info.csv' in the simulation folder.
Simulation settings and parameters may be altered in the text files within the simulation folder (e.g. bubble initial conditions in constant/transportProperties or simulation end time and field output interval in system/controlDict).
The simulation geometry may be altered by changing the file "system/blockMeshDict".
