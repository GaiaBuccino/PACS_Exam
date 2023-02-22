# PACS_Exam
Repository with the code associated to the PACS exam

## 0. Prerequisites: compatible version of OpenFOAM

## 1. fork and clone the repository
   ```
   git clone {url_my_repo}
   ```

## 2. SIMULATION OF THE W-Psi MODEL on the Vortex Merger test case

   ```
   cd SOLVER_and_TEST_CASE
   cd vorticity-streamFoam
   wclean
   wmake
   cd ..
   cd vortexMergerWPsi
   blockMesh
   setExprFields -dict system/setFieldsVorticity
   vorticity-streamFoam
   ```

   The results can be viewed in paraFoam
   ```
   paraFoam
   ```

## 3. VALIDATION: 
   simulation of the same case with a usual solver (velocity-pressure)
   ### 3.1. generation of the initial condition for the model U-p (using W-Psi solver on a small time step)
   ```
   cd ..
   cd vortexMergerIC_U
   blockMesh
   vorticity-streamFoam
   ```
   The results can be viewed in paraFoam
   ```
   paraFoam
   ```
   this code provides the initial condition for the u-p formulation
   ### 3.2 run Vortex Merger test case with the model U-p
   ```
   cd ..
   cd icoFoamUp
   wclean
   wmake
   cd ..
   cd vortexMergerUp
   blockMesh
   icoFoamUp 
   ```                  
   it takes a while
   The results can be viewed in paraFoam
   ```
   paraFoam
   ```

## 4. ITHACA-FV modifications in the context of the project
   This code is only to show the modifications made to ITHACA-FV library, but they have not been validated yet.

   ```
   cd ITHACA-FV
   git submodule update --init
   source etc/bashrc
   ./Allmake
   ```

   in order to see the code implemented 
   ```
   cd src/ITHACA_FOMPROBLEMS/unsteadyNSWPsi
   ```
   OR
   ```
   cd src/ITHACA_FOMPROBLEMS/steadyNSWPsi
   ```
   of course also the Make/files has been modified

## 5. Animations folder
   
In the SOLVER_and_TEST_CASE folder we have collect the animations of the simulations for each variable simulated in both formulations.



   





