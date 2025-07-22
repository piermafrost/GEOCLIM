# GEOCLIM7 version for TICE (360Ma) experiments

This current branch is modified from GEOCLIM7 (release "7.0") of "main" branch, to be used specifically for the simulations
presented in Zhang et al. "Enhanced Continental Weathering Initiates Late Paleozoic Ice Age".

This README only concerns the above-mentioned simulations.
For a complete description of GEOCLIM and of the whole GitHub repository, please refer to the README of "main" branch.


## Code modifications

* Implementation of an uplift ("source/uplift.f90"): a mask file is read at the beginning of the run (in "source/dynsoil_read_input.f90"),
  and the slope field is modified during the run according to the uplift mask and the shape of the uplift pulse.
* Implementation of lithium cycle: add two main GEOCLIM variables "Li concentration" and "d7Li", with computation of Li and 7Li budgets
  in "source/creades.f", implementation of Cave Rugenstein et al. (2019)'s model in "source/dynsoil_lithium.f90" and route Li fluxes
  in "source/cont_weath.f".


## Stored data

### netCDF inputs
This branch contains the "default" netCDF outputs used by `make_test` (ERA5, IPSL, ...), the FOAM pre-industrial outputs,
and the FOAM 360Ma outputs. Only the variables needed by GEOCLIM are present.
The other netCDF files concern the slope fields and the modern lithology field.

### Restart files
The COMBINE restart files for FOAM pre-industrial and 360 simulations are present in `restart/geoclim/`.
The DynSoil restart files are not stored, due to memory constraint.

### COMBINE boundary condition files
The COMBINE boundary condition files for FOAM pre-industrial and 360 simulations are stored in "INPUT/COMBINE/FOAM-CTRL" and
"INPUT/COMBINE/FOAM-360Ma-split" (respectively). They were generated using `preproc/BC/make_FOAM_BC.py` from raw FOAM outputs.

### Configuration files
The configuration files for the FOAM pre-industrial run, the 360Ma spin-up run, and the 360Ma "uplift" run, are stored in "config/templates/".


## HOW TO REPRODUCE THE SIMULATIONS

Here are the steps to follow to reproduce the FOAM pre-industrial and 360 simulations.

##### Make sure the code of GEOCLIM compiles and runs properly

For this purpose, the bash script `make_test`, generating test-runs, is available on the branch "main" of current GitHub repository.

##### FOAM pre-instrial simulation.

The default source and configuration files are for this run.
The only needed modification is to set "scaling_factor = 1d-3" in "source/dynsoil_physical_parameters.f90".
Compile the code with `./build_GEOCLIM --res 1,128,128 --compset default --mode optim` (and other options if needed).
Then, simply run the executable `executable/geoclim_B-10_R-1x128x128x6_DS-10lev.exe 0 1 3 0 0`

##### FOAM 360Ma spin-up equilibrium run

1. Fetch the source files for 360Ma configuration:  
    * `cp source/templates/constante.f90.FOAM-360Ma source/constante.f90`  
    * `cp source/templates/basin_geometry.f.FOAM-360Ma source/basin_geometry.f`
2. Set "scaling_factor = 1d-3" in "source/dynsoil_physical_parameters.f90"
3. Compile the code with `./build_GEOCLIM --res 7,128,128 --compset default --mode optim` (and other options if needed).
4. Fetch the configuration files for the 360Ma case:  
    * `cp config/templates/IO_CONDITIONS.FOAM-360Ma-split-lithunif-Liauth config/IO_CONDITIONS`  
    * `cp config/templates/GCM_input_conditions.FOAM-360Ma config/GCM_input_conditions`  
    * `cp config/templates/cond_p20.FOAM-360Ma-split-lithunif.dat config/cond_p20.dat`
5. Run the executable `executable/geoclim_B-10_R-7x128x128x6_DS-10lev.exe 0 1 3 0 0`

##### FOAM 360Ma with uplift perturbation

1. Follow all the instructions for "FOAM 360Ma spin-up equilibrium run" (without running the executable),
   except "scaling_factor" that must be put back to "1"  
   NOTE: if the "360Ma spin-up" configuration was already set, it is still necessary to recompile the code
   (with `build_GEOCLIM` or `make`) because of the parameter "scaling_factor" is changed.
2. Uncommented the lines "uplift_mask_fname = ..." and "uplift_mask_vname = ..." in "config/IO_CONDITIONS".
3. Use the "360Ma spin-up" COMBINE restart file ("&INIT_INFO" block in "config/IO_CONDITIONS"), set DynSoil in "restart" mode
   and use the "360Ma spin-up" DynSoil restart file  ("&DYNSOIL_INIT_INFO" in "config/IO_CONDITIONS").
4. Change the run name
5. In "source/cond_p20.dat", follow the instruction of the "NOTE ON TIME-STEPS" black:
   divide the timestep by 2, multiply the asyncrhonicity and printing factors by 2, and set the acceleration factors back to 1.
6. Run the executable `executable/geoclim_B-10_R-7x128x128x6_DS-10lev.exe 0 1 3 0 0`

