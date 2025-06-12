# GEOCLIM7 version for Turonian (90Ma) experiments with Milankovitch cycles - v1.0

This current branch is modified from GEOCLIM7 (release "7.0") of "main" branch, to be used specifically for the simulations
presented in Maffre et al., "GEOCLIM7, an Earth System Model for multi-million years evolution of the geochemical cycles and climate.",
submitted to GMD (November 2024).

This README only concerns the above-mentioned simulations (hereafter called "90Ma simulations").
For a complete description of GEOCLIM and of the whole GitHub repository, please refer to the README of "main" branch.

## Data storage

### Raw netCDF inputs
The only inputs stored in the current repository are the IPSL-CM5A2 pre-industrial control IPSL-CM5A2, the slope fields (both modern
and 90Ma) and the modern lithology field.
Because they require too much memory, the climate fields of the 90Ma IPSL-CM5A2 simulations (conducted for all orbital configurations
and 2 CO2 levels) are stored in the Zenodo archive https://doi.org/10.5281/zenodo.14228131

### Restart files
The COMBINE restart files for all 90Ma simulations are present in `restart/geoclim/`.
The DynSoil restart files require too much memory, and are stored in the Zenodo archive https://doi.org/10.5281/zenodo.14228131

### COMBINE boundary condition files
The COMBINE boundary condition files are generated using `preproc/BC/BC_generator.py`, from IPSL-CMA2 simulation outputs.
The script `preproc/BC/IPSL_90Ma_all.py` was specifically written to generate all boundary conditions needed for the 90Ma simulations.
Those boundary condition files are all stored in `INPUT/COMBINE/` (one subdirectory per case). However, the IPSL outputs needed to
generate them are stored in the Zenodo archive https://doi.org/10.5281/zenodo.14228131

### Configuration files
The configuration files for all GEOCLIM simulations presented in Maffre et al. (submitted to GMD), plus additional simulations,
are stored in `config/90Ma_templates/`.
Each subdirectory corresponds to a simulation case (see following section "Summary of GEOCLIM simulations").
Only two configuration files are needed to replicate the simulation: "IO\_CONDITION" and "cond\_p20.dat".

### GEOCLIM simulation outputs
The outputs of all 90Ma simulations are not stored in this repository because of they require too much memory.
They are available on the Zenodo archive https://doi.org/10.5281/zenodo.14228131


## GEOCLIM simulations

### Nomenclature
For each individual simulation, there are 3 GEOCLIM output files (all written in `OUTPUT/`):
`geoclim_output$RUN_NAME$.nc`, `geographic_output$RUN_NAME$.nc`, and `dynsoil_output$RUN_NAME$.nc`, with "$RUN\_NAME$" being the name
of the simulation.
For the simulations presented here, all run names are in the form `.90Ma-$config$.$run_type$-$other$`, where:
* "$config$ is the main GEOCLIM configuration (details of boxes definition and specificities).
  "-AveOrb" means that climate inputs for that run are the average of all orbital configurations.
* "$run-type$" can be "2X.equil" (or simply "equil") for equilibrium spin-up simulations with pCO2 fixed at 2xpre-industrial,
  "deg.equil" for equilibrium spin-up simulations with free pCO2 (force by degassing), "Laskar" (or "Lsk") for simulations run with
  95Ma-85Ma time-series fo orbital parameters from Laskar (2004).
* "$other" indicates additional sensitivity tests.

All details can be found in the last section "Summary of all GEOCLIM simulations".

### HOW TO REPRODUCE THE SIMULATIONS

Here are the steps to follow to reproduce the 90Ma simulations presented in Maffre et al. (submitted to GMD)

0. Make sure the code of GEOCLIM compiles and runs properly

    For this purpose, the bash script `make_test`, generating test-runs, is available on the branch "main" of current GitHub repository.

1. Gather the needed inputs from the Zenodo archive https://doi.org/10.5281/zenodo.14228131

    * IPSL-CM5A2 netCDF inputs: the directory`"90Ma_Laugie/` from Zenodo archive must be put in `INPUT/IPSL/`
    * DynSoil restarts: all DynSoil restart files (`dynsoil_restart.*.nc`) from Zenodo archive must be put in `restart/dynsoil/`

    If desired, the COMBINE boundary conditions files can be remade with `preproc/BC/IPSL_90Ma_all.py`

2. Modify the source code if needed.

    The needed source code modifications are indicated, for each simulation case, in following section "Summary of all GEOCLIM simulations".

    * Spin-up equilibrium runs: the `scaling_factor` (defined in "dynsoil\_physical\_parameters.f90") must be reduced to 1d-3 for the regolith model
      to reach steady-state in due time. (Note: the acceleration of oxygen and sulfur cycles are already set in the configuration files).
    * Sensitivity experiments: one line of `geoclim_mainprog.f` should be commented/uncommented to deactivate a given process.

    **Important:**
    * When compiling the code for new simulation case, do not forget to undo the source code modification you have made.
    * Whenever the source code is modified, or put back as original, the compilation command must be re-run, and will erase the former executable.

3. Compile the code

    Use the compilation command indicated in following section "Summary of all GEOCLIM simulations" for each specific simulation case.

4. Configure GEOCLIM for the simulation

    Copy the 2 configuration files (one "IO_CONDITION..." and one "cond_p20...") of the desired simulation case, in the corresponding
    subdirectory in `config/90Ma_templates/` (see following section "Summary of all GEOCLIM simulations"), and **replace** the files
    `config/IO_CONDITIONS` and `config/cond_p20.dat`

5. Run GEOCLIM

    Run the executable file **corresponding to the desired simulation case** (i.e., the executable generated by the compilation command
    indicated in following section "Summary of all GEOCLIM simulations").

    The calculation performance is about 30 minutes per million years of simulation.
    The 10-Myr long 90Ma simulations may take up to 24 hours to complete.


### Summary of all GEOCLIM simulations

Here follows a summarized description of all the simulations corresponding to the configuration cases (directories) available in `config/90Ma_templates/`

* "globoce.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | "Classic" 10-boxes configuration, fixed CO2 (2X), equilibrium run          |
  | Config name in Maffre et al. (Table 6) | -                                                       |
  | GEOCLIM run name                       | .90Ma-globoce-AveOrb.equil                              |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --mode optim`             |
* "splitepic.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | 14-boxes config (split coastal boxes), fixed CO2 (2X), equilibrium run     |
  | Config name in Maffre et al. (Table 6) | o13                                                     |
  | GEOCLIM run name                       | .90Ma-splitepic-AveOrb.equil                            |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 14 --mode optim` |
* "AtlPac.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Split Atl/rest of mid-lat + coastal (23-boxes), fixed CO2 (2X), equil. run |
  | Config name in Maffre et al. (Table 6) | o22                                                     |
  | GEOCLIM run name                       | .90Ma-AtlPac-AveOrb.equil                               |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 23 --mode optim` |
* "AtlPac.Arct-3.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Same as previous case with 3 vertical levels in Arctic (24-boxes)          |
  | Config name in Maffre et al. (Table 6) | -                                                       |
  | GEOCLIM run name                       | .90Ma-AtlPac-Arct3-AveOrb.equil                         |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 24 --mode optim` |
* "AtlPac.Arct-3.epol.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Same as previous case with bioproduct. reduced in high-lat coastal boxes   |
  | Config name in Maffre et al. (Table 6) | -                                                       |
  | GEOCLIM run name                       | .90Ma-AtlPac-Arct3-epol-AveOrb.equil                    |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 24 --mode optim` |
* "3bas.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Split N Atl/S Atl/rest + coastal (28-boxes), fixed CO2 (2X), equil. run    |
  | Config name in Maffre et al. (Table 6) | o27                                                     |
  | GEOCLIM run name                       | .90Ma-3bas-AveOrb.equil                                 |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 28 --mode optim` |
* "3bas.Arct-3.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Same as previous case with 3 vertical levels in Arctic (29-boxes)          |
  | Config name in Maffre et al. (Table 6) | o28                                                     |
  | GEOCLIM run name                       | .90Ma-3bas-Arct3-AveOrb.equil                           |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 29 --mode optim` |
* "3bas.Arct-3.AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Same as previous case with free CO2 (imposed degassing: 5 Tmol/yr)         |
  | Config name in Maffre et al. (Table 6) | o28                                                     |
  | GEOCLIM run name                       | .90Ma-3bas-Arct3-AveOrb.deg.equil                       |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --mode optim` |
* "3bas.Arct-3.Laskar/" directory
  |                     |                                                                                                                      |
  | ------------------- | -------------------------------------------------------------------------------------------------------------------- |
  | Description         | Same configuration than "3bas.Arct-3.AveOrb.equil/" forced with "Laskar" 95Ma--85 Ma time-series of orbital cycles   |
  | Config name in Maffre et al. (Table 6) | o28                                                                                               |
  | GEOCLIM run name                       | .90Ma-3bas-Arct3.Laskar                                                                           |
  | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --clim-param 2,4,2 --param-periods ,360., --mode optim` |
* "3bas.Arct-3.epol.2X-AveOrb.equil/" directory
  * "IO\_CONDITIONS" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as "3bas.Arct-3.2X-AveOrb.equil" case with biop. reduc. in HL cst box |
    | Config name in Maffre et al. (Table 6) | o28'                                                    |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-AveOrb.2X.equil                   |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 29 --mode optim` |
  * "IO\_CONDITIONS.APx0.5" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case with N Alt/Pac exchanges reduced by 50%              |
    | Config name in Maffre et al. (Table 6) | o28'-APx0.5                                             |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.5-AveOrb.2X.equil            |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0.25" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 75%          |
    | Config name in Maffre et al. (Table 6) | o28'-APx0.25                                            |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25-AveOrb.2X.equil           |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 100%         |
    | Config name in Maffre et al. (Table 6) | o28'-APx0                                               |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0-AveOrb.2X.equil              |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
* "3bas.Arct-3.epol.AveOrb.equil/" directory
  * "IO\_CONDITIONS" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as "3bas.Arct-3.epol.2X-AveOrb.equil/" case with free CO2             |
    | Config name in Maffre et al. (Table 6) | o28'                                                    |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-AveOrb.deg.equil                  |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --mode optim` |
  * "IO\_CONDITIONS.APx0.5" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case with N Alt/Pac exchanges reduced by 50%              |
    | Config name in Maffre et al. (Table 6) | o28'-APx0.5                                             |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.5-AveOrb.deg.equil           |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0.25" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 75%          |
    | Config name in Maffre et al. (Table 6) | o28'-APx0.25                                            |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25-AveOrb.deg.equil          |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 100%         |
    | Config name in Maffre et al. (Table 6) | o28'-APx0                                               |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0-AveOrb.deg.equil             |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
* "3bas.Arct-3.epol.Laskar/" directory
  * "IO\_CONDITIONS" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same configuration than "3bas.Arct-3.epol.AveOrb.equil/" forced with "Laskar" 95Ma--85 Ma time-series of orbital cycles   |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol.Laskar                                                                           |
    | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --clim-param 2,4,2 --param-periods ,360., --mode optim`      |
  * "IO\_CONDITIONS.APx0.5" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous case with N Alt/Pac exchanges reduced by 50%                                                             |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.5.Laskar                                                                    |
    | Compilation command | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 75%                                                         |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | Exper. name in Maffre et al. (Table 8) | "All processes"                                                                                        |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.Laskar                                                                      |
    | Compilation command | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 100%                                                        |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25.Laskar                                                                   |
    | Compilation command | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceXchg" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial oceanic circulation                                                 |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | Exper. name in Maffre et al. (Table 8) | -                                                                                                      |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceXchg                                                           |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 755 `call creades_noX(time)`     |
    | Compilation command | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTemp" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial oceanic temperature                                                 |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | Exper. name in Maffre et al. (Table 8) | -                                                                                                      |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTemp                                                           |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 754 `call creades_noT(time)`     |
    | Compilation command | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteCntFlx" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial continental fluxes                                                  |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | Exper. name in Maffre et al. (Table 8) | -                                                                                                      |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-csteCntFlx                                                           |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 753 `call creades_noWth(time)`   |
    | Compilation command | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTX" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial oceanic temperature and circulation                                 |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | Exper. name in Maffre et al. (Table 8) | "Cont. fluxes"                                                                                         |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTX                                                             |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 756 `call creades_noWth(time)`   |
    | Compilation command | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceXCntFlx" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial continental fluxes and oceanic circulation                          |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | Exper. name in Maffre et al. (Table 8) | "Oce. temperature"                                                                                     |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceXchgCntFlx                                                     |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 757 `call creades_noXWth(time)`  |
    | Compilation command | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTCntFlx" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial continental fluxes and oceanic temperature                          |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | Exper. name in Maffre et al. (Table 8) | "Oce. exchanges"                                                                                       |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTempCntFlx                                                     |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 758 `call creades_noTWth(time)`  |
    | Compilation command | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTXCntFlx" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial continental fluxes, oceanic temperature and circulation             |
    | Config name in Maffre et al. (Table 6) | o28'                                                                                                   |
    | Exper. name in Maffre et al. (Table 8) | "None"                                                                                                 |
    | GEOCLIM run name                       | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTXCntFlx                                                       |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 760 `call creades_noTXWth(time)` |
    | Compilation command | *same as previous case*                                                                                                   |

