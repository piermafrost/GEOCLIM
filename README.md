# GEOCLIM7 version for Turonian (90Ma) experiments with Milankovitch cycles - v2.0

This version 2.0 is an extension of "Mil-90Ma" original branch with new experiments designed for simulating the OAE2 as a consequence
of magmatic CO2 emission. These new experiments were conducted for the study "A shift in Earth system feedback strength explains
the disappearance of OAEs in the Cenozoic" (Maffre et al., submitted to Nature Geoscience).
They were build upon the original experiments from v1.0 and published in Maffre et al. (2025) (doi: 10.5194/gmd-18-6367-2025)
This branch originally forked from GEOCLIM7 main branch (release "7.0").

## Data storage

### Raw netCDF inputs
The only inputs stored in the current repository are the IPSL-CM5A2 pre-industrial control IPSL-CM5A2, the slope fields (both modern
and 90Ma) and the modern lithology field.
Because they require too much memory, the climate fields of the 90Ma IPSL-CM5A2 simulations (conducted for all orbital configurations
and 2 CO2 levels) and pre-insudrial "no ice" simulations are stored in the Zenodo archive https://doi.org/10.5281/zenodo.14228131
(updated since first publication).

### Restart files
The COMBINE restart files for all simulations are present in `restart/geoclim/`.
The DynSoil restart files require too much memory, and are stored in the Zenodo archive https://doi.org/10.5281/zenodo.14228131
(updated since first publication).

### COMBINE boundary condition files
The COMBINE boundary condition files are generated using `preproc/BC/BC_generator.py`, from IPSL-CMA2 simulation outputs.
The script `preproc/BC/IPSL_90Ma_all.py` was specifically written to generate all boundary conditions needed for the 90Ma simulations
published in v1.0, and still used for v2.0. The script `preproc/BC/IPSL_90Ma_HR.py` was written for the new configuration (46 boxes)
added in v2.0.
Those boundary condition files are all stored in `INPUT/COMBINE/` (one subdirectory per case). However, the IPSL outputs needed to
generate them are stored in the Zenodo archive https://doi.org/10.5281/zenodo.14228131 (updated since first publication)

### Configuration files
The configuration files for all GEOCLIM simulations presented in Maffre et al. (2025) (doi: 10.5194/gmd-18-6367-2025) and Maffre et
al. (submitted to NatureGeosience), and some additional simulations, are stored in `config/90Ma_templates/`.
Each subdirectory corresponds to a simulation case (see following section "List of GEOCLIM simulations").
Only two configuration files are needed to replicate the simulation: "IO\_CONDITION" and "cond\_p20.dat".

### GEOCLIM simulation outputs
The outputs of all GEOCLIM simulations are not stored in this repository because of they require too much memory.
They are available on the Zenodo archive https://doi.org/10.5281/zenodo.14228131 (updated since first publication).


## GEOCLIM simulations

### Nomenclature
For each individual simulation, there are 3 GEOCLIM output files (all written in `OUTPUT/`):
`geoclim_output$RUN_NAME$.nc`, `geographic_output$RUN_NAME$.nc`, and `dynsoil_output$RUN_NAME$.nc`, with "$RUN\_NAME$" being the name
of the simulation.
For the simulations presented here, all run names are in the form `.$age-$config$.$run_type$-$other$`, where:
* "$age" is either "PI-NOICE" (for pre-industrial boundary conditions) or "90Ma" (for Cenomano-Turonian boundary conditions).
* "$config$ is the main GEOCLIM configuration (details of boxes definition and specificities).
  "-AveOrb" means that climate inputs for that run are the average of all orbital configurations.
* "$run-type$" can be "2X.equil" (or simply "equil") for equilibrium spin-up simulations with pCO2 fixed at 2xpre-industrial,
  "deg.equil" (sometimes simply "equil") for equilibrium spin-up simulations with free pCO2 (forced by degassing), "Laskar"
  (or "Lsk") for simulations run with95Ma-85Ma time-series fo orbital parameters from Laskar (2004).
* "$other" indicates additional sensitivity tests.

Further technical details can be found in the last section "List of all GEOCLIM simulations".

### Oceanic basins configuration
Several configuration (i.e., discretization of oceanic basins) are present in this repository. Their description is provided in the
section "List of all GEOCLIM simulations". In addition, the full discretization (# and name of all GEOCLIM boxes) is provided for each
case in the COMBINE input subrepertory `INPUT/COMBINE/*/BOXES_NUMBER_AND_NAMES.txt`.


### HOW TO REPRODUCE THE SIMULATIONS (AND BOUNDARY CONDITIONS)

Here are the steps to follow to reproduce the "v1.0" simulations presented in Maffre et al. (2025) (doi: 10.5194/gmd-18-6367-2025) and
the "v2.0" simulations presented in Maffre et al. (submitted to NatureGeosience)

0. Make sure the code of GEOCLIM compiles and runs properly

    For this purpose, the bash script `make_test`, generating test-runs, is available on the branch "main" of current GitHub repository.

1. Gather the needed inputs from the Zenodo archive https://doi.org/10.5281/zenodo.14228131

    * IPSL-CM5A2 netCDF inputs: the directory`"90Ma_Laugie/` from Zenodo archive must be put in `INPUT/IPSL/`
    * DynSoil restarts: all DynSoil restart files (`dynsoil_restart.*.nc`) from Zenodo archive must be put in `restart/dynsoil/`

    One can re-generate all boundary conditions files from the raw climate netCDF ouputs, with the following scripts:
    * `preproc/BC/IPSL_90Ma_all.py` and `preproc/BC/IPSL_90Ma_HR.py` for COMBINE boundary condition files (toolkit #1)
    * `preproc/BC/basinmap/IPSL_90Ma_allroutingmaps.py` and `preproc/BC/basinmap/IPSL_90Ma_HR_basinmap.py` for the land-to-ocean routing
       map (toolkit #2). Warning: the scripts in preproc/BC/ must be called first.
    * `preproc/generate_slope_field.py` for the non-parametric slope field reconstruction (toolkit #3)

2. Modify the source code if needed.

    The needed source code modifications are indicated, for each simulation case, in following section "List of all GEOCLIM simulations".

    * Spin-up equilibrium runs: the `scaling_factor` (defined in "dynsoil\_physical\_parameters.f90") must be reduced to 1d-3 for the regolith model
      to reach steady-state in due time. (Note: the acceleration of oxygen and sulfur cycles are already set in the configuration files).
    * Sensitivity experiments: one line of `geoclim_mainprog.f` should be commented/uncommented to deactivate a given process.
    * Configuration "46bx-TSS" (46 boxes withprimary productivity deactivated in trans-Sahara seaway): replace the file `source/Makefile` by
      `source/templates/Makefile.noSahara`. This new Makefile **will only work with 46 boxes**.

    **Important:**
    * When compiling the code for new simulation case, do not forget to undo the source code modification you have made.
    * Whenever the source code is modified, or put back as original, the compilation command must be re-run, and will erase the former executable.

3. Compile the code

    Use the compilation command indicated in following section "List of all GEOCLIM simulations" for each specific simulation case.

4. Configure GEOCLIM for the simulation

    Copy the 2 configuration files (one "IO_CONDITION..." and one "cond_p20...") of the desired simulation case, in the corresponding
    subdirectory in `config/90Ma_templates/` (see following section "List of all GEOCLIM simulations"), and **replace** the files
    `config/IO_CONDITIONS` and `config/cond_p20.dat`

5. Run GEOCLIM

    Run the executable file **corresponding to the desired simulation case** (i.e., the executable generated by the compilation command
    indicated in following section "List of all GEOCLIM simulations").

    The calculation performance ranges from ~15 minutes per million years of simulation for the smallest configurations (10 boxes) to ~6 hours
    per million yeras of simulation forthe largest configurations (46 boxes).


### List of all GEOCLIM simulations

Here follows a summarized description of all the simulations corresponding to the configuration cases (directories) available in `config/90Ma_templates/`

#### Spin-up "equilibrium" simulations

##### Pre-industrial boundary conditions

* "PI-NOICE.equil/" directory
  * "IO\_CONDITIONS" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | "Classic" 10-boxes configuration, free CO2 (imposed degassing: 5Tmol/yr)   |
    | Used in version     | v2.0                                                                       |
    | Config name in Maffre et al. (sub.) [Tab. 1] | modern 10bx                                       |
    | GEOCLIM run name                             | .PI-NOICE.deg.equil                               |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --mode optim`             |
    | COMBINE input dir   | "INPUT/COMBINE/IPSL-PI-NOICE/"                                             |
  * "IO\_CONDITIONS.lithunif" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case with uniform lithology                               |
    | Used in version     | v2.0                                                                       |
    | GEOCLIM run name                             | .PI-NOICE-lithunif.eq                             |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
    | COMBINE input dir   | *same as previous case*                                                    |

##### 90Ma boundary conditions

* "globoce.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | "Classic" 10-boxes configuration, fixed CO2 (2X), equilibrium run          |
  | Used in version     | v1.0                                                                       |
  | GEOCLIM run name                             | .90Ma-globoce-AveOrb.equil                        |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --mode optim`             |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-2X/"                                              |
* "globoce.AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Same as previous case with free CO2 (imposed degassing: 5 Tmol/yr)         |
  | Used in version     | v2.0                                                                       |
  | Config name in Maffre et al. (sub.) [Tab. 1] | CT 10bx                                           |
  | GEOCLIM run name                             | .90Ma-globoce-AveOrb.equil                        |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | *same as previous case*                                                    |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma/"                                                 |
* "splitepic.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | 14-boxes config (split coastal boxes), fixed CO2 (2X), equilibrium run     |
  | Used in version     | v1.0                                                                       |
  | Config name in Maffre et al. (2025) [Tab. 6] | o13                                               |
  | GEOCLIM run name                             | .90Ma-splitepic-AveOrb.equil                      |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 14 --mode optim` |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-splitepic-2X/"                                    |
* "AtlPac.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Split Atl/rest of mid-lat + coastal (23-boxes), fixed CO2 (2X), equil. run |
  | Used in version     | v1.0                                                                       |
  | Config name in Maffre et al. (2025) [Tab. 6] | o22                                               |
  | GEOCLIM run name                             | .90Ma-AtlPac-AveOrb.equil                         |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 23 --mode optim` |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-AtlPac-2X/"                                       |
* "AtlPac.Arct-3.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Same as previous case with 3 vertical levels in Arctic (24-boxes)          |
  | Used in version     |                                                                            |
  | Config name in Maffre et al. (2025) [Tab. 6] | -                                                 |
  | GEOCLIM run name                             | .90Ma-AtlPac-Arct3-AveOrb.equil                   |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 24 --mode optim` |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-AtlPac-Arct3-2X/"                                 |
* "AtlPac.Arct-3.epol.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Same as previous case with bioproduct. reduced in high-lat coastal boxes   |
  | Used in version     | v1.0                                                                       |
  | Config name in Maffre et al. (2025) [Tab. 6] | -                                                 |
  | GEOCLIM run name                             | .90Ma-AtlPac-Arct3-epol-AveOrb.equil              |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 24 --mode optim` |
  | COMBINE input dir   | *same as previous case*                                                    |
* "3bas.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Split N Atl/S Atl/rest + coastal (28-boxes), fixed CO2 (2X), equil. run    |
  | Used in version     | v1.0                                                                       |
  | Config name in Maffre et al. (2025) [Tab. 6] | o27                                               |
  | GEOCLIM run name                             | .90Ma-3bas-AveOrb.equil                           |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 28 --mode optim` |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-3bas-2X/"                                         |
* "3bas.Arct-3.2X-AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Same as previous case with 3 vertical levels in Arctic (29-boxes)          |
  | Used in version     | v1.0                                                                       |
  | Config name in Maffre et al. (2025) [Tab. 6] | o28                                               |
  | GEOCLIM run name                             | .90Ma-3bas-Arct3-AveOrb.equil                     |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 29 --mode optim` |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-3bas-Arct3-2X/"                                   |
* "3bas.Arct-3.AveOrb.equil/" directory
  |                     |                                                                            |
  | ------------------- | -------------------------------------------------------------------------- |
  | Description         | Same as previous case with free CO2 (imposed degassing: 5 Tmol/yr)         |
  | Used in version     | v1.0                                                                       |
  | Config name in Maffre et al. (2025) [Tab. 6] | o28                                               |
  | GEOCLIM run name                             | .90Ma-3bas-Arct3-AveOrb.deg.equil                 |
  | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --mode optim` |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-3bas-Arct3/"                                      |
* "3bas.Arct-3.epol.2X-AveOrb.equil/" directory
  * "IO\_CONDITIONS" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as "3bas.Arct-3.2X-AveOrb.equil" case with biop. reduc. in HL cst box |
    | Used in version     | v1.0                                                                       |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                              |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-AveOrb.2X.equil             |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 29 --mode optim` |
    | COMBINE input dir   | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0.5" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case with N Alt/Pac exchanges reduced by 50%              |
    | Used in version     | v1.0                                                                       |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'-APx0.5                                       |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.5-AveOrb.2X.equil      |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
    | COMBINE input dir   | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0.25" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 75%          |
    | Used in version     | v1.0                                                                       |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'-APx0.25                                      |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25-AveOrb.2X.equil     |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
    | COMBINE input dir   | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 100%         |
    | Used in version     | v1.0                                                                       |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'-APx0                                         |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0-AveOrb.2X.equil        |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
    | COMBINE input dir   | *same as previous case*                                                    |
* "3bas.Arct-3.epol.AveOrb.equil/" directory
  * "IO\_CONDITIONS" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as "3bas.Arct-3.epol.2X-AveOrb.equil/" case with free CO2             |
    | Used in version     | v1.0                                                                       |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                              |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-AveOrb.deg.equil            |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --mode optim` |
    | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-3bas-Arct3/"                                      |
  * "IO\_CONDITIONS.APx0.5" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case with N Alt/Pac exchanges reduced by 50%              |
    | Used in version     | v1.0                                                                       |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'-APx0.5                                       |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.5-AveOrb.deg.equil     |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
    | COMBINE input dir   | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0.25" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 75%          |
    | Used in version     | v1.0 and v2.0                                                              |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'-APx0.25                                      |
    | Config name in Maffre et al. (sub.) [Tab. 1] | CT 29bx                                           |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25-AveOrb.deg.equil    |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
    | COMBINE input dir   | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 100%         |
    | Used in version     | v1.0                                                                       |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'-APx0                                         |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0-AveOrb.deg.equil       |
    | Code modification   | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
    | COMBINE input dir   | *same as previous case*                                                    |
* "HR.AveOrb.equil/" directory ("IO\_CONDITIONS.APTx0.25" file)
  * "cond\_p2O.dat" file
    * default case
      |                     |                                                                            |
      | ------------------- | -------------------------------------------------------------------------- |
      | Description         | Mid-latitudes split in N. Atl., S. Atl., Tethys and Pacific;
                              S. high latitudes split in Atl.-facing and Pac.-facing;
                              Arctic, N. Atl. and Tethys coastal boxes split by continent;
                              N. Atl/Pacific and N. Atl./Tethys exchanged reduced by 75%;
                              free CO2 (imposed degassing: 5 Tmol/yr)                                    |
      | Used in version     | v2.0                                                                       |
      | Config name in Maffre et al. (sub.) [Tab. 1] | 46 bx                                             |
      | GEOCLIM run name                             | .90Ma-HR-APTx0.25-AveOrb.deg.equil                |
      | Code modification   | Replace `source/Makefile` by `source/templates/Makefile.noSahara`;
                              Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
      | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 46 --mode optim` |
      | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-HR/"                                              |
    * "no trans-Sahara seaway (TSS)" case
      |                     |                                                                            |
      | ------------------- | -------------------------------------------------------------------------- |
      | Description         | Same as previous case with primary productivity deaactivated in TSS        |
      | Used in version     | v2.0                                                                       |
      | Config name in Maffre et al. (sub.) [Tab. 1] | 46 bx-TSS                                         |
      | GEOCLIM run name                             | .90Ma-HR-APTx0.25-AveOrb.noSah.equil              |
      | Code modification   | Replace `source/Makefile` by `source/templates/Makefile.noSahara`;
                              Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
      | Compilation command | *same as previous case*                                                    |
      | COMBINE input dir   | *same as previous case*                                                    |
  * "cond\_p2O-warmer.dat" file
    |                     |                                                                            |
    | ------------------- | -------------------------------------------------------------------------- |
    | Description         | Same as previous case (with TSS deactivation) with  degass. of 6.6 Tmol/yr |
    | Used in version     | v2.0                                                                       |
    | Config name in Maffre et al. (sub.) [Tab. 1] | 46 bx-TSS warm                                    |
    | GEOCLIM run name                             | .90Ma-HR-APTx0.25-AveOrb.noSah.warmer.equil       |
    | Code modification   | Replace `source/Makefile` by `source/templates/Makefile.noSahara`;
                            Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command | *same as previous case*                                                    |
    | COMBINE input dir   | *same as previous case*                                                    |

#### Simulations with "perturbation"

##### Pre-industrial boundary conditions - LIP degassing

* "PI-LIPOAE/" directory
  * "IO\_CONDITIONS" file
    * default case
      |                     |                                                                                                                                                    |
      | ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
      | Description         | "Classic" 10-boxes configuration, free CO2, LIP degassing scenario                                                                                 |
      | Used in version     | v2.0                                                                                                                                               |
      | Config name in Maffre et al. (sub.) [Tab. 1] | modern 10bx                                                                                                               |
      | GEOCLIM run name (Clarkson degas. scenario)  | .PI-NOICE.LIPOAE                                                                                                          |
      | GEOCLIM run name (Papadom. degas. scenario)  | .PI-NOICE.LIPNina                                                                                                         |
      | Code modification (Clarkson degas. scenario) | Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                           |
      | Code modification (Papadom. degas. scenario) | Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPNina.inc`                                          |
      | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --mode optim`                                                                                     |
      | COMBINE input dir   | "INPUT/COMBINE/IPSL-PI-NOICE/"                                                                                                                     |
    * with phosphorus ositive feedback
      |                     |                                                                                                                                                    |
      | ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
      | Description         | Same as previous case with phosphorus positive feedback (no inorganic phosphorus sink)                                                             |
      | Used in version     | v2.0                                                                                                                                               |
      | Config name in Maffre et al. (sub.) [Tab. 1] | N/A                                                                                                                       |
      | GEOCLIM run name (Clarkson degas. scenario)  | .PI-NOICE.LIPOAE.noPHydPhos                                                                                               |
      | Code modification (Clarkson degas. scenario) | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 755 `call creades_noX(time)`     |
                                                       Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                           |
      | Compilation command | *same as previous case*                                                                                                                            |
      | COMBINE input dir   | *same as previous case*                                                                                                                            |
  * "IO\_CONDITIONS.lithunif" file
    |                     |                                                                                                                                                    |
    | ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as two cases above with uniform lithology                                                                                                     |
    | Used in version     | v2.0                                                                                                                                               |
    | Config name in Maffre et al. (sub.) [Tab. 1] | N/A                                                                                                                       |
    | GEOCLIM run name (Clarkson degas. scenario)  | .PI-NOICE-lithunif.LIPOAE                                                                                                 |
    | Code modification (Clarkson degas. scenario) | Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                           |
    | Compilation command | *same as previous case*                                                                                                                            |
    | COMBINE input dir   | *same as previous case*                                                                                                                            |

##### 90Ma boundary conditions - orbital forcing

* "3bas.Arct-3.Laskar/" directory
  |                     |                                                                                                                      |
  | ------------------- | -------------------------------------------------------------------------------------------------------------------- |
  | Description         | Same configuration than "3bas.Arct-3.AveOrb.equil/" forced with "Laskar" 95Ma--85 Ma time-series of orbital cycles   |
  | Used in version     | v1.0                                                                                                                 |
  | Config name in Maffre et al. (2025) [Tab. 6] | o28                                                                                         |
  | GEOCLIM run name                             | .90Ma-3bas-Arct3.Laskar                                                                     |
  | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --clim-param 2,4,2 --param-periods ,360., --mode optim` |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-3bas-Arct3/"                                                                                |
* "3bas.Arct-3.epol.Laskar/" directory
  * "IO\_CONDITIONS" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same configuration than "3bas.Arct-3.epol.AveOrb.equil/" forced with "Laskar" 95Ma--85 Ma time-series of orbital cycles   |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol.Laskar                                                                     |
    | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --clim-param 2,4,2 --param-periods ,360., --mode optim`      |
    | COMBINE input dir   | *same as previous case*                                                                                              |
  * "IO\_CONDITIONS.APx0.5" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous case with N Alt/Pac exchanges reduced by 50%                                                             |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.5.Laskar                                                              |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 75%                                                         |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | Exper. name in Maffre et al. (2025) [Tab. 8] | "All processes"                                                                                  |
    | GEOCLIM run name                           | .90Ma-3bas-Arct3-epol-APx0.Laskar                                                                  |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous case but N Alt/Pac exchanges area reduced by 100%                                                        |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25.Laskar                                                             |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceXchg" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial oceanic circulation                                                 |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | Exper. name in Maffre et al. (2025) [Tab. 8] | -                                                                                                |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceXchg                                                     |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 755 `call creades_noX(time)`     |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTemp" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial oceanic temperature                                                 |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | Exper. name in Maffre et al. (2025) [Tab. 8] | -                                                                                                |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTemp                                                     |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 754 `call creades_noT(time)`     |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteCntFlx" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial continental fluxes                                                  |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | Exper. name in Maffre et al. (2025) [Tab. 8] | -                                                                                                |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-csteCntFlx                                                     |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 753 `call creades_noWth(time)`   |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTX" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial oceanic temperature and circulation                                 |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | Exper. name in Maffre et al. (2025) [Tab. 8] | "Cont. fluxes"                                                                                   |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTX                                                       |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 756 `call creades_noWth(time)`   |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceXCntFlx" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial continental fluxes and oceanic circulation                          |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | Exper. name in Maffre et al. (2025) [Tab. 8] | "Oce. temperature"                                                                               |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceXchgCntFlx                                               |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 757 `call creades_noXWth(time)`  |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTCntFlx" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial continental fluxes and oceanic temperature                          |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | Exper. name in Maffre et al. (2025) [Tab. 8] | "Oce. exchanges"                                                                                 |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTempCntFlx                                               |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 758 `call creades_noTWth(time)`  |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTXCntFlx" file
    |                     |                                                                                                                           |
    | ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous "APx0.25" case with constant initial continental fluxes, oceanic temperature and circulation             |
    | Used in version     | v1.0                                                                                                                      |
    | Config name in Maffre et al. (2025) [Tab. 6] | o28'                                                                                             |
    | Exper. name in Maffre et al. (2025) [Tab. 8] | "None"                                                                                           |
    | GEOCLIM run name                             | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTXCntFlx                                                 |
    | Code modification   | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 760 `call creades_noTXWth(time)` |
    | Compilation command | *same as previous case*                                                                                                   |
    | COMBINE input dir   | *same as previous case*                                                                                                   |

##### 90Ma boundary conditions - LIP degassing

* "LIPOAE2\_globoce.AveOrb/" directory
  |                     |                                                                                                                                                       |
  | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
  | Description         | Same configuration than "globoce.AveOrb.equil/" forced with LIP degassing scenario                                                                    |
  | Used in version     | v2.0                                                                                                                                                  |
  | Config name in Maffre et al. (sub.) [Tab. 1] | modern 10bx                                                                                                                  |
  | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-globoce-AveOrb.LIPOAE                                                                                                  |
  | Code modification (Clarkson degas. scenario) | Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                              |
  | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --mode optim`                                                                                        |
  | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma/"                                                                                                                            |

* "LIPOAE2\_3bas.Arct-3.epol.AveOrb/" directory
  * default case
    |                     |                                                                                                                                                       |
    | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same configuration than "3bas.Arct-3.epol.AveOrb.equil/" (APx0.25) forced with LIP degassing scenario                                                 |
    | Used in version     | v2.0                                                                                                                                                  |
    | Config name in Maffre et al. (sub.) [Tab. 1] | CT 29bx                                                                                                                      |
    | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-3bas-Arct3-epol-APx0.25-AveOrb.LIPOAE                                                                                  |
    | Code modification (Clarkson degas. scenario) | Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                              |
    | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --mode optim`                                                                            |
    | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-3bas-Arct3/"                                                                                                                 |
  * Case with LIP CO2 d13C = -5 permil
    |                     |                                                                                                                                                       |
    | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as previous case with LIP CO2 d13C signature of-5 permil instead of 0 permil                                                                     |
    | Used in version     | v2.0                                                                                                                                                  |
    | Config name in Maffre et al. (sub.) [Tab. 1] | CT 29bx                                                                                                                      |
    | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-3bas-Arct3-epol-APx0.25-AveOrb.LIPOAE_-5pm                                                                             |
    | Code modification (Clarkson degas. scenario) | Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE_-5pm.inc`                                         |
    | Compilation command | *same as previous case*                                                                                                                               |
    | COMBINE input dir   | *same as previous case*                                                                                                                               |
  * Case with phosphorus positive feedback
    |                     |                                                                                                                                                       |
    | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same as two cases above with phosphorus positive feedback (no inorganic P sink)                                                                       |
    | Used in version     | v2.0                                                                                                                                                  |
    | Config name in Maffre et al. (sub.) [Tab. 1] | N/A                                                                                                                          |
    | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-3bas-Arct3-epol-APx0.25-AveOrb.noPHydPhos.LIPOAE                                                                       |
    | Code modification (Clarkson degas. scenario) | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 752 `call creades_noPHydPhos(time)`;
                                                     Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE_-5pm.inc`                                         |
    | Compilation command | *same as previous case*                                                                                                                               |
    | COMBINE input dir   | *same as previous case*                                                                                                                               |

* "LIPOAE2\_HR.AveOrb/" directory
  * "cond\_p20.dat" file
    * default case
      |                     |                                                                                                                                                       |
      | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
      | Description         | Same configuration than "HR.AveOrb.equil/" forced with LIP degassing scenario                                                                         |
      | Used in version     | v2.0                                                                                                                                                  |
      | Config name in Maffre et al. (sub.) [Tab. 1] | CT 46bx                                                                                                                      |
      | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-HR-APTx0.25-AveOrb.LIPOAE                                                                                              |
      | Code modification (Clarkson degas. scenario) | Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                              |
      | Compilation command | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 46 --mode optim`                                                                            |
      | COMBINE input dir   | "INPUT/COMBINE/IPSL-90Ma-HR/"                                                                                                                         |
    * case with phosphorus positive feedback                                  
      |                     |                                                                                                                                                       |
      | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
      | Description         | Same as previous case with phosphorus positive feedback (no inorganic phosphorus sink)                                                                |
      | Used in version     | v2.0                                                                                                                                                  |
      | Config name in Maffre et al. (sub.) [Tab. 1] | N/A                                                                                                                          |
      | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-HR-APTx0.25-AveOrb.noPHydPhos.LIPOAE                                                                                   |
      | Code modification (Clarkson degas. scenario) | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 752 `call creades_noPHydPhos(time)`;
                                                       Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                              |
      | Compilation command | *same as previous case*                                                                                                                               |
      | COMBINE input dir   | *same as previous case*                                                                                                                               |
    * case with primary productivity deactived in trans-Sahara seaway (TSS)
      |                     |                                                                                                                                                       |
      | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
      | Description         | Same as "default" case with primary productivity deactived in TSS                                                                                     |
      | Used in version     | v2.0                                                                                                                                                  |
      | Config name in Maffre et al. (sub.) [Tab. 1] | CT 46bx-TSS                                                                                                                  |
      | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-HR-APTx0.25-AveOrb.noSah.LIPOAE                                                                                        |
      | Code modification (Clarkson degas. scenario) | Replace `source/Makefile` by `source/templates/Makefile.noSahara`;
                                                       Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                              |
      | Compilation command | *same as previous case*                                                                                                                               |
      | COMBINE input dir   | *same as previous case*                                                                                                                               |
    * case with primary productivity deactived in TSS and phosphorus positive feedback
      |                     |                                                                                                                                                       |
      | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
      | Description         | Same as previous case with phosphorus positive feedback (no inorganic phosphorus sink)                                                                |
      | Used in version     | v2.0                                                                                                                                                  |
      | Config name in Maffre et al. (sub.) [Tab. 1] | CT 46bx-TSS P-fbk                                                                                                            |
      | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-HR-APTx0.25-AveOrb.noSah.noPHydPhos.LIPOAE                                                                             |
      | Code modification (Clarkson degas. scenario) | Replace `source/Makefile` by `source/templates/Makefile.noSahara`;
                                                       In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 752 `call creades_noPHydPhos(time)`;
                                                       Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                              |
      | GEOCLIM run name (Papadom. degas. scenario)  | .90Ma-HR-APTx0.25-AveOrb.noSah.noPHydPhos.LIPNina                                                                            |
      | Code modification (Papadom. degas. scenario) | Same as "Clarkson" scenario, but replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPNina.inc`            |
      | Compilation command | *same as previous case*                                                                                                                               |
      | COMBINE input dir   | *same as previous case*                                                                                                                               |
    * case with primary productivity deactived in TSS and LIP CO2 d13C = -5 permil       
      |                     |                                                                                                                                                       |
      | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
      | Description         | Same as two cases above case with d13C signature of LIP emitted CO2 set to -5 permil (instead of 0 pm)                                                |
      | Used in version     | v2.0                                                                                                                                                  |
      | Config name in Maffre et al. (sub.) [Tab. 1] | CT 46bx-TSS                                                                                                                  |
      | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-HR-APTx0.25-AveOrb.noSah.LIPOAE_-5pm                                                                                   |
      | Code modification (Clarkson degas. scenario) | Replace `source/Makefile` by `source/templates/Makefile.noSahara`;
                                                       Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE_-5pm.inc`                                         |
      | GEOCLIM run name (Papadom. degas. scenario)  | .90Ma-HR-APTx0.25-AveOrb.noSah.LIPNina_-5pm                                                                                  |
      | Code modification (Papadom. degas. scenario) | Same as "Clarkson" scenario, but replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPNina_-5pm.inc`       |
      | Compilation command | *same as previous case*                                                                                                                               |
      | COMBINE input dir   | *same as previous case*                                                                                                                               |
    * case with primary productivity deactived in TSS, phosphorus positive feedback, and LIP CO2 d13C = -5 permil       
      |                     |                                                                                                                                                       |
      | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
      | Description         | Same as two cases above case with d13C signature of LIP emitted CO2 set to -5 permil (instead of 0 pm)                                                |
      | Used in version     | v2.0                                                                                                                                                  |
      | Config name in Maffre et al. (sub.) [Tab. 1] | CT 46bx-TSS P-fbk                                                                                                            |
      | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-HR-APTx0.25-AveOrb.noSah.noPHydPhos.LIPOAE_-5pm                                                                        |
      | Code modification (Clarkson degas. scenario) | Replace `source/Makefile` by `source/templates/Makefile.noSahara`;
                                                       In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 752 `call creades_noPHydPhos(time)`;
                                                       Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE_-5pm.inc`                                         |
      | Compilation command | *same as previous case*                                                                                                                               |
      | COMBINE input dir   | *same as previous case*                                                                                                                               |
  * "cond\_p20-warmer.dat" file
    |                     |                                                                                                                                                       |
    | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
    | Description         | Same configuration than "HR.AveOrb.equil/" with degass. set at 6.6 Tmol/yr, primary product. deactivated in TSS, and forced with LIP degass. scenario |
    | Used in version     | v2.0                                                                                                                                                  |
    | Config name in Maffre et al. (sub.) [Tab. 1] | CT 46bx-TSS warmer                                                                                                           |
    | GEOCLIM run name (Clarkson degas. scenario)  | .90Ma-HR-APTx0.25-AveOrb.noSah.warmer.LIPOAE                                                                                 |
    | Code modification (Clarkson degas. scenario) | Replace `source/Makefile` by `source/templates/Makefile.noSahara`;
                                                     Replace `source/combine_foam.inc` by `source/templates/combine_foam_LIPOAE.inc`                                              |
    | Compilation command | *same as previous case*                                                                                                                               |
    | COMBINE input dir   | *same as previous case*                                                                                                                               |

##### 90Ma boundary conditions - orbital forcing and LIP degassing

Those simulations were not used in any of the scientific articles.
Their configuration files are "LIPOAE2\_3bas.Arct-3.epol.Laskar/" and "LIPOAE2\_HR.Laskar/" directories

