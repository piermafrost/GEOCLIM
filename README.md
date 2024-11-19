# GEOCLIM7.0 version for Turonian (90Ma) experiments with Milankovitch cycles

This current branch is modified from GEOCLIM7.0 (release ...) of "main" branch, to be used specifically for the
simulations presented in Maffre et al., submitted to GMD (2024).

This README only concerns the above-mentioned simulations (hereafter called "90Ma simulations").
For a complete description of GEOCLIM and of the whole Github repository, please refer to the README of "main" branch.

## Data storage

### Raw netCDF inputs
A few inputs are stored here...

### Restart files
COMBINE restart files for all 90Ma simulations are present in `restart/geoclim`.
DynSoil restart files require too much memory, and are stored in the PANGAEA archive ... 

### COMBINE boundary condition files
The COMBINE boundary condition files are generated using `preproc/BC/BC_generator.py`, from IPSL-CMA2 simulation outputs.
The script `preproc/BC/IPSL_90Ma_all.py` was specifically written to generate all boundary conditions needed for the 90Ma simulations.
Those boundary condition files are all stored in `INPUT/COMBINE/` (one subdirectory per case). However, the IPSL outputs needed to
generate them are stored in the PANGAEA archive ...

### Configuration files
The configuration files for all GEOCLIM simulations presented in Maffre et al. (submitted to GMD), plus additional simulations,
are stored in `config/90Ma_templates/`.
Each subdirectory corresponds to a simulation case (see following section "Summary of GEOCLIM simulations").
Only two configuration files are needed to replicate the simulation: "IO\_CONDITION" and "cond\_p20.dat".

### GEOCLIM simulation sutputs
The outputs of all 90Ma simulations are not stored in this repository because of they require too much memory.
They are available on the PANGAEA archive ...


## GEOCLIM simulations

### Run nomenclature
...

### Summary of GEOCLIM simulations

* "globoce.2X-AveOrb.equil/" directory
  |                       |                                                                            |
  | --------------------- | -------------------------------------------------------------------------- |
  | Article's config name | -                                                                          |
  | GEOCLIM run name      | .90Ma-globoce-AveOrb.equil                                                 |
  | Description           | ...                                                                        |
  | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command   | `./build_GEOCLIM --compset default --res 1,96,96 --mode optim`             |
* "splitepic.2X-AveOrb.equil/" directory
  |                       |                                                                            |
  | --------------------- | -------------------------------------------------------------------------- |
  | Article's config name | o13                                                                        |
  | GEOCLIM run name      | .90Ma-splitepic-AveOrb.equil                                               |
  | Description           | ...                                                                        |
  | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command   | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 14 --mode optim` |
* "AtlPac.2X-AveOrb.equil/" directory
  |                       |                                                                            |
  | --------------------- | -------------------------------------------------------------------------- |
  | Article's config name | o22                                                                        |
  | GEOCLIM run name      | .90Ma-AtlPac-AveOrb.equil                                                  |
  | Description           | ...                                                                        |
  | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command   | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 23 --mode optim` |
* "AtlPac.Arct-3.2X-AveOrb.equil/" directory
  |                       |                                                                            |
  | --------------------- | -------------------------------------------------------------------------- |
  | Article's config name | -                                                                          |
  | GEOCLIM run name      | .90Ma-AtlPac-Arct3-AveOrb.equil                                            |
  | Description           | ...                                                                        |
  | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command   | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 24 --mode optim` |
* "AtlPac.Arct-3.epol.2X-AveOrb.equil/" directory
  |                       |                                                                            |
  | --------------------- | -------------------------------------------------------------------------- |
  | Article's config name | -                                                                          |
  | GEOCLIM run name      | .90Ma-AtlPac-Arct3-epol-AveOrb.equil                                       |
  | Description           | ...                                                                        |
  | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command   | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 24 --mode optim` |
* "3bas.2X-AveOrb.equil/" directory
  |                       |                                                                            |
  | --------------------- | -------------------------------------------------------------------------- |
  | Article's config name | o27                                                                        |
  | GEOCLIM run name      | .90Ma-3bas-AveOrb.equil                                                    |
  | Description           | ...                                                                        |
  | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command   | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 28 --mode optim` |
* "3bas.Arct-3.2X-AveOrb.equil/" directory
  |                       |                                                                            |
  | --------------------- | -------------------------------------------------------------------------- |
  | Article's config name | o28                                                                        |
  | GEOCLIM run name      | .90Ma-3bas-Arct3-AveOrb.equil                                              |
  | Description           | ...                                                                        |
  | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command   | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 29 --mode optim` |
* "3bas.Arct-3.AveOrb.equil/" directory
  |                       |                                                                            |
  | --------------------- | -------------------------------------------------------------------------- |
  | Article's config name | o28                                                                        |
  | GEOCLIM run name      | .90Ma-3bas-Arct3-AveOrb.deg.equil                                          |
  | Description           | ...                                                                        |
  | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
  | Compilation command   | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --mode optim` |
* "3bas.Arct-3.Laskar/" directory
  |                       |                                                                                                                      |
  | --------------------- | -------------------------------------------------------------------------------------------------------------------- |
  | Article's config name | o28                                                                                                                  |
  | GEOCLIM run name      | .90Ma-3bas-Arct3.Laskar                                                                                              |
  | Description           | ...                                                                                                                  |
  | Compilation command   | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --clim-param 2,4,2 --param-periods ,360., --mode optim` |
* "3bas.Arct-3.epol.2X-AveOrb.equil/" directory
  * "IO\_CONDITIONS" file
    |                       |                                                                            |
    | --------------------- | -------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                       |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-AveOrb.2X.equil                                      |
    | Description           | ...                                                                        |
    | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command   | `./build_GEOCLIM --compset default --res 1,96,96 --nbasin 29 --mode optim` |
  * "IO\_CONDITIONS.APx0.5" file
    |                       |                                                                            |
    | --------------------- | -------------------------------------------------------------------------- |
    | Article's config name | o28'-APx0.5                                                                |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.5-AveOrb.2X.equil                               |
    | Description           | ...                                                                        |
    | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command   | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0.25" file
    |                       |                                                                            |
    | --------------------- | -------------------------------------------------------------------------- |
    | Article's config name | o28'-APx0.25                                                               |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25-AveOrb.2X.equil                              |
    | Description           | ...                                                                        |
    | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command   | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0" file
    |                       |                                                                            |
    | --------------------- | -------------------------------------------------------------------------- |
    | Article's config name | o28'-APx0                                                                  |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0-AveOrb.2X.equil                                 |
    | Description           | ...                                                                        |
    | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command   | *same as previous case*                                                    |
* "3bas.Arct-3.epol.AveOrb.equil/" directory
  * "IO\_CONDITIONS" file
    |                       |                                                                            |
    | --------------------- | -------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                       |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-AveOrb.deg.equil                                     |
    | Description           | ...                                                                        |
    | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command   | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --mode optim` |
  * "IO\_CONDITIONS.APx0.5" file
    |                       |                                                                            |
    | --------------------- | -------------------------------------------------------------------------- |
    | Article's config name | o28'-APx0.5                                                                |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.5-AveOrb.deg.equil                              |
    | Description           | ...                                                                        |
    | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command   | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0.25" file
    |                       |                                                                            |
    | --------------------- | -------------------------------------------------------------------------- |
    | Article's config name | o28'-APx0.25                                                               |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25-AveOrb.deg.equil                             |
    | Description           | ...                                                                        |
    | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command   | *same as previous case*                                                    |
  * "IO\_CONDITIONS.APx0" file
    |                       |                                                                            |
    | --------------------- | -------------------------------------------------------------------------- |
    | Article's config name | o28'-APx0                                                                  |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0-AveOrb.deg.equil                                |
    | Description           | ...                                                                        |
    | Additional modif      | Set `scaling_factor = 1d-3` in `source/dynsoil_physical_parameters.f90`    | 
    | Compilation command   | *same as previous case*                                                    |
* "3bas.Arct-3.epol.Laskar/" directory
  * "IO\_CONDITIONS" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol.Laskar                                                                                              |
    | Description           | ...                                                                                                                       |
    | Compilation command   | `./build_GEOCLIM --compset default --res 2,96,96 --nbasin 29 --clim-param 2,4,2 --param-periods ,360., --mode optim`      |
  * "IO\_CONDITIONS.APx0.5" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.5.Laskar                                                                                       |
    | Description           | ...                                                                                                                       |
    | Compilation command   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | Article Table 8 name  | "All processes"                                                                                                           |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.Laskar                                                                                         |
    | Description           | ...                                                                                                                       |
    | Compilation command   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25.Laskar                                                                                      |
    | Description           | ...                                                                                                                       |
    | Compilation command   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceXchg" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | Article Table 8 name  | -                                                                                                                         |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceXchg                                                                              |
    | Description           | ...                                                                                                                       |
    | Additional modif      | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 755 `call creades_noX(time)`     |
    | Compilation command   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTemp" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | Article Table 8 name  | -                                                                                                                         |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTemp                                                                              |
    | Description           | ...                                                                                                                       |
    | Additional modif      | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 754 `call creades_noT(time)`     |
    | Compilation command   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteCntFlx" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | Article Table 8 name  | -                                                                                                                         |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-csteCntFlx                                                                              |
    | Description           | ...                                                                                                                       |
    | Additional modif      | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 753 `call creades_noWth(time)`   |
    | Compilation command   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTX" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | Article Table 8 name  | "Cont. fluxes"                                                                                                            |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTX                                                                                |
    | Description           | ...                                                                                                                       |
    | Additional modif      | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 756 `call creades_noWth(time)`   |
    | Compilation command   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceXCntFlx" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | Article Table 8 name  | "Oce. temperature"                                                                                                        |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceXchgCntFlx                                                                        |
    | Description           | ...                                                                                                                       |
    | Additional modif      | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 757 `call creades_noXWth(time)`  |
    | Compilation command   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTCntFlx" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | Article Table 8 name  | "Oce. exchanges"                                                                                                          |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTempCntFlx                                                                        |
    | Description           | ...                                                                                                                       |
    | Additional modif      | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 758 `call creades_noTWth(time)`  |
    | Compilation command   | *same as previous case*                                                                                                   |
  * "IO\_CONDITIONS.APx0.25.csteOceTXCntFlx" file
    |                       |                                                                                                                           |
    | --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
    | Article's config name | o28'                                                                                                                      |
    | Article Table 8 name  | "None"                                                                                                                    |
    | GEOCLIM run name      | .90Ma-3bas-Arct3-epol-APx0.25.Lsk-cstOceTXCntFlx                                                                          |
    | Description           | ...                                                                                                                       |
    | Additional modif      | In `source/geoclim_mainprog.f`: comment line 751 `call creades(time)` and uncomment line 760 `call creades_noTXWth(time)` |
    | Compilation command   | *same as previous case*                                                                                                   |


### **How to reproduce the simulations**

Here are the steps to follow to reproduce the 90Ma simulations presented in Maffre et al. (submitted to GMD)

0. Make sure the code of GEOCLIM compiles and runs properly

    FOr this purpose, the bash script `make_test`, generating test-runs, is available on the branch "main" of current Github repository.

1. Gather the needed inputs from the PANGAEA archive

    * IPSL-CM5A2 netCDF inputs: ...
    * DynSoil restarts: ...

    If desired, the COMBINE boundary conditions files can be remade with `preproc/BC/IPSL_90Ma_all.py`

2. Modify the source code if needed.

    The needed source code modifications are indicated, for each simulation case, in previous section "Summary of GEOCLIM simulations".

    * Spin-up equilibrium run
    * Sensitivity experiments

    **Important:**
    * When compiling the code for new simulation case, do not forget to undo the source code modification you have made.
    * Whenever the source code is modified, or put back as original, the compilation command must be re-run, and will erase the former executable.

3. Compile the code

    Use the compilation command indicated in previous section "Summary of GEOCLIM simulations" for each specific simulation case.

4. Configure GEOCLIM for the simulation

    Copy the 2 configuration files (one "IO_CONDITION..." and one "cond_p20...") of the desired simulation case, in the corresponding
    subdirectory in `config/90Ma_templates/` (see previous section "Summary of GEOCLIM simulations"), and **replace** the files
    `config/IO_CONDITIONS` and `config/cond_p20.dat`

5. Run GEOCLIM

    Run the executable file **corresponding to the desired simulation case** (i.e., the executable generated by the compilation command
    indicated in previous section "Summary of GEOCLIM simulations").

    The calculation performance is about 30 minutes per million years of simulation.
    The 10-Myr long 90Ma simulations may take up to 24 hours to complete.

