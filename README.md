# GEOCLIM7.0




## GEOCLIM model in a nutshell
GEOCLIM is cluster of models, more or less adaptable, computing geochemical cycles of several species (C, O...) at geological
timescale.
The core of the model is an ocean-atmosphere chemistry model (advection-reaction) COMBINE. Ocean an atmosphere are discretized
in 10 (or more) reservoirs (boxes).
This core is closely associated to an early diagenesis module computing the "output" fluxes (burial of elements in marine sediments)
for each box.
It is also associated to an continental weathering module computing the "input" fluxes. This weathering module is spatially-resolved
(using a geographic mesh grid), its resolution is adaptable, and several options exist for the silicate weathering part.
This triplet is indirectly coupled to a climate model (GCM).
Climate simulations must be run before using GEOCLIM, for a range of CO2 levels. Any climate model can be used, as long as it computes
surface air temperature and continental runoff. Oceanic temperature from the climate model can be used, or parameterized if not
available. Climate fields are then interpolated on the "CO2 dimension", at the current atmospheric CO2 computed by GEOCLIM, and
used to compute continental weathering, and oceanic boxes temperature.
The resolution (ie, the spatial grid) of the continental weathering module must be the same than the GCM.

##### How to run the model:
After downloading the present repository, type `./make_test testname` (testname being one of "ERA5", "GFDL", "CESM", "paleo" and
"ascii"). This command will compile and execute a short GEOCLIM run and compare the output to a reference template.
This allows to verify that the compilation and execution of the model are performed without error, and yield the same
results than reference runs. If not, the command should tell what type of error was encountered (see also section *Frequent issues*
at the end of this file). You could try the 5 tests to make sure everything works as excepted.

If the tests are conclusive, follow those step to create your run:
1. **Gather the input files you need**.  
The standard input files are **climate model netCDF outputs: ATMOSPHERE, LAND, and OCEAN files**.
The alternative way ("old" GEOCLIM setup) is to use ASCII input files, but not all the recent options are compatible with ASCII inputs.
2. **Generate the boundary conditions files with "preproc/BC/BC_generator.py" and "preproc/BC/basinmap/basinmap_editor.py"**.
The latter is only needed for explicit land-to-ocean routing.  
This script generates *from climate model netCDF OCEAN files* a bunch of text files with all the information needed by GEOCLIM 
for the ocean box model. See "preproc/BC/examples\_oceanic\_input.py".
This step is not strictly necessary, since those text files can be created manually, and they are options to use simplified configurations,
like parametric ocean temperature and circulation.
3. **Configure the run by editing the files "config/IO_CONDITIONS", "config/GCM_input_conditions" and "config/cond_p20.dat"**.  
Indicate the name of the run, initial conditions, input files (ATMOSPHERE and LAND netCDF files, oceanic text files created by
"BC\_generator.py"), solver parameters...
4. **Compile the code with `build_GEOCLIM`**.  
Specify the model set of components, resolution... Try `./build_GEOCLIM --help` for more information.
This command create an executable file "geoclim\_$$$.exe" in "executable/"
5. **Run the model with `executable/geoclim_$$$.exe`** (name of the executable created at the compilation step).  
Alternatively, and if you are running the model on a cluster, you can use the files in "job/" to submit the run as a batch process.
"job/submit\_chain.sh" is a script designed for submitting a series of runs, each new one starting from the end of the previous one.
See also section *Multiple runs and job submission* in *Run GEOCLIM*.

#### Required software
* Fortran compiler (note: the Makefile is configured for gfortran, ifort or pgfortran)
* netCDF-Fortran library (see section "How to install netCDF-Fortran library")

#### Required software, except for "old" GEOCLIM setup
* Python. The main boundary condition generating scripts (highly recommended) is written in Python.
Some scripts for pre-processing and visualization scripts are in Python.
* "netcdf\_visualization" package (https://github.com/piermafrost/netcdf_visualization).
Only needed by the specific boundary condition generating script "preproc/BC/basinmap/basinmap\_editor.py".
Those scripts use the following packages:
    * numpy
    * netCDF4
    * matplotlib (plotting scripts)
* 'make' software (installed by default on all UNIX/Mac OS)




## Updates from GEOCLIM6:

... TO BE WRITTEN ...

### Summary of former updates (from GEOCLIM5 to GEOCLIM6):
* Flexibilization of inputs/outputs (netCDF), units and axis recognition (for inputs), use of "namelist" syntax (Fortran).
* Other code simplifications, remove redundant (multiple locations) or hard-coded definitions.
* Implementation of climatic parameters, multilinear interpolation, and simplified sulfur cycle.




## Templates
A couple of template GEOCLIM runs are defined. They consist of a set of input data and corresponding configuration files.
The script `make_test` does the complete configuration (pre- and post-compilation), compiles and runs the desired template
(the pre-existing configuration and source files will be restored when the run is complete).
The expected GEOCLIM outputs for those templates are available in "OUTPUT/templates/", so you can check that you obtain
the same results.
Type `./make_test` for more information.




## Configure GEOCLIM

GEOCLIM is designed to be adaptable the largest possible range of paleogeographic configurations, and to get the information it needs
(forcing fields) from climate model (GCM) outputs.

Hence, the first step of GEOCLIM configuration is to define:
* The number of ocean boxes, and their characteristics: temperature, water fluxes between them, positions (with respect to continents,
latitude, depth), volumes and horizontal areas.
* The continental resolution. Continents are represented on a 2D grid, traditionnally (but not necessarily) longitude x latitude.
* The routing of continental fluxes to ocean boxes, in case there is more than 1 "coastal" box (1 coastal box is the traditional
  configuration).
* The routing of marine sediments on seafloor (meant to represent turbidite fluxes).
* The number of lithological classes, and their distribution on continents.
* The parameters ("dimensions") along which climate fields can vary, i.e., CO2 (mandatory) and optional climatic parameters (e.g., orbital
  parameters).
  Climate fields are read by the code from netCDF or text files, and the code expect, for each climate variable (temperature, runoff...),
  **one field per combination of CO2 and climatic parameters**.

NOTE: the script 'configure.sh' has not yet been updated to include the new features of GEOCLIM.

##### Note on climatology average
The climate fields (land temperature and runoff, ocean temperature and circulation) considered by GEOCLIM model, are supposed to be
**annual mean** variables, **averaged over many years** (*climatology average*. eg, 30 years).
The model does not take into account the seasonal cycle, or any other "short-term" climatic variability.


### Input data

#### Ocean boundary conditions

##### Description
The ocean boundary conditions files are conventionally stored in individual repetories in "INPUT/COMBINE/".
There are 13 files (11 mandatory).
Their standard names are (with "N" the number of GEOCLIM boxes, and "K" the number of CO2 and climatic parameters combinations):
* "indice\_epicont.dat"     : N lines of 1 element ("1" or "0"). 1 line per GEOCLIM box, indicating if each box is "coastal" (1) or
                              "open-ocean" (0)
* "indice\_polar.dat"       : Similarly, indicating if each box is "polar" (1) or not (0)
* "indice\_surface.dat"     : Similarly, indicating if each box is at top of water column (1) or not (0)
* "thermocline.dat"         : Similarly, indicating if each box is in the middle of water column (1) or not (0)
* "indice\_deep.dat"        : Similarly, indicating if each box is at the bottom of water column (1) or not (0)
* "apport\_ct.dat"          : Similarly, indicating if each box receives continental fluxes (1) or not (0)
* "oce\_surf.dat"           : Similarly, indicating the horizontal area of each box (in 1e15m2)
* "surf\_sedi.dat"          : Similarly, indicating for each box its horizontal area (in 1e15m2) intercepting the seafloor
* "oce\_vol.dat"            : Similarly, indicating the volume of each box (in 1e15m3, and *always* 1e-15 for the atmosphere!)
* "press\_box.dat"          : N lines (1 per box) of 2 elements, indicating the depth of the bottom of each box (in km), and its mean
                              pressure (in bar)
* "box\_sf\_bnd\_len.dat"   : [OPTIONAL] N lines of N elements. The j-th column and i-th line indicates the length of the horizontal
                              boundary between the j-th and and i-th boxes (is symmetric by construction). Could be in any length units,
                              since it will be normalized.
* "GCM\_oceanic\_temp.dat"  : [OPTIONAL] K lines of N elements. Each line corresponds to 1 CO2/climatic parameter combination.
                              The 1st column is null (was historically CO2 level), the other columns are the mean temperature of the N-1
                              oceanic boxes (in °C). The atmospheric (last) box is omitted!
* "exchange\_2.dat"         : K blocks of N lines of N elements. Each block corresponds to 1 CO2/climatic parameter combination.
                              On each block, the j-th column and i-th line indicates the flux of water from the j-th to the i-th box (in Sv).
                              The blocks are not symmetric, but *the sum of each line must equals the sum of each corresponding column*. 
                              Alternatively, there can be a single N*N block, that the code interprets as an unvarying oceanic circulation.

###### Implicit rules of GEOCLIM boxes
The horizontal definition of boxes is entirely customizable, but the vertical definition is subject to some rules:
* The sequence of the GEOCLIM boxes must be organized in "blocks" of water columns: i.e., the k-th box is either
  * the box immediately below the k-th box in the water column
  * the surface box of a different column, if the k-1-th box is at the bottom of the water column (or if k = 1).
* As a consequence, there cannot be any "branching" within a water column: only 1 box is located below a given box
  (or 0 if it's the bottom of the water column).
* There is no retriction *a priori* on the number of vertical levels within each water column, but there has to be
**exactly 4 different vertical levels** that intercept the seafloor. 2 for the coastal boxes, 2 for the open oceans.
* There can only be 1 atmospheric box, and it must be the last box.

##### Generate ocean boundary conditions from climate model netCDF outputs
The Python function "oce\_output\_2\_GEOCLIM\_BC" in the file "preproc/BC/BC\_generator.py" is designed for generating the 13 GEOCLIM ocean
boundary condition files **from netCDF "ocean" outputs of any climate model**.
In other words, the function converts "raw" 3-D ocean fields from climate model outputs into GEOCLIM box-integrated values.

"oce\_output\_2\_GEOCLIM\_BC" takes as input:
* The variables defining the ocean grid: longitude, latitude, depth, depth bounds (top and bottom of cells),
  x-weighting, y-weighting, area-weighting, z-weighting and bathymetry.
  Not all these variables are required, e.g., depth-weighting can be computed if depth bounds are given.
* 3-D oceanic temperature field, and water velocity or flux (U, V and W).
  Note that the standard way for W is to deduce it from U,V divergence. 
* Information about how to split the "ocean": cutting depths and latitudes, and (optional) horizontal mask to
  define regional basins (e.g., Atlantic/Pacific), given as a netCDF variable.

"oce\_output\_2\_GEOCLIM\_BC" can read inputs from multiple netCDF files, up to 5 different (for grid definition, temperature, U, V and W).
Since temperature and oceanic circulation are climate-dependent, "oce\_output\_2\_GEOCLIM\_BC" actually expects *lists* of netCDF files,
each corresponding to a combination of CO2 and climatic parameters (see previous section "Description").
If temperature, U, V and W are stored in different files, the 4 lists of files must be given in the same order, order in which the BC files
"GCM\_oceanic\_temp.dat" and "exchange\_2.dat" will be written.  
"oce\_output\_2\_GEOCLIM\_BC" accepts time-dependent inputs. In other words, it can be given 12-month fields, and will compute the annual
mean.

"oce\_output\_2\_GEOCLIM\_BC" can also save in a netCDF file the mask of the actually defined GEOCLIM oceanic boxes *on the "native" 3D
oceanic grid* (option `export_mask_fields=True`).
This mask (saved in this way) can be used to generate the land-to-ocean routing map ("preproc/BC/basinmap/basinmap\_editor.py", see next
section). 

More information about the function can be found in its docstring (`help(oce_output_2_GEOCLIM_BC)` in a Python console).
Several examples of use can also be found in "preproc/BC/examples\_oceanic\_input.py".

#### Land-to-ocean fluxes routing
The explicit routing of continental fluxes to oceanic boxes *is only necessary if there are more than 1 epicontinental boxes* (which is not
the case on the "traditional" GEOCLIM set-up).
Even in case of multiple epicontinental boxes, in the absence of information about the routing, there is the option `uniform_routing=.true.`
in "config/IO\_CONDITIONS", to uniformaly distribute the continental fluxes to all epicontinental boxes, proportionally to the area of
these boxes.
Otherwise, GEOCLIM expect, as input, a land-to-ocean routing map. I.e., a map of drainage basin, indicating, on each point of the
continental grid, the # of the GEOCLIM box in which the pixel's fluxes will be delivered.

The Python function "make_routing_map", in the file "preproc/BC/basinmap/basinmap\_editor.py", generates that land-to-ocean routing map,
directly from climate model netCDF outputs, and from the ocean basin mask definition (on the oceanic grid), that is genetated by
"preproc/BC/BC\_generator.py". It is computed in 2 steps:
1. For each continental point, determining the *closest oceanic point*, and assigning its oceanic basin # to the continental point
2. Launching the interactive editor (from "netcdf\_visualization" package) for the user to manually edit the generated drainge map

This interactive edition is almost always necessary for two reasons:
* The land routing must be towards GEOCLIM *epicontinental* boxes. However, on the oceanic-grid basin mask, the grid cells at the edges of
  continents can often belong to *open-ocean* boxes (not epicontinental), if the seafloor drops directly below the epicontinental depth limit,
  at the border of some continent. So, by finding the closest oceanic point, the algorithm will route land fluxes toward both open-ocean and
  epicontinental boxes. This need to be manually corrected.
* The routing algorithm cannot use the information about pre-defined river basins (in the climate model), since there is no generic way to
  indicate which river flows in which oceanic basin.
  For this reason, this function can takes as input argument the map of pre-defined river basins, and river outlets map), for purpose of
  plotting and editing. The user should then "hand-pick" each basin and edit the routing map. This can be easily done with the selecting
  option "click on value" from the netCDF editor. The selected area is kept when switching from one variable to another in the editor.

More practical information can be found in the docstrings of "preproc/BC/basinmap/basinmap\_editor.py".

#### Other boundary conditions

##### Atmospheric and continental climate

GEOCLIM Fortran code directly reads climate model's netCDF outputs from 'atmosphere' module (grid area and world-wide 2-meter temperature)
and 'land' module (land area, temperature, and runoff).
Those 2-D variables need to be *annual mean* climatology (e.g., not month-by-month climatology).
It is required that the 'atmosphere' and 'land' variables are defined on the same grid, that can be (and often is) different from the
'ocean' grid.

##### Lithology

The lithology map should be given as a netCDF 3-D variable indicating the fraction of each lithological class on each cell of the
**2-D land/atmosphere grid**.
The first two dimensions of that variable are therefore the horizontal dimensions of the land/atmosphere grid, the third dimension is
the "lithology" dimension.

Alternatively, a uniform lithology distribution can be specified, without providing a map (see next section "Input/Output user interface")

###### Note on lithological classes
The number of lithological classes is defined by the parameter 'nlitho' in "source/shape.inc", and is conventionally configured with
`build_GEOCLIM`.
The model is traditionally parameterized for 6 classes. This can be modified, as long as it stays consistent with **two implicit rules**:
* The last lithological class is "carbonate".
* There is one (and only one) lithological class corresponding to "basalts", indicated by the parameter 'BASALT\_LITHO\_NUM' in
  "source/constante.f90".

If you change the number of lithological classes, you also need to update the lithology-dependent parameters in "source/constante.f90"
and "source/dynsoil\_physical\_parameters.f90".
Note that even though 'CaMg\_rock' and DynSoil parameters are *not used* if GEOCLIM is not coupled with DynSoil, they must be defined,
consistently with 'nlitho', or the compilation will fail. You may simply fill them with zeros.

##### Slope

Topographic slope is only used if DynSoil module is activated, to compute erosion rate.
It should be given as a netCDF 2-D variable, defined  on  **the same land/atmosphere grid**.
It represents the average of the absolute slope of the high-resolution topography on each grid cell.

Even though this input field is critical (it largely determines the erosion rates, that controls many other geochemical fluxes), there is
no pre-defined way to generate the slope field of paleogeographic reconstruction.
What is often done is to use empirical correlation between topography and mean absolute slope.

##### Time-series of climatic parameters

A basic text file, each line contains the values of the user-defined climatic parameters (up to 5 are possible) at a given time step.
The time interval at which the lines of this file are read is indicated in "config.cond\_p20.dat" (see next section "Input/Output user
interface").

#### Note on physical units

**In all the netCDF files given to either boundary condition generating Python scripts, or GEOCLIM Fortan code,
the physical units of every variable should be indicated in its netCDF attribute "units"**.
Specific functions are implemented to identify physical units and perform the needed conversions.

If the units are not recognized, the default units will be assumed, but the user will be asked interactively to validate it.
New units and conversion factors can easily be added in "source/physical\_units.f90" and "preproc/BC/units.py"
(see section *defining new physical units* in *Further information*).

#### Tips for specific climate models

##### IPSL
* "Total" continental runoff (= precipipation minus evaporation) is the sum of 'land' (sechiba) outputs variables "runoff" and "drainage".
* The map of continental river basins is "basinmap", and the map of river outlets is "nbrivers" (both in 'land' sechiba outputs).
* Grid cell area in the 'atmosphere' outputs is the variable "aire", whose true units is "m2".
  **Its "units" attribute "-" in the netCDF file needs to be manually modified to "m2"** ("-" is interpreted as "cell fraction" by GEOCLIM).
* The grid of 'ocean' outputs has some peculiarities: duplicated cells, connection of North Pole fold...
  They are taken into account by specifying "special\_wrap='ORCA...'" in the function "oce\_output\_2\_GEOCLIM\_BC" of
  "preproc/BC/BC\_generator.py". One must also specify if the North Pole pivot is on a "T" point or "F" point.

##### FOAM
* The total continental runoff variable is "RNF", stored in the "coupler" output netCDF file.
  Its units ("m") means "meter accumulated during a time step", which is 30 minutes.
  Therefore **this variable needs to be manually converted in units understandable by GEOCLIM, and its "units" attribute properly updated**.
* The land fraction of grid cells ("land-sea mask") is in variable "ORO" of "coupler" netCDF output, but mixed with sea-ice fraction
  (values between 0 and 2). Make sure to get rid of sea-ice values to retrieve the "true" land-sea mask information.


### Initial conditions
GEOCLIM can automatically generate initial conditions ("coldstart", to indicate in "config/IO\_CONDITION").
When a run is over, it writes "restart" files, that can be used as initial conditions for branching new runs.


### Input/Output user interface

The GEOCLIM input/output interface is managed with 3 files (2 mandatory):
* "config/IO\_CONDITIONS":
  Main IO interface, provides the paths of all the boundary condition files (except those in "GCM\_input\_conditions"),
  initial condition files, and the name of the output files and which variables to output (under which name).
* "config/GCM\_input\_conditions":
  [OPTIONAL] You can indicate in "config/IO\_CONDITIONS" to read the land/atmosphere boundary conditions in this file.
* "config/cond\_p20.dat":
  States the physical and numerical parameters: solver time steps, output writing frequency, duration of run, when to generate restarts,
  acceleration parameters, and volcanic degassing (CO2, SO4, Trapp setting...)

Those files are self-describing. Examples can be found in config/templates/

##### Information given in "IO\_CONDITIONS"

###### "MAIN IO CONDITIONS"
Name of current run, directory outputs will be written, name of the 2nd config file ("cond\_p20.dat"), and killing file ("deathnote.txt")

###### "INITIALIZATION FILES (ie, restarts)"
Path of initial condition files (if needed), and, for DynSoil initial condition, name of netCDF variables.  
Here, one may specify `COMB_init_mode='coldstart'` and `init_mode='startup:eq'` (DynSoil) to use internally-generated initial conditions.

###### "INPUT FILES"
* "MAIN CONTINENTAL INPUTS":
  Name of "GCM\_input\_condition" file.
  Alternatively, the option `cont_input_mode='ascii'` allows to read area, land area, land temperature and runoff from ascii files.
  In that case, areas must be in 1e6 km2, temperature in °C and runoff in cm/y, there is no possibility to define global temperature,
  and to use climatic parameters (other than CO2). See ascii input examples in INPUT/FOAM/CTRL\_48x40/.
* "LITHOLOGY:
  Path of lithology netCDF file + name of variable, or prescribed uniform lithology fractions.
* "CONTINENTAL FLUX ROUTING":
  Routing scheme and path of (potential) land routing map netCDF file + name of variable.
  Alternative option: `uniform_routing=.true.` (to distribute continental fluxes uniformly on all epicontinental boxes, proportionally to
  their area).
* "SLOPE:
  Path of (optional) slope netCDF file + name of variable
* "OTHER CONTINENTAL INPUTS":
  Vegetation map file (archaism) and climatic parameter time-series text file (optional)
* "OCEANIC INPUTS":
  The 13 ocean boundary condition files described in section "Ocean boundary conditions" +
  the input mode for the 2 optional ones ('ocean\_temp' and 'sedim\_transport').  
  `ocean_temp_mode='parametric'` if for using the parameterization of oceanic box temperature (if that information is missing), while
  `ocean_temp_mode='file'` is for reading it in the input file.   
  `sedim_transport_mode='uniform'` is for uniformly routing seafloor sediments from box to box (but only from shallow to deep), while
  `sedim_transport_mode='neighbour'` is for routing sediment from one box to the connected boxes (still from shallow to deep),
   proportionally to the length of the boundary between the boxes. In this latter case, the input file "box\_sf\_bnd\_len.dat" is needed.

###### "OUTPUT CONDITIONS"
The name of netCDF files, name and definition of variables (dimensions, attributes...) of GEOCLIM ouputs, split in 4 modules:
* "COMBINE OUTPUTS": main geochemical variables
* "GEOGRAPHIC OUTPUTS": continental variables (weathering fields)
* "DYNSOIL OUTPUTS": continental variables of DynSoil module (erosion, regolith variables...)

This section usually doesn't need to be modified, save for indicating which variable to write in the ouputs, or not to (`writevar(*)=FALSE`)

###### "names of restart files"
Definition of restart offline files written by GEOCLIM, usually at the end of a run.
For a regular use of GEOCLIM, this section does not need to be modified. 

##### Information given in "GCM\_input\_conditions"

###### "Climatic parameters"
The list of combinations of CO2 and climatic parameters, **entered in the same order than the oceanic boundary condition files
"GCM_oceanic_temp.dat" and "exchange_2.dat"**

###### "CONTINENTAL ANNUAL CLIMATOLOGY"
Path of netCDF files, name of dimension and variables for:
* land/atmosphere grid cell area
* land area or fraction of grid cells
* land temperature
* land runoff
* world-wide 2-m temperature (usually in 'atmosphere' file). OPTIONAL, used to compute global mean temperature.

1 netCDF file is expected per combination of CO2/climatic parameters, both for the 'land' and 'atmosphere' variables.
The order in which those files are entered must be consistent with the given list of combinations of CO2/climatic parameters.

##### Information given in "cond\_p20.dat"
Values of solid Earth degassing fluxes, flags for several GEOCLIM modules, solver parameters (starting and ending time, timestep),
acceleration factors, time interval for reading the climatic parameters time-series, time intervals for writting the outputs, and
time at which writting the restart files.


### "Manual" compilation configuration
This part is unnecessary if you compile GEOCLIM with `build_GEOCLIM` (which is the recommended way).

If you wish not to do so, they are additional configuration steps (that `build_GEOCLIM` automatically perform) you should follow:
* Indicate, in "source/coupler.inc"
  * which module to use: 'coupling\_dynsoil', 'use\_dynsoil\_steady\_state', 'coupling\_ecogeo' (.true. or .false.)
  * the type of CO2 interpolation ("log" or "linear")
  * if using climatic parameters, whether or not loop back to the beginning ('climparam\_loop') or stop the run ('climparam\_kill') when
    reaching the end of the climatic parameters time-series
* Indicate, in "source/shape.inc":
  * the total number of ocean boxes + atmosphere 'nbasin'
  * the number of CO2 levels 'nclimber'
  * the geographic resolution 'nlon' and 'nlat'
  * the number of lithological classes 'nlitho'
  * the number of DynSoil vertical levels 'nDSlev'. It must always be defined, but if "coupling\_dynsoil == .false.", its value is not used.
  * if using climatic parameters, their number 'nclimparam', the length of the dimensions associated with each climatic parameter
    'len\_p\*', and their potential value's period ('p\*\_period', keep it as "0" if the parameter is not periodic). 
    Note that if a climatic parameter is periodic, its dimension length 'len\_p\*' must here be increased by 1 with respect to the number
    of inputs.

IMPORTANT: all those parameters must be consistent with the initialization and forcing files (initial ocean chemistry, climate files...)




## Compile GEOCLIM

The command `build_GEOCLIM` performs the "compilation configuration" (i.e., editing the needed source files) and compile the code,
using the Makefile in "source/".
This command supports a large number of options, including:
* The set of components (modules) to use
* The model's resolution (number of ocean and atmosphere boxes, continental resolution, number of lithological classes,
  DynSoil vertical resolution, and the definition of climatic parameters).

IMPORTANT: all the information passed to `build_GEOCLIM` must be consistent with the input/output interface files (see above).  
Type `./build_GEOCLIM --help` to get detailed information on how to use it.

If you wish to compile without `build_GEOCLIM` (to recompile the code, for instance), here are the steps to follow,
**assuming that the "compilation configuration" is done** (either `build_GEOCLIM` was previously invoked, or it was manually done,
see previous section "'Manual' compilation configuration").

#### With Makefile:

Note that if you use the command `build_GEOCLIM`, **it will tell you the "make" command that was used** (make + arguments).
This command is also **saved in the file "GEOCLIM_environment"**.

Go to 'source/'.

If you use gfortran compiler (default one), just type `make`, the executable 'geoclim.exe' will be created.
If you use a different compiler, you must specify it by typing `make FC=your_compiler`. Note that the configuration is
only made for gfortran, ifort and pgfortran. If you use a different one, the compilation options will not be defined.
See following paragraph.

###### Details of Makefile:
The "standard" options to customize the Makefile are:

`make [FC=...] [MODE=...] [NCPATH=...] [execut=...] [main_flags="..."] [NETCDF_FLAGS="..."] [FFLAGS="..."]`

All those arguments are optional (hence the []).

* FC=...: Fortran compiler. The Makefile is configured for 'gfortran', 'ifort' or 'pgfortran'
  Default: gfortran
* MODE=...: sets the compilation options. 3 options are accepted
    * 'standard': (default), standard check options.
    * 'debug': extra debugging options
    * 'optim': with optimization flags (and less debugging options, like traceback)
* execut=...: sets the name of the created executable file. Default is 'geoclim.exe'
* NCPATH=...: States the path of the directory where the netCDF-Fortran library is installed. This assumes that the only
needed options are "-I$NCPATH/lib", "-L$NCPATH/include" and "-lnetcdf -lnetcdff". If not, use NETCDF\_FLAGS option.
By default, NCPATH is '/usr', but you may have installed your netCDF library elsewhere (/usr/local, /usr/local/netcdf, ...)
It can be inquired with `nc-config --prefix`. The command `build_GEOCLIM` first tries this to get it. 
* NETCDF\_FLAGS="...": Override the netCDF flags.
* main\_flags="...": Override the main compilation flags (all but the netCDF flags).
  The variable 'MODE' becomes useless. Useful if you use a different compiler whose options are not configured in the Makefile.
* FFLAGS="...": Override all compilation flags. The variables 'MODE' and 'ncpath' become useless. Useful if you use a
  different compiler whose options are not configured in the Makefile and who does not support "-I", "-L", "-lnetcdf" or
  "-lnetcdff" options.

In any case, you can check the compilation command that will be used by doing:

`make echo [all the options you want]`

#### Without Makefile
`make` is not necessary to compile the code. If you do not use it, make sure that:
* The file 'source/path.inc' exists and contains the line `character(len=*), parameter:: geoclim_path = "..."`
  Where '...' is the path of the GEOCLIM root directory (this file is edited automatically by the Makefile)
* Your compilation command uses the netCDF options. Usually, it must have the options `-I/usr/include -L/usr/lib -lnetcdf -lnetcdff`
  (/usr/lib and /usr/include are the 2 directories where the netCDF library is commonly installed, but it may be elsewhere).
* You use the 'gnu' or 'Fortran 2003' standard for *all* source files (usually `-free` or `-ffree-from`, possibly `-std=gnu` or `-std=f2003`)
* To make sure that the executable is up to date with the source files, do `rm -f *.o *.mod *__gen_mod*` before compiling (or `make clean`)
* Compile 4 times (there are 4 levels of nested subroutine calls). It is normal to get many errors on the first 3 times this is executed




## Run GEOCLIM


### Executable files
`build_GEOCLIM` put the executable file in "executable/". The other methods let it in "source/". It can be run from any
directory, as all the paths in the code are absolute paths.


### Inputs error handling
GEOCLIM performs tests on the input files before the "main" execution. They are 5 types:

1. axis mismatch between the input files (for instance, shifted longitude)
2. missing values on continental pixels (continental pixels are defined by the "land area" input variable)
3. invalid value for runoff (negative) and slope (negative or null)
4. sum of all lithology classes differs from 1
5. units not recognized in netCDF inputs

By default, the executable interactively asks the user what to do when an error is encountered. This can be problematic when run
as a batch process on a cluster (with no interactive interface). It is possible to pass 5 arguments to the executable, as follows:

`./executable i1 i2 i3 i4 i5`

where i1...i5 are integer numbers, between -1 and 3, and correspond respectively to the 5 kind of errors above-mentioned.
* -1 means 'ask the user interactively' (default)
*  0 means 'abort the execution'
*  1 means 'remove the problematic pixels' (not possible for axis mismatch or units not recognized)
*  2 means 'ignore the issue and continue execution without any change'
*  3 means 'replace the invalid value' (only possible for runoff and slope)

For instance, I recommend `./geoclim.exe 0 1 3 0 0`, or if you are sure of your axis, lithology mask and units `./geoclim.exe 2 1 3 2 2`

NOTE: these tests cannot be done with ascii input files (except the negative runoff test), use that format at your own risk!

Also, errors exist that the code cannot handle, for instance, if one of the input file specified in config/IO\_CONDITIONS
does not exist, or if the shape of the netCDF variables does not match the one specified in the code (source/shape.inc).
These will cause the run to crash. The code can, however, handle transposed 2D (x-y) variables in GCM input files, as well as
degenerated (size-1) extra dimensions.


### Output
GEOCLIM writes outputs in netCDF format, in 3 files (see previous section "OUTPUT CONDITIONS")
It is possible to automatically convert those outputs in ascii format, with the option `convert2ascii=.true.` in config/cond\_p20.dat

The frequency at which outputs are written is specified in "config/cond\_p20.dat". There are 3 frequencies: one for COMBINE
output, one for geographic outputs (ie, continental variables) and one for specific DynSoil outputs.


### Restart
Restarts are created in the output directory. It is a good habit to move them to "restart/.../". Automatic launching script,
like `submit_chain.sh` (in job/) will automatically move the restart files in that directory.

The time for restart generation is specified in "config/cond\_p20.dat" (`ageYprint`). By default, this time is set to "-1.", which is
interpreted as "when the run is complete". Note that the restarts are generated 1 time only.


### Killing a run
Sometimes one may want to end a run and create restart files precociously. For instance, if a run has reached the steady-state
sooner than expected and one wants to launch "perturbation runs" from that steady-state.

To do so, simply write the name of the run (as specified at the first uncommented line of "config/IO\_CONDITIONS") in the file
"deathnote.txt". It will cause the run to stop and generate restart files (if they were not already created).
You can put as many run names as you want in the deathnote, one by line. Don't forget to erase the names afterward!


### Multiple runs and job submission
**! THESE SCRIPTS HAVE NOT BEEN UPDATED FROM GEOCLIM6.1 AND ARE CURRENTLY NOT WORKING !**

`submit_chain.sh` (in the repertory job/) is a bash script for automatically launching a series of GEOCLIM runs.
A series of runs are runs that have exactly the same configuration (except for their timesteps, starting and stopping times),
each one starts from the end of the previous one (the very first one starts from the initial condition given by the user).
This is useful for runs with an initial perturbation, requiring a short timestep, but whose long-term evolution (after
the adaptation to the perturbation) can be computed with a longer timestep.
In addition, it offers to possibility to submit the GEOCLIM run as batch processes (jobs), which is required on clusters.
Clusters usually have a time limitation for jobs, which makes the automatic resubmission (series of run) helpful.
Finally, this script provides a security for conflicting access to the configuration files, that is helpful for running
several independent runs in parallel.

Practically, the script `submit_chain.sh` works in pair with a second script (usually, `run_geoclim.sh`). The main script
(submit\_chain.sh) "submits" the second one (either executes it, or submits it with the cluster submission command), that
actually run the geoclim model, and call the first script back when the run is completed.
The main script does all the configuration, and move the restarts. The second is only for running the GEOCLIM executable,
but must be configured for the current cluster (whereas the main one is a bash script meant to be executed directly).

When using 'submit\_chain.sh', the pre-compilation configuration must be done (and the code compiled). If you wants to launch
in parallel several runs that need different pre-compilation configurations, save as many different GEOCLIM executable files.
Here is the list of options that can be customized with 'submit\_chain.sh':
* The name of the run (for a series of runs, suffix '\_1', '\_2'... are automatically added).
* The submission command (cluster-dependent) and the name of the running script (usually, 'run\_geoclim.sh').
* The name of the GEOCLIM executable file.
* GEOCLIM (COMBINE) and DynSoil initial condition.
* Stopping (and restarting) times. Note: The "first" starting time is given by config/cond\_p20.dat, and should normally be 0.
* The different model timesteps and printing timesteps.
* The job log file.
* The name of GEOCLIM main configuration file (normally, config/IO\_CONDITIONS. In case extra configuration customization
is needed. Usually, keeping the default one is sufficient).

The script is designed for parallel runs. It edits the configuration files and ensures there is no conflict.
Once you have submitted one run (series of run), you can safely edit the file 'submit\_chain.sh' and submit a second run (series of
runs). The script will tell you if a run is waiting to access the configuration files.
If you need to do configuration modifications not available in 'submit\_chain.sh', the safest way is to create a new config file
"IO\_CONDITIONS" and to tell 'submit\_chain.sh' to use it (note: the name of the other config files, like cond\_p20.dat, are stated
in the main one). Remember that if you edit any of the configuration files, it will impact **all** the series of run that have
been launched. When all the runs are completed, the original configuration files will be reinstated.


### Special runs

#### Equilibrium (accelerated) run
An "equilibrium run", or accelerated run, is a run whose transient evolution is of no interest because one only wants to get
the geochemical steady-state (for instance, to start perturbation from that steady-state).
In that case, a couple of things can be modified to shorten the time needed to reach the steady-state, without modifying it.

###### Before compilation: 
Only if you are using DynSoil module in its dynamic version, should you decrease the value of 'scaling\_factor' in
"source/dynsoil\_physical\_parameters.f90". The scaling factor controls the inertia of the regolith, and does not affect
its steady-state. 1 is for a normal regolith. In some places, regolith can take millions of years to reach its steady-state.
To shorten that time, set it to 1d-3. You can put a value as close to zero as you want, it will not generate
any numerical instability. However, it will become useless if the evolution time-scale of the regolith is lower than the
model time step.

Do not forget to put the 'scaling\_factor' back at 1 after the run is complete!

Alternatively, you can use the "steady-state" version of Dynsoil. The code will directly compute the analytical steady-state.
This has a lower computation cost, but it has no visible effect, if the asynchronous time step of continental weathering is
high enough.
Note that there will be a slight difference between the analytical steady-state and the numerical one (reached with the
"dynamic" version of DynSoil) simply because of the vertical discretization.

###### After compilation
* Oxygen cycle acceleration: O2 has a residence time of ~8 Myr, so it requires around 20 Myr to reach equilibrium.
An acceleration coefficient can be tuned in 'config/cond\_p20.dat'. Setting it to 100 is enough to bring the equilibration
time down to Carbon residence time. An excessively value will cause the model to crash.
* Sulfur cycle acceleration: Similarly to oxygen, an acceleration coefficient can be tuned in 'config/cond\_p20.dat' to reduce
the time needed for sulfur cycle to reach equilibrium (the residence time of sulfur is ~30 Myr). 100 is a good value.
* Asynchronous coupling with continental weathering: The standard time-step for continental weathering is 25 years, and
100 years for DynSoil module (if activated).
The model spends a significant amount of time on the continental computation, especially at high resolution (1° or less)
and when coupled to DynSoil.
Increasing that time step to 250 years, or 1000 years will hasten the run, only degrading the quality of the transient
evolution. If it is too high, however, it can increase the model time needed to reach the steady-state.
Moreover, the steady-state weathering flux of DynSoil module (in its "dynamic" version) are actually dependent of the DynSoil
timestep, because of numerical accuracy. A longer timestep will result in slightly lower weathering flux (ie, higher equilibrium
CO2). For instance, increasing the timestep from 100 year (default) to 10000 years cause the CO2 to rise by a ~10 ppmv.

#### Fixed CO2 run
This is a special case of model configuration. If there is only 1 CO2 level (nclimber=1), it will be run in "fixed CO2 mode".
This means the amount of CO2 in the atmospheric reservoir will be held constant at the value of the unique CO2 level, whatever
the carbon fluxes. The concentration of the various forms of carbon in the other reservoirs will adjust freely to the atmospheric
concentration (ocean-atmosphere diffusion) and to the carbon fluxes.
In other words, in fixed CO2 mode, the mass balance is not respected for carbon.

This is useful for calibration runs where one wants to hold the atmospheric CO2 constant and adjust the degassing to balance the
silicate weathering flux.

#### Run with locked geochemical cycles
This can be useful if one wants to investigate the behavior of inorganic carbon cycle only, without the feedback of the
other cycles, while still respecting mass-balance. 2 geochemical cycles can be "locked": the sulfur cycle, and the oxygen cycle.
The model "lock" a cycle be imposing that the sources balance the sink at each timestep. More specifically, the sinks (oceanic
processes) are computed freely, and the sources (continental weathering) are force to match the sinks.
* For the sulfur cycle, the sulfuric silicate weathering is still compute freely, and the sulfuric carbonate weathering is adjusted
so that the sum of the two matches the sulfate reduction (the release of H2SO4 is set to 0).
* For the oxygen cycle, the kerogen weathering is adjusted so that when added to sulfide weathering, it matches the organic carbon
burial (whether or not the sulfur cycle is locked).

The mass balance is still respected, which means those modified fluxes affect the other geochemical species (carbon, alkalinity...)

To lock one or several cycle, set the value `.true.` of the corresponding parameters in source/coupler.inc (*before compilation*),
or use the options `--lock OS` (O for oxygen, S for sulfur, a single one works as well) in `build_GEOCLIM`.




## Further information


### Calibration procedure

**This section needs to be improved!**

To be fully consistent, the model should be recalibrated for each new set of boundary conditions.
The current calibration is done with ERA5 reanalysis fields for temperature and runoff, SRTM slope, and
Hartmann et al. 2013 lithology mask, all at a resolution of 0.5°. A second calibration is available for the GFDL boundary
conditions.
The climate fields of a General Circulation Model will inevitably differ from ERA5 fields, and differ from one GCM to
another, yielding many unique geochemical steady-states. The spatial resolution may also affect the steady-state.

Here are suggested steps to properly recalibrate the model:
* Run a Pre-Industrial simulation (1850 boundary conditions, 1xCO2) with the GCM you intend to use, preferably with
the same set of components and resolution. Retrieve the equilibrium annual climatology.
* Configure GEOCLIM at the given geographic resolution and with 1 CO2 level. It will set the model in fixed CO2 mode.
If you are using the 'GCM' input mode, you will need to remove (or comment) the lines stating the netCDF inputs that are
not at 1xCO2 in config/GCM\_input\_conditions. If you are using the 'ascii' input mode, you will need to remove the not-1xCO2
inputs *in the ascii files* (or create new ascii files with only 1xCO2 inputs).
* For the lithology mask, it is preferable to keep it the same way than you intend to use it for the paleo runs
(uniform or spatially-resolved, with same number of classes), while being consistent with present-day lithology.
* Do a first short run with the pre-industrial forcings to retrieve the silicate weathering flux. If you used DynSoil
in "dynamic" mode, run the model from "startup:eq" with acceleration tuning (see previous paragraph) during 1000-10000 yr.
If you use another set of components, only 1 model time step is needed.
* Get the "total silicate weathering" flux and the "total sulfuric silicate weathering flux" from the outputs (traditionnaly,
thos variables are named 'sil\_wth\_C\_flux' and 'sil\_sulfwth\_Ca\_flux'). The sum of the 2 must be used as CO2 degassing flux.
The degassing flux value should around 2-6 Tmol/yr.
As today's degassing flux is not well constrained, it is better to tuned it and keep the weathering parameters unchanged.
Note that the degassing flux (specified in "config/cond\_p20.dat") is split in 2: Volcanic (continental) and MOR (oceanic).
Doing so, the equilibrium CO2 with Pre-Industrial boundary conditions will be 1 PAL.
* Although it is not strictly necessary, you may want to adjust the parameters of Phosphorus and carbonate weathering to get the
desired flux. Phosphorus weathering will impact the oxygen levels. Carbonate weathering has no impact on equilibrium CO2, and
virtually no impact on O2 (though it may affect the biological pump). However, it directly impacts the oceanic DIC, Calcium and
alkalinity. Phosphorus weathering parameters (ie, P amount in source rocks: P\_rock, P2C\_ker and P2C\_carb) are defined in
source/constante.f90. Carbonate weathering should be modified directly in cont\_weath.f. In both cases, a simple cross-multiplication
is sufficient to get the right flux.
* Finally, re-do as many runs as needed to get 1 PAL of atmospheric O2 and 29 mol/m3 of mean oceanic sulfate (using acceleration
coefficients will help). There is no other way than to manually run the model, adjust the parameters if O2 and sulfate are too
low or too high, re-compile, re-run, and so on. I recommend tuning the value of the parameter 'OC\_in\_rocks' (in source/constante.f90)
that corresponds to "silicate sediments", because it is the most poorly constrained parameter. This parameter specifies the mass
fraction of petrogenic carbon in each rock type. A higher value will result in higher kerogen weathering, and thus less oxygen (and vice
versa). A standard value for silicate sediment is ~1%, though is highly depends on the type of sediment. For the sulfur cycle, the parameter
'Sulf\_rock' (amount of sulfide in source rocks, still in source/cont\_weath.f) is controlled by the S:C ratio, and determine the
sulfide weathering flux. With the acceleration parameters at 100, the model should be run for 2-5 Myr to have an idea of the equilibrium
O2 and SO4^2-
* To be more accurate, the oceanic alkalinity, DIC and the O2 gradient can be checked. If needed, they can be adjusted by modifying the
parameters controlling ocean chemistry, that are defined in source/constante.f90.


### Defining new physical units
When reading netCDF inputs (in "GCM" input mode), the code read the attribute "units" (a string) of the netCDF variables. If the "units"
string matches a defined ones, the corresponding conversion is performed to set the variable into the model's reference units.
It may be needed to add new unit string if the code do not recognize a given netCDF input file (for instance, in my experience, there are
as many runoff units as there are GCM). The string has to match exactly to be recognized (space and case sensitive).

Physical units are defined in the source code in "source/physical\_units.f90". To define a new one, go to the desired variable (eg,
'temperature\_units'), increment the variable 'naccepted' by 1, then go to the last line of 'accept\_unit' definition, and add a line
defining 'accept\_unit(n)%string' (your unit string) and 'accept\_unit(n)%conversion' (conversion factor and offset), 'n' being the number
of your newly defined units. The rule for converting variable into the reference unit is 'ref\_unit\_var = factor\*var + offset'.

Physical units of netCDF inputs are also interpreted by Python code "preproc/BC/BC\_generator.py".
The definition of units, and conversion parameters, are coded in "preproc/BC/units.py" (at the end of the file), in a similar way than in
"source/physical\_units.f90"


### Basic code modification

##### Model parameters
Most of the empirical parameters of GEOCLIM are defined in 'source/constante.f90'. Those parameters concern the oceanic and diagenesis
components (COMBINE module) and the continental weathering parameters, with the exception of DynSoil module parameters.

The chemical equilibrium constants, for oceanic chemistry (like carbonate speciation), are computed dynamically in 'source/eqcte.f90',
because they depends on temperature, pressure and salinity. Those relationships are constrained by thermodynamics, and less susceptible
to be modified.

All the parameters of DynSoil module are defined in 'source/dynsoil\_physical\_parameters.f90'

##### Output additional variables
There are a certain number of predefined output variables that the user can choose to output or not simply by editing the main
configuration file 'config/IO\_CONDITIONS' ("writevar(\*)" option, see previous section "Run GEOCLIM", "Output").
However, to output a variable that is not in that predefined list, here are the steps to follow.

Firstly, you need to know the name of the variable *in the Fortran source code*, or the way to compute it.

Secondly, depending on the type of the variable, it should be outputted in a different file:
* COMBINE variables. i.e., oceanic variables, that have a value for each COMBINE box (like salinity), or a single value (for instance,
atmospheric variable, or continental flux).
* Geographic variables. i.e., 2D geographic fields (for instance, runoff), or 3D if lithology-dependent (like weathering fluxes)
* DynSoil variables. Similar to geographic variables, but outputted only if DynSoil module is activated, and can be defined on lithology
and/or DynSoil vertical dimension.

The next steps are:

* In 'source/'output\_size.inc, increment by 1 the parameter defining the number of variable *of the corresponding outputfile*
('nCOMBoutvar' for COMBINE, 'nGEOGoutvar' for geographic, 'nDYNSoutvar' for DynSoil). 
* In the main configuration file (config/IO\_CONDITIONS), in the section "OUTPUT CONDITIONS" and corresponding block ("COMBINE OUTPUTS",
"GEOGRAPHIC OUTPUTS", or "DYNSOIL OUTPUTS"), add a line in the namelist (respectively, "&CMB\_OUTPUT\_VAR", "&GEO\_OUTPUT\_VAR",
or "&DYN\_OUTPUT\_VAR") – e.g., copy and paste the last line – stating the name of the variable *in the netCDF output file* ("vname(\*)"),
its units ("units(\*)"), under which dimension it is defined ("defdim(\*,:)", must be consistent with the source code!), its fill-value,
and "long\_name" description. "\*" is the output variable number, i.e., the number incremented in the previous step.
* In the corresponding source file "source/...\_write\_output.f90" ("..." being, respectively, "geoclim", "geographic", or "dynsoil"),
at the end of the section "write output variables", add a block of lines (e.g., copy and paste the following block) that looks
like:
> i = 106  
>   if (COMB_outvar_info(i)%writevar) &  
>     call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fO2_odc_tot), &  
>                  stt=(/nt/), cnt=(/1/))
> !

Or, for multi-dimensional variables (like geographic or DynSoil variable):
> i = 12  
>   if (GEOG_outvar_info(i)%writevar) &  
>     call put_var(fid, varname=GEOG_outvar_info(i)%vname, &  
> var_real3D=real(reshape(wth_litho_wgh, shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/)  
> !  

> i = 9 ! Z  
> if (DYNS_outvar_info(i)%writevar) &  
>   call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real4D=real(reshape(z, shape=(/nlon,nlat,nDSlev,nlitho/), &  
>                                            order=(/4,3,1,2/))), stt=(/1,1,1,1,nt/), cnt=(/nlon,nlat,nlitho,nDSlev,1/))  
> !  

"i = " is the output variable number, the one defined in the output namelist in IO\_CONDITIONS. "\*\_outvar\_info" is the namelist
variable defined in IO\_CONDITIONS.

Geographic variables have the particularity that they need a "reshape", as the two horizontal dimensions are unfolded in 1 dimension,
the lithology dimension is before the horizontal dimensions (and needs to be put at the end in the netCDF outputs), and for DynSoil
vertical variables, the vertical dimension ("z") is the first (before lithology and horizontal dimensions), and must be the last
in the netCDF outputs. This is the meaning of "order=..."


### Advanced customization

##### Add a new geochemical species
The main *prognostic* variables (for which advection and chemical reactions are computed) are stored in 3 Fortran variables: "var\_diss"
(for dissolved species, e.g., DIC), "var\_part" (for particulate species, e.g., PIC), and "var\_isot" (for isotopic variables, e.g., d13C).  
The variables are expressed in *concentration* (mol/m3) in the oceanic boxes and in *amount* (mol) in the atmospheric box, except the
isotopic variables, that are ratios (or "deltas") in all the boxes.

Here are the instruction to add a new main ocean-atmosphere variable:

* Determine if it is an dissolved, particulate, or isotopic variable.
* Update the corresponding number of variables (parameters 'nvar\_diss', 'nvar\_part', or 'nvar\_isot' in "source/combine\_foam.inc").
* Add a section in 'source/creades.f' (note that it is organized in blocks of 'diss', 'part', and 'isot' variables). It is the subroutine
  defining the source and sink fluxes due to chemical reaction, continental input, sinking (for particulate variables only) and
  sedimentation on seafloor. This net local source-sink rate is the code variable 'R'.
* Update whatever routine needed to compute the geochemical fluxes of that new variable.
* Update the section "Made-up semi-realistic concentrations" in "souece/geoclim\_mainprog.f (lines 347-449), that generates the initial
  conditions. This way is easier than manually editing the restart files. Add a block, with initial values suitable for the newly added
  variable. Blocks look like:
> j0 = j0+nbasin  
> ! Calcium  
> y(j0+1 : j0+nbasin-1) = 9d0  
> y(j0+nbasin) = 0 ! atmosphere

  Warning: the order of restart values must be: 1. all the 'diss' variables, 2. all the 'part' variables, 3. all the 'isot' variables.
  Make sure to write the initial values of the new variable in the right place!
* At the time of the first GEOCLIM run with this added variable, use the "internally-generated initial condition" option
  (`COMB_init_mode='coldstart'` in "config/IO\_CONDITIONS"),

Do note that this newly added "main" variable will not be outputted. To do so, follow the procedure described in
*output additional variables*, *COMBINE variable* section.
To compute (and output) the oceanic average of the variable, update the lines 33-106 of 'source/geoclim\_write\_output.f90'.
Note that the oceanic average of isotopic ratios is different than concentration variables.




## Visualization
A couple of external scripts are designed for the visualization of GEOCLIM output, in 'visualization/python/' (Python scripts) and
'visualization/jnl/' (Ferret scripts)

'visualization/python/plot\_final\_state.py' draw the main oceanic profiles and geochemical fluxes at the end of a run and save them
in two pdf files 'final\_fluxes--\*.pdf' and 'final\_ocean\_chemistry--\*.pdf'.
Usage: `python plot_final_state.py geoclim_output_file_path`.



## Technical notes

### Install ifort compiler on Mac OS
Note: the following instructions worked in June 2018, on Mac OS High Sierra 10.13.5

* download the student version
    * https://software.intel.com/en-us/qualify-for-free-software/student
* click `macOS (Fortran)`
* follow sign up
* download and install
    * `source /opt/intel/bin/compilervars.sh intel64`
* documentation
    * file:///opt/intel/documentation_2018/en/ps2018/getstart_comp_mf.htm

### How to install netCDF-Fortran library

#### With ifort compiler, on Mac OS
Check the following before completing the steps outlined below:
* check if gfortran is installed by typing `gfortran` onto the command line
    * there may be conflict issues between ifort and gfortran (these instructions are meant for ifort compiler)
    * if the command is recognized, then gfortran is installed, and may need to be removed
* XCode is (probably) required - make sure it is installed, up to date, and the license has been accepted
    * to accept the license: `sudo xcodebuild -license`
* if something went wrong with the installation, make sure that the source directories (zlib, hdf5, netcdf) are "fresh"
    * don't use a source directory that has been used before
    * instead, delete the old source directory, redownload and unzip, and attempt installation again
* restarting the computer after step 1 and step 4 may be helpful

1. netCDF-C
    * download zlib:
        * http://www.zlib.net
    * download hdf5
        * https://www.hdfgroup.org/downloads/hdf5/source-code/
    * download netCDF
        * https://github.com/Unidata/netcdf-c/releases/v4.6.1
    * unzip and `cd` into the zlib directory
        * `export ZDIR="/usr/local"`
        * `./configure --prefix=${ZDIR}`
        * `make check`
        * `make install`
    * unzip and `cd` into the hdf5 directory
        * `export H5DIR="/usr/local"`
        * `./configure --with-zlib=${ZDIR} --prefix=${H5DIR} --enable-hl`
        * `make check`
        * `make install`
    * unzip and `cd` into the netCDF directory
        * `export CPPFLAGS="-I${H5DIR}/include"`
        * `export LDFLAGS="-L${H5DIR}/lib"`
        * `./configure`
        * `make check`
        * `make install`
2. netCDF-F
    * download netCDF-F:
        * https://www.unidata.ucar.edu/downloads/netcdf/index.jsp
    * unzip and `cd` into the netCDF directory
        * `export NCDIR="/usr/local"`
        * `export NFDIR="/usr/local"`
        * `export CPPFLAGS="-I${NCDIR}/include"`
        * `export LDFLAGS="-L${NCDIR}/lib"`
        * `./configure --prefix=${NFDIR}`
        * `make check`
        * `make install`

#### With Linux OS
The simplest way is with apt: `apt install libnetcdff`
This will ensure the compatibility with the installed Fortran compiler.
Issues may however rise if several Fortran compiler are installed.

### Non-required software instructions for installation (Mac OS)
* pyFerret
    * download for python 3.6
        * https://github.com/NOAA-PMEL/PyFerret/releases
    * follow instructions here:
        * https://github.com/NOAA-PMEL/PyFerret#installation-from-prebuilt-targz-file
    * every time a new terminal is opened, you have to run `. ferret_paths.sh` (note the space)
        * recommend adding to `.bash_profile`:
            * `alias enable_pyferret=‘. /path/to/dir/pyferret-7.4-MacOSX-Python-3.6/ferret_paths.sh’`




## Frequent issues

**Note: this section has not been updated. Source file names and line numbers in error examples may not be up-to-date**


### Errors during compilation

#### incompatibility with previously compiled files
It happens when some files were previously compiled with another compiler, or other options, and are not compatible
with your new compilation command. Type `make clc` in the "source/" directory (or use option `--reset` in `build_GEOCLIM`)
and see if the error still persists.

#### netCDF library
This is by far the most frequent source of error, and they can be hard to detect and solve. They are basically 2 possibilities:

###### library not found
If the netCDF library is not found, the compilation fails when it reaches the file "io\_module.f90". The error should look like:

> io_module.f90:193:5:  
>   
>   use netcdf  
>      1  
> Fatal Error: Can't open module file ‘netcdf.mod’ for reading at (1): No such file or directory  
> compilation terminated.  
> make: \*\*\* [Makefile:229: io_module.o] Error 1  

To be sure, you can try and compile only the file "netcdf\_io\_module.f90" (that does not use any other source file),
for instance: `make netcdf_io_module.o [your potential Make options]`
The error message should look like:

> netcdf_io_module.f90:8:5:
>
>   use netcdf  
>      1  
> Fatal Error: Can't open module file ‘netcdf.mod’ for reading at (1): No such file or directory

This error occurs when no netCDF library exists in the specified path. Check that your compilation command
contains `-I.../include -l.../lib` (for instance, try `make echo`). If it does contain it, the paths `.../include'
and '.../lib' probably do not have netCDF library. Try to find where the library is installed.
This information can normally be obtained with `nc-config --prefix` (note that it is what `build_GEOCLIM` uses).

###### library not recognized
This error may be harder to detect. Sometimes, the compiler gives a specific indication, for instance, with gfortran:

> io_module.f90:193:5:  
>  
>   193 |  use netcdf  
>       |     1  
> Fatal Error: Cannot read module file ‘netcdf.mod’ opened at (1), because it was created by a different version of GNU Fortran  
> compilation terminated.  
> Makefile:229: recipe for target 'io_module.o' failed  
> make: *** [io_module.o] Error 1

Sometimes, the compilation of the modules and objects is successful, but the error comes while creating the executable, and the
compiler returns plenty of error messages, which does not make the task any easier, most of them looking like:

> /usr/bin/ld: /tmp/cc7LF9wU.o: in function \`\_\_netcdf_io_module_MOD_put_att_int':  
> netcdf_io_module.f90:(.text+0x8148): undefined reference to \`\_\_netcdf_MOD_nf90_put_att_one_fourbyteint'

This indicates that the compiler failed to used the netCDF library. 
It can happen for several reasons:
* Some compilation options are missing. Compilers generally need a specific option to use netCDF, `-lnetcdff` (gfortran),
`-lnetcdf` (ifort). The Makefile put those 2 options among the compilation flags (unless you override the flags with `FFLAGS=...`).
Make sure that you have those 2 options in your compilation command.
Sometimes (notably with ifort on Mac OS) those options need to be *at the very end* of the compilation command, which is what the
Makefile normally does.
* Incompatible library: another possibility is that the netCDF library is not compatible with the compiler. It may be a version
issue, or the fact that it is installed for the wrong compiler (try `nc-config --fc` to check the compiler the library is
configured for). In that case, there is no better solution than reinstalling the netCDF library (or using another Fortran compiler).

#### Fortran fixed-format interpretation
Normally, this type of error should not happen with the Makefile. With plain compilation command, however, one must be careful to
specify free format interpretation (`-ffree-form` with gfortran, `-free -132` with ifort...), as by default. Some compilers consider
files with extension '.f' (or all files) as Fortran 77 format.
This ".f" file interpretation is also an implicit rule in `make`, but this rule is overridden in the current Makefile.


### Error during execution

#### With netCDF input format
ie, input mode = 'GCM'.

With that format, the code is able to do many compliance checks and detect most error sources (and notify the user).
The code does not check invalid value for area (or land area) and temperature.
Be careful for instance if you define land fraction as a difference of variables, not to generate negative values.

The netCDF-Fortran library may not be able to read file in the most recent netCDF versions, like netCDF4.
The code is meant to read netCDF3 "classic" format input files.
This should not be an issue for GCM outputs, but be careful if you export data (like slope and lithology) in netCDF
format. In python with netCDF4 package, specify "format=NETCDF3\_CLASSIC".

#### With ascii input format
ie, input mode = 'ascii'

The code does not perform any checks besides negative runoff and slope. There can be several sources of errors.
The run will not necessarily crash, it sometimes continues with NaN or Infinity values. If that happens, recompile the code
with debug options (use `MODE=debug` with the Makefile, or `--mode debug --reset` with `build_GEOCLIM`) and re-run it.
It will tell you where the first error happened.

The "standard" GEOCLIM ascii format for geographic fields is:
* values unraveled with increasing latitude and longitude, longitude being the most rapidly varying axis.
* for total area and continental area, one pixel value per line, or values separated by commas or space (works as well)
* for climatic fields (temperature and runoff), **the first value must be the current CO2 level**, then, all the pixels
  values of the field (similarly unraveled), then the next CO2 level, and so on.
  Usually, each line is for one CO2 level, and the values are separated by comma or blank, but it works just as well with
  line breaks, or even one value per line, as long as the order is respected. 
  CO2 levels must be in **decreasing order**.

Slope and lithology mask must be in netCDF format. Note however that you can specify a uniform lithology directly in
"config/IO\_CONDITIONS", and that slope file is only needed if DynSoil module is activated.

###### Missing values
It is possible that on some continental points (ie, points with area > 0), climatic fields have missing value (runoff notably).
Note that the code will always check if there are points with negative runoff.

Re-compile the code with debugging option to determine where the first error occured (with `build_GEOCLIM`, add options
`--mode debug --reset`, with Makefile, do `make clc`, then `make MODE=debug [your personal options]`)
If the missing value is far enough from the valid range (like -9d33, or 1d36), you should get an error message like:

> Program received signal SIGFPE: Floating-point exception - erroneous arithmetic operation.

In one of those files:

> 0x55a5bfc51603 in eqcte\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5\_sulf/source/eqcte.f:16

> 0x55f31b341c98 in phfunc\_  
>	at /home/piermafrost/GitHub/GEOCLIM5\_sulf/source/phfunc.f:34

> 0x55fc6d79cb98 in bio\_frac\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5\_sulf/source/bio_frac.f:8

> 0x5578b564025d in carbo\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5\_sulf/source/carbo.f:9

> 0x55e19454c99d in newton\_  
>	at /home/piermafrost/GitHub/GEOCLIM5\_sulf/source/newton.f:23

> 0x55a479b87442 in ocean\_atm\_flu\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5\_sulf/source/ocean_atm_flu.f:16

If the missing value is close from the valid range, there may be no error, simply wrong continental fluxes.

###### Shifted grid
If the formatting of climatic or slope fields is different than the one of area (eg, flipped latitude axis, or longitude
starting at -180° instead of 0°), the code will likely read missing values, and the same errors than previous paragraph
will happen. However, it may read regular temperature values out of continents, or null runoff, depending on how the input
ascii file handle non-continental points. In that case, you will simply have wrong continental weathering fluxes, which may
be difficult to identify.

###### Units
Ascii files carry no information on variable unit, so the code cannot check it.  
The unit assumed by the code are:
* area: m2
* temperature: °C
* runoff: cm/yr

If the temperature is in Kelvin, you should receive (with debugging compilation options) an error message like:

> Program received signal SIGFPE: Floating-point exception - erroneous arithmetic operation.
>
> 0x55a479b87442 in ocean\_atm\_flu\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5\_sulf/source/ocean_atm_flu.f:16

Wrong runoff or area units will generally not trigger a crash, but will generate aberrant continental fluxes:
too high by a factor ~10 (runoff in mm/yr), or too low by a factor that is often 1/10, 1/100, 1d-6, 1d-12.

The output variable "discharge" (water discharge) is a good indicator of wrong runoff or area units, as its order
of magnitude is normally ~4d13 m3/yr.
You can also simply check the climatic and area variables in the geographic output file.

#### Error with COMBINE input data, or initial state
Combine input data (size of oceanic basins, seawater temperature...) is expected to be different for each paleo configuration.
It is possible that a new input dataset may make the model crash, because of error in its generation, or because it is not
compatible with the initial condition.

To see if the errors come from Combine input data, try and re-run the model with the reference dataset (files in
"INPUT/COMBINE/ref/"), keeping your COMBINE initial condition.

If the COMBINE initial condition is solely responsible for the crash, reducing the time step just for the time to the
ocean mixing to dissipate aberrant concentrations (100-1000 years) may solve the problem. The model could then be
run normally from the new restart.




## Notes and acknowledgements
Reference for ERA5 climate dataset:
Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I.,
Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2019): ERA5 monthly averaged data on single" levels from 1979 to present.
Copernicus Climate Change Service (C3S). Climate Data Store (CDS). (accessed on 19 Feb 2020)
https://doi.org/10.24381/cds.f17050d7
distributed under Copernicus Products license: https://cds.climate.copernicus.eu/api/v2/terms/static/license-to-use-copernicus-products.pdf




## Contact
pierre.maffre@normalesup.org

godderis was here

