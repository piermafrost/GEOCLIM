module dynsoil_physical_parameters

    implicit none


    !################################################################################################################!
    !# Vertical streching factor. "Artificial" parameter that modulate the inertia (i.e., response time) of         #!
    !# the regolith without affecting the steady-state weathering rate.                                             #!
    !# Keep at 1 for realistic model run. Modify to explore the sensitivity of transient state to regolith inertia, #!
    !# Decrease to accelerate the equilibrium-reaching in fast runs (e.g., down to 1d-3)                            #!
    !################################################################################################################!
    double precision, parameter:: scaling_factor = 1.


    ! UNIVERSAL PARAMETERS
    ! --------------------
    !
    ! Ideal gas constant (J/K/mol)
    double precision, parameter:: Rgas = 8.314


    ! REGOLITH (SOIL) PHYSICAL PARAMETERS
    ! -----------------------------------
    !
    ! Bedrock density (kg/m3)
    double precision, parameter:: rho_ucc = 2500


    include 'shape.inc'
    ! -> to get number of lithological classes: 'nlitho', number of vertical levels 'nDSlevs',
    !    also contains 'nlon' and 'nlat'

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !%  The following parameters are designed for 6 lithological classes: %!
    !%  #1: Metamorphic                                                   %!
    !%  #2: Felsic                                                        %!
    !%  #3: Intermediate                                                  %!
    !%  #4: Mafic                                                         %!
    !%  #5: Siliclastic sediments                                         %!
    !%  #6: Carbonates (ignored by DynSoil module)                        %!
    !%                                                                    %!
    !% WARNING: Carbonate **must be** the last lithological class         %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


    ! CONSTANTS FOR SOIL OPT. PRODUCTION RATE (m/y)
    ! ---------------------------------------------
    !
    ! Regolith production constante (-) (equivalent to water-to-rock ratio: optimal reg. prod. rate at T0 / runoff)
    double precision, dimension(nlitho), parameter:: krp   = (/ 1d-2, 1d-2, 1d-2, 1d-2, 1d-2, 0d0 /)
    ! Reference temperature (K)
    double precision, dimension(nlitho), parameter:: T0_rp = (/ 286.,   286.,   286.,   286.,   286.,   0.  /)
    ! Apparent activation energy at reference temperature T0 (J/mol)
    double precision, dimension(nlitho), parameter:: Ea_rp = (/ 42d3,   42d3,   42d3,   42d3,   42d3,   0d0 /)


    ! CONSTANTS FOR SOIL PRODUCTION FUNCTION
    ! --------------------------------------
    !
    ! Characteristic depth of decay (m)
    !    - for exponential SPF:
    double precision, parameter:: h0 = 2.73 * scaling_factor
    !!    - for inverse SPF:
    !double precision, parameter:: h0 = 1.21 * scaling_factor !(m) characteristic depth of decay
    !!
    !! SPF parameters used for "humped" SPF only:
    !double precision, parameter:: k1 = 0.95
    !double precision, parameter:: d1 = 10 * scaling_factor   !(m)
    !double precision, parameter:: d2 = 0.12 * scaling_factor !(m)
    !!    Maximum of non-normalized humped SPF:
    !double precision, parameter:: knorm = (k1*d1/d2)**(-d2/(d1-d2)) - k1*((k1*d1/d2)**(-d1/(d1-d2)))


    ! CONSTANTS FOR REGOLITH EROSION LAW (m/y)
    ! ----------------------------------------
    !
    ! Erodibility constant (m^(1-a)/yr^(1-a))
    double precision, dimension(nlitho), parameter:: ke = (/ 0.0030713, 0.0030713, 0.0030713, 0.0030713, 0.0030713, 0. /)
    ! Runoff exponent
    double precision, parameter:: a = 0.5
    ! Slope exponent
    double precision, parameter:: b = 1

    ! CONSTANTS FOR ROCK DISSOLUTION LAW
    ! ----------------------------------
    !
    ! Runoff sensitivity parameter (m^-1.y)
    double precision, dimension(nlitho), parameter:: kw    = (/ 1.,   1.,   1.,   1.,   1.,   0.  /)
    ! Age exponent for dissolution rate constant (-)
    double precision, dimension(nlitho), parameter:: sigma = (/ -0.4, -0.4, -0.4, -0.4, -0.4, 0.  /)
    ! Main dissolution constant (y^(-1-sigma))
    double precision, dimension(nlitho), parameter:: Kwest = (1/(scaling_factor**(sigma+1))) * &
                                                             (/ 5d-4, 5d-4, 5d-4, 5d-4, 5d-4, 0d0 /) ! standard value: 6e-5
    ! Reference temperature (K)
    double precision, dimension(nlitho), parameter:: T0    = (/ 286., 286., 286., 286., 286., 0.  /)
    ! Apparent activation energy at reference temperature T0 (J/mol)
    double precision, dimension(nlitho), parameter:: Ea    = (/ 42d3, 42d3, 42d3, 42d3, 42d3, 0d0 /)


    ! CONSTANTS FOR NUMERICAL INTEGRATION OF SOIL PROFILE
    ! ---------------------------------------------------
    !
    ! Threshold for new x-point creatiom (-)
    double precision, parameter:: epsl = 0.01
    ! Useless for dynsoil_steeady_state, only used for the dynamic version


    ! CONSTANTS FOR SECONDARY PHASES AND LITHIUM CYCLE
    ! ------------------------------------------------
    !
    ! Lithium abundance in bedrock (mol/kg)
    double precision, parameter:: Li_UCC = 0.011 !77ppm ! 176ppm !22ppm !
    ! Lithium delta isotopic in bedrock (-)
    double precision, parameter:: dLi_UCC = 0.0017
    ! Temperature coefficients for fractionation parameter: Delta = a/T^2 + b
    double precision, parameter:: aDland = 1.63e+3 , bDland = -2.04e-3


end module
