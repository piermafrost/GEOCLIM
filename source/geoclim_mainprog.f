!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   + COupled Model of BIogeochemical cycles aNd climatE +
!   + COMBINE                                            +
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   Calibrated with ERA5 reanalysis (continental temperature and runoff),
!   STRM (slope) and Hartmann & Moosdorf 2012 (lithology map),
!   resolution: 0.5deg x 0.5deg


program geoclim



    use omp_lib
    use constante,                     only: PI_n_CO2_atm
    use GCM_io_module,                 only: load_climatology
    use geoclim_create_output_mod,     only: geoclim_create_output
    use geographic_create_output_mod,  only: geographic_create_output
    use dynsoil_read_input_mod,        only: dynsoil_read_input
    use dynsoil_create_output_mod,     only: dynsoil_create_output
    use dynsoil_initialization_mod,    only: dynsoil_initialization
    use io_module,                     only: netcdf_output_var, check_invalid
    use output_netcdf2ascii_mod,       only: output_netcdf2ascii
    use read_oceanic_temperature_mod,  only: read_oceanic_temperature
    use utils,                         only: set_error_handling_option
    use climatic_parameters,           only: get_clim_param
    implicit none
    include 'combine_foam.inc'
    integer:: ierr
    ! Declare structures for netCDF output info (cannot be declared in combine_foam.inc)
    type(netcdf_output_var), dimension(nCOMBoutvar) :: COMB_outvar_info
    type(netcdf_output_var), dimension(nGEOGoutvar) :: GEOG_outvar_info
    type(netcdf_output_var), dimension(nDYNSoutvar) :: DYNS_outvar_info





!   ####################################################################################
    open(unit=0, status='old', action='read', file=geoclim_path//'config/IO_CONDITIONS')
!   ####################################################################################


    ! Set the variable ERROR_HANDLING_OPTION telling the code what to do when
    ! some errors are encountered
    call set_error_handling_option()


    ! read run name, output directory and open cond_p20:
    call READ_MAIN_IO_CONDITIONS(0)


    ! Display run name
    print *
    print *
    n = len(trim(run_name))
    do k = 1,n+13
        write(*, fmt='(A1)', advance='no') '#'
    end do
    print *
    write(*, fmt='(A)') '##  RUN: '//trim(run_name)//'  ##'
    do k = 1,n+13
        write(*, fmt='(A1)', advance='no') '#'
    end do
    print *


    call read_conditions
    tbegin=tsta
    tend=tfi
    stept=ts
!    check=ijump_print
    time=tbegin
!    compteur=ijump_print+1
    htry=ts
    eps=1.d-8
    icount_geogprint = ijump_geogprint
    icount_veget = ijump_veget
    if (.not. coupling_veget) then
        !no modification
        veget_factor = 1
        veget_eros_factor = 1
    end if
    icount_climparam = ijump_climparam


!   ############ FILE OPENNING ############ !
    call OPEN_ASCII_FILES(0)
!   ####################################### !


!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   Read climatic inputs (global, land, ocean) and other land data:

    if (cont_input_mode == 'ascii') then

        if (nclimparam > 0) then
            print *, 'ERROR: use of climatic parameters (other than CO2) is not supported in "ascii" input mode'
            stop
        end if

        ref_x_axis = 0d0
        ref_y_axis = 0d0
        ! => tells the code to create them later

        do i=1,npixel  !S ---> N
            read(7,*)areaclimber(i)  !continental area
        end do
        close(unit=7)

        do i=1,npixel  !S ---> N
            read(301,*)areaEarth(i)  !continental and oceanic area
        end do
        close(unit=301)

        ! Note: both areas are expected to be in m2

        do j=1,nclimber
            read(30,*)co2climber(nclimber+1-j), &
                    (Tairclimber(i,nclimber+1-j,1,1,1,1,1),i=1,npixel)
            read(31,*)dummy,(Runclimber(i,nclimber+1-j,1,1,1,1,1),i=1,npixel)
        end do
        close(unit=30)
        close(unit=31)

        GMSTclimber = -1d99 ! no definition of GMST in inputs

        ! Units conversion: CO2: ppmv -> PAL
        co2climber = co2climber/280.

        ! temperature of oceanic basins
        call read_oceanic_temperature(32, Toceclimber, co2_axis=co2climber)


    else ! cont_input_mode == 'GCM'

        call load_climatology(333)
        ! Get global variables:
        !    'co2climber'
        !    'ref_x_axis'
        !    'ref_y_axis'
        !    'areaEarth'
        !    'areaclimber'
        !    'Tairclimber'
        !    'Runclimber'
        !    'Toceclimber'
        !    'GMSTclimber'

    end if

    ! Units conversion: area: m2 -> 1e6km2:
    areaclimber = areaclimber/1e12
    areaEarth = areaEarth/1e12

    ! get total land area
    areatot = sum(areaclimber)

    if (.not. uniform_lithology) call load_lithology(0)
    ! get global variable 'litho_frac' and check consistency with continental grid

!   Check for negative runoff
    print *
    call check_invalid('runoff', areaclimber, ERROR_HANDLING_OPTION(3), var7D=Runclimber, axis=2)

    ! Store GMST in last box (atmosphere) of Toceclimber:
    Toceclimber(nbasin,:,:,:,:,:,:) = GMSTclimber


!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        


    ! Read other oceanic inputs

    ! Obsolete: hypso file
    !do j=1,13
    !    read(4,*)deephyp(j),shyp(j)
    !end do
    !close(unit=4)

    areal=0.

    call basin_geometry()

    close(unit=333)

!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        




!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        

!   <><><><><><><><><><><><>
!   read initial conditions
!   <><><><><><><><><><><><>


    ! Rewind main IO file, because earlier namelists need to be read
    rewind(unit=0)


!   CLIMATIC PARAMETERS: (potentially needed to generate DynSoil initial condition) 
!   --------------------

    ! Lower bounds of climatic parameters (needed for climatic parameters with periodic ranges)
    cp_lowest_value = (/clim_param_1(1), clim_param_2(1), clim_param_3(1), clim_param_4(1), clim_param_5(1)/)

    call get_clim_param(time, tend, ijump_climparam, icount_climparam, prev_cpvec, next_cpvec, cpvec, cp_lowest_value, &
                        initialization=.true.)


!   COMBINE:
!   --------

    do j=1,nvb
        read(2,*)y(j)
    enddo
    close(unit=2)

    ! convert extensive variables: concentration => chemical amount
    ! (note: isotopic variables have vol(j)=1)
    y = y*vol

    if (nclimber==1) then
        ! fixed atm CO2 run => set atmospheric CO2 as defined by the unique CO2 level
        ! = co2_level (PAL) * Pre-Ind moles of CO2 in atm (0.0508d18)
        ! (atm PCO2 is 12th variable of last basin (nbasin))
        y(7*nbasin) = co2climber(1)*PI_n_CO2_atm
    endif

    call varset(y)


!   ECOLOGY MODULE:
!   ---------------

    if (coupling_ecogeo) then
        print *, 'ECOGEO module not available'
        stop -1
    end if


!    DynSoil
!    -------

    if (coupling_dynsoil) then

        ! read init file info or create init cond
        call dynsoil_read_input( 0,                                                               &
                                 ref_x_axis,ref_y_axis,areaclimber, xlevs,                        &
                                 reg_thick,reg_x_surf,reg_tau_surf,reg_z_prof,reg_tau_prof,slope, &
                                 DS_missingpoints, DYNS_restart_dim, DYNS_restart_var             )

        ! initialisation of offline dynsoil variables (ktop and P_vol):
        call dynsoil_initialization( xlevs, reg_thick, reg_x_surf, reg_z_prof, DS_missingpoints,  &
                                     reg_P_vol, reg_ktop                                          )

        DS_timestep = ijump_DS_integration * ijump_cont_weath * ts
        icount_DS_int = ijump_DS_integration
        icount_DS_pri = ijump_DS_print

    end if

!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        




!   Create geographic coordinate (longitude and latitude) if they were not given in input data
    if (all(ref_x_axis==0d0) .and. all(ref_y_axis==0d0)) then
        ref_x_axis = (/ ( ((dble(i)-0.5)/nlon) * 360 - 180 + 360/nlon/2   , i=1,nlon ) /)
        ref_y_axis = (/ ( ((dble(j)-0.5)/nlat) * 180 - 90  + 180/nlat/2   , j=1,nlat ) /)
    end if

!   Get the list of continental pixels (areaclimber > 0)
    call get_cont_pixel()




!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        

!   <><><><><><><><><><>
!   Create output files
!   <><><><><><><><><><>

!   read combine (oceanic) output conditions and create output file(s):
    call geoclim_create_output(0, output_directory, run_name, box_vol, oce_surf, surf_sedi, &
                               COMB_ofile_name, COMB_time_dimname, COMB_outvar_info     )
    ! Special case of GMST (optional input, output variable #93)
    where(GMSTclimber==-1d99)
        GMSTclimber                     = COMB_outvar_info(93)%fillval
        Toceclimber(nbasin,:,:,:,:,:,:) = COMB_outvar_info(93)%fillval
    end where

!   read geographic output conditions and create output file(s):
    call geographic_create_output(0, output_directory, run_name, ref_x_axis, ref_y_axis, areaclimber, litho_frac, slope, &
                                  GEOG_ofile_name, GEOG_time_dimname, GEOG_outvar_info                                   )

!   read dynsoil output conditions and create output file(s):
    if (coupling_dynsoil) then
        call dynsoil_create_output(0, output_directory, run_name, xlevs, ref_x_axis, ref_y_axis, areaclimber, litho_frac, slope, &
                                   DYNS_ofile_name, DYNS_time_dimname, DYNS_outvar_info                                          )
    end if

!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        


!   Restart file not yet created:
    fog = 0


!   Close IO condition file:
!   ------------------------
!   ########
    close(0)
!   ########


!   <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   !
!   indicate that the current run is done with its configuration files,
!   they can therefore be safely modified by other GEOCLIM program running.
!   <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   !
    open(unit=0, file=geoclim_path//'.config-queue', status='old', action='readwrite', iostat=ierr)
    if (ierr/=0)  open(unit=0, file=geoclim_path//'.config-queue', status='replace', action='readwrite')
    ierr=0
    do while (ierr==0) ! read until end of file
        read(unit=0, fmt=*, iostat=ierr)
    end do
    if (ierr>0) read(unit=0, fmt=*) ! if error other than end-of-file raised, re-execute reading action to display it
    backspace(unit=0)
    write(unit=0, fmt='(A)') 'free (run: '//trim(run_name)//')'
    close(unit=0)
!   <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   !




!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Put missing-value on non-continental points (area==0) of time-varying geographic output variables
!   (fillvalues are read during creation of output files, time-constant variables already written)

    do j = 1,npixel
        if (areaclimber(j)==0) then
            Tclim(j)           = GEOG_outvar_info(8)%fillval
            Runclim(j)         = GEOG_outvar_info(9)%fillval
            wth_allsil(j)      = GEOG_outvar_info(11)%fillval
            wth_litho_wgh(:,j) = GEOG_outvar_info(12)%fillval
            wth_litho(:,j)     = GEOG_outvar_info(13)%fillval
            fker(j)            = GEOG_outvar_info(14)%fillval
            POC_export_rate(j) = GEOG_outvar_info(15)%fillval
            fp(j)              = GEOG_outvar_info(16)%fillval
        end if
    end do

    ! DynSoil output variables
    if (coupling_dynsoil) then
        do j = 1,npixel
            ! missing-value for all variables on the carbonate lithology (#nlitho)
            reg_thick(nlitho,j)      = DYNS_outvar_info(6)%fillval
            reg_x_surf(nlitho,j)     = DYNS_outvar_info(7)%fillval
            reg_tau_surf(nlitho,j)   = DYNS_outvar_info(8)%fillval
            reg_z_prof(:,nlitho,j)   = DYNS_outvar_info(9)%fillval
            reg_tau_prof(:,nlitho,j) = DYNS_outvar_info(10)%fillval
            reg_x_mean(nlitho,j)     = DYNS_outvar_info(16)%fillval
            reg_mean_age(nlitho,j)   = DYNS_outvar_info(17)%fillval
            reg_Li_Friv(nlitho,j)    = DYNS_outvar_info(18)%fillval
            reg_Li_Fsp(nlitho,j)     = DYNS_outvar_info(19)%fillval
            reg_Li_driv(nlitho,j)    = DYNS_outvar_info(20)%fillval
            ! missing-value for initial fluxes on ALL pixels (can only be computed between two time steps)
            reg_prod(:,j)        = DYNS_outvar_info(11)%fillval
            reg_eros(:,j)        = DYNS_outvar_info(12)%fillval
            reg_P_diss(:,j)      = DYNS_outvar_info(13)%fillval
            reg_P_eros(:,j)      = DYNS_outvar_info(14)%fillval
            reg_x_surf_eros(:,j) = DYNS_outvar_info(15)%fillval
            ! missing-value on non-continental pixels for other output variables
            if (areaclimber(j)==0) then
                reg_thick(:,j)      = DYNS_outvar_info(6)%fillval
                reg_x_surf(:,j)     = DYNS_outvar_info(7)%fillval
                reg_tau_surf(:,j)   = DYNS_outvar_info(8)%fillval
                reg_z_prof(:,:,j)   = DYNS_outvar_info(9)%fillval
                reg_tau_prof(:,:,j) = DYNS_outvar_info(10)%fillval
                reg_x_mean(:,j)     = DYNS_outvar_info(16)%fillval
                reg_mean_age(:,j)   = DYNS_outvar_info(17)%fillval
                reg_Li_Friv(:,j)    = DYNS_outvar_info(18)%fillval
                reg_Li_Fsp(:,j)     = DYNS_outvar_info(19)%fillval
                reg_Li_driv(:,j)    = DYNS_outvar_info(20)%fillval
            end if
        end do
    end if

!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        

!   End of initialization:


    htot=0.
    icount=ijump_print+1
    icount_cont_weath = ijump_cont_weath

!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        




!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   !
!   dummy integration (do not actualized
!   prognostic variables) to get initiale fluxes

    print *
    print *
    print *, 'Compute initial values of fluxes (for outputs only)'

    ! SAVE INITIAL VAR:
    y_0 = y
    ! dynsoil prognostic variables:
    if (coupling_dynsoil) then
      reg_thick_0    = reg_thick
      reg_x_surf_0   = reg_x_surf
      reg_tau_surf_0 = reg_tau_surf
      reg_P_vol_0    = reg_P_vol
      reg_z_prof_0   = reg_z_prof
      reg_tau_prof_0 = reg_tau_prof
      reg_ktop_0     = reg_ktop
    end if

    call read_veget()
    ! no need to call "get_clim_param" => 'cpvec' variable already assigned during initialization
    call creades(time)
    call rk4(y,dydx,nvb,time,htry,yout)

    do j=1,nvb
      y(j)=yout(j)
    end do


    ! RESET VARIABLES WITH INITIAL VALUES:
    y = y_0
    if (coupling_dynsoil) then
      reg_thick    = reg_thick_0
      reg_x_surf   = reg_x_surf_0
      reg_tau_surf = reg_tau_surf_0
      reg_P_vol    = reg_P_vol_0
      reg_z_prof   = reg_z_prof_0
      reg_tau_prof = reg_tau_prof_0
      reg_ktop     = reg_ktop_0
      icount_DS_int = ijump_DS_integration
    end if
    call varset(y)
    icount_cont_weath = ijump_cont_weath
    icount_veget = ijump_veget
    ! 'icount_climparam' not modified during dummy integration

    call printf(time, icount, y, COMB_outvar_info, GEOG_outvar_info, DYNS_outvar_info)

    ! For not printing the initial state a second time but still counting it:
    if (coupling_dynsoil) icount_DS_pri = -1
    icount_geogprint = -1
!    compteur=-1
    icount=-1

!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   !




!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   %%         TIME LOOP         %%
!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    print *
    print *
    print *, 'Start run'
    call date_and_time(values=computer_time)
    write(*, fmt='(I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') &
        computer_time(1),'/',computer_time(2),'/',computer_time(3),' ', &
        computer_time(5),':',computer_time(6),':',computer_time(7),':',computer_time(8)


!   &&&&&&&&&&  start iteration  &&&&&&&&&&
    do while (time.le.tend)

      call read_veget()
      call get_clim_param(time, tend, ijump_climparam, icount_climparam, prev_cpvec, next_cpvec, cpvec, cp_lowest_value)
      call creades(time)
      call rk4(y,dydx,nvb,time,htry,yout)
      do j=1,nvb
        y(j)=yout(j)
      end do

      call varset(y)

      call printf(time, icount, y, COMB_outvar_info, GEOG_outvar_info, DYNS_outvar_info)
      time=time+htry

    end do
!   &&&&&&&&&&  end iteration  &&&&&&&&&&


    ! call 'printf' again only for creating restart files (in case asked to do at end of run)
    icount = -1
    icount_geogprint = -1
    icount_DS_pri = -1
    ! => do not print anything
    ageYprint = time - 2*ts
    ! => create restart, if not already created
    call printf(time, icount, y, COMB_outvar_info, GEOG_outvar_info, DYNS_outvar_info)


    print *
    print *, 'End of run'
    call date_and_time(values=computer_time)
    write(*, fmt='(I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') &
        computer_time(1),'/',computer_time(2),'/',computer_time(3),' ', &
        computer_time(5),':',computer_time(6),':',computer_time(7),':',computer_time(8)


    ! CONVERT NETCDF OUTPUTS TO ASCII OUTPUTS:
    if (convert2ascii) call output_netcdf2ascii(COMB_outvar_info)





end program

