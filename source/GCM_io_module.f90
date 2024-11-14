module GCM_io_module
! contains subroutines handle to GCM outputs (ie, annual climatologies)
!   1. read the paths of configuration variables in the config file (config/GCM_io_condition)
!   2. read the continental variables (area, landfrac, temperature and runoff) from GCM outputs

    implicit none
    integer, parameter, private:: fmaxlen=500, vmaxlen=200
    double precision, parameter, private:: axis_relat_accuracy=1d-5
    ! maximum allowed relative difference between axis (lon, lat) from different files

    contains



    ! ###################################################################################### !
    ! Main subroutine of module: read all GCM input and put them in GEOCLIM common variables !
    ! ###################################################################################### !

    subroutine load_climatology(fID)

        use netcdf
        use io_module, only: UNDEFINED_VALUE_CHAR
        use read_oceanic_temperature_mod, only: read_oceanic_temperature

        include 'combine_foam.inc'

        integer, parameter:: paramspace_size = nclimber*len_p1*len_p2*len_p3*len_p4*len_p5

        integer, intent(in):: fID
        character(len=fmaxlen), dimension(2):: area_file
        character(len=fmaxlen), dimension(paramspace_size):: climo_file, gmst_climo_file
        character(len=vmaxlen), dimension(2):: area_x_varname, area_y_varname, area_varname
        character(len=vmaxlen):: x_varname, y_varname, temp_varname, runf_varname, var_units, vname
        character(len=vmaxlen):: gmst_x_varname, gmst_y_varname, glob_temp_varname
        double precision:: fillvalue

        integer, dimension(paramspace_size, 6):: paramspace_filling_order

        double precision, dimension(2+2*paramspace_size, nlon):: all_x
        double precision, dimension(2+2*paramspace_size, nlat):: all_y
        double precision, dimension(npixel, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5):: glob_temperature
        integer:: ierr, nax, i1, i2, i3, i4, i5, red_pspace_size


        print *
        print *
        print *, 'Read GCM conditions file and load land inputs from GCM climatology files + oceanic temperature from ascii file'
        print *, '--------------------------------------------------------------------------------------------------------------'

        ! get netCDF file info
        ! --------------------
        call read_GCM_condition(fID, &
        !                       <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                                co2climber, clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5, &
                                red_pspace_size, paramspace_filling_order, &
        !                       <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                                area_file(1), area_x_varname(1), area_y_varname(1), area_varname(1), &
                                area_file(2), area_x_varname(2), area_y_varname(2), area_varname(2), &
                                climo_file, x_varname, y_varname, temp_varname, runf_varname, &
                                gmst_climo_file, gmst_x_varname, gmst_y_varname, glob_temp_varname)

        ! Units conversion: CO2: ppmv -> PAL
        co2climber = co2climber/280.


        ! get total area
        ! --------------
        call load_variable('area', area_varname(1), area_x_varname(1), area_y_varname(1), &
                           !                                        ,,,,,,,,,
                           single_input_file=area_file(1), varout1D=areaEarth, x=all_x(1,:), y=all_y(1,:), &
                           !                                        ^^^^^^^^^
                           fill_missing=.true.)

        ! get land fraction/area
        ! ----------------------
        call load_variable('landarea', area_varname(2), area_x_varname(2), area_y_varname(2), &
                           !                                        ,,,,,,,,,,,
                           single_input_file=area_file(2), varout1D=areaclimber, x=all_x(2,:), y=all_y(2,:), &
                           !                                        ^^^^^^^^^^^
                           totarea=areaEarth, fill_missing=.true.)

        ! get temperature
        ! ---------------
        call load_variable('temperature', temp_varname, x_varname, y_varname, &
                           !                                                           ,,,,,,,,,,,
                           multiple_input_file=climo_file(1:red_pspace_size), varout7D=Tairclimber, &
                           !                                                           ^^^^^^^^^^^
                           xvec=all_x(3:2+red_pspace_size,:), yvec=all_y(3:2+red_pspace_size,:), &
                           landarea=areaclimber, paramspace_filling_order=paramspace_filling_order, &
                           fill_missing=.false., fillval_handling_option=ERROR_HANDLING_OPTION(2))
        ! note: ERROR_HANDLING_OPTION(2) tells the code what to do if it found missing-value where area > 0

        ! get runoff
        ! ----------
        call load_variable('runoff', runf_varname, x_varname, y_varname, &
                           !                                                           ,,,,,,,,,,
                           multiple_input_file=climo_file(1:red_pspace_size), varout7D=Runclimber, &
                           !                                                           ^^^^^^^^^^
                           landarea=areaclimber, paramspace_filling_order=paramspace_filling_order, &
                           fill_missing=.false., fillval_handling_option=ERROR_HANDLING_OPTION(2))

        ! get optional global temperature
        ! -------------------------------
        if (glob_temp_varname /= UNDEFINED_VALUE_CHAR) then
            call load_variable('temperature', glob_temp_varname, gmst_x_varname, gmst_y_varname, &
                               multiple_input_file=gmst_climo_file(1:red_pspace_size), varout7D=glob_temperature, &
                               xvec=all_x(3+red_pspace_size:2+2*red_pspace_size, :), &
                               yvec=all_y(3+red_pspace_size:2+2*red_pspace_size, :), &
                               fill_missing=.true., paramspace_filling_order=paramspace_filling_order)
            ! number of axis to check
            nax = 2+2*red_pspace_size
        else
            ! number of axis to check
            nax = 2+red_pspace_size
        end if

        ! Check that all axis are the same
        ! --------------------------------
        call check_axis(all_x(1:nax,:), all_y(1:nax,:), ERROR_HANDLING_OPTION(1))
        ! note: ERROR_HANDLING_OPTION(1) tells the code what to do if it find axis mismatch


        ! Get oceanic temperature (from ascii file #32)
        ! ---------------------------------------------
        !                                 ,,,,,,,,,,,
        call read_oceanic_temperature(32, Toceclimber, co2_axis=co2climber, order=paramspace_filling_order, ndata=red_pspace_size)
        !                                 ^^^^^^^^^^^

        ! Get water exchange fluxes between basins (from ascii file #42)
        ! --------------------------------------------------------------
        !                                                                                    ,,,,,,,,,
        call load_oce_ascii_variable(42, red_pspace_size, paramspace_filling_order, varout8D=Fxch_clim, singleclim=variable_watxch)
        !                                                                                    ^^^^^^^^^
        variable_watxch = (.not. variable_watxch) ! single-clim => NOT variable water exchanges
        !
        if (.not. variable_watxch) then
            ! In case of unmutable oceanic circulation (water exchanges between basins) the common "Fxch_clim" is ignored,
            ! the common variable "F" is drectly used
            F = reshape(Fxch_clim(:,1,1,1,1,1,1), (/nbasin,nbasin/)) ! basin-to-basin ravelling
            F = 31.536d12*F ! conversion Sv => m3/yr
        end if


        ! Reference axis (GEOCLIM global variables)
        ! =========================================

        ref_x_axis = all_x(1,:)
        ref_y_axis = all_y(1,:)


        ! Fill last hyperslabs of climate arrays for climatic parameters with periodic ranges
        ! ===================================================================================

        if (p1_period /= 0d0) then
            Tairclimber(:,:,len_p1,:,:,:,:) = Tairclimber(:,:,1,:,:,:,:)
            Runclimber( :,:,len_p1,:,:,:,:) = Runclimber( :,:,1,:,:,:,:)
            Toceclimber(:,:,len_p1,:,:,:,:) = Toceclimber(:,:,1,:,:,:,:)
            if (glob_temp_varname/=UNDEFINED_VALUE_CHAR) glob_temperature(:,:,len_p1,:,:,:,:) = glob_temperature(:,:,1,:,:,:,:)
            if (variable_watxch)                                Fxch_clim(:,:,len_p1,:,:,:,:) =        Fxch_clim(:,:,1,:,:,:,:)
        end if
        if (p2_period /= 0d0) then
            Tairclimber(:,:,:,len_p2,:,:,:) = Tairclimber(:,:,:,1,:,:,:)
            Runclimber( :,:,:,len_p2,:,:,:) = Runclimber( :,:,:,1,:,:,:)
            Toceclimber(:,:,:,len_p2,:,:,:) = Toceclimber(:,:,:,1,:,:,:)
            if (glob_temp_varname/=UNDEFINED_VALUE_CHAR) glob_temperature(:,:,:,len_p2,:,:,:) = glob_temperature(:,:,:,1,:,:,:)
            if (variable_watxch)                                Fxch_clim(:,:,:,len_p2,:,:,:) =        Fxch_clim(:,:,:,1,:,:,:)
        end if
        if (p3_period /= 0d0) then
            Tairclimber(:,:,:,:,len_p3,:,:) = Tairclimber(:,:,:,:,1,:,:)
            Runclimber( :,:,:,:,len_p3,:,:) = Runclimber( :,:,:,:,1,:,:)
            Toceclimber(:,:,:,:,len_p3,:,:) = Toceclimber(:,:,:,:,1,:,:)
            if (glob_temp_varname/=UNDEFINED_VALUE_CHAR) glob_temperature(:,:,:,:,len_p3,:,:) = glob_temperature(:,:,:,:,1,:,:)
            if (variable_watxch)                                Fxch_clim(:,:,:,:,len_p3,:,:) =        Fxch_clim(:,:,:,:,1,:,:)
        end if
        if (p4_period /= 0d0) then
            Tairclimber(:,:,:,:,:,len_p4,:) = Tairclimber(:,:,:,:,:,1,:)
            Runclimber( :,:,:,:,:,len_p4,:) = Runclimber( :,:,:,:,:,1,:)
            Toceclimber(:,:,:,:,:,len_p4,:) = Toceclimber(:,:,:,:,:,1,:)
            if (glob_temp_varname/=UNDEFINED_VALUE_CHAR) glob_temperature(:,:,:,:,:,len_p4,:) = glob_temperature(:,:,:,:,:,1,:)
            if (variable_watxch)                                Fxch_clim(:,:,:,:,:,len_p4,:) =        Fxch_clim(:,:,:,:,:,1,:)
        end if
        if (p5_period /= 0d0) then
            Tairclimber(:,:,:,:,:,:,len_p5) = Tairclimber(:,:,:,:,:,:,1)
            Runclimber( :,:,:,:,:,:,len_p5) = Runclimber( :,:,:,:,:,:,1)
            Toceclimber(:,:,:,:,:,:,len_p5) = Toceclimber(:,:,:,:,:,:,1)
            if (glob_temp_varname/=UNDEFINED_VALUE_CHAR) glob_temperature(:,:,:,:,:,:,len_p5) = glob_temperature(:,:,:,:,:,:,1)
            if (variable_watxch)                                Fxch_clim(:,:,:,:,:,:,len_p5) =        Fxch_clim(:,:,:,:,:,:,1)
        end if


        ! conversion of water exchange fluxes
        if (variable_watxch)  Fxch_clim = 31.536d12*Fxch_clim ! Sv => m3/yr


        ! Compute GMST (GEOCLIM global variable)
        ! ======================================

        if (glob_temp_varname == UNDEFINED_VALUE_CHAR) then
            !,,,,,,,,,,
            GMSTclimber = -1d99 ! "internal" missing-value -> will be replaced by "output" missing-value
            !^^^^^^^^^^
        else
            totarea = sum(areaEarth)
            do i5 = 1,len_p5
                do i4 = 1,len_p4
                    do i3 = 1,len_p3
                        do i2 = 1,len_p2
                            do i1 = 1,len_p1
                                do k=1,nclimber
                                    !,,,,,,,,,,,,,,,,,,,,,,,,,,,,
                                    GMSTclimber(k,i1,i2,i3,i4,i5) = sum(glob_temperature(:,k,i1,i2,i3,i4,i5)*areaEarth) / totarea
                                    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end if


    end subroutine




    !======================================================================!




    subroutine read_GCM_condition(fID, &
                                  CO2_levels, clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5, &
                                  red_pspace_size, paramspace_filling_order, &
                                  area_file, area_x_dim_name, area_y_dim_name, area_var_name, &
                                  landarea_file, landarea_x_dim_name, landarea_y_dim_name, landarea_var_name, &
                                  land_list_file, land_x_dim_name, land_y_dim_name, land_temp_var_name, runoff_var_name, &
                                  glob_list_file, glob_x_dim_name, glob_y_dim_name, glob_temp_var_name)

        use io_module, only: check_namelist_def, UNDEFINED_VALUE_CHAR, UNDEFINED_VALUE_DBLE
        use climatic_parameters, only: retrieve_clim_param_space
        use utils, only: add_path

        include 'shape.inc'
        integer, parameter:: paramspace_size = nclimber*len_p1*len_p2*len_p3*len_p4*len_p5

        integer, intent(in):: fID
        double precision, intent(out):: CO2_levels(nclimber), clim_param_1(len_p1), clim_param_2(len_p2), clim_param_3(len_p3), &
                                        clim_param_4(len_p4), clim_param_5(len_p5)
        integer, intent(out):: red_pspace_size
        integer, dimension(paramspace_size, 6), intent(out):: paramspace_filling_order
        character(len=fmaxlen), intent(out):: area_file, landarea_file
        character(len=vmaxlen), intent(out):: area_x_dim_name, area_y_dim_name, area_var_name, &
                                              landarea_x_dim_name, landarea_y_dim_name, landarea_var_name
        character(len=fmaxlen), dimension(paramspace_size), intent(out):: land_list_file, glob_list_file
        character(len=vmaxlen), intent(out):: land_x_dim_name, land_y_dim_name, land_temp_var_name, runoff_var_name
        character(len=vmaxlen), intent(out):: glob_x_dim_name, glob_y_dim_name, glob_temp_var_name
        double precision, dimension(paramspace_size,6):: all_climparam
        logical, dimension(nclimber,len_p1,len_p2,len_p3,len_p4,len_p5):: filled_space
        integer:: k, n_undef, n_def, alenp1, alenp2, alenp3, alenp4, alenp5

        ! Namelist declaration
        namelist /CLIM_PARAMS/    all_climparam
        namelist /AREA_INFO/      area_file, area_x_dim_name, area_y_dim_name, area_var_name
        namelist /LAND_AREA_INFO/ landarea_file, landarea_x_dim_name, landarea_y_dim_name, landarea_var_name
        namelist /LAND_CLIM_INFO/ land_x_dim_name, land_y_dim_name, land_temp_var_name, runoff_var_name, land_list_file
        namelist /GLOB_CLIM_INFO/ glob_x_dim_name, glob_y_dim_name, glob_temp_var_name, glob_list_file

        ! Default values of namelist variables
        all_climparam       = UNDEFINED_VALUE_DBLE
        area_file           = UNDEFINED_VALUE_CHAR
        area_x_dim_name     = UNDEFINED_VALUE_CHAR
        area_y_dim_name     = UNDEFINED_VALUE_CHAR
        area_var_name       = UNDEFINED_VALUE_CHAR
        landarea_file       = UNDEFINED_VALUE_CHAR
        landarea_x_dim_name = UNDEFINED_VALUE_CHAR
        landarea_y_dim_name = UNDEFINED_VALUE_CHAR
        landarea_var_name   = UNDEFINED_VALUE_CHAR
        land_x_dim_name     = UNDEFINED_VALUE_CHAR
        land_y_dim_name     = UNDEFINED_VALUE_CHAR
        land_temp_var_name  = UNDEFINED_VALUE_CHAR
        runoff_var_name     = UNDEFINED_VALUE_CHAR
        land_list_file      = UNDEFINED_VALUE_CHAR
        glob_x_dim_name     = UNDEFINED_VALUE_CHAR
        glob_y_dim_name     = UNDEFINED_VALUE_CHAR
        glob_temp_var_name  = UNDEFINED_VALUE_CHAR
        glob_list_file      = UNDEFINED_VALUE_CHAR


        ! CO2 and climatic parameters
        ! ---------------------------

        ! Read CO2 and climatic parameters values
        !<><><><><><><><><><><><><><><>!
        read(unit=fID, nml=CLIM_PARAMS)
        !<><><><><><><><><><><><><><><>!

        ! Get individual vectors of CO2 and climatic parameter, and the order in which the whole parameter space if filled 
        call retrieve_clim_param_space(all_climparam, &
                                       CO2_levels, clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5, &
                                       alenp1, alenp2, alenp3, alenp4, alenp5, &
                                       paramspace_filling_order, is_paramspace_filled=filled_space)
        ! "actual" parameter space size (removing 1 value for all periodic ranges)
        red_pspace_size = nclimber*alenp1*alenp2*alenp3*alenp4*alenp5

        ! Check that enough parameter combinations were read (<=> no undefined values in "all_climparam" array)
        n_undef = count(all_climparam(1:red_pspace_size, 1:1+nclimparam) == UNDEFINED_VALUE_DBLE)
        if (n_undef > 0) then
            print *
            print *, 'Error: not enough values of CO2 and climatic parameters were given in namelist "CLIM_PARAMS"'
            print *, 'in GCM configuration file'
            print *, n_undef, ' entries of CO2/climatic parameters are missing out of ', red_pspace_size
            stop 89
        end if

        ! Check that parameter space is filled
        n_def = count(filled_space(1:nclimber, 1:alenp1, 1:alenp2, 1:alenp3, 1:alenp4, 1:alenp5))
        if (n_def /= red_pspace_size) then
            print *
            print *, 'ERROR: underfilled parameter space (ie, CO2 and climatic parameters)'
            print *, 'The expected number of parameter combinations was loaded, which means there are duplicated combinations'
            print *, 'Filled parameter combinations: ', n_def, '/', red_pspace_size
            stop 89
        end if


        ! Area / land area / land climatology / global temperature
        ! --------------------------------------------------------

        !<><><><><><><><><><><><><><>!
        read(unit=fID, nml=AREA_INFO)
        !<><><><><><><><><><><><><><>!
        call check_namelist_def('Error: "area_file" not defined in "AREA_INFO" namelist of GCM input config file', &
                                char_var=area_file)
        call check_namelist_def('Error: "area_x_dim_name" not defined in "AREA_INFO" namelist of GCM input config file', &
                                char_var=area_x_dim_name)
        call check_namelist_def('Error: "area_y_dim_name" not defined in "AREA_INFO" namelist of GCM input config file', &
                                char_var=area_y_dim_name)
        call check_namelist_def('Error: "area_var_name" not defined in "AREA_INFO" namelist of GCM input config file', &
                                char_var=area_var_name)
        call add_path(area_file)

        !<><><><><><><><><><><><><><><><>!
        read(unit=fID, nml=LAND_AREA_INFO)
        !<><><><><><><><><><><><><><><><>!
        call check_namelist_def('Error: "landarea_file" not defined in "AREA_INFO" namelist of GCM input config file', &
                                char_var=landarea_file)
        call check_namelist_def('Error: "landarea_x_dim_name" not defined in "AREA_INFO" namelist of GCM input config file', &
                                char_var=landarea_x_dim_name)
        call check_namelist_def('Error: "landarea_y_dim_name" not defined in "AREA_INFO" namelist of GCM input config file', &
                                char_var=landarea_y_dim_name)
        call check_namelist_def('Error: "landarea_var_name" not defined in "AREA_INFO" namelist of GCM input config file', &
                                char_var=landarea_var_name)
        call add_path(landarea_file)

        !<><><><><><><><><><><><><><><><>!
        read(unit=fID, nml=LAND_CLIM_INFO)
        !<><><><><><><><><><><><><><><><>!
        call check_namelist_def('Error: "land_x_dim_name" not defined in "LAND_CLIM_INFO" namelist of GCM input config file', &
                                char_var=land_x_dim_name)
        call check_namelist_def('Error: "land_y_dim_name" not defined in "LAND_CLIM_INFO" namelist of GCM input config file', &
                                char_var=land_y_dim_name)
        call check_namelist_def('Error: "land_temp_var_name" not defined in "LAND_CLIM_INFO" namelist of GCM input config file', &
                                char_var=land_temp_var_name)
        call check_namelist_def('Error: "runoff_var_name" not defined in "LAND_CLIM_INFO" namelist of GCM input config file', &
                                char_var=runoff_var_name)
        n_undef = count(land_list_file(1:red_pspace_size) == UNDEFINED_VALUE_CHAR)
        if (n_undef > 0) then
            print *
            print *, 'Error: not enough file names were given in "LAND_CLIM_INFO" namelist (variable "land_list_file") of GCM'// &
                     'input config file'
            print *, n_undef, ' entries missing out of ', red_pspace_size
            stop 89
        end if
        do k = 1, red_pspace_size
            call add_path(land_list_file(k))
        end do

        !<><><><><><><><><><><><><><><><>!
        read(unit=fID, nml=GLOB_CLIM_INFO)
        !<><><><><><><><><><><><><><><><>!
        if (glob_temp_var_name == '') then
            glob_temp_var_name = UNDEFINED_VALUE_CHAR
        elseif (glob_temp_var_name /= UNDEFINED_VALUE_CHAR) then
            if (glob_list_file(1) == '' .or. glob_list_file(1) == UNDEFINED_VALUE_CHAR) then
                glob_x_dim_name = land_x_dim_name
                glob_y_dim_name = land_y_dim_name
                glob_list_file  = land_list_file
            else
                call check_namelist_def('Error: "glob_x_dim_name" not defined in "GLOB_CLIM_INFO" namelist of GCM input file', &
                                        char_var=glob_x_dim_name)
                call check_namelist_def('Error: "glob_y_dim_name" not defined in "GLOB_CLIM_INFO" namelist of GCM input file', &
                                        char_var=glob_y_dim_name)
                n_undef = count(glob_list_file(1:red_pspace_size) == UNDEFINED_VALUE_CHAR)
                if (n_undef > 0) then
                    print *
                    print *, 'Error: not enough file names were given in "GLOB_CLIM_INFO" namelist (variable "GLOB_list_file")'// &
                             ' of GCM input config file'
                    print *, n_undef, ' entries missing out of ', red_pspace_size
                    stop 89
                end if
                do k = 1,red_pspace_size
                    call add_path(glob_list_file(k))
                end do
            end if
        end if


    end subroutine


    !======================================================================!


    subroutine load_variable(internal_varname, varname, x_varname, y_varname,               &
                             single_input_file, multiple_input_file,                        &
                             varout1D, varout7D, x, y, xvec, yvec, landarea, totarea,       &
                             paramspace_filling_order, fill_missing, fillval_handling_option)
        use netcdf
        include 'shape.inc'
        integer, parameter:: npixel = nlon*nlat
        integer, parameter:: paramspace_size = nclimber*len_p1*len_p2*len_p3*len_p4*len_p5

        character(len=*), intent(in):: internal_varname
        character(len=vmaxlen), intent(in):: varname, x_varname, y_varname
        character(len=fmaxlen), intent(in), optional:: single_input_file
        character(len=fmaxlen), dimension(:), intent(in), optional:: multiple_input_file
        double precision, dimension(npixel), intent(inout), optional:: landarea
        double precision, dimension(npixel), intent(in), optional:: totarea
        double precision, dimension(npixel), intent(out), optional:: varout1D
        double precision, dimension(npixel, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5), intent(out), optional:: varout7D
        double precision, intent(out), optional:: x(nlon), y(nlat), xvec(:,:), yvec(:,:)
        integer, dimension(paramspace_size, 6), intent(in), optional:: paramspace_filling_order
        integer, intent(in), optional:: fillval_handling_option
        logical, intent(in), optional:: fill_missing

        character(len=fmaxlen):: loc_input_file
        character(len=8):: load_status
        double precision, dimension(npixel):: dummyvar, dummyvarout
        double precision:: loc_x(nlon), loc_y(nlat)
        double precision:: fillvalue
        character(len=vmaxlen):: var_units, vname
        character(len=1):: oper
        logical:: loc_fill_missing, axis_is_present
        integer:: ierr, loc_fillval_handling_option
        integer:: i1, i2, i3, i4, i5, i6, k, n, noperations, n_input


        ! How to handle missing-value
        if (present(fill_missing)) then
            loc_fill_missing = fill_missing
        else
            loc_fill_missing = .false. ! default behaviour
        end if
        !
        ! Note: checking missing-points is done only if fill_missing==.false. and "landarea" argument is given
        if (present(fillval_handling_option)) then
            loc_fillval_handling_option = fillval_handling_option
        else
            loc_fillval_handling_option = -1 ! default value => ask user
        end if


        ! Check consistency between 1D variable or 7D (ie, multiple CO2 levels/climatic parameters) variable
        ! => the right "input_file" and "var" arguments must be given.
        !
        if (present(varout1D) .and. present(single_input_file) .and. (.not. present(varout7D)) .and. &
            (.not. present(multiple_input_file)) .and. (.not. present(paramspace_filling_order))) then
            load_status = 'single  '
            n_input = 1
            axis_is_present = (present(x) .and. present(y))
        !
        elseif (present(varout7D) .and. present(multiple_input_file) .and. present(paramspace_filling_order) .and. &
               (.not. present(varout1D)) .and. (.not. present(single_input_file))) then
            load_status = 'multiple'
            n_input = size(multiple_input_file)
            if (present(xvec) .and. present(yvec)) then
                axis_is_present = .true.
                if (size(xvec,1) /= n_input) then
                    print *
                    print *, 'INTERNAL ERROR: size of "xvec" variable passed to subroutine "load_variable" in module'
                    print *, '"GMC_io_module" is inconsistent with number of input file ("multiple_input_file")'
                    stop 421
                end if
                if (size(yvec,1) /= n_input) then
                    print *
                    print *, 'INTERNAL ERROR: size of "yvec" variable passed to subroutine "load_variable" in module'
                    print *, '"GMC_io_module" is inconsistent with number of input file ("multiple_input_file")'
                    stop 421
                end if
                if (size(xvec,2) /= nlon) then
                    print *
                    print *, 'INTERNAL ERROR: length of 2nd dimension of "xvec" variable passed to subroutine "load_variable"'
                    print *, 'in module "GMC_io_module" inconsistent with length of x axis'
                    stop 421
                end if
                if (size(yvec,2) /= nlat) then
                    print *
                    print *, 'INTERNAL ERROR: length of 2nd dimension of "yvec" variable passed to subroutine "load_variable"'
                    print *, 'in module "GMC_io_module" inconsistent with length of y axis'
                    stop 421
                end if
            else
                axis_is_present = .false.
            end if
        !
        else
            print *
            print *, 'INTERNAL ERROR: inconsistent set of optional variable passed to the subroutine'
            print *, '"load_variable", in module "GCM_io_module"'
            stop 421
        end if


        ! Instead of reading just 1 variable name, they may be several variables,
        ! with arithmetic operation between them (addition or substraction)
        noperations = get_arithmetic_operations(varname)
        ! => write the scratch files:
        !        unit=334: list of individual variables
        !        unit=335: list of mathematics operators (+ or -)


        ! Loop on (potential) input files (1 per CO2/climatic parameters combination(
        do k = 1,n_input

            if (load_status == 'single  ') then
                loc_input_file = single_input_file
            else ! (load_status == multiple')
                loc_input_file = multiple_input_file(k)
            end if


            ! Initialization
            dummyvarout = 0d0

            ! Loop to get variables and perform arithmetic operations:
            do n = 1,noperations

                ! load current variable name (expect one name per operation (addition, subtraction))
                read(unit=334, fmt=*) vname

                if (loc_fill_missing) then
                    ! load variable
                    if (axis_is_present) then
                        call load_netcdf_dble_horiz(loc_input_file, x_varname, y_varname, vname, dummyvar, &
                                                    x=loc_x, y=loc_y, varunits=var_units, fillval=fillvalue, fillval_iostat=ierr)
                    else
                        call load_netcdf_dble_horiz(loc_input_file, x_varname, y_varname, vname, dummyvar, &
                                                    varunits=var_units, fillval=fillvalue, fillval_iostat=ierr)
                    end if
                    !
                    ! set var=0 on "missing" cells
                    if (ierr==NF90_NOERR) where (dummyvar==fillvalue) dummyvar = 0d0
                    !
                    ! check variable units
                    if (present(totarea)) then
                        call check_units(internal_varname, vname, dummyvar, var_units, totarea=totarea)
                    else
                        call check_units(internal_varname, vname, dummyvar, var_units)
                    end if
                    !
                else
                    ! load variable
                    if (axis_is_present) then
                        call load_netcdf_dble_horiz(loc_input_file, x_varname, y_varname, vname, dummyvar, &
                                                    x=loc_x, y=loc_y, varunits=var_units, fillval=fillvalue)
                    else
                        call load_netcdf_dble_horiz(loc_input_file, x_varname, y_varname, vname, dummyvar, &
                                                    varunits=var_units, fillval=fillvalue)
                    end if
                    !
                    ! check missingpoint
                    if (present(landarea)) then
                        call check_missingpoints(vname, loc_input_file, dummyvar, fillvalue, landarea, &
                                                 loc_fillval_handling_option)
                    end if
                    !
                    ! check variable units
                    if (present(totarea)) then
                        call check_units(internal_varname, vname, dummyvar, var_units, fillvalue=fillvalue, totarea=totarea)
                    else
                        call check_units(internal_varname, vname, dummyvar, var_units, fillvalue=fillvalue)
                    end if
                end if

                ! perform arithmetic operation
                read(unit=335, fmt=*) oper
                select case (oper)
                !
                    case ("+")
                        if (present(landarea)) then
                            where (landarea>0) dummyvarout = dummyvarout + dummyvar
                        else
                            dummyvarout = dummyvarout + dummyvar
                        end if
                !
                    case ("-")
                        if (present(landarea)) then
                            where (landarea>0) dummyvarout = dummyvarout - dummyvar
                        else
                            dummyvarout = dummyvarout - dummyvar
                        end if
                !
                end select

            end do


            ! put currently loaded variables in subroutine's output variable !
            ! -------------------------------------------------------------- !

            if (load_status == 'single  ') then

                varout1D = dummyvarout

                if (axis_is_present) then
                    x = loc_x
                    y = loc_y
                end if

            else !(load_status == 'multiple')

                i1 = paramspace_filling_order(k,1)
                i2 = paramspace_filling_order(k,2)
                i3 = paramspace_filling_order(k,3)
                i4 = paramspace_filling_order(k,4)
                i5 = paramspace_filling_order(k,5)
                i6 = paramspace_filling_order(k,6)

                varout7D(:, i1, i2, i3, i4, i5, i6) = dummyvarout

                if (axis_is_present) then
                    xvec(k,:) = loc_x
                    yvec(k,:) = loc_y
                end if
                
            end if

            ! -------------------------------------------------------------- !


            ! rewind scratch files (they need to be re-loaded for each climatic parameters combination)
            rewind(unit=334)
            rewind(unit=335)

        end do


        ! close scratch files
        close(unit=334)
        close(unit=335)


    end subroutine


    !======================================================================!


    function get_arithmetic_operations(varstring)
    ! => create scratch files 334 and 335 containing (respectively) the lists of variables and operators
    ! Note: the first operators (corresponding to the first variable) is automatically '+'
    ! return the number of variables/operations
        character(len=*), intent(in):: varstring
        integer:: get_arithmetic_operations
        character(len=1), dimension(2), parameter:: operators = (/'+', '-'/) ! Legal arithmetic operators
        integer:: i0, i1
        open(unit=334, status='scratch') ! individual variables
        open(unit=335, status='scratch') ! operators
        i0 = 1
        i1 = 1
        get_arithmetic_operations = 1
        write(unit=335, fmt='(A)') '+'
        do i1 = 1,len(varstring)
            if (any(operators == varstring(i1:i1))) then
                write(unit=334, fmt='(A)') trim(adjustl(varstring(i0:i1-1)))
                write(unit=335, fmt='(A)') varstring(i1:i1)
                i0 = i1+1
                get_arithmetic_operations = get_arithmetic_operations + 1
            end if
        end do
        write(unit=334, fmt='(A)') trim(adjustl(varstring(i0:)))
        rewind(unit=334)
        rewind(unit=335)
    end function


    !======================================================================!


    subroutine load_netcdf_dble_horiz(fname, x_varname, y_varname, varname, var, x, y, varunits, fillval, fillval_iostat)

        use netcdf
        use netcdf_io_module, only: nf90_check
        include 'shape.inc'
        integer, parameter:: npixel = nlon*nlat
        integer, parameter:: nmax = max(nlon, nlat)

        character(len=*), intent(in):: fname, x_varname, y_varname, varname
        double precision, dimension(npixel), intent(out):: var
        double precision, intent(out), optional:: x(nlon), y(nlat)
        character(len=*), intent(out), optional:: varunits
        double precision, intent(out), optional:: fillval
        integer, intent(out), optional:: fillval_iostat
        integer:: nx, ny, ndim, k, n
        integer:: ierr, fid, varid, xdimid, ydimid, shp_ix, shp_iy
        integer, parameter:: ndim_max = 30 ! note: assume no more than 30 dimensions for a given variable
        integer, dimension(ndim_max):: dimids, shp
        logical:: transp
        double precision, dimension(nmax,nmax):: loc_var

        ! open netCDF file
        ierr = nf90_open(fname, NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//trim(fname))

        ! Get axis data
        if (present(x)) then
            call load_axis(fname, fid, x_varname, xdimid, x)
        else
            call load_axis(fname, fid, x_varname, xdimid)
        end if
        if (present(y)) then
            call load_axis(fname, fid, y_varname, ydimid, y)
        else
            call load_axis(fname, fid, y_varname, ydimid)
        end if

        ! Get variable ID
        ierr = nf90_inq_varid(fid, varname, varid)
        call nf90_check(ierr, 'Error while inquiring variable "'//trim(varname)//'" ID in file "'//trim(fname)//'"')

        ! Get variable number of dimension
        ierr = nf90_inquire_variable(fid, varid, ndims=ndim)
        call nf90_check(ierr, 'Error while inquiring number of dim. of variable "'//trim(varname)//'" in file "'//trim(fname)//'"')

        if (ndim > ndim_max) then
            print *, 'ERROR: too may dimensions in variable "'//trim(varname)//'" of file "'//trim(fname)//'"'
            print *, 'Maximum supported number of dimensions is: ', ndim_max
            stop 17
        end if

        ! Get variable dimension IDs
        ierr = nf90_inquire_variable(fid, varid, dimids=dimids(1:ndim))
        call nf90_check(ierr, 'Error while inquiring dimensions IDs of variable "'//trim(varname)//'" in file "'//trim(fname)//'"')

        ! Get variable shape (inquire length of all dimensions)
        do k = 1,ndim
            ierr = nf90_inquire_dimension(fid, dimids(k), len=shp(k))
            call nf90_check(ierr, 'Error while inquiring length of dimension in file "'//trim(fname)//'"')
        end do

        ! Check variable shape (must have exactly 2 non degenerated dimensions)
        if (count(shp(1:ndim) > 1) /= 2) then
            print *
            print *, 'Error: variable "'//trim(varname)//'" of file "'//trim(fname)// &
                     '" must have exactly 2 non-degenerated (ie, size-1) dimensions'
            stop 17
        else
            k = 1
            do while (shp(k)==1)
                k = k + 1
            end do
            shp_ix = k
            nx = shp(k)
            k = k + 1
            do while (shp(k)==1)
                k = k + 1
            end do
            shp_iy = k
            ny = shp(k)
        end if

        ! Check that variable is defined on the given dimensions
        if (dimids(shp_ix)==xdimid .and. dimids(shp_iy)==ydimid) then
            transp = .false.
        elseif (dimids(shp_ix)==ydimid .and. dimids(shp_iy)==xdimid) then
            transp = .true.
        else
            print *
            print *
            print *, 'Error: variable "'//trim(varname)//'" of file "'//trim(fname)// &
            '" is not defined on the given dimensions "'//trim(x_varname)//'" and "'//trim(y_varname)//'"'
            stop 17
        end if

        ! load variable
        if (ndim>2) then
            print *
            print *,                 'Variable "'//trim(varname)//'" of file "'//trim(fname)//'"'
            write(*, fmt='(A,I2,A)') '    Ignore ', ndim-2, ' degenerated (size-1) dimension(s).'
        end if
        ! >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< !
        ierr = nf90_get_var(fid, varid, loc_var(1:nx,1:ny))
        call nf90_check(ierr, 'Error while getting variable "'//trim(varname)//'" of file "'//trim(fname)//'"')
        ! >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< !
        ! >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< !
        var = geographic_ravelling(loc_var(1:nx,1:ny), transp=transp)
        ! >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< >< !

        ! Get variable units (if asked)
        if (present(varunits)) then
            varunits = ''
            ierr = nf90_get_att(fid, varid, 'units', varunits)
            call nf90_check(ierr, 'Warning: unable to get attribute "units" of variable "'//trim(varname)//'" in file "' &
                                   //trim(fname)//'". Variable assumed to be dimensionless.', kill=.false.)
        end if

        ! Get variale fill-value (if asked)
        if (present(fillval)) then
            ierr = nf90_get_att(fid, varid, '_FillValue', fillval)
            if (present(fillval_iostat)) then
                fillval_iostat = ierr
            else
                call nf90_check(ierr, 'Error while getting attribute "_FillValue" of variable "'//trim(varname)//'" in file "' &
                                                                                                          //trim(fname)//'"')
            end if
        end if

        ! close netCDF file
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//trim(fname), kill=.false.)

    end subroutine


    !======================================================================!


    function geographic_ravelling(varin2D, transp)
        include 'shape.inc'
        integer, parameter:: npixel = nlon*nlat
        double precision, dimension(:,:), intent(in):: varin2D
        logical, intent(in):: transp
        !
        double precision, dimension(npixel):: geographic_ravelling
        !
        ! unravel 2D-arrays in 1D-array, the first dimension (x) begin the most rapidly varying => "natural" order
        if (transp) then
            geographic_ravelling = reshape(transpose(varin2D), shape=(/npixel/))
        else
            geographic_ravelling = reshape(varin2D, shape=(/npixel/))
        end if
        !
    end function


    !======================================================================!


    subroutine load_axis(fname, fid, axname, axdimid, ax)
        use netcdf
        use netcdf_io_module, only: nf90_check
        character(len=*), intent(in):: fname, axname
        integer, intent(in):: fid
        integer, intent(out):: axdimid
        double precision, dimension(:), intent(out), optional:: ax
        integer:: varid, ierr

        ! Inquire dim ID
        ierr = nf90_inq_dimid(fid, axname, axdimid)
        call nf90_check(ierr, 'Error while inquiring dimension "'//trim(axname)//'" ID in file "'//trim(fname)//'"')
        ! Inquire var ID
        ierr = nf90_inq_varid(fid, axname, varid)
        call nf90_check(ierr, 'Error while inquiring variable "'//trim(axname)//'" ID in file "'//trim(fname)//'"')
        ! Get variable
        if (present(ax)) then
            ierr = nf90_get_var(fid, varid, ax)
            call nf90_check(ierr, 'Error while getting variable "'//trim(axname)//'" in file "'//trim(fname)//'"')
        end if

    end subroutine


    !======================================================================!


    subroutine compare_units(units_string, known_units, conversion_factor, conversion_offset, passed)

        use physical_units, only: units

        character(len=*), intent(in):: units_string
        type(units), intent(in):: known_units
        double precision, intent(out):: conversion_factor, conversion_offset
        logical, intent(out):: passed
        integer:: k

        conversion_factor = 1d0
        conversion_offset = 0d0

        passed = .false.

        if (units_string == known_units%reference) then
            passed = .true.
        else
            do k = 1,known_units%naccepted
                if (units_string == known_units%accepted(k)%string) then
                    passed = .true.
                    conversion_factor = known_units%accepted(k)%conversion(1)
                    conversion_offset = known_units%accepted(k)%conversion(2)
                end if
            end do
        end if

    end subroutine

    !======================================================================!

    subroutine check_units(which_variable, varname, var, varunits, fillvalue, totarea, error_handling)

        use physical_units, only: units, area_units, fraction_units, temperature_units, runoff_units

        include 'shape.inc'
        integer, parameter:: npixel = nlon*nlat

        character(len=*), intent(in):: which_variable, varname, varunits
        double precision, dimension(npixel), intent(inout):: var
        double precision, dimension(npixel), intent(in), optional:: totarea
        double precision, intent(in), optional:: fillvalue
        integer, intent(in), optional:: error_handling
        type(units):: known_units
        double precision:: factor, offset
        logical:: multiply_by_area, passed, loop
        integer:: k, loc_error_handling, answer


        if (present(error_handling)) then
            loc_error_handling = error_handling
        else
            loc_error_handling = -1 ! default behaviour: ask user
        end if


        multiply_by_area = .false.

        ! case-dependent statement
        select case (which_variable)

            case ("area")
                known_units = area_units()

            case ("landarea")
                known_units = area_units()

            case ("temperature")
                known_units = temperature_units()

            case ("runoff")
                known_units = runoff_units()

            case default
                print *
                print *, 'INTERNAL ERROR in function "check_units" of module "GCM_io_module": unkown variable case "' &
                                                                                            //trim(which_variable)//'"'
                stop 421

        end select


        ! Compare given units to reference units, and get conversion values if not equals
        call compare_units(varunits, known_units, factor, offset, passed)
        ! in the case "landarea", also try "fraction" units
        if (which_variable == "landarea" .and. (.not. passed)) then
            call compare_units(varunits, fraction_units(), factor, offset, passed)
            multiply_by_area = passed
        end if



        ! Print message
        if (passed) then
            if (factor/=1 .or. offset/=0 .or. multiply_by_area) then
                print *
                print *, 'Automatic conversion of variable "'//trim(varname)//'":'
                print *, '    "'//trim(varunits)//'" => "'//trim(known_units%reference)//'"'
            end if
        else
            print *
            print *, 'WARNING: unkown units of variable "'//trim(varname)//'".'
            print *, '    got units:      "'//trim(varunits)//'"'
            print *, '    expected units: "'//trim(known_units%reference)//'"'
            ! Error handling:
            select case (loc_error_handling)
                case (-1) ! ask user interactively
                    print *
                    loop = .true.
                    do while (loop)
                        loop = .false.
                        print *, 'Enter one of the following options:'
                        print *, '    0: abort the program'
                        print *, '    2: ignore the issue and continue the execution'
                        read *, answer
                        select case (answer)
                            case (0); stop 2
                            case (2) ! do nothing
                            case default; loop=.true.
                        end select
                    end do

                case (0) ! abort the program
                    stop 2

                case (2) ! do nothing

                case default
                    print *
                    print *, 'INTERNAL ERROR: illegal error handling option:', error_handling
                    stop 421
            end select
        end if


        ! Conversion
        if (present(fillvalue)) then
            where (var/=fillvalue)  var = factor*(var + offset)
        else
            var = factor*(var + offset)
        end if

        ! Multiply value by cell area
        if (multiply_by_area) then
            if (present(totarea)) then
                print *, 'Variable "'//trim(varname)//'" multiplied by total area'
                if (present(fillvalue)) then
                    where (var/=fillvalue) var = var * totarea
                else
                    var = var * totarea
                end if
            else
                print *
                print *, 'INTERNAL ERROR: total area not passed to the function "check_units" in module "GCM_io_module".'
                print *, 'You have the right to be mad at the developer.'
                stop 421
            end if
        end if


    end subroutine


    !======================================================================!


    subroutine check_missingpoints(varname, filename, var, fillvalue, landarea, error_handling)

        include 'shape.inc'
        integer, parameter:: npixel = nlon*nlat
        !
        character(len=*), intent(in):: varname, filename
        double precision, dimension(npixel), intent(in):: var
        double precision, intent(in):: fillvalue
        double precision, dimension(npixel), intent(inout):: landarea
        integer, intent(in):: error_handling
        logical, dimension(npixel):: errormask
        double precision:: area_err
        integer:: nerr, answer
        logical:: loop

        errormask = (landarea>0 .and. var==fillvalue)
        nerr = count(errormask)
        area_err = sum(landarea, mask=errormask)

        ! If error found:
        if (nerr > 0) then

            print *
            write(*,'(A)')         ' WARNING: found missing-value on continental points (land area > 0) in variable "' &
                                    //trim(varname)//'" of file "'//trim(filename)//'"'
            write(*,'(A,I0,A,I0)') '     Number of problematic continent cells:  ', nerr, ' / ', count(landarea>0)
            write(*,'(A,E14.7)')   '     Total area of those cells (km2):        ', 1d-6*area_err
            write(*,'(A,E14.7)')   '     Which is a fraction of total land area: ', area_err/sum(landarea)

            ! Error handling:
            select case (error_handling)
                case (-1) ! ask user interactively
                    print *
                    loop = .true.
                    do while (loop)
                        loop = .false.
                        print *, 'Enter one of the following options:'
                        print *, '    0: abort the program'
                        print *, '    1: remove erratic points (set area=0)'
                        print *, '    2: ignore the issue and continue the execution'
                        read *, answer
                        select case (answer)
                            case (0); stop 2
                            case (1); where (errormask) landarea=0
                            case (2) ! do nothing
                            case default; loop=.true.
                        end select
                    end do

                case (0) ! abort the program
                    stop 2

                case (1) ! automatic correction
                    print *
                    print *, 'Automatic correction => remove erratic points (set area=0)'
                    where (errormask) landarea=0

                case (2) ! do nothing

                case default
                    print *
                    print *, 'INTERNAL ERROR: illegal error handling option:', error_handling
                    stop 421

            end select

        end if

    end subroutine


    !======================================================================!


    subroutine check_axis(all_x, all_y, error_handling)
        use io_module, only: raise_axis_error
        double precision, dimension(:,:), intent(in):: all_x, all_y ! dim #1: files (temp, runoff, ...), dim #2: axis
        integer, intent(in):: error_handling
        double precision:: daxis, max_mismatch
        logical, dimension(:), allocatable:: errormask
        integer:: n, k, axlen, nerr

        n = size(all_x, 1) ! also = size(all_y, 1)


        ! x axis
        ! ------

        axlen = size(all_x, 2)
        allocate(errormask(axlen))
        errormask = .false.

        daxis = abs(all_x(1,2) - all_x(1,1))
        max_mismatch = 0d0

        do k = 2,n
            errormask = (errormask  .or.  maxval(abs(all_x(k,:) - all_x(1,:)))/daxis > axis_relat_accuracy)
            max_mismatch = max(max_mismatch, maxval(abs(all_x(k,:) - all_x(1,:))))
        end do
        nerr = count(errormask)
        if (nerr > 0) then
            call raise_axis_error('x', nerr, axlen, max_mismatch, error_handling, ref_axis_message='between GCM input files')
        end if

        deallocate(errormask)


        ! y axis
        ! ------

        axlen = size(all_y, 2)
        allocate(errormask(axlen))
        errormask = .false.

        daxis = abs(all_y(1,2) - all_y(1,1))
        max_mismatch = 0d0

        do k = 2,n
            errormask = (errormask  .or.  maxval(abs(all_y(k,:) - all_y(1,:)))/daxis > axis_relat_accuracy)
            max_mismatch = max(max_mismatch, maxval(abs(all_y(k,:) - all_y(1,:))))
        end do
        nerr = count(errormask)
        if (nerr > 0) then
            call raise_axis_error('y', nerr, axlen, max_mismatch, error_handling, ref_axis_message='between GCM input files')
        end if

        deallocate(errormask)

    end subroutine


    !======================================================================!


    subroutine load_oce_ascii_variable(id, red_pspace_size, paramspace_filling_order, varout7D, varout8D, singleclim)

        include 'shape.inc'
        integer, parameter:: npixel = nlon*nlat
        integer, parameter:: paramspace_size = nclimber*len_p1*len_p2*len_p3*len_p4*len_p5

        integer, intent(in):: id, red_pspace_size
        integer, dimension(paramspace_size, 6), intent(in):: paramspace_filling_order
        double precision, dimension(nbasin, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5), intent(out), optional:: varout7D
        double precision, dimension(nbasin**2, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5), intent(out), optional:: varout8D
        logical, intent(out), optional:: singleclim

        character(len=fmaxlen):: loc_input_file
        double precision, dimension(nbasin**2):: dummyvar
        integer:: i1, i2, i3, i4, i5, i6, k, npts, ierr

        if (present(varout7D)) then
            npts = nbasin
        elseif (present(varout8D)) then
            npts = nbasin**2
        else
            print *
            print *, 'INTERNAL ERROR: in subroutine "load_oce_ascii_variable" of "GCM_io_module":'
            print *, 'either "varout7D" or "varout8D" output arguments must be present.'
            stop 421
        end if

        ! first reading (single climatic parameter)
        read(unit=id, fmt=*) dummyvar(1:npts)

        if (red_pspace_size > 1) then
            ! second reading (to determine if unique or multiple climatic paraters values are in file)
            read(unit=id, fmt=*, iostat=ierr) dummyvar(1:npts)
        else
            ierr = -1
        end if

        if (ierr<0) then ! end-of-file error => single value

            ! rewind file to load first (unique) value
            rewind(unit=id)
            if (present(varout7D)) then
                read(unit=id, fmt=*) varout7D(:,1,1,1,1,1,1)
            else
                read(unit=id, fmt=*) varout8D(:,1,1,1,1,1,1)
            end if

            if (present(singleclim)) singleclim = .true.

            print *
            print *
            print *, 'Loaded constant oceanic circulation'

        else

            rewind(unit=id)

            ! load entire file (if ierr was >0, the error will be displayed)
            do k = 1,red_pspace_size

                read(unit=id, fmt=*) dummyvar(1:npts)

                i1 = paramspace_filling_order(k,1)
                i2 = paramspace_filling_order(k,2)
                i3 = paramspace_filling_order(k,3)
                i4 = paramspace_filling_order(k,4)
                i5 = paramspace_filling_order(k,5)
                i6 = paramspace_filling_order(k,6)

                if (present(varout7D)) then
                    varout7D(:, i1, i2, i3, i4, i5, i6) = dummyvar(1:npts)
                else
                    varout8D(:, i1, i2, i3, i4, i5, i6) = dummyvar(1:npts)
                end if

            end do

            if (present(singleclim)) singleclim = .false.

            print *
            print *
            print *, 'Loaded climate-variable oceanic circulation'

        end if

    end subroutine


end module

