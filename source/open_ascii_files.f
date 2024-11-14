subroutine open_ascii_files(ID)
    use utils, only: add_path, read_comment
    use io_module, only: UNDEFINED_VALUE_CHAR, UNDEFINED_VALUE_DBLE, check_namelist_def
    implicit none
    integer, intent(in):: ID
    include 'combine_foam.inc'

    ! Local variables
    character(len=30):: ocean_temp_mode, sedim_transport_mode
    character(len=100):: var_name, fillval_name ! note : unused, but needed because present in LITHO_INFO namelist
    character(len=500):: file_name, COMB_var_init_file, grid_area_file, cont_area_file, &
                         temperature_file, runoff_file, interp_T_factor, interp_R_factor, GCM_input_condition_file, &
                         ocean_temp_file, box_volume_file,box_surf_file, box_surf_sedi_file, deep_box_file, &
                         appcont_box_file, thermocline_file, surface_box_file, epicont_box_file, polar_box_file, box_press_file, &
                         exchange_file, boundary_length_file, COMBINE_restart_file
    integer:: ierr
    double precision, dimension(nlitho):: singlepixel_lithofrac


    ! Namelists declaration
    ! ---------------------
    namelist /INIT_INFO/ COMB_init_mode, COMB_var_init_file
    namelist /CONT_INFO/ cont_input_mode, grid_area_file, cont_area_file, temperature_file, runoff_file, &
                         interp_T_factor, interp_R_factor, GCM_input_condition_file
    namelist /LITHO_INFO/ file_name, var_name, fillval_name, singlepixel_lithofrac
    namelist /VEGET_INFO/ file_name
    namelist /CLIMPARAM_INFO/ file_name
    namelist /COMBINE_INFO/ ocean_temp_mode, sedim_transport_mode, ocean_temp_file, &
                            box_volume_file, box_surf_file, box_surf_sedi_file, &
                            deep_box_file, appcont_box_file, thermocline_file, surface_box_file, epicont_box_file, &
                            polar_box_file, box_press_file, exchange_file, boundary_length_file
    namelist /RESTART_INFO/ COMBINE_restart_file, DynSoil_restart_file


    ! Default values of namelist variables
    ! ------------------------------------
    COMB_init_mode           = 'restart'
    COMB_var_init_file       = UNDEFINED_VALUE_CHAR
    cont_input_mode          = UNDEFINED_VALUE_CHAR
    grid_area_file           = UNDEFINED_VALUE_CHAR
    cont_area_file           = UNDEFINED_VALUE_CHAR
    temperature_file         = UNDEFINED_VALUE_CHAR
    runoff_file              = UNDEFINED_VALUE_CHAR
    interp_T_factor          = UNDEFINED_VALUE_CHAR
    interp_R_factor          = UNDEFINED_VALUE_CHAR
    GCM_input_condition_file = UNDEFINED_VALUE_CHAR
    singlepixel_lithofrac    = UNDEFINED_VALUE_DBLE
    ocean_temp_mode          = UNDEFINED_VALUE_CHAR
    sedim_transport_mode     = UNDEFINED_VALUE_CHAR
    ocean_temp_file          = UNDEFINED_VALUE_CHAR
    box_volume_file          = UNDEFINED_VALUE_CHAR
    box_surf_file            = UNDEFINED_VALUE_CHAR
    box_surf_sedi_file       = UNDEFINED_VALUE_CHAR
    deep_box_file            = UNDEFINED_VALUE_CHAR
    appcont_box_file         = UNDEFINED_VALUE_CHAR
    thermocline_file         = UNDEFINED_VALUE_CHAR
    surface_box_file         = UNDEFINED_VALUE_CHAR
    epicont_box_file         = UNDEFINED_VALUE_CHAR
    polar_box_file           = UNDEFINED_VALUE_CHAR
    box_press_file           = UNDEFINED_VALUE_CHAR
    exchange_file            = UNDEFINED_VALUE_CHAR
    boundary_length_file     = UNDEFINED_VALUE_CHAR
    COMBINE_restart_file     = UNDEFINED_VALUE_CHAR
    DynSoil_restart_file     = UNDEFINED_VALUE_CHAR



    !! =========================================================== !!
    !!                INPUT AND INITIALIZATION FILES               !!
    !! =========================================================== !!


    ! <><><><><><><><><><><><> !
    read(unit=ID, nml=INIT_INFO)
    ! <><><><><><><><><><><><> !

    if (COMB_init_mode == "restart") then
        call check_namelist_def('Error - in open_ascii_files.f: variable "COMB_var_init_file" was not given in'// &
                                'config/IO_CONDITIONS', char_var=COMB_var_init_file)
        call add_path(COMB_var_init_file)
        open (2  , status='old', action='read', file=COMB_var_init_file)
        !    ***********************************************************
    elseif (COMB_init_mode /= "coldstart") then
        print *
        print *, 'Error: illegal value for "COMB_init_mode" (read in config/IO_CONDITIONS): ', COMB_init_mode
        print *, 'Expect "restart" or "coldstart"'
        stop 31
    end if


    !--------!


    ! <><><><><><><><><><><><> !
    read(unit=ID, nml=CONT_INFO)
    ! <><><><><><><><><><><><> !

    call check_namelist_def('Error - in open_ascii_files.f: variable "cont_input_mode" was not given in config/IO_CONDITIONS', &
                            char_var=cont_input_mode)

    select case (cont_input_mode)

        case ("ascii")

            call check_namelist_def('Error - in open_ascii_files.f: variable "grid_area_file" was not given in '// &
                                    'config/IO_CONDITIONS', char_var=grid_area_file)
            call add_path(grid_area_file)
            open (301, status='old', action='read', file=grid_area_file)
            !    *******************************************************

            call check_namelist_def('Error - in open_ascii_files.f: variable "cont_area_file" was not given in '// &
                                    'config/IO_CONDITIONS', char_var=cont_area_file)
            call add_path(cont_area_file)
            open (7  , status='old', action='read', file=cont_area_file)
            !    *******************************************************

            !ccccccc +++++++++++++++++++++++++++++++++ 
            !       Climate under several CO2 levels
            !
            call check_namelist_def('Error - in open_ascii_files.f: variable "temperature_file" was not given in '// &
                                    'config/IO_CONDITIONS', char_var=temperature_file)
            call add_path(temperature_file)
            open (30 , status='old', action='read', file=temperature_file)
            !    *********************************************************
            !
            call check_namelist_def('Error - in open_ascii_files.f: variable "runoff_file" was not given in '// &
                                    'config/IO_CONDITIONS', char_var=runoff_file)
            call add_path(runoff_file)
            open (31 , status='old', action='read', file=runoff_file)
            !    ****************************************************
            !
            !! ALTERNATIVE: log(CO2)-interpolation on each grid point:
            !
            !call check_namelist_def('Error - in open_ascii_files.f: variable "interp_T_factor" was not given in'// &
            !                        'config/IO_CONDITIONS', char_var=interp_T_factor)
            !call add_path(interp_T_factor)
            !open (302 , status='old', action='read', file=interp_T_factor)
            !!    *********************************************************
            !
            !call check_namelist_def('Error - in open_ascii_files.f: variable "interp_R_factor" was not given in '// &
            !                        'config/IO_CONDITIONS', char_var=interp_R_factor)
            !call add_path(interp_R_factor)
            !open (303 , status='old', action='read', file=interp_R_factor)
            !!    *********************************************************

        case ("GCM")

            call check_namelist_def('Error - in open_ascii_files.f: variable "GCM_input_condition_file" was not given in'// &
                                    ' config/IO_CONDITIONS', char_var=GCM_input_condition_file)
            call add_path(GCM_input_condition_file)
            open (333, status='old', action='read', file=GCM_input_condition_file)
            !    *****************************************************************

        case default

            print *
            print *, 'Error: illegal value for "cont_input_mode" (read in config/IO_CONDITIONS): ', cont_input_mode
            print *, 'Expect "GCM" or "ascii"'
            stop 31

    end select


    !--------!


    ! Generic variables default value
    file_name = UNDEFINED_VALUE_CHAR

    ! <><><><><><><><><><><><><> !
    read(unit=ID, nml=LITHO_INFO)
    ! <><><><><><><><><><><><><> !

    if (file_name=='' .or. file_name==UNDEFINED_VALUE_CHAR) then

        call check_namelist_def('Error - in open_ascii_files.f: variable "singlepixel_lithofrac" was not given in'// &
                                ' config/IO_CONDITIONS', dble_var=singlepixel_lithofrac(1))
        !
        print *
        print *, '    Uniform lithology  => loaded from main config file'
        !
        do k = 1,nlitho
            litho_frac(k,:) = singlepixel_lithofrac(k)
        end do
        !
        uniform_lithology = .true.

    else

        uniform_lithology = .false.
        ! => lithology data will be load later
        
    end if


    !--------!


    if (coupling_veget) then

        ! Generic variables default value
        file_name = UNDEFINED_VALUE_CHAR

        ! <><><><><><><><><><><><><> !
        read(unit=ID, nml=VEGET_INFO)
        ! <><><><><><><><><><><><><> !

        call check_namelist_def('Error - in open_ascii_files.f: variable "file_name" from namelist "VEGET_INFO" was not given'// &
                                ' in config/IO_CONDITIONS', char_var=file_name)
        call add_path(file_name)
        open (48 , status='old', action='read', file=file_name)
        !    **************************************************

    end if


    !--------!


    if (nclimparam > 0) then

        ! Generic variables default value
        file_name = UNDEFINED_VALUE_CHAR

        ! <><><><><><><><><><><><><><><> !
        read(unit=ID, nml=CLIMPARAM_INFO)
        ! <><><><><><><><><><><><><><><> !

        call check_namelist_def('Error - in open_ascii_files.f: variable "file_name" from namelist "CLIMPARAM_INFO" was not'// &
                                ' given in config/IO_CONDITIONS', char_var=file_name)
        call add_path(file_name)
        open (49 , status='old', action='read', file=file_name)
        !    **************************************************
        call read_comment(49) ! only accept comments at beginning of climatic parameters file

    end if



    !! =========================================================== !!
    !!                     COMBINE INPUT FILES                     !!
    !! =========================================================== !!


    ! <><><><><><><><><><><><><><> !
    read(unit=ID, nml=COMBINE_INFO)
    ! <><><><><><><><><><><><><><> !

    if (ocean_temp_mode=='parametric') then

        ! parametric ocean temperature => write "parametric" in scratch input file (that will be read later)
        open(unit=32, status='scratch', action='readwrite')
        !    **********************************************
        write(unit=32, fmt=*) 'parametric'
        rewind(unit=32)
        !
        print *
        print *, '    Parametric ocean temperature => internally-defined relationship'

    else

        call check_namelist_def('Error - in open_ascii_files.f: variable "ocean_temp_file" was not given in config/IO_CONDITIONS', &
                                char_var=ocean_temp_file)
        call add_path(ocean_temp_file)
        open (32 , status='old', action='read', file=ocean_temp_file)
        !    ********************************************************

    end if

    !ccccccc +++++++++++++++++++++++++++++++++ 
    !       Earth physical dimensions

    call check_namelist_def('Error - in open_ascii_files.f: variable "box_volume_file" was not given in config/IO_CONDITIONS', &
                            char_var=box_volume_file)
    call add_path(box_volume_file)
    open (33 , status='old', action='read', file=box_volume_file)
    !    ********************************************************

    call check_namelist_def('Error - in open_ascii_files.f: variable "box_surf_file" was not given in config/IO_CONDITIONS', &
                            char_var=box_surf_file)
    call add_path(box_surf_file)
    open (34 , status='old', action='read', file=box_surf_file)
    !    ******************************************************

    call check_namelist_def('Error - in open_ascii_files.f: variable "box_surf_sedi_file" was not given in config/IO_CONDITIONS', &
                            char_var=box_surf_sedi_file)
    call add_path(box_surf_sedi_file)
    open (37 , status='old', action='read', file=box_surf_sedi_file)
    !    ***********************************************************


    !ccccccc+++++++++++++++++++++++++++++++++
    !       definition of boxes

    call check_namelist_def('Error - in open_ascii_files.f: variable "deep_box_file" was not given in config/IO_CONDITIONS', &
                            char_var=deep_box_file)
    call add_path(deep_box_file)
    open (35 , status='old', action='read', file=deep_box_file)
    !    ******************************************************

    call check_namelist_def('Error - in open_ascii_files.f: variable "box_press_file" was not given in config/IO_CONDITIONS', &
                            char_var=box_press_file)
    call add_path(box_press_file)
    open (38 , status='old', action='read', file=box_press_file)
    !    *******************************************************

    call check_namelist_def('Error - in open_ascii_files.f: variable "appcont_box_file" was not given in config/IO_CONDITIONS', &
                            char_var=appcont_box_file)
    call add_path(appcont_box_file)
    open (39 , status='old', action='read', file=appcont_box_file)
    !    *********************************************************

    call check_namelist_def('Error - in open_ascii_files.f: variable "thermocline_file" was not given in config/IO_CONDITIONS', &
                            char_var=thermocline_file)
    call add_path(thermocline_file)
    open (40 , status='old', action='read', file=thermocline_file)
    !    *********************************************************

    call check_namelist_def('Error - in open_ascii_files.f: variable "surface_box_file" was not given in config/IO_CONDITIONS', &
                            char_var=surface_box_file)
    call add_path(surface_box_file)
    open (41 , status='old', action='read', file=surface_box_file)
    !    *********************************************************

    call check_namelist_def('Error - in open_ascii_files.f: variable "exchange_file" was not given in config/IO_CONDITIONS', &
                            char_var=exchange_file)
    call add_path(exchange_file)
    open (42 , status='old', action='read', file=exchange_file)
    !    ******************************************************

    call check_namelist_def('Error - in open_ascii_files.f: variable "epicont_box_file" was not given in config/IO_CONDITIONS', &
                            char_var=epicont_box_file)
    call add_path(epicont_box_file)
    open (46 , status='old', action='read', file=epicont_box_file)
    !    *********************************************************

    call check_namelist_def('Error - in open_ascii_files.f: variable "polar_box_file" was not given in config/IO_CONDITIONS', &
                            char_var=polar_box_file)
    call add_path(polar_box_file)
    open (47 , status='old', action='read', file=polar_box_file)
    !    *******************************************************


    !ccccccc+++++++++++++++++++++++++++++++++
    !       choice of sediment transport scheme

    call check_namelist_def('Error - in open_ascii_files.f: variable "sedim_transport_mode" was not given in config/IO_CONDITIONS',&
                            char_var=sedim_transport_mode)

    if (sedim_transport_mode == 'uniform') then

        sed_routing_flag = .false.

    elseif (sedim_transport_mode == 'neighbour') then

        call check_namelist_def('Error - in open_ascii_files.f: variable "boundary_length_file" (for the sediment transport '// &
                                'scheme) was not given in config/IO_CONDITIONS. Consider using `sedim_transport_mode="uniform"`.', &
                                char_var=boundary_length_file)
        call add_path(boundary_length_file)
        open (43 , status='old', action='read', file=boundary_length_file)
        !    *************************************************************

        sed_routing_flag = .true.

    else

        print *
        print *, 'Error: illegal value for "sedim_transport_mode" (read in config/IO_CONDITIONS): ', sedim_transport_mode
        print *, 'Expect "uniform" or "neighbour"'
        stop 31

    end if



    !! =========================================================== !!
    !!                  RESTART AND OFFLINE FILES                  !!
    !! =========================================================== !!


    ! <><><><><><><><><><><><><><> !
    read(unit=ID, nml=RESTART_INFO)
    ! <><><><><><><><><><><><><><> !


    !ccccccc+++++++++++++++++++++++++++++++++ 
    !       geoclim restart file:

    call check_namelist_def('Error - in open_ascii_files.f: variable "COMBINE_restart_file" was not given in config/IO_CONDITIONS',&
                            char_var=COMBINE_restart_file)
    open (10 , status='REPLACE', action='write', file=trim(output_directory)//trim(COMBINE_restart_file)//trim(run_name))
    !    ****************************************************************************************************************


    if (coupling_dynsoil) then

        !ccccccc+++++++++++++++++++++++++++++++++ 
        !        DynSoil restart file:

        call check_namelist_def('Error - in open_input_files.f: variable "DynSoil_restart_file" was not given in'// &
                                ' config/IO_CONDITIONS', char_var=DynSoil_restart_file)
        ! note: "DynSoil_restart_file" is a global variable

    end if




end subroutine
