!=================================================================================!
!------------------ SUBROUTINE FOR LOADING DYNSOIL INTPUT DATA -------------------!
!=================================================================================!


module dynsoil_read_input_mod
implicit none

contains

subroutine dynsoil_read_input(IOcondID, &
                              lon, lat, area, xlevs, reg_z_s, reg_x_s, reg_tau_s, reg_z, reg_tau, slope, missingpoints, &
                              DYNS_restart_dim, DYNS_restart_var)

    use io_module, only: load_netcdf_dimvar, load_netcdf, load_netcdf_3D, load_netcdf_4D, check_axis, &
                         check_landcells, check_landcells_generic, check_invalid, check_namelist_def, UNDEFINED_VALUE_CHAR
    use dynsoil_create_init_condition, only: create_init_condition
    use utils, only: add_path

    use dynsoil_physical_parameters, only: nlon, nlat, nlitho, nDSlev, scaling_factor
    integer, parameter:: npxl = nlon*nlat
    integer, intent(in):: IOcondID
    double precision, intent(inout), dimension(nlon):: lon
    double precision, intent(inout), dimension(nlat):: lat
    double precision, intent(inout), dimension(npxl):: area
    double precision, intent(out), dimension(nDSlev):: xlevs
    double precision, intent(out), dimension(npxl):: slope
    double precision, intent(out), dimension(nlitho,npxl):: reg_z_s,reg_x_s,reg_tau_s
    double precision, intent(out), dimension(nDSlev,nlitho,npxl):: reg_z,reg_tau
    logical, intent(out), dimension(npxl):: missingpoints
    character(len=100), intent(out):: DYNS_restart_dim(4), DYNS_restart_var(5)
    !
    double precision, dimension(nlon):: loc_lon
    double precision, dimension(nlat):: loc_lat
    !
    integer, parameter:: flg=200, vlg=50, ulg=30, mlg=20
    integer, parameter:: nvar=9
    character(len=30) :: init_mode
    character(len=flg):: init_file, file_name
    character(len=vlg):: var_name
    character(len=ulg):: units(nvar)
    character(len=mlg):: fillval_name
    double precision  :: missval(nvar)
    !
    character(len=2):: num
    integer:: i, nlev
    !
    integer, dimension(5):: ERROR_HANDLING_OPTION
    common /error/ ERROR_HANDLING_OPTION

    ! Default name and value of missing-values (if unspecified in text input files)
    character(len=*), parameter:: defmissvalname = '_FillValue'
    double precision, parameter:: defmissval = 9.96921e+36


    ! Namelist declaration
    namelist /DYNSOIL_INIT_INFO/ init_mode, init_file, DYNS_restart_dim, DYNS_restart_var
    namelist /SLOPE_INFO/ file_name, var_name, fillval_name


    ! Set default fill_value (in case cannot be read in netCDF files)
    missval = defmissval



    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !%   Read restart information from main config file:   %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

    print *
    print *
    print *, 'Read DynSoil conditions'
    print *, '-----------------------'

    ! Default values of namelist variables
    init_mode           = UNDEFINED_VALUE_CHAR
    init_file           = UNDEFINED_VALUE_CHAR
    DYNS_restart_dim    = UNDEFINED_VALUE_CHAR
    DYNS_restart_var    = UNDEFINED_VALUE_CHAR

    ! <><><><><><><><><><><><><><><><><><><> !
    read(unit=IOcondID, nml=DYNSOIL_INIT_INFO)
    ! <><><><><><><><><><><><><><><><><><><> !

    call check_namelist_def('Error - in dynsoil_read_input.f90: variable "init_mode" from namelist "DYNSOIL_INIT_INFO" was not'// &
                            ' given in config/IO_CONDITIONS', char_var=init_mode)

    if (init_mode(1:8) == 'startup:') then

        print *
        print *, '    => "cold" restart (internally-generated)'
        missingpoints = .false.

    elseif(init_mode == 'restart') then

        print *
        print *, '  > Load initialization files'

        call check_namelist_def('Error - in dynsoil_read_input.f90: variable "init_file" from namelist'// &
                                ' "DYNSOIL_INIT_INFO" was not given in config/IO_CONDITIONS', char_var=init_file)
        do i = 1,4
            write(num, fmt="(I0)") i
            call check_namelist_def('Error - in dynsoil_read_input.f90: "DYNS_restart_dim('//trim(num)//')" from namelist'// &
                                    ' "DYNSOIL_INIT_INFO" was not given in config/IO_CONDITIONS', char_var=DYNS_restart_dim(i))
        end do
        do i = 1,5
            write(num, fmt="(I0)") i
            call check_namelist_def('Error - in dynsoil_read_input.f90: "DYNS_restart_var('//trim(num)//')" from namelist'// &
                                    ' "DYNSOIL_INIT_INFO" was not given in config/IO_CONDITIONS', char_var=DYNS_restart_var(i))
        end do

        ! Add GEOCLIM root path if needed
        call add_path(init_file)


        ! ------------------------------ !
        !          loading data:         !
        ! ------------------------------ !

        call load_netcdf_dimvar(init_file, varnameX=DYNS_restart_dim(1), varnameY=DYNS_restart_dim(2), &
                                varnameZ=DYNS_restart_dim(4), varX=loc_lon, varY=loc_lat, varZ=xlevs)
        ! if reference axis do not already exist, create then
        if (all(lon==0d0) .and. all(lat==0d0)) then
            lon = loc_lon
            lat = loc_lat
        else ! otherwise, check axis consistency
            call check_axis('DynSoil initialization file', loc_lon, loc_lat, lon, lat, ERROR_HANDLING_OPTION(1))
        end if
        call load_netcdf_3D(init_file, DYNS_restart_var(1), reg_z_s,   units(4), defmissvalname, missval(4), expected_units='m')
        call load_netcdf_3D(init_file, DYNS_restart_var(2), reg_x_s,   units(5), defmissvalname, missval(5), &
                            expected_units='dimensionless')
        call load_netcdf_3D(init_file, DYNS_restart_var(3), reg_tau_s, units(6), defmissvalname, missval(6), expected_units='yr')
        call load_netcdf_4D(init_file, DYNS_restart_var(4), reg_z,     units(7), defmissvalname, missval(7), expected_units='m')
        call load_netcdf_4D(init_file, DYNS_restart_var(5), reg_tau,   units(8), defmissvalname, missval(8), expected_units='yr')


        ! Multiply by scaling factor
        ! --------------------------

        where (reg_z_s/=missval(4))
            reg_z_s   = reg_z_s * scaling_factor
        end where
        where (reg_tau_s/=missval(6))
            reg_tau_s = reg_tau_s * scaling_factor
        end where
        where (reg_z/=missval(7))
            reg_z = reg_z * scaling_factor
        end where
        where (reg_tau/=missval(8))
            reg_tau = reg_tau * scaling_factor
        end where


        ! Get missing points:
        ! -------------------

        where( reg_z_s(1,:)   == missval(4) &
          .or. reg_x_s(1,:)   == missval(5) &
          .or. reg_tau_s(1,:) == missval(6) &
          .or. reg_z(1,1,:)   == missval(7) &
          .or. reg_tau(1,1,:) == missval(8) )
            missingpoints = .true.
        elsewhere
            missingpoints = .false.
        end where


    else

        print *
        print *, 'ERROR: invalid DynSoil initialization mode "'//trim(init_mode)//'".'
        print *, 'Legal ones are "restart", "startup:null" and "startup:eq".'
        stop 47

    end if



    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !%             Slope            %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

    print *
    print *
    print *, 'Read Slope info and load data'
    print *, '-----------------------------'
    print *

    ! Default namelist values
    file_name    = UNDEFINED_VALUE_CHAR
    var_name     = UNDEFINED_VALUE_CHAR
    fillval_name = defmissvalname

    ! <><><><><><><><><><><><><><><><> !
    read(unit=IOcondID, nml=SLOPE_INFO)
    ! <><><><><><><><><><><><><><><><> !

    call check_namelist_def('Error - in dynsoil_read_input.f90: variable "file_name" from namelist "SLOPE_INFO" was not'// &
                            ' given in config/IO_CONDITIONS', char_var=file_name)
    call check_namelist_def('Error - in dynsoil_read_input.f90: variable "var_name" from namelist "SLOPE_INFO" was not'// &
                            ' given in config/IO_CONDITIONS', char_var=var_name)

    ! Add GEOCLIM root path if needed
    call add_path(file_name)

    ! Load data
    if (all(lon==0d0) .and. all(lat==0d0)) then ! no reference axis
        call load_netcdf(file_name, var_name, slope, units(9), fillval_name, missval(9), expected_units='dimensionless')
    else ! check consistency with reference axis
        call load_netcdf(file_name, var_name, slope, units(9), fillval_name, missval(9), &
                         x_axis_ref=lon, y_axis_ref=lat, expected_units='dimensionless')
    end if
    nlev = size(xlevs)



    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !%       Input data checks      %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

    ! Check consistency between missing-points in GEOCLIM inputs (ie: continental area = 0)
    ! in dynsoil inputs (missingpoint = .true.) and in slope

    call check_landcells('slope', missval(9), area, ERROR_HANDLING_OPTION(2), var1D=slope)
    call check_landcells_generic('WARNING: found missing value on continental cells of DynSoil init. variables (h, x_s or tau_s)', &
                                 missingpoints, area, ERROR_HANDLING_OPTION(2))

    ! Check for null slope
    call check_invalid('slope', area, ERROR_HANDLING_OPTION(3), var1D=slope)

    where (area==0) missingpoints=.true.



    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% Create init condition if "startup" mode %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

    if (init_mode(1:7)=='startup') then
        call create_init_condition(init_mode(9:), missval(9), xlevs,reg_z_s,reg_x_s,reg_tau_s,reg_z,reg_tau, missingpoints)
        ! Notes:
        !   * remove "startup:" form 'init_mode' variable
        !   * scaling factor is already integrated in dynsoil physical parameters
        !     (for 'eq' initial conditions), and not needed for 'null' initial conditions
    end if



end subroutine

end module
