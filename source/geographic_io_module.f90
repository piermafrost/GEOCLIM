module geographic_io_module

! This module contains subrouting that read the IO information in main config
! file (IO_CONDITION) and load the corresponding boundary conditions netCDF
! files that are not GCM outputs (basin routing map, lithology, ...)


use netcdf
use netcdf_io_module, only: nf90_check
use io_module, only: check_landcells, check_axis, check_invalid, UNDEFINED_VALUE_CHAR, check_namelist_def
use utils, only: add_path
implicit none


contains



    subroutine load_lithology(ID, ref_x_axis, ref_y_axis, areaclimber, litho_frac)

        !include 'combine_foam.inc'
        include 'shape.inc'
        integer, parameter :: npixel=nlon*nlat

        integer, intent(in):: ID
        double precision, dimension(nlon), intent(in):: ref_x_axis
        double precision, dimension(nlat), intent(in):: ref_y_axis
        double precision, dimension(npixel), intent(inout) :: areaclimber
        double precision, dimension(nlitho,npixel), intent(out) :: litho_frac

        character(len=500):: file_name
        character(len=100):: var_name, fillval_name
        integer:: fileid, varid, ierr
        integer, dimension(3):: dimids, dimlen
        character(len=50), dimension(3):: dimname
        double precision:: fillval
        logical:: got_fillval
        double precision, dimension(nlon):: loc_x_axis
        double precision, dimension(nlat):: loc_y_axis
        double precision, dimension(nlitho):: singlepixel_lithofrac ! unused variable, but present in namelist
        integer :: i, j, k

        ! Common block (must be the same than in combine_foam.inc)
        integer, dimension(5) :: ERROR_HANDLING_OPTION
        common /error/ ERROR_HANDLING_OPTION

        ! Namelist declaration
        namelist /LITHO_INFO/ file_name, var_name, fillval_name, singlepixel_lithofrac



        !----------------------------------------------!
        ! Get input information from main config file: !
        !----------------------------------------------!

        ! variables default value
        file_name    = UNDEFINED_VALUE_CHAR
        var_name     = UNDEFINED_VALUE_CHAR
        fillval_name = UNDEFINED_VALUE_CHAR

        rewind(unit=ID)
        ! <><><><><><><><><><><><><> !
        read(unit=ID, nml=LITHO_INFO)
        ! <><><><><><><><><><><><><> !

        call check_namelist_def('Error - in geographic_io_module.f90: variable "file_name" from namelist "LITHO_INFO" was'// &
                                ' not given in config/IO_CONDITIONS', char_var=file_name)
        call check_namelist_def('Error - in geographic_io_module.f90: variable "var_name" from namelist "LITHO_INFO" was not'// &
                                ' given in config/IO_CONDITIONS', char_var=var_name)
        call check_namelist_def('Error - in geographic_io_module.f90: variable "fillval_name" from namelist "LITHO_INFO" was not'//&
                                ' given in config/IO_CONDITIONS', char_var=fillval_name)

        call add_path(file_name)


        ! ========== !


        print *
        print *
        print *, 'Load lithology map'
        print *, '------------------'
        print *


        !----------------!
        ! File openning: !
        !----------------!

        ierr = nf90_open(file_name, NF90_NOWRITE, fileid)
        call nf90_check(ierr, 'Error while openning input file '//trim(file_name))


        !-------------------!
        ! Loading variable: !
        !-------------------!

        ! get variable ID and number of dimension:
        ierr = nf90_inq_varid( fileid, var_name, varid )
        call nf90_check(ierr, 'Error while getting identifiers of variable '//trim(var_name))

        ! get variable fillvalue
        ierr = nf90_get_att( fileid, varid, fillval_name, fillval )
        call nf90_check(ierr, 'Warning: unable to get attribute "'//trim(fillval_name)//'" of variable '//trim(var_name), &
                        kill=.false.)
        got_fillval = (ierr==NF90_NOERR)

        ! get variable shape:
        ierr = nf90_inquire_variable( fileid, varid, dimids=dimids )
        call nf90_check(ierr, 'Error while getting dimensions identifiers of variable '//trim(var_name))
        do i = 1,3
            ierr = nf90_inquire_dimension( fileid, dimids(i), dimname(i), dimlen(i) )
            call nf90_check(ierr, 'Error while getting dimensions lengths and names of variable '//trim(var_name))
        end do

        ! print report:
        print *, 'netCDF variable: '//trim(var_name)
        print *, '  - Loaded dimension:'
        print *, '      lon:   ',trim(dimname(1))
        print *, '      lat:   ',trim(dimname(2))
        print *, '      litho: ',trim(dimname(3))
        print *

        ! variable loading:
        if (dimlen(1)/=nlon .or. dimlen(2)/=nlat .or. dimlen(3)/=nlitho) then
            print *, 'Error while loading variable '//trim(var_name)
            print *, 'Inconsistent shape, get: ',dimlen(1),'x',dimlen(2),'x',dimlen(3)
            print *, 'while expected:          ',nlon,'x',nlat,'x',nlitho
            stop 13
        else
            do j=1,nlat
                do k=1,nlitho
                    ierr = nf90_get_var(fileid, varid, litho_frac(k , 1+nlon*(j-1) : nlon*j), start=(/1,j,k/), count=(/nlon,1,1/))
                    !         dimension unravelling ------------------^^^^^^^^^^^^^^^^^^^^^
                    call nf90_check(ierr, 'Error while getting variable '//trim(var_name))
                end do
            end do
        end if


        !----------------------------------------!
        ! Check consistency with reference axis: !
        !----------------------------------------!

        if (.not. (all(ref_x_axis==0d0) .and. all(ref_y_axis==0d0))) then ! if reference axis are defined
            ! Assume: 1st dim = x axis, 2nd dim = y axis
            ierr = nf90_inq_varid(fileid, dimname(1), varid)
            call nf90_check(ierr, 'Warning: unable to get identifier of variable '//trim(dimname(1)), kill=.false.)
            if (ierr==NF90_NOERR) then
                ierr = nf90_get_var(fileid, varid, loc_x_axis)
                call nf90_check(ierr, 'Error while getting variable '//trim(dimname(1)))
                ierr = nf90_inq_varid(fileid, dimname(2), varid)
                call nf90_check(ierr, 'Warning: unable to get identifier of variable '//trim(dimname(2)), kill=.false.)
                if (ierr==NF90_NOERR) then
                    ierr = nf90_get_var(fileid, varid, loc_y_axis)
                    call nf90_check(ierr, 'Error while getting variable '//trim(dimname(2)))
                    call check_axis('lithology file', loc_x_axis, loc_y_axis, ref_x_axis, ref_y_axis, ERROR_HANDLING_OPTION(1))
                else
                    print *, 'Cannot check axis consistency'
                end if
            else
                print *, 'Cannot check axis consistency'
            end if
        end if


        !---------------!
        ! File closing: !
        !---------------!

        ierr = nf90_close(fileid)
        call nf90_check(ierr, 'Error while closing input file '//trim(file_name))


        !---------!
        ! Checks: !
        !---------!

        ! Missingpoints
        if (got_fillval) then
            call check_landcells('lithology fraction', fillval, areaclimber, ERROR_HANDLING_OPTION(2), var2D=litho_frac, axis=1)
        end if

        ! Fraction consistency
        call check_invalid('lithology fraction', areaclimber, ERROR_HANDLING_OPTION(4), var2D=litho_frac, axis=1)


    end subroutine



    ! ================================================================================================================ !



    subroutine load_landbasinmap(ID, ref_x_axis, ref_y_axis, areaclimber, indice_appcont, jbox_appcont, basinmap, uniform_routing)

        !include 'combine_foam.inc'
        include 'shape.inc'
        integer, parameter :: npixel=nlon*nlat

        integer, intent(in):: ID
        double precision, dimension(nlon), intent(in):: ref_x_axis
        double precision, dimension(nlat), intent(in):: ref_y_axis
        double precision, dimension(npixel), intent(inout) :: areaclimber
        integer, dimension(nbasin), intent(in) :: indice_appcont, jbox_appcont
        integer, dimension(npixel), intent(out) :: basinmap
        logical, intent(out):: uniform_routing

        character(len=500):: file_name
        character(len=100):: var_name, fillval_name
        integer:: fileid, varid, ierr
        integer, dimension(2):: dimids, dimlen
        character(len=50), dimension(2):: dimname
        integer:: fillval
        logical:: got_fillval, kill
        double precision, dimension(nlon):: loc_x_axis
        double precision, dimension(nlat):: loc_y_axis
        integer :: i, j, idx(nbasin)

        ! Common block (must be the same than in combine_foam.inc)
        integer, dimension(5) :: ERROR_HANDLING_OPTION
        common /error/ ERROR_HANDLING_OPTION

        ! Namelist declaration
        namelist /BASINMAP_INFO/ file_name, var_name, fillval_name, uniform_routing



        !----------------------------------------------!
        ! Get input information from main config file: !
        !----------------------------------------------!

        ! variables default value
        file_name    = UNDEFINED_VALUE_CHAR
        var_name     = UNDEFINED_VALUE_CHAR
        fillval_name = UNDEFINED_VALUE_CHAR
        !
        uniform_routing = .false.

        rewind(unit=ID)
        ! <><><><><><><><><><><><><><><> !
        read(unit=ID, nml=BASINMAP_INFO)
        ! <><><><><><><><><><><><><><><> !


        if (uniform_routing) then

            basinmap(:) = jbox_appcont(1)
            ! put all the land fluxes in the 1st "appcont" box
            ! => the global continental fluxes will after be redistributed in all epicontinental basins proportionnaly to their areas


        else

            call check_namelist_def('Error - in geographic_io_module.f90: variable "file_name" from namelist "BASINMAP_INFO" was'//&
                                    ' not given in in config/IO_CONDITIONS', char_var=file_name)
            call check_namelist_def('Error - in geographic_io_module.f90: variable "var_name" from namelist "BASINMAP_INFO" was'// &
                                    ' not given in config/IO_CONDITIONS', char_var=var_name)

            call add_path(file_name)


            ! ========== !


            print *
            print *
            print *, 'Load land basins map'
            print *, '--------------------'
            print *


            !----------------!
            ! File openning: !
            !----------------!

            ierr = nf90_open(file_name, NF90_NOWRITE, fileid)
            call nf90_check(ierr, 'Error while openning input file '//trim(file_name))


            !-------------------!
            ! Loading variable: !
            !-------------------!

            ! get variable ID and number of dimension:
            ierr = nf90_inq_varid( fileid, var_name, varid )
            call nf90_check(ierr, 'Error while getting identifiers of variable '//trim(var_name))

            ! get variable fillvalue
            if (fillval_name == UNDEFINED_VALUE_CHAR) then
                got_fillval = .false.
            else
                ierr = nf90_get_att( fileid, varid, fillval_name, fillval )
                call nf90_check(ierr, 'Warning: unable to get attribute "'//trim(fillval_name)//'" of variable '//trim(var_name), &
                                kill=.false.)
                got_fillval = (ierr==NF90_NOERR)
            end if

            ! get variable shape:
            ierr = nf90_inquire_variable( fileid, varid, dimids=dimids )
            call nf90_check(ierr, 'Error while getting dimensions identifiers of variable '//trim(var_name))
            do i = 1,2
                ierr = nf90_inquire_dimension( fileid, dimids(i), dimname(i), dimlen(i) )
                call nf90_check(ierr, 'Error while getting dimensions lengths and names of variable '//trim(var_name))
            end do

            ! print report:
            print *, 'netCDF variable: '//trim(var_name)
            print *, '  - Loaded dimension:'
            print *, '      lon:   ',trim(dimname(1))
            print *, '      lat:   ',trim(dimname(2))
            print *

            ! variable loading:
            if (dimlen(1)/=nlon .or. dimlen(2)/=nlat) then
                print *, 'Error while loading variable '//trim(var_name)
                print *, 'Inconsistent shape, get: ',dimlen(1),'x',dimlen(2)
                print *, 'while expected:          ',nlon,'x',nlat
                stop 13
            else
                do j=1,nlat
                    ierr = nf90_get_var(fileid, varid, basinmap(1+nlon*(j-1) : nlon*j), start=(/1,j/), count=(/nlon,1/))
                    !   dimension unravelling ------------------^^^^^^^^^^^^^^^^^^^^^
                    call nf90_check(ierr, 'Error while getting variable '//trim(var_name))
                end do
            end if


            !----------------------------------------!
            ! Check consistency with reference axis: !
            !----------------------------------------!

            if (.not. (all(ref_x_axis==0d0) .and. all(ref_y_axis==0d0))) then ! if reference axis are defined
                ! Assume: 1st dim = x axis, 2nd dim = y axis
                ierr = nf90_inq_varid(fileid, dimname(1), varid)
                call nf90_check(ierr, 'Warning: unable to get identifier of variable '//trim(dimname(1)), kill=.false.)
                if (ierr==NF90_NOERR) then
                    ierr = nf90_get_var(fileid, varid, loc_x_axis)
                    call nf90_check(ierr, 'Error while getting variable '//trim(dimname(1)))
                    ierr = nf90_inq_varid(fileid, dimname(2), varid)
                    call nf90_check(ierr, 'Warning: unable to get identifier of variable '//trim(dimname(2)), kill=.false.)
                    if (ierr==NF90_NOERR) then
                        ierr = nf90_get_var(fileid, varid, loc_y_axis)
                        call nf90_check(ierr, 'Error while getting variable '//trim(dimname(2)))
                        call check_axis('basinmap file', loc_x_axis, loc_y_axis, ref_x_axis, ref_y_axis, ERROR_HANDLING_OPTION(1))
                    else
                        print *, 'Cannot check axis consistency'
                    end if
                else
                    print *, 'Cannot check axis consistency'
                end if
            end if


            !---------------!
            ! File closing: !
            !---------------!

            ierr = nf90_close(fileid)
            call nf90_check(ierr, 'Error while closing input file '//trim(file_name))


            !---------!
            ! Checks: !
            !---------!

            ! Missingpoints
            if (got_fillval) then
                call check_landcells('land basins map',dble(fillval), areaclimber, ERROR_HANDLING_OPTION(2), var1D=dble(basinmap))
            end if

            ! Consistency with COMBINE boxes with continental inputs ("appcont" boxes):
            ! =========================================================================

            idx(:) = 0
            kill = .false.

            do j = 1,npixel
                if (areaclimber(j) /= 0d0) then
                    if (basinmap(j)<1 .or. basinmap(j)>nbasin-1) then
                        write(*,'(A,I0)')   ' Error: found out-of-range value in loaded land basin map: ',basinmap(j)
                        write(*,'(A,I0,A)') ' (current GEOCLIM oceanic basins range from #1 to #',nbasin-1,')'
                        stop 18
                    else
                        idx(basinmap(j)) = 1
                    end if
                end if
            end do

            if (any(idx/=indice_appcont)) then
                if (any(indice_appcont==1 .and. idx==0)) then
                    print *, 'Warning: there are oceanic basins indicated as "land input" ("appcont") basins'
                    print *, 'by current COMBINE input files, that do not receive continental fluxes according'
                    print *, 'the loaded land basin map.'
                    print *
                end if
                if (any(indice_appcont==0 .and. idx==1)) then
                    print *, 'Error: inconsistency between COMBINE inputs and continental routing input.'
                    print *, 'The loaded land basin map routes some continental fluxes into basins that are'
                    print *, 'indicated as "no land input" ("noappcont") by current COMBINE input files.'
                    print *
                    kill = .true.
                end if
                print *, '  "land input" ("appcont") basins according to COMBINE input files (1|0):'
                write(*, fmt='(A5)', advance='no') '     '
                do i = 1,nbasin-1
                    write(*, fmt='(I0,A1)', advance='no') indice_appcont(i), ' '
                end do
                print *
                print *, '  basins where continental fluxes are routed according to "basinmap" input file:'
                write(*, fmt='(A5)', advance='no') '     '
                do i = 1,nbasin-1
                    write(*, fmt='(I0,A1)', advance='no') idx(i), ' '
                end do
                print *
                print *
                if (kill) stop 18
            end if


        end if


    end subroutine



end module

