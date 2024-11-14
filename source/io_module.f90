module io_module
implicit none


! Structure-type to store information about netCDF output variables.
! ------------------------------------------------------------------
! <><><><><><><><><><> !
type netcdf_output_var
    integer :: id
    character(len=30) :: vartype
    character(len=100) :: vname, units
    character(len=500) :: long_name
    integer :: ndim_tot               ! -> Total number of dimension of the output file (not of the variable!)
    logical, dimension(16) :: def_dim ! -> Indicate for each dimension if the variable is define on the dimension
    double precision :: fillval       ! -> Type will be converted according to "vartype"
    logical :: writevar
end type
! <><><><><><><><><><> !


! Default values for undefined namelist variables (ie, values that are set before trying to read variable in namelists)
character(len=17), parameter :: UNDEFINED_VALUE_CHAR = '!#\_UNDEFINED_/#!'
integer, parameter ::           UNDEFINED_VALUE_INT  = -142857103
real, parameter ::              UNDEFINED_VALUE_REAL = -0.142857e33
double precision, parameter ::  UNDEFINED_VALUE_DBLE = -0.142857103693d99


! Default values of netCDF-related variables
character(len=*), parameter:: DEFAULT_FILLVAL_NAME = '_FillValue'
double precision, parameter:: DEFAULT_FILLVAL = 9.96921e+36
character(len=*), parameter:: DEFAULT_VAR_TYPE = 'real'


contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 subroutine load_netcdf(filename,varname,VAR, units, missval_name,missval, x_axis_ref, y_axis_ref, expected_units)
 ! load variables only on the first two dimensions, the other dimensions will be loaded
 ! at the end of the dimension. ie: start = (/1,1,n3,n4,.../) ; stop = (/n1,n2,n3,n4,.../)
 ! The 2D-variable is then unravelled in a 1D variable, the first dimension
 ! being the inner-loop unravelled one.
 use netcdf
 use netcdf_io_module, only: nf90_check

  character(len=*), intent(in):: filename, varname
  double precision, intent(out), dimension(:):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in),  optional:: missval_name
  double precision, intent(out), optional:: missval
  double precision, dimension(:), intent(in), optional:: x_axis_ref, y_axis_ref ! reference axis that current file should match
  character(len=*), intent(in), optional:: expected_units

  character(len=50):: missval_n
  integer:: fileid, varid, ierr, i, l, n
  integer, dimension(:), allocatable:: dimids, dimlen, stt, cnt
  character(len=50), dimension(:), allocatable:: dimname
  double precision, dimension(:), allocatable:: x, y

  integer, dimension(5):: ERROR_HANDLING_OPTION
  common /error/ ERROR_HANDLING_OPTION


!----------------!
! File openning: !
!----------------!

  ierr = nf90_open(filename, NF90_NOWRITE, fileid)
  call nf90_check(ierr,'Error while openning input file '//filename)


!-------------------!
! Loading variable: !
!-------------------!

  ! get variable ID and number of dimension:
  ierr = nf90_inq_varid( fileid, varname, varid )
  call nf90_check(ierr,'Error while getting identifiers of variable '//varname)
  ierr = nf90_inquire_variable( fileid, varid, ndims=n )
  call nf90_check(ierr,'Error while getting number of dimensions of variable '//varname)

  if (n<2) then
    print *
    print *, 'Error while loading variable'//varname
    print *, 'Two dimension at least are expected'
    stop 17
  end if

  allocate(dimids(n))
  allocate(dimlen(n))
  allocate(stt(n))
  allocate(cnt(n))
  allocate(dimname(n))

  ! get variable shape:
  ierr = nf90_inquire_variable( fileid, varid, dimids=dimids )
  call nf90_check(ierr,'Error while getting dimensions identifiers of variable '//varname)
  do i = 1,n
    ierr = nf90_inquire_dimension( fileid, dimids(i), dimname(i), dimlen(i) )
    call nf90_check(ierr,'Error while getting dimensions lengths and names of variable '//varname)
  end do

  ! print report:
  print *
  print *, 'netCDF variable: '//varname
  print *, '  - Loaded dimension:'
  print *, '      x1: ',dimname(1)
  print *, '      x2: ',dimname(2)
  if (n>3) print *, '  - unloaded dimension (get end of dim):'
  do i = 3,n
    print *, '        ', dimname(i)
  end do

  ! select data range (except 2nd dimension):
  stt(1) = 1
  stt(3:n) = dimlen(3:n)
  cnt(1) = dimlen(1)
  cnt(2:n)=1

  ! variable loading:
  if (dimlen(1)*dimlen(2) /= size(VAR)) then
    print *
    write(*,'(A)') &
      ' Error while loading variable '//varname
    write(*,'(A,I0,A,I0,A,I0)') &
      ' Inconsistent shape, get: ',dimlen(1),'x',dimlen(2),' = ',dimlen(1)*dimlen(2)
    write(*,'(A,I0)') &
      ' while expected:          ',size(VAR)
    stop 13
  else
    do i=1,dimlen(2) ! loop on second dimension
      stt(2) = i
      ierr = nf90_get_var( fileid, varid, VAR( 1+dimlen(1)*(i-1) : dimlen(1)*i ), start=stt, count=cnt )
      !         dimension unravelling ---------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    call nf90_check(ierr,'Error while getting variable '//varname)
    end do
  end if


!----------------------------!
! Loading attribute "units": !
!----------------------------!

  if (present(units)) then
    ierr = nf90_get_att( fileid, varid, 'units', units)
    call nf90_check(ierr,'Error while getting variable "'//trim(varname)//'" attribute "units"')
  end if


!-------------------------------------!
! Loading attributes "missing-value": !
!-------------------------------------!

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_FillValue'
    end if
    ierr = nf90_get_att( fileid, varid, missval_n(1:l), missval)
    call nf90_check(ierr,'Error while getting variable "'//trim(varname)//'" attribute '//missval_n(1:l))
  end if


!--------------------------------------------!
! Check axis consistency with reference axis !
!--------------------------------------------!

  if (present(x_axis_ref) .and. present(y_axis_ref)) then
      ! Assume: 1st dim = x axis, 2nd dim = y axis
      allocate(x(dimlen(1)))
      allocate(y(dimlen(2)))
      ierr = nf90_inq_varid(fileid, dimname(1), varid)
      call nf90_check(ierr,'Warning: unable to get identifier of variable '//dimname(1), kill=.false.)
      if (ierr==NF90_NOERR) then
          ierr = nf90_get_var(fileid, varid, x)
          call nf90_check(ierr,'Error while getting variable '//dimname(1))
          ierr = nf90_inq_varid(fileid, dimname(2), varid)
          call nf90_check(ierr,'Warning: unable to get identifier of variable '//dimname(2), kill=.false.)
          if (ierr==NF90_NOERR) then
              ierr = nf90_get_var(fileid, varid, y)
              call nf90_check(ierr,'Error while getting variable '//dimname(2))
              call check_axis('', x, y, x_axis_ref, y_axis_ref, ERROR_HANDLING_OPTION(1))
          else
              print *, 'Cannot check axis consistency'
          end if
      else
          print *, 'Cannot check axis consistency'
      end if
      deallocate(x)
      deallocate(y)
  end if


!---------------!
! File closing: !
!---------------!

  ierr = nf90_close(fileid)
  call nf90_check(ierr,'Error while closing input file '//filename)


  ! print report:
  if (present(expected_units)) then
      print *, '  - units: '//trim(units)//' (expected: '//trim(expected_units)//')'
  else
      print *, '  - units: '//trim(units)
  end if
  print *


 end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 subroutine load_netcdf_3D(filename,varname,VAR, units, missval_name,missval, expected_units)
 ! load variables only on the first three dimensions, the other dimensions will be loaded
 ! at the end of the dimension. ie: start = (/1,1,1,n4,.../) ; stop = (/n1,n2,n3,n4,.../)
 ! The first two dimensions of the 3D-variable are then unravelled, leading to a 2D variable.
 ! The first dimension being the inner-loop unravelled one.
 ! Then, the order of the 2 remaining dimensions is reversed.
 use netcdf
 use netcdf_io_module, only: nf90_check

  character(len=*), intent(in):: filename, varname
  double precision, intent(out), dimension(:,:):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in),  optional:: missval_name
  double precision, intent(out), optional:: missval
  character(len=*), intent(in), optional:: expected_units

  character(len=50):: missval_n
  integer:: fileid, varid, ierr, i, j, l, n
  integer, dimension(:), allocatable:: dimids, dimlen, stt, cnt
  character(len=50), dimension(:), allocatable:: dimname


!----------------!
! File openning: !
!----------------!

  ierr = nf90_open(filename, NF90_NOWRITE, fileid)
  call nf90_check(ierr,'Error while openning input file '//filename)


!-------------------!
! Loading variable: !
!-------------------!

  ! get variable ID and number of dimension:
  ierr = nf90_inq_varid( fileid, varname, varid )
  call nf90_check(ierr,'Error while getting identifiers of variable '//varname)
  ierr = nf90_inquire_variable( fileid, varid, ndims=n )
  call nf90_check(ierr,'Error while getting number of dimensions of variable '//varname)

  if (n<3) then
    print *
    print *, 'Error while loading variable '//varname
    print *, 'Three dimension at least are expected'
    stop 17
  end if

  allocate(dimids(n))
  allocate(dimlen(n))
  allocate(dimname(n))
  allocate(stt(n))
  allocate(cnt(n))

  ! get variable shape:
  ierr = nf90_inquire_variable( fileid, varid, dimids=dimids )
  call nf90_check(ierr,'Error while getting dimensions identifiers of variable '//varname)
  do i = 1,n
    ierr = nf90_inquire_dimension( fileid, dimids(i), dimname(i), dimlen(i) )
    call nf90_check(ierr,'Error while getting dimensions lengths and names of variable '//varname)
  end do

  ! print report:
  print *, 'netCDF variable: '//varname
  print *, '  - Loaded dimension:'
  print *, '      x1: ',dimname(1)
  print *, '      x2: ',dimname(2)
  print *, '      x3: ',dimname(3)
  if (n>4) print *, '  - unloaded dimension (get end of dim):'
  do i = 4,n
    print *, '        ', dimname(i)
  end do

  ! select data range (except "start" of 2nd and 3rd dimension):
  stt(1) = 1
  cnt(1) = dimlen(1)
  cnt(2:n) = 1
  stt(4:n) = dimlen(4:n)

  ! variable loading:
  if (dimlen(1)*dimlen(2) /= size(VAR,2) .or. dimlen(3) /= size(VAR,1) ) then
    print *
    write(*,'(A)') &
      ' Error while loading variable '//varname
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') &
      ' Inconsistent shape, get: ',dimlen(1),'x',dimlen(2),'x',dimlen(3),' => ',dimlen(1)*dimlen(2),'x',dimlen(3)
    write(*,'(A,I0,A,I0)') &
      ' while expected:          ',size(VAR,2),'x',size(VAR,1)
    stop 13
  else
    do i=1,dimlen(2)
      do j=1,dimlen(3)
        stt(2) = i
        stt(3) = j
        !   3rd dimension becomes 1st dimension: v
        ierr = nf90_get_var( fileid, varid, VAR( j , 1+dimlen(1)*(i-1) : dimlen(1)*i ), start=stt, count=cnt )
        !     1st-2nd dimension unravelling ---------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        call nf90_check(ierr,'Error while getting variable '//varname)
      end do
    end do
  end if


!----------------------------!
! Loading attribute "units": !
!----------------------------!

  if (present(units)) then
    ierr = nf90_get_att( fileid, varid, 'units', units)
    call nf90_check(ierr,'Error while getting variable "'//trim(varname)//'" attribute "units"')
  end if


!-------------------------------------!
! Loading attributes "missing-value": !
!-------------------------------------!

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_FillValue'
    end if
    ierr = nf90_get_att( fileid, varid, missval_n(1:l), missval)
    call nf90_check(ierr,'Error while getting variable "'//trim(varname)//'" attribute '//missval_n(1:l))
  end if


!---------------!
! File closing: !
!---------------!

  ierr = nf90_close(fileid)
  call nf90_check(ierr,'Error while closing input file '//filename)


  ! print report:
  if (present(expected_units)) then
      print *, '  - units: '//trim(units)//' (expected: '//trim(expected_units)//')'
  else
      print *, '  - units: '//trim(units)
  end if
  print *


 end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 subroutine load_netcdf_4D(filename,varname,VAR, units, missval_name,missval, expected_units)
 ! load variables only on the first four dimensions, the other dimensions will be loaded
 ! at the end of the dimension. ie: start = (/1,1,1,1,n5,.../) ; stop = (/n1,n2,n3,n4,n5,.../)
 ! The first two dimensions of the 4D-variable are then unravelled, leading to a 3D variable.
 ! The first dimension being the inner-loop unravelled one.
 ! Then, the order of the 3 remaining dimensions is reversed.
 use netcdf
 use netcdf_io_module, only: nf90_check

  character(len=*), intent(in):: filename, varname
  double precision, intent(out), dimension(:,:,:):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in),  optional:: missval_name
  double precision, intent(out), optional:: missval
  character(len=*), intent(in), optional:: expected_units

  character(len=50):: missval_n
  integer:: fileid, varid, ierr, i, j, k, l, n
  integer, dimension(:), allocatable:: dimids, dimlen, stt, cnt
  character(len=50), dimension(:), allocatable:: dimname


!----------------!
! File openning: !
!----------------!

  ierr = nf90_open(filename, NF90_NOWRITE, fileid)
  call nf90_check(ierr,'Error while openning input file '//filename)


!-------------------!
! Loading variable: !
!-------------------!

  ! get variable ID and number of dimension:
  ierr = nf90_inq_varid( fileid, varname, varid )
  call nf90_check(ierr,'Error while getting identifiers of variable '//varname)
  ierr = nf90_inquire_variable( fileid, varid, ndims=n )
  call nf90_check(ierr,'Error while getting number of dimensions of variable '//varname)

  if (n<4) then
    print *
    print *, 'Error while loading variable '//varname
    print *, 'Four dimension at least are expected'
    stop 17
  end if

  allocate(dimids(n))
  allocate(dimlen(n))
  allocate(dimname(n))
  allocate(stt(n))
  allocate(cnt(n))

  ! get variable shape:
  ierr = nf90_inquire_variable( fileid, varid, dimids=dimids )
  call nf90_check(ierr,'Error while getting dimensions identifiers of variable '//varname)
  do i = 1,n
    ierr = nf90_inquire_dimension( fileid, dimids(i), dimname(i), dimlen(i) )
    call nf90_check(ierr,'Error while getting dimensions lengths and names of variable '//varname)
  end do

  ! print report:
  print *, 'netCDF variable: '//varname
  print *, '  - Loaded dimension:'
  print *, '      x1: ',dimname(1)
  print *, '      x2: ',dimname(2)
  print *, '      x3: ',dimname(3)
  print *, '      x4: ',dimname(4)
  if (n>5) print *, '  - unloaded dimension (get end of dim):'
  do i = 5,n
    print *, '        ', dimname(i)
  end do

  ! select data range (except "start" of 2nd, 3rd and 4th dimension):
  stt(1) = 1
  cnt(1) = dimlen(1)
  cnt(2:n) = 1
  stt(5:n) = dimlen(5:n)

  ! variable loading:
  if (dimlen(1)*dimlen(2) /= size(VAR,3) .or. dimlen(3) /= size(VAR,2) .or. dimlen(4) /= size(VAR,1) ) then
    print *
    write(*,'(A)') &
      ' Error while loading variable '//varname
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
      ' Inconsistent shape, get: ',dimlen(1),'x',dimlen(2),'x',dimlen(3),'x',dimlen(4),&
      ' => ',dimlen(1)*dimlen(2),'x',dimlen(3),'x',dimlen(4)
    write(*,'(A,I0,A,I0,A,I0)') &
      ' while expected:          ',size(VAR,3),'x',size(VAR,2),'x',size(VAR,1)
    stop 13
  else
    do i=1,dimlen(2)
      do j=1,dimlen(3)
        do k=1,dimlen(4)
          stt(2) = i
          stt(3) = j
          stt(4) = k
          !   3rd dimension becomes 2nd dimension: ----|
          !   4th dimension becomes 1st dimension: v   v
          ierr = nf90_get_var( fileid, varid, VAR( k , j , 1+dimlen(1)*(i-1) : dimlen(1)*i ), start=stt, count=cnt )
          !         1st-2nd dimension unravelling ---------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          call nf90_check(ierr,'Error while getting variable '//varname)
        end do
      end do
    end do
  end if


!----------------------------!
! Loading attribute "units": !
!----------------------------!

  if (present(units)) then
    ierr = nf90_get_att( fileid, varid, 'units', units)
    call nf90_check(ierr,'Error while getting variable "'//trim(varname)//'" attribute "units"')
  end if


!-------------------------------------!
! Loading attributes "missing-value": !
!-------------------------------------!

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_FillValue'
    end if
    ierr = nf90_get_att( fileid, varid, missval_n(1:l), missval)
    call nf90_check(ierr,'Error while getting variable "'//trim(varname)//'" attribute '//missval_n(1:l))
  end if


!---------------!
! File closing: !
!---------------!

  ierr = nf90_close(fileid)
  call nf90_check(ierr,'Error while closing input file '//filename)


  ! print report:
  if (present(expected_units)) then
      print *, '  - units: '//trim(units)//' (expected: '//trim(expected_units)//')'
  else
      print *, '  - units: '//trim(units)
  end if
  print *


 end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 subroutine load_netcdf_dimvar(filename,varnameX,varnameY,varnameZ,varX,varY,varZ,unitsX,unitsY,unitsZ)

 use netcdf
 use netcdf_io_module, only: nf90_check

  character(len=*), intent(in):: filename
  character(len=*), intent(in), optional:: varnameX, varnameY, varnameZ
  double precision, intent(out), dimension(:), optional:: varX, varY, varZ
  character(len=*), intent(out), optional:: unitsX, unitsY, unitsZ
  integer:: fileid, varid, ierr


! File openning:

  ierr = nf90_open(filename, NF90_NOWRITE, fileid)
  call nf90_check(ierr,'Error while openning input file '//filename)


!-----------------------------------!
! Loading variables and attributes: !
!-----------------------------------!

  if (present(varnameX)) then
    ierr = nf90_inq_varid( fileid, varnameX, varid )
    call nf90_check(ierr,'Error while getting variable "'//trim(varnameX)//'" identifiant')
    if (present(varX)) then
      ierr = nf90_get_var( fileid, varid, varX )
      call nf90_check(ierr,'Error while getting variable '//varnameX)
    end if
    if (present(unitsX)) then
      ierr = nf90_get_att( fileid, varid, 'units', unitsX)
      call nf90_check(ierr,'Error while getting variable "'//trim(varnameX)//'" attribute "units"')
    end if
  end if
  if (present(varnameY)) then
    ierr = nf90_inq_varid( fileid, varnameY, varid )
    call nf90_check(ierr,'Error while getting  variable "'//trim(varnameY)//'" identifiant')
    if (present(varY)) then
      ierr = nf90_get_var( fileid, varid, varY )
      call nf90_check(ierr,'Error while getting variable '//varnameY)
    end if
    if (present(unitsY)) then
      ierr = nf90_get_att( fileid, varid, 'units', unitsY)
      call nf90_check(ierr,'Error while getting variable "'//trim(varnameY)//'" attribute "units"')
    end if
  end if
  if (present(varnameZ)) then
    ierr = nf90_inq_varid( fileid, varnameZ, varid )
    call nf90_check(ierr,'Error while getting  variable "'//trim(varnameZ)//'" identifiant')
    if (present(varZ)) then
      ierr = nf90_get_var( fileid, varid, varZ )
      call nf90_check(ierr,'Error while getting variable '//varnameZ)
    end if
    if (present(unitsZ)) then
      ierr = nf90_get_att( fileid, varid, 'units', unitsZ)
      call nf90_check(ierr,'Error while getting variable "'//trim(varnameZ)//'" attribute "units"')
    end if
  end if


!---------------!
! File closing: !
!---------------!

  ierr = nf90_close(fileid)
  call nf90_check(ierr,'Error while closing input file '//filename)



end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine check_landcells_generic(message, erratic_points, landarea, error_handling)

    character(len=*), intent(in):: message
    logical, dimension(:), intent(in):: erratic_points
    double precision, dimension(:), intent(inout):: landarea
    integer, intent(in):: error_handling
    logical, dimension(:), allocatable:: errormask
    double precision:: area_err, tot_landarea
    integer:: nerr, answer
    logical:: loop

    tot_landarea = sum(landarea)

    allocate(errormask(size(landarea)))

    errormask = (landarea>0 .and. erratic_points)

    nerr = count(errormask)
    area_err = sum(landarea, mask=errormask)

    ! If error found:
    if (nerr > 0) then

        print *
        write(*,'(A)')         ' '//trim(message)
        write(*,'(A,I0,A,I0)') '     Number of problematic continent cells:  ', nerr, ' / ', count(landarea>0)
        write(*,'(A,E14.7)')   '     Total area of those cells (km2):        ', 1d6*area_err ! Note: expect area from geoclim in 1e6km2
        write(*,'(A,E14.7)')   '     Which is a fraction of total land area: ', area_err/tot_landarea

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

    deallocate(errormask)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine check_landcells(varname, fillval, landarea, error_handling, var1D, var2D, axis)

    character(len=*), intent(in):: varname
    double precision, intent(in):: fillval
    double precision, dimension(:), intent(inout):: landarea
    integer, intent(in):: error_handling
    double precision, dimension(:), intent(in), optional:: var1D
    double precision, dimension(:,:), intent(in), optional:: var2D
    integer, intent(in), optional:: axis
    integer:: loc_axis


    ! Default value for non-pixel axis (CO2 levels, lithology, ...)
    if (present(axis)) then
        loc_axis = axis
    else
        loc_axis = 2
        ! Default value => Climate variable convention: 1st axis is pixels, 2nd is CO2 level
    end if

    if (present(var1D)) then
        call check_landcells_generic('WARNING: found missing values on continental cells of variable "'//trim(varname)//'"', &
                                    (var1D==fillval), landarea, error_handling)
    end if

    if (present(var2D)) then
        call check_landcells_generic('WARNING: found missing values on continental cells of variable "'//trim(varname)//'"', &
                                    any(var2D==fillval, dim=loc_axis), landarea, error_handling)
    end if

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine check_invalid(which_var, landarea, error_handling, var1D, var2D, var7D, axis)

    character(len=*), intent(in):: which_var
    double precision, dimension(:), intent(inout):: landarea
    integer, intent(in):: error_handling
    double precision, dimension(:), intent(inout), optional:: var1D
    double precision, dimension(:,:), intent(inout), optional:: var2D
    double precision, dimension(:,:,:,:,:,:,:), intent(inout), optional:: var7D
    integer, intent(in), optional:: axis
    logical, dimension(:), allocatable:: errormask
    double precision:: area_err, tot_landarea, example_invalid, replacement_value
    integer:: k, i1, i2, i3, i4, i5, nerr, answer, loc_axis
    logical:: loop

    if (.not. (present(var1D) .or. (present(var2D)) .or. (present(var7D)))) then
        print *
        print *, 'INTERNAL ERROR in subroutine "check_invalid" of module "io_module".'
        print *, 'No variable to check was given to the subroutine.'
        stop 421
    end if

    ! Default value for non-pixel axis (CO2 levels, lithology, ...)
    if (present(axis)) then
        loc_axis = axis
    else
        loc_axis = 2
        ! Default value => Climate variable convention: 1st axis is pixels, 2nd is CO2 level
    end if

    allocate(errormask(size(landarea)))

    select case (which_var)

        case ('runoff')
            ! Check if negative runoff
            if (present(var1D)) then
                errormask = (landarea>0 .and. var1D<0)
                example_invalid = minval(var1D, mask=errormask)
            elseif (present(var2D)) then
                errormask = (landarea>0 .and. any(var2D<0, dim=loc_axis))
                example_invalid = minval(minval(var2D, dim=loc_axis), mask=errormask)
            else ! => var7D. !!ASSUME GEOGRAPHIC AXIS IS #1!!
                errormask = (landarea>0 .and. any(any(any(any(any(any(var7D<0, dim=2), dim=2), dim=2), dim=2), dim=2), dim=2))
                example_invalid = minval(minval(minval(minval(minval(minval(minval(var7D,dim=2),dim=2),dim=2),dim=2),dim=2),dim=2),&
                                         mask=errormask)
            end if
            replacement_value = 0d0

        case ('slope')
            ! Check if negative or null slope
            if (present(var1D)) then
                errormask = (landarea>0 .and. var1D<=0)
                example_invalid = minval(var1D, mask=errormask)
                replacement_value = minval(var1D, mask=(landarea>0 .and. var1D>0))
            else ! var2D
                errormask = (landarea>0 .and. any(var2D<=0, dim=loc_axis))
                example_invalid = minval(minval(var2D, dim=loc_axis), mask=errormask)
                replacement_value = minval(minval(var2D, dim=loc_axis), mask=(landarea>0 .and. all(var2D>0, dim=loc_axis)))
            end if

        case ('lithology fraction')
            ! Check if sum of lithology is 1 (with 1d-6 accuracy)
            if (present(var1D)) then
                print *
                print *, 'INTERNAL ERROR in subroutine "check_invalid" of module "io_module".'
                print *, 'Unexpected 1-dimension lithology variable (should be nlitho x npixel)'
                stop 421
            else ! var2D
                errormask = (landarea>0 .and. abs(sum(var2D, dim=loc_axis) - 1d0)>1d-6)
                example_invalid = maxval(abs(sum(var2D, dim=loc_axis) - 1d0), mask=errormask)
            end if

        case default
            print *
            print *, 'INTERNAL ERROR in subroutine "check_invalid" of module "io_module".'
            print *, 'Unkown variable case "'//trim(which_var)//'"'
            stop 421

        end select

    tot_landarea = sum(landarea)

    nerr = count(errormask)
    area_err = sum(landarea, mask=errormask)

    ! If error found:
    if (nerr > 0) then

        print *
        select case (which_var)
            case ('runoff')
                write(*,'(A)')         ' WARNING: found negative runoff on continental cells.'
            case ('slope')
                write(*,'(A)')         ' WARNING: found null or negative slope on continental cells.'
            case ('lithology fraction')
                write(*,'(A)')         ' WARNING: found continental cells where sum of all lithology fractions is not 1.'
        end select
        write(*,'(A,I0,A,I0)') '     Number of problematic continent cells:  ', nerr, ' / ', count(landarea>0)
        write(*,'(A,E14.7)')   '     Total area of those cells (km2):        ', 1d6*area_err ! Note: expect area from geoclim in 1e6km2
        write(*,'(A,E14.7)')   '     Which is a fraction of total land area: ', area_err/tot_landarea
        if (which_var=='lithology fraction') then
            write(*,'(A,E14.7)')   '     Maximum deviation from 1 found: ', example_invalid
        else
            write(*,'(A,E14.7)')   '     Minimum '//trim(which_var)//' value found: ', example_invalid
        end if

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
                    select case (which_var)
                        case ('runoff')
                            print *, '    3: replace invalid values by 0'
                        case ('slope')
                            print *, '    3: replace invalid values by minimum positive slope:', replacement_value
                    end select
                    read *, answer
                    select case (answer)
                        case (0)
                            stop 2
                        case (1)
                            where (errormask) landarea=0
                        case (2)
                            ! do nothing
                        case (3)
                            if (which_var=='lithology fraction') then
                                loop=.true.
                            else
                                if (present(var1D)) then
                                    where (errormask) var1D=replacement_value
                                ! note: Runclimber convention -> 1st dimension is pixels, 2nd is CO2 levels, others are climatic parameters
                                elseif (present(var2D)) then
                                    do k = 1,size(var2D,2)
                                        where (errormask) var2D(:,k)=replacement_value
                                    end do
                                else ! => var7D
                                    do i5 = 1,size(var7D,7)
                                        do i4 = 1,size(var7D,6)
                                            do i3 = 1,size(var7D,5)
                                                do i2 = 1,size(var7D,4)
                                                    do i1 = 1,size(var7D,3)
                                                        do k = 1,size(var7D,2) ! Runclimber convention: 1st dimension is pixels, 2nd is CO2 levels
                                                            where (errormask) var7D(:,k,i1,i2,i3,i4,i5) = replacement_value
                                                        end do
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end if
                            end if
                        case default
                            loop=.true.
                    end select
                end do

            case (0) ! abort the program
                stop 2

            case (1) ! automatic correction
                print *
                print *, 'Automatic correction => remove erratic points (set area=0)'
                where (errormask) landarea=0

            case (2) ! do nothing

            case (3)
                print *
                select case (which_var)
                    case ('runoff')
                        print *, 'Automatic correction => replace invalid values by 0'
                    case ('slope')
                        print *, 'Automatic correction => replace invalid values minimum positive slope'
                    case ('lithology fraction')
                        print *, 'Bad error handling option "3": No replacement value for invalid lithology fractions sum'
                        stop 59
                end select
                if (present(var1D)) then
                    where (errormask) var1D=replacement_value
                end if
                if (present(var2D)) then
                    do k = 1,size(var2D,2) ! Runclimber convention: 1st dimension is pixels, 2nd is CO2 levels
                        where (errormask) var2D(:,k)=replacement_value
                    end do
                end if

            case default
                print *
                print *, 'INTERNAL ERROR: illegal error handling option:', error_handling
                stop 421

        end select

    end if

    deallocate(errormask)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine raise_axis_error(which_axis, nerr, axis_len, max_mismatch, error_handling, ref_axis_message)
    character(len=*), intent(in):: which_axis
    integer, intent(in):: nerr, axis_len
    double precision, intent(in):: max_mismatch
    integer, intent(in):: error_handling
    character(len=*), intent(in), optional:: ref_axis_message
    integer:: answer
    logical:: loop

    print *
    if (present(ref_axis_message)) then
        write(*,'(A)')         ' WARNING: found mismatch of '//which_axis//' axis with reference axis'
    else
        write(*,'(A)')         ' WARNING: found mismatch of '//which_axis//' axis '//trim(ref_axis_message)
    end if
        write(*,'(A,I0,A,I0)') '     Number of points concerned / axis length:  ', nerr, ' / ', axis_len
        write(*,'(A,E14.7)')   '     Maximum mismatch found (in axis units):    ', max_mismatch

    ! Error handling:
    select case (error_handling)
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

        case (1) ! asked for automatic correction
            print *
            print *, 'Bad error handling option "1": cannot remove points with axis mismatch.'
            stop 59

        case (2) ! do nothing

        case default
            print *
            print *, 'INTERNAL ERROR: illegal error handling option:', error_handling
            stop 421

    end select

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine check_axis(filename, x, y, xref, yref, error_handling)
    character(len=*), intent(in):: filename
    double precision, dimension(:), intent(in):: x, y, xref, yref
    integer, intent(in):: error_handling
    double precision:: daxis, max_mismatch
    integer:: nerr
    !
    double precision, parameter:: axis_relat_accuracy = 1d-5


    ! x axis
    daxis = abs(x(2) - x(1))
    nerr = count( abs(x - xref)  >  daxis*axis_relat_accuracy )
    max_mismatch = maxval(abs(x - xref))
    if (nerr > 0) then
        call raise_axis_error(trim(filename)//' x', nerr, size(x), max_mismatch, error_handling)
    end if

    ! y axis
    daxis = abs(y(2) - y(1))
    nerr = count( abs(y - yref)  >  daxis*axis_relat_accuracy )
    max_mismatch = maxval(abs(y - yref))
    if (nerr > 0) then
        call raise_axis_error(trim(filename)//' y', nerr, size(y), max_mismatch, error_handling)
    end if

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! Functions and subroutines specific to output_var structures and namelists !
! ========================================================================= !


  subroutine set_default_nml_values(vname, units, defdim, writevar, long_name, fillval, vartype)
    ! Set as "undefined value" the outvar info variables that will be read
    ! in a namelist in the config file

    character(len=30), intent(out), optional:: vartype 
    character(len=100), intent(out), optional:: vname, units
    character(len=500), intent(out), optional:: long_name
    integer, dimension(:), intent(out), optional:: defdim
    logical, intent(out), optional:: writevar
    double precision, intent(out), optional:: fillval

    if (present(vname))     vname     = UNDEFINED_VALUE_CHAR
    if (present(units))     units     = UNDEFINED_VALUE_CHAR
    if (present(defdim))    defdim    = UNDEFINED_VALUE_INT
    if (present(writevar))  writevar  = .true.
    if (present(long_name)) long_name = UNDEFINED_VALUE_CHAR
    if (present(fillval))   fillval   = DEFAULT_FILLVAL
    if (present(vartype))   vartype   = DEFAULT_VAR_TYPE

  end subroutine


  ! ---------- !


  function set_outvar_info(vname, units, def_dim, writevar, long_name, fillval, vartype)
    !
    ! put the output variable information read in the namelist in the netCDF output var structure,

    character(len=30), intent(in):: vartype 
    character(len=100), intent(in):: vname, units
    character(len=500), intent(in):: long_name
    integer, dimension(:), intent(in):: def_dim
    logical, intent(in):: writevar
    double precision, intent(in):: fillval
    !
    type(netcdf_output_var):: set_outvar_info
    !
    integer:: ndim

    set_outvar_info%writevar = writevar

    set_outvar_info%vname = vname
    set_outvar_info%units = units
    ndim = size(def_dim)
    set_outvar_info%ndim_tot = ndim
    set_outvar_info%def_dim(1:ndim) = (def_dim == 1) ! convert "1|0" in ".true.|.false."

    if (long_name/=UNDEFINED_VALUE_CHAR) then
      set_outvar_info%long_name = long_name
    else
      set_outvar_info%long_name = ''
    end if

    if (vartype/=UNDEFINED_VALUE_CHAR) then
      set_outvar_info%vartype = vartype
    else
      set_outvar_info%vartype = DEFAULT_VAR_TYPE
    end if

    if (fillval/=UNDEFINED_VALUE_DBLE) then
      set_outvar_info%fillval = fillval
    else
      set_outvar_info%fillval = DEFAULT_FILLVAL
    end if

  end function



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine check_namelist_def(message, char_var, int_var, real_var, dble_var, kill)

    character(len=*), intent(in):: message
    character(len=*), intent(in), optional:: char_var
    integer, intent(in), optional::          int_var
    real, intent(in), optional::             real_var
    double precision, intent(in), optional:: dble_var
    logical, optional, intent(in):: kill
    logical:: loc_kill

    if (present(kill)) then
      loc_kill = kill
    else
      loc_kill = .true.
    end if

    if (present(char_var)) then
      if (char_var == UNDEFINED_VALUE_CHAR) then
        print *
        print *, message
        if (loc_kill) stop 191
      end if
    end if

    if (present(int_var)) then
      if (int_var == UNDEFINED_VALUE_INT) then
        print *
        print *, message
        if (loc_kill) stop 191
      end if
    end if

    if (present(real_var)) then
      if (real_var == UNDEFINED_VALUE_REAL) then
        print *
        print *, message
        if (loc_kill) stop 191
      end if
    end if

    if (present(dble_var)) then
      if (dble_var == UNDEFINED_VALUE_DBLE) then
        print *
        print *, message
        if (loc_kill) stop 191
      end if
    end if

  end subroutine




end module
