module local_functions
implicit none
contains

  subroutine mult_noise(x,ampl)
  !
    real,intent(out):: x
    real,intent(in):: ampl
    real:: logampl
  !
    logampl = log(ampl)
    call random_number(x)
    x = (2*x - 1)  *  logampl
    x = exp(x)
  !
  end subroutine mult_noise


  subroutine init_random_seed()
  !
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed
  !
    call random_seed(size = n)
    allocate(seed(n))
  !
    call system_clock(count=clock)
  !
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
  !
    deallocate(seed)
  !
  end subroutine init_random_seed


  subroutine check(ierr)
    use netcdf
    integer, intent(in):: ierr
    if (ierr/=NF90_NOERR) then
      print *, nf90_strerror(ierr)
      !stop
    end if
  end subroutine

end module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program create_slope
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! create a random slope matching geographic forcings from the area (or       !!
!! landfrac) variable in a "geographic" storing netCDF file                   !! 
!! The generated slope is pseudo-uniform: mean value * multiplicative noise   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use netcdf
use local_functions, only: mult_noise, init_random_seed, check
implicit none


! ================================================================================================ !
!                                          USER PARAMETERS                                         !
! ================================================================================================ !
  integer, parameter:: nx=96, ny=96
  character(len=*), parameter:: geofile='../INPUT/IPSL/CTRL_CM5A2/piControl_SE_2750_2849_1Y_sechiba.nc'
  character(len=*), parameter:: xvar='lon', yvar='lat', areavar='Areas'
  character(len=*), parameter:: fileout='test_uniform_slope.nc'
! ================================================================================================ !
  real, parameter:: NOISE_AMPLITUDE = 5.
  real, parameter:: MEANSLOPE = 0.03354 & ! 0.03354 is present-day Earth mean slope
                               *  2.*log(NOISE_AMPLITUDE)/(NOISE_AMPLITUDE-1./NOISE_AMPLITUDE)
! ================================================================================================ !


  integer:: i,j, ierr, gfid,ofid,gdimid(2),odimid(2),gvarid(3),ovarid(3)
  real:: lon(nx), lat(ny)
  real, dimension(nx,ny):: area, slope
  real:: noise, missval
  
  character(len=*), parameter:: missvalname='_FillValue'


  ierr = nf90_open(   geofile , NF90_NOWRITE , gfid )  ;  call check(ierr)
  ierr = nf90_create( fileout , NF90_CLOBBER , ofid )  ;  call check(ierr)

  ierr = nf90_inq_varid( gfid , xvar    , gvarid(1) )  ;  call check(ierr)
  ierr = nf90_inq_varid( gfid , yvar    , gvarid(2) )  ;  call check(ierr)
  ierr = nf90_inq_varid( gfid , areavar , gvarid(3) )  ;  call check(ierr)

  ierr = nf90_get_var( gfid , gvarid(1) , lon )   ;  call check(ierr)
  ierr = nf90_get_var( gfid , gvarid(2) , lat )   ;  call check(ierr)
  ierr = nf90_get_var( gfid , gvarid(3) , area )  ;  call check(ierr)
  !
  ! Default missing-value:
  missval = 9.96921e+36
  ierr = nf90_get_att( gfid , gvarid(3) , missvalname  , missval )  ;  call check(ierr)

  !---!

  ierr = nf90_def_dim( ofid , xvar , nx , odimid(1) )  ;  call check(ierr)
  ierr = nf90_def_dim( ofid , yvar , ny , odimid(2) )  ;  call check(ierr)

  ierr = nf90_def_var( ofid , xvar , NF90_REAL , odimid(1) , ovarid(1) )       ;  call check(ierr)
  ierr = nf90_copy_att( gfid , gvarid(1) , 'axis'       , ofid , ovarid(1) )    ;  call check(ierr)
  ierr = nf90_copy_att( gfid , gvarid(1) , 'nav_model'  , ofid , ovarid(1) )    ;  call check(ierr)
  ierr = nf90_copy_att( gfid , gvarid(1) , 'name'       , ofid , ovarid(1) )    ;  call check(ierr)
  ierr = nf90_copy_att( gfid , gvarid(1) , 'units'      , ofid , ovarid(1) )    ;  call check(ierr)
  !
  ierr = nf90_def_var( ofid , yvar , NF90_REAL , odimid(2) , ovarid(2) )       ;  call check(ierr)
  ierr = nf90_copy_att( gfid , gvarid(2) , 'axis'       , ofid , ovarid(2) )    ;  call check(ierr)
  ierr = nf90_copy_att( gfid , gvarid(2) , 'nav_model'  , ofid , ovarid(2) )    ;  call check(ierr)
  ierr = nf90_copy_att( gfid , gvarid(2) , 'name'       , ofid , ovarid(2) )    ;  call check(ierr)
  ierr = nf90_copy_att( gfid , gvarid(2) , 'units'      , ofid , ovarid(2) )    ;  call check(ierr)
  !
  ierr = nf90_def_var( ofid , 'slope' , NF90_REAL , odimid  , ovarid(3) )       ;  call check(ierr)
  ierr = nf90_put_att( gfid , ovarid(3) , 'name'  , 'slope' )                   ;  call check(ierr)
  ierr = nf90_put_att( ofid , ovarid(3) , 'units' , 'm/m' )                     ;  call check(ierr)
  ierr = nf90_put_att( ofid , ovarid(3) , missvalname  , missval )              ;  call check(ierr)

  ierr = nf90_close( gfid )

  ierr = nf90_enddef( ofid )

  ierr = nf90_put_var( ofid , ovarid(1) , lon )   ;  call check(ierr)
  ierr = nf90_put_var( ofid , ovarid(2) , lat )   ;  call check(ierr)

  !---------------------------------------------------------------------!

  call init_random_seed()

  do j = 1,ny
    do i = 1,nx

      if (area(i,j)/=missval .and. area(i,j)/=0.) then
         call mult_noise(noise,NOISE_AMPLITUDE)
!        ******************************
         slope(i,j) = MEANSLOPE * noise
!        ******************************
      else
        slope(i,j) = missval
      end if

    end do
  end do

  !---------------------------------------------------------------------!

  ierr = nf90_put_var( ofid , ovarid(3) , slope )  ;  call check(ierr)

  ierr = nf90_close( ofid )  ;  call check(ierr)


end program

