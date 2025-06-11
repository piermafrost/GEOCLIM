module nf90check_mod
implicit none
contains
  subroutine nf90check(ierr, message)
  use netcdf
    integer, intent(in):: ierr
    character(len=*), intent(in), optional:: message
    if (ierr/=0) then
      if (present(message)) print *, message
      print *, nf90_strerror(ierr)
      stop
    end if
  end subroutine
end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program create_uniform_litho
! Create a netCDF file with uniform lithologial fraction of the prescribed
! classes. Dimension could be specified by "--size n,m" where n is the size of
! longitude dimension and m the size of latitude dimension. In that case, the
! lower-left corner is assumed to be {-180°E,-90°N}.
! Or, it can be specified by "--match fname". It will then match the first 2
! dimensions of the file 'fname'. The format of 'fname" must be netCDF.
! The remaining arguments are interpreted as fraction of given lithological
! classes. The program return an error if the sum of all litholgy proportion is
! not 1.
!
! Examples of use:
!
!    ./create_uniform_litho --size 48,40 -o "litho_frac_48x40.nc" 0.1 0.7 0.2
!
!    --> create a 48x40 file (longitude x latitude) with 3 lithological classes
!        whose proportions are (respectively) 0.1, 0.7 and 0.2
!        the lower-left corner is assumed to be {-180°E,-90°N}
!        the name of the output file is "litho_frac_48x40.nc"
!
!
!    ./create_uniform_litho --match ../INPUT/temperature_piCTRL.nc 0.5 0.5
!
!    ---> create a file matching the shape of ../INPUT/temperature_piCTRL.nc
!         with 2 lithological classes whose proportions are (respectively) 0.5
!         and 0.5
!         the default output file name is used
!
!

use netcdf
use nf90check_mod, only: nf90check
implicit none

character(len=*), parameter:: def_ofname='litho_frac.nc'

double precision, dimension(:), allocatable:: lon, lat, litho_frac
character(len=500):: ofname
character(len=100):: arg, fname, varname
logical:: got_shape
integer:: Narg, n, i,j,k, nlon,nlat, nlitho
integer:: ierr, fid, dimid(3), varid(4)


! Scanning arguments:
!--------------------

open(unit=11,status='scratch',action='readwrite')

got_shape = .false.
nlon = 0
nlat = 0

Narg = command_argument_count()

ofname = def_ofname
nlitho = 0
n = 0
do while (n < Narg)

  n = n+1

  call get_command_argument(n, arg)


  if (arg(1:7)=='--match') then

    if (got_shape) then
      print *, 'WARNING: conflicting prescribed shape. The shape of the last command will be considered'
      deallocate(lon)
      deallocate(lat)
    end if
    got_shape = .true.
    n = n+1
    call get_command_argument(n, fname)
    ierr = nf90_open(fname, NF90_NOWRITE, fid)
    call nf90check(ierr,'Error while openning file '//trim(fname))
    ierr = nf90_inquire_dimension(fid, 1, name=varname, len=nlon)
    call nf90check(ierr,'Error while inquiring 1st dimension (longitude) of '//trim(fname))
    allocate(lon(nlon))
    ierr = nf90_inq_varid(fid, varname, varid(1))
    call nf90check(ierr,'Error while inquiring ID of longitude variable '//trim(varname))
    ierr = nf90_get_var(fid, varid(1), lon)
    call nf90check(ierr,'Error while getting longitude variable '//trim(varname))
    ierr = nf90_inquire_dimension(fid, 2, name=varname, len=nlat)
    call nf90check(ierr,'Error while inquiring 1st dimension (latitude) of '//trim(fname))
    allocate(lat(nlat))
    ierr = nf90_inq_varid(fid, varname, varid(1))
    call nf90check(ierr,'Error while inquiring ID of latitude variable '//trim(varname))
    ierr = nf90_get_var(fid, varid(1), lat)
    call nf90check(ierr,'Error while getting latitude variable '//trim(varname))
    ierr = nf90_close(fid)
    call nf90check(ierr,'Error while closing file '//trim(fname))
    

  elseif (arg(1:6)=='--size') then

    if (got_shape) then
      print *, 'WARNING: conflicting prescribed shape. The shape of the last command will be considered'
      deallocate(lon)
      deallocate(lat)
    end if
    got_shape = .true.
    n = n+1
    call get_command_argument(n, fname)
    read(fname,*) nlon,nlat
    allocate(lon(nlon))
    allocate(lat(nlat))
    lon = (/ ( -180 + 360*(dble(i)-0.5)/nlon , i=1,nlon ) /)
    lat = (/ ( -90  + 180*(dble(j)-0.5)/nlat , j=1,nlat ) /)


  elseif (arg(1:2)=='-o') then
    n = n+1
    call get_command_argument(n, ofname)


  else

    nlitho = nlitho + 1
    write(unit=11,fmt=*) arg

  end if



end do


if (.not. got_shape) then
  print *, 'ERROR: no rule for variable shape (longitude and latitude dimension lengths)'
  stop
end if


allocate(litho_frac(nlitho))

rewind(unit=11)
k = 1
do while (k<=nlitho)
  read(unit=11,fmt=*) arg
  read(arg,fmt=*,iostat=ierr) litho_frac(k)
  if (ierr/=0) then
    print *, 'WARNING: unexpected argument ',arg
    print *, 'Will be ignored.'
    nlitho = nlitho - 1
  else
    if (litho_frac(k)<0 .or. litho_frac(k)>1) then
      print *, 'ERROR: lithology fraction must be between 0 and 1'
      stop
    end if
    k = k+1
  end if
end do

write(*,fmt='(A43,I3)') 'Number of lithological classe(s) detected: ',nlitho
do k = 1,nlitho
  write(*,fmt='(A21,I2.2,A1,F10.5)') '  Fraction of class #',k,':',litho_frac(k)
end do

if ( abs(sum(litho_frac(1:nlitho))-1) > 1d-12 ) then
  print *, 'ERROR: the sum of all lithology fraction must be 1'
  stop
end if



!!!!!!!!!!
! OUTPUT !
!!!!!!!!!!

! Output file creation:
ierr = nf90_create(ofname, NF90_CLOBBER, fid)
call nf90check(ierr,'Error while output file creation')

! Dimensions:
ierr = nf90_def_dim( fid, 'lon',     nlon,   dimid(1) )
call nf90check(ierr,'Error while defining dimension "lon"')
ierr = nf90_def_dim( fid, 'lat',     nlat,   dimid(2) )
call nf90check(ierr,'Error while defining dimension "lat"')
ierr = nf90_def_dim( fid, 'litho',   nlitho, dimid(3) )
call nf90check(ierr,'Error while defining dimension "litho"')

! Variables:
ierr = nf90_def_var(fid, 'lon',        NF90_FLOAT, dimid(1), varid(1))
call nf90check(ierr,'Error while defining variable "lon"')
ierr = nf90_def_var(fid, 'lat',        NF90_FLOAT, dimid(2), varid(2))
call nf90check(ierr,'Error while defining variable "lat"')
ierr = nf90_def_var(fid, 'litho',      NF90_INT,   dimid(3), varid(3))
call nf90check(ierr,'Error while defining variable "litho"')
ierr = nf90_def_var(fid, 'litho_frac', NF90_FLOAT, dimid,    varid(4))
call nf90check(ierr,'Error while defining variable "litho_frac"')

! Variables attributes:
! Names:
ierr = nf90_put_att(fid, varid(1) , 'name', 'longitude')
call nf90check(ierr,'Error while putting variable "lon" attribute "name"')
ierr = nf90_put_att(fid, varid(2) , 'name', 'latitude')
call nf90check(ierr,'Error while putting variable "lat" attribute "name"')
ierr = nf90_put_att(fid, varid(3) , 'name', 'litho')
call nf90check(ierr,'Error while putting variable "litho" attribute "name"')
! Axis:
ierr = nf90_put_att(fid, varid(1) , 'axis',      'X')
call nf90check(ierr,'Error while putting variable "lon" attribute "axis" ')
ierr = nf90_put_att(fid, varid(1) , 'nav_model', 'Default grid')
ierr = nf90_put_att(fid, varid(1) , 'long_name', 'Longitude')
ierr = nf90_put_att(fid, varid(2) , 'axis',      'Y')
call nf90check(ierr,'Error while putting variable "lat" attribute "axis" ')
ierr = nf90_put_att(fid, varid(2) , 'nav_model', 'Default grid')
ierr = nf90_put_att(fid, varid(2) , 'long_name', 'Latitude')
ierr = nf90_put_att(fid, varid(3) , 'long_name', 'Lithology')
! Units:
ierr = nf90_put_att(fid, varid(1) , 'units', 'degrees_east')
call nf90check(ierr,'Error while putting variable "lon" attribute "units" ')
ierr = nf90_put_att(fid, varid(2) , 'units', 'degrees_north')
call nf90check(ierr,'Error while putting variable "lat" attribute "units" ')
ierr = nf90_put_att(fid, varid(3) , 'units', '-')
call nf90check(ierr,'Error while putting variable "litho" attribute "-" ')
ierr = nf90_put_att(fid, varid(4) , 'units', '-')

! End of definition
ierr = nf90_enddef(fid)
call nf90check(ierr,'Error while end of definition of file'//trim(ofname))

! Writing:
ierr = nf90_put_var(fid, varid(1), real(lon))
call nf90check(ierr,'Error while writing variable "lon"')
ierr = nf90_put_var(fid, varid(2), real(lat))
call nf90check(ierr,'Error while writing variable "lat"')
ierr = nf90_put_var(fid, varid(3), (/(k,k=1,nlitho)/))
call nf90check(ierr,'Error while writing variable "litho"')
do j = 1,nlat
  do i = 1,nlon
    ierr = nf90_put_var(fid, varid(4), real(litho_frac(1:nlitho)), start=(/i,j,1/), count=(/1,1,nlitho/) )
    call nf90check(ierr,'Error while writing variable "litho_frac"')
  end do
end do

! Close file
ierr = nf90_close(fid)
call nf90check(ierr,'Error while closing file'//trim(ofname))


print *
print *, 'Sucessfully created netCDF file '//trim(ofname)


end program
