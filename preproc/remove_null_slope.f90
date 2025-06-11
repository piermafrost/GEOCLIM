! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Use: `./remove_null_slope slope_file_name`
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program remove_null_slope
use netcdf

integer:: nx, ny, narg, fid, vid, dimids(2), ierr
character(len=500):: filename
double precision, dimension(:,:), allocatable:: slope
double precision:: minslp

narg = command_argument_count()
if (narg >= 1) then
    call get_command_argument(1, filename)
    if (narg > 1) write(unit=*, fmt='(A,I0,A)') 'warning: ignore ', narg-1, ' extra argument(s) passed to the script'
else
    write(unit=*, fmt='(A)', advance='no') 'Name of "slope" netCDF file: '
    read(unit=*, fmt='(A)') filename
end if

! ============================================================================ !

ierr = nf90_open(filename, NF90_WRITE , fid)
if (ierr/=NF90_NOERR) print *, nf90_strerror(ierr)

ierr = nf90_inq_varid(fid, 'slope', vid)
if (ierr/=NF90_NOERR) print *, nf90_strerror(ierr)

ierr = nf90_inquire_variable(fid, vid, dimids=dimids)
if (ierr/=NF90_NOERR) print *, nf90_strerror(ierr)
ierr = nf90_inquire_dimension(fid, dimids(1), len=nx)
if (ierr/=NF90_NOERR) print *, nf90_strerror(ierr)
ierr = nf90_inquire_dimension(fid, dimids(2), len=ny)
if (ierr/=NF90_NOERR) print *, nf90_strerror(ierr)

!++++++++++++++++++++!
allocate(slope(nx,ny))
!++++++++++++++++++++!

ierr = nf90_get_var(fid, vid, slope)
if (ierr/=NF90_NOERR) print *, nf90_strerror(ierr)

!++++++++++++++++++++++++++++++++++++++++++++++++!
minslp =  minval(minval(slope, 1, slope>0), 1)
print *
print *, 'minimum non-zero slope value: ', minslp
print *
where (slope==0) slope = minslp
!++++++++++++++++++++++++++++++++++++++++++++++++!

ierr = nf90_put_var(fid, vid, slope)
if (ierr/=NF90_NOERR) print *, nf90_strerror(ierr)
ierr = nf90_close(fid)
if (ierr/=NF90_NOERR) print *, nf90_strerror(ierr)

deallocate(slope)

end program
