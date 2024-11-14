module geographic_write_output_mod
implicit none

contains

subroutine geographic_write_output(ofile_name, time_dimname, GEOG_outvar_info, t, cpvec, &
                                   temp, runoff, wth_allsil, wth_litho_wgh, wth_litho, ker_wth, bioC_exp, phos_wth, pyr_wth)
  use netcdf
  use io_module,        only: netcdf_output_var
  use netcdf_io_module, only: open_file, close_file, inquire_dim, put_var
  include 'shape.inc' ! => nlon, nlat, nlitho
  include 'output_size.inc' ! => nGEOGoutvar
  character(len=*):: ofile_name, time_dimname
  type(netcdf_output_var), dimension(nGEOGoutvar), intent(in):: GEOG_outvar_info
  double precision, intent(in):: t
  double precision, dimension(5), intent(in):: cpvec
  double precision, dimension(nlon*nlat), intent(in):: temp,runoff, wth_allsil, ker_wth, bioC_exp, phos_wth, pyr_wth 
  double precision, dimension(nlitho,nlon*nlat), intent(in):: wth_litho_wgh, wth_litho
  integer:: i, ierr, nt, fid, timevarid, dimid

!=================================================================================================!
!------------------------------------- OPENNING AND WRITING: -------------------------------------!
!=================================================================================================!


    ! Open file, get size of time dimension, and put current time
    ! -----------------------------------------------------------
    call open_file(ofile_name, fid, mode=NF90_WRITE)
    call inquire_dim(fid, time_dimname, dimid)
    ierr = nf90_inquire_dimension(fid, dimid, len=nt)
    nt = nt + 1
    call put_var(fid, varname=time_dimname, var_real0D=real(t), stt=(/nt/), cnt=(/1/))



    !<><><><><><><><><><><><><><><><><>!
    !<> %%%%%%%%%%%%%%%%%%%%%%%%%%%% <>!
    !<> %  write output variables  % <>!
    !<> %%%%%%%%%%%%%%%%%%%%%%%%%%%% <>!
    !<><><><><><><><><><><><><><><><><>!

    ! Note: the following list of blocks needs to be updated if one wants to add new output variables
    ! ***********************************************************************************************

    do i = 4,8 ! => climatic parameters # 1 to 5
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real0D=real(cpvec(i-3)), stt=(/nt/), cnt=(/1/))
    end do
    !
    i = 9 ! Temperature
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(temp, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 10 ! Runoff
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(runoff, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 12 ! Silicate weathering (sum of all lihtology)
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(wth_allsil, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 13 ! Silicate weathering by lithology (fraction-weighted)
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
 var_real3D=real(reshape(wth_litho_wgh, shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 14 ! Silicate weathering by lithology (unweighted)
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
     var_real3D=real(reshape(wth_litho, shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 15 ! Kerogen weathering
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(ker_wth, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 16 ! Biogenic organic carbon export
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(bioC_exp, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 17 ! Phosphorus weathering
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(phos_wth, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 18 ! Pyrite weathering
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(pyr_wth, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))


    ! Close output file
    ! -----------------

    call close_file(fid)



end subroutine


end module
