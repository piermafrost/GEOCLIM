module output_netcdf2ascii_mod
implicit none
contains


subroutine output_netcdf2ascii(COMB_outvar_info)
use netcdf
use netcdf_io_module, only: nf90_check
use io_module, only: netcdf_output_var, DEFAULT_FILLVAL
use constante, only: PI_n_CO2_atm, PI_n_O2_atm, PI_n_rest_of_atm

  double precision, parameter:: PI_atm_moles = PI_n_O2_atm+PI_n_CO2_atm+PI_n_rest_of_atm ! molecular mass of the atmosphere
! ********************************************************************************
  include 'combine_foam.inc'
! ********************************************************************************
  character(len=200) :: dummystr
  integer, dimension(nCOMBoutvar):: varid
  integer:: fid, Bfid
  integer:: ierr, ntime, Bntime, tfid, tvarid, Btvarid
  double precision:: vector(1), PIC_snk_flx(nbasin)
  logical:: open_COMB_file, got_var(nCOMBoutvar)
  ! Declare structures for netCDF output info (cannot be declared in combine_foam.inc)
  type(netcdf_output_var), dimension(nCOMBoutvar), intent(in) :: COMB_outvar_info


  print *
  print *
  print *, '***********************************************'
  print *, '* conversion of netCDF output in ASCII format *'
  print *, '***********************************************'
  print *


  ! create ascii output:
  !----------------------

  open (7,file=trim(output_directory)//'var'//run_name,status='REPLACE')
  open (11,file=trim(output_directory)//'box1'//run_name,status='REPLACE')
  open (12,file=trim(output_directory)//'box2'//run_name,status='REPLACE')
  open (13,file=trim(output_directory)//'box3'//run_name,status='REPLACE')
  open (14,file=trim(output_directory)//'box4'//run_name,status='REPLACE')
  open (15,file=trim(output_directory)//'box5'//run_name,status='REPLACE')
  open (16,file=trim(output_directory)//'box6'//run_name,status='REPLACE')
  open (17,file=trim(output_directory)//'box7'//run_name,status='REPLACE')
  open (18,file=trim(output_directory)//'box8'//run_name,status='REPLACE')
  open (19,file=trim(output_directory)//'box9'//run_name,status='REPLACE')
  open (20,file=trim(output_directory)//'box10'//run_name,status='REPLACE')
  open (21,file=trim(output_directory)//'cflux'//run_name,status='REPLACE')
  open (22,file=trim(output_directory)//'chimie'//run_name,status='REPLACE')
  open (23,file=trim(output_directory)//'lysocline'//run_name,status='REPLACE')
  open (24,file=trim(output_directory)//'trap1'//run_name,status='REPLACE')
  open (25,file=trim(output_directory)//'trap2'//run_name,status='REPLACE')
  open (26,file=trim(output_directory)//'seacarb_diss'//run_name,status='REPLACE')
  open (27,file=trim(output_directory)//'speciation'//run_name,status='REPLACE')
  open (28,file=trim(output_directory)//'speciation_c13'//run_name,status='REPLACE')
  open (29,file=trim(output_directory)//'forcing'//run_name,status='REPLACE')
  open (205,file=trim(output_directory)//'lithium'//run_name,status='REPLACE')



  ! open geoclim netcdf output:
  !----------------------------

  ierr = nf90_open(COMB_ofile_name, NF90_NOWRITE, fid)
  call nf90_check(ierr, 'Error while openning file "'//trim(COMB_ofile_name)// &
                  '". Cannot convert COMBINE netCDF outputs in ASCII format', kill=.false.)

  if (ierr==NF90_NOERR) then
    open_COMB_file = .true.
    ! get time length:
    ierr = nf90_inq_dimid(fid, COMB_time_dimname, tvarid)
    call nf90_check(ierr, 'Error while getting ID of dimension "'//trim(COMB_time_dimname)//'"')
    ierr = nf90_inquire_dimension(fid, tvarid, len=ntime)
    call nf90_check(ierr, 'Error while inquiring length of dimension "'//trim(COMB_time_dimname)//'"')
    ierr = nf90_inq_varid(fid, COMB_time_dimname, tvarid)
    call nf90_check(ierr, 'Error while getting ID of variable "'//trim(COMB_time_dimname)//'"')
    ! get variables identifier:
    do k = 14,nCOMBoutvar ! Skip first 13 variables (time-invariant fields)
      if (COMB_outvar_info(k)%writevar) then
        ierr = nf90_inq_varid(fid, COMB_outvar_info(k)%vname, varid(k))
        call nf90_check(ierr, 'Error while getting ID of variable "'//trim(COMB_outvar_info(k)%vname)//'"', kill=.false.)
        got_var(k) = .true.
      else
        got_var(k) = .false. ! indicate that variable was not saved
      end if
    end do
  else
    open_COMB_file = .false. ! indicate that netCDF file couldn't be open
  end if



  ! ========================================================================== !

  ! Fill-value on all variables (in case they were not saved)
  t                    = DEFAULT_FILLVAL
  var_diss             = DEFAULT_FILLVAL
  var_part             = DEFAULT_FILLVAL
  var_isot             = DEFAULT_FILLVAL
  h2co3(:)             = DEFAULT_FILLVAL
  hco3(:)              = DEFAULT_FILLVAL
  co3(:)               = DEFAULT_FILLVAL
  dh2co3(:)            = DEFAULT_FILLVAL
  dhco3(:)             = DEFAULT_FILLVAL
  dco3(:)              = DEFAULT_FILLVAL
  ph(:)                = DEFAULT_FILLVAL
  omega_0(:)           = DEFAULT_FILLVAL
  temp_box(:)          = DEFAULT_FILLVAL
  salin(:)             = DEFAULT_FILLVAL
  dplysc(:)            = DEFAULT_FILLVAL
  dplysa(:)            = DEFAULT_FILLVAL
  vartot_diss          = DEFAULT_FILLVAL
  vartot_part          = DEFAULT_FILLVAL
  vartot_isot          = DEFAULT_FILLVAL
  ph_tot               = DEFAULT_FILLVAL
  temp_tot             = DEFAULT_FILLVAL
  salin_tot            = DEFAULT_FILLVAL
  pO2_atm              = DEFAULT_FILLVAL
  pCO2_atm             = DEFAULT_FILLVAL
  xPOPexport           = DEFAULT_FILLVAL
  fanthros             = DEFAULT_FILLVAL
  fco2atm_ocean(:)     = DEFAULT_FILLVAL
  fco2atm_ocean_tot    = DEFAULT_FILLVAL
  finorgC(:)           = DEFAULT_FILLVAL
  fdissol_carb(:)      = DEFAULT_FILLVAL
  fsilw(1)             = DEFAULT_FILLVAL
  fbasw(1)             = DEFAULT_FILLVAL
  fkerw(1)             = DEFAULT_FILLVAL
  freef(:)             = DEFAULT_FILLVAL
  freef_tot            = DEFAULT_FILLVAL
  fcarbdep_tot         = DEFAULT_FILLVAL
  fodc(:)              = DEFAULT_FILLVAL
  fodc_tot             = DEFAULT_FILLVAL
  PIC_snk_flx(:)       = DEFAULT_FILLVAL
  fpw(1)               = DEFAULT_FILLVAL
  fbioC(:)             = DEFAULT_FILLVAL
  F_seafloor_cdiss(:)  = DEFAULT_FILLVAL
  F_seafloor_cdiss_tot = DEFAULT_FILLVAL
  ftrap                = DEFAULT_FILLVAL
  fCO2_crust(:)        = DEFAULT_FILLVAL
  fSO4_basin(:)        = DEFAULT_FILLVAL
  fSO4_crust(:)        = DEFAULT_FILLVAL
  FrivLi(1)            = DEFAULT_FILLVAL
  dLiriv(1)            = DEFAULT_FILLVAL


  do k = 1,ntime

    ! read geoclim netcdf output:
    !----------------------------

    if (open_COMB_file) then
      ierr = nf90_get_var( fid, tvarid , vector , start=(/k/), count=(/1/) )
      t = vector(1)
      !
      do i = 20,27
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  var_diss(i-19,:)     , start=(/1,k/), count=(/nbasin,1/)  )
      end do
      !
      do i = 28,32
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  var_part(i-27,:)     , start=(/1,k/), count=(/nbasin,1/)  )
      end do
      !
      do i = 33,38
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  var_isot(i-32,:)     , start=(/1,k/), count=(/nbasin,1/)  )
      end do
      !
      i = 39
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  h2co3(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 40
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  hco3(:)              , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 41
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  co3(:)               , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 42
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dh2co3(:)            , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 43
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dhco3(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 44
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dco3(:)              , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 45
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  ph(:)                , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 46
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  omega_0(:)           , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 47
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  temp_box(:)          , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 48
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  salin(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 49
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dplysc(:)            , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 51
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dplysa(:)            , start=(/1,k/), count=(/nbasin,1/)  )
      !
      do i = 54,59
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vartot_diss(i-53:i-53) , start=(/k/), count=(/1/)          )
      end do
      !
      i = 60
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vartot_diss(8:8) ,      start=(/k/), count=(/1/)  )
      !
      do i = 61,65
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vartot_part(i-60:i-60) , start=(/k/), count=(/1/)          )
      end do
      !
      do i = 66,68
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vartot_isot(i-65:i-65) , start=(/k/), count=(/1/)          )
      end do
      !
      do i = 69,70
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vartot_isot(i-64:i-64) , start=(/k/), count=(/1/)          )
      end do
      !
      i = 71
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         ph_tot = vector(1)
      !
      i = 72
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         temp_tot = vector(1)
      !
      i = 73
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         salin_tot = vector(1)
      !
      i = 75
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         pO2_atm = vector(1)
      !
      i = 77
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         pCO2_atm = vector(1)
      !
      i = 80
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         xPOPexport = vector(1)
      !
      i = 81
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fanthros = vector(1)
      !
      i = 84
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fco2atm_ocean(:)     , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 110
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fco2atm_ocean_tot = vector(1)
      !
      i = 85
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  finorgC(:)           , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 89
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fdissol_carb(:)      , start=(/1,k/), count=(/nbasin,1/) )
      !
      i = 121
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fsilw(1) = vector(1)
      !
      i = 122
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fbasw(1) = vector(1)
      !
      i = 123
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fcarbw(1) = vector(1)
      !
      i = 124
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fkerw(1) = vector(1)
      !
      i = 100
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  freef(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 111
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         freef_tot = vector(1)
      !
      i = 112
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fcarbdep_tot = vector(1)
      !
      i = 102
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fodc(:)              , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 113
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fodc_tot = vector(1)
      !
      i = 87
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  PIC_snk_flx(:)       , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 128
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fpw(1) = vector(1)
      !
      i = 86
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fbioC(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 109
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  F_seafloor_cdiss(:)  , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 117
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         F_seafloor_cdiss_tot = vector(1)
      !
      i = 82
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         ftrap = vector(1)
      !
      i = 141
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fCO2_crust(:)        , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 142
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fSO4_basin(:)        , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 143
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fSO4_crust(:)        , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 145
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         FrivLi(1) = vector(1)
      !
      i = 144
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         dLiriv(1) = vector(1)
      i = 129
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         total_cont_POC_export(1) = vector(1)
    end if



    ! write ascii output:
    !--------------------

    write(7,9) t, pO2_atm*PI_atm_moles/PI_n_O2_atm, 1e-6*pCO2_atm*PI_atm_moles/PI_n_CO2_atm, &
                  fanthros,fco2atm_ocean(1),fco2atm_ocean(3), &
                  fco2atm_ocean(6),fco2atm_ocean(8),xPOPexport, &
                  finorgC(3),fdissol_carb(3),dplysa(3)

    do j = 1,10
      write(10+j,11) t, (var_diss(i,j),i=1,5), (var_part(i,j),i=1,5), (var_diss(i,j),i=6,7), (var_isot(i,j),i=1,6), var_diss(8,j)
    end do

    write(21,12)t,fsilw(1),fbasw(1),fcarbw(1),fkerw(1),freef_tot, &
                    fcarbdep_tot,fodc_tot, &
                    (freef(j),j=1,9), &
                    (fodc(j),j=1,9), &
             (PIC_snk_flx(j),j=1,9),fpw(1), &
             (fbioC(j),j=1,9),total_cont_POC_export(1)

    write(205,12)t,FrivLi(1),dLiriv(1)

    write(23,10)t,(dplysc(i),i=1,nbasin-1)

    write(22,15)t,vartot_diss(1),vartot_diss(2),ph_tot, &
                    (pH(j),j=1,9),(omega_0(j),j=1,9)

    write(24,12)t,ftrap,(fCO2_crust(j),j=1,9),(fSO4_basin(j),j=1,9),(fSO4_crust(j),j=1,9)
    write(26,15)t,(F_seafloor_cdiss(j),j=1,nbasin-1)

    write(27,12)t,(h2co3(j),j=1,nbasin-1),(hco3(j),j=1,nbasin-1),(co3(j),j=1,nbasin-1)
    write(28,12)t,(dh2co3(j),j=1,nbasin-1),(dhco3(j),j=1,nbasin-1),(dco3(j),j=1,nbasin-1)

    write(29,17)t/1e6,var_diss(7,nbasin)/PI_n_CO2_atm,temp_box(1),temp_box(3), &
                  temp_box(6),temp_box(8),salin(3),salin(4),salin(6),dco3(3)*1d3,dco3(6)*1d3


   9    format(1f15.5,1x,2(f15.6,1x),20(es15.6e3,1x))
   11   format(21(es15.6e3,1x))
   12   format(500(es15.6e3,1x))
   10   format(10(es15.6e3,1x))
   15   format(23(es15.6e3,1x))
   16   format(1921(es15.6e3,1x))
   17   format(30(f10.4,1x))
   18   format(1920(es15.6e2,1x))
   19   format(1f15.5,1x,4i5)



  end do



  ! close netcdf output:
  !---------------------

  if (open_COMB_file) ierr = nf90_close(fid)



  ! close ascii output:
  !--------------------

  close(7)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)
  close(23)
  close(24)
  close(25)
  close(26)
  close(27)
  close(28)
  close(29)
  close(205)



end subroutine


end module
