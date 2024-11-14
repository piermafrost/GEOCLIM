module geoclim_write_output_mod
implicit none

contains

subroutine geoclim_write_output(t, COMB_outvar_info)
  use netcdf
  use io_module,        only: netcdf_output_var
  use netcdf_io_module, only: open_file, close_file, inquire_dim, put_var
  use water_chemistry,  only: omega_with_pressure
  use constante,        only: sink_veloc, rho_sed, PI_n_O2_atm, PI_n_CO2_atm, PI_n_rest_of_atm
  !
  double precision, parameter:: PI_atm_moles = PI_n_O2_atm+PI_n_CO2_atm+PI_n_rest_of_atm ! molecular mass of the atmosphere
  !
  include 'combine_foam.inc'
  integer:: ierr, fid, timevarid, dimid, nt
  double precision, dimension(nbasin):: omega, fcarb_dep, diag_dplysc, diag_dplysa
  ! Declare structures for netCDF output info (cannot be declared in combine_foam.inc)
  type(netcdf_output_var), dimension(nCOMBoutvar), intent(in) :: COMB_outvar_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!           offline computations:           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xPOPexport=0.
    do j0=1,nsurface
        j = jbox_surface(j0)
        xPOPexport=xPOPexport + sink_veloc*var_part(2,j)*oce_surf(j)
        !totalCO2_exch=totalCO2_exch+fco2atm_ocean(j)
    enddo

    vol_tot=0.
    vartot_diss(:)=0.
    vartot_part(:)=0.
    vartot_isot(:)=0.
    ph_tot=0.
    temp_tot=0.
    salin_tot=0.
    freef_tot=0.
    fodc_tot=0.
    fcarbdep_tot=0.
    fco2atm_ocean_tot=0.
    F_seafloor_cdiss_tot=0.
    fSulfRed_tot=0.
    fO2_odc_tot=0.
    do j=1,nbasin-1
        vol_tot              = vol_tot + box_vol(j)
        vartot_diss(:)       = vartot_diss(:) + var_diss(:,j)*box_vol(j) ! dissolved variables
        vartot_part(:)       = vartot_part(:) + var_part(:,j)*box_vol(j) ! particulate variables
        vartot_isot(1)       = vartot_isot(1) + var_isot(1,j)*var_diss(1,j)*box_vol(j) ! DIC d13C
        vartot_isot(2)       = vartot_isot(2) + var_isot(2,j)*var_part(5,j)*box_vol(j) ! PIC d13C
        vartot_isot(3)       = vartot_isot(3) + var_isot(3,j)*var_part(4,j)*box_vol(j) ! POC d13C
        vartot_isot(5)       = vartot_isot(5) + var_isot(5,j)*var_diss(5,j)*box_vol(j) ! diss Sr isotopic ration
        vartot_isot(6)       = vartot_isot(6) + var_isot(6,j)*var_part(1,j)*box_vol(j) ! PIC Sr isotopic ration
        !vartot_isot()       = vartot_isot() + var_isot(,j)*var_diss(,j)*box_vol(j) ! d7Li
        ph_tot               = ph_tot + pH(j)*box_vol(j)
        temp_tot             = temp_tot + temp_box(j)*box_vol(j)
        salin_tot            = salin_tot + salin(j)*box_vol(j)
        freef_tot            = freef_tot+freef(j)
        fodc_tot             = fodc_tot+fodc(j)
        fco2atm_ocean_tot    = fco2atm_ocean_tot + fco2atm_ocean(j)
        F_seafloor_cdiss_tot = F_seafloor_cdiss_tot + F_seafloor_cdiss(j)
        fSulfRed_tot         = fSulfRed_tot + fSulfRed(j)
        fO2_odc_tot          = fO2_odc_tot + fO2_odc(j)
    end do

    ! Carbonate deposition on sedimenting oceanic basins: all "raining" PIC is buried
    fcarb_dep = frain_PIC
    do j0=1,nsedi
        j = jbox_sedi(j0)
        fcarbdep_tot = fcarbdep_tot + fcarb_dep(j)
    end do

    ! Normalize isotopic variables
    vartot_isot(1) = vartot_isot(1) / vartot_diss(1) ! DIC d13C
    !
    if (vartot_part(5) == 0d0) then
        vartot_isot(2) = COMB_outvar_info(14)%fillval
    else
        vartot_isot(2) = vartot_isot(2) / vartot_diss(5) ! PIC d13C
    end if
    !
    vartot_isot(3) = vartot_isot(3) / vartot_part(4) ! POC d13C
    !
    if (vartot_diss(5) == 0d0) then
        vartot_isot(5) = COMB_outvar_info(17)%fillval
    else
        vartot_isot(5) = vartot_isot(5) / vartot_diss(5) ! diss Sr iso ration
    end if
    !
    if (vartot_part(1) == 0d0) then
        vartot_isot(6) = COMB_outvar_info(18)%fillval
    else
        vartot_isot(6) = vartot_isot(6) / vartot_part(1) ! PIC Sr iso ration
    end if
    !
    !vartot_isot()    = vartot_isot() / vartot_diss() ! d7Li

    ! Normalize mean intensive variables (e.g., concentration)
    vartot_diss    = vartot_diss / vol_tot ! dissolved variables
    vartot_part    = vartot_part / vol_tot ! particulate variables

    ph_tot         = ph_tot / vol_tot
    temp_tot       = temp_tot / vol_tot
    salin_tot      = salin_tot / vol_tot

    ! Partial pressure of O2 and CO2 (in units of PI atmosphere)
    pO2_atm  = var_diss(6,nbasin) / PI_atm_moles
    pCO2_atm = 1e6 * var_diss(7,nbasin) / PI_atm_moles

    ! convert mol to volume fraction (= molar fraction) for 02 and CO2
    O2_atm_conc   =  ( var_diss(6,nbasin) / ( PI_n_rest_of_atm + var_diss(6,nbasin) + var_diss(7,nbasin) ) )*100 ! (%)
    CO2_atm_conc  =  ( var_diss(7,nbasin) / ( PI_n_rest_of_atm + var_diss(6,nbasin) + var_diss(7,nbasin) ) )*1e6 ! (ppm)

    ! calcite saturation at actual boxes pressure
    do j0 = 1,nsurface
        j = jbox_surface(j0)
        omega(j) = omega_0(j)
    end do
    do j0 = 1,nnosurface
        j = jbox_nosurface(j0)
        omega(j) = omega_with_pressure(omega_0(j), press_box(j), temp_box(j), aragonite=.false.)
    end do
    omega(nbasin) = COMB_outvar_info(31)%fillval

    ! Interpretable lysocline depth (only if outputs are written)
    if (COMB_outvar_info(50)%writevar .or. COMB_outvar_info(74)%writevar) then
       !                 ^^ diag calcite lysocline depth    ^^ area-averaged diag calcite lysocline depth
        call diagnose_lysocline_depth(nsurface, jbox_surface, depth_box, depth_top_box, dplysc, diag_dplysc, &
                                      fillval=COMB_outvar_info(50)%fillval)
    end if
    if (COMB_outvar_info(52)%writevar) then
       !                 ^^ diag aragonite lysocline depth
        call diagnose_lysocline_depth(nsurface, jbox_surface, depth_box, depth_top_box, dplysa, diag_dplysa, &
                                      fillval=COMB_outvar_info(52)%fillval)
    end if



!=================================================================================================!
!------------------------------------- OPENNING AND WRITING: -------------------------------------!
!=================================================================================================!


    ! Open file, get size of time dimension, and put current time
    ! -----------------------------------------------------------
    call open_file(COMB_ofile_name, fid, mode=NF90_WRITE)
    call inquire_dim(fid, COMB_time_dimname, dimid)
    ierr = nf90_inquire_dimension(fid, dimid, len=nt)
    nt = nt + 1
    call put_var(fid, varname=COMB_time_dimname, var_real0D=real(t), stt=(/nt/), cnt=(/1/))



    !<><><><><><><><><><><><><><><><><>!
    !<> %%%%%%%%%%%%%%%%%%%%%%%%%%%% <>!
    !<> %  write output variables  % <>!
    !<> %%%%%%%%%%%%%%%%%%%%%%%%%%%% <>!
    !<><><><><><><><><><><><><><><><><>!

    ! Note: the following list of blocks needs to be updated if one wants to add new output variables
    ! ***********************************************************************************************

    do i = 14,18 ! => climatic parameter # 1 to 5
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(cpvec(i-13)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    i = 19 ! Water exchange
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real2D=real(F(:,:)/31.536d12), &
                     stt=(/1,1,nt/), cnt=(/nbasin,nbasin,1/))
    !
    do i = 20,27 ! => 8 main COMBINE dissolved variables
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(var_diss(i-19,:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    end do
    !
    do i = 28,32 ! => 5 main COMBINE particulate variables
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(var_part(i-27,:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    end do
    !
    do i = 33,38 ! => 6 main COMBINE isotopic variables
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(var_isot(i-32,:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    end do
    !
    i = 39
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(h2co3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 40
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(hco3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 41
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(co3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 42
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dh2co3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 43
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dhco3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 44
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dco3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 45
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(ph(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 46
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(omega(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 47
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(temp_box(1:nbasin-1)-273.15), &
                     stt=(/1,nt/), cnt=(/nbasin-1,1/))
    !
    i = 48
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(salin(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 49 ! Calcite lysocline depth
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dplysc(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 50
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(diag_dplysc(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 51 ! Aragonite lysocline depth
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dplysa(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 52
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(diag_dplysa(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 53 ! org C fractionation parameter EPSILON
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(epsiC(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    do i = 54,59 ! => box-average values of the 6 first main COMBINE dissolved variables
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_diss(i-53)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    i = 60 ! => box-average value of the main COMBINE dissolved variable #8 (SO4^2-)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_diss(8)), &
                     stt=(/nt/), cnt=(/1/))
    !
    do i = 61,65 ! => box-average values of the 5 main COMBINE particulate variables
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_part(i-60)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    do i = 66,68 ! => box-average values of the main COMBINE isotopic variables #1 to #3
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_isot(i-65)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    do i = 69,70 ! => box-average values of the main COMBINE isotopic variables #5 to #6
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_isot(i-64)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    i = 71
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(ph_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 72
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(temp_tot-273.15), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 73
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(salin_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 74 ! Area-averaged calcite lysocline depth
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real( &
                         area_average_lysdep(nsurfnoappcont, jbox_surfnoappcont, oce_surf, diag_dplysc)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 75
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(pO2_atm), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 76
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(O2_atm_conc), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 77
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(pCO2_atm), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 78
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(CO2_atm_conc), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 79 ! Global Mean Surface Temperature
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(temp_box(nbasin)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 80
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(xPOPexport), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 81 ! Anthropogenic CO2 degassing
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fanthros), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 82 ! Trapp CO2 degassing
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(ftrap), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 83 ! MOR CO2 degassing
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fmor(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 84 ! Atmosphere -> Ocean CO2 flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fco2atm_ocean(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 85 ! Net carbonate bioproduction
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(finorgC(:)) ,&
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 86 ! Net organic bioproduction
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fbioC(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 87 ! PIC sink flux: sink_veloc * [PIC] * horizontal area
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(sink_veloc*var_part(5,:)*oce_surf(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 88 ! POC sink flux: sink_veloc * [POC] * horizontal area
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(sink_veloc*var_part(4,:)*oce_surf(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 89 ! PIC re-dissolution flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fdissol_carb(:)*var_part(5,:)*box_vol(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 90 ! POC remineralization flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(roxyd(:)*var_part(4,:)*box_vol(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 91 ! Bulk particles rain flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(rho_sed*frain_sed(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 92 ! Dowslope sediment export flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real2D=real(rho_sed*fexch_sed(:,:)), &
                     stt=(/1,1,nt/), cnt=(/nbasin,nbasin,1/))
    !
    i = 93 ! Sedimentation flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(rho_sed*fin_sed(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 94 ! Sedimentation rate
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(ws(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 95 ! PIC rain flux on seafloor: sink_veloc * [PIC] * sedimenting area
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(frain_PIC(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 96 ! POC rain flux on seafloor: sink_veloc * [POC] * sedimenting area
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(frain_POC(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 97 ! Neritic carbonate deposition flux before redissolution in sediment
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(freef(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 98 ! Pelagic carbonate deposition flux before redissolution in sediment
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fcarb_dep(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 99 ! Organic C deposition flux before remineralization in sediment
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fin_POC(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 100 ! Neritic carbonate burial flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(freef(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 101 ! Pelagic carbonate burial flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fcarb_dep(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 102 ! Organic carbon burial flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fodc(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 103 ! Burial efficiency
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(Corg_BE(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 104 ! P burial flux orgC-bound
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fodp(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 105 ! P burial flux in form of phosphorite
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fphos(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 106 ! P burial flux hydrothermal Fe-bound
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fhydP(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 107 ! Sulfate-reduction
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fSulfRed(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 108 ! O2 flux due to org C burial
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fO2_odc(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 109 ! Seafloor weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(F_seafloor_cdiss(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 110 ! Total atmosphere -> Ocean CO2 flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fco2atm_ocean_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 111 ! Carbonate neritic total burial flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(freef_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 112 ! Carbonate pelagic total burial flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fcarbdep_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 113 ! Organic carbon total burial flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fodc_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 114 ! Phosphorus total burial flux (on all forms)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(fodp+fphos+fhydP)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 115 ! Total sulfate-reduction flux
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fSulfRed_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 116 ! Total O2 flux due to organic C burial
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fO2_odc_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 117 ! Total seafloor weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(F_seafloor_cdiss_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 118 ! Total CO2 degassing (all forms)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fvol+sum(fmor)+ftrap+fanthros), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 119 ! Freshwater discharge
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(discharge)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 120 ! Sediment discharge
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(tss)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 121 ! Silicate weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(fsilw+fbasw)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 122 ! Basalt weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(fbasw)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 123 ! Carbonate weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(fcarbw)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 124 ! Kerogen weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(fkerw)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 125 ! Pyrite weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(fcarbsulfw+fsilsulfw+fH2SO4sulfw)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 126 ! Silicate sulfuric weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(fsilsulfw)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 127 ! Carbonate sulfuric weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(fcarbsulfw)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 128 ! Phosphorus weathering
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(fpw)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 129 ! Continental biogenic C export
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(total_cont_POC_export)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 130 ! Freshwater discharge in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(discharge(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 131 ! Sediment discharge in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(tss(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 132 ! Silicate weathering in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fsilw(:)+fbasw(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 133 ! Basalt weathering in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fbasw(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 134 ! Carbonate weathering in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fcarbw(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 135 ! Kerogen weathering in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fkerw(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 136 ! Pyrite weathering in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fcarbsulfw(:)+fsilsulfw(:)+fH2SO4sulfw(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 137 ! Silicate sulfuric weathering in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fsilsulfw(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 138 ! Carbonate sulfuric weathering in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fcarbsulfw(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 139 ! Phosphorus weathering in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fpw(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 140 ! Continental biogenic C export in each oceanic box
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(total_cont_POC_export(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 141
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fCO2_crust(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 142
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fSO4_basin(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 143
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fSO4_crust(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 144
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(dLiriv*FrivLi)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 145
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(sum(FrivLi)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 146
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dLiriv), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 147
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(FrivLi), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 148
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(epsiCont), &
                     stt=(/nt/), cnt=(/1/))


  ! Close output file
  ! -----------------

  call close_file(fid)



end subroutine


! ============================== !


subroutine diagnose_lysocline_depth(nsurface, jbox_surface, box_depth, box_top_depth, lysdep, diag_lysdep, fillval)
    include 'shape.inc'
    integer, intent(in):: nsurface, jbox_surface(nbasin)
    double precision, dimension(nbasin), intent(in):: box_depth, box_top_depth, lysdep
    double precision, dimension(nbasin), intent(out):: diag_lysdep
    double precision, intent(in), optional:: fillval
    integer:: j, j0
    double precision:: h
    logical:: go_on

    do j0 = 1,nsurface

        go_on = .true.
        do j = jbox_surface(j0), jbox_surface(j0+1)-1
            if (go_on) then
                if (lysdep(j) <= box_depth(j)) then
                    go_on = .false.
                    h = max(lysdep(j), box_top_depth(j))
                end if
            end if
        end do
        if (go_on) h = lysdep(j-1) ! = last value of the loop

        diag_lysdep(jbox_surface(j0) : jbox_surface(j0+1)-1) = h

    end do

    ! fillvalue on atmospheric box
    if (present(fillval)) diag_lysdep(nbasin) = fillval

end subroutine


! ----- !


function area_average_lysdep(nsurfnoappcont, jbox_surfnoappcont, oce_surf, lysdep)
    include 'shape.inc'
    integer, intent(in):: nsurfnoappcont, jbox_surfnoappcont(nbasin)
    double precision, dimension(nbasin), intent(in):: oce_surf, lysdep
    double precision:: area_average_lysdep, tot_ocesurf
    integer:: j0, j
    !
    area_average_lysdep = 0
    tot_ocesurf = 0
    !
    do j0 = 1,nsurfnoappcont
        ! => exclude epicontinental boxes for the computation of average lysocline depth
        !    because the estimation of lysocline in those boxes is not realistic
        j = jbox_surfnoappcont(j0)
        area_average_lysdep = area_average_lysdep + lysdep(j)*oce_surf(j)
        tot_ocesurf = tot_ocesurf + oce_surf(j)
    end do
    !
    area_average_lysdep = area_average_lysdep/tot_ocesurf
    !
end function



end module

