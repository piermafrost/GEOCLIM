!   Technical parameters and variables:
!   ==================================================================

    include 'path.inc' ! => get the variable 'geoclim_path'
    character(len=100)    :: run_name, check_run_name, killing_file
    character(len=500)    :: output_directory
    character(len=30)     :: COMB_init_mode, cont_input_mode
    integer, dimension(5) :: ERROR_HANDLING_OPTION

    include 'coupler.inc' ! coupling parameters: 'coupling_dynsoil' and 'coupling_ecogeo' plus 'CO2_interpolation' 

    include 'shape.inc'
    ! => define: number of boxes ('nbasin')
    !            number of CO2 levels ('nclimber')
    !            geographic dimensions ('nlon' and 'nlat') -> number of pixels
    !            number of lithology classes ('nlitho'),
    !            number of DynSoil levels ('nDSlev') 
    !            number of additional climatic parameters ('len_p1', 'len_p2', 'len_p3', 'len_p4', and 'len_p5')
    integer, parameter :: nvar_diss=8, nvar_part=5, nvar_isot=6
    integer, parameter :: nvb = (nvar_diss + nvar_part + nvar_isot)*nbasin
    integer, parameter :: npixel=nlon*nlat

    integer, parameter :: ndegas=1000000,nsres=24
    include 'output_size.inc' ! => define nCOMBoutvar, nGEOGoutvar and nDYNSoutvar
    integer, parameter :: n_peaks = 5





!   Variables:
!   ==================================================================
    logical decide,check
    integer i,j,j0,k,k0,kl,icount,n,ishelfal,k1,k2,ncount,kbasin
    integer clock
    integer ipeak
    logical:: convert2ascii
    integer ijump_climate,icount_climate, ijump_geogprint,icount_geogprint
    integer, dimension(:), allocatable :: seed
    integer, dimension(8):: computer_time

    double precision tstart_deg2,tend_deg2
    double precision out,O2,shelfal,xkill,sluggish
    double precision DOA(1:nbasin),cp_burial(1:nbasin),t
    double precision fCO2atm_ocean(1:nbasin),po2,fC13atm_ocean(1:nbasin), &
                     fC13ocean_atm(1:nbasin)
    double precision oce_surf(nbasin), oce_surf_tot, temp_box(nbasin), depth_box(nbasin), depth_top_box(nbasin), box_thick(nbasin),&
                     press_box(nbasin), surf_sedi(nbasin), sedbastype_fract(nbasin),accumul_capacity(nbasin)
    integer, dimension(nbasin) :: indice_deep, indice_sedi, app_cont,indice_thermo, indice_surface, indice_epicont, indice_polar
    double precision, dimension(nbasin) :: appcont_fract
    integer :: ndeep, nnodeep, nsedi, nnosedi, nthermo, nnothermo, nsurface, nnosurface, &
               nepicont, nnoepicont, npolar, nnopolar, nappcont, nnoappcont, nsediepicsurf, &
               nsediepicnosurf, nsediintermed, nsedideep, nsurfnoappcont, nnosurfnoappcont
    integer, dimension(nbasin) :: jbox_deep, jbox_nodeep, jbox_sedi, jbox_nosedi, jbox_thermo, jbox_nothermo, &
                                  jbox_surface, jbox_nosurface, jbox_epicont, jbox_noepicont, jbox_polar, jbox_nopolar, &
                                  jbox_appcont, jbox_noappcont, jbox_sediepicsurf, jbox_sediepicnosurf, &
                                  jbox_sediintermed, jbox_sedideep, jbox_surfnoappcont, jbox_nosurfnoappcont
    double precision carb_ratio(1:nbasin),fbioP(1:nbasin),fbioC(1:nbasin),finorgC(1:nbasin), &
                     finorgP(1:nbasin)
    double precision fdownt,fdowntP
    double precision dc13PTAS,dc13PTAD,dc13PTACS,dc13PTACD,dc13PTAT &
                     ,dco3(1:nbasin),dhco3(1:nbasin),dh2co3(1:nbasin) &
                     ,dc13atm,dctrap
    double precision epsiC(1:nbasin)
    double precision epsiCont
    double precision Toceclimber(nbasin, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5)
    double precision Tairclimber(npixel, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5)
    double precision Runclimber(npixel, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5)
    double precision GMSTclimber(nclimber, len_p1, len_p2, len_p3, len_p4, len_p5)
    double precision areaclimber(1:npixel),co2climber(1:nclimber),areatot
    double precision areaEarth(1:npixel),area_ocean(1:nbasin)
    double precision hdeepC
    double precision fvol,fmor(1:nbasin),ftrap,fSO4_atmos,fSO4_ocean,fSO4_deep,fSO4_basin(1:nbasin) &
                     ,fSO4_crust(1:nbasin),fCO2_atmos,fCO2_crust(1:nbasin),fCO2_deep
    double precision clo,phosss
    double precision var_diss(nvar_diss,nbasin), var_part(nvar_part,nbasin), var_isot(nvar_isot,nbasin)
    double precision volin,xMORin,xlumin,berner,fvegin,oblin,shells &
                     ,xnoorg,xdrift,xsolar,temcondi,async,ageYprint &
                     ,oxy_acc_fact,sulf_acc_fact
    character(len=30) :: solver
    double precision fanthros
    double precision areaPTAS,areaPTACS,dcdiss,aconv,fevapnow &
                     ,fog,ak1soil,ak2soil,akspsoil,betasoil
    double precision freef(1:nbasin),fodc(1:nbasin),fodp(1:nbasin) &
                     ,freefP(1:nbasin),f_diss_c(1:nbasin),f_diss_a(1:nbasin) &
                     ,fcarbP_dep(1:nbasin) &
                     ,fdissol_carb(1:nbasin),fdissol_carbP(1:nbasin) &
                     ,ws(nbasin),frain_sed(nbasin),frain_POC(nbasin),frain_PIC(nbasin),frain_POP(nbasin) &
                     ,fin_sed(nbasin),fin_POC(nbasin),fin_POP(nbasin), &
                     Xsed(nbasin,nbasin),fexch_sed(nbasin,nbasin),Corg_BE(nbasin)
    logical sed_routing_flag
    double precision F_seafloor_cdiss(1:nbasin)
    double precision tend  
    !double precision deephyp(13),shyp(13) ! obsolete
    double precision dplysc(1:nbasin),dplysa(1:nbasin)
    double precision value_keep1(1:npixel),value_keep2(1:npixel)
    double precision fphos(1:nbasin)
    double precision fhydP(1:nbasin)
    double precision pco2soil,Hsoil,Pmax
    double precision rC_Corg(1:nbasin),reff(1:nbasin)
    double precision roxyd(1:nbasin),ranoxyd(1:nbasin)
    double precision runclim(1:npixel), discharge(nbasin)
    double precision tss(nbasin)
    double precision Tclim(1:npixel)
    double precision salin(1:nbasin)
    double precision omega_0(1:nbasin),omega_ara_0(1:nbasin)
    double precision St,S
    double precision h2co3(1:nbasin),hco3(1:nbasin),co3(1:nbasin), &
                     ch(1:nbasin),ak1(1:nbasin),ak2(1:nbasin), &
                     akb(1:nbasin),akc_0(1:nbasin),bco2(1:nbasin),ph(1:nbasin) &
                     ,pco2_dissous(1:nbasin)
    double precision ts,tsta,tfi,ifog
    integer :: ijump_print
    double precision rSrdep(1:nbasin)
    double precision SO4_conc(1:npixel),SO4soil
    double precision fodt,ept
    double precision tdata_deg(1:ndegas),data_deg(1:ndegas)
    double precision fracveg
    double precision box_vol(1:nbasin), vol(1:nvb)
    double precision, dimension(nbasin) :: fsilw,fcarbw,fkerw,fpw,fbasw,fcarbsulfw,fsilsulfw,fH2SO4sulfw ! add carbonate and silicate weathering by sulphuric acid
    ! Sulfate Reduction:
    double precision:: fodc_noSR_m2, fSulfRed_m2, fO2_odc_m2
    double precision, dimension(nbasin):: fSulfRed, fO2_odc
    double precision total_cont_POC_export(nbasin), eros_galy_unit
    double precision F(1:nbasin,1:nbasin), Fxch_clim(nbasin**2, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5)
    logical:: variable_watxch
    double precision R_diss(nvar_diss,nbasin), R_part(nvar_part,nbasin), R_isot(nvar_isot,nbasin)
    double precision ATe(1:npixel),BTe(1:npixel)
    double precision ARu(1:npixel),BRu(1:npixel)
    double precision litho_frac(1:nlitho,1:npixel)
    logical uniform_lithology
    double precision weighted_rsw(nbasin)
    double precision ref_x_axis(nlon), ref_y_axis(nlat)

    ! Climatic parameters
    double precision:: clim_param_1(len_p1), clim_param_2(len_p2), clim_param_3(len_p3), clim_param_4(len_p4), clim_param_5(len_p5)
    double precision:: cpvec(5), prev_cpvec(5), next_cpvec(5)
    double precision:: cp_lowest_value(5)
    integer:: ijump_climparam, icount_climparam

    double precision y(nvb),y_0(nvb),dydx(nvb),yout(nvb)
    double precision Corg_hml(1:nbasin),Corg_hsr(1:nbasin), &
                     fodcm2
    double precision area_cal(1:nbasin),area_ara(1:nbasin)
    double precision fsilj,fbasj,fcarbj,fker(1:npixel)
    double precision wth_allsil(1:npixel),wth_litho_wgh(1:nlitho,1:npixel),wth_litho(1:nlitho,1:npixel),pyr_wth(1:npixel)
    double precision POC_export_rate(1:npixel)
    double precision fp(1:npixel),fcarbw_pla(1:npixel)
    double precision Runinterd(1:nclimber),Tinterd(1:nclimber)
    double precision edb(1:nbasin),ebc(1:nbasin)
    double precision beta(1:nbasin), diff_diss(nvar_diss,nbasin), diff_part(nvar_part,nbasin), diff_isot(nvar_isot,nbasin), &
                     bilan_water(1:nbasin)

    double precision runoff,Patm,temperature,temper
    double precision x1,x2,xacc,xn,fx,dfx
    double precision h2co3soil,hco3soil,co3soil,calcium,xmass_carbonate

    double precision p,totarea,pco2,fco2,fe
    double precision facgra,facbas,facker,carb_weath_conc
    double precision x,xnum,xden,somme,xlogK,dum,RT0
    double precision pco2pal,fco2v,fev,fco2nv,fenv,time,pmin,prmax,aco2
    double precision htot,areal,dummy,eps,stept,tbegin,htry!,compteur
    double precision salsoil,prsoil,akbsoil,akcsoil,aksp,xo,fsoil,dfsoil
    double precision pr,fph,xv,dfph
    double precision xPOPexport,totalCO2_exch,freef_tot,fodc_tot,fcarbdep_tot,fco2atm_ocean_tot,F_seafloor_cdiss_tot,&
                     fSulfRed_tot,fO2_odc_tot
    double precision vartot_diss(nvar_diss), vartot_part(nvar_part), vartot_isot(nvar_isot), ph_tot,vol_tot,temp_tot,salin_tot,&
                     pO2_atm,pCO2_atm,O2_atm_conc,CO2_atm_conc

    double precision tka
    double precision ftrap1a,ftrap1b,ftrap2,ftrap3,ftrap4

    ! Lithium cycle
    double precision, dimension(nbasin) :: FrivLi, dLiriv
    double precision, dimension(nbasin) :: FhydLi

    character(len=500) :: COMB_ofile_name, GEOG_ofile_name
    character(len=100) :: COMB_time_dimname, GEOG_time_dimname

    ! continental variables
    logical:: uniform_routing
    integer, dimension(npixel):: list_cont_pixel, cont_basin_map
    integer:: ncontpxl

    ! dynsoil module variables:
    integer :: ijump_DS_integration, ijump_DS_print, icount_DS_int, icount_DS_pri
    double precision :: DS_timestep, xlevs(nDSlev)
    double precision, dimension(nlitho,npixel) :: reg_thick,reg_x_surf,reg_tau_surf,reg_P_vol,reg_mean_age,reg_x_mean
    double precision, dimension(nlitho,npixel) :: reg_thick_0,reg_x_surf_0,reg_tau_surf_0,reg_P_vol_0
    double precision, dimension(nlitho,npixel) :: reg_prod,reg_eros,reg_P_diss,reg_P_eros,reg_x_surf_eros
    double precision, dimension(nlitho,npixel) :: reg_Li_Friv, reg_Li_Fsp, reg_Li_driv
    double precision, dimension(nDSlev,nlitho,npixel) :: reg_z_prof,reg_tau_prof
    double precision, dimension(nDSlev,nlitho,npixel) :: reg_z_prof_0,reg_tau_prof_0
    double precision, dimension(npixel) :: slope, reg_eros_lithmean
    logical, dimension(npixel) :: DS_missingpoints
    integer, dimension(nlitho,npixel) :: reg_ktop, reg_ktop_0
    character(len=500) :: DYNS_ofile_name
    character(len=100) :: DYNS_time_dimname, DYNS_restart_dim(4), DYNS_restart_var(5)
    character(len=200) :: DynSoil_restart_file


    !LIP degassing
    double precision, dimension(n_peaks) :: pulse_start,pulse_end
    double precision, dimension(n_peaks) :: dc13_peak,amount_peak
    double precision, dimension(n_peaks) :: peak_duration

    ! vegetation variables:
    logical :: coupling_veget
    integer :: icount_veget, ijump_veget, veget_occup(npixel)
    double precision :: veget_factor(npixel), veget_eros_factor(npixel)



!   time constraint for CAMP degassing
    data pulse_start    /0, 175e3, 195e3, 265e3, 623e3/  !yrs
    data pulse_end      /50e3, 195e3, 265e3, 285e3, 643e3/  !yrs
    data dc13_peak      /-13e-3, -41e-3, -5e-3, -41e-3, -5e-3/
    data amount_peak    /8800, 4800, 4000, 2500, 4000/  !units=Gt/C


!   Global variables:
!   ==================================================================    

    common /IOcond/ run_name, output_directory, killing_file, COMB_init_mode, cont_input_mode
    common /error/ ERROR_HANDLING_OPTION
    common /convert/ convert2ascii
    common /jump/ ijump_climate,icount_climate, ijump_geogprint,icount_geogprint
    common /anoxia/ DOA,cp_burial
    common /axis/ ref_x_axis, ref_y_axis
    common /atmosphere/ fCO2atm_ocean,po2,fC13atm_ocean,fC13ocean_atm
    common /basin_parameter/  oce_surf, oce_surf_tot, temp_box, depth_box, depth_top_box, box_thick, press_box, &
                              surf_sedi, sedbastype_fract, accumul_capacity, indice_deep, indice_sedi, app_cont, &
                              indice_thermo, indice_surface, indice_epicont, indice_polar
    common /basin_paramdble/ appcont_fract
    common /basin_ndex/ ndeep, nnodeep, nsedi, nnosedi, nthermo, nnothermo, nsurface, nnosurface, &
                        nepicont, nnoepicont, npolar, nnopolar, nappcont, nnoappcont, nsediepicsurf, &
                        nsediepicnosurf, nsediintermed, nsedideep, nsurfnoappcont, nnosurfnoappcont
    common /basin_jbox/ jbox_deep, jbox_nodeep, jbox_sedi, jbox_nosedi, jbox_thermo, jbox_nothermo, &
                        jbox_surface, jbox_nosurface, jbox_epicont, jbox_noepicont, jbox_polar, jbox_nopolar, &
                        jbox_appcont, jbox_noappcont, jbox_sediepicsurf, jbox_sediepicnosurf, jbox_sediintermed, &
                        jbox_sedideep, jbox_surfnoappcont
    common /biopump/ carb_ratio,fbioP,fbioC,finorgC,finorgP
    common /calcite_diss_therm/ fdownt,fdowntP
    common /carb13/ dco3,dhco3,dh2co3 &
                    ,dc13atm,dctrap,edb
    common /Cfrac/ epsiCont,epsiC
    common /CLIMBER1/ Toceclimber
    common /CLIMBER2/ Tairclimber,Runclimber,GMSTclimber
    common /CLIMBER3/ areaclimber,co2climber,areaEarth,areatot
    common /coastal_depth/ hdeepC
    common /CO2_degas/ fvol,fmor,ftrap,fSO4_atmos,fSO4_ocean,fSO4_deep,fSO4_basin &
                       ,fSO4_crust,fCO2_atmos,fCO2_crust,fCO2_deep &
                       ,tstart_deg2,tend_deg2,ipeak
    common /closed_system/ clo,phosss
    common /concent/ var_diss, var_part, var_isot
    common /shell_alone/ shells,oxy_acc_fact,sulf_acc_fact
    common /condi/ volin,xMORin,xlumin,berner,fvegin,oblin &
                   ,xnoorg,xdrift,xsolar,temcondi,async,ageYprint
    common /condi_char/ solver
    common /CO2_em/ fanthros
    common /terrestrial_POC/ total_cont_POC_export, POC_export_rate
    common /cste/ areaPTAS,areaPTACS,dcdiss,aconv,fevapnow &
                  ,fog,ak1soil,ak2soil,akspsoil,betasoil
    common /deposition/ freef,fodc,fodp &
                        ,freefP,f_diss_c,f_diss_a &
                        ,fcarbP_dep,fdissol_carb,fdissol_carbP,ws,frain_sed,frain_POC,frain_PIC,frain_POP,fin_sed,fin_POC,fin_POP &
                        ,Xsed,fexch_sed,Corg_BE
    common /depos_bool/ sed_routing_flag
    common /dissolution/ F_seafloor_cdiss
    common /final_t/ tend  
    !common /hypsometry/ deephyp,shyp ! obsolete
    common /interclimate/ ATe,BTe,ARu,BRu
    common /litho1/ litho_frac
    common /litho2/ uniform_lithology
    common /lithium/ FrivLi,dLiriv,FhydLi
    common /lysodepth/ dplysc,dplysa
    common /xnewton_keep/ value_keep1,value_keep2
    common /phosphor/ fphos
    common /Phydrot/ fhydP
    common /pressure/ pco2soil,Hsoil,Pmax
    common /ratioC_Corg/ rC_Corg,reff
    common /recycling/ roxyd,ranoxyd
    common /runoffCLIMBER/ runclim,discharge
    common /erosion/ tss
    common /temperatureCLIMBER/ Tclim
    common /saline/ salin
    common /saturation/ omega_0,omega_ara_0
    common /sedirate/ St,S
    common /species/ h2co3,hco3,co3, &
                     ch,ak1,ak2, &
                     akb,akc_0,bco2,ph &
                     ,pco2_dissous
    common /solver_condi/ ts,tsta,tfi,ifog
    common /solver_condi_int/ ijump_print
    common /str_no_cste/ rSrdep,weighted_rsw
    common /sulfate/ SO4_conc,SO4soil
    common /terorg/ fodt,ept
    common /trap/ tdata_deg,data_deg
    common /vegetal_frac/ fracveg
    common /volumes/ box_vol, vol
    common /weathfl/ fsilw,fcarbw,fkerw,fpw,fbasw,fcarbsulfw,fsilsulfw,fH2SO4sulfw ! add carbonate and silicate weathering by sulphuric acid
    common /weathrt/ wth_allsil, wth_litho_wgh, wth_litho, pyr_wth, fker, fp
    common /sulfred/ fSulfRed, fO2_odc ! variables for sulfate reduction
    common /water/ F, Fxch_clim
    common /water_bool/ variable_watxch
    common /creation/ R_diss, R_part, R_isot
    common /shellsbug/ ishelfal
    common /outputs1/ COMB_ofile_name, GEOG_ofile_name
    common /outputs2/ COMB_time_dimname, GEOG_time_dimname
    common /contpxl/ list_cont_pixel, cont_basin_map, ncontpxl
    common /routing/ uniform_routing
    common /dynsoil_intg0D/   ijump_DS_integration, ijump_DS_print, icount_DS_int, icount_DS_pri
    common /dynsoil_dim/      DS_timestep, xlevs
    common /dynsoil_dble2D_1/ reg_thick,reg_x_surf,reg_tau_surf,reg_P_vol,reg_mean_age,reg_x_mean
    common /dynsoil_dble2D_2/ reg_prod,reg_eros,reg_P_diss,reg_P_eros,reg_x_surf_eros
    common /dynsoil_dble2D_3/ reg_Li_Friv, reg_Li_Fsp, reg_Li_driv
    common /dynsoil_dble3D/   reg_z_prof,reg_tau_prof
    common /dynsoil_dble1D/   slope
    common /dynsoil_bool1D/   DS_missingpoints
    common /dynsoil_intg2D/   reg_ktop
    common /dynsoil_charvar1/ DYNS_ofile_name
    common /dynsoil_charvar2/ DYNS_time_dimname, DYNS_restart_dim, DYNS_restart_var
    common /dynsoil_charvar3/ DynSoil_restart_file
    ! vegetation:
    common /veget_log/ coupling_veget
    common /veget_int/ icount_veget, ijump_veget, veget_occup
    common /veget/     veget_factor,veget_eros_factor
    ! climatic parameters
    common /climparam_int/ ijump_climparam, icount_climparam
    common /climparam/     clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5, cpvec, prev_cpvec, next_cpvec
