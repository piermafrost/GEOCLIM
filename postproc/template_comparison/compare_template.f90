program compare_template
!
! Load given run and compare it (for several key variables) to the corresponding template run.
! The run name is read in GEOCLIM's main file: "config/IO_CONDITIONS"

use netcdf
use netcdf_io_module, only: open_file, close_file, inquire_dim, nf90_check, get_var
use local_functions, only: compare_var1D, compare_var2D
implicit none

include 'path.inc' ! => get main path: 'geoclim_path'
include 'shape.inc' ! => get 'nbasin'
include 'output_size.inc' ! => get 'nCOMBoutvar'

! namelist variables
character(len=100):: run_name, check_run_name, phys_cond_file, killing_file
character(len=500):: output_directory, file_name
integer, parameter:: ndim=3
character(len=30):: vartype(nCOMBoutvar)
character(len=100):: dname(ndim), vname(nCOMBoutvar), units(nCOMBoutvar)
character(len=500):: long_name(nCOMBoutvar)
integer, dimension(nCOMBoutvar,ndim):: defdim
logical, dimension(nCOMBoutvar):: writevar
double precision, dimension(nCOMBoutvar):: fillval

! main variables
character(len=1000):: fname, fname_0
real, dimension(:), allocatable:: CO2_level,   O2_level,   sil_wth,   bas_wth,   carb_wth,   ker_wth,   pyr_wth,   phos_wth, &
                                  terr_bio_C_exp,   time
real, dimension(:), allocatable:: CO2_level_0, O2_level_0, sil_wth_0, bas_wth_0, carb_wth_0, ker_wth_0, pyr_wth_0, phos_wth_0, &
                                  terr_bio_C_exp_0, time_0
real, dimension(:,:), allocatable:: temp,   O2,   DIC,   alk,   pH,   Ca,   PO4,   PIC,   POC,   POP,   lysdpt,   &
                                    Sr,   Sr_ratio,   d13C,   Qpbox,   tsspbox,   silwpbox,   sedim_flux,   sedim_rate,   &
                                    bioprod,   carbprod,   ner_carb_bur,   pel_carb_bur,   org_C_bur,   bur_eff,   &
                                    P_org_bur,   P_phos_bur,   P_hydro_bur
real, dimension(:,:), allocatable:: temp_0, O2_0, DIC_0, alk_0, pH_0, Ca_0, PO4_0, PIC_0, POC_0, POP_0, lysdpt_0,  &
                                    Sr_0, Sr_ratio_0, d13C_0, Qpbox_0, tsspbox_0, silwpbox_0, sedim_flux_0, sedim_rate_0, &
                                    bioprod_0, carbprod_0, ner_carb_bur_0, pel_carb_bur_0, org_C_bur_0, bur_eff_0, &
                                    P_org_bur_0, P_phos_bur_0, P_hydro_bur_0
real:: delta, xref
integer:: fid, fid_0, dimid, ierr
integer:: nt, nt0, k, j
logical, dimension(38):: passed

! Namelist declaration
namelist /MAIN_INFO/ run_name, output_directory, phys_cond_file, killing_file
namelist /CMB_OUTPUT_FILE/ file_name
namelist /CMB_OUTPUT_DIM/ dname, units, vartype
namelist /CMB_OUTPUT_VAR/ vname, units, defdim, writevar, long_name, fillval, vartype




print *
print *, 'GEOCLIM template-reference comparison check'
print *, '###########################################'
print *
print *




! Load GEOCLIM OUTPUT information
! ===============================


! Open IO_CONDITION file
open(unit=1, status='old', action='read', file=geoclim_path//'config/IO_CONDITIONS')

! Read namelists
read(unit=1, nml=MAIN_INFO)
read(unit=1, nml=CMB_OUTPUT_FILE)
read(unit=1, nml=CMB_OUTPUT_DIM)
read(unit=1, nml=CMB_OUTPUT_VAR)

! close IO_CONDITION file
close(unit=1)



! Load variables
! ==============


! TIME
!
!     - new output file
fname = geoclim_path//'OUTPUT/'//trim(file_name)//trim(run_name)//'.nc'
call open_file(fname, fid)
call inquire_dim(fid, trim(dname(2)), dimid)
ierr = nf90_inquire_dimension(fid, dimid, len=nt)
call nf90_check(ierr, 'Error while inquiring ID of variable "'//trim(dname(2))//'" in file "'//trim(fname)//'"')
!
!     - Reference file:
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(file_name)//trim(run_name)//'.nc'
call open_file(fname_0, fid_0)
call inquire_dim(fid_0, trim(dname(2)), dimid)
ierr = nf90_inquire_dimension(fid_0, dimid, len=nt0)
call nf90_check(ierr, 'Error while inquiring ID of variable "'//trim(dname(2))//'" in file "'//trim(fname_0)//'"')
!
! Compare length
if (nt/=nt0) then
    print *
    print *, 'TEST NOT PASSED: time dimensions do not have the same length'
    stop
end if
!
!################################!
! ALLOCATE VARIABLES
allocate(time(nt))
allocate(time_0(nt))
allocate(CO2_level(nt))
allocate(CO2_level_0(nt))
allocate(O2_level(nt))
allocate(O2_level_0(nt))
allocate(sil_wth(nt))
allocate(sil_wth_0(nt))
allocate(bas_wth(nt))
allocate(bas_wth_0(nt))
allocate(carb_wth(nt))
allocate(carb_wth_0(nt))
allocate(ker_wth(nt))
allocate(ker_wth_0(nt))
allocate(pyr_wth(nt))
allocate(pyr_wth_0(nt))
allocate(phos_wth(nt))
allocate(phos_wth_0(nt))
allocate(terr_bio_C_exp(nt))
allocate(terr_bio_C_exp_0(nt))
!
allocate(temp(nbasin,nt))
allocate(temp_0(nbasin,nt))
allocate(O2(nbasin,nt))
allocate(O2_0(nbasin,nt))
allocate(DIC(nbasin,nt))
allocate(DIC_0(nbasin,nt))
allocate(alk(nbasin,nt))
allocate(alk_0(nbasin,nt))
allocate(pH(nbasin,nt))
allocate(pH_0(nbasin,nt))
allocate(Ca(nbasin,nt))
allocate(Ca_0(nbasin,nt))
allocate(PO4(nbasin,nt))
allocate(PO4_0(nbasin,nt))
allocate(PIC(nbasin,nt))
allocate(PIC_0(nbasin,nt))
allocate(POC(nbasin,nt))
allocate(POC_0(nbasin,nt))
allocate(POP(nbasin,nt))
allocate(POP_0(nbasin,nt))
allocate(lysdpt(nbasin,nt))
allocate(lysdpt_0(nbasin,nt))
allocate(Sr(nbasin,nt))
allocate(Sr_0(nbasin,nt))
allocate(Sr_ratio(nbasin,nt))
allocate(Sr_ratio_0(nbasin,nt))
allocate(d13C(nbasin,nt))
allocate(d13C_0(nbasin,nt))
allocate(Qpbox(nbasin,nt))
allocate(Qpbox_0(nbasin,nt))
allocate(tsspbox(nbasin,nt))
allocate(tsspbox_0(nbasin,nt))
allocate(silwpbox(nbasin,nt))
allocate(silwpbox_0(nbasin,nt))
allocate(sedim_flux(nbasin,nt))
allocate(sedim_flux_0(nbasin,nt))
allocate(sedim_rate(nbasin,nt))
allocate(sedim_rate_0(nbasin,nt))
allocate(bioprod(nbasin,nt))
allocate(bioprod_0(nbasin,nt))
allocate(carbprod(nbasin,nt))
allocate(carbprod_0(nbasin,nt))
allocate(ner_carb_bur(nbasin,nt))
allocate(ner_carb_bur_0(nbasin,nt))
allocate(pel_carb_bur(nbasin,nt))
allocate(pel_carb_bur_0(nbasin,nt))
allocate(org_C_bur(nbasin,nt))
allocate(org_C_bur_0(nbasin,nt))
allocate(bur_eff(nbasin,nt))
allocate(bur_eff_0(nbasin,nt))
allocate(P_org_bur(nbasin,nt))
allocate(P_org_bur_0(nbasin,nt))
allocate(P_phos_bur(nbasin,nt))
allocate(P_phos_bur_0(nbasin,nt))
allocate(P_hydro_bur(nbasin,nt))
allocate(P_hydro_bur_0(nbasin,nt))
!################################!

! load time variable
call get_var(fid,   dname(3), var_real1D=time)
call get_var(fid_0, dname(3), var_real1D=time_0)

! Atmospheric CO2
call get_var(fid,   vname(77), var_real1D=CO2_level)
call get_var(fid_0, vname(77), var_real1D=CO2_level_0)

! Atmospheric O2
call get_var(fid,   vname(75), var_real1D=O2_level)
call get_var(fid_0, vname(75), var_real1D=O2_level_0)

! Silicate weathering (all litho)
call get_var(fid,   vname(121), var_real1D=sil_wth)
call get_var(fid_0, vname(121), var_real1D=sil_wth_0)

! Basalt weathering
call get_var(fid,   vname(122), var_real1D=bas_wth)
call get_var(fid_0, vname(122), var_real1D=bas_wth_0)

! Carbonate weathering
call get_var(fid,   vname(123), var_real1D=carb_wth)
call get_var(fid_0, vname(123), var_real1D=carb_wth_0)

! Kerogen weathering
call get_var(fid,   vname(124), var_real1D=ker_wth)
call get_var(fid_0, vname(124), var_real1D=ker_wth_0)

! Pyrite weathering
call get_var(fid,   vname(125), var_real1D=pyr_wth)
call get_var(fid_0, vname(125), var_real1D=pyr_wth_0)

! Phosphorus weathering
call get_var(fid,   vname(128), var_real1D=phos_wth)
call get_var(fid_0, vname(128), var_real1D=phos_wth_0)

! Terrestrial biospheric organic C export
call get_var(fid,   vname(129), var_real1D=terr_bio_C_exp)
call get_var(fid_0, vname(129), var_real1D=terr_bio_C_exp_0)

! Oceanic temperature
call get_var(fid,   vname(47), var_real2D=temp)
call get_var(fid_0, vname(47), var_real2D=temp_0)

! Dissolved oxygen
call get_var(fid,   vname(25), var_real2D=O2)
call get_var(fid_0, vname(25), var_real2D=O2_0)

! Dissolved Inorganic Carbon
call get_var(fid,   vname(20), var_real2D=DIC)
call get_var(fid_0, vname(20), var_real2D=DIC_0)

! Alkalinity
call get_var(fid,   vname(21), var_real2D=alk)
call get_var(fid_0, vname(21), var_real2D=alk_0)

! PH
call get_var(fid,   vname(45), var_real2D=pH)
call get_var(fid_0, vname(45), var_real2D=pH_0)

! Dissolved calcium
call get_var(fid,   vname(23), var_real2D=Ca)
call get_var(fid_0, vname(23), var_real2D=Ca_0)

! Dissolved phosphorus
call get_var(fid,   vname(22), var_real2D=PO4)
call get_var(fid_0, vname(22), var_real2D=PO4_0)

! Particulate Inorganic Carbon
call get_var(fid,   vname(32), var_real2D=PIC)
call get_var(fid_0, vname(32), var_real2D=PIC_0)

! Particulate Organic Carbon
call get_var(fid,   vname(31), var_real2D=POC)
call get_var(fid_0, vname(31), var_real2D=POC_0)

! Particulate Organic Phosphorus
call get_var(fid,   vname(29), var_real2D=POP)
call get_var(fid_0, vname(29), var_real2D=POP_0)

! Calcite lysocline depth
call get_var(fid,   vname(49), var_real2D=lysdpt)
call get_var(fid_0, vname(49), var_real2D=lysdpt_0)

! Dissolved strontium
call get_var(fid,   vname(24), var_real2D=Sr)
call get_var(fid_0, vname(24), var_real2D=Sr_0)

! Dissolved strontium 87Sr/86Sr ratio
call get_var(fid,   vname(37), var_real2D=Sr_ratio)
call get_var(fid_0, vname(37), var_real2D=Sr_ratio_0)

! DIC delta 13C
call get_var(fid,   vname(33), var_real2D=d13C)
call get_var(fid_0, vname(33), var_real2D=d13C_0)

! Water discharge per boxes
call get_var(fid,   vname(130), var_real2D=Qpbox)
call get_var(fid_0, vname(130), var_real2D=Qpbox_0)

! Sediment flux (Total Suspende Solid) per box
call get_var(fid,   vname(131), var_real2D=tsspbox)
call get_var(fid_0, vname(131), var_real2D=tsspbox_0)

! Silicate weathering per box
call get_var(fid,   vname(132), var_real2D=silwpbox)
call get_var(fid_0, vname(132), var_real2D=silwpbox_0)

! Seaflor sedimentation flux
call get_var(fid,   vname(93), var_real2D=sedim_flux)
call get_var(fid_0, vname(93), var_real2D=sedim_flux_0)

! Seaflor sedimentation rate
call get_var(fid,   vname(94), var_real2D=sedim_rate)
call get_var(fid_0, vname(94), var_real2D=sedim_rate_0)

! Organic primary productivity
call get_var(fid,   vname(86), var_real2D=bioprod)
call get_var(fid_0, vname(86), var_real2D=bioprod_0)

! Carbonate primary productivity
call get_var(fid,   vname(85), var_real2D=carbprod)
call get_var(fid_0, vname(85), var_real2D=carbprod_0)

! Neritic carboante burial flux
call get_var(fid,   vname(100), var_real2D=ner_carb_bur)
call get_var(fid_0, vname(100), var_real2D=ner_carb_bur_0)

! Pelagic carbonate burial flux
call get_var(fid,   vname(101), var_real2D=pel_carb_bur)
call get_var(fid_0, vname(101), var_real2D=pel_carb_bur_0)

! Organic C burial flux
call get_var(fid,   vname(102), var_real2D=org_C_bur)
call get_var(fid_0, vname(102), var_real2D=org_C_bur_0)

! Burial efficiency
call get_var(fid,   vname(103), var_real2D=bur_eff)
call get_var(fid_0, vname(103), var_real2D=bur_eff_0)

! P burial associated with POC
call get_var(fid,   vname(104), var_real2D=P_org_bur)
call get_var(fid_0, vname(104), var_real2D=P_org_bur_0)

! P burial in form of phosphorite
call get_var(fid,   vname(105), var_real2D=P_phos_bur)
call get_var(fid_0, vname(105), var_real2D=P_phos_bur_0)

! P burial associated with hydothermal Fe
call get_var(fid,   vname(106), var_real2D=P_hydro_bur)
call get_var(fid_0, vname(106), var_real2D=P_hydro_bur_0)


! Close netCDF files
call close_file(fid)
call close_file(fid_0)




! Compare variables and print report
! ==================================

print *
print *

passed(1)  = compare_var1D('Time',                           time,           time_0)
passed(2)  = compare_var1D('Atmospheric pCO2',               CO2_level,      CO2_level_0)
passed(3)  = compare_var1D('Atmospheric pO2',                O2_level,       O2_level_0)
passed(4)  = compare_var1D('Silicate weathering (all)',      sil_wth,        sil_wth_0)
passed(5)  = compare_var1D('Basalt weathering',              bas_wth,        bas_wth_0)
passed(6)  = compare_var1D('Carbonate weathering',           carb_wth,       carb_wth_0)
passed(7)  = compare_var1D('Kerogen weathering',             ker_wth,        ker_wth_0)
passed(8)  = compare_var1D('Pyrite weathering',              pyr_wth,        pyr_wth_0)
passed(9)  = compare_var1D('Phosphorus weathering',          phos_wth,       phos_wth_0)
passed(10) = compare_var1D('Land biospheric POC export',     terr_bio_C_exp, terr_bio_C_exp_0)
passed(11) = compare_var2D('Oceanic temperature',            temp,           temp_0)
passed(12) = compare_var2D('Dissolved O2',                   O2,             O2_0)
passed(13) = compare_var2D('DIC',                            DIC,            DIC_0)
passed(14) = compare_var2D('Alkalinity',                     alk,            alk_0)
passed(15) = compare_var2D('PH',                             pH,             pH_0)
passed(16) = compare_var2D('Dissolved Ca',                   Ca,             Ca_0)
passed(17) = compare_var2D('Dissolved P',                    PO4,            PO4_0)
passed(18) = compare_var2D('PIC',                            PIC,            PIC_0)
passed(19) = compare_var2D('POC',                            POC,            POC_0)
passed(20) = compare_var2D('POP',                            POP,            POP_0)
passed(21) = compare_var2D('Calcite lysocline depth',        lysdpt,         lysdpt_0)
passed(22) = compare_var2D('Dissolved Sr',                   Sr,             Sr_0)
passed(23) = compare_var2D('Dissolved Sr - 87Sr/86Sr ratio', Sr_ratio,       Sr_ratio_0)
passed(24) = compare_var2D('DIC delta 13C',                  d13C,           d13C_0)
passed(25) = compare_var2D('Water discharge per boxes',      Qpbox,          Qpbox_0)
passed(26) = compare_var2D('Sediment TSS per box',           tsspbox,        tsspbox_0)
passed(27) = compare_var2D('Silicate weathering per box',    silwpbox,       silwpbox_0)
passed(28) = compare_var2D('Seafloor sedimentation flux',    sedim_flux,     sedim_flux_0)
passed(29) = compare_var2D('Seafloor sedimentation rate',    sedim_rate,     sedim_rate_0)
passed(30) = compare_var2D('Organic primary productivity',   bioprod,        bioprod_0)
passed(31) = compare_var2D('Carbonate primary productivity', carbprod,       carbprod_0)
passed(32) = compare_var2D('Neritic carboante burial flux',  ner_carb_bur,   ner_carb_bur_0)
passed(33) = compare_var2D('Pelagic carbonate burial flux',  pel_carb_bur,   pel_carb_bur_0)
passed(34) = compare_var2D('Organic C burial flux',          org_C_bur,      org_C_bur_0)
passed(35) = compare_var2D('Burial efficiency',              bur_eff,        bur_eff_0)
passed(36) = compare_var2D('P burial assoc. with POC',       P_org_bur,      P_org_bur_0)
passed(37) = compare_var2D('Phosphorite burial',             P_phos_bur,     P_phos_bur_0)
passed(38) = compare_var2D('P burial assooc. with hyd. Fe',  P_hydro_bur,    P_hydro_bur_0)




print *
print *
if (all(passed)) then
    print *, 'Template comparison check PASSED'
else
    print *, 'Template comparison check FAILED'
end if



end program
