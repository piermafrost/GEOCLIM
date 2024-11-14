module geoclim_create_output_mod
implicit none

contains

subroutine geoclim_create_output(ID, output_path, run_name, &
                                 indice_surface, indice_thermo, indice_deep, indice_polar, indice_epicont, app_cont, indice_sedi, &
                                 vol, surf, surf_sedi, depth_box, depth_top_box, accumul_capacity, &
                                 file_name, time_dimname, outvar_info)

  use io_module, only: netcdf_output_var, set_default_nml_values, set_outvar_info, &
                       check_namelist_def, UNDEFINED_VALUE_INT, UNDEFINED_VALUE_CHAR, DEFAULT_FILLVAL_NAME
  use netcdf_io_module, only: create, def_dim, def_var, put_att, enddef, put_var, close_file
  use netcdf

  include 'coupler.inc'
  include 'shape.inc' ! => nbasin
  include 'output_size.inc' ! => nCOMBouvar

  integer, intent(in):: ID
  character(len=*), intent(in):: output_path, run_name
  integer, dimension(nbasin), intent(in):: indice_surface,indice_thermo,indice_deep,indice_polar,indice_epicont,app_cont,indice_sedi
  double precision, dimension(nbasin), intent(in):: vol, surf, surf_sedi, depth_box, depth_top_box, accumul_capacity
  !
  character(len=500), intent(out):: file_name
  character(len=100), intent(out):: time_dimname
  type(netcdf_output_var), dimension(nCOMBoutvar), intent(out):: outvar_info

  integer, parameter:: ndim=3
  character(len=30):: vartype(nCOMBoutvar)
  character(len=100):: dname(ndim), vname(nCOMBoutvar), units(nCOMBoutvar)
  character(len=500):: long_name(nCOMBoutvar)
  integer, dimension(nCOMBoutvar,ndim):: defdim
  logical, dimension(nCOMBoutvar):: writevar
  double precision, dimension(nCOMBoutvar):: fillval

  character(len=3):: num
  integer:: fid, dimvarid(ndim), dimid(ndim)
  integer:: i

  ! Namelists declaration
  ! ---------------------
  namelist /CMB_OUTPUT_FILE/ file_name
  namelist /CMB_OUTPUT_DIM/ dname, units, vartype
  namelist /CMB_OUTPUT_VAR/ vname, units, defdim, writevar, long_name, fillval, vartype



  !%%%%%%%%%%%%%%%%%!
  !%  output file  %!
  !%%%%%%%%%%%%%%%%%!

  file_name = UNDEFINED_VALUE_CHAR
  !<><><><><><><><><><><><><><><><>!
  read(unit=ID, nml=CMB_OUTPUT_FILE)
  !<><><><><><><><><><><><><><><><>!
  call check_namelist_def('Error: "file_name" not defined in "CMB_OUTPUT_FILE" namelist of main config file', char_var=file_name)

  file_name = trim(output_path)//trim(file_name)//trim(run_name)//'.nc'

  ! Create file
  call create(file_name, fid)

  ! Global attributes
  call put_att(    fid, NF90_GLOBAL, 'title',                     attribute_text='COMBINE module ouputs')
  call put_att(    fid, NF90_GLOBAL, 'run_name',                  attribute_text=trim(run_name))
  if (coupling_dynsoil) then
    call put_att(  fid, NF90_GLOBAL, 'silicate_weathering_model', attribute_text='DynSoil')
    if (use_dynsoil_steady_state) then
      call put_att(fid, NF90_GLOBAL, 'DynSoil_mode',              attribute_text='steady-state')
    else
      call put_att(fid, NF90_GLOBAL, 'DynSoil_mode',              attribute_text='dynamic')
    end if
  else
    call put_att(  fid, NF90_GLOBAL, 'silicate_weathering_model', attribute_text='old')
  endif
  call put_att(    fid, NF90_GLOBAL, 'CO2_interpolation',         attribute_text=CO2_interpolation)



  !%%%%%%%%%%%%%%%%%%%%%%%!
  !%  output dimensions  %!
  !%%%%%%%%%%%%%%%%%%%%%%%!


  ! Put default "undefined" values in namelist variables that will be read
  do i = 1,ndim
    call set_default_nml_values(dname(i), units(i), vartype=vartype(i))
  end do

  ! Read all information
  !<><><><><><><><><><><><><><><><>!
  read(unit=ID, nml=CMB_OUTPUT_DIM)
  !<><><><><><><><><><><><><><><><>!

  ! dim #1: boxes
  i = 1
  !
  write(num, fmt="(I0)") i
  call check_namelist_def('Error: "dname('//trim(num)//')" not defined in "CMB_OUTPUT_DIM" namelist of main config file', &
                          char_var=dname(i))
  call check_namelist_def('Error: "units('//trim(num)//')" not defined in "CMB_OUTPUT_DIM" namelist of main config file', &
                          char_var=units(i))
  !
  call def_dim(fid, dname(i), nbasin, dimid(i))
  call def_var(fid, dname(i), vartype(i), dimid(i:i), (/.true./), dimvarid(i))
  ! dimension variable attributes:
  call put_att(fid, dimvarid(i), 'axis',   attribute_text='X')
  call put_att(fid, dimvarid(i), 'name',   attribute_text=dname(i))
  call put_att(fid, dimvarid(i), 'units',  attribute_text=units(i))

  ! dim #1: boxes-2
  i = 2
  !
  write(num, fmt="(I0)") i
  call check_namelist_def('Error: "dname('//trim(num)//')" not defined in "CMB_OUTPUT_DIM" namelist of main config file', &
                          char_var=dname(i))
  call check_namelist_def('Error: "units('//trim(num)//')" not defined in "CMB_OUTPUT_DIM" namelist of main config file', &
                          char_var=units(i))
  !
  call def_dim(fid, dname(i), nbasin, dimid(i))
  call def_var(fid, dname(i), vartype(i), dimid(i:i), (/.true./), dimvarid(i))
  ! dimension variable attributes:
  call put_att(fid, dimvarid(i), 'axis',   attribute_text='Y')
  call put_att(fid, dimvarid(i), 'name',   attribute_text=dname(i))
  call put_att(fid, dimvarid(i), 'units',  attribute_text=units(i))

  ! dim #3: time
  i = 3
  !
  write(num, fmt="(I0)") i
  call check_namelist_def('Error: "dname('//trim(num)//')" not defined in "CMB_OUTPUT_DIM" namelist of main config file', &
                          char_var=dname(i))
  call check_namelist_def('Error: "units('//trim(num)//')" not defined in "CMB_OUTPUT_DIM" namelist of main config file', &
                          char_var=units(i))
  !
  call def_dim(fid, dname(i), NF90_UNLIMITED, dimid(i))
  call def_var(fid, dname(i), vartype(i), dimid(i:i), (/.true./), dimvarid(i))
  ! dimension variable attributes:
  call put_att(fid, dimvarid(i), 'axis',  attribute_text='T')
  call put_att(fid, dimvarid(i), 'name',  attribute_text=dname(i))
  call put_att(fid, dimvarid(i), 'units', attribute_text=units(i))
  !,,,,,,,,,,,,,,,,,,,,,!
  time_dimname = dname(i)
  !'''''''''''''''''''''!



  !%%%%%%%%%%%%%%%%%%%%%%!
  !%  output variables  %!
  !%%%%%%%%%%%%%%%%%%%%%%!


  ! Put default "undefined" values in namelist variables that will be read
  do i = 1,nCOMBoutvar
    call set_default_nml_values(vname(i), units(i), defdim(i,:), writevar(i), long_name(i), fillval(i), vartype(i))
  end do

  ! read all information
  !<><><><><><><><><><><><><><><><>!
  read(unit=ID, nml=CMB_OUTPUT_VAR)
  !<><><><><><><><><><><><><><><><>!

  do i = 1,nCOMBoutvar

    ! Put output variable information in "outvar" structure
    outvar_info(i) = set_outvar_info(vname(i), units(i), defdim(i,:), writevar(i), long_name(i), fillval(i), vartype(i))

    if (writevar(i)) then

      write(num, fmt="(I0)") i
      call check_namelist_def('Error: "vname('//trim(num)//')" not defined in CMB_OUTPUT_VAR namelist of config/IO_CONDITIONS', &
                               char_var=vname(i))
      call check_namelist_def('Error: "units('//trim(num)//')" not defined in CMB_OUTPUT_VAR namelist of config/IO_CONDITIONS', &
                               char_var=units(i))
      call check_namelist_def('Error: "defdim('//trim(num)//',:)" not defined in CMB_OUTPUT_VAR namelist of config/IO_CONDITIONS', &
                               int_var=defdim(i,1))

      ! define netcdf variable
      call def_var(fid, outvar_info(i)%vname, outvar_info(i)%vartype, &
                   dimid, outvar_info(i)%def_dim(1:outvar_info(i)%ndim_tot), &
                   outvar_info(i)%id)

      ! put all netcdf attributes
      call put_att(fid, outvar_info(i)%id, 'name',               attribute_text=outvar_info(i)%vname)
      call put_att(fid, outvar_info(i)%id, 'units',              attribute_text=outvar_info(i)%units)
      call put_att(fid, outvar_info(i)%id, 'long_name',          attribute_text=outvar_info(i)%long_name)
      call put_att(fid, outvar_info(i)%id, DEFAULT_FILLVAL_NAME, attribute_numeric=outvar_info(i)%fillval, &
                                                                      convert2type=outvar_info(i)%vartype  )

    end if
  end do



  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% end of output file definition %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


  call enddef(fid)


  ! write dimension variables and time-independent variables (volume and surfaces):
  !--------------------------------------------------------------------------------

  call put_var(fid, dimvarid(1), var_int1D=(/(i,i=1,nbasin)/))
  call put_var(fid, dimvarid(2), var_int1D=(/(i,i=1,nbasin)/))
  if (outvar_info(1)%writevar)  call put_var(fid, outvar_info(1)%id,  var_int1D=indice_surface)
  if (outvar_info(2)%writevar)  call put_var(fid, outvar_info(2)%id,  var_int1D=indice_thermo)
  if (outvar_info(3)%writevar)  call put_var(fid, outvar_info(3)%id,  var_int1D=indice_deep)
  if (outvar_info(4)%writevar)  call put_var(fid, outvar_info(4)%id,  var_int1D=indice_polar)
  if (outvar_info(5)%writevar)  call put_var(fid, outvar_info(5)%id,  var_int1D=indice_epicont)
  if (outvar_info(6)%writevar)  call put_var(fid, outvar_info(6)%id,  var_int1D=app_cont)
  if (outvar_info(7)%writevar)  call put_var(fid, outvar_info(7)%id,  var_int1D=indice_sedi)
  if (outvar_info(8)%writevar)  call put_var(fid, outvar_info(8)%id,  var_real1D=real(vol))
  if (outvar_info(9)%writevar)  call put_var(fid, outvar_info(9)%id,  var_real1D=real(surf))
  if (outvar_info(10)%writevar) call put_var(fid, outvar_info(10)%id, var_real1D=real(surf_sedi))
  if (outvar_info(11)%writevar) call put_var(fid, outvar_info(11)%id, var_real1D=real(depth_box))
  if (outvar_info(12)%writevar) call put_var(fid, outvar_info(12)%id, var_real1D=real(depth_top_box))
  if (outvar_info(13)%writevar) call put_var(fid, outvar_info(13)%id, var_real1D=real(accumul_capacity))


  ! close output file:
  !-------------------

  call close_file(fid)


  ! Erase variable IDs
  ! (they are meaningless once the netCDF output file is closed)
  ! ------------------

  do i = 1,nCOMBoutvar
    outvar_info(i)%id = UNDEFINED_VALUE_INT
  end do



 end subroutine

end module
