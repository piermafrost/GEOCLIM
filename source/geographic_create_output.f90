module geographic_create_output_mod
implicit none

contains

subroutine geographic_create_output(ID, output_path, run_name, lon, lat, area, litho_frac, slope, cont_basin_map, &
                                    file_name, time_dimname, outvar_info                                          )

  use io_module, only: netcdf_output_var, set_default_nml_values, set_outvar_info, &
                       check_namelist_def, UNDEFINED_VALUE_INT, UNDEFINED_VALUE_CHAR, DEFAULT_FILLVAL_NAME
  use netcdf_io_module, only: create, def_dim, def_var, put_att, enddef, put_var, close_file
  use netcdf

  include 'coupler.inc'
  include 'shape.inc' ! => nlon, nlat, nlitho
  include 'output_size.inc' ! => nGEOGoutvar

  integer, intent(in):: ID
  character(len=*), intent(in):: output_path, run_name
  double precision, intent(in), dimension(nlon):: lon
  double precision, intent(in), dimension(nlat):: lat
  double precision, intent(in), dimension(nlon*nlat):: area, slope
  double precision, intent(in), dimension(nlitho,nlon*nlat):: litho_frac
  integer, intent(inout), dimension(nlon*nlat):: cont_basin_map
  !
  character(len=500), intent(out):: file_name
  character(len=100), intent(out):: time_dimname
  type(netcdf_output_var), dimension(nGEOGoutvar), intent(out):: outvar_info

  integer, parameter:: ndim=4
  character(len=30):: vartype(nGEOGoutvar)
  character(len=100):: dname(ndim), vname(nGEOGoutvar), units(nGEOGoutvar)
  character(len=500):: long_name(nGEOGoutvar)
  integer, dimension(nGEOGoutvar,ndim):: defdim
  logical, dimension(nGEOGoutvar):: writevar
  double precision, dimension(nGEOGoutvar):: fillval

  character(len=3):: num
  integer:: fid, dimvarid(ndim), dimid(ndim)
  integer:: i

  ! Namelists declaration
  ! ---------------------
  namelist /GEO_OUTPUT_FILE/ file_name
  namelist /GEO_OUTPUT_DIM/ dname, units, vartype
  namelist /GEO_OUTPUT_VAR/ vname, units, defdim, writevar, long_name, fillval, vartype



  !%%%%%%%%%%%%%%%%%!
  !%  output file  %!
  !%%%%%%%%%%%%%%%%%!

  file_name = UNDEFINED_VALUE_CHAR
  !<><><><><><><><><><><><><><><><>!
  read(unit=ID, nml=GEO_OUTPUT_FILE)
  !<><><><><><><><><><><><><><><><>!
  call check_namelist_def('Error: "file_name" not defined in "GEO_OUTPUT_FILE" namelist of main config file', char_var=file_name)

  file_name = trim(output_path)//trim(file_name)//trim(run_name)//'.nc'

  ! Create file
  call create(file_name, fid)

  ! Global attributes
  call put_att(    fid, NF90_GLOBAL, 'title',                     attribute_text='Main continental ouputs')
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
  read(unit=ID, nml=GEO_OUTPUT_DIM)
  !<><><><><><><><><><><><><><><><>!

  ! dim #1: horizontal "x" axis (e.g., longitude)
  i = 1
  !
  write(num, fmt="(I0)") i
  call check_namelist_def('Error: "dname('//trim(num)//')" not defined in "GEO_OUTPUT_DIM" namelist of main config file', &
                          char_var=dname(i))
  call check_namelist_def('Error: "units('//trim(num)//')" not defined in "GEO_OUTPUT_DIM" namelist of main config file', &
                          char_var=units(i))
  !
  call def_dim(fid, dname(i), nlon, dimid(i))
  call def_var(fid, dname(i), vartype(i), dimid(i:i), (/.true./), dimvarid(i))
  ! dimension variable attributes:
  call put_att(fid, dimvarid(i), 'axis',      attribute_text='X')
  call put_att(fid, dimvarid(i), 'nav_model', attribute_text='Default grid')
  call put_att(fid, dimvarid(i), 'name',      attribute_text=dname(i))
  call put_att(fid, dimvarid(i), 'units',     attribute_text=units(i))

  ! dim #2: horizontal "y" axis (e.g., latitude)
  i = 2
  !
  write(num, fmt="(I0)") i
  call check_namelist_def('Error: "dname('//trim(num)//')" not defined in "GEO_OUTPUT_DIM" namelist of main config file', &
                          char_var=dname(i))
  call check_namelist_def('Error: "units('//trim(num)//')" not defined in "GEO_OUTPUT_DIM" namelist of main config file', &
                          char_var=units(i))
  !
  call def_dim(fid, dname(i), nlat, dimid(i))
  call def_var(fid, dname(i), vartype(i), dimid(i:i), (/.true./), dimvarid(i))
  ! dimension variable attributes:
  call put_att(fid, dimvarid(i), 'axis',      attribute_text='Y')
  call put_att(fid, dimvarid(i), 'nav_model', attribute_text='Default grid')
  call put_att(fid, dimvarid(i), 'name',      attribute_text=dname(i))
  call put_att(fid, dimvarid(i), 'units',     attribute_text=units(i))

  ! dim #3: lithology classes
  i = 3
  !
  write(num, fmt="(I0)") i
  call check_namelist_def('Error: "dname('//trim(num)//')" not defined in "GEO_OUTPUT_DIM" namelist of main config file', &
                          char_var=dname(i))
  call check_namelist_def('Error: "units('//trim(num)//')" not defined in "GEO_OUTPUT_DIM" namelist of main config file', &
                          char_var=units(i))
  !
  call def_dim(fid, dname(i), nlitho, dimid(i))
  call def_var(fid, dname(i), vartype(i), dimid(i:i), (/.true./), dimvarid(i))
  ! dimension variable attributes:
  call put_att(fid, dimvarid(i), 'name',  attribute_text=dname(i))
  call put_att(fid, dimvarid(i), 'units', attribute_text=units(i))

  ! dim #4: time
  i = 4
  !
  write(num, fmt="(I0)") i
  call check_namelist_def('Error: "dname('//trim(num)//')" not defined in "GEO_OUTPUT_DIM" namelist of main config file', &
                          char_var=dname(i))
  call check_namelist_def('Error: "units('//trim(num)//')" not defined in "GEO_OUTPUT_DIM" namelist of main config file', &
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
  do i = 1,nGEOGoutvar
    call set_default_nml_values(vname(i), units(i), defdim(i,:), writevar(i), long_name(i), fillval(i), vartype(i))
  end do

  ! read all information
  !<><><><><><><><><><><><><><><><>!
  read(unit=ID, nml=GEO_OUTPUT_VAR)
  !<><><><><><><><><><><><><><><><>!

  !   .   SPECIAL CASE FOR slope VARIABLE (#11): 
  !  /!\  -> do not write it if DynSoil module is not activated
  ! '''''
  ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - !
  if (.not. coupling_dynsoil) writevar(11) = .false.  !
  ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - !

  do i = 1,nGEOGoutvar

    ! Put output variable information in "outvar" structure
    outvar_info(i) = set_outvar_info(vname(i), units(i), defdim(i,:), writevar(i), long_name(i), fillval(i), vartype(i))

    if (writevar(i)) then

      write(num, fmt="(I0)") i
      call check_namelist_def('Error: "vname('//trim(num)//')" not defined in GEO_OUTPUT_VAR namelist of config/IO_CONDITIONS', &
                               char_var=vname(i))
      call check_namelist_def('Error: "units('//trim(num)//')" not defined in GEO_OUTPUT_VAR namelist of config/IO_CONDITIONS', &
                               char_var=units(i))
      call check_namelist_def('Error: "defdim('//trim(num)//',:)" not defined in GEO_OUTPUT_VAR namelist of config/IO_CONDITIONS', &
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


  ! write dimension variables and time-independent variables (area, litho_frac and slope):
  !---------------------------------------------------------------------------------------

  call put_var(fid, dimvarid(1), var_real1D=real(lon)       )
  call put_var(fid, dimvarid(2), var_real1D=real(lat)       )
  call put_var(fid, dimvarid(3),  var_int1D=(/(i,i=1,nlitho)/))
  if (outvar_info(1)%writevar)  call put_var(fid, outvar_info(1)%id,  var_real2D=real(1e12*reshape(area, shape=(/nlon,nlat/))))
                                                                                    ! ^^^^ --> convert 1e6km2 => m2
  if (outvar_info(2)%writevar)  call put_var(fid, outvar_info(2)%id,  var_real3D=real(reshape(litho_frac, &
                                                                                      shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))))
  ! put fill-value on "basinmap" before writing it
  where (area==0d0) cont_basin_map = int(outvar_info(3)%fillval)
  if (outvar_info(3)%writevar)  call put_var(fid, outvar_info(3)%id,  var_int2D=reshape(cont_basin_map, shape=(/nlon,nlat/)))
  if (outvar_info(11)%writevar) call put_var(fid, outvar_info(11)%id, var_real2D=real(reshape(slope, shape=(/nlon,nlat/))))


  ! close output file:
  !-------------------

  call close_file(fid)


  ! Erase variable IDs
  ! (they are meaningless once the netCDF output file is closed)
  ! ------------------

  do i = 1,nGEOGoutvar
    outvar_info(i)%id = UNDEFINED_VALUE_INT
  end do



 end subroutine

end module

