module netcdf_io_module
implicit none

contains


! ================================================================================================ !
! ================================================================================================ !


! Subroutines for netCDF file openning, closing, creation, end of definiiton...
! -----------------------------------------------------------------------------


  subroutine create(filename, ID, mode)
   use netcdf

    character(len=*), intent(in):: filename
    integer, intent(out):: ID
    integer, intent(in), optional:: mode
    integer:: ierr

    if (present(mode)) then
      ierr = nf90_create(filename, mode, ID)
    else
      ierr = nf90_create(filename, NF90_CLOBBER, ID)
    end if
    call nf90_check(ierr, 'Error while creating output file '//filename)

  end subroutine


  ! ========================================================================== !


  subroutine open_file(filename, ID, mode)
  use netcdf

    character(len=*), intent(in):: filename
    integer, intent(out):: ID
    integer, intent(in), optional:: mode
    integer:: ierr

    if (present(mode)) then
      ierr = nf90_open(filename, mode, ID)
    else
      ierr = nf90_open(filename, NF90_NOWRITE, ID)
    end if
    call nf90_check(ierr, 'Error while openning file '//filename)

  end subroutine


  ! ========================================================================== !


  subroutine enddef(ID)
  use netcdf

    integer, intent(in):: ID
    character(len=30):: num
    integer:: ierr

    ierr = nf90_enddef(ID)
    write(num,fmt="(I0)") ID
    call nf90_check(ierr, 'Error while end of definition of output file #'//num)

  end subroutine


  ! ========================================================================== !


  subroutine redef(ID)
  use netcdf

    integer, intent(in):: ID
    character(len=30):: num
    integer:: ierr

    ierr = nf90_redef(ID)
    write(num,fmt="(I0)") ID
    call nf90_check(ierr, 'Error while returning in definition mode of output file #'//num)

  end subroutine


  ! ========================================================================== !


  subroutine close_file(ID)
  use netcdf

    integer, intent(in):: ID
    character(len=30):: num
    integer:: ierr

    ierr = nf90_close(ID)
    write(num,fmt="(I0)") ID
    call nf90_check(ierr, 'Error while closing file #'//num)

  end subroutine


! ================================================================================================ !
! ================================================================================================ !


! Subroutines for defining dimension, variables and attributes
! ------------------------------------------------------------


  subroutine def_dim(ID, dimname, length, dimid)
    use netcdf

    integer, intent(in):: ID ! file identifier
    character(len=*), intent(in):: dimname
    integer, intent(in):: length
    integer, intent(out):: dimid
    integer:: ierr

    ierr = nf90_def_dim(ID, trim(dimname), length, dimid)
    call nf90_check(ierr, 'Error while defining dimension '//dimname)

  end subroutine


  ! ========================================================================== !


  subroutine get_dimid(alldimid, defdim, defdimid, ndim)
    integer, dimension(:), intent(in):: alldimid
    logical, dimension(:), intent(in):: defdim
    integer, dimension(:), intent(out):: defdimid
    integer, intent(out):: ndim
    integer:: k
    !
    if (size(defdim) > size(alldimid)) then
      print *, 'ERROR - in function "get_dimid" of file "netcdf_io_module.f90": "defdim" cannot have more elements than "alldimid"'
      stop 997
    end if
    !
    if (size(alldimid) > size(defdimid)) then
      print *,'ERROR - in function "get_dimid" of file "netcdf_io_module.f90": "alldimid" cannot have more elements than "defdimid"'
      stop 997
    end if
    !
    ndim = 0
    do k = 1,size(defdim)
      if (defdim(k)) then
        ndim = ndim + 1
        defdimid(ndim) = alldimid(k)
      end if
    end do
  end subroutine

  ! ---------- !

  subroutine def_var(ID, varname, vartype, alldimid, defdim, varid)
  use netcdf

    integer, intent(in):: ID ! file identifier
    character(len=*), intent(in):: varname, vartype
    integer, dimension(:), intent(in):: alldimid
    logical, dimension(:), intent(in):: defdim
    integer, intent(out):: varid
    integer, dimension(16):: loc_dimid
    integer:: nctype, ierr, ndim, i, nd

    ndim = size(alldimid) ! => determine total number of dimension in the netCDF file

    ! convert character "vartype" in NF90 integer "vartype"
    if (vartype=='integer' .or. vartype=='int' .or. vartype=='INTEGER' .or. vartype=='INT') then
      nctype = NF90_INT
    elseif (vartype=='single precision' .or. vartype=='real' .or. vartype=='float' .or. vartype=='float32' .or. vartype=='REAL'.or.&
            vartype=='FLOAT32' ) then
      nctype = NF90_FLOAT
    elseif (vartype=='double precision' .or. vartype=='double' .or. vartype=='dble' .or. vartype=='float64' .or. vartype=='DOUBLE' &
            .or. vartype=='DBLE'.or. vartype=='FLOAT64' ) then
      nctype = NF90_DOUBLE
    else
      print *, 'Error - in "def_var" in "netcdf_io_module.f90": unknown variable type "'//trim(vartype)//'"'
      stop 997
    end if

    ! There are 2 possible cases to interpret "defdim"
    !   - case #1 size(defdim)==1 => variable is defined on all dimension or none
    !   - case #2 size(defdim)==ndim => variable is defined on its own set of dimension

    if (size(defdim) == 1) then ! case #1
      if (defdim(1)) then
        nd = ndim
        loc_dimid(1:nd) = alldimid
      else
        nd = 0
      end if
    elseif (size(defdim) == ndim) then ! case #2
      call get_dimid(alldimid, defdim, loc_dimid, nd)
    else
      print *, 'Error - in subroutine "def_var" in "netcdf_io_module.f90": unexpected size of "defdim"'
    end if

    ierr=nf90_def_var(ID, varname, nctype, loc_dimid(1:nd), varid)
    call nf90_check(ierr, 'Error while defining variable "'//trim(varname)//'"')

  end subroutine


  ! ========================================================================== !


  subroutine put_att(ID, varid, attname, attribute_text, attribute_numeric, convert2type)
  use netcdf

    integer, intent(in):: ID ! file identifier
    integer, intent(in):: varid
    character(len=*), intent(in):: attname
    character(len=*), intent(in), optional:: attribute_text
    double precision, intent(in), optional:: attribute_numeric
    character(len=*), intent(in), optional:: convert2type
    integer:: ierr
    character(len=30):: num

    if (present(attribute_text)) then

      if (present(attribute_numeric)) then
        print *, 'Warning: choose "attribute_text" over "attribute_numeric" in "put_att" in "netcdf_io_module.f90'
      end if
      if (attribute_text=='') then
        ierr = NF90_NOERR
      else
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - !
        ierr = nf90_put_att(ID, varid, attname, attribute_text)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - !
      end if
  
    elseif (present(attribute_numeric)) then

      if (present(convert2type)) then
        if (convert2type=='integer' .or. convert2type=='int' .or. convert2type=='INTEGER' .or. convert2type=='INT') then
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
          ierr = nf90_put_att(ID, varid, attname, int(attribute_numeric))
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
        elseif (convert2type=='single precision' .or. convert2type=='real' .or. convert2type=='float' .or. &
                convert2type=='float32' .or. convert2type=='REAL'.or. convert2type=='FLOAT32' ) then
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
          ierr = nf90_put_att(ID, varid, attname, real(attribute_numeric))
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
        elseif (convert2type=='double precision' .or. convert2type=='double' .or. convert2type=='dble' .or. &
                convert2type=='float64' .or. convert2type=='DOUBLE' .or. convert2type=='DBLE'.or. convert2type=='FLOAT64' ) then
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
          ierr = nf90_put_att(ID, varid, attname, attribute_numeric)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
        else
          print *, 'Error - in "put_att" in "netcdf_io_module.f90": unknown variable type "'//trim(convert2type)//'"'
          stop 997
        end if
      else
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
        ierr = nf90_put_att(ID, varid, attname, attribute_numeric)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      end if

    else

      print *, 'Warning: no attribute given in "put_att" in "netcdf_io_module.f90'
      ierr = NF90_NOERR

    end if

    write(num, fmt="(I0)") varid
    call nf90_check(ierr,  'Error while putting variable #'//trim(num)//' attribute: "'//attname//'"')


  end subroutine


! ================================================================================================ !
! ================================================================================================ !


! Subroutines load variable from netCDF file
! (with one specific case for each variable type and dimension)
! -------------------------------------------------------------


  subroutine get_var(ID, varname, &
                     var_real0D, var_real1D, var_real2D, var_real3D, var_real4D, &
                     var_dble0D, var_dble1D, var_dble2D, var_dble3D, var_dble4d, &
                     units, fillval_name, fillval_real, fillval_dble, stt, cnt)
  use netcdf

    integer, intent(in):: ID ! file identifier
    character(len=*), intent(in):: varname
    real, intent(out), optional:: var_real0D
    real, dimension(:), intent(out), optional:: var_real1D
    real, dimension(:,:), intent(out), optional:: var_real2D
    real, dimension(:,:,:), intent(out), optional:: var_real3D
    real, dimension(:,:,:,:), intent(out), optional:: var_real4D
    double precision, intent(out), optional:: var_dble0D
    double precision, dimension(:), intent(out), optional:: var_dble1D
    double precision, dimension(:,:), intent(out), optional:: var_dble2D
    double precision, dimension(:,:,:), intent(out), optional:: var_dble3D
    double precision, dimension(:,:,:,:), intent(out), optional:: var_dble4D
    character(len=*), intent(out), optional:: units
    character(len=*), intent(in), optional:: fillval_name
    real, intent(out), optional:: fillval_real
    real, intent(out), optional:: fillval_dble
    integer, intent(in), dimension(:), optional:: stt, cnt
    double precision, dimension(1):: dumdble
    real, dimension(1):: dumreal

    character(len=50):: loc_fillval_name
    integer:: varid, ierr

    ! Getting variable ID
    ierr = nf90_inq_varid(ID, trim(varname), varid)
    call nf90_check(ierr, 'Error while getting "'//trim(varname)//'" variable identifiers')

    ! Loading variable
    if (present(var_real0D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, dumreal, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, dumreal)
      end if
      var_real0D = dumreal(1)
    elseif (present(var_real1D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, var_real1D, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, var_real1D)
      end if
    elseif (present(var_real2D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, var_real2D, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, var_real2D)
      end if
    elseif (present(var_real3D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, var_real3D, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, var_real3D)
      end if
    elseif (present(var_real4D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, var_real4D, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, var_real4D)
      end if
    elseif (present(var_dble0D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, dumdble, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, dumdble)
      end if
      var_dble0D = dumdble(1)
    elseif (present(var_dble1D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, var_dble1D, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, var_dble1D)
      end if
    elseif (present(var_dble2D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, var_dble2D, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, var_dble2D)
      end if
    elseif (present(var_dble3D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, var_dble3D, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, var_dble3D)
      end if
    elseif (present(var_dble4D)) then
      if (present(stt) .and. present(cnt)) then
        ierr = nf90_get_var(ID, varid, var_dble4D, start=stt, count=cnt)
      else
        ierr = nf90_get_var(ID, varid, var_dble4D)
      end if
    else
      print *, 'Error - in "get_var" in "netcdf_io_module.f90": no "intent(out)" array were given to the subroutine'
      stop 997
    end if
    call nf90_check(ierr, 'Error while getting variable "'//trim(varname)//'"')

    ! Getting optional attribute "units"
    if (present(units)) then
      ierr = nf90_get_att(ID, varid, 'units', units)
      call nf90_check(ierr, 'Error while getting variable "'//trim(varname)//'" attribute "units"')
    end if

    ! Getting optional attribute fillvalue
    if (present(fillval_real) .or. present(fillval_dble)) then
      if (present(fillval_name)) then
        loc_fillval_name = trim(fillval_name)
      else
        loc_fillval_name = '_Fillvalue'
      end if
      if (present(fillval_real)) then
        if (present(fillval_dble)) then
          print *, 'Warning - in "get_var" in "netcdf_io_module.f90": choose "fillval_real" over "fillval_dble" (both were given)'
        end if
        ierr = nf90_get_att(ID, varid, trim(loc_fillval_name), fillval_real)
      else ! present(fillval_dble)
        ierr = nf90_get_att(ID, varid, trim(loc_fillval_name), fillval_dble)
      end if
      call nf90_check(ierr, 'Error while getting variable "'//trim(varname)//'" attribute "'//trim(loc_fillval_name)//'"')
    end if

  end subroutine


! ================================================================================================ !
! ================================================================================================ !


! Subroutines writing variable in netCDF file
! (with one specific case for each variable type and dimension)
! -------------------------------------------------------------


  subroutine put_var(ID, varid, varname, &
                    var_real0D, var_real1D, var_real2D, var_real3D, var_real4D, &
                    var_dble0D, var_dble1D, var_dble2D, var_dble3D, var_dble4D, &
                    var_int0D, var_int1D, var_int2D, var_int3D, var_int4D, &
                    stt, cnt)
  use netcdf

    integer, intent(in):: ID
    integer, intent(in), optional:: varid
    character(len=*), intent(in), optional:: varname
    real, intent(in), optional:: var_real0D
    real, dimension(:), intent(in), optional:: var_real1D
    real, dimension(:,:), intent(in), optional:: var_real2D
    real, dimension(:,:,:), intent(in), optional:: var_real3D
    real, dimension(:,:,:,:), intent(in), optional:: var_real4D
    double precision, intent(in), optional:: var_dble0D
    double precision, dimension(:), intent(in), optional:: var_dble1D
    double precision, dimension(:,:), intent(in), optional:: var_dble2D
    double precision, dimension(:,:,:), intent(in), optional:: var_dble3D
    double precision, dimension(:,:,:,:), intent(in), optional:: var_dble4D
    integer, intent(in), optional:: var_int0D
    integer, dimension(:), intent(in), optional:: var_int1D
    integer, dimension(:,:), intent(in), optional:: var_int2D
    integer, dimension(:,:,:), intent(in), optional:: var_int3D
    integer, dimension(:,:,:,:), intent(in), optional:: var_int4D
    integer, dimension(:), optional:: stt, cnt
    integer:: ierr, loc_vid
    character(len=30):: num

    if (present(varid)) then
      loc_vid = varid
    elseif (present(varname)) then
      call inquire_var(ID, varname, loc_vid)
    end if

    if (present(var_real0D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, (/var_real0D/), start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, (/var_real0D/))
      end if
    elseif (present(var_real1D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_real1D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_real1D)
      end if
    elseif (present(var_real2D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_real2D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_real2D)
      end if
    elseif (present(var_real3D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_real3D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_real3D)
      end if
    elseif (present(var_real4D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_real4D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_real4D)
      end if
    elseif (present(var_dble0D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, (/var_dble0D/), start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, (/var_dble0D/))
      end if
    elseif (present(var_dble1D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_dble1D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_dble1D)
      end if
    elseif (present(var_dble2D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_dble2D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_dble2D)
      end if
    elseif (present(var_dble3D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_dble3D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_dble3D)
      end if
    elseif (present(var_dble4D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_dble4D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_dble4D)
      end if
    elseif (present(var_int0D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, (/var_int0D/), start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, (/var_int0D/))
      end if
    elseif (present(var_int1D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_int1D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_int1D)
      end if
    elseif (present(var_int2D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_int2D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_int2D)
      end if
    elseif (present(var_int3D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_int3D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_int3D)
      end if
    elseif (present(var_int4D)) then
      if ( present(stt) .and. present(cnt) ) then
        ierr=nf90_put_var(ID, loc_vid, var_int4D, start=stt, count=cnt)
      else
        ierr=nf90_put_var(ID, loc_vid, var_int4D)
      end if
    else
      print *, 'Error - in "put_var" in "netcdf_io_module.f90": no "intent(in)" array were given to the subroutine'
      stop 997
    end if

    write(num,fmt="(I0)") loc_vid
    call nf90_check(ierr, 'Error while putting variable #'//num)

 end subroutine


! ================================================================================================ !
! ================================================================================================ !


! Subroutines for getting variable or dimension ID in netCDF file openning
! + checking that netCDF fillvalue and code internal fillvalue matche
! ------------------------------------------------------------------------


  subroutine inquire_dim(ID, dimname, dimid)
  use netcdf

    integer, intent(in):: ID ! file identifier
    character(len=*), intent(in):: dimname
    integer, intent(out):: dimid
    integer:: ierr

    ierr = nf90_inq_dimid(ID, dimname, dimid)
    call nf90_check(ierr, 'Error while inquiring dimension "'//dimname//'" identifier')

  end subroutine


  ! ========================================================================== !


  subroutine inquire_var(ID, varname, varid)
  use netcdf

    integer, intent(in):: ID ! file identifier
    character(len=*), intent(in):: varname
    integer, intent(out):: varid
    integer:: ierr

    ierr=nf90_inq_varid(ID, varname, varid)
    call nf90_check(ierr, 'Error while inquiring variable "'//varname//'" identifier')

  end subroutine


  ! ========================================================================== !


  subroutine check_fillvalue( ID, varid, fillvalname, fillval )
  use netcdf

    integer, intent(in):: ID ! file identifier
    integer, dimension(:), intent(in):: varid
    character(len=*), dimension(:), intent(in):: fillvalname
    real, dimension(:), intent(inout):: fillval
    real:: fillvalue
    integer:: ierr, N, i
    character(len=30):: num

    N = size(varid)

    do i=1,N
      ierr=nf90_get_att( ID, varid(i), fillvalname(i), fillvalue )
      write(num,fmt="(I0)") varid(i)
      call nf90_check(ierr, 'Error while getting variable #'//trim(num)//' attribute "'//trim(fillvalname(i))//'"')
      if (fillvalue /= fillval(i)) then
        print *, 'WARNING: variable #'//trim(num)//' missing-value is not consistent with'
        print *, 'existing file. It will be changed to keep existing file value' 
        fillval(i) = fillvalue
      end if
    end do

  end subroutine


! ================================================================================================ !
! ================================================================================================ !


! Generic subroutine to check IO status and display error message
! ---------------------------------------------------------------

  subroutine nf90_check(ierr, message, kill)

    use netcdf
    integer, intent(in):: ierr
    character(len=*), intent(in):: message
    logical, optional, intent(in):: kill
    logical:: loc_kill

    if (present(kill)) then
      loc_kill = kill
    else
      loc_kill = .true.
    end if

    if (ierr/=NF90_NOERR) then
      print *
      print *, trim(message)
      print *, nf90_strerror(ierr)
      if (loc_kill) stop 997
    end if

  end subroutine



end module
