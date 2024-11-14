subroutine read_conditions
! ************************

    use io_module, only: UNDEFINED_VALUE_INT, UNDEFINED_VALUE_CHAR, UNDEFINED_VALUE_DBLE, check_namelist_def
    implicit none
    include 'combine_foam.inc'
    ! Local variables
    double precision :: tot_degass_in, fract_MORin

    ! Namelist declaration
    namelist /COND_P20/ tot_degass_in, volin, xMORin, fract_MORin, fSO4_atmos, fSO4_ocean, fSO4_deep, fCO2_atmos, fCO2_deep, &
                        ipeak, shells, ishelfal, clo, phosss, xnoorg, temcondi, &
                        coupling_veget, &
                        solver, ts, tsta, tfi, &
                        oxy_acc_fact, sulf_acc_fact, &
                        ijump_climate, ijump_DS_integration, ijump_climparam, ijump_veget, &
                        ijump_print, ijump_geogprint, ijump_DS_print, ageYprint, convert2ascii

    ! Default values of namelist variables
    tot_degass_in        = UNDEFINED_VALUE_DBLE
    volin                = UNDEFINED_VALUE_DBLE
    xMORin               = UNDEFINED_VALUE_DBLE
    fract_MORin          = UNDEFINED_VALUE_DBLE
    fSO4_atmos           = UNDEFINED_VALUE_DBLE
    fSO4_ocean           = UNDEFINED_VALUE_DBLE
    fSO4_deep            = UNDEFINED_VALUE_DBLE
    fCO2_atmos           = UNDEFINED_VALUE_DBLE
    fCO2_deep            = UNDEFINED_VALUE_DBLE
    ipeak                = UNDEFINED_VALUE_INT
    shells               = UNDEFINED_VALUE_DBLE
    ishelfal             = UNDEFINED_VALUE_INT
    clo                  = UNDEFINED_VALUE_DBLE
    phosss               = UNDEFINED_VALUE_DBLE
    xnoorg               = UNDEFINED_VALUE_DBLE
    temcondi             = UNDEFINED_VALUE_DBLE
    coupling_veget       = .false. ! Default behaviour: no coupling with vegetation if nothing is specified in config file
    solver               = UNDEFINED_VALUE_CHAR
    ts                   = UNDEFINED_VALUE_DBLE
    tsta                 = UNDEFINED_VALUE_DBLE
    tfi                  = UNDEFINED_VALUE_DBLE
    oxy_acc_fact         = UNDEFINED_VALUE_DBLE
    sulf_acc_fact        = UNDEFINED_VALUE_DBLE
    ijump_climate        = UNDEFINED_VALUE_INT
    ijump_DS_integration = UNDEFINED_VALUE_INT
    ijump_climparam      = UNDEFINED_VALUE_INT
    ijump_veget          = UNDEFINED_VALUE_INT
    ijump_print          = UNDEFINED_VALUE_INT
    ijump_geogprint      = UNDEFINED_VALUE_INT
    ijump_DS_print       = UNDEFINED_VALUE_INT
    ageYprint            = UNDEFINED_VALUE_DBLE
    convert2ascii        = .false.! Default behaviour: no output conversion in ascii format if nothing is specified in config file


    ! Read namelist variables in physical condition file
    ! <><><><><><><><><><><><><><><><><><><><><><><><> !
    read(unit=3, nml=COND_P20)
    ! <><><><><><><><><><><><><><><><><><><><><><><><> !


    ! Check that all the needed variables were defined
    ! ------------------------------------------------

    if (fract_MORin == UNDEFINED_VALUE_DBLE) then
        !
        call check_namelist_def('Error - in read_conditions.f: neither "fract_MORin" nor "xMORin" were given in physical'// &
                                'conditions config file', dble_var=xMORin)
        !
        if (tot_degass_in == UNDEFINED_VALUE_DBLE) then
            call check_namelist_def('Error - in read_conditions.f: neither "tot_degass_in" nor "volin" were given in physical'// &
                                    'conditions config file', dble_var=volin)
        else
            ! Define "volin" by subtracting MORB degassing from total degassing
            if (volin /= UNDEFINED_VALUE_DBLE) then
                print *
                print *, 'Warning: gotten value of "volin" variable will be overriden by "tot_degass_in"'
            end if
            volin = tot_degass_in - xMORin
        end if
        !
    else
        !
        if (xMORin /= UNDEFINED_VALUE_DBLE) then
            print *
            print *, 'Warning: gotten value of "xMORin" variable will be overriden by "fract_MORin"'
        end if
        if (tot_degass_in == UNDEFINED_VALUE_DBLE) then
            call check_namelist_def('Error - in read_conditions.f: neither "tot_degass_in" nor "volin" were given in physical'// &
                                    'conditions config file', dble_var=volin)
            ! Define "xMORin" from the fraction of total degassing (xMORin + volin)
            xMORin = volin * fract_MORin/(1-fract_MORin)
        else
            if (volin /= UNDEFINED_VALUE_DBLE) then
                print *
                print *, 'Warning: gotten value of "volin" variable will be overriden by "tot_degass_in"'
            end if
            ! Define "xMORin" from the fraction of total degassing and "volin" by subtracting MORB degassing from total degassing
            xMORin = fract_MORIN*tot_degass_in
            volin = tot_degass_in - xMORin
        end if
        !
    end if
    !
    call check_namelist_def('Error - in read_conditions.f: variable "fSO4_atmos" was not given in physical conditions config file',&
                            dble_var=fSO4_atmos)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "fSO4_ocean" was not given in physical conditions config file',&
                            dble_var=fSO4_ocean)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "fSO4_deep" was not given in physical conditions config file', &
                            dble_var=fSO4_deep)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "fCO2_atmos" was not given in physical conditions config file',&
                            dble_var=fCO2_atmos)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "fCO2_deep" was not given in physical conditions config file', &
                            dble_var=fCO2_deep)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ipeak" was not given in physical conditions config file', &
                            int_var=ipeak)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "shells" was not given in physical conditions config file', &
                            dble_var=shells)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ishelfal" was not given in physical conditions config file', &
                            int_var=ishelfal)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "clo" was not given in physical conditions config file', &
                            dble_var=clo)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "phosss" was not given in physical conditions config file', &
                            dble_var=phosss)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "xnoorg" was not given in physical conditions config file', &
                            dble_var=xnoorg)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "temcondi" was not given in physical conditions config file', &
                            dble_var=temcondi)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "solver" was not given in physical conditions config file', &
                            char_var=solver)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ts" was not given in physical conditions config file', &
                            dble_var=ts)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "tsta" was not given in physical conditions config file', &
                            dble_var=tsta)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "tfi" was not given in physical conditions config file', &
                            dble_var=tfi)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "oxy_acc_fact" was not given in physical conditions'// &
                            ' config file', dble_var=oxy_acc_fact)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "sulf_acc_fact" was not given in physical conditions'// &
                            ' config file', dble_var=sulf_acc_fact)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ijump_climate" was not given in physical conditions'// &
                            ' config file', int_var=ijump_climate)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ijump_DS_integration" was not given in physical conditions'// &
                            ' config file', int_var=ijump_DS_integration)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ijump_climparam" was not given in physical conditions'// &
                            ' config file', int_var=ijump_climparam)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ijump_veget" was not given in physical conditions'// &
                            ' config file', int_var=ijump_veget)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ijump_print" was not given in physical conditions'// &
                            ' config file', int_var=ijump_print)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ijump_geogprint" was not given in physical conditions'// &
                            ' config file', int_var=ijump_geogprint)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ijump_DS_print" was not given in physical conditions'// &
                            ' config file', int_var=ijump_DS_print)
    !
    call check_namelist_def('Error - in read_conditions.f: variable "ageYprint" was not given in physical conditions config file', &
                            dble_var=ageYprint)


    ! generate restart at end of run if asked (0 or negative time)
    if (ageYprint <= 0) ageYprint = tfi + ts


end subroutine

