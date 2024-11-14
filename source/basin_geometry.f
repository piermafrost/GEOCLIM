    subroutine basin_geometry()
!    **************************
    use constante, only: ksed
    use io_module, only: DEFAULT_FILLVAL
    implicit none
    include 'combine_foam.inc'
    ! local variables:
    double precision:: area_appcont, area_sediepicsurf, area_sediepicnosurf, area_sediintermed, area_sedideep

!   volumes in cubic meters (read in 1e6 km3):
    do k=1,nbasin
        read(33,*) box_vol(k)
        box_vol(k) = box_vol(k)*1.d+15 ! 1e6 km3 => m3
    end do

!   Define a volume for each variable (with special case for isotopic variables!)
!   COMBINE variable order: all dissolved, then all particulate, then all isotopic
    !
    i = 0
    do j=1,nvar_diss ! dissolved variables
        do k = 1,nbasin
            i = i+1
            vol(i) = box_vol(k)
        end do
    end do
    do j=1,nvar_part ! particulate variables
        do k = 1,nbasin
            i = i+1
            vol(i) = box_vol(k)
        end do
    end do
    do j=1,nvar_isot ! isotopic variables do not need volume => will be multiplied by 1
        do k = 1,nbasin
            i = i+1
            vol(i) = 1
        end do
    end do

    ndeep            = 0
    nnodeep          = 0
    nsedi            = 0
    nnosedi          = 0
    nthermo          = 0
    nnothermo        = 0
    nsurface         = 0
    nnosurface       = 0
    nepicont         = 0
    nnoepicont       = 0
    npolar           = 0
    nnopolar         = 0
    nappcont         = 0
    nnoappcont       = 0
    nsediepicsurf    = 0
    nsediepicnosurf  = 0
    nsediintermed    = 0
    nsedideep        = 0
    nsurfnoappcont   = 0

    do i=1,nbasin-1

        read(34,*)oce_surf(i)
        read(35,*)indice_deep(i) ! =1 si oce profond sinon 0
        read(37,*)surf_sedi(i)
        read(38,*) depth_box(i), press_box(i)
        read(39,*)app_cont(i)
        read(40,*)indice_thermo(i)
        read(41,*)indice_surface(i)
        read(46,*)indice_epicont(i)
        read(47,*)indice_polar(i)

        if (surf_sedi(i)==0d0) then
            indice_sedi(i) = 0
        else
            indice_sedi(i) = 1
        end if

        if(indice_deep(i)==1) then
            ndeep = ndeep + 1
            jbox_deep(ndeep) = i
        else
            nnodeep = nnodeep + 1
            jbox_nodeep(nnodeep) = i
        end if
        if(indice_sedi(i)==1) then
            nsedi = nsedi + 1
            jbox_sedi(nsedi) = i
            if (indice_epicont(i).eq.1) then
                if (indice_surface(i).eq.1) then
                    nsediepicsurf = nsediepicsurf + 1
                    jbox_sediepicsurf(nsediepicsurf) = i
                else
                    nsediepicnosurf = nsediepicnosurf + 1
                    jbox_sediepicnosurf(nsediepicnosurf) = i
                end if
            elseif (indice_deep(i)==1) then
                nsedideep = nsedideep + 1
                jbox_sedideep(nsedideep) = i
            else ! all other sedimenting boxes are assumed to be intermediate depth (0-1000 m or 100-1000 m)
                nsediintermed = nsediintermed + 1
                jbox_sediintermed(nsediintermed) = i
            endif

        else
            nnosedi = nnosedi + 1
            jbox_nosedi(nnosedi) = i
        end if
        if(indice_thermo(i)==1) then
            nthermo = nthermo + 1
            jbox_thermo(nthermo) = i
        else
            nnothermo = nnothermo + 1
            jbox_nothermo(nnothermo) = i
        end if
        if(indice_surface(i)==1) then
            nsurface = nsurface + 1
            jbox_surface(nsurface) = i
            depth_top_box(i) = 0
        else
            nnosurface = nnosurface + 1
            jbox_nosurface(nnosurface) = i
            depth_top_box(i) = depth_box(i-1)
        end if
        if(indice_epicont(i)==1) then
            nepicont = nepicont + 1
            jbox_epicont(nepicont) = i
        else
            nnoepicont = nnoepicont + 1
            jbox_noepicont(nnoepicont) = i
        end if
        if(indice_polar(i)==1) then
            npolar = npolar + 1
            jbox_polar(npolar) = i
        else
            nnopolar = nnopolar + 1
            jbox_nopolar(nnopolar) = i
        end if
        if(app_cont(i)==1) then
            nappcont = nappcont + 1
            jbox_appcont(nappcont) = i
        else
            nnoappcont = nnoappcont + 1
            jbox_noappcont(nnoappcont) = i
            if(indice_surface(i)==1) then
                nsurfnoappcont = nsurfnoappcont + 1
                jbox_surfnoappcont(nsurfnoappcont) = i
            end if
        end if

    end do


    ! Consistency checks
    ! ------------------

    ! * for seafloor sedimentation scheme
    n = min(nsediepicsurf, nappcont)
    if ((nsediepicsurf/=nappcont) .or. (.not. all(jbox_sediepicsurf(1:n)==jbox_appcont(1:n)))) then
        print *, ''
        print *, 'ERROR: inconsistent boundary conditions.'
        print *, '  GEOCLIM seafloor sedimentation scheme (defined in seafloor_sedimentation.f90)'
        print *, '  requires that the sedimenting epicontinental surface basins MUST CORRESPOND to'
        print *, '  the basins with continental inputs.'
        print *, '  This condition is not satisfied by the current COMBINE input files'
        print *, '  (see files "indice_epicont.dat", "indice_surface.dat", "surf_sedi.dat" and'
        print *, '  "apport_ct.dat").'
        stop 18
    end if


    ! Trick to perform "vertical" loops for each series of "stacked" basins
    jbox_surface(nsurface+1) = nbasin

    ! One more "read" for atmosphere box
    read(34,*)oce_surf(nbasin)
    read(35,*)indice_deep(nbasin)
    read(37,*)surf_sedi(nbasin)
    read(38,*)depth_box(nbasin), press_box(nbasin)
    read(39,*)app_cont(nbasin)
    read(40,*)indice_thermo(nbasin)
    read(41,*)indice_surface(nbasin)
    read(46,*)indice_epicont(nbasin)
    read(47,*)indice_polar(nbasin)

    ! thickness of boxes
    box_thick = depth_box - depth_top_box

    ! area conversion: 1e9 km2 => m2
    oce_surf = 1d15*oce_surf
    surf_sedi = 1d15*surf_sedi

    oce_surf_tot=oce_surf(nbasin)


    ! Area fraction of coastal basins:
    ! --------------------------------

    ! Total epicontinental surf area
    area_appcont = 0d0
    do j0 = 1,nappcont
        j = jbox_appcont(j0)
        area_appcont = area_appcont + oce_surf(j) 
    end do

    ! area fraction
    appcont_fract = 0d0
    do j0 = 1,nappcont
        j = jbox_appcont(j0)
        appcont_fract(j) = oce_surf(j)/area_appcont
    end do


    ! Sedimentation capacity of bottom bassins + area fraction of coastal basins:
    ! ---------------------------------------------------------------------------

    ! Total sedimentation area
    area_sediepicsurf   = 0d0
    area_sediepicnosurf = 0d0
    area_sedideep       = 0d0
    area_sediintermed   = 0d0
    do j0 = 1,nsediepicsurf
        j = jbox_sediepicsurf(j0)
        area_sediepicsurf = area_sediepicsurf + surf_sedi(j) 
    end do
    do j0 = 1,nsediepicnosurf
        j = jbox_sediepicnosurf(j0)
        area_sediepicnosurf = area_sediepicnosurf + surf_sedi(j) 
    end do
    do j0 = 1,nsediintermed
        j = jbox_sediintermed(j0)
        area_sediintermed = area_sediintermed + surf_sedi(j) 
    end do
    do j0 = 1,nsedideep
        j = jbox_sedideep(j0)
        area_sedideep = area_sedideep + surf_sedi(j) 
    end do

    ! fraction of total area of bassin type (-), and sedimentation capacity (in m3/yr)
    sedbastype_fract = 0d0
    accumul_capacity = 0d0
    do j0 = 1,nsediepicsurf
        j = jbox_sediepicsurf(j0)
        sedbastype_fract(j) = surf_sedi(j)/area_sediepicsurf
        accumul_capacity(j) = sedbastype_fract(j) * ksed*area_sediepicsurf**1.5
    end do
    do j0 = 1,nsediepicnosurf
        j = jbox_sediepicnosurf(j0)
        sedbastype_fract(j) = surf_sedi(j)/area_sediepicnosurf
        accumul_capacity(j) = sedbastype_fract(j) * ksed*area_sediepicnosurf**1.5
    end do
    do j0 = 1,nsediintermed
        j = jbox_sediintermed(j0)
        sedbastype_fract(j) = surf_sedi(j)/area_sediintermed
        accumul_capacity(j) = sedbastype_fract(j) * ksed*area_sediintermed**1.5
    end do
    do j0 = 1,nsedideep
        j = jbox_sedideep(j0)
        sedbastype_fract(j) = surf_sedi(j)/area_sedideep
        accumul_capacity(j) = DEFAULT_FILLVAL ! infinity
    end do



    return
    end
