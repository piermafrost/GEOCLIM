module dynsoil_create_init_condition
    implicit none

    contains

    subroutine create_init_condition(which_initcond, loc_slopemissval, &
                                     loc_xlevs, loc_reg_thick, loc_x_P_surf, loc_tau_surf, loc_z_prof, loc_tau_prof, &
                                     loc_missingpoints)

        use constante, only: PI_n_CO2_atm
        use multidimensional_interpolation, only: climate_interpolation

        include 'combine_foam.inc'

        character(len=*), intent(in):: which_initcond
        double precision, intent(in):: loc_slopemissval
        double precision, dimension(nDSlev), intent(out):: loc_xlevs
        double precision, dimension(nlitho,npixel), intent(out):: loc_reg_thick, loc_x_P_surf, loc_tau_surf
        double precision, dimension(nDSlev,nlitho,npixel), intent(out):: loc_z_prof, loc_tau_prof
        logical, dimension(npixel), intent(out):: loc_missingpoints

        loc_missingpoints = (areaclimber==0 .or. slope==loc_slopemissval)

        loc_xlevs = (/( dble(k)/dble(nDSlev) , k=nDSlev,1,-1 )/)

        select case (which_initcond)

            case ("null")
                do j = 1,npixel
                    if (.not. loc_missingpoints(j)) then
                            loc_reg_thick(:,j)  = 0
                            loc_x_P_surf(:,j)   = 1
                            loc_tau_surf(:,j)   = 0
                            loc_z_prof(1,:,j)   = 0
                            loc_tau_prof(1,:,j) = 0
                        end if
                end do

            case ("eq")
                call get_cont_pixel() ! => set global variable 'list_cont_pixel' needed for climate interpolation
                p=var_diss(7,nbasin)/PI_n_CO2_atm
                call climate_interpolation(co2climber, clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5,    &
                                           p, cpvec, list_cont_pixel=list_cont_pixel, ncontpxl=ncontpxl,                        &
                                           temp_array=Tairclimber, runf_array=Runclimber, interp_temp=Tclim, interp_runf=runclim)
                call equilibrium_values(loc_missingpoints, Tclim, Runclim, slope, loc_xlevs, &
                                        loc_reg_thick, loc_x_P_surf, loc_tau_surf, loc_z_prof, loc_tau_prof)

            case default
                print *
                print *, 'ERROR: invalid DynSoil startup case "'//trim(which_initcond)//'"'
                print *, 'Only "eq" or "null" are allowed.'
                stop 47

        end select

        print *
        print *, 'Created "'//trim(which_initcond)//'" DynSoil initial conditions.'
        print *

    end subroutine


    !=============================================================================================!


    subroutine equilibrium_values(missingpoints, temp, runoff, slope, x, reg_thick, x_P_surf, tau_surf, z_prof, tau_prof)

        use dynsoil_empirical_laws, only: reg_prod_opt, erosion, dissolution_constant
        use dynsoil_physical_parameters, only: nDSlev, nlon, nlat, nlitho, h0, sigma, epsl

        integer, parameter:: npixel = nlon*nlat
        logical, dimension(npixel), intent(in):: missingpoints
        double precision, dimension(npixel), intent(in):: temp, runoff, slope
        double precision, dimension(nDSlev), intent(in):: x
        double precision, dimension(nlitho,npixel), intent(out):: reg_thick, x_P_surf, tau_surf
        double precision, dimension(nDSlev,nlitho,npixel), intent(out):: z_prof, tau_prof

        double precision:: RPopt, E, Kmain
        integer:: i, j, k

        do i = 1,npixel
            if (.not. missingpoints(i)) then

                if (runoff(i) > 0) then

                    do j = 1,nlitho-1 ! skip carbonate

                        ! optimal regolith production rate:
                        RPopt  = reg_prod_opt(temp(i),runoff(i),j)
                        ! regolith erosion rate:
                        E = erosion(temp(i),runoff(i),slope(i),j)
                        ! dissolution constant:
                        Kmain = dissolution_constant(temp(i),runoff(i),j)

                        ! TOP OF REGOLITH
                        reg_thick(j,i) = h0 * log(RPopt/E)
                        if (reg_thick(j,i)<0) reg_thick(j,i) = 0
                        tau_surf(j,i) = reg_thick(j,i)/E
                        x_P_surf(j,i) = exp( -1*Kmain * ((reg_thick(j,i)/E)**(sigma(j)+1)) / (sigma(j)+1) )

                        ! INNER REGOLITH (z and tau): !
                        z_prof(1,j,i)   = 0
                        tau_prof(1,j,i) = 0
                        k = 1
                        do while ( k < nDSlev )
                            ! test if next point is outside the regolith:
                            if ( x(k+1) < x_P_surf(j,i)+epsl ) then
                                k = nDSlev
                            else
                                tau_prof(k+1,j,i) = (-1*(sigma(j)+1)*log(x(k+1))/Kmain)**(1/(sigma(j)+1))
                                z_prof(k+1,j,i) = E*tau_prof(k+1,j,i)
                                k = k + 1
                            end if
                        end do

                    end do

                else ! if null runoff

                    reg_thick(:,i)  = 0
                    x_P_surf(:,i)   = 1
                    tau_surf(:,i)   = 0
                    z_prof(1,:,i)   = 0
                    tau_prof(1,:,i) = 0

                end if

            end if
        end do

    end subroutine

end module
