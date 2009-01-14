
 
subroutine ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,annual_time&
                    ,writing_dail,writing_mont,history_time,the_end)
  
  use ed_state_vars,only:edgrid_g

  use grid_coms, only: ngrids,nzg  ! INTENT(IN)

  use ed_node_coms,only : mynum,nnodetot



  use misc_coms, only: dtlsm, current_time, &
       idoutput, &
       imoutput, &
       iyoutput, &
       isoutput, &
       ifoutput, &
       iprintpolys, &
       frqsum

  implicit none
    
  logical, intent(in)  :: the_end,analysis_time,dail_analy_time
  logical, intent(in)  :: writing_dail,writing_mont
  logical, intent(in) :: mont_analy_time,history_time,new_day,annual_time
  integer :: ifm

  if(analysis_time .or. history_time .or. (new_day .and. (writing_dail .or. writing_mont))) then
     do ifm=1,ngrids
        call normalize_averaged_vars_ar(edgrid_g(ifm),frqsum,dtlsm)
     enddo

     !  Perform averaging and data preparation
     call spatial_averages
     
     if (writing_dail .or. writing_mont) then
        do ifm=1,ngrids
           call integrate_ed_daily_output_flux(edgrid_g(ifm))
        end do
     end if
  endif
  
  
  if (analysis_time) then
     
    
     !  Write out analysis fields - mostly polygon averages
     if (ifoutput.eq.3) then
        call h5_output('INST')
     endif
     
     ! If printpolys is on then print this info to
     ! the screen
     
     if (iprintpolys.eq.1) then
        do ifm=1,ngrids
           call print_array(ifm,edgrid_g(ifm))
        enddo
     endif

  endif

  ! Daily analysis output and monthly integration
  if (new_day .and. (writing_dail .or. writing_mont)) then

     do ifm=1,ngrids
        call normalize_ed_daily_output_vars(edgrid_g(ifm))
        if (writing_mont) call integrate_ed_monthly_output_vars(edgrid_g(ifm))
     end do

     if (dail_analy_time) call h5_output('DAIL')

     do ifm=1,ngrids
       call zero_ed_daily_output_vars(edgrid_g(ifm))
     end do

  end if

  ! Monthly analysis output
  if (mont_analy_time) then

     do ifm=1,ngrids
        call normalize_ed_monthly_output_vars(edgrid_g(ifm))
     end do

     call h5_output('MONT')

     do ifm=1,ngrids
        call zero_ed_monthly_output_vars(edgrid_g(ifm))
     end do
  end if

  if (annual_time) then

     call h5_output('YEAR')

  endif

  ! History files should only be output at a frequency which
  ! divides by frqanl, thus the integrated fast-time variables
  ! are valid, but representative of the last frqanl period, not
  ! the last frqhist period.


  if(history_time) then
     
     call h5_output('HIST')
     
  endif



  return
end subroutine ed_output
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine spatial_averages
  
  ! ----------------------------------------------------------------------------
  ! The following subroutine performs several spatial averaging functions
  ! and temporal integrations.  Specifically, it area averages patch level
  ! quantities to the site, and site level quantities to the polygon
  ! -----------------------------------------------------------------------------

  use ed_state_vars,only:edtype,polygontype,sitetype,patchtype,edgrid_g
  use grid_coms, only : ngrids,nzg,nzs
  use canopy_radiation_coms, only: lai_min
  use consts_coms, only : alvl
  use misc_coms, only: frqsum

  implicit none
  
  type(edtype),pointer :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer :: csite
  type(patchtype),pointer :: cpatch
  integer :: igr,ipy,isi,ipa
  integer :: k
  integer :: lai_index
  real :: lai_sum,site_area_i,poly_area_i
  real :: frqsumi
  real :: snow_min = 0.0000001
  real,dimension(3) :: area_sum

  frqsumi = 1.0 / frqsum
  do igr=1,ngrids
     cgrid => edgrid_g(igr)

     cgrid%avg_lai_ebalvars = 0.0

     do ipy=1,cgrid%npolygons
        cpoly => cgrid%polygon(ipy)

        area_sum=0.0

        cgrid%avg_balive(ipy)      = 0.0
        cgrid%avg_bdead(ipy)       = 0.0

        cgrid%lai(ipy)             = 0.0

        cgrid%avg_gpp(ipy)         = 0.0
        cgrid%avg_leaf_resp(ipy)   = 0.0
        cgrid%avg_root_resp(ipy)   = 0.0
        cgrid%avg_plant_resp(ipy)  = 0.0
        cgrid%avg_htroph_resp(ipy) = 0.0

        cgrid%max_veg_temp(ipy)    = -10.0
        cgrid%min_veg_temp(ipy)    = 390.0
        cgrid%max_soil_temp(ipy)    = -10.0
        cgrid%min_soil_temp(ipy)    = 390.0

        poly_area_i = 1./sum(cpoly%area)

        do isi=1,cpoly%nsites
           csite => cpoly%site(isi)
           
           if (csite%npatches>0) then
              site_area_i=1./sum(csite%area)

           cpoly%lai(isi)                = sum(csite%lai                * csite%area ) * site_area_i


           ! Average Fast Time Flux Dynamics Over Sites
           cpoly%avg_vapor_vc(isi)       = sum(csite%avg_vapor_vc       * csite%area ) * site_area_i
           cpoly%avg_dew_cg(isi)         = sum(csite%avg_dew_cg         * csite%area ) * site_area_i
           cpoly%avg_vapor_gc(isi)       = sum(csite%avg_vapor_gc       * csite%area ) * site_area_i
           cpoly%avg_wshed_vg(isi)       = sum(csite%avg_wshed_vg       * csite%area ) * site_area_i
           cpoly%avg_vapor_ac(isi)       = sum(csite%avg_vapor_ac       * csite%area ) * site_area_i
           cpoly%avg_transp(isi)         = sum(csite%avg_transp         * csite%area ) * site_area_i
           cpoly%avg_evap(isi)           = sum(csite%avg_evap           * csite%area ) * site_area_i
           cpoly%aux(isi)                = sum(csite%aux                * csite%area ) * site_area_i
           cpoly%avg_sensible_vc(isi)    = sum(csite%avg_sensible_vc    * csite%area ) * site_area_i
           cpoly%avg_sensible_2cas(isi)  = sum(csite%avg_sensible_2cas  * csite%area ) * site_area_i
           cpoly%avg_qwshed_vg(isi)      = sum(csite%avg_qwshed_vg      * csite%area ) * site_area_i
           cpoly%avg_sensible_gc(isi)    = sum(csite%avg_sensible_gc    * csite%area ) * site_area_i
           cpoly%avg_sensible_ac(isi)    = sum(csite%avg_sensible_ac    * csite%area ) * site_area_i
           cpoly%avg_sensible_tot(isi)   = sum(csite%avg_sensible_tot   * csite%area ) * site_area_i

           !! for NACP intercomparision (MCD)
           cpoly%avg_fsc(isi)            = sum(csite%fast_soil_C        * csite%area ) * site_area_i
           cpoly%avg_ssc(isi)            = sum(csite%slow_soil_C        * csite%area ) * site_area_i
           cpoly%avg_stsc(isi)           = sum(csite%structural_soil_C  * csite%area ) * site_area_i
           cpoly%avg_co2can(isi)         = sum(csite%can_co2            * csite%area ) * site_area_i
           
           cpoly%avg_snowdepth(isi)   = 0.0
           cpoly%avg_snowmass(isi)    = 0.0
           cpoly%avg_snowtempk(isi)   = 0.0
           cpoly%avg_snowfracliq(isi) = 0.0
           do k=1,csite%npatches
              if(csite%nlev_sfcwater(k) > 0 .and. csite%sfcwater_mass(1,k) > snow_min) then
                 !! if snow is present, sum mass and depth over layers
                 cpoly%avg_snowdepth(isi)   = cpoly%avg_snowdepth(isi)   + &
                      &sum(csite%sfcwater_depth(1:csite%nlev_sfcwater(k)  ,k)) * csite%area(k)
                 cpoly%avg_snowmass(isi)    = cpoly%avg_snowmass(isi)    + &
                      &sum(csite%sfcwater_mass(1:csite%nlev_sfcwater(k)   ,k)) * csite%area(k)
                 !! average temp and frac liq weighted by **mass**
                 cpoly%avg_snowtempk(isi)   = cpoly%avg_snowtempk(isi)   + &
                      &sum(csite%sfcwater_tempk(1:csite%nlev_sfcwater(k),k)*csite%sfcwater_mass(1:csite%nlev_sfcwater(k),k))&
                      &/sum(csite%sfcwater_mass(1:csite%nlev_sfcwater(k)   ,k))* csite%area(k)
                 cpoly%avg_snowfracliq(isi) = cpoly%avg_snowfracliq(isi) + &
                      & sum(csite%sfcwater_fracliq(1:csite%nlev_sfcwater(k),k)*csite%sfcwater_mass(1:csite%nlev_sfcwater(k),k)) &
                      &/sum(csite%sfcwater_mass(1:csite%nlev_sfcwater(k)   ,k)) * csite%area(k)

              else
                 !! assume a thin layer in equilibrium with the top soil layer
                 cpoly%avg_snowtempk(isi)   = cpoly%avg_snowtempk(isi) + csite%soil_tempk(nzg,k) * csite%area(k)
                 cpoly%avg_snowtempk(isi)   = cpoly%avg_snowtempk(isi) + csite%soil_fracliq(nzg,k) * csite%area(k)
              endif
           enddo
           cpoly%avg_snowfracliq(isi) = cpoly%avg_snowfracliq(isi) * site_area_i
           cpoly%avg_snowdepth(isi)   = cpoly%avg_snowdepth(isi)   * site_area_i
           cpoly%avg_snowmass(isi)    = cpoly%avg_snowmass(isi)    * site_area_i
           cpoly%avg_snowtempk(isi)   = cpoly%avg_snowtempk(isi)   * site_area_i


           !!--------------------

           do k=cpoly%lsl(isi),nzg

              cpoly%avg_sensible_gg(k,isi)  = sum(csite%avg_sensible_gg(k,:)  * csite%area ) * site_area_i
              cpoly%avg_smoist_gg(k,isi)    = sum(csite%avg_smoist_gg(k,:)    * csite%area ) * site_area_i
              cpoly%avg_smoist_gc(k,isi)    = sum(csite%avg_smoist_gc(k,:)    * csite%area ) * site_area_i
              cpoly%aux_s(k,isi)            = sum(csite%aux_s(k,:)            * csite%area ) * site_area_i

              cpoly%avg_soil_energy(k,isi)  = sum(csite%soil_energy(k,:)      * csite%area ) * site_area_i
              cpoly%avg_soil_water(k,isi)   = real(sum(csite%soil_water(k,:)       * dble(csite%area) )) * site_area_i
              cpoly%avg_soil_temp(k,isi)    = sum(csite%soil_tempk(k,:)       * csite%area ) * site_area_i
              cpoly%avg_soil_fracliq(k,isi) = sum(csite%soil_fracliq(k,:)     * csite%area ) * site_area_i

           enddo
           
           ! Average over patches
           
           do ipa=1,csite%npatches
              cpatch => csite%patch(ipa)
              
              if (cpatch%ncohorts>0) then
                 
                 
                 lai_sum = max(lai_min,sum(cpatch%lai, cpatch%lai > lai_min))
                 csite%avg_veg_energy(ipa)  = sum(cpatch%veg_energy * cpatch%lai,cpatch%lai > lai_min)   / lai_sum
                 csite%avg_veg_temp(ipa)  = sum(cpatch%veg_temp * cpatch%lai,cpatch%lai > lai_min)   / lai_sum
                 csite%avg_veg_water(ipa) = sum(cpatch%veg_water * cpatch%lai,cpatch%lai > lai_min)  / lai_sum
                 
                 if (lai_sum > lai_min) then
                    csite%laiarea(ipa) = csite%area(ipa)
                 else
                    csite%laiarea(ipa) = 0.0
                 end if

                 cgrid%avg_gpp(ipy)       = cgrid%avg_gpp(ipy)       + &
                      csite%area(ipa)*cpoly%area(isi)*sum(cpatch%mean_gpp)
                 
                 cgrid%avg_leaf_resp(ipy) = cgrid%avg_leaf_resp(ipy) + &
                      csite%area(ipa)*cpoly%area(isi)*sum(cpatch%mean_leaf_resp)
                 cgrid%avg_root_resp(ipy) = cgrid%avg_root_resp(ipy) + &
                      csite%area(ipa)*cpoly%area(isi)*sum(cpatch%mean_root_resp)
                 cgrid%avg_balive(ipy)    = cgrid%avg_balive(ipy) + &
                      csite%area(ipa)*cpoly%area(isi)*sum(cpatch%balive*cpatch%nplant)
                 cgrid%avg_bdead(ipy)     = cgrid%avg_bdead(ipy) + &
                      csite%area(ipa)*cpoly%area(isi)*sum(cpatch%bdead*cpatch%nplant)

                 ! Check the extremes
                 
                 if (maxval(cpatch%veg_temp) > cgrid%max_veg_temp(ipy)) then
                    cgrid%max_veg_temp(ipy) = maxval(cpatch%veg_temp)
                 endif

                 if (minval(cpatch%veg_temp) < cgrid%min_veg_temp(ipy)) then
                    cgrid%min_veg_temp(ipy) = minval(cpatch%veg_temp)
                 endif

                 

                 
              else
                 ! Set veg-temp to air temp
                 csite%avg_veg_energy(ipa)  = 0.0          
                 csite%avg_veg_temp(ipa)  = csite%can_temp(ipa)
                 csite%avg_veg_water(ipa) = 0.0

              endif

              cgrid%avg_plant_resp(ipy)  = cgrid%avg_plant_resp(ipy)  + &
                   csite%area(ipa)*cpoly%area(isi)*csite%co2budget_plresp(ipa)*frqsumi
              cgrid%avg_htroph_resp(ipy) = cgrid%avg_htroph_resp(ipy) + &
                   csite%area(ipa)*cpoly%area(isi)*csite%co2budget_rh(ipa)    *frqsumi
              

              lai_index = min(3,max(1, floor(csite%lai(ipa)/2.0) + 1)  )
              
              area_sum(lai_index) = area_sum(lai_index) + csite%area(ipa)


              ! Net radiation
              cgrid%avg_lai_ebalvars(lai_index,1,ipy) = &
                   cgrid%avg_lai_ebalvars(lai_index,1,ipy) + csite%avg_netrad(ipa)*csite%area(ipa)

              ! Latent heat flux
              cgrid%avg_lai_ebalvars(lai_index,2,ipy) = &
                   cgrid%avg_lai_ebalvars(lai_index,2,ipy) + csite%avg_vapor_ac(ipa)*csite%area(ipa)

              ! Sensible heat flux
              cgrid%avg_lai_ebalvars(lai_index,3,ipy) = &
                   cgrid%avg_lai_ebalvars(lai_index,3,ipy) + csite%avg_sensible_ac(ipa)*csite%area(ipa)

              ! Canopy temperature
              cgrid%avg_lai_ebalvars(lai_index,4,ipy) = &
                   cgrid%avg_lai_ebalvars(lai_index,4,ipy) + csite%can_temp(ipa)*csite%area(ipa)


              
              do k=1,nzg
                 if ( csite%soil_tempk(k,ipa) > cgrid%max_soil_temp(ipy)) then
                    cgrid%max_soil_temp(ipy) = csite%soil_tempk(k,ipa)
                 endif
                 
                 if ( csite%soil_tempk(k,ipa) < cgrid%min_soil_temp(ipy)) then
                    cgrid%min_soil_temp(ipy) = csite%soil_tempk(k,ipa)
                 endif
              enddo
                 

           enddo
           
           csite%laiarea = csite%laiarea / max(sum(csite%laiarea),1.0)

           cpoly%avg_veg_energy(isi) = sum(csite%avg_veg_energy   * csite%laiarea)
           cpoly%avg_veg_temp(isi)   = sum(csite%avg_veg_temp     * csite%laiarea)
           cpoly%avg_veg_water(isi)  = sum(csite%avg_veg_water    * csite%laiarea)
           cpoly%avg_can_temp(isi)   = sum(csite%can_temp         * csite%area)
           cpoly%avg_can_shv(isi)    = sum(csite%can_shv          * csite%area)


        else
           call fatal_error('No patches in this site, impossible','spatial_averages','edio.f90')
        endif

        

     enddo
     
     ! Normalize the lai specific quantities
     if (area_sum(1)>0.) then
        cgrid%avg_lai_ebalvars(1,1,ipy) = cgrid%avg_lai_ebalvars(1,1,ipy)/area_sum(1)
        cgrid%avg_lai_ebalvars(1,2,ipy) = cgrid%avg_lai_ebalvars(1,2,ipy)/area_sum(1)
        cgrid%avg_lai_ebalvars(1,3,ipy) = cgrid%avg_lai_ebalvars(1,3,ipy)/area_sum(1)  
        cgrid%avg_lai_ebalvars(1,4,ipy) = cgrid%avg_lai_ebalvars(1,4,ipy)/area_sum(1)
     else
        cgrid%avg_lai_ebalvars(1,:,ipy) = -9999.0
     endif
     if (area_sum(2)>0.) then
        cgrid%avg_lai_ebalvars(2,1,ipy) = cgrid%avg_lai_ebalvars(2,1,ipy)/area_sum(2)
        cgrid%avg_lai_ebalvars(2,2,ipy) = cgrid%avg_lai_ebalvars(2,2,ipy)/area_sum(2)
        cgrid%avg_lai_ebalvars(2,3,ipy) = cgrid%avg_lai_ebalvars(2,3,ipy)/area_sum(2)
        cgrid%avg_lai_ebalvars(2,4,ipy) = cgrid%avg_lai_ebalvars(2,4,ipy)/area_sum(2)
     else
        cgrid%avg_lai_ebalvars(2,:,ipy) = -9999.0
     endif
     if (area_sum(3)>0.) then
        cgrid%avg_lai_ebalvars(3,1,ipy) = cgrid%avg_lai_ebalvars(3,1,ipy)/area_sum(3)
        cgrid%avg_lai_ebalvars(3,2,ipy) = cgrid%avg_lai_ebalvars(3,2,ipy)/area_sum(3)
        cgrid%avg_lai_ebalvars(3,3,ipy) = cgrid%avg_lai_ebalvars(3,3,ipy)/area_sum(3)
        cgrid%avg_lai_ebalvars(3,4,ipy) = cgrid%avg_lai_ebalvars(3,4,ipy)/area_sum(3)
     else
        cgrid%avg_lai_ebalvars(3,:,ipy) = -9999.0
     endif


        cgrid%lai(ipy)                = sum(cpoly%lai                * cpoly%area ) * poly_area_i
        
        ! Average Fast Time Flux Dynamics Over Polygons
        cgrid%avg_vapor_vc(ipy)       = sum(cpoly%avg_vapor_vc       * cpoly%area ) * poly_area_i
        cgrid%avg_dew_cg(ipy)         = sum(cpoly%avg_dew_cg         * cpoly%area ) * poly_area_i
        cgrid%avg_vapor_gc(ipy)       = sum(cpoly%avg_vapor_gc       * cpoly%area ) * poly_area_i
        cgrid%avg_wshed_vg(ipy)       = sum(cpoly%avg_wshed_vg       * cpoly%area ) * poly_area_i
        cgrid%avg_co2can(ipy)         = sum(cpoly%avg_co2can         * cpoly%area ) * poly_area_i
        cgrid%avg_fsc(ipy)            = sum(cpoly%avg_fsc            * cpoly%area ) * poly_area_i
        cgrid%avg_stsc(ipy)           = sum(cpoly%avg_stsc           * cpoly%area ) * poly_area_i
        cgrid%avg_ssc(ipy)            = sum(cpoly%avg_ssc            * cpoly%area ) * poly_area_i
        cgrid%avg_runoff_heat(ipy)    = sum(cpoly%avg_runoff_heat    * cpoly%area ) * poly_area_i
        cgrid%avg_runoff(ipy)         = sum(cpoly%avg_runoff         * cpoly%area ) * poly_area_i
        cgrid%avg_snowmass(ipy)       = sum(cpoly%avg_snowmass       * cpoly%area ) * poly_area_i
        cgrid%avg_snowtempk(ipy)      = sum(cpoly%avg_snowtempk      * cpoly%area ) * poly_area_i
        cgrid%avg_snowfracliq(ipy)    = sum(cpoly%avg_snowfracliq    * cpoly%area ) * poly_area_i
        cgrid%avg_snowdepth(ipy)      = sum(cpoly%avg_snowdepth      * cpoly%area ) * poly_area_i

        ! I REMOVED THE FOLLOWING SCALING BY ALVL. THIS VARIABLE IS ALREADY IN ENERGY UNITS, SEE
        ! RK4_DERIVS.F90. RGK-6-30-08

        cgrid%avg_vapor_ac(ipy)       = sum(cpoly%avg_vapor_ac       * cpoly%area ) * poly_area_i !* alvl
        cgrid%avg_transp(ipy)         = sum(cpoly%avg_transp         * cpoly%area ) * poly_area_i
        cgrid%avg_evap(ipy)           = sum(cpoly%avg_evap           * cpoly%area ) * poly_area_i
        cgrid%aux(ipy)                = sum(cpoly%aux                * cpoly%area ) * poly_area_i
        cgrid%avg_sensible_vc(ipy)    = sum(cpoly%avg_sensible_vc    * cpoly%area ) * poly_area_i
        cgrid%avg_sensible_2cas(ipy)  = sum(cpoly%avg_sensible_2cas  * cpoly%area ) * poly_area_i
        cgrid%avg_qwshed_vg(ipy)      = sum(cpoly%avg_qwshed_vg      * cpoly%area ) * poly_area_i
        cgrid%avg_sensible_gc(ipy)    = sum(cpoly%avg_sensible_gc    * cpoly%area ) * poly_area_i
        cgrid%avg_sensible_ac(ipy)    = sum(cpoly%avg_sensible_ac    * cpoly%area ) * poly_area_i
        cgrid%avg_sensible_tot(ipy)   = sum(cpoly%avg_sensible_tot   * cpoly%area ) * poly_area_i


        cgrid%avg_veg_energy(ipy)     = sum(cpoly%avg_veg_energy     * cpoly%area ) * poly_area_i
        cgrid%avg_veg_temp(ipy)       = sum(cpoly%avg_veg_temp       * cpoly%area ) * poly_area_i
        cgrid%avg_veg_water(ipy)      = sum(cpoly%avg_veg_water      * cpoly%area ) * poly_area_i
        cgrid%avg_can_temp(ipy)       = sum(cpoly%avg_can_temp       * cpoly%area ) * poly_area_i
        cgrid%avg_can_shv(ipy)        = sum(cpoly%avg_can_shv        * cpoly%area ) * poly_area_i

        do k=cgrid%lsl(ipy),nzg

           cgrid%avg_sensible_gg(k,ipy)  = sum(cpoly%avg_sensible_gg(k,:)  * cpoly%area ) * poly_area_i
           cgrid%avg_smoist_gg(k,ipy)    = sum(cpoly%avg_smoist_gg(k,:)    * cpoly%area ) * poly_area_i
           cgrid%avg_smoist_gc(k,ipy)    = sum(cpoly%avg_smoist_gc(k,:)    * cpoly%area ) * poly_area_i
           cgrid%aux_s(k,ipy)            = sum(cpoly%aux_s(k,:)            * cpoly%area ) * poly_area_i
           cgrid%avg_soil_energy(k,ipy)  = sum(cpoly%avg_soil_energy(k,:)  * cpoly%area ) * poly_area_i
           cgrid%avg_soil_water(k,ipy)   = sum(cpoly%avg_soil_water(k,:)   * cpoly%area ) * poly_area_i
           cgrid%avg_soil_temp(k,ipy)    = sum(cpoly%avg_soil_temp(k,:)    * cpoly%area ) * poly_area_i
           cgrid%avg_soil_fracliq(k,ipy) = sum(cpoly%avg_soil_fracliq(k,:) * cpoly%area ) * poly_area_i

        enddo
        
     enddo

  enddo

  return
end subroutine spatial_averages

! ==============================

subroutine get2d(m1,m2,temp2,in_ptr)
  
  implicit none
  integer :: m1,m2,i,j
  real,dimension(m1,m2) :: temp2
  real,dimension(m1,m2) :: in_ptr
  do i = 1,m1
     do j = 1,m2
        temp2(i,j) = in_ptr(i,j)
     enddo
  enddo
  return
end subroutine get2d

! =============================

subroutine get3d(m1,m2,m3,temp3,in_ptr)
  
  implicit none
  integer :: m1,m2,m3,i,j,k
  real,dimension(m1,m2,m3) :: temp3
  real,dimension(m1,m2,m3) :: in_ptr
  do i = 1,m1
     do j = 1,m2
        do k = 1,m3
           temp3(i,j,k) = in_ptr(i,j,k)
        enddo
     enddo
  enddo
  return
end subroutine get3d

! =================================================

subroutine print_array(ifm,cgrid)
  
  !------------------------------------------------------
  ! PRINT OUT FIELDS OF INTEREST
  !
  ! This subroutine prints patch, cohort, polygon or site level
  ! data, upscales that data to the site level and stores
  ! the data in its spatial array coordinate.
  ! The data is then printed to the screen, based on a
  ! specified window of data.
  ! Note, this may be printing windows on various nodes,
  ! or this may be called from a master process.  Be
  ! conscious of this; as it will dictate what part of the
  ! domain you are accessing variables from, and whether
  ! or not the variable of interest is stored in memory
  ! at that time.  For instance, many variables are stored
  ! only on the slave nodes, and need not be passed back
  ! to the master.  Likewise, many slave node data will
  ! accumulate after each lsm call, until they are passed back
  ! to the master, where they are normalized. These variables
  ! will be immediately zeroed on the slaves after being
  ! sent to the master.
  ! Dont forget to adjust the number precision on the 
  ! format string at the end.  The X.Xf
  !--------------------------------------------------------
  
  use ed_node_coms,only: mynum,nnodetot,sendnum,recvnum,master_num,machs
  use ed_state_vars,only: edtype,polygontype
  use misc_coms, only: &
            printvars,  &
            ipmax,      &
            ipmin,      &
            pfmtstr,    &
            iprintpolys
  
  use var_tables_array,only:vt_info,num_var


  implicit none

  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE) :: status
  integer  :: ifm,nv,np,i,ip
  integer  :: ping
  integer  :: npolys
  integer  :: g_idmin,g_idmax,l_idmin,l_idmax,ierr
  integer  :: node_idmin,node_idmax
  integer  :: mast_idmin,mast_idmax
  integer  :: g_id,g_ln,nm
  integer  :: ncols,row,maxrows,col
  integer,parameter :: maxcols = 10

  real,pointer,dimension(:) :: pvar_l
  real,pointer,dimension(:) :: pvar_g
  
  character(len=30)     :: fmtstr
  character(len=32)     :: pvar_name
  
  ! Linked structures
  type(edtype),target     :: cgrid
  
  logical :: pvartrue
  logical :: ptr_recv
  logical :: ptr_send

  ping = 8675309

  
  npolys = ipmax - ipmin + 1
  
  ! Adjust the format string according to the chosen
  ! Variables
  ! ------------------------------------------------

  ! CHeck the window size
  ! ---------------------
  
  if (ipmax.gt.cgrid%npolygons_global) then
     print*,"========================================="
     print*,"You have specified a print index"
     print*,"greater than the total number of polygons"
     print*,"You must reduce this number. Stopping"
     stop
  end if


  ! Allocate the print and scratch vector
  allocate(pvar_l(npolys))

  if (mynum .eq. nnodetot .or. nnodetot .eq. 1) allocate(pvar_g(npolys))

  
  ! Loop through the printvar entries from the namelist

  ip = 0
  
  count_pvars: do

     ip = ip+1
     pvar_name = printvars(ip)
     if (len_trim(pvar_name).eq.32) then
        
        exit count_pvars
     endif

!     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     pvartrue = .false.
     do nv = 1,num_var(ifm)

        if(trim(vt_info(nv,ifm)%name) .eq. trim(pvar_name)) then
           pvartrue = .true.
           
           ! If this is root, then collect the sends, keep reading to find out what it is
           ! receiving
           if (nnodetot.gt.1) then
              
              if (mynum .eq. nnodetot) then
                 
                 pvar_g = -99.9
                 ! Loop through the variable table to match the print variable
                 print*,""
                 print*," ============ ",trim(pvar_name)," =================="
                 print*,""

                 do nm = 1,nnodetot-1

                    call MPI_Recv(ptr_recv,1,MPI_LOGICAL,machs(nm),31,MPI_COMM_WORLD,status,ierr)

                    if (ptr_recv) then
                       call MPI_Recv(mast_idmin,1,MPI_INTEGER,machs(nm),32,MPI_COMM_WORLD,status,ierr)
                       call MPI_Recv(mast_idmax,1,MPI_INTEGER,machs(nm),33,MPI_COMM_WORLD,status,ierr)
                       call MPI_Recv(pvar_g(mast_idmin:mast_idmax),mast_idmax-mast_idmin+1,MPI_REAL,&
                            machs(nm),34,MPI_COMM_WORLD,status,ierr)
                    end if
                 enddo
              endif
           else
              
              ! Loop through the variable table to match the print variable
              print*,""
              print*," ============ ",trim(pvar_name)," =================="
              print*,""
              pvar_g = -99.9
           endif
        
              
           ! The namelist print entry has been matched with the var_table
           ! entry.  Now lets cycle through our machines and determine
           ! if those machines hold data that should be printed.  If so, 
           ! then send that data to the master (machine 0).

           ! This is scratch space that all machines will use

           pvar_l =-99.9
           
           ! Set the blocking recieve to allow ordering, start with machine 1
           if (mynum /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,212,MPI_COMM_WORLD,status,ierr)
           
           ! Cycle through this node's pointers for the current variable.  If the index
           ! falls within the printable range. Save that pointer to a local array.
           ! Once all the pointers have been cycled, send the local array to the 
           ! master to populate a global array and print
           
           ptr_send = .false.
           node_idmin = -1
           node_idmax = -1

           do np = 1,vt_info(nv,ifm)%nptrs
              g_id = vt_info(nv,ifm)%vt_vector(np)%globid+1
              g_ln = vt_info(nv,ifm)%vt_vector(np)%varlen
              
              ! Determine if any of this segment falls within the 
              ! range desired for output, globid+1 is the global index

              if (g_id .le. ipmax .and. g_id+g_ln-1 .ge. ipmin ) then

                 ! OK, this segment is good, set the send flag
                 ptr_send = .true.

                 ! These are the indices of the data in the current segment to use
                 ! and the indices in the global array they will be sent to
                 if (g_id >= ipmin) then
                    l_idmin = 1
                    g_idmin = g_id - ipmin + 1
                 else
                    l_idmin = ipmin - g_id + 1
                    g_idmin = 1
                 endif

                 if (g_id+g_ln-1 < ipmax) then
                    l_idmax = g_ln
                    g_idmax = g_id + g_ln - ipmin
                 else
                    l_idmax = ipmax - g_id + 1
                    g_idmax = ipmax - ipmin + 1
                 endif
                 
                 ! These should be the same size, if not...
                 if (l_idmax - l_idmin .ne. g_idmax - g_idmin ) then
                    print*,"NOT THE SAME LENGTHS"
                    print*,l_idmax - l_idmin,g_idmax - g_idmin
                    print*,l_idmax,l_idmin,g_idmax,g_idmin
                    stop
                 endif

                 ! Shift the global dataset so that it is applied to the
                 ! first index

                 call fillvar_l(pvar_l,vt_info(nv,ifm)%vt_vector(np)%var_rp,npolys,g_ln,g_idmin,l_idmin,l_idmax)

                 ! Determine the minimum and maximum indices that will be sent
                 if (g_idmin < node_idmin .or. node_idmin.eq.-1) node_idmin = g_idmin
                 if (g_idmax > node_idmax .or. node_idmax.eq.-1) node_idmax = g_idmax

              end if
              
           enddo
           
           if (nnodetot.gt.1) then
              
              ! The local array for this machine has been created. Send it off to the master

              if (mynum /= nnodetot) then

                 call MPI_Send(ptr_send,1,MPI_LOGICAL,machs(nnodetot),31,MPI_COMM_WORLD,ierr)
                 if (ptr_send) then
                    call MPI_Send(node_idmin,1,MPI_INTEGER,machs(nnodetot),32,MPI_COMM_WORLD,ierr)
                    call MPI_Send(node_idmax,1,MPI_INTEGER,machs(nnodetot),33,MPI_COMM_WORLD,ierr)
                    call MPI_Send(pvar_l(node_idmin:node_idmax),node_idmax-node_idmin+1, &
                         MPI_REAL,machs(nnodetot),34,MPI_COMM_WORLD,ierr)
                 end if
                 
              
                 ! When this node is finished, send the blocking MPI_Send to the next machine

                 call MPI_Send(ping,1,MPI_INTEGER,sendnum,212,MPI_COMM_WORLD,ierr)
                 
                 ! If this is root, then just copy the array to the global
              else
                            
                 if (ptr_send) then

                    pvar_g(node_idmin:node_idmax) = pvar_l(node_idmin:node_idmax)
           
                 end if
                 
              endif


              
           else
              
              pvar_g = pvar_l
              
           endif
           
           
           ! The data over the desired range of indices have been collected
           ! if this is the only machine or the root machine, then print it 
           ! to standard output
           
           if (mynum .eq. nnodetot .or. nnodetot .eq. 1) then

              ! Print out a maximum of 10 variables per row...

              maxrows = ceiling(real(npolys)/real(maxcols))

              do row = 1,maxrows
                 
                 ncols = min( maxcols,npolys-((row-1)*maxcols)   )
                 col   = ( (row-1)*maxcols)+1
                 
                 write(fmtstr,'(i3)')ncols
                 fmtstr = '('// trim(fmtstr)  // '(2x,' // trim(pfmtstr(ip)) // '))'
                 
                 print(trim(fmtstr)),(pvar_g(i),i=col,col+ncols-1)
              enddo
              print*,""
              print*,""
              
           endif
           
        endif
        
     enddo

     ! Check to see if we matched the variable
     if(.not.pvartrue) then
        print*,"The diagnostic variable named:",trim(pvar_name)
        print*,"does not match any of the var_table variables"
        print*,"Check you namelist entries, and the variable "
        print*,"registry and/or remove this"
        print*,"diagnostic variable."
        stop
     endif

  enddo count_pvars

  ! Dont proceed until everything is written out
  ! if this is not done, then there will be other writing
  ! contaminating the output, and thats icky

!  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  
  return
end subroutine print_array

! =======================================================

subroutine fillvar_l(pvar_l,vt_ptr,npts_out,npts_in,out1,in1,in2)
  
  implicit none
  real,dimension(npts_out)  :: pvar_l
  real,dimension(npts_in)   :: vt_ptr
  integer,intent(in)      :: npts_in,npts_out
  integer,intent(in)      :: out1,in1,in2
  integer :: i,j
  
  j = out1
  do i = in1,in2     
     pvar_l(j) = vt_ptr(i)
     j = j + 1
  enddo
  return
end subroutine fillvar_l

