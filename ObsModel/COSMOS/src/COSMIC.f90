!========================================================================================
! COsmic-ray Soil Moisture Interaction Code (COSMIC) - Version 1.5
!
! W. James Shuttleworth and Rafael Rosolem - January/2012
! Additional support: Marek Zreda, Trenton Franz, Xubin Zeng, and Christopher Zweck
! Fortran code developed by Rafael Rosolem
! 
! Shuttleworth, J., R. Rosolem, M. Zreda, and T. Franz (2013), The COsmic-ray Soil
!               Moisture Interaction Code (COSMIC) for use in data assimilation,
!               Hydrol. Earth Syst. Sci., 17(8), 3205–3217, doi:10.5194/hess-17-3205-2013.
! Rosolem, R., et al. (2014). "Translating aboveground cosmic-ray neutron intensity 
!               to high-frequency soil moisture profiles at sub-kilometer scale." 
!               Hydrology and Earth System Sciences 18(11): 4363-4379.
! COSMIC has been developed under the COsmic-ray Soil Moisture Observing System
! (COSMOS) project. The COSMOS project is funded by the Atmospheric Science, Hydrology,
! and Ecology Programs of the US National Science Foundation (grant ATM-0838491).
! 
! Please contact Rafael Rosolem at rafael.rosolem@bristol.ac.uk for any questions and
! support.
! 
!
!========================================================================================
! Updates:
! 2012/01/20 - Version 1.0: * Original version
! 2012/01/28 - Version 1.1: * Some parameters are re-defined for better physical realism
! 2012/02/17 - Version 1.2: * After contribution from all angles are taken, need to
!                             multiply fastflux by 2/pi
!                           * Angle increments can be specified here (ideg)
! 2012/02/29 - Version 1.3: * Reduced number of parameters based on relationship of ns
!                             and nw (now given as alpha = ns/nw)
! 2012/04/03 - Version 1.4: * Soil thickness (i.e., input.dat file) needs to be specified
!                             at finer resolution (i.e., 0.1 cm)
! 2012/04/04 - Version 1.5: * Now the contributions to soil and water densities/mass are
!                             taken to be at the center of a given soil layer
! 2012/10/08 - Version 1.5: * Includes modifications made in DART (obs_def_COSMOS_mod.f90) 
!========================================================================================

subroutine cosmic(nlyr, soil_moisture, layerz, bd_in, lattwat_in, N_in, alpha_in, L1_in, L2_in, L3_in, L4_in, nsteps, nlevels, nthreads, val, soil_moisture_weighted, sensor_depth_weighted, istatus)
!----------------------------------------------------------------------
! Uses a weighting function calculated by COSMIC (COsmic-ray Soil Moisture Interaction Code)
use omp_lib

real(8), parameter :: PI = 3.14159265358979323846

integer,             intent(in)  :: nlyr ! Total number of soil layers - ( default 3000 each 1mm for 3m)
real(8),             intent(in)  :: soil_moisture(nsteps,nlevels) ! soil_moisture vector
real(8),             intent(in)  :: layerz(nsteps,nlevels) ! original soil layer depths, COSMIC needs them to be positive and in centimeters
real(8), 			 intent(in)  :: bd_in(nsteps), lattwat_in(nsteps), N_in(nsteps), alpha_in(nsteps), L1_in(nsteps), L2_in(nsteps), L3_in(nsteps), L4_in(nsteps)
integer,             intent(in)  :: nsteps  ! num of time steps
integer,             intent(in)  :: nlevels  ! num of soil layers
integer,             intent(in)  :: nthreads  ! num of nthreads
real(8),             intent(out) :: val(nsteps),soil_moisture_weighted(nsteps),sensor_depth_weighted(nsteps)      ! value of obs
integer,             intent(out) :: istatus  ! status of the calculation
character(len=256) :: string1, string2

!=================================================================================
! COSMIC: Variables list
!=================================================================================

real(8) :: bd     = 0.0 ! Dry soil bulk density (g/m3)
real(8) :: vwclat = 0.0 ! Volumetric "lattice" water content (m3/m3)
real(8) :: N      = 0.0 ! High energy neutron flux (-)
real(8) :: alpha  = 0.0 ! Ratio of Fast Neutron Creation Factor (Soil to Water), alpha (-)
real(8) :: L1     = 0.0 ! High Energy Soil   Attenuation Length (g/cm2)
real(8) :: L2     = 0.0 ! High Energy Water  Attenuation Length (g/cm2)
real(8) :: L3     = 0.0 ! Fast Neutron Soil  Attenuation Length (g/cm2)
real(8) :: L4     = 0.0 ! Fast Neutron Water Attenuation Length (g/cm2)
real(8) :: zdeg
real(8) :: zrad
real(8) :: ideg
real(8) :: costheta
real(8) :: dtheta
real(8) :: totflux     ! Total flux of above-ground fast neutrons

real(8), dimension(:), allocatable :: dz          ! Soil layers (cm)
real(8), dimension(:), allocatable :: zthick      ! Soil layer thickness (cm)
real(8), dimension(:), allocatable :: vwc         ! Volumetric Water Content (m3/m3)
real(8), dimension(:), allocatable :: isoimass    ! Integrated dry soil mass above layer (g)
real(8), dimension(:), allocatable :: iwatmass    ! Integrated water mass above layer (g)
real(8), dimension(:), allocatable :: hiflux      ! High energy neutron flux
real(8), dimension(:), allocatable :: fastpot     ! Fast neutron source strength of layer
real(8), dimension(:), allocatable :: h2oeffdens  ! "Effective" density of water in layer (g/cm3)
real(8), dimension(:), allocatable :: idegrad     ! Integrated neutron degradation factor (-)
real(8), dimension(:), allocatable :: fastflux    ! Contribution to above-ground neutron flux
real(8), dimension(:), allocatable :: vwc_w    	  ! soil moisture weighted

!rr: Not needed for DART
!rr: real(8), dimension(:), allocatable :: normfast ! Normalized contribution to neutron flux (-) [weighting factors]

real(8), parameter   :: h2odens = 1000.0 ! Density of water (g/cm3)

integer,  parameter   :: maxlayers = 1000 ! more than the maximum # of model soil layers
!real(8), allocatable :: layerz(:)        ! original soil layer depths
!real(8), allocatable :: soil_moisture(:) ! original soil layer moistures

integer  :: angle, angledz, maxangle  ! loop indices for an integration interval
integer  :: i, zi, step

CALL OMP_SET_NUM_THREADS(nthreads)

!$OMP PARALLEL SHARED(val) PRIVATE(step,bd,vwclat,N,alpha,L1,L2,L3,L4,i,dz,zthick,vwc,hiflux,fastpot,h2oeffdens,idegrad,fastflux,isoimass,iwatmass,vwc_w,totflux,zi,ideg,angledz,maxangle,dtheta,string1)
!$OMP DO SCHEDULE(DYNAMIC,1)	
do step = 1,nsteps
	!=================================================================================
	
	val(step) = 0.0 ! set return value early
	soil_moisture_weighted(step) = 0.0
	
	!=================================================================================
	! COSMIC: Site specific-parameters come from the observation metadata
	!=================================================================================

	bd     = bd_in(step)
	vwclat = lattwat_in(step)
	N      = N_in(step)
	alpha  = alpha_in(step)
	L1     = L1_in(step)
	L2     = L2_in(step)
	L3     = L3_in(step)
	L4     = L4_in(step)
	!print*,bd,vwclat,N,alpha,L1,L2,L3,L4
	!=================================================================================
	! COSMIC: Allocate arrays and initialize variables
	!=================================================================================

	allocate(dz(nlyr), zthick(nlyr), vwc(nlyr), &
			 hiflux(nlyr), fastpot(nlyr), h2oeffdens(nlyr), &
			idegrad(nlyr), fastflux(nlyr), &
			isoimass(nlyr), iwatmass(nlyr), vwc_w(nlyr))

	totflux = 0.0
	
	do i = 1,nlyr

	   dz(i)         = (real(i,8)/10.0) * (3000.0/nlyr)	! default 0.1 cm intervals
	   zthick(i)     = 0.0
	   h2oeffdens(i) = 0.0
	   vwc(i)        = 0.0
	   isoimass(i)   = 0.0
	   iwatmass(i)   = 0.0
	   hiflux(i)     = 0.0
	   fastpot(i)    = 0.0
	   idegrad(i)    = 0.0
	   fastflux(i)   = 0.0
	   vwc_w(i)		 = 0.0
		!rr: Not needed for DART
		!rr: normfast(i)    = 0.0

	enddo
	
	!=================================================================================
	! Get soil moisture from individual model layers and assign them to
	! 1 mm intervals (down to 3 meters)
	
	do i = 1,nlyr

	   SOIL : do zi = 1,nlevels
		  if (dz(i) >= layerz(step,nlevels)) then
			 vwc(i) = soil_moisture(step,nlevels)
			 exit SOIL
		  elseif (dz(i) <= layerz(step,zi)) then
			 vwc(i) = soil_moisture(step,zi) ! soil moisture (m3/m3)
			 exit SOIL
		  endif
	   enddo SOIL

	enddo
	
	!=================================================================================
	! COSMIC: Neutron flux calculation
	!=================================================================================

	! At some point, you might want to tinker around with non-uniform
	! soil layer thicknesses.
	zthick(1) = dz(1) - 0.0 ! Surface layer
	do i = 2,nlyr
	   zthick(i) = dz(i) - dz(i-1) ! Remaining layers
	enddo
	
	! Angle distribution parameters (HARDWIRED)
	!rr: Using 0.5 deg angle intervals appears to be sufficient
	!rr: (smaller angles increase the computing time for COSMIC)
	ideg     = 0.5                   ! ideg ultimately controls the number of trips through
	angledz  = nint(ideg*10.0)       ! the ANGLE loop. Make sure the 10.0 is enough
	maxangle = 900 - angledz            ! to create integers with no remainder
	dtheta   = ideg*(PI/180.0)
	
	if ( real(angledz,8) /= ideg*10.0 ) then
	   write(string1,*) 'ideg*10.0 must result in an integer - it results in ',ideg*10.0
	   !call error_handler(E_ERR,'get_expected_neutron_intensity',string1,source,revision,revdate)
	endif
	
	do i = 1,nlyr

	   ! High energy neutron downward flux
	   ! The integration is now performed at the node of each layer (i.e., center of the layer)

	   h2oeffdens(i) = ((vwc(i)+vwclat)*h2odens)/1000.0
	   !print*,h2oeffdens(i),vwc(i),vwclat
	   if(i > 1) then
		  ! Assuming an area of 1 cm2
		  isoimass(i) = isoimass(i-1) + bd*(0.5*zthick(i-1))*1.0 + &
										bd*(0.5*zthick(i  ))*1.0
		  ! Assuming an area of 1 cm2
		  iwatmass(i) = iwatmass(i-1) + h2oeffdens(i-1)*(0.5*zthick(i-1))*1.0 + &
										h2oeffdens(i  )*(0.5*zthick(i  ))*1.0
	   else
		  isoimass(i) =            bd*(0.5*zthick(i))*1.0 ! Assuming an area of 1 cm2
		  iwatmass(i) = h2oeffdens(i)*(0.5*zthick(i))*1.0 ! Assuming an area of 1 cm2
	   endif

	   hiflux( i) = N*exp(-(isoimass(i)/L1 + iwatmass(i)/L2) )
	   fastpot(i) = zthick(i)*hiflux(i)*(alpha*bd + h2oeffdens(i))
	   
	   ! This second loop needs to be done for the distribution of angles for fast neutron release
	   ! the intent is to loop from 0 to 89.5 by 0.5 degrees - or similar.
	   ! Because Fortran loop indices are integers, we have to divide the indices by 10 - you get the idea.

	   do angle=0,maxangle,angledz
		  zdeg     = real(angle,8)/10.0   ! 0.0  0.5  1.0  1.5 ...
		  zrad     = (zdeg*PI)/180.0
		  costheta = cos(zrad)

		  ! Angle-dependent low energy (fast) neutron upward flux
		  fastflux(i) = fastflux(i) + fastpot(i)*exp(-(isoimass(i)/L3 + iwatmass(i)/L4)/costheta)*dtheta
	   enddo
	   
	   ! After contribution from all directions are taken into account,
	   ! need to multiply fastflux by 2/PI

	   fastflux(i) = (2.0/PI)*fastflux(i)
	   
	   ! Low energy (fast) neutron upward flux
	   totflux = totflux + fastflux(i)
	   
	enddo
	
	! Calculate the weighted soil moisture
	vwc_w = fastflux / totflux
	soil_moisture_weighted(step) = sum(vwc * vwc_w)
	!print*,soil_moisture_weighted(step),sum(vwc_w)

	do i = 2,nlyr
	    vwc_w(i) = vwc_w(i) + vwc_w(i-1)
	    if (vwc_w(i) >= 0.86) then
	        sensor_depth_weighted(step) = i
	        exit
	    endif
	enddo

	!print*,sensor_depth_weighted(step)
	
	deallocate(dz, zthick, vwc, hiflux, fastpot, &
			   h2oeffdens, idegrad, fastflux, isoimass, iwatmass, vwc_w)

	!=================================================================================
	! ... and finally set the return the neutron intensity (i.e. totflux)
	
	val(step)     = totflux
	istatus = 0        ! assume all is well if we get this far

enddo
!$OMP END PARALLEL

return

end subroutine cosmic
