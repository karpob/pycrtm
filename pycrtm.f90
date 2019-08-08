module pycrtm 
contains
subroutine wrap_forward( coefficientPath, sensor_id_in, & 
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, &
                        nChan, N_Profiles, N_LAYERS, &
                        pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers, & 
                        co2ConcLayers, & 
                        aerosolEffectiveRadius, aerosolConcentration, aerosolType, & 
                        cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m, n_absorbers, & 
                        landType, soilType, vegType, waterType, snowType, iceType, nthreads, &  
                        outTb, outTransmission, & 
                        emissivity )      

  ! ============================================================================
  ! STEP 1. **** ENVIRONMENT SETUP FOR CRTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================
  ! variables for interface
  character(len=*), intent(in) :: coefficientPath
  character(len=*), intent(in) :: sensor_id_in
  ! The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  integer, intent(in) :: nChan, N_Profiles, N_Layers
  real(kind=8), intent(in) :: zenithAngle(n_profiles), scanAngle(n_profiles) 
  real(kind=8), intent(in) :: azimuthAngle(n_profiles), solarAngle(2,n_profiles)
  real(kind=8), intent(in) :: pressureLevels(N_LAYERS+1, N_Profiles)
  real(kind=8), intent(in) :: pressureLayers(N_LAYERS, N_Profiles), temperatureLayers(N_LAYERS, N_Profiles)
  real(kind=8), intent(in) ::humidityLayers(N_LAYERS,N_Profiles), ozoneConcLayers(N_LAYERS, N_Profiles)
  real(kind=8), intent(in) :: co2ConcLayers(N_LAYERS, N_Profiles)
  real(kind=8), intent(in) :: aerosolEffectiveRadius(N_LAYERS, N_Profiles), aerosolConcentration(N_LAYERS, N_Profiles)
  real(kind=8), intent(in) :: cloudEffectiveRadius(N_LAYERS, N_Profiles), cloudConcentration(N_LAYERS, N_Profiles) 
  real(kind=8), intent(in) :: cloudFraction(N_LAYERS, N_Profiles)
  integer, intent(in) :: aerosolType(N_Profiles), cloudType(N_Profiles)
  integer, intent(in) :: n_absorbers(N_Profiles), climatology(N_profiles)
  real(kind=8), intent(in) :: surfaceTemperatures(4, N_Profiles), surfaceFractions(4, N_Profiles)
  real(kind=8), intent(in) :: LAI(N_Profiles), salinity(N_Profiles),  windSpeed10m(N_Profiles), windDirection10m(N_Profiles)
  integer, intent(in) ::  landType(N_Profiles), soilType(N_Profiles), vegType(N_Profiles), waterType(N_Profiles)
  integer, intent(in) ::  snowType(N_Profiles), iceType(N_Profiles) 
  integer, intent(in) :: nthreads
  real(kind=8), intent(out) :: outTb(nChan, N_Profiles), emissivity(nChan, N_Profiles)
  real(kind=8), intent(out) :: outTransmission(nChan, N_LAYERS, N_Profiles)
  character(len=256), dimension(1) :: sensor_id

  ! --------------------------
  ! Some non-CRTM-y Parameters
  ! --------------------------
  CHARACTER(*), PARAMETER :: SUBROUTINE_NAME   = 'wrap_forward'
  CHARACTER(*), PARAMETER :: PROGRAM_VERSION_ID = '0.01'



  ! ============================================================================
  ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
  !
  ! Directory location of coefficients

  ! Profile dimensions
  INTEGER :: N_CLOUDS
  INTEGER :: N_AEROSOLS
  
  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 1
  ! ============================================================================
  


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: message, version
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels
  INTEGER :: i, l, m, n, nc, ll,mm, nn
  real, dimension(n_layers,nchan) :: outTau
  logical :: cloudsOn, aerosolsOn

  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(1)
  TYPE(CRTM_Geometry_type)                :: geo(1)


  ! 3b. Define the FORWARD variables
  ! --------------------------------
  TYPE(CRTM_Atmosphere_type)              :: atm(1)
  TYPE(CRTM_Surface_type)                 :: sfc(1)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
  type(crtm_options_type)                :: options
  sensor_id(1) = sensor_id_in
  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( SUBROUTINE_NAME, &
    'Running simulation.', &
    'CRTM Version: '//TRIM(Version) )

  ! ============================================================================
  ! STEP 4. **** INITIALIZE THE CRTM ****
  !
  ! 4a. Initialise all the sensors at once
  ! --------------------------------------
  ! Karpowicz addition... if we have less than -9999 for all concentrations, don't turn on aerosols/clouds.
  if( all(aerosolConcentration < -9999.0) ) then
    N_AEROSOLS = 0
    aerosolsOn = .False.
  else
    N_AEROSOLS = 1
    aerosolsOn = .True. 
  endif

  if( all(cloudConcentration < -9999.0) ) then
    N_CLOUDS = 0
    cloudsOn = .False.
  else
    N_CLOUDS = 1
    cloudsOn = .True. 
  endif
  ! End Karpowicz change to CRTM-style interface.

   err_stat = CRTM_Init( sensor_id,  chinfo, &
                        File_Path=coefficientPath, &
                        Load_CloudCoeff = cloudsOn, &  
                        Load_AerosolCoeff = aerosolsOn, &  
                        Quiet=.True. )

  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error initializing CRTM'
    CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
    STOP
  END IF


  WRITE( *,'(/5x,"Initializing the CRTM...")' )

  ! if aerosols or cloud concentrations specified as < -9999, don't load cloud/aerosol coefficients.


  ! 4b. Output some channel information
  ! -----------------------------------
  n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels




  ! ============================================================================

    ! ==========================================================================
    ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
    !
    ! 5a. Determine the number of channels
    !     for the current sensor
    ! ------------------------------------
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))
    
    ! 5b. Allocate the ARRAYS
    ! -----------------------
    ALLOCATE( rts( n_channels, 1), STAT = alloc_stat )

    IF ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF
 

  ! Begin loop over profile
  ! ----------------------
  !$ call omp_set_num_threads(nthreads)
  !$omp parallel do default(firstprivate) shared(chinfo,emissivity,outTb)& 
  !$omp& shared(zenithAngle, scanAngle, azimuthAngle, solarAngle)&
  !$omp& shared(nChan, N_Profiles, N_LAYERS)&
  !$omp& shared(pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers)& 
  !$omp& shared(co2ConcLayers)& 
  !$omp& shared(aerosolEffectiveRadius, aerosolConcentration, aerosolType)& 
  !$omp& shared(cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology)& 
  !$omp& shared(surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m, n_absorbers)& 
  !$omp& shared(landType, soilType, vegType, waterType, snowType, iceType)&
  !$omp& num_threads(nthreads) 
  Profile_Loop: DO n = 1, N_Profiles
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))

    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS(n), N_CLOUDS, N_AEROSOLS )
    IF ( ANY (.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      !STOP
    END IF
  

    ! ==========================================================================
    ! STEP 6. **** ASSIGN INPUT DATA ****
    !
    ! 6a. Atmosphere and Surface input
    !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
    !     by which the atmosphere and surface data are loaded in to their
    !     respective structures below was done purely to keep the step-by-step
    !     instructions in this program relatively "clean".
    ! ------------------------------------------------------------------------
    atm(1)%Absorber_Id(1:2)    = (/ H2O_ID                 , O3_ID /)
    atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    ! ...Profile data
    atm(1)%Climatology = climatology(n)
    atm(1)%Level_Pressure = pressureLevels(:,n)
    atm(1)%Pressure = pressureLayers(:,n)
    atm(1)%Temperature = temperatureLayers(:,n)
    atm(1)%Absorber(:,1) = humidityLayers(:,n)
    atm(1)%Absorber(:,2) = ozoneConcLayers(:,n)

    if( aerosolsOn )  then
      atm(1)%Aerosol(1)%Type                = aerosolType(n)
      atm(1)%Aerosol(1)%Effective_Radius(:) = aerosolEffectiveRadius(:,n)
      atm(1)%Aerosol(1)%Concentration(:)    = aerosolConcentration(:,n)
    endif

    if( cloudsOn ) then
      atm(1)%Cloud(1)%Type                = cloudType(n)
      atm(1)%Cloud(1)%Effective_Radius(:) = cloudEffectiveRadius(:,n)
      atm(1)%Cloud(1)%Water_Content(:)    = cloudConcentration(:,n)
      atm(1)%Cloud_Fraction(:)            = cloudFraction(:,n)
    endif

    if(n_absorbers(n) >2) then 
      atm(1)%Absorber(:,3)     = co2ConcLayers(:,n)
      atm(1)%Absorber_Id(3)    = CO2_ID
      atm(1)%absorber_units(3) = VOLUME_MIXING_RATIO_UNITS
    endif


    ! 6b. Geometry input
    ! ------------------
    CALL CRTM_Geometry_SetValue( geo, &
                                 Sensor_Zenith_Angle  = zenithAngle(n),   &
                                 Sensor_Scan_Angle    = scanAngle(n),     & 
                                 Sensor_Azimuth_Angle = azimuthAngle(n),  &  
                                 Source_Zenith_Angle  = solarAngle(1,n),  & 
                                 Source_Azimuth_Angle = solarAngle(2,n) )
    !print *, atm(1)%Temperature
    !geo%Sensor_Zenith_Angle = zenithAngle(n)
    !geo%Sensor_Scan_Angle = scanAngle(n)
    !geo%Sensor_Azimuth_Angle = azimuthAngle(n) 
    ! ==========================================================================
    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    sfc%Land_Coverage     = surfaceFractions(1,n)
    sfc%Land_Type         = landType(n) !TUNDRA_SURFACE_TYPE
    sfc%Land_Temperature  = surfaceTemperatures(1,n)
    sfc%Lai               = LAI(n)
    sfc%Soil_Type         = soilType(n) !COARSE_SOIL_TYPE
    sfc%Vegetation_Type   = vegType(n) !GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc%Water_Coverage    = surfaceFractions(2,n)
    sfc%Water_Type        = waterType(n) !SEA_WATER_TYPE
    sfc%Water_Temperature = surfaceTemperatures(2,n)

    ! ...Snow coverage characteristics
    sfc%Snow_Coverage    = surfaceFractions(3,n)
    sfc%Snow_Type        = snowType(n) !FRESH_SNOW_TYPE
    sfc%Snow_Temperature = surfaceTemperatures(3,n)
    ! ...Ice surface characteristics
    sfc%Ice_Coverage    = surfaceFractions(4,n)
    sfc%Ice_Type        = iceType(n) !FRESH_ICE_TYPE
    sfc%Ice_Temperature = surfaceTemperatures(4,n)

   
    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    !WRITE( *, '( /5x, "Calling the CRTM functions for ",a,"..." )' ) TRIM(sensor_id(1))
    
    ! 8a. The forward model
    ! ---------------------

    ! Need this to get transmission out of solution, otherwise won't be allocated !!!
    !call crtm_rtsolution_create( rts, n_layers )
    !if ( any(.not. crtm_rtsolution_associated( rts )) ) then
    !    call display_message( subroutine_name, 'error allocating rts', err_stat)
    !    return
    !end if

    ! why did I put this here????
    !call crtm_options_create( options, nChan )
    !if ( any(.not. crtm_options_associated( options )) ) then 
    !    call display_message( subroutine_name, 'error allocating options', FAILURE)  
    !    return
    !endif

    !options%Use_Emissivity = .false.
    !options%Use_Direct_Reflectivity = .false.

    !! Karpowicz addition to CRTM-style interface, switch to force scattering off when aerosols and clouds turned off.
    !if (.not. cloudsOn .and. .not. aerosolsOn) then 
    !    options%Include_Scattering = .false.
    !end if 

    err_stat = CRTM_Forward( atm        , &  ! Input
                             sfc        , &  ! Input
                             geo        , &  ! Input
                             chinfo, &  ! Input
                             rts ) !,    & ! Output
    !                        options = options ) 
    IF ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM Forward Model for '//TRIM(sensor_id(1))
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF
    ! ============================================================================
    ! 8c. **** OUTPUT THE RESULTS TO SCREEN **** (Or transfer it into a series of arrays out of this thing!)
    !
    ! User should read the user guide or the source code of the routine
    ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
    ! select the needed variables for outputs.  These variables are contained
    ! in the structure RTSolution.
    ! do l=1,nChan
    !    outTransmission(l,1:n_layers, n) = rts(l,1)%Layer_Optical_Depth
    !enddo
    emissivity(:,n) = rts(:,1)%Surface_Emissivity 
    outTb(:,n) = rts(:,1)%Brightness_Temperature 
    
    ! ==========================================================================
    ! STEP 9. **** CLEAN UP FOR NEXT PROFILE ****
    !
    ! 9a. Deallocate the structures
    ! -----------------------------


    ! 9b. Deallocate the arrays
    ! -------------------------
    ! ==========================================================================
    CALL CRTM_Atmosphere_Destroy(atm)
  END DO Profile_Loop
  !$omp end parallel do

  
    DEALLOCATE(rts, STAT = alloc_stat)
  
  ! ==========================================================================
  ! 10. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
    STOP
  END IF
  ! ==========================================================================
end subroutine wrap_forward

subroutine wrap_k_matrix( coefficientPath, sensor_id_in, & 
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, &  
                        nChan, N_profiles, N_LAYERS, & 
                        pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers, & 
                        co2ConcLayers, & 
                        aerosolEffectiveRadius, aerosolConcentration, aerosolType, & 
                        cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, salinity, windSpeed10m, windDirection10m, n_absorbers, & 
                        landType, soilType, vegType, waterType, snowType, iceType, &  
                        outTb, outTransmission, & 
                        temperatureJacobian, humidityJacobian, ozoneJacobian, emissivity, nthreads )      

  ! ============================================================================
  ! STEP 1. **** ENVIRONMENT SETUP FOR CRTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================
  


  ! --------------------------
  ! Some non-CRTM-y Parameters
  ! --------------------------
  CHARACTER(*), PARAMETER :: SUBROUTINE_NAME   = 'wrap_k_matrix'
  CHARACTER(*), PARAMETER :: PROGRAM_VERSION_ID = '0.01'

  ! variables for interface
  character(len=1024), intent(in) :: coefficientPath
  character(len=*), intent(in) :: sensor_id_in
  integer, intent(in) :: nChan, N_profiles, N_Layers 
  ! The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  real(kind=8), intent(in) :: zenithAngle(N_profiles), scanAngle(N_profiles)
  real(kind=8), intent(in) :: azimuthAngle(N_profiles), solarAngle(2,N_profiles)
  real(kind=8), intent(in) :: pressureLevels(N_LAYERS+1, N_profiles)
  real(kind=8), intent(in) :: pressureLayers(N_LAYERS, N_profiles), temperatureLayers(N_LAYERS, N_profiles)
  real(kind=8), intent(in) :: humidityLayers(N_LAYERS, N_profiles)
  real(kind=8), intent(in) :: ozoneConcLayers(N_LAYERS, N_profiles)
  real(kind=8), intent(in) :: co2ConcLayers(N_LAYERS, N_profiles)
  real(kind=8), intent(in) :: aerosolEffectiveRadius(N_LAYERS, N_profiles), aerosolConcentration(N_LAYERS, N_profiles)
  real(kind=8), intent(in) :: cloudEffectiveRadius(N_LAYERS, N_profiles) 
  real(kind=8), intent(in) :: cloudConcentration(N_LAYERS, N_profiles), cloudFraction(N_LAYERS,N_profiles)
  integer, intent(in) :: aerosolType(N_profiles), cloudType(N_profiles), n_absorbers(N_profiles), climatology(N_profiles)
  real(kind=8), intent(in) :: surfaceTemperatures(4, N_profiles), surfaceFractions(4, N_profiles), LAI(N_profiles) 
  real(kind=8), intent(in) :: salinity(N_profiles), windSpeed10m(N_profiles), windDirection10m(N_profiles)
  integer, intent(in) :: landType(N_profiles), soilType(N_profiles), vegType(N_profiles), waterType(N_profiles) 
  integer, intent(in) :: snowType(N_profiles), iceType(N_profiles) 
  real(kind=8), intent(out) :: outTb(nChan, N_profiles), emissivity(nChan, N_profiles)
  real(kind=8), intent(out) :: outTransmission(nChan, N_LAYERS, N_profiles) 
  real(kind=8), intent(out) :: temperatureJacobian(nChan,N_LAYERS, N_profiles)
  real(kind=8), intent(out) ::  humidityJacobian(nChan, N_LAYERS, N_profiles)
  real(kind=8), intent(out) :: ozoneJacobian(nChan, N_LAYERS, N_profiles)
  integer, intent(in) :: nthreads
  

  character(len=256) :: sensor_id(1)
  ! ============================================================================
  ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
  !
  ! Directory location of coefficients

  ! Profile dimensions
  INTEGER :: N_CLOUDS 
  INTEGER :: N_AEROSOLS
  
  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 1
  ! ============================================================================
  


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: message, version
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels
  INTEGER :: l, m, n, nc
  Logical :: cloudsOn, aerosolsOn


  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(1)
  TYPE(CRTM_Geometry_type)                :: geo(1)

  ! 3b. Define the FORWARD variables
  ! --------------------------------
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm(:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfc(:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
 
  ! 3c. Define the K-MATRIX variables
  ! ---------------------------------
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_K(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: sfc_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_K(:,:)
  ! ============================================================================


  sensor_id(1) = sensor_id_in
  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( SUBROUTINE_NAME, &
    'Running simulation.', &
    'CRTM Version: '//TRIM(Version) )



  ! ============================================================================
  ! STEP 4. **** INITIALIZE THE CRTM ****
  !
  ! 4a. Initialise all the sensors at once
  ! --------------------------------------
  if( all(aerosolConcentration < -9999.0) ) then
    N_AEROSOLS = 0
    aerosolsOn = .False.
  else
    N_AEROSOLS = 1
    aerosolsOn = .True. 
  endif

  if( all(cloudConcentration < -9999.0) ) then
    N_CLOUDS = 0
    cloudsOn = .False.
  else
    N_CLOUDS = 1
    cloudsOn = .True. 
  endif

  WRITE( *,'(/5x,"Initializing the CRTM...")' )
 
 ! if aerosols or cloud concentrations specified as < -9999, don't load cloud/aerosol coefficients.

  err_stat = CRTM_Init( sensor_id,  chinfo, &
                        File_Path=coefficientPath, &
                        Load_CloudCoeff = cloudsOn, &  
                        Load_AerosolCoeff = aerosolsOn, &  
                        Quiet=.True. )


  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error initializing CRTM'
    CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
    STOP
  END IF

  ! 4b. Output some channel information
  ! -----------------------------------
  n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels
  WRITE( *,'(7x,i0," from ",a)' )  CRTM_ChannelInfo_n_Channels(chinfo(1)), TRIM(sensor_id(1))
  ! ============================================================================

    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure

!    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS(1), N_CLOUDS, N_AEROSOLS )
!    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
!      message = 'Error allocating CRTM Forward Atmosphere structure'
!      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
!      STOP
!    END IF
 
  ! Begin loop over sensors
  ! ----------------------
  ! openmp won't work for K-matrix. Should be a fix in future release of CRTM. 
  !!$ call omp_set_num_threads(nthreads)
  !!$omp parallel do default(firstprivate) shared(chinfo,emissivity,outTb,atm_k,rts_k)& 
  !!$omp& shared(temperatureJacobian,humidityJacobian,ozoneJacobian,outTransmission)&
  !!$omp& shared(zenithAngle, scanAngle, azimuthAngle, solarAngle)&
  !!$omp& shared(nChan, N_Profiles, N_LAYERS)&
  !!$omp& shared(pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers)& 
  !!$omp& shared(co2ConcLayers)& 
  !!$omp& shared(aerosolEffectiveRadius, aerosolConcentration, aerosolType)& 
  !!$omp& shared(cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology)& 
  !!$omp& shared(surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m, n_absorbers)& 
  !!$omp& shared(landType, soilType, vegType, waterType, snowType, iceType)&
  !!$omp& num_threads(nthreads) 
 
  Profile_Loop: DO n = 1, N_profiles

  
    ! ==========================================================================
    ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
    !
    ! 5a. Determine the number of channels
    !     for the current sensor
    ! ------------------------------------
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))

    
    ! 5b. Allocate the ARRAYS
    ! -----------------------
    allocate(atm(1), sfc(1), STAT = alloc_stat)
    IF (alloc_stat /= 0 ) THEN
      message = 'Error allocating atm, sfc'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    else 
    END IF


    ALLOCATE( rts( n_channels, 1 ), STAT = alloc_stat )
    IF (alloc_stat /= 0 ) THEN
      message = 'Error allocating rts'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    else 
    END IF


    ALLOCATE( atm_K( n_channels, 1 ), STAT = alloc_stat )
    IF (alloc_stat /= 0 ) THEN
      message = 'Error allocating atm_K'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF



    ALLOCATE( sfc_K( n_channels, 1 ), STAT = alloc_stat )
    IF (alloc_stat /= 0 ) THEN
      message = 'Error allocating sfc_K'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF

    ALLOCATE( rts_K( n_channels, 1 ), STAT = alloc_stat )
    IF (alloc_stat /= 0 ) THEN
      message = 'Error allocating rts_K'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF

    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS(n), N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF

    ! The output K-MATRIX structure
    CALL CRTM_Atmosphere_Create( atm_K, N_LAYERS, N_ABSORBERS(n), N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF

    call crtm_rtsolution_create( rts, n_layers )
    if ( any(.not. crtm_rtsolution_associated( rts )) ) then
        call display_message( subroutine_name, 'error allocating rts', err_stat)
        !return
    end if
    
    call crtm_rtsolution_create( rts_k, n_layers )
    if ( any(.not. crtm_rtsolution_associated( rts_k )) ) then
        call display_message( subroutine_name, 'error allocating rts_k', err_stat)
        !return
    end if

    ! ==========================================================================
    ! STEP 6. **** ASSIGN INPUT DATA ****
    !
    ! 6a. Atmosphere and Surface input
    !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
    !     by which the atmosphere and surface data are loaded in to their
    !     respective structures below was done purely to keep the step-by-step
    !     instructions in this program relatively "clean".
    ! ------------------------------------------------------------------------
    atm(1)%Climatology         = climatology(n)  
    atm(1)%Absorber_Id(1:2)    = (/ H2O_ID                 , O3_ID /)
    atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    ! ...Profile data
    atm(1)%Level_Pressure = pressureLevels(:,n)
    atm(1)%Pressure = pressureLayers(:,n)
    atm(1)%Temperature = temperatureLayers(:,n)
    atm(1)%Absorber(:,1) = humidityLayers(:,n)
    atm(1)%Absorber(:,2) = ozoneConcLayers(:,n)
    if( aerosolsOn )  then
        atm(1)%Aerosol(1)%Type                = aerosolType(n)
        atm(1)%Aerosol(1)%Effective_Radius(:) = aerosolEffectiveRadius(:, n)
        atm(1)%Aerosol(1)%Concentration(:)    = aerosolConcentration(:, n)
    endif

    if( cloudsOn ) then
        atm(1)%Cloud(1)%Type                = cloudType(n)
        atm(1)%Cloud(1)%Effective_Radius(:) = cloudEffectiveRadius(:, n)
        atm(1)%Cloud(1)%Water_Content(:)    = cloudConcentration(:, n)
        atm(1)%Cloud_Fraction(:)            = cloudFraction(:, n)
    endif

    if( n_absorbers(n) > 2 ) then 
        atm(1)%Absorber(:,3)      = co2ConcLayers(:,n)
        atm(1)%Absorber_Id(3)     = CO2_ID
        atm(1)%absorber_units(3)  = VOLUME_MIXING_RATIO_UNITS
    endif

    ! 6b. Geometry input
    ! ------------------
    ! All profiles are given the same value
    CALL CRTM_Geometry_SetValue( geo, &
                                 Sensor_Zenith_Angle  = zenithAngle(n),   &
                                 Sensor_Scan_Angle    = scanAngle(n),     & 
                                 Sensor_Azimuth_Angle = azimuthAngle(n),  &  
                                 Source_Zenith_Angle  = solarAngle(1,n),  & 
                                 Source_Azimuth_Angle = solarAngle(2,n) )
 
    ! ==========================================================================




    ! ==========================================================================
    ! STEP 7. **** INITIALIZE THE K-MATRIX ARGUMENTS ****
    !
    ! 7a. Zero the K-matrix OUTPUT structures
    ! ---------------------------------------
    CALL CRTM_Atmosphere_Zero( atm_K )
    CALL CRTM_Surface_Zero( sfc_K )


    ! 7b. Inintialize the K-matrix INPUT so
    !     that the results are dTb/dx
    ! -------------------------------------
    rts_K%Radiance               = ZERO
    rts_K%Brightness_Temperature = ONE
    ! ==========================================================================
    
    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    sfc%Land_Coverage     = surfaceFractions(1,n)
    sfc%Land_Type         = landType(n) !TUNDRA_SURFACE_TYPE
    sfc%Land_Temperature  = surfaceTemperatures(1,n)
    sfc%Lai               = LAI(n)
    sfc%Soil_Type         = soilType(n) !COARSE_SOIL_TYPE
    sfc%Vegetation_Type   = vegType(n) !GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc%Water_Coverage    = surfaceFractions(2,n)
    sfc%Water_Type        = waterType(n) !SEA_WATER_TYPE
    sfc%Water_Temperature = surfaceTemperatures(2,n)

    ! ...Snow coverage characteristics
    sfc%Snow_Coverage    = surfaceFractions(3,n)
    sfc%Snow_Type        = snowType(n) !FRESH_SNOW_TYPE
    sfc%Snow_Temperature = surfaceTemperatures(3,n)
    ! ...Ice surface characteristics
    sfc%Ice_Coverage    = surfaceFractions(4,n)
    sfc%Ice_Type        = iceType(n) !FRESH_ICE_TYPE
    sfc%Ice_Temperature = surfaceTemperatures(4,n)
    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    ! 8b. The K-matrix model
    ! ----------------------
!    err_stat = CRTM_Forward( atm        , &  ! Input
!                             sfc        , &  ! Input
!                             geo        , &  ! Input
!                             chinfo, &  ! Input
!                             rts ) !,    & ! Output
!    !                        options = options ) 
! 
    err_stat = CRTM_K_Matrix( atm        , &  ! FORWARD  Input
                              sfc        , &  ! FORWARD  Input
                              rts_K      , &  ! K-MATRIX Input
                              geo        , &  ! Input
                              chinfo     , &  ! Input
                              atm_K      , &  ! K-MATRIX Output
                              sfc_K      , &  ! K-MATRIX Output
                              rts          )  ! FORWARD  Output
    
    IF ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM K-Matrix Model'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF

    ! ==========================================================================
    ! STEP 9. **** CLEAN UP FOR NEXT SENSOR ****
    !
    ! 9a. Deallocate the structures
    ! -----------------------------
    

    ! 9b. Deallocate the arrays
    ! -------------------------
    ! transfer jacobians out
    do l=1,nChan
        temperatureJacobian(l, 1:n_layers, n) = atm_k(l, 1)%Temperature(1:n_layers)
        humidityJacobian(l, 1:n_layers, n) = atm_k(l, 1)%Absorber(1:n_layers, 1)
        ozoneJacobian(l, 1:n_layers, n) = atm_k(l, 1)%Absorber(1:n_layers, 2)
        outTransmission(l, 1:n_layers, n) = rts(l, 1)%Layer_Optical_Depth
    enddo
    outTb(:,n) = rts(:,1)%Brightness_Temperature 
    emissivity(:,n) = rts(:,1)%Surface_Emissivity
    CALL CRTM_Atmosphere_Destroy(atm)
    DEALLOCATE(atm_k, STAT = alloc_stat)
    DEALLOCATE(rts_K, sfc_k, STAT = alloc_stat)
    DEALLOCATE(rts, STAT = alloc_stat)
    DEALLOCATE(atm, sfc, STAT = alloc_stat)
    ! ==========================================================================

  END DO Profile_Loop
  !!$omp end parallel do
  ! ==========================================================================
  ! 10. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
    STOP
  END IF
  ! ==========================================================================
end subroutine wrap_k_matrix

end module pycrtm
