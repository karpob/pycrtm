module pycrtm 
contains
subroutine wrap_forward( coefficientPath, sensor_id_in, & 
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, &
                        year, month, day, & 
                        nChan, N_Profiles, N_LAYERS, N_aerosols, N_clouds, &
                        pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers, & 
                        co2ConcLayers, & 
                        aerosolEffectiveRadius, aerosolConcentration, aerosolType, & 
                        cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m, n_absorbers, & 
                        landType, soilType, vegType, waterType, snowType, iceType, nthreads, &  
                        outTb, outTransmission, & 
                        emissivityReflectivity )      

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
  integer, intent(in) :: nChan, N_Profiles, N_Layers, N_aerosols, N_clouds
  real(kind=8), intent(in) :: zenithAngle(n_profiles), scanAngle(n_profiles) 
  real(kind=8), intent(in) :: azimuthAngle(n_profiles), solarAngle(n_profiles,2)
  integer, intent(in) :: year(n_profiles), month(n_profiles), day(n_profiles) 
  real(kind=8), intent(in) :: pressureLevels(N_profiles, N_LAYERS+1)
  real(kind=8), intent(in) :: pressureLayers(N_profiles, N_LAYERS), temperatureLayers(N_Profiles,N_Layers)
  real(kind=8), intent(in) ::humidityLayers(N_profiles,N_LAYERS), ozoneConcLayers(N_profiles,N_LAYERS)
  real(kind=8), intent(in) :: co2ConcLayers(N_Profiles,N_layers)
  real(kind=8), intent(in) :: aerosolEffectiveRadius(N_Profiles,N_layers, N_aerosols)
  real(kind=8), intent(in) :: aerosolConcentration(N_profiles,N_layers, N_aerosols)
  real(kind=8), intent(in) :: cloudEffectiveRadius(N_Profiles,N_layers, N_clouds)
  real(kind=8), intent(in) :: cloudConcentration(N_profiles, N_LAYERS, N_clouds) 
  real(kind=8), intent(in) :: cloudFraction(N_Profiles, N_layers)
  integer, intent(in) :: aerosolType(N_Profiles, N_aerosols), cloudType(N_Profiles, N_clouds)
  integer, intent(in) :: n_absorbers(N_Profiles), climatology(N_profiles)
  real(kind=8), intent(in) :: surfaceTemperatures(N_Profiles,4), surfaceFractions(N_profiles, 4)
  real(kind=8), intent(in) :: LAI(N_Profiles), salinity(N_Profiles),  windSpeed10m(N_Profiles), windDirection10m(N_Profiles)
  integer, intent(in) ::  landType(N_Profiles), soilType(N_Profiles), vegType(N_Profiles), waterType(N_Profiles)
  integer, intent(in) ::  snowType(N_Profiles), iceType(N_Profiles) 
  integer, intent(in) :: nthreads
  real(kind=8), intent(out) :: outTb(N_Profiles,nChan), emissivityReflectivity(2,N_Profiles,nChan)
  real(kind=8), intent(out) :: outTransmission(N_profiles, nChan, N_Layers)
  character(len=256), dimension(1) :: sensor_id
  ! --------------------------
  ! Some non-CRTM-y Parameters
  ! --------------------------
  CHARACTER(*), PARAMETER :: SUBROUTINE_NAME   = 'wrap_forward'
  CHARACTER(*), PARAMETER :: PROGRAM_VERSION_ID = '0.01'



  ! ============================================================================
  ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
  !
  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 1
  ! ============================================================================
  


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: message, version
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels, N_clouds_crtm, N_aerosols_crtm
  INTEGER :: i, l, m, n, nc, ll,mm, nn, species
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
  TYPE(CRTM_Atmosphere_type),allocatable  :: atm(:)
  TYPE(CRTM_Surface_type), allocatable    :: sfc(:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
  type(crtm_options_type)                 :: options
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
  ! Karpowicz addition... if we have less than 0 for aerosol/cloud, don't turn on aerosols/clouds.
  if( all(aerosolType < 0) ) then
    N_AEROSOLS_crtm = 0
    aerosolsOn = .False.
  else
    N_AEROSOLS_crtm = N_aerosols
    aerosolsOn = .True. 
  endif

  if( all(cloudType < 0) ) then
    N_CLOUDS_crtm = 0
    cloudsOn = .False.
  else
    N_CLOUDS_crtm = N_clouds
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
  ! Begin loop over profile
  ! ----------------------
  !$ call omp_set_num_threads(nthreads)
  !$omp parallel do default(private) shared(chinfo,emissivityReflectivity,outTb,outTransmission)& 
  !$omp& shared(zenithAngle, scanAngle, azimuthAngle, solarAngle)&
  !$omp& shared(nChan, N_Profiles, N_LAYERS, N_Clouds_crtm, N_aerosols_crtm)&
  !$omp& shared(pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers)& 
  !$omp& shared(co2ConcLayers)& 
  !$omp& shared(aerosolEffectiveRadius, aerosolConcentration, aerosolType, cloudsOn, aerosolsOn)& 
  !$omp& shared(cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology)& 
  !$omp& shared(surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m, n_absorbers)& 
  !$omp& shared(landType, soilType, vegType, waterType, snowType, iceType, year, month, day)&
  !$omp& num_threads(nthreads) 
  Profile_Loop: DO n = 1, N_Profiles
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))

    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
     ALLOCATE( rts( n_channels, 1), STAT = alloc_stat )

    IF ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF
    ALLOCATE( atm( 1), STAT = alloc_stat )

    IF ( alloc_stat /= 0 ) THEN
      message = 'Error allocating atm arrays'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF
  ALLOCATE( sfc(1), STAT = alloc_stat )

    IF ( alloc_stat /= 0 ) THEN
      message = 'Error allocating sfc arrays'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF



    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS(n), N_CLOUDS_crtm, N_AEROSOLS_crtm )
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
    atm(1)%Level_Pressure = pressureLevels(n,:)
    atm(1)%Pressure = pressureLayers(n,:)
    atm(1)%Temperature = temperatureLayers(n,:)
    atm(1)%Absorber(:,1) = humidityLayers(n,:)
    atm(1)%Absorber(:,2) = ozoneConcLayers(n,:)
    if( aerosolsOn )  then
      do species = 1, N_aerosols_crtm
        atm(1)%Aerosol(species)%Type                = aerosolType(n, species)
        atm(1)%Aerosol(species)%Effective_Radius(:) = aerosolEffectiveRadius(n,:, species)
        atm(1)%Aerosol(species)%Concentration(:)    = aerosolConcentration(n,:, species)
      enddo
    endif
    if( cloudsOn ) then
      do species = 1, N_clouds_crtm
        atm(1)%Cloud(species)%Type                = cloudType(n, species)
        atm(1)%Cloud(species)%Effective_Radius(:) = cloudEffectiveRadius(n,:, species)
        atm(1)%Cloud(species)%Water_Content(:)    = cloudConcentration(n,:, species)
      enddo
      atm(1)%Cloud_Fraction(:)            = cloudFraction(n,:)
    endif

    if(n_absorbers(n) >2) then 
      atm(1)%Absorber(:,3)     = co2ConcLayers(n,:)
      atm(1)%Absorber_Id(3)    = CO2_ID
      atm(1)%absorber_units(3) = VOLUME_MIXING_RATIO_UNITS
    endif


    ! 6b. Geometry input
    ! ------------------
    CALL CRTM_Geometry_SetValue( geo, &
                                 year = year(n), & 
                                 month = month(n), & 
                                 day = day(n), & 
                                 Sensor_Zenith_Angle  = zenithAngle(n),   &
                                 Sensor_Scan_Angle    = scanAngle(n),     & 
                                 Sensor_Azimuth_Angle = azimuthAngle(n),  &  
                                 Source_Zenith_Angle  = solarAngle(n,1),  & 
                                 Source_Azimuth_Angle = solarAngle(n,2) )
 
    ! ==========================================================================
    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    sfc%Land_Coverage     = surfaceFractions(n,1)
    sfc%Land_Type         = landType(n) !TUNDRA_SURFACE_TYPE
    sfc%Land_Temperature  = surfaceTemperatures(n,1)
    sfc%Lai               = LAI(n)
    sfc%Soil_Type         = soilType(n) !COARSE_SOIL_TYPE
    sfc%Vegetation_Type   = vegType(n) !GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc%Water_Coverage    = surfaceFractions(n,2)
    sfc%Water_Type        = waterType(n) !SEA_WATER_TYPE
    sfc%Water_Temperature = surfaceTemperatures(n,2)

    ! ...Snow coverage characteristics
    sfc%Snow_Coverage    = surfaceFractions(n,3)
    sfc%Snow_Type        = snowType(n) !FRESH_SNOW_TYPE
    sfc%Snow_Temperature = surfaceTemperatures(n,3)
    ! ...Ice surface characteristics
    sfc%Ice_Coverage    = surfaceFractions(n,4)
    sfc%Ice_Type        = iceType(n) !FRESH_ICE_TYPE
    sfc%Ice_Temperature = surfaceTemperatures(n,4)

    sfc%Wind_Speed = windSpeed10m(n)
    sfc%Wind_Direction = windDirection10m(n)
    sfc%Salinity = salinity(n)
   
    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    !WRITE( *, '( /5x, "Calling the CRTM functions for ",a,"..." )' ) TRIM(sensor_id(1))
    
    ! 8a. The forward model
    ! ---------------------

    ! Need this to get transmission out of solution, otherwise won't be allocated !!!
    call crtm_rtsolution_create( rts, n_layers )
    if ( any(.not. crtm_rtsolution_associated( rts )) ) then
        call display_message( subroutine_name, 'error allocating rts', err_stat)
        !return
    end if

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
    do l=1,nChan
       outTransmission(n, l,1:n_layers) = rts(l,1)%Layer_Optical_Depth
    enddo
    emissivityReflectivity(1,n,:) = rts(:,1)%Surface_Emissivity 
    emissivityReflectivity(2,n,:) = rts(:,1)%Surface_Reflectivity 
    outTb(n,:) = rts(:,1)%Brightness_Temperature 
    
    ! ==========================================================================
    ! STEP 9. **** CLEAN UP FOR NEXT PROFILE ****
    !
    ! 9a. Deallocate the structures
    ! -----------------------------


    ! 9b. Deallocate the arrays
    ! -------------------------
    ! ==========================================================================
    CALL CRTM_Atmosphere_Destroy(atm)
    call crtm_rtsolution_destroy(rts)
    deallocate(atm,sfc)
    DEALLOCATE(rts, STAT = alloc_stat)
  END DO Profile_Loop
  !$omp end parallel do

  
  
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
                        year, month, day, & 
                        nChan, N_profiles, N_LAYERS, N_aerosols, N_clouds, & 
                        pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers, & 
                        co2ConcLayers, & 
                        aerosolEffectiveRadius, aerosolConcentration, aerosolType, & 
                        cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, salinity, windSpeed10m, windDirection10m, n_absorbers, & 
                        landType, soilType, vegType, waterType, snowType, iceType, &  
                        outTb, outTransmission, & 
                        temperatureJacobian, humidityJacobian, ozoneJacobian, co2Jacobian, emissivityReflectivity, nthreads )      

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
  integer, intent(in) :: nChan, N_profiles, N_Layers, N_aerosols, N_clouds 
  ! The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  real(kind=8), intent(in) :: zenithAngle(N_profiles), scanAngle(N_profiles)
  real(kind=8), intent(in) :: azimuthAngle(N_profiles), solarAngle(N_profiles,2)
  integer, intent(in) :: year(n_profiles), month(n_profiles), day(n_profiles)
  real(kind=8), intent(in) :: pressureLevels(N_profiles, N_Layers+1)
  real(kind=8), intent(in) :: pressureLayers(N_profiles, N_layers), temperatureLayers(N_profiles, N_layers)
  real(kind=8), intent(in) :: humidityLayers(N_profiles, N_layers)
  real(kind=8), intent(in) :: ozoneConcLayers(N_profiles, N_layers)
  real(kind=8), intent(in) :: co2ConcLayers(N_profiles, N_layers)
  real(kind=8), intent(in) :: aerosolEffectiveRadius(N_profiles,N_layers, N_aerosols)
  real(kind=8), intent(in) :: aerosolConcentration(N_profiles,N_layers, N_aerosols)
  real(kind=8), intent(in) :: cloudEffectiveRadius(N_profiles,N_layers, N_clouds) 
  real(kind=8), intent(in) :: cloudConcentration(N_profiles, N_layers, N_clouds), cloudFraction(N_profiles,N_layers)
  integer, intent(in) :: aerosolType(N_profiles, N_aerosols), cloudType(N_profiles, N_clouds)
  integer, intent(in) :: n_absorbers(N_profiles), climatology(N_profiles)
  real(kind=8), intent(in) :: surfaceTemperatures(N_profiles,4), surfaceFractions(N_profiles,4), LAI(N_profiles) 
  real(kind=8), intent(in) :: salinity(N_profiles), windSpeed10m(N_profiles), windDirection10m(N_profiles)
  integer, intent(in) :: landType(N_profiles), soilType(N_profiles), vegType(N_profiles), waterType(N_profiles) 
  integer, intent(in) :: snowType(N_profiles), iceType(N_profiles) 
  real(kind=8), intent(out) :: outTb(N_profiles,nChan), emissivityReflectivity(2,N_profiles,nChan)
  real(kind=8), intent(out) :: outTransmission(N_profiles, nChan, N_LAYERS) 
  real(kind=8), intent(out) :: temperatureJacobian(N_profiles, nChan, N_LAYERS)
  real(kind=8), intent(out) ::  humidityJacobian(N_profiles, nChan, N_LAYERS)
  real(kind=8), intent(out) :: ozoneJacobian(N_profiles, nChan, N_LAYERS)
  real(kind=8), intent(out) :: co2Jacobian(N_profiles, nChan, N_LAYERS)
  integer, intent(in) :: nthreads

  character(len=256) :: sensor_id(1)
  ! ============================================================================
  ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
  !
  
  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 1
  ! ============================================================================
  


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: message, version
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels, N_aerosols_crtm, N_clouds_crtm
  INTEGER :: l, m, n, nc, species
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
  if( all(aerosolType < 0) ) then
    N_AEROSOLS_crtm = 0
    aerosolsOn = .False.
  else
    N_AEROSOLS_crtm = N_aerosols
    aerosolsOn = .True. 
  endif

  if( all(cloudType < 0) ) then
    N_CLOUDS_crtm = 0
    cloudsOn = .False.
  else
    N_CLOUDS_crtm = N_clouds
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
  !$ call omp_set_num_threads(nthreads)
  !$omp parallel do default(private) shared(emissivityReflectivity,outTb)&
  !$omp& shared(temperatureJacobian,humidityJacobian)& 
  !$omp& shared(ozoneJacobian,co2Jacobian,outTransmission)&
  !$omp& shared(nChan, N_Layers,N_Absorbers,N_CLOUDS_crtm, N_AEROSOLS_crtm)&
  !$omp& shared(pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers)& 
  !$omp& shared(co2ConcLayers, cloudsOn, aerosolsOn, zenithAngle,scanAngle,azimuthAngle,solarAngle)& 
  !$omp& shared(aerosolEffectiveRadius, aerosolConcentration, aerosolType)& 
  !$omp& shared(cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology)& 
  !$omp& shared(surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m)& 
  !$omp& shared(landType, soilType, vegType, waterType, snowType, iceType)&
  !$omp& shared(sensor_id,coefficientPath, chinfo, year, month, day)&
 
  !$omp& num_threads(nthreads) 
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
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS(n), N_CLOUDS_crtm, N_AEROSOLS_crtm )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF

    ! The output K-MATRIX structure
    CALL CRTM_Atmosphere_Create( atm_K, N_LAYERS, N_ABSORBERS(n), N_CLOUDS_crtm, N_AEROSOLS_crtm )
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
    atm(1)%Level_Pressure = pressureLevels(n,:)
    atm(1)%Pressure = pressureLayers(n,:)
    atm(1)%Temperature = temperatureLayers(n,:)
    atm(1)%Absorber(:,1) = humidityLayers(n,:)
    atm(1)%Absorber(:,2) = ozoneConcLayers(n,:)
    if( aerosolsOn )  then
        do species = 1,N_aerosols_crtm
            atm(1)%Aerosol(species)%Type                = aerosolType(n,species)
            atm(1)%Aerosol(species)%Effective_Radius(:) = aerosolEffectiveRadius( n, :, species)
            atm(1)%Aerosol(species)%Concentration(:)    = aerosolConcentration(n, :, species)
        enddo
    endif
    if( cloudsOn ) then
        do species = 1,N_clouds_crtm
            atm(1)%Cloud(species)%Type                = cloudType(n, species)
            atm(1)%Cloud(species)%Effective_Radius(:) = cloudEffectiveRadius(n,:, species)
            atm(1)%Cloud(species)%Water_Content(:)    = cloudConcentration(n,:, species)
        enddo
        atm(1)%Cloud_Fraction(:)            = cloudFraction(n, :)
    
    endif

    if( n_absorbers(n) > 2 ) then 
        atm(1)%Absorber(:,3)      = co2ConcLayers(n,:)
        atm(1)%Absorber_Id(3)     = CO2_ID
        atm(1)%absorber_units(3)  = VOLUME_MIXING_RATIO_UNITS
    endif
    ! 6b. Geometry input
    ! ------------------
    ! All profiles are given the same value
    CALL CRTM_Geometry_SetValue( geo, &
                                 year = year(n), & 
                                 month = month(n), & 
                                 day = day(n), & 
                                 Sensor_Zenith_Angle  = zenithAngle(n),   &
                                 Sensor_Scan_Angle    = scanAngle(n),     & 
                                 Sensor_Azimuth_Angle = azimuthAngle(n),  &  
                                 Source_Zenith_Angle  = solarAngle(n,1),  & 
                                 Source_Azimuth_Angle = solarAngle(n,2) )
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
    sfc%Land_Coverage     = surfaceFractions(n,1)
    sfc%Land_Type         = landType(n) !TUNDRA_SURFACE_TYPE
    sfc%Land_Temperature  = surfaceTemperatures(n,1)
    sfc%Lai               = LAI(n)
    sfc%Soil_Type         = soilType(n) !COARSE_SOIL_TYPE
    sfc%Vegetation_Type   = vegType(n) !GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc%Water_Coverage    = surfaceFractions(n,2)
    sfc%Water_Type        = waterType(n) !SEA_WATER_TYPE
    sfc%Water_Temperature = surfaceTemperatures(n,2)

    ! ...Snow coverage characteristics
    sfc%Snow_Coverage    = surfaceFractions(n,3)
    sfc%Snow_Type        = snowType(n) !FRESH_SNOW_TYPE
    sfc%Snow_Temperature = surfaceTemperatures(n,3)
    ! ...Ice surface characteristics
    sfc%Ice_Coverage    = surfaceFractions(n,4)
    sfc%Ice_Type        = iceType(n) !FRESH_ICE_TYPE
    sfc%Ice_Temperature = surfaceTemperatures(n,4)
    
    sfc%Wind_Speed = windSpeed10m(n)
    sfc%Wind_Direction = windDirection10m(n)
    sfc%Salinity = salinity(n)
    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    ! 8b. The K-matrix model
    ! ----------------------

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
        temperatureJacobian(n, l, 1:n_layers) = atm_k(l, 1)%Temperature(1:n_layers)
        humidityJacobian(n, l, 1:n_layers) = atm_k(l, 1)%Absorber(1:n_layers, 1)
        ozoneJacobian(n, l, 1:n_layers) = atm_k(l, 1)%Absorber(1:n_layers, 2)
        if(n_absorbers(n)>2) then
            co2Jacobian(n,l, 1:n_layers) = atm_k(l,1)%Absorber(1:n_layers,3)
        endif
        outTransmission(n, l, 1:n_layers) = rts(l, 1)%Layer_Optical_Depth
    enddo
    outTb(n,:) = rts(:,1)%Brightness_Temperature 
    emissivityReflectivity(1,n,:) = rts(:,1)%Surface_Emissivity
    emissivityReflectivity(2,n,:) = rts(:,1)%Surface_Reflectivity
    CALL CRTM_Atmosphere_Destroy(atm)
      if(alloc_stat /= SUCCESS) then
        print*, 'atm destroy failed'
        STOP
    endif
    CALL CRTM_Atmosphere_Destroy(atm_k)
    DEALLOCATE(atm_k, STAT = alloc_stat)
    if(alloc_stat /= SUCCESS) then
        print*, 'atm_k destroy failed'
        STOP
    endif
  DEALLOCATE(rts_K, sfc_k, STAT = alloc_stat)
   if(alloc_stat /= SUCCESS) then
        print*, 'rts_k dealloc failed'
        STOP
    endif
   DEALLOCATE(rts, STAT = alloc_stat)
   if(alloc_stat /= SUCCESS) then
        print*, 'rts dealloc failed'
        STOP
    endif
   DEALLOCATE(atm, sfc, STAT = alloc_stat)
    if(alloc_stat /= SUCCESS) then
        print*, 'dealloc atm,sfc failed'
        STOP
    endif
    ! ==========================================================================

  END DO Profile_Loop
  !$omp end parallel do
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
!--------------------------------------------------------------------------------
!
! NAME:
!       LayerAvg
!
! PURPOSE:
!    Given px1 (output domain) and px2 (input domain) that are ascending arrays,
!    it computes the accumulated weighting factors for interpolation from input
!    to output domainsan, and the index array, interp_index, such that output
!    variable
!      to(jo) = sum(pz(interp_index(1,jo):interp_index(2,jo),jo) &
!                   *(ti(interp_index(1,jo):interp_index(2,jo))))

!
! CALLING SEQUENCE:
!       CALL LayerAvg(PX1,PX2,PZ,Interp_index)
!
! INPUT ARGUMENTS:
!       PX1:          The abscissa values for the target data (output domain) and
!                     they must be monotonically ascending (e.g. lnP; in increasing values).
!                     UNITS:      N/A
!                     TYPE:       fp
!                     DIMENSION:  rank-1 (KN1)
!                     ATTRIBUTES: INTENT(IN)
!       PX2:          The abscissa values for the source data (input domain) and
!                     they must be monotonically ascending (e.g. lnP; in increasing values).
!                     UNITS:      N/A
!                     TYPE:       fp
!                     DIMENSION:  rank-1 (KN2)
!                     ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       PZ:           Resultant accumulated weighting factors for
!                     interpolation from input to output domains.
!                     UNITS:      N/A
!                     TYPE:       fp
!                     DIMENSION:  rank-2 (KN2,KN1)
!                     ATTRIBUTES: INTENT(OUT)
!    interp_index:    The index array of dimension (2 x KN1) for
!                     Start index for relevant PZ row array segment and
!                     End index for relevant PZ row array segment, where KN1 = SIZE(PX1)
!                     UNITS:      N/A
!                     TYPE:       Integer
!                     DIMENSION:  rank-2
!                     ATTRIBUTES: INTENT(IN)
!
!     Comments:
!
!     1) PX1(i)<PX1(i+1) & PX2(i)<PX2(i+1)
!
! RESTRICTIONS:
!     To be efficient, this routine does not check that x and u are both
!     monotonically ascending and the index bounds.
!
!     - Journal reference:
!       Rochon, Y.J., L. Garand, D.S. Turner, and S. Polavarapu.
!       Jacobian mapping between vertical coordinate systems in data assimilation,
!       Q. J. of the Royal Met. Soc., 133, 1547-1558 (2007) DOI: 10.1002/qj.117
!
! CREATION HISTORY:
!       Written by:     Yong Chen, 08-May-2009
!-----------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  SUBROUTINE LayerAvg( PX1,PX2,PZ,interp_index)
    REAL(kind=8),     DIMENSION(:),       INTENT(IN)     :: PX1(:)
    REAL(kind=8),     DIMENSION(:),       INTENT(IN)     :: PX2(:)
    REAL(kind=8),     DIMENSION(:,:),     INTENT(OUT)    :: PZ(:, :)
    INTEGER,      DIMENSION(:,:),     INTENT(OUT)    :: interp_index

    ! Local variables
    INTEGER  ::  KN1,KN2, ibot,itop,ii,istart,iend, J,IC,ISKIP,KI
    REAL(kind=8) ::  z1,z2,z3,zw1,zw2,zsum
    REAL(kind=8) ::  y1,y2,d,w10,w20,dz,dx,dy,dzd,dxd

  ! for using in SUBROUTINE LayerAvg
    REAL(kind=8), PARAMETER :: SMALLDIFF  = 1.0E-20
    REAL(kind=8), PARAMETER :: ONE =1.0d0
    REAL(kind=8), PARAMETER :: ZERO =0.0d0
    KN1 = SIZE(PX1)
    KN2 = SIZE(PX2)

    istart=1
    iend=kn1
    DO KI=1,KN1
      z2=px1(ki)
!
      if (ki == 1) then
         z1=2.0*z2-px1(ki+1)
      else
         z1=px1(ki-1)
      endif
!
      if (ki == kn1) then
         z3=2.0*z2-z1
      else
         z3=px1(ki+1)
      endif
      if (z3 > px2(kn2)) z3=px2(kn2)
!
      iskip=0
      if (z2 >= px2(kn2)) then
         z3=px2(kn2)
         z2=px2(kn2)
         iskip=1
      endif

! --- Determine forward interpolator
!
      pz(1:kn2,ki)=ZERO
      ic=0
      do j=istart,kn2-1
        if (px2(j) > z3) go to 1000
!
        if (px2(j) <= z2 .and. px2(j+1) > z1) then
          itop=0
          ibot=0
          if (z1 < z3) then
             y1=z1
             if (px2(j) > z1) then
                y1=px2(j)
                itop=1
             endif
             y2=z2
             if (px2(j+1) < z2) then
                y2=px2(j+1)
                ibot=1
             endif
          else
             y1=z2
             if (px2(j) > z2) then
                y1=px2(j)
                itop=1
             endif
             y2=z1
             if (px2(j+1) < z1) then
                y2=px2(j+1)
                ibot=1
             endif
          endif
!
! ---     Set weights for forward interpolator
!
          dy=y2-y1
          dz=z1-z2
          if (abs(dz) < SMALLDIFF) then
             write(6,*) 'SUBLAYER: ERROR: dz is <=0. dz = ',dz
             write(6,*) 'z1,z2,z3 = ',z1,z2,z3
             write(6,*) 'px2(j),px2(j+1)    = ',px2(j),px2(j+1)
             return
          else
             dzd=ONE/dz
          endif
          zw1=(z1-y1)*dzd*dy
          zw2=(z1-y2)*dzd*dy
          w10=zw1
          w20=zw2
          dx=(px2(j+1)-px2(j))
          if (abs(dx) < SMALLDIFF) then
             write(6,*) 'SUBLAYER: ERROR: dx is <=0. dx = ',dx
             write(6,*) 'z1,z2,z3 = ',z1,z2,z3
             write(6,*) 'px2(j),px2(j+1)    = ',px2(j),px2(j+1)
             return
          else
             dxd=ONE/dx
          endif
!
          d=(px2(j+1)-z2)*dxd
          if (z1 < z3 .and. ibot == 0) then
             zw1=zw1+zw2*d
             zw2=zw2*(ONE-d)
          else if (z1 > z3 .and. itop == 0) then
             zw2=zw2+zw1*(ONE-d)
             zw1=zw1*d
          end if
          pz(j,ki)=pz(j,ki)+zw1
          pz(j+1,ki)=pz(j+1,ki)+zw2
          ic=1
        endif
!
        if (px2(j) < z3 .and. px2(j+1) >= z2 .and. iskip == 0) then
          itop=0
          ibot=0
          if (z3 < z1) then
             y1=z3
             if (px2(j) > z3) then
                y1=px2(j)
                itop=1
             endif
             y2=z2
             if (px2(j+1) < z2) then
                y2=px2(j+1)
                ibot=1
             endif
          else
             y1=z2
             if (px2(j) > z2) then
                y1=px2(j)
                itop=1
             endif
             y2=z3
             if (px2(j+1) < z3) then
                y2=px2(j+1)
                ibot=1
             endif
          endif
!
! ---     Set weights for forward interpolator
!
          dy=y2-y1
          dz=z3-z2
          if (abs(dz) < SMALLDIFF) then
             write(6,*) 'SUBLAYER: ERROR: dz is <=0. dz = ',dz
             write(6,*) 'z3,z2,z1 = ',z3,z2,z1
             write(6,*) 'px2(j),px2(j+1)    = ',px2(j),px2(j+1)
             return
          else
             dzd=ONE/dz
          endif
          zw1=(z3-y1)*dzd*dy
          zw2=(z3-y2)*dzd*dy
          w10=zw1
          w20=zw2
          dx=(px2(j+1)-px2(j))
          if (abs(dx) < SMALLDIFF) then
             write(6,*) 'SUBLAYER: ERROR: dx is <=0. dx = ',dx
             write(6,*) 'z3,z2,z1 = ',z3,z2,z1
             write(6,*) 'px2(j),px2(j+1)    = ',px2(j),px2(j+1)
             return
          else
             dxd=ONE/dx
          endif
!
          d=(px2(j+1)-z2)*dxd
          if (z3 < z1 .and. ibot == 0) then
             zw1=zw1+zw2*d
             zw2=zw2*(ONE-d)
          else if (z3 > z1 .and. itop == 0) then
             zw2=zw2+zw1*(ONE-d)
             zw1=zw1*d
          end if
          pz(j,ki)=pz(j,ki)+zw1
          pz(j+1,ki)=pz(j+1,ki)+zw2
          ic=1
        endif
      enddo
      j=kn2
 1000 continue
      if (ic == 0) pz(j,ki)=ONE
!
!     Normalize sum to unity (instead of calculating and dividing by
!     weighting denominator)
!
      do ii=istart,kn2
        if(PZ(ii,ki) /= ZERO)then
          interp_index(1, ki)=ii
          exit
        endif
      enddo
      istart=interp_index(1, ki)

      interp_index(2,ki)=kn2

      do ii=interp_index(1, ki)+1,kn2
        if(PZ(ii,ki) == ZERO)then
          interp_index(2,ki)=ii-1
          exit
        endif
      enddo
      iend=interp_index(2,ki)

      zsum=sum(pz(interp_index(1, ki):interp_index(2,ki),ki))
      pz(interp_index(1,ki):interp_index(2,ki),ki)= &
                 pz(interp_index(1,ki):interp_index(2, ki),ki)/zsum
    ENDDO


  END SUBROUTINE LayerAvg

  subroutine applyAvg( Ref_LnPressure, User_LnPressure, Nref, Nuser, Xin, Xout ) 
    integer,         intent(in) :: Nref, Nuser
    REAL(kind=8),   INTENT(IN)     :: Ref_LnPressure(Nref)
    REAL(kind=8),   INTENT(IN)     :: User_LnPressure(Nuser)
    REAL(kind=8),   INTENT(IN)     :: Xin(Nuser)
    REAL(kind=8),   INTENT(OUT)    :: Xout(Nref)
    !locals   
    INTEGER      :: k, interp_index(2,Nref)
    REAL(kind=8) :: Acc_Weighting(Nuser,Nref)
 
    CALL LayerAvg( Ref_LnPressure   , &
                   User_LnPressure  , &
                   Acc_Weighting    , &
                   interp_index)



    DO k = 1, Nref
       Xout(k) = SUM(Acc_Weighting(interp_index(1,k):interp_index(2,k), k)  &
                           * Xin(interp_index(1,k):interp_index(2,k)) )
    END DO
  end subroutine applyAvg

end module pycrtm
