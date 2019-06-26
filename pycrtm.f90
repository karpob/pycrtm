module pycrtm 
contains
subroutine wrap_forward( coefficientPath, sensor_id_in, & 
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, nChan, &
                        N_LAYERS, pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers, & 
                        co2ConcLayers, & 
                        aerosolEffectiveRadius, aerosolConcentration, aerosolType, & 
                        cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, windSpeed10m, windDirection10m, n_absorbers, & 
                        landType, soilType, vegType, waterType, snowType, iceType, &  
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
  real(kind=8), intent(in) :: zenithAngle, scanAngle, azimuthAngle, solarAngle
  integer, intent(in) :: nChan, N_Layers
  real(kind=8), intent(in) :: pressureLevels(N_LAYERS+1)
  real(kind=8), intent(in) :: pressureLayers(N_LAYERS), temperatureLayers(N_LAYERS), humidityLayers(N_LAYERS)
  real(kind=8), intent(in) :: ozoneConcLayers(N_LAYERS)
  real(kind=8), intent(in) :: co2ConcLayers(N_LAYERS)
  real(kind=8), intent(in) :: aerosolEffectiveRadius(N_LAYERS), aerosolConcentration(N_LAYERS)
  real(kind=8), intent(in) :: cloudEffectiveRadius(N_LAYERS), cloudConcentration(N_LAYERS), cloudFraction(N_LAYERS)
  integer, intent(in) :: aerosolType, cloudType, n_absorbers, climatology
  real(kind=8), intent(in) :: surfaceTemperatures(4), surfaceFractions(4), LAI, windSpeed10m, windDirection10m
  integer, intent(in) ::  landType, soilType, vegType, waterType, snowType, iceType 
  real(kind=8), intent(out) :: outTb(nChan), emissivity(nChan)
  real(kind=8), intent(out) :: outTransmission(nChan,N_LAYERS)
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
  INTEGER, PARAMETER :: N_PROFILES  = 1
  INTEGER, PARAMETER :: N_CLOUDS    = 1 
  INTEGER, PARAMETER :: N_AEROSOLS  = 1
  
  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 1
  !CHARACTER(len=20) :: sensor_id 
  ! ============================================================================
  


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: message, version
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels
  INTEGER :: l, m, n, nc, ll,mm, nn
  real, dimension(n_layers,nchan) :: outTau


  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: geo(N_PROFILES)


  ! 3b. Define the FORWARD variables
  ! --------------------------------
  TYPE(CRTM_Atmosphere_type)              :: atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
  type(crtm_options_type)     , dimension(N_PROFILES)        :: options
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
  WRITE( *,'(/5x,"Initializing the CRTM...")' )
  err_stat = CRTM_Init( sensor_id,  chinfo, File_Path=coefficientPath, Quiet=.TRUE.)

  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error initializing CRTM'
    CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
    STOP
  END IF

  ! 4b. Output some channel information
  ! -----------------------------------
  n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels
  DO n = 1, N_SENSORS
    WRITE( *,'(7x,i0," from ",a)' ) &
      CRTM_ChannelInfo_n_Channels(chinfo(n)), TRIM(sensor_id(n))
  END DO
  ! ============================================================================



  ! Begin loop over sensors
  ! ----------------------
  Sensor_Loop: DO n = 1, N_SENSORS

  
    ! ==========================================================================
    ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
    !
    ! 5a. Determine the number of channels
    !     for the current sensor
    ! ------------------------------------
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(n))

    
    ! 5b. Allocate the ARRAYS
    ! -----------------------
    ALLOCATE( rts( n_channels, N_PROFILES ), STAT = alloc_stat )

    IF ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF


    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
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
    atm(1)%Climatology = climatology
    atm(1)%Level_Pressure = pressureLevels
    atm(1)%Pressure = pressureLayers
    atm(1)%Temperature = temperatureLayers
    atm(1)%Absorber(:,1) = humidityLayers
    atm(1)%Absorber(:,2) = ozoneConcLayers

    atm(1)%Aerosol(1)%Type = aerosolType
    atm(1)%Aerosol(1)%Effective_Radius = aerosolEffectiveRadius
    atm(1)%Aerosol(1)%Concentration = aerosolConcentration

    atm(1)%Cloud(1)%Type = cloudType
    atm(1)%Cloud(1)%Effective_Radius = cloudEffectiveRadius
    atm(1)%Cloud(1)%Water_Content = cloudConcentration
    atm(1)%Absorber(:,3)     = co2ConcLayers
    atm(1)%Cloud_Fraction = cloudFraction

    ! 6b. Geometry input
    ! ------------------
    ! All profiles are given the same value
    !  The Sensor_SCAN_ANGLE is optional.  !! BMK- Oh? this would be nice. Not sure if that's true though. Think you need it for FastEm?
    CALL CRTM_Geometry_SetValue( geo, &
                                 Sensor_Zenith_Angle = zenithAngle, &
                                 Sensor_Scan_Angle   = scanAngle,   & 
                                 Sensor_Azimuth_Angle = azimuthAngle )
    ! ==========================================================================
    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    sfc%Land_Coverage     = surfaceFractions(1)
    sfc%Land_Type         = landType !TUNDRA_SURFACE_TYPE
    sfc%Land_Temperature  = surfaceTemperatures(1)
    sfc%Lai               = LAI
    sfc%Soil_Type         = soilType !COARSE_SOIL_TYPE
    sfc%Vegetation_Type   = vegType !GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc%Water_Coverage    = surfaceFractions(2)
    sfc%Water_Type        = waterType !SEA_WATER_TYPE
    sfc%Water_Temperature = surfaceTemperatures(2)
    !Sfc%Wind_Direction = windDirection10m
    !Sfc%Wind_Speed = windSpeed10m
    !Sfc%Salinity = 0.0_fp   



    ! ...Snow coverage characteristics
    sfc%Snow_Coverage    = surfaceFractions(3)
    sfc%Snow_Type        = snowType !FRESH_SNOW_TYPE
    sfc%Snow_Temperature = surfaceTemperatures(3)
    ! ...Ice surface characteristics
    sfc%Ice_Coverage    = surfaceFractions(4)
    sfc%Ice_Type        = iceType !FRESH_ICE_TYPE
    sfc%Ice_Temperature = surfaceTemperatures(4)

   
    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    WRITE( *, '( /5x, "Calling the CRTM functions for ",a,"..." )' ) TRIM(sensor_id(n))
    
    ! 8a. The forward model
    ! ---------------------

    call crtm_options_create( options, nChan )
    if ( any(.not. crtm_options_associated( options )) ) then 
        call display_message( subroutine_name, 'error allocating options', FAILURE)  
        return
    endif
    ! Need this to get transmission out of solution, otherwise won't be allocated !!!
    call crtm_rtsolution_create( rts, n_layers )
    if ( any(.not. crtm_rtsolution_associated( rts )) ) then
        call display_message( subroutine_name, 'error allocating rts', err_stat)
        return
    end if

    options%Use_Emissivity = .false.
    options%Use_Direct_Reflectivity = .false.

    err_stat = CRTM_Forward( atm        , &  ! Input
                             sfc        , &  ! Input
                             geo        , &  ! Input
                             chinfo, &  ! Input
                             rts,    & ! Output
                            options = options ) 
    IF ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM Forward Model for '//TRIM(sensor_id(n))
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
        outTransmission(l,1:n_layers) = rts(l,1)%Layer_Optical_Depth
    enddo
    emissivity = rts(:,1)%Surface_Emissivity 
    outTb = rts(:,1)%Brightness_Temperature 
    
    ! ==========================================================================
    ! STEP 9. **** CLEAN UP FOR NEXT SENSOR ****
    !
    ! 9a. Deallocate the structures
    ! -----------------------------
    CALL CRTM_Atmosphere_Destroy(atm)


    ! 9b. Deallocate the arrays
    ! -------------------------
    DEALLOCATE(rts, STAT = alloc_stat)
    ! ==========================================================================

  END DO Sensor_Loop


  
  
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
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, nChan, &
                        N_LAYERS, pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers, & 
                        co2ConcLayers, & 
                        aerosolEffectiveRadius, aerosolConcentration, aerosolType, & 
                        cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, windSpeed10m, windDirection10m, n_absorbers, & 
                        landType, soilType, vegType, waterType, snowType, iceType, &  
                        outTb, outTransmission, & 
                        temperatureJacobian, humidityJacobian, ozoneJacobian, emissivity )      

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
  ! The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  real(kind=8), intent(in) :: zenithAngle, scanAngle, azimuthAngle, solarAngle
  integer, intent(in) :: nChan, N_Layers 
  real(kind=8), intent(in) :: pressureLevels(N_LAYERS+1)
  real(kind=8), intent(in) :: pressureLayers(N_LAYERS), temperatureLayers(N_LAYERS), humidityLayers(N_LAYERS)
  real(kind=8), intent(in) :: ozoneConcLayers(N_LAYERS)
  real(kind=8), intent(in) :: co2ConcLayers(N_LAYERS)
  real(kind=8), intent(in) :: aerosolEffectiveRadius(N_LAYERS), aerosolConcentration(N_LAYERS)
  real(kind=8), intent(in) :: cloudEffectiveRadius(N_LAYERS), cloudConcentration(N_LAYERS), cloudFraction(N_LAYERS)
  integer, intent(in) :: aerosolType, cloudType, n_absorbers, climatology
  real(kind=8), intent(in) :: surfaceTemperatures(4), surfaceFractions(4), LAI, windSpeed10m, windDirection10m
  integer, intent(in) :: landType, soilType, vegType, waterType, snowType, iceType 
  real(kind=8), intent(out) :: outTb(nChan), emissivity(nChan)
  real(kind=8), intent(out) :: outTransmission(nChan,N_LAYERS), temperatureJacobian(nChan,N_LAYERS)
  real(kind=8), intent(out) ::  humidityJacobian(nChan,N_LAYERS), ozoneJacobian(nChan, N_LAYERS)
  INTEGER, PARAMETER :: TUNDRA_SURFACE_TYPE         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
  INTEGER, PARAMETER :: SCRUB_SURFACE_TYPE          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
  INTEGER, PARAMETER :: COARSE_SOIL_TYPE            =  1  ! Soil type                for MW land SfcOptics
  INTEGER, PARAMETER :: GROUNDCOVER_VEGETATION_TYPE =  7  ! Vegetation type          for MW Land SfcOptics
  INTEGER, PARAMETER :: BARE_SOIL_VEGETATION_TYPE   = 11  ! Vegetation type          for MW Land SfcOptics
  INTEGER, PARAMETER :: SEA_WATER_TYPE              =  1  ! Water type               for all SfcOptics
  INTEGER, PARAMETER :: FRESH_SNOW_TYPE             =  2  ! NPOESS Snow type         for IR/VIS SfcOptics
  INTEGER, PARAMETER :: FRESH_ICE_TYPE              =  1  ! NPOESS Ice type          for IR/VIS SfcOptics

  character(len=256) :: sensor_id(1)
  ! ============================================================================
  ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
  !
  ! Directory location of coefficients

  ! Profile dimensions
  INTEGER, PARAMETER :: N_PROFILES  = 1
  INTEGER, PARAMETER :: N_CLOUDS    = 1
  INTEGER, PARAMETER :: N_AEROSOLS  = 1
  
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



  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: geo(N_PROFILES)

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
  allocate(atm(N_PROFILES), sfc(N_PROFILES))
  WRITE( *,'(/5x,"Initializing the CRTM...")' )
  err_stat = CRTM_Init( sensor_id, &
                        chinfo, &
                        File_Path=coefficientPath, &
                        !MWwaterCoeff_File = 'FASTEM5.MWwater.EmisCoeff.bin', & 
                        Quiet=.False.)
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error initializing CRTM'
    CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
    STOP
  END IF

  ! 4b. Output some channel information
  ! -----------------------------------
  n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels
  DO n = 1, N_SENSORS
    WRITE( *,'(7x,i0," from ",a)' ) &
      CRTM_ChannelInfo_n_Channels(chinfo(n)), TRIM(sensor_id(n))
  END DO
  ! ============================================================================

    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF
 

  ! Begin loop over sensors
  ! ----------------------
  Sensor_Loop: DO n = 1, N_SENSORS

  
    ! ==========================================================================
    ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
    !
    ! 5a. Determine the number of channels
    !     for the current sensor
    ! ------------------------------------
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(n))

    
    ! 5b. Allocate the ARRAYS
    ! -----------------------
    ALLOCATE( rts( n_channels, N_PROFILES),    & 
              atm_K( n_channels, N_PROFILES ), &
              sfc_K( n_channels, N_PROFILES ), &
              rts_K( n_channels, N_PROFILES ), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF


    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF

    ! The output K-MATRIX structure
    CALL CRTM_Atmosphere_Create( atm_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure'
      CALL Display_Message( SUBROUTINE_NAME, message, FAILURE )
      STOP
    END IF

    call crtm_rtsolution_create( rts, n_layers )
    if ( any(.not. crtm_rtsolution_associated( rts )) ) then
        call display_message( subroutine_name, 'error allocating rts', err_stat)
        return
    end if
    
    call crtm_rtsolution_create( rts_k, n_layers )
    if ( any(.not. crtm_rtsolution_associated( rts_k )) ) then
        call display_message( subroutine_name, 'error allocating rts_k', err_stat)
        return
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
    atm(1)%Climatology         = climatology  
    atm(1)%Absorber_Id(1:3)    = (/ H2O_ID                 , O3_ID, CO2_ID /)
    atm(1)%Absorber_Units(1:3) = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    ! ...Profile data
    atm(1)%Level_Pressure = pressureLevels
    atm(1)%Pressure = pressureLayers
    atm(1)%Temperature = temperatureLayers
    atm(1)%Absorber(:,1) = humidityLayers
    atm(1)%Absorber(:,2) = ozoneConcLayers
    atm(1)%Aerosol(1)%Type = aerosolType
    atm(1)%Aerosol(1)%Effective_Radius = aerosolEffectiveRadius
    atm(1)%Aerosol(1)%Concentration = aerosolConcentration

    atm(1)%Cloud(1)%Type = cloudType
    atm(1)%Cloud(1)%Effective_Radius = cloudEffectiveRadius
    atm(1)%Cloud(1)%Water_Content = cloudConcentration
    atm(1)%Cloud_Fraction = cloudFraction
    atm(1)%Absorber(:,3)     = co2ConcLayers
    ! 6b. Geometry input
    ! ------------------
    ! All profiles are given the same value
    !  The Sensor_SCAN_ANGLE is optional.  !! BMK- Oh? this would be nice. Not sure if that's true though.
    CALL CRTM_Geometry_SetValue( geo, &
                                 Sensor_Zenith_Angle = zenithAngle, &
                                 Sensor_Scan_Angle   = scanAngle ,  &
                                 Sensor_Azimuth_Angle = azimuthAngle )
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
    sfc%Land_Coverage     = surfaceFractions(1)
    sfc%Land_Type         = landType !TUNDRA_SURFACE_TYPE
    sfc%Land_Temperature  = surfaceTemperatures(1)
    sfc%Lai               = LAI
    sfc%Soil_Type         = soilType !COARSE_SOIL_TYPE
    sfc%Vegetation_Type   = vegType !GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc%Water_Coverage    = surfaceFractions(2)
    sfc%Water_Type        = waterType !SEA_WATER_TYPE
    sfc%Water_Temperature = surfaceTemperatures(2)
    !Sfc%Wind_Direction = windDirection10m
    !Sfc%Wind_Speed = windSpeed10m
    !Sfc%Salinity = 0.0_fp   



    ! ...Snow coverage characteristics
    sfc%Snow_Coverage    = surfaceFractions(3)
    sfc%Snow_Type        = snowType !FRESH_SNOW_TYPE
    sfc%Snow_Temperature = surfaceTemperatures(3)
    ! ...Ice surface characteristics
    sfc%Ice_Coverage    = surfaceFractions(4)
    sfc%Ice_Type        = iceType !FRESH_ICE_TYPE
    sfc%Ice_Temperature = surfaceTemperatures(4)

    
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
      message = 'Error calling CRTM K-Matrix Model for '//TRIM(SENSOR_ID(n))
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
        temperatureJacobian(l,1:n_layers) = atm_k(l,1)%Temperature(1:n_layers)
        humidityJacobian(l,1:n_layers) = atm_k(l,1)%Absorber(1:n_layers,1)
        ozoneJacobian(l,1:n_layers) = atm_k(l,1)%Absorber(1:n_layers,2)
        outTransmission(l,1:n_layers) = rts(l,1)%Layer_Optical_Depth
    enddo
    outTb = rts(:,1)%Brightness_Temperature 
    emissivity = rts(:,1)%Surface_Emissivity
    CALL CRTM_Atmosphere_Destroy(atm)
    DEALLOCATE(atm_k, STAT = alloc_stat)
    DEALLOCATE(rts_K, sfc_k, STAT = alloc_stat)
    ! ==========================================================================

  END DO Sensor_Loop
  
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
