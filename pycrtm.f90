module pycrtm 
contains
subroutine wrap_forward( coefficientPath, sensor_id_in, & 
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, nChan, &
                        N_LAYERS, pressureLevels, pressureLayers, temperatureLayers, humidityLayers, ozoneConcLayers, & 
                        co2ConcLayers, & 
                        aerosolEffectiveRadius, aerosolConcentration, aerosolType, & 
                        cloudEffectiveRadius, cloudConcentration, cloudType, & 
                        surfaceType, surfaceTemperature, windSpeed10m, windDirection10m, & 
                        outTb, outTransmission )      

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
  real, intent(in) :: pressureLevels(N_LAYERS+1)
  real, intent(in) :: pressureLayers(N_LAYERS), temperatureLayers(N_LAYERS), humidityLayers(N_LAYERS)
  real, intent(in) :: ozoneConcLayers(N_LAYERS)
  real, intent(in) :: co2ConcLayers(N_LAYERS)
  real, intent(in) :: aerosolEffectiveRadius(N_LAYERS), aerosolConcentration(N_LAYERS)
  real, intent(in) :: cloudEffectiveRadius(N_LAYERS), cloudConcentration(N_LAYERS)
  integer, intent(in) :: surfaceType, aerosolType, cloudType
  real, intent(in) :: surfaceTemperature, windSpeed10m, windDirection10m
  real, intent(out) :: outTb(nChan)
  real, intent(out) :: outTransmission(nChan,N_LAYERS)
  character(len=256), dimension(1) :: sensor_id
  INTEGER, PARAMETER :: TUNDRA_SURFACE_TYPE         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
  INTEGER, PARAMETER :: SCRUB_SURFACE_TYPE          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
  INTEGER, PARAMETER :: COARSE_SOIL_TYPE            =  1  ! Soil type                for MW land SfcOptics
  INTEGER, PARAMETER :: GROUNDCOVER_VEGETATION_TYPE =  7  ! Vegetation type          for MW Land SfcOptics
  INTEGER, PARAMETER :: BARE_SOIL_VEGETATION_TYPE   = 11  ! Vegetation type          for MW Land SfcOptics
  INTEGER, PARAMETER :: SEA_WATER_TYPE              =  1  ! Water type               for all SfcOptics
  INTEGER, PARAMETER :: FRESH_SNOW_TYPE             =  2  ! NPOESS Snow type         for IR/VIS SfcOptics
  INTEGER, PARAMETER :: FRESH_ICE_TYPE              =  1  ! NPOESS Ice type          for IR/VIS SfcOptics



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
  INTEGER, PARAMETER :: N_ABSORBERS = 2
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
    ! 6b. Geometry input
    ! ------------------
    ! All profiles are given the same value
    !  The Sensor_SCAN_ANGLE is optional.  !! BMK- Oh? this would be nice. Not sure if that's true though. Think you need it for FastEm?
    CALL CRTM_Geometry_SetValue( geo, &
                                 Sensor_Zenith_Angle = zenithAngle, &
                                 Sensor_Scan_Angle   = scanAngle ) !,   & 
                                 !Sensor_Azimuth_Angle = dble(azimuthAngle) )
    ! ==========================================================================

    ! surface stuff! need to put something more advanced here!
    !Sfc%Water_Coverage = 1
    Sfc%Water_Temperature = surfaceTemperature
    Sfc%Wind_Direction = windDirection10m
    Sfc%Wind_Speed = windSpeed10m
    Sfc%Salinity = 35.0 
    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    sfc%Land_Coverage     = 0.1_fp
    sfc%Land_Type         = TUNDRA_SURFACE_TYPE
    sfc%Land_Temperature  = 272.0_fp
    sfc%Lai               = 0.17_fp
    sfc%Soil_Type         = COARSE_SOIL_TYPE
    sfc%Vegetation_Type   = GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc%Water_Coverage    = 0.5_fp
    sfc%Water_Type        = SEA_WATER_TYPE
    sfc%Water_Temperature = 275.0_fp
    ! ...Snow coverage characteristics
    sfc%Snow_Coverage    = 0.25_fp
    sfc%Snow_Type        = FRESH_SNOW_TYPE
    sfc%Snow_Temperature = 265.0_fp
    ! ...Ice surface characteristics
    sfc%Ice_Coverage    = 0.15_fp
    sfc%Ice_Type        = FRESH_ICE_TYPE
    sfc%Ice_Temperature = 269.0_fp

   
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
    print *, 'Emissivity Crtm',rts(:,1)%Surface_Emissivity 
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
                        cloudEffectiveRadius, cloudConcentration, cloudType, & 
                        surfaceType, surfaceTemperature, windSpeed10m, windDirection10m, & 
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
  real(kind=8), intent(in) :: cloudEffectiveRadius(N_LAYERS), cloudConcentration(N_LAYERS)
  integer, intent(in) :: surfaceType, aerosolType, cloudType
  real(kind=8), intent(in) :: surfaceTemperature, windSpeed10m, windDirection10m
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
  INTEGER, PARAMETER :: N_ABSORBERS = 2
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
    atm(1)%Climatology         = US_STANDARD_ATMOSPHERE  
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

    atm(1)%Absorber(:,3)     = co2ConcLayers
    WHERE(atm(1)%Cloud(1)%Water_Content > 0.0_fp) atm(1)%Cloud_Fraction = 0.0_fp
    ! 6b. Geometry input
    ! ------------------
    ! All profiles are given the same value
    !  The Sensor_SCAN_ANGLE is optional.  !! BMK- Oh? this would be nice. Not sure if that's true though.
    CALL CRTM_Geometry_SetValue( geo, &
                                 Sensor_Zenith_Angle = 30.0_fp, &
                                 Sensor_Scan_Angle   = 26.37293341421_fp)!,   &
    !                             Sensor_Azimuth_Angle = dble(azimuthAngle) )
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
    ! surface stuff! need to put something more advanced here!
    !Sfc%Water_Coverage = 1
    !Sfc%Water_Temperature = surfaceTemperature
    !Sfc%Wind_Direction = windDirection10m
    !Sfc%Wind_Speed = windSpeed10m
    !Sfc%Salinity = 35.0    


    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    sfc%Land_Coverage     = 0.1_fp
    sfc%Land_Type         = TUNDRA_SURFACE_TYPE
    sfc%Land_Temperature  = 272.0_fp
    sfc%Lai               = 0.17_fp
    sfc%Soil_Type         = COARSE_SOIL_TYPE
    sfc%Vegetation_Type   = GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc%Water_Coverage    = 0.5_fp
    sfc%Water_Type        = SEA_WATER_TYPE
    sfc%Water_Temperature = 275.0_fp
    ! ...Snow coverage characteristics
    sfc%Snow_Coverage    = 0.25_fp
    sfc%Snow_Type        = FRESH_SNOW_TYPE
    sfc%Snow_Temperature = 265.0_fp
    ! ...Ice surface characteristics
    sfc%Ice_Coverage    = 0.15_fp
    sfc%Ice_Type        = FRESH_ICE_TYPE
    sfc%Ice_Temperature = 269.0_fp

    
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
  SUBROUTINE test_data_us_std(Level_Pressure, Pressure, Temperature, water_vapor, ozone, &
                              aerosolEffectiveRadius, aerosolConcentration, aerosolType, &
                              cloudEffectiveRadius, cloudConcentration, cloudType, co2 )
    
    USE CRTM_Module
    IMPLICIT NONE 
    real(kind=8), intent(out), dimension(93) :: Level_Pressure 
    real(kind=8), intent(out), dimension(92) :: Pressure, Temperature, water_vapor, ozone, aerosolEffectiveRadius, &
                                        aerosolConcentration,  cloudEffectiveRadius, cloudConcentration, co2
    integer, intent(out) :: aerosolType, cloudType
    
    ! ...Profile data
    Level_Pressure = &
    (/0.714_fp,   0.975_fp,   1.297_fp,   1.687_fp,   2.153_fp,   2.701_fp,   3.340_fp,   4.077_fp, &
      4.920_fp,   5.878_fp,   6.957_fp,   8.165_fp,   9.512_fp,  11.004_fp,  12.649_fp,  14.456_fp, &
     16.432_fp,  18.585_fp,  20.922_fp,  23.453_fp,  26.183_fp,  29.121_fp,  32.274_fp,  35.650_fp, &
     39.257_fp,  43.100_fp,  47.188_fp,  51.528_fp,  56.126_fp,  60.990_fp,  66.125_fp,  71.540_fp, &
     77.240_fp,  83.231_fp,  89.520_fp,  96.114_fp, 103.017_fp, 110.237_fp, 117.777_fp, 125.646_fp, &
    133.846_fp, 142.385_fp, 151.266_fp, 160.496_fp, 170.078_fp, 180.018_fp, 190.320_fp, 200.989_fp, &
    212.028_fp, 223.441_fp, 235.234_fp, 247.409_fp, 259.969_fp, 272.919_fp, 286.262_fp, 300.000_fp, &
    314.137_fp, 328.675_fp, 343.618_fp, 358.967_fp, 374.724_fp, 390.893_fp, 407.474_fp, 424.470_fp, &
    441.882_fp, 459.712_fp, 477.961_fp, 496.630_fp, 515.720_fp, 535.232_fp, 555.167_fp, 575.525_fp, &
    596.306_fp, 617.511_fp, 639.140_fp, 661.192_fp, 683.667_fp, 706.565_fp, 729.886_fp, 753.627_fp, &
    777.790_fp, 802.371_fp, 827.371_fp, 852.788_fp, 878.620_fp, 904.866_fp, 931.524_fp, 958.591_fp, &
    986.067_fp,1013.948_fp,1042.232_fp,1070.917_fp,1100.000_fp/)

    Pressure = &
    (/0.838_fp,   1.129_fp,   1.484_fp,   1.910_fp,   2.416_fp,   3.009_fp,   3.696_fp,   4.485_fp, &
      5.385_fp,   6.402_fp,   7.545_fp,   8.822_fp,  10.240_fp,  11.807_fp,  13.532_fp,  15.423_fp, &
     17.486_fp,  19.730_fp,  22.163_fp,  24.793_fp,  27.626_fp,  30.671_fp,  33.934_fp,  37.425_fp, &
     41.148_fp,  45.113_fp,  49.326_fp,  53.794_fp,  58.524_fp,  63.523_fp,  68.797_fp,  74.353_fp, &
     80.198_fp,  86.338_fp,  92.778_fp,  99.526_fp, 106.586_fp, 113.965_fp, 121.669_fp, 129.703_fp, &
    138.072_fp, 146.781_fp, 155.836_fp, 165.241_fp, 175.001_fp, 185.121_fp, 195.606_fp, 206.459_fp, &
    217.685_fp, 229.287_fp, 241.270_fp, 253.637_fp, 266.392_fp, 279.537_fp, 293.077_fp, 307.014_fp, &
    321.351_fp, 336.091_fp, 351.236_fp, 366.789_fp, 382.751_fp, 399.126_fp, 415.914_fp, 433.118_fp, &
    450.738_fp, 468.777_fp, 487.236_fp, 506.115_fp, 525.416_fp, 545.139_fp, 565.285_fp, 585.854_fp, &
    606.847_fp, 628.263_fp, 650.104_fp, 672.367_fp, 695.054_fp, 718.163_fp, 741.693_fp, 765.645_fp, &
    790.017_fp, 814.807_fp, 840.016_fp, 865.640_fp, 891.679_fp, 918.130_fp, 944.993_fp, 972.264_fp, &
    999.942_fp,1028.025_fp,1056.510_fp,1085.394_fp/)

    Temperature = &
    (/256.186_fp, 252.608_fp, 247.762_fp, 243.314_fp, 239.018_fp, 235.282_fp, 233.777_fp, 234.909_fp, &
      237.889_fp, 241.238_fp, 243.194_fp, 243.304_fp, 242.977_fp, 243.133_fp, 242.920_fp, 242.026_fp, &
      240.695_fp, 239.379_fp, 238.252_fp, 236.928_fp, 235.452_fp, 234.561_fp, 234.192_fp, 233.774_fp, &
      233.305_fp, 233.053_fp, 233.103_fp, 233.307_fp, 233.702_fp, 234.219_fp, 234.959_fp, 235.940_fp, &
      236.744_fp, 237.155_fp, 237.374_fp, 238.244_fp, 239.736_fp, 240.672_fp, 240.688_fp, 240.318_fp, &
      239.888_fp, 239.411_fp, 238.512_fp, 237.048_fp, 235.388_fp, 233.551_fp, 231.620_fp, 230.418_fp, &
      229.927_fp, 229.511_fp, 229.197_fp, 228.947_fp, 228.772_fp, 228.649_fp, 228.567_fp, 228.517_fp, &
      228.614_fp, 228.861_fp, 229.376_fp, 230.223_fp, 231.291_fp, 232.591_fp, 234.013_fp, 235.508_fp, &
      237.041_fp, 238.589_fp, 240.165_fp, 241.781_fp, 243.399_fp, 244.985_fp, 246.495_fp, 247.918_fp, &
      249.073_fp, 250.026_fp, 251.113_fp, 252.321_fp, 253.550_fp, 254.741_fp, 256.089_fp, 257.692_fp, &
      259.358_fp, 261.010_fp, 262.779_fp, 264.702_fp, 266.711_fp, 268.863_fp, 271.103_fp, 272.793_fp, &
      273.356_fp, 273.356_fp, 273.356_fp, 273.356_fp/)

    water_vapor = &
    (/4.187E-03_fp,4.401E-03_fp,4.250E-03_fp,3.688E-03_fp,3.516E-03_fp,3.739E-03_fp,3.694E-03_fp,3.449E-03_fp, &
      3.228E-03_fp,3.212E-03_fp,3.245E-03_fp,3.067E-03_fp,2.886E-03_fp,2.796E-03_fp,2.704E-03_fp,2.617E-03_fp, &
      2.568E-03_fp,2.536E-03_fp,2.506E-03_fp,2.468E-03_fp,2.427E-03_fp,2.438E-03_fp,2.493E-03_fp,2.543E-03_fp, &
      2.586E-03_fp,2.632E-03_fp,2.681E-03_fp,2.703E-03_fp,2.636E-03_fp,2.512E-03_fp,2.453E-03_fp,2.463E-03_fp, &
      2.480E-03_fp,2.499E-03_fp,2.526E-03_fp,2.881E-03_fp,3.547E-03_fp,4.023E-03_fp,4.188E-03_fp,4.223E-03_fp, &
      4.252E-03_fp,4.275E-03_fp,4.105E-03_fp,3.675E-03_fp,3.196E-03_fp,2.753E-03_fp,2.338E-03_fp,2.347E-03_fp, &
      2.768E-03_fp,3.299E-03_fp,3.988E-03_fp,4.531E-03_fp,4.625E-03_fp,4.488E-03_fp,4.493E-03_fp,4.614E-03_fp, &
      7.523E-03_fp,1.329E-02_fp,2.468E-02_fp,4.302E-02_fp,6.688E-02_fp,9.692E-02_fp,1.318E-01_fp,1.714E-01_fp, &
      2.149E-01_fp,2.622E-01_fp,3.145E-01_fp,3.726E-01_fp,4.351E-01_fp,5.002E-01_fp,5.719E-01_fp,6.507E-01_fp, &
      7.110E-01_fp,7.552E-01_fp,8.127E-01_fp,8.854E-01_fp,9.663E-01_fp,1.050E+00_fp,1.162E+00_fp,1.316E+00_fp, &
      1.494E+00_fp,1.690E+00_fp,1.931E+00_fp,2.226E+00_fp,2.574E+00_fp,2.939E+00_fp,3.187E+00_fp,3.331E+00_fp, &
      3.352E+00_fp,3.260E+00_fp,3.172E+00_fp,3.087E+00_fp/)
    ozone = &
    (/3.035E+00_fp,3.943E+00_fp,4.889E+00_fp,5.812E+00_fp,6.654E+00_fp,7.308E+00_fp,7.660E+00_fp,7.745E+00_fp, &
      7.696E+00_fp,7.573E+00_fp,7.413E+00_fp,7.246E+00_fp,7.097E+00_fp,6.959E+00_fp,6.797E+00_fp,6.593E+00_fp, &
      6.359E+00_fp,6.110E+00_fp,5.860E+00_fp,5.573E+00_fp,5.253E+00_fp,4.937E+00_fp,4.625E+00_fp,4.308E+00_fp, &
      3.986E+00_fp,3.642E+00_fp,3.261E+00_fp,2.874E+00_fp,2.486E+00_fp,2.102E+00_fp,1.755E+00_fp,1.450E+00_fp, &
      1.208E+00_fp,1.087E+00_fp,1.030E+00_fp,1.005E+00_fp,1.010E+00_fp,1.028E+00_fp,1.068E+00_fp,1.109E+00_fp, &
      1.108E+00_fp,1.071E+00_fp,9.928E-01_fp,8.595E-01_fp,7.155E-01_fp,5.778E-01_fp,4.452E-01_fp,3.372E-01_fp, &
      2.532E-01_fp,1.833E-01_fp,1.328E-01_fp,9.394E-02_fp,6.803E-02_fp,5.152E-02_fp,4.569E-02_fp,4.855E-02_fp, &
      5.461E-02_fp,6.398E-02_fp,7.205E-02_fp,7.839E-02_fp,8.256E-02_fp,8.401E-02_fp,8.412E-02_fp,8.353E-02_fp, &
      8.269E-02_fp,8.196E-02_fp,8.103E-02_fp,7.963E-02_fp,7.741E-02_fp,7.425E-02_fp,7.067E-02_fp,6.702E-02_fp, &
      6.368E-02_fp,6.070E-02_fp,5.778E-02_fp,5.481E-02_fp,5.181E-02_fp,4.920E-02_fp,4.700E-02_fp,4.478E-02_fp, &
      4.207E-02_fp,3.771E-02_fp,3.012E-02_fp,1.941E-02_fp,9.076E-03_fp,2.980E-03_fp,5.117E-03_fp,1.160E-02_fp, &
      1.428E-02_fp,1.428E-02_fp,1.428E-02_fp,1.428E-02_fp/)
    !AerosolType = DUST_AEROSOL
    aerosolType = 1
    aerosolEffectiveRadius = & ! microns
      (/0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 5.305110E-16_fp, &
        7.340409E-16_fp, 1.037097E-15_fp, 1.496791E-15_fp, 2.207471E-15_fp, 3.327732E-15_fp, &
        5.128933E-15_fp, 8.083748E-15_fp, 1.303055E-14_fp, 2.148368E-14_fp, 3.622890E-14_fp, &
        6.248544E-14_fp, 1.102117E-13_fp, 1.987557E-13_fp, 3.663884E-13_fp, 6.901587E-13_fp, &
        1.327896E-12_fp, 2.608405E-12_fp, 5.228012E-12_fp, 1.068482E-11_fp, 2.225098E-11_fp, &
        4.717675E-11_fp, 1.017447E-10_fp, 2.229819E-10_fp, 4.960579E-10_fp, 1.118899E-09_fp, &
        2.555617E-09_fp, 5.902789E-09_fp, 1.376717E-08_fp, 3.237321E-08_fp, 7.662427E-08_fp, &
        1.822344E-07_fp, 4.346896E-07_fp, 1.037940E-06_fp, 2.475858E-06_fp, 5.887266E-06_fp, &
        1.392410E-05_fp, 3.267943E-05_fp, 7.592447E-05_fp, 1.741777E-04_fp, 3.935216E-04_fp, &
        8.732308E-04_fp, 1.897808E-03_fp, 4.027868E-03_fp, 8.323272E-03_fp, 1.669418E-02_fp, &
        3.239702E-02_fp, 6.063055E-02_fp, 1.090596E-01_fp, 1.878990E-01_fp, 3.089856E-01_fp, &
        4.832092E-01_fp, 7.159947E-01_fp, 1.001436E+00_fp, 1.317052E+00_fp, 1.622354E+00_fp, &
        1.864304E+00_fp, 1.990457E+00_fp, 1.966354E+00_fp, 1.789883E+00_fp, 1.494849E+00_fp, &
        1.140542E+00_fp, 7.915451E-01_fp, 4.974823E-01_fp, 2.818937E-01_fp, 1.433668E-01_fp, &
        6.514795E-02_fp, 2.633057E-02_fp, 9.421763E-03_fp, 2.971053E-03_fp, 8.218245E-04_fp/)
 
    aerosolConcentration = & ! kg/m^2
        (/0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, &
        0.000000E+00_fp, 0.000000E+00_fp, 0.000000E+00_fp, 2.458105E-18_fp, 1.983430E-16_fp, &
        1.191432E-14_fp, 5.276880E-13_fp, 1.710270E-11_fp, 4.035105E-10_fp, 6.911389E-09_fp, &
        8.594215E-08_fp, 7.781797E-07_fp, 5.162773E-06_fp, 2.534018E-05_fp, 9.325154E-05_fp, &
        2.617738E-04_fp, 5.727150E-04_fp, 1.002153E-03_fp, 1.446048E-03_fp, 1.782757E-03_fp, &
        1.955759E-03_fp, 1.999206E-03_fp, 1.994698E-03_fp, 1.913109E-03_fp, 1.656122E-03_fp, &
        1.206328E-03_fp, 6.847261E-04_fp, 2.785695E-04_fp, 7.418821E-05_fp, 1.172680E-05_fp, &
        9.900895E-07_fp, 3.987399E-08_fp, 6.786932E-10_fp, 4.291151E-12_fp, 8.785440E-15_fp/)


        cloudType = 1
        cloudEffectiveRadius(:)=0.0_fp
        cloudConcentration(:)=0.0_fp

        cloudEffectiveRadius(75:79) = 20.0_fp ! microns
        cloudConcentration(75:79)    = 5.0_fp  ! kg/m^2
        co2(:) = 380.0_fp
 
  end subroutine test_data_us_std

  SUBROUTINE test_data_tropical(Level_Pressure, Pressure, Temperature, water_vapor, ozone)
    real, intent(out), dimension(93) :: Level_Pressure 
    real, intent(out), dimension(92) :: Pressure, Temperature, water_vapor, ozone
    ! ...Profile data
    Level_Pressure = &
    (/0.714,   0.975,   1.297,   1.687,   2.153,   2.701,   3.340,   4.077, &
      4.920,   5.878,   6.957,   8.165,   9.512,  11.004,  12.649,  14.456, &
     16.432,  18.585,  20.922,  23.453,  26.183,  29.121,  32.274,  35.650, &
     39.257,  43.100,  47.188,  51.528,  56.126,  60.990,  66.125,  71.540, &
     77.240,  83.231,  89.520,  96.114, 103.017, 110.237, 117.777, 125.646, &
    133.846, 142.385, 151.266, 160.496, 170.078, 180.018, 190.320, 200.989, &
    212.028, 223.441, 235.234, 247.409, 259.969, 272.919, 286.262, 300.000, &
    314.137, 328.675, 343.618, 358.967, 374.724, 390.893, 407.474, 424.470, &
    441.882, 459.712, 477.961, 496.630, 515.720, 535.232, 555.167, 575.525, &
    596.306, 617.511, 639.140, 661.192, 683.667, 706.565, 729.886, 753.627, &
    777.790, 802.371, 827.371, 852.788, 878.620, 904.866, 931.524, 958.591, &
    986.067,1013.948,1042.232,1070.917,1100.000/)

    Pressure = &
    (/0.838,   1.129,   1.484,   1.910,   2.416,   3.009,   3.696,   4.485, &
      5.385,   6.402,   7.545,   8.822,  10.240,  11.807,  13.532,  15.423, &
     17.486,  19.730,  22.163,  24.793,  27.626,  30.671,  33.934,  37.425, &
     41.148,  45.113,  49.326,  53.794,  58.524,  63.523,  68.797,  74.353, &
     80.198,  86.338,  92.778,  99.526, 106.586, 113.965, 121.669, 129.703, &
    138.072, 146.781, 155.836, 165.241, 175.001, 185.121, 195.606, 206.459, &
    217.685, 229.287, 241.270, 253.637, 266.392, 279.537, 293.077, 307.014, &
    321.351, 336.091, 351.236, 366.789, 382.751, 399.126, 415.914, 433.118, &
    450.738, 468.777, 487.236, 506.115, 525.416, 545.139, 565.285, 585.854, &
    606.847, 628.263, 650.104, 672.367, 695.054, 718.163, 741.693, 765.645, &
    790.017, 814.807, 840.016, 865.640, 891.679, 918.130, 944.993, 972.264, &
    999.942,1028.025,1056.510,1085.394/)

    Temperature = &
    (/266.536, 269.608, 270.203, 264.526, 251.578, 240.264, 235.095, 232.959, &
      233.017, 233.897, 234.385, 233.681, 232.436, 231.607, 231.192, 230.808, &
      230.088, 228.603, 226.407, 223.654, 220.525, 218.226, 216.668, 215.107, &
      213.538, 212.006, 210.507, 208.883, 206.793, 204.415, 202.058, 199.718, &
      197.668, 196.169, 194.993, 194.835, 195.648, 196.879, 198.830, 201.091, &
      203.558, 206.190, 208.900, 211.736, 214.601, 217.522, 220.457, 223.334, &
      226.156, 228.901, 231.557, 234.173, 236.788, 239.410, 242.140, 244.953, &
      247.793, 250.665, 253.216, 255.367, 257.018, 258.034, 258.778, 259.454, &
      260.225, 261.251, 262.672, 264.614, 266.854, 269.159, 271.448, 273.673, &
      275.955, 278.341, 280.822, 283.349, 285.826, 288.288, 290.721, 293.135, &
      295.609, 298.173, 300.787, 303.379, 305.960, 308.521, 310.916, 313.647, &
      315.244, 315.244, 315.244, 315.244/)

    water_vapor = &
    (/3.887E-03,3.593E-03,3.055E-03,2.856E-03,2.921E-03,2.555E-03,2.392E-03,2.605E-03, &
      2.573E-03,2.368E-03,2.354E-03,2.333E-03,2.312E-03,2.297E-03,2.287E-03,2.283E-03, &
      2.282E-03,2.286E-03,2.296E-03,2.309E-03,2.324E-03,2.333E-03,2.335E-03,2.335E-03, &
      2.333E-03,2.340E-03,2.361E-03,2.388E-03,2.421E-03,2.458E-03,2.492E-03,2.523E-03, &
      2.574E-03,2.670E-03,2.789E-03,2.944E-03,3.135E-03,3.329E-03,3.530E-03,3.759E-03, &
      4.165E-03,4.718E-03,5.352E-03,6.099E-03,6.845E-03,7.524E-03,8.154E-03,8.381E-03, &
      8.214E-03,8.570E-03,9.672E-03,1.246E-02,1.880E-02,2.720E-02,3.583E-02,4.462E-02, &
      4.548E-02,3.811E-02,3.697E-02,4.440E-02,2.130E-01,6.332E-01,9.945E-01,1.073E+00, &
      1.196E+00,1.674E+00,2.323E+00,2.950E+00,3.557E+00,4.148E+00,4.666E+00,5.092E+00, &
      5.487E+00,5.852E+00,6.137E+00,6.297E+00,6.338E+00,6.234E+00,5.906E+00,5.476E+00, &
      5.176E+00,4.994E+00,4.884E+00,4.832E+00,4.791E+00,4.760E+00,4.736E+00,6.368E+00, &
      7.897E+00,7.673E+00,7.458E+00,7.252E+00/)

    ozone = &
    (/2.742E+00,3.386E+00,4.164E+00,5.159E+00,6.357E+00,7.430E+00,8.174E+00,8.657E+00, &
      8.930E+00,9.056E+00,9.077E+00,8.988E+00,8.778E+00,8.480E+00,8.123E+00,7.694E+00, &
      7.207E+00,6.654E+00,6.060E+00,5.464E+00,4.874E+00,4.299E+00,3.739E+00,3.202E+00, &
      2.688E+00,2.191E+00,1.710E+00,1.261E+00,8.835E-01,5.551E-01,3.243E-01,1.975E-01, &
      1.071E-01,7.026E-02,6.153E-02,5.869E-02,6.146E-02,6.426E-02,6.714E-02,6.989E-02, &
      7.170E-02,7.272E-02,7.346E-02,7.383E-02,7.406E-02,7.418E-02,7.424E-02,7.411E-02, &
      7.379E-02,7.346E-02,7.312E-02,7.284E-02,7.274E-02,7.273E-02,7.272E-02,7.270E-02, &
      7.257E-02,7.233E-02,7.167E-02,7.047E-02,6.920E-02,6.803E-02,6.729E-02,6.729E-02, &
      6.753E-02,6.756E-02,6.717E-02,6.615E-02,6.510E-02,6.452E-02,6.440E-02,6.463E-02, &
      6.484E-02,6.487E-02,6.461E-02,6.417E-02,6.382E-02,6.378E-02,6.417E-02,6.482E-02, &
      6.559E-02,6.638E-02,6.722E-02,6.841E-02,6.944E-02,6.720E-02,6.046E-02,4.124E-02, &
      2.624E-02,2.623E-02,2.622E-02,2.622E-02/)

  end subroutine test_data_tropical
end module pycrtm
