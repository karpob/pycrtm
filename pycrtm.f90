module pycrtm 
real(kind=8), allocatable :: outTransmission(:, :, :) ! outTransmission(N_profiles, nChan, N_Layers)

real(kind=8), allocatable :: aerosolEffectiveRadius(:,:,:) !(N_Profiles,N_layers, N_aerosols)
real(kind=8), allocatable :: aerosolConcentration(:,:,:)   !(N_profiles,N_layers, N_aerosols)
integer, allocatable :: aerosolType(:,:)                   !(N_Profiles, N_aerosols)


real(kind=8), allocatable :: cloudEffectiveRadius(:,:,:) !(N_Profiles,N_layers, N_clouds)
real(kind=8), allocatable :: cloudConcentration(:,:,:)   !(N_profiles,N_layers, N_clouds)
real(kind=8), allocatable :: cloudFraction(:,:)          !(N_profiles,N_layers)
integer, allocatable :: cloudType(:,:)                   !(N_Profiles, N_clouds)

contains
subroutine wrap_forward( coefficientPath, sensor_id_in, IRwaterCoeff_File, MWwaterCoeff_File, & 
                        output_tb_flag, output_transmission_flag,  zenithAngle, scanAngle, azimuthAngle, solarAngle, &
                        year, month, day, & 
                        nChan, N_Profiles, N_LAYERS, N_trace, &
                        pressureLevels, pressureLayers, temperatureLayers, & 
                        traceConcLayers, trace_IDs, & 
                        climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m, & 
                        landType, soilType, vegType, waterType, snowType, iceType, nthreads, &  
                        outTb, & 
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
  character(len=*), intent(in) :: IRwaterCoeff_File
  character(len=*), intent(in) :: MWwaterCoeff_File
  logical, intent(in) :: output_tb_flag, output_transmission_flag 
  ! The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  integer, intent(in) :: nChan, N_Profiles, N_Layers, N_trace
  real(kind=8), intent(in) :: zenithAngle(n_profiles), scanAngle(n_profiles) 
  real(kind=8), intent(in) :: azimuthAngle(n_profiles), solarAngle(n_profiles,2)
  integer, intent(in) :: year(n_profiles), month(n_profiles), day(n_profiles) 
  real(kind=8), intent(in) :: pressureLevels(N_profiles, N_LAYERS+1)
  real(kind=8), intent(in) :: pressureLayers(N_profiles, N_LAYERS), temperatureLayers(N_Profiles,N_Layers)
  real(kind=8), intent(in) :: traceConcLayers(N_Profiles,N_layers,N_trace)
  integer, intent(in) :: trace_IDs(N_trace)
  integer, intent(in) :: climatology(N_profiles)
  real(kind=8), intent(in) :: surfaceTemperatures(N_Profiles,4), surfaceFractions(N_profiles, 4)
  real(kind=8), intent(in) :: LAI(N_Profiles), salinity(N_Profiles),  windSpeed10m(N_Profiles), windDirection10m(N_Profiles)
  integer, intent(in) ::  landType(N_Profiles), soilType(N_Profiles), vegType(N_Profiles), waterType(N_Profiles)
  integer, intent(in) ::  snowType(N_Profiles), iceType(N_Profiles) 
  integer, intent(in) :: nthreads
  real(kind=8), intent(out) :: outTb(N_Profiles,nChan) 
  real(kind=8), intent(inout) :: emissivityReflectivity(2,N_Profiles,nChan)
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
  INTEGER :: i, l, m, n, nc, ll,mm, nn, species,i_abs
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
  type(crtm_options_type)                 :: options(1)
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
  ! allocate globals in the module based upon user selection through interface.
  call check_and_allocate_globals(output_transmission_flag, N_Profiles, nChan, N_layers)
  ! Figure out what needs allocating for Clouds and Aerosols. Are they on?
  call aerosols_and_clouds_on(N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)
  
  err_stat = CRTM_Init( sensor_id,  chinfo, &
                        File_Path=coefficientPath, &
                        Load_CloudCoeff = cloudsOn, &  
                        Load_AerosolCoeff = aerosolsOn, &
                        IRwaterCoeff_File = IRwaterCoeff_File, & 
                        MWwaterCoeff_File = MWwaterCoeff_File, & 
                        Quiet=.True. )
  call check_allocate_status(err_stat,'Error Initializing CRTM')

  WRITE( *,'(/5x,"Initializing the CRTM...")' )

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
  !$omp& shared(zenithAngle, scanAngle, azimuthAngle, solarAngle, output_tb_flag)&
  !$omp& shared(nChan, N_Profiles, N_LAYERS, N_Clouds_crtm, N_aerosols_crtm, N_trace)&
  !$omp& shared(pressureLevels, pressureLayers, temperatureLayers)& 
  !$omp& shared(traceConcLayers,trace_IDs, output_transmission_flag)& 
  !$omp& shared(aerosolEffectiveRadius, aerosolConcentration, aerosolType, cloudsOn, aerosolsOn)& 
  !$omp& shared(cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology)& 
  !$omp& shared(surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m)& 
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
    call check_allocate_status(alloc_stat, "Error allocating Solution rts(n_channels,1).")
    ! allocate 1 profile
    ALLOCATE( atm(1), STAT = alloc_stat )
    call check_allocate_status(alloc_stat, "Error allocating atm.")
    ! allocate 1 surface
    ALLOCATE( sfc(1), STAT = alloc_stat )
    call check_allocate_status(alloc_stat, "Error allocating sfc.")

    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_trace, N_CLOUDS_crtm, N_AEROSOLS_crtm )
    call check_logical_status(ANY(.not. CRTM_Atmosphere_Associated(atm) ), "Failed in CRTM_Atmopsphere_Create")

    ! ==========================================================================
    ! STEP 6. **** ASSIGN INPUT DATA ****
    !
    ! 6a. Atmosphere and Surface input
    !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
    !     by which the atmosphere and surface data are loaded in to their
    !     respective structures below was done purely to keep the step-by-step
    !     instructions in this program relatively "clean".
    ! ------------------------------------------------------------------------
    ! ...Profile data
    call set_profile(atm, n, climatology(n), pressureLevels(n,:), pressureLayers(n,:), temperatureLayers(n,:),&
                         traceConcLayers(n,:,:), trace_IDs(:), &
                         N_trace, N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)

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
    ! set the surface properties for the profile.
    call set_surface(sfc, surfaceFractions(n,:), landType(n), surfaceTemperatures(n,:), LAI(n), & 
                           soilType(n), vegType(n), waterType(n), snowType(n), iceType(n), &
                           windSpeed10m(n), windDirection10m(n), salinity(n))
    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    ! 8a. The forward model
    ! ---------------------

    ! Need this to get transmission out of solution, otherwise won't be allocated !!!
    call crtm_rtsolution_create( rts, n_layers )
    call check_logical_status( any(.not. crtm_rtsolution_associated( rts ) ),'rts failed to create.') 

    call crtm_options_create( options, nChan )
    call check_logical_status( any(.not. crtm_options_associated( options ) ),'options failed to create.' )
    call set_emissivity(options, emissivityReflectivity(1,n,:), emissivityReflectivity(2,n,:))

    err_stat = CRTM_Forward( atm        , &  ! Input
                             sfc        , &  ! Input
                             geo        , &  ! Input
                             chinfo     , &  ! Input
                             rts        , &  ! Output
                             options = options ) 

    call check_allocate_status(err_stat, "Error Calling CRTM_Forward.")

    ! ============================================================================
    ! 8c. **** OUTPUT THE RESULTS TO SCREEN **** (Or transfer it into a series of arrays out of this thing!)
    !
    ! User should read the user guide or the source code of the routine
    ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
    ! select the needed variables for outputs.  These variables are contained
    ! in the structure RTSolution.
    if (output_transmission_flag) then 
        do l=1,nChan
            outTransmission(n, l,1:n_layers) = & 
             dexp(-1.0*cumsum( rts(l,1)%Layer_Optical_Depth ) )
        enddo
    endif
    emissivityReflectivity(1,n,:) = rts(:,1)%Surface_Emissivity 
    emissivityReflectivity(2,n,:) = rts(:,1)%Surface_Reflectivity

    if (output_tb_flag) then 
        outTb(n,:) = rts(:,1)%Brightness_Temperature
    else
        outTb(n,:) = rts(:,1)%Radiance
    endif 
 
    
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
    call crtm_options_destroy(options)
    deallocate(atm,stat=alloc_stat)
    call check_allocate_status(alloc_stat,"Atm failed to deallocate.")
    deallocate(sfc, stat=alloc_stat)
    call check_allocate_status(alloc_stat,"Sfc failed to deallocate.")
    DEALLOCATE(rts, STAT = alloc_stat)
    call check_allocate_status(alloc_stat,"Rts failed to deallocate.")
  END DO Profile_Loop
  !$omp end parallel do
  
  ! ==========================================================================
  ! 10. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  call check_allocate_status(err_stat, 'Error Destroying the CRTM')
  write(*,*)'wrap_forward done!'

end subroutine wrap_forward

subroutine wrap_k_matrix( coefficientPath, sensor_id_in, IRwaterCoeff_File, MWwaterCoeff_File, & 
                        output_tb_flag, output_transmission_flag, zenithAngle, scanAngle, azimuthAngle, solarAngle, &  
                        year, month, day, & 
                        nChan, N_profiles, N_LAYERS, N_trace, & 
                        pressureLevels, pressureLayers, temperatureLayers, & 
                        traceConcLayers, trace_IDs, & 
                        climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, salinity, windSpeed10m, windDirection10m, & 
                        landType, soilType, vegType, waterType, snowType, iceType, &  
                        nthreads, outTb, & 
                        temperatureJacobian, traceJacobian, skinK, emisK, reflK, &
                        windSpeedK, windDirectionK,  emissivityReflectivity )      

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
  character(len=*), intent(in) :: coefficientPath
  character(len=*), intent(in) :: sensor_id_in
  character(len=*), intent(in) :: IRwaterCoeff_File
  character(len=*), intent(in) :: MWwaterCoeff_File
  logical, intent(in) :: output_tb_flag, output_transmission_flag
  integer, intent(in) :: nChan, N_profiles, N_Layers, N_trace 
  ! The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  real(kind=8), intent(in) :: zenithAngle(N_profiles), scanAngle(N_profiles)
  real(kind=8), intent(in) :: azimuthAngle(N_profiles), solarAngle(N_profiles,2)
  integer, intent(in) :: year(n_profiles), month(n_profiles), day(n_profiles)
  real(kind=8), intent(in) :: pressureLevels(N_profiles, N_Layers+1)
  real(kind=8), intent(in) :: pressureLayers(N_profiles, N_layers), temperatureLayers(N_profiles, N_layers)
  real(kind=8), intent(in) :: traceConcLayers(N_profiles, N_layers, N_trace)
  integer, intent(in) :: trace_IDs(N_trace)
  integer, intent(in) ::  climatology(N_profiles)
  real(kind=8), intent(in) :: surfaceTemperatures(N_profiles,4), surfaceFractions(N_profiles,4), LAI(N_profiles) 
  real(kind=8), intent(in) :: salinity(N_profiles), windSpeed10m(N_profiles), windDirection10m(N_profiles)
  integer, intent(in) :: landType(N_profiles), soilType(N_profiles), vegType(N_profiles), waterType(N_profiles) 
  integer, intent(in) :: snowType(N_profiles), iceType(N_profiles) 
  real(kind=8), intent(out) :: outTb(N_profiles,nChan)
  real(kind=8), intent(inout):: emissivityReflectivity(2,N_profiles,nChan)
  real(kind=8), intent(out) :: skinK(N_profiles,nChan,4), emisK(N_profiles,nChan), reflK(N_profiles,nChan)
  real(kind=8), intent(out) :: windSpeedK(N_profiles,nChan), windDirectionK(N_profiles,nChan)
  real(kind=8), intent(out) :: temperatureJacobian(N_profiles, nChan, N_LAYERS)
  real(kind=8), intent(out) :: traceJacobian(N_profiles, nChan, N_LAYERS, N_trace)
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
  INTEGER :: l, m, n, nc, species, i_abs
  Logical :: cloudsOn, aerosolsOn


  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(1)
  TYPE(CRTM_Geometry_type)                :: geo(1)
  type(crtm_options_type)                 :: options(1)
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
  ! allocate globals in the module based upon user selection through interface.
  call check_and_allocate_globals(output_transmission_flag, N_Profiles, nChan, N_layers)

  ! figure out how to allocate aerosols/clouds and are the even turned on by the user?
  call aerosols_and_clouds_on( N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)

  WRITE( *,'(/5x,"Initializing the CRTM...")' )

  err_stat = CRTM_Init( sensor_id,  chinfo, &
                        File_Path=coefficientPath, &
                        Load_CloudCoeff = cloudsOn, &  
                        Load_AerosolCoeff = aerosolsOn, &
                        IRwaterCoeff_File = IRwaterCoeff_File, & 
                        MWwaterCoeff_File = MWwaterCoeff_File, &   
                        Quiet=.True. )
  call check_allocate_status(err_stat, 'Error initializing CRTM')

  ! 4b. Output some channel information
  ! -----------------------------------
  n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels
  WRITE( *,'(7x,i0," from ",a)' )  CRTM_ChannelInfo_n_Channels(chinfo(1)), TRIM(sensor_id(1))

  ! Begin loop over sensors
  ! ----------------------
  !$ call omp_set_num_threads(nthreads)
  !$omp parallel do default(private) shared(emissivityReflectivity,outTb)&
  !$omp& shared(temperatureJacobian, output_tb_flag, output_transmission_flag)& 
  !$omp& shared(traceJacobian,outTransmission)&
  !$omp& shared(nChan, N_Layers,N_CLOUDS_crtm, N_AEROSOLS_crtm, N_trace)&
  !$omp& shared(pressureLevels, pressureLayers, temperatureLayers)& 
  !$omp& shared(traceConcLayers, trace_IDs, cloudsOn, aerosolsOn, zenithAngle,scanAngle,azimuthAngle,solarAngle)& 
  !$omp& shared(aerosolEffectiveRadius, aerosolConcentration, aerosolType)& 
  !$omp& shared(cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology)& 
  !$omp& shared(surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m)& 
  !$omp& shared(skinK, emisK, reflK, windSpeedK, windDirectionK)& 
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
    allocate(atm(1), STAT = alloc_stat)
    call check_allocate_status(alloc_stat,'Error allocating atm(1)')

    allocate(sfc(1), STAT = alloc_stat)
    call check_allocate_status(alloc_stat,'Error allocating sfc(1)')

    ALLOCATE( rts( n_channels, 1 ), STAT = alloc_stat )
    call check_allocate_status(alloc_stat,'Error allocating rts(n_channels,1)')

    ALLOCATE( atm_K( n_channels, 1 ), STAT = alloc_stat )
    call check_allocate_status(alloc_stat,'Error allocating atm_k')

    ALLOCATE( sfc_K( n_channels, 1 ), STAT = alloc_stat )
    call check_allocate_status(alloc_stat,'Error allocating sfc_k')

    ALLOCATE( rts_K( n_channels, 1 ), STAT = alloc_stat )
    call check_allocate_status(alloc_stat,'Error allocating rts_k')

    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_trace, N_CLOUDS_crtm, N_AEROSOLS_crtm )
    call check_logical_status(ANY(.NOT. CRTM_Atmosphere_Associated(atm)), 'Error in CRTM_Atmosphere_Create Atm()')

    ! The output K-MATRIX structure
    CALL CRTM_Atmosphere_Create( atm_K, N_LAYERS, N_trace, N_CLOUDS_crtm, N_AEROSOLS_crtm )
    call check_logical_status(ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)),  'Error in CRTM_Atmosphere_Create Atm_k()')

    call crtm_rtsolution_create( rts, n_layers )
    call check_logical_status(any(.not. crtm_rtsolution_associated( rts )),  'Error in crtm_rtsolution_create rts()')
    
    call crtm_rtsolution_create( rts_k, n_layers )
    call check_logical_status(any(.not. crtm_rtsolution_associated( rts_k )),  'Error in crtm_rtsolution_create rts_k()')

    ! ==========================================================================
    ! STEP 6. **** ASSIGN INPUT DATA ****
    !
    ! 6a. Atmosphere and Surface input
    !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
    !     by which the atmosphere and surface data are loaded in to their
    !     respective structures below was done purely to keep the step-by-step
    !     instructions in this program relatively "clean".
    ! ------------------------------------------------------------------------
    ! ...Profile data
    call set_profile(atm, n, climatology(n), pressureLevels(n,:), pressureLayers(n,:), temperatureLayers(n,:),&
                         traceConcLayers(n,:,:), trace_IDs(:), &
                         N_trace, N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)
 

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
    if (output_tb_flag) then
        rts_K%Radiance               = ZERO
        rts_K%Brightness_Temperature = ONE
    else 
        rts_K%Radiance               = ONE
        rts_K%Brightness_Temperature = ZERO
    endif
    ! ==========================================================================
    
    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
     call set_surface(sfc, surfaceFractions(n,:), landType(n), surfaceTemperatures(n,:), LAI(n), & 
                           soilType(n), vegType(n), waterType(n), snowType(n), iceType(n), &
                           windSpeed10m(n), windDirection10m(n), salinity(n))
    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    ! 8b. The K-matrix model
    ! ----------------------
    call crtm_options_create( options, nChan )
    call check_logical_status( any(.not. crtm_options_associated( options ) ),'options failed to create' )
    call set_emissivity(options, emissivityReflectivity(1,n,:), emissivityReflectivity(2,n,:))

    err_stat = CRTM_K_Matrix( atm        , &  ! FORWARD  Input
                              sfc        , &  ! FORWARD  Input
                              rts_K      , &  ! K-MATRIX Input
                              geo        , &  ! Input
                              chinfo     , &  ! Input
                              atm_K      , &  ! K-MATRIX Output
                              sfc_K      , &  ! K-MATRIX Output
                              rts        , &  ! FORWARD  Output
                              Options=options)
    call check_allocate_status(err_stat,'Error calling the CRTM K-Matrix Model')    

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
        !jacobians of H2O, O3, etc... will be determined by the order in which they were assigned in atm. 
        do i_abs=1,N_trace
            traceJacobian(n,l, 1:n_layers,i_abs) = atm_k(l,1)%Absorber(1:n_layers,i_abs)
        enddo
        skinK(n,l,1) = sfc_K(l,1)%Land_Temperature
        skinK(n,l,2) = sfc_K(l,1)%Water_Temperature
        skinK(n,l,3) = sfc_K(l,1)%Ice_Temperature
        skinK(n,l,4) = sfc_K(l,1)%Snow_Temperature
        windSpeedK(n,l) = sfc_K(l,1)%Wind_Speed
        windDirectionK(n,l) = sfc_K(l,1)%Wind_Direction
        emisK(n,l) = RTS_K(l,1)%Surface_Emissivity
        reflK(n,l) = RTS_K(l,1)%Surface_Reflectivity
        if(output_transmission_flag) then 
             outTransmission(n, l,1:n_layers) = & 
             dexp(-1.0* cumsum( rts(l,1)%Layer_Optical_Depth ) ) 
        endif
    enddo
    if (output_tb_flag) then 
        outTb(n,:) = rts(:,1)%Brightness_Temperature
    else
        outTb(n,:) = rts(:,1)%Radiance
    endif 
    emissivityReflectivity(1,n,:) = rts(:,1)%Surface_Emissivity
    emissivityReflectivity(2,n,:) = rts(:,1)%Surface_Reflectivity
    CALL CRTM_Atmosphere_Destroy(atm)
    CALL CRTM_Atmosphere_Destroy(atm_k)
    call crtm_options_destroy(options)
    DEALLOCATE(atm_k, STAT = alloc_stat)
    call check_allocate_status(alloc_stat, 'Atm_k deallocate failed')

    DEALLOCATE(rts_K, STAT = alloc_stat)
    call check_allocate_status(alloc_stat, 'rts_k deallocate failed')

    DEALLOCATE(sfc_k, STAT = alloc_stat)
    call check_allocate_status(alloc_stat, 'sfc_k deallocate failed')

    DEALLOCATE(rts, STAT = alloc_stat)
    call check_allocate_status(alloc_stat, 'rts deallocate failed')

    DEALLOCATE(atm, STAT = alloc_stat)
    call check_allocate_status(alloc_stat, 'atm deallocate failed')

    DEALLOCATE(sfc, STAT = alloc_stat)
    call check_allocate_status(alloc_stat, 'sfc deallocate failed')
    ! ==========================================================================

  END DO Profile_Loop
  !$omp end parallel do
  ! ==========================================================================
  ! 10. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  call check_allocate_status(err_stat, 'Error destroying the CRTM.')
  ! ==========================================================================
end subroutine wrap_k_matrix

  subroutine check_and_allocate_globals(output_transmission_flag, N_Profiles, nChan, N_layers)
  logical, intent(in) :: output_transmission_flag
  integer, intent(in) :: N_profiles, nChan, N_layers
  if(output_transmission_flag) then
    if ( allocated(outTransmission) ) deallocate(outTransmission)
    allocate( outTransmission(N_Profiles, nChan, N_layers) ) 
  endif
  end subroutine check_and_allocate_globals

  subroutine applyAvg( Ref_LnPressure, User_LnPressure, Nref, Nuser, Xin, Xout ) 
    USE ODPS_CoordinateMapping
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

  subroutine aerosols_and_clouds_on(N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)

  integer, intent(out) :: N_AEROSOLS_crtm, N_CLOUDS_crtm
  logical, intent(out) :: aerosolsOn, cloudsOn
  integer :: shp(2)
  if( .not. allocated(aerosolType) ) then
    N_AEROSOLS_crtm = 0
    aerosolsOn = .False.
  else
    shp = shape(aerosolType)
    N_AEROSOLS_crtm = shp(2)
    aerosolsOn = .True. 
  endif

  if( .not. allocated(cloudType) ) then
    N_CLOUDS_crtm = 0
    cloudsOn = .False.
  else
    shp = shape(cloudType)
    N_CLOUDS_crtm = shp(2)  
    cloudsOn = .True. 
  endif
 
  end subroutine aerosols_and_clouds_on

  subroutine set_emissivity(options, emiss, refl)
    use crtm_module
    implicit none
    type(crtm_options_type), intent(inout) :: options(1)
    real(kind=8) :: emiss(:), refl(:)

    if ( all( emiss < -0.9999 ) ) then
        Options(1)%Use_Emissivity = .false.   ! compute it
    else
        Options(1)%Use_Emissivity = .true.    ! user supplied
        Options(1)%Emissivity(:)=emiss(:)
    endif 

    if ( all( refl < -0.9999) ) then
        Options(1)%Use_Direct_Reflectivity = .false.
    else
        Options(1)%Use_Direct_Reflectivity = .true.  ! 1: User-supplied
        Options(1)%Direct_Reflectivity(:)=refl(:) 
    endif
  end subroutine set_emissivity

  subroutine check_allocate_status(alloc_stat,message)
    integer :: alloc_stat
    character(len=*) :: message
    IF ( alloc_stat /= 0 ) THEN
      write(*,*) message
      STOP
    END IF
  end subroutine check_allocate_status

  subroutine check_logical_status(stat,message)
    logical :: stat
    character(len=*) :: message
    IF ( stat ) THEN
      write(*,*) message
      STOP
    END IF
  end subroutine check_logical_status

  subroutine set_profile(atm, n, climatology, pressureLevels, pressureLayers, temperatureLayers,&
                         traceConcLayers, trace_IDs, & 
                         N_trace, N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)
  USE CRTM_module
  TYPE(CRTM_Atmosphere_type), intent(inout) :: atm(:)
  integer :: n,climatology
  real(kind=8) :: pressureLevels(:), pressureLayers(:), temperatureLayers(:), traceConcLayers(:,:)
  integer :: trace_IDs(:) 
  integer :: N_trace, N_aerosols_crtm, N_clouds_crtm
  logical :: aerosolsOn, cloudsOn
  integer :: i_abs,species  
 
    atm(1)%Climatology = climatology
    atm(1)%Level_Pressure = pressureLevels(:)
    atm(1)%Pressure = pressureLayers(:)
    atm(1)%Temperature = temperatureLayers(:)
   
    do i_abs = 1,N_trace 
      atm(1)%Absorber(:,i_abs)      = traceConcLayers(:,i_abs)
      atm(1)%Absorber_Id(i_abs)     = trace_IDs(i_abs)
      if( trace_IDs(i_abs) == H2O_ID ) then 
        atm(1)%absorber_units(i_abs) = MASS_MIXING_RATIO_UNITS
      else 
        atm(1)%absorber_units(i_abs)  = VOLUME_MIXING_RATIO_UNITS
      endif
    enddo

    if( aerosolsOn )  then
      do species = 1, N_aerosols_crtm
        atm(1)%Aerosol(species)%Type                = aerosolType(n, species)
        atm(1)%Aerosol(species)%Effective_Radius(:) = aerosolEffectiveRadius(n, :, species)
        atm(1)%Aerosol(species)%Concentration(:)    = aerosolConcentration(n, :, species)
      enddo
    endif
    if( cloudsOn ) then
      do species = 1, N_clouds_crtm
        atm(1)%Cloud(species)%Type                = cloudType(n, species)
        atm(1)%Cloud(species)%Effective_Radius(:) = cloudEffectiveRadius(n, :, species)
        atm(1)%Cloud(species)%Water_Content(:)    = cloudConcentration(n, :, species)
      enddo
      atm(1)%Cloud_Fraction(:)            = cloudFraction(n,:)
    endif

  end subroutine set_profile
  subroutine set_surface(sfc, surfaceFractions, landType, surfaceTemperatures, LAI, soilType, & 
                         vegType, waterType, snowType, iceType, windSpeed10m, windDirection10m, & 
                         salinity)
    USE CRTM_module  
    TYPE(CRTM_Surface_type) :: sfc(:)
    real(kind=8) :: surfaceFractions(:), surfaceTemperatures(:), LAI, windSpeed10m, windDirection10m, salinity
    integer :: landType, soilType, vegType, waterType, snowType, iceType

    sfc%Land_Coverage     = surfaceFractions(1)
    sfc%Land_Type         = landType 
    sfc%Land_Temperature  = surfaceTemperatures(1)
    sfc%Lai               = LAI
    sfc%Soil_Type         = soilType 
    sfc%Vegetation_Type   = vegType 
    ! ...Water surface characteristics
    sfc%Water_Coverage    = surfaceFractions(2)
    sfc%Water_Type        = waterType 
    sfc%Water_Temperature = surfaceTemperatures(2)

    ! ...Snow coverage characteristics
    sfc%Snow_Coverage    = surfaceFractions(3)
    sfc%Snow_Type        = snowType 
    sfc%Snow_Temperature = surfaceTemperatures(3)
    ! ...Ice surface characteristics
    sfc%Ice_Coverage    = surfaceFractions(4)
    sfc%Ice_Type        = iceType 
    sfc%Ice_Temperature = surfaceTemperatures(4)

    sfc%Wind_Speed = windSpeed10m
    sfc%Wind_Direction = windDirection10m
    sfc%Salinity = salinity
  end subroutine set_surface
  function cumsum(x) result(xout)
    real(kind=8), DIMENSION(:), INTENT(IN) :: x
    real(kind=8), DIMENSION(size(x)) :: xout
    integer :: n,j
    n = size(x)
    do j=2,n
        xout(j) = xout(j-1) + x(j)
    end do
  end function cumsum

end module pycrtm
