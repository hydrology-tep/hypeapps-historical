!> \file assimilation_variables.f90

!> Contain module assimilation_variables, with variables needed to run the (EnKF) data assimilation routines within HYSS
!> Author: D.Gustafsson (SMHI)
!> Versions: 
!> 2012.05.10 Original adapted to HYPE ver 3.6
!> 2013.09.19 Adapted to HYPE ver 4.5.x
!> 2013.12.03 Renamed to assim_data.f90, to open up for more DA methods in addition to EnKF
!> 2014.11.09 Reformated these comments to be read  by doxygen
!> 2014.12.10 New data for saving ensemble data to direct access binary files instead of keeping in memory.
!> 2015.06.22 Clean-up version (not finalized)
!> 2015.09.18 Quick and dirty clean-up for the HYPE course 2015-09-25
!> 2016.10.20 File and module name changed to assimilation_variables(.f90), plus adaption to HYPE coding style (FORTRAN names in UPPER CASE and variable/routine names in lower case.
  
MODULE ASSIMILATION_VARIABLES
!Copyright 2016 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------------------

USE random_routines, ONLY: randxy_data

IMPLICIT NONE
! VARIABLES declared in this file
!
! name          content
! ----------    ---------------------
! EnkfInfo      General EnKF Settings
!
! EnkfXinfo     Info about model state variables    (prognostics)
! EnkfAinfo     -"-  -"-   -"-   auxiliary variables (diagnostics; HYPE outvars')
! EnkfPinfo     -"-  -"-   -"-   parameters
! EnkfFinfo     -"-  -"-   -"-   forcing variables 
! EnkfDinfo     -"-  -"-   -"-   observation variables 
!
! EnkfX         model state ensemble
! EnkfA         ensembles of auxiliary model variables used for output (outvars')
! EnkfP         parameter ensemble
! EnkfF         forcing ensemble
! EnkfD         observation ensemble  (size will change from time to time, depending on missing data)
! EnkfHX        model "Observations"  (size will change from time to time, depending on missing data)
! EnkfR         obs.error covariance (diagonal) matrix  (size will change from time to time, depending on missing data)
!---------------------------------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------------------------------
! ASSIMILATION APPLICATION INFO VARIABLE (application specific settings)
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_flag_type
  LOGICAL           :: add           !<assimilation status
  CHARACTER(LEN=50) :: category = '' !<name of category or statetype
  CHARACTER(LEN=50) :: varname  = '' !<name of variable, empty string means apply to whole category
ENDTYPE assim_flag_type

TYPE assim_info_type
  
  !Model structure variables
  INTEGER nX           !number of model state variables       (dimension of assimData%X vector)
  INTEGER nF           !number of model forcing variables     (dimension of assimData%F vector)
  INTEGER nA           !number of model outvar variables      (dimension of assimData%A vector)
  INTEGER nA2          !number of model outvarbasin variables (dimension of assimData%A2 vector)
  INTEGER nObs         !number of observation variables       (dimension of assimData%Obs vector)
  INTEGER nD           !number of observations (rows) in D and HX ensemble matrices for the current Analysis
  !INTEGER nP           !number of parameters in P ensemble  (rows in the P ensemble)
  !INTEGER nDT          !number of timesteps in assimilation time window
  
  INTEGER nLoc         !number of localization matrices     (depend on the number of spatial domains represented by model state variables, 
                       !                                     in HYPE there are several potential domains: subbasins, aquifers, parameter regions,
                       !                                     as well as non-spatial states such as general and land use parameters)
  INTEGER ncoord       !number of spatial domains           (length of EnkfCoordinates(:) vector, as above)
  
  !ENKF settings
  INTEGER nE           !number of ensemble members
  INTEGER locFunc      !localization function (0 none, 1 exponential (xy, z)
  
  !LOGICAL flags, may be modified...
  LOGICAL FA           !include auxiliaries in kalman filter (general switch on/off)
  LOGICAL FP           !include parameters  in kalman filter (general switch on/off) 
  LOGICAL FF           !include forcing     in kalman filter (general switch on/off)
  LOGICAL meanout      !print ensemble mean (.true.) or median (.false.) in output files
  
  !assimilation flags for different variable categories
  !the categories are defined (hardcoded) in initialize_assim_categories_MODEL in the interface
  TYPE(assim_flag_type),ALLOCATABLE :: assim_flag(:)
  INTEGER :: nFlag  !number of data in assim_flag structure
  INTEGER :: nCat   !size of assim_categories (number of possible state variables)
  CHARACTER(LEN=20), ALLOCATABLE :: assim_categories(:)
  !assimilation flags for different (HYPE model specific) variable categories, CP kan tas bort nu
  LOGICAL a_snow
  LOGICAL a_glacier
  LOGICAL a_lakeice
  LOGICAL a_riverice
  LOGICAL a_soil
  LOGICAL a_riverwt
  LOGICAL a_lakewt
  LOGICAL a_misc
  LOGICAL a_aquifer
    
  !some parameters
  REAL :: moradkhani_delta ! coef.[0-1] to retain variance in parameter ensemble (Moradkhani et al, 2004), value typical around 0.95
  
  !random number generation seed
  INTEGER :: seed
  
  !missing value
  REAL :: missing
  
  !minimum allowed sigma in random number generation (general)
  !REAL :: ensgen_minsigma
  
  !option to read/write ensemble data to direct access binary files instead of saving in allocated memory
  !LOGICAL useBinFiles
  INTEGER useBinFilesX   !0=save in memory, 1=save to one bin-fil, 2=save to several bin-files
  INTEGER useBinFilesFA  !0=save in memory, 1=save to one bin-fil, 2=save to several bin-files
  INTEGER nBinFiles
  
  !options regarding outputs
  INTEGER,ALLOCATABLE :: outvarflag(:), parameterflag(:)  !CP161207 These are not used today
  INTEGER :: nstatout   !Number of statistical extra outputs of simulation (e.g. 2 means min and max)
  
  !LOcalization parameters
  REAL xy_scalefac,z_scalefac

END TYPE assim_info_type

!---------------------------------------------------------------------------------------------------------------------------
! ENSEMBLE DATA MATRIX TYPE 
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_ensemble_type
  !number of rows, columns, and the ensemble matrix (y)
  INTEGER           :: nvar           !number of variables, ie n:o "model units" (for instance, number of sub-basins)
  INTEGER           :: nens           !number of ensemble members
  !ensemble data matrix (x)
  REAL, ALLOCATABLE :: x(:,:)         !ensemble data [nvar x nens]
  !variables used for min/max checks, etc
  REAL              :: minimum        !min allowed 
  REAL              :: maximum        !max allowed 
  !variables for output statistics
  REAL, ALLOCATABLE :: outmean(:), outquant(:,:), outmin(:), outmax(:), outsigma(:)
  !indices used for I/O
  INTEGER           :: fileID         !file id used for read/write to direct access binary file (if needed)
  INTEGER           :: rec            !record of binary file, CP161201
  !indices to link to spatial information
  INTEGER           :: coordID     !spatial coordinates ID
  INTEGER           :: locID       !localization ID
  !Flag for assimilation (true) or re-initialization(false)
  LOGICAL           :: assimilate
END TYPE assim_ensemble_type

!---------------------------------------------------------------------------------------------------------------------------
! ENSEMBLE ENKF2MODEL INTERFACE DATA TYPE 
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_interface_type
  !Indices used to interface between ENKF data and Model data
  CHARACTER(LEN=30)    :: varName     !character string for model variables (used for debugging, and possibly file names and outputs)
  INTEGER              :: varID       !variable ID (id number used by interface for linking to model variables)  (in HYPE it can be an outvar index, or the order in which the state variables are considered by interface)
  INTEGER              :: modID       !model ID,  link to the corresponding variables used for H(X)              (in HYPE: outvar index)
  INTEGER              :: selID       !special index, to select 1st or 2nd set of varID:s and modID:s (ie. outvars or outvarbasinoutvars)
  INTEGER              :: nSubDim     !number of sub-dimensions (if needed, for instance lateral sub-units, vertical layers, substances, etc, in HYPE for instance SLC, substances, landuse, or slc, etc)
  INTEGER, ALLOCATABLE :: subDimID(:) !index in the sub-dimensions
END TYPE assim_interface_type
  
!---------------------------------------------------------------------------------------------------------------------------
! ENSEMBLE GENERATION INFORMATION DATA TYPE 
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_generation_type
  !number of variables
  INTEGER           :: nvar      !number of variables, ie n:o "model units" (for instance, number of sub-basins)
  INTEGER           :: ndata     !number of non-missing data
  !standard deviation (sigma) and mean (mean) used in ensemble generation
  REAL, ALLOCATABLE :: sigma(:)  !standard deviation [nvar x 1]
  REAL, ALLOCATABLE :: mean(:)   !mean               [nvar x 1]
  !parameter used for ensemble generation
  INTEGER           :: ensgen   ! type of ensemble generation        (0 none, 1 unrestricted, 2 [min,+inf], 3 [-inf,max], 4 restricted [min,max])   
  REAL              :: fixsigma ! fixed standard deviation           (ensgen=1)
  REAL              :: semimeta ! relative sigma for semi-restricted (ensgen=2,3,4, following Turner et al 2008)
  REAL              :: restmeta ! relative sigma for restricted      (ensgen=2,3,4)
  REAL              :: minsigma ! minimum sigma                      (ensgen=2,3,4)
  LOGICAL           :: ECMWF
  !variable for generation of spatially correlated random data, if needed
  LOGICAL           :: dorandxy
  TYPE(randxy_data) :: myrandxy_data
END TYPE assim_generation_type

!---------------------------------------------------------------------------------------------------------------------------
! STATE VARIABLE ENSEMBLE TYPE (one variable, all model units)
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_state_ensemble_type
  !ensemble data
  TYPE(assim_ensemble_type) :: x
  !interface data
  TYPE(assim_interface_type)  :: info
END TYPE assim_state_ensemble_type

!---------------------------------------------------------------------------------------------------------------------------
! ENSEMBLE TYPE USED FOR FORCING, OBSERVATION, and PARAMETER ENSEMBLES 
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_input_ensemble_type
  !ensemble matrix
  TYPE(assim_ensemble_type) :: x
  !interface data
  TYPE(assim_interface_type)  :: info
  !generation data
  TYPE(assim_generation_type) :: gen
END TYPE assim_input_ensemble_type

!---------------------------------------------------------------------------------------------------------------------------
! SPATIAL COORDINATES MATRIX VARIABLES TYPE 
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_coordinate_type
  INTEGER n                 !number of model positions
  REAL, ALLOCATABLE :: x(:) !x-coordinate
  REAL, ALLOCATABLE :: y(:) !y-coordinate
  REAL, ALLOCATABLE :: z(:) !z-coordinate
END TYPE assim_coordinate_type

!---------------------------------------------------------------------------------------------------------------------------
! LOCALIZATION MATRIX VARIABLES TYPE !CP161214 needed several locCXY-matrixes for different coordinate systems
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_twoDmatrix_type
  INTEGER n1,n2               !dimensions
  REAL, ALLOCATABLE :: x(:,:) !matrix
END TYPE assim_twoDmatrix_type

!---------------------------------------------------------------------------------------------------------------------------
! ASSIMILATION DATA TOP LEVEL DATA TYPE
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_data_type
  !general information about the assimilation application
  TYPE(assim_info_type)                           :: info
  !states, forcing, auxiliaries, parameters and observations ensembles
  TYPE(assim_state_ensemble_type), ALLOCATABLE :: X(:)     !model state     ensemble (vector over variables)
  TYPE(assim_input_ensemble_type), ALLOCATABLE :: F(:)     !model forcing   -"-       -"-    -"-  -"-
  TYPE(assim_state_ensemble_type), ALLOCATABLE :: A(:)     !model auxiliary -"-       -"-    -"-  -"-
  TYPE(assim_state_ensemble_type), ALLOCATABLE :: A2(:)     !model auxiliary -"-       -"-    -"-  -"-
  !TYPE(assim_input_ensemble_type), ALLOCATABLE :: P(:)     !model parameter -"-       -"-    -"-  -"-
  TYPE(assim_input_ensemble_type), ALLOCATABLE :: Obs(:)   !observation     -"-       -"-    -"-  -"-
  
  !observations, predicted observations, observation error covariance
     !TYPE(assim_ensemble_type)                    :: D        !observations    -"-      (incl all variables, used for next analysis, may span several timesteps)
     !TYPE(assim_ensemble_type)                    :: HX       !predictions(HX) -"-      (incl all variables, used for next analysis, may span several timesteps)
     !TYPE(assim_ensemble_type)                    :: R        !obs.error covar matrix   (incl all variables, used for next analysis, may span several timesteps)
  REAL, ALLOCATABLE :: D(:,:)
  REAL, ALLOCATABLE :: HX(:,:)
  REAL, ALLOCATABLE :: R(:)
  
 !spatial coordinates and localization data
 TYPE(assim_coordinate_type),     ALLOCATABLE :: Coordinates(:)
     !TYPE(assim_ensemble_type)                    :: LocCYY
     !TYPE(assim_ensemble_type),       ALLOCATABLE :: LocCXY(:)
  REAL, ALLOCATABLE :: LocCYY(:,:)
!  REAL, ALLOCATABLE :: LocCXY(:,:) !CP161214 need one LocCXY for every Coordinate system
  TYPE(assim_twoDmatrix_type),    ALLOCATABLE :: LocCXY(:)

END TYPE assim_data_type

!------------------------------------------------------------------------------
! 2 Declaration of the variables needed for the assimilation application
!------------------------------------------------------------------------------
!myAssimData with all necessary ensemble data etc
TYPE(assim_data_type), SAVE :: myAssimData
!---------------------------------------------------
!file id for input/output files
INTEGER fid_assim_info, fid_assim_min, fid_assim_max, fid_assim_median, fid_assim_final
INTEGER fid_assim_mean, fid_assim_min95, fid_assim_max95, fid_assim_std
!file id for binary direct access files
INTEGER, ALLOCATABLE :: fid_assim_bin(:)
INTEGER              :: fid_assim_bin_base
!defalt min and max values
REAL, PARAMETER :: assim_defmin = -1.1E38
REAL, PARAMETER :: assim_defmax = 1.1E38
REAL, PARAMETER :: assim_minsigma = 1.E-6
!
! Below, the old data structures
!

!!---------------------------------------------------------------------------------------------------------------------------
!! ENSEMBLE DATA TYPE
!!---------------------------------------------------------------------------------------------------------------------------
!TYPE enkf_ensemble_data
!  
!  !some of the fields are unnecessary for some type of ensembles, and will not be allocated
!  
!  INTEGER :: ne                 ! n:o ensembles (columns)
!  INTEGER :: n                  ! n:o variables (rows)
!  INTEGER :: ntypes             ! n:o types, only used for forcing and observation ensembles
!  INTEGER, ALLOCATABLE :: types(:) ! type ID
!  LOGICAL :: ECMWF              ! switch to select the ECMWF mean and std ensembles instead of PObs and TObs
!  
!  character(LEN=1) :: etype     ! Ensemble type: X, A, P, F, D, HX
!  
!  REAL, ALLOCATABLE :: y(:,:)   ! [n x ne] matrix
!  
!  INTEGER, ALLOCATABLE :: mi(:)            ! id number in relation to the model data, defined by the user in "interface"
!  INTEGER, ALLOCATABLE :: xi(:)            ! id number in ensemble matrix X (0 if not included)
!  INTEGER, ALLOCATABLE :: qt(:)            ! quality type (0 non-missing, -1 <=, -2 <, 1 >=, 2 >)
!  REAL,    ALLOCATABLE :: ql(:)            ! quality limit (if qt ~= 0, used for quality screening of observations)
!  INTEGER, ALLOCATABLE :: id(:)            ! additional id used by H operator to locate subid i 
!  
!  INTEGER, ALLOCATABLE :: kalman(:)        ! 0 no, 1 yes (switch Kalman filter on/off)
!  
!  INTEGER, ALLOCATABLE :: ensgen(:)        ! type of ensemble (0 none, 1 unrestricted, 2 semi-min, 3 semi-max, 4 restricted (min-max))   
!  REAL, ALLOCATABLE    :: mean(:)          ! mean value (used for initialization and generation)
!  REAL, ALLOCATABLE    :: minimum(:)       ! min allowed (for check)
!  REAL, ALLOCATABLE    :: maximum(:)       ! max allowed (for check)
!  
!  REAL, ALLOCATABLE    :: sigma(:)         ! stdev (used for generation)
!  REAL, ALLOCATABLE    :: semimeta(:)      ! relative sigma for semi-restricted (for generation a la Turner2008)
!  REAL, ALLOCATABLE    :: restmeta(:)      ! relative sigma for restricted (for generation a la Turner2008)
!  REAL, ALLOCATABLE    :: minsigma(:)      ! minimum sigma, used for variables with relative sigma (ensgen = 2,3, and 4)
!  REAL, ALLOCATABLE    :: covar(:,:)       ! covar (used for generation) - not needed!!! Only needed at assimilation - keep unallocated as long as possible
!    
!  !INTEGER, ALLOCATABLE :: statout(:)       ! yes(1) no(0) - print out statistics for this variable
!  !INTEGER, ALLOCATABLE :: ensout(:)        ! yes(1) no(0) - print out entire ensemble for this variable
!  
!  !some variables for output statistics:
!  REAL, ALLOCATABLE   :: outmean(:), outquant(:,:), outmin(:), outmax(:), outsigma(:) !, outcovar(:,:), outcorr(:,:)
!      
!  ! NEW variable for generation of spatially correlated random data, if needed
!  TYPE(randxy_data), ALLOCATABLE :: myrandxy_data(:)
!   
!end type enkf_ensemble_data
!
!! Ensemble declarations: X states, A additional model variables, P PARAMETERs, F forcing, D observations, HX model observ.
!type(enkf_ensemble_data)  :: EnkfX, EnkfA, EnkfP, EnkfF, EnkfD, EnkfHX
!type(enkf_ensemble_data)  :: EnkfXinfo, EnkfAinfo, EnkfPinfo, EnkfFinfo, EnkfDinfo
!
!!---------------------------------------------------------------------------------------------------------------------------
!! ENKF INFO: general information, settings, etc that could be modified by the user
!!---------------------------------------------------------------------------------------------------------------------------
!TYPE enkf_info
!    INTEGER nE         ! number of ensembles (columns in ensemble matrix X)
!    INTEGER nX         ! number of state variables in ensemble matrix X (excluding augmented variables)
!    INTEGER nA         ! number of additional model variables in A ensemble (used for output)
!    INTEGER nF         ! number of input in F ensemble
!    INTEGER nD         ! number of observations in D and HX ensembles (maximum size)
!    INTEGER nDA        ! number of observations in D and HX ensembles in current Assimilation Task
!    INTEGER nDT        ! number of observations per subbasin (HYPE specific)
!    INTEGER nP         ! number of parameters in P ensemble
!    INTEGER nAX        ! number of A in X matrix
!    INTEGER nFX        !  -"- of F in X
!    INTEGER nPX        !  -"- of P in X
!    
!    ! LOGICAL flags, may be modified...
!    LOGICAL FQ          ! include fluxes in kalman filter     (general switch on/off)
!    LOGICAL FA          ! include auxil. in kalman filter     (general switch on/off)
!    LOGICAL FP          ! include parameters in kalman filter (general switch on/off) 
!    LOGICAL FF          ! include inputs in kalman filter     (general switch on/off)
!    LOGICAL XS          ! include X states in ensemble statistics output (probably not...)
!    LOGICAL EC          ! ECMWF forecast switch
!    LOGICAL meanout     ! print ensemble mean (.true.) or median (.false.) in output files 
!    
!    ! some parameters
!    REAL :: moradkhani_delta ! coef.[0-1] to retain variance in parameter ensemble (Moradkhani et al, 2004), value typical around 0.95 
!
!    ! random number generation seed
!    INTEGER :: seed
!
!    ! missing value 
!    REAL :: missing
!
!    ! minimum allowed sigma in random number generation
!    REAL :: enkf_minsigma
!    
!    ! flags to indicate if outvars and parameters are included in ensembles (included in the enkf filer???)
!    REAL, ALLOCATABLE :: outvarflag(:), parameterflag(:)
!    
!    !option to read/write ensemble data to direct access binary files instead of saving in allocated memory
!    LOGICAL useBinFiles
!    
!end type enkf_info
!
!TYPE(enkf_info) EnkfInfo
!
!TYPE enkf_localization
!  REAL, ALLOCATABLE :: locCYY(:,:)
!  REAL, ALLOCATABLE :: locCXY(:,:)
!END TYPE enkf_localization
!
!TYPE(enkf_localization) EnkfLoc

! end of module
END MODULE ASSIMILATION_VARIABLES