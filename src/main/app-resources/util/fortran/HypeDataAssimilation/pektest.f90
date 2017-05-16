  
  MODULE VARS
    PUBLIC
    INTEGER,PARAMETER :: nsub=2
    INTEGER,PARAMETER :: nclass=4
  ENDMODULE VARS

  MODULE TYPES
    IMPLICIT NONE
    PUBLIC
    TYPE STATETYPE
      REAL,ALLOCATABLE :: glacvol(:)
      REAL,ALLOCATABLE :: snow(:,:)
      REAL,ALLOCATABLE :: csnow(:)
    END TYPE STATETYPE
    TYPE HELPSTATETYPE
      LOGICAL :: allok
      INTEGER :: ndim
      INTEGER :: dims(4)
      CHARACTER(LEN=10) :: name
    ENDTYPE
  
  CONTAINS
    SUBROUTINE point_to_current_state1(ivar,mypoint,states)
      INTEGER,INTENT(IN) :: ivar  !index of state variable e.g. snow
      REAL,POINTER,INTENT(INOUT) :: mypoint(:)
      TYPE(STATETYPE),TARGET,INTENT(IN) :: states
      
      IF(ivar==1)THEN
        mypoint => states%glacvol
      ELSEIF(ivar==2)THEN
        CONTINUE !ivar 2 is not 1dimensional
      ELSEIF(ivar==3)THEN
        mypoint => states%csnow
      ENDIF
    ENDSUBROUTINE point_to_current_state1
    
    SUBROUTINE point_to_current_state2(ivar,mypoint,states)
      INTEGER,INTENT(IN) :: ivar  !index of state variable e.g. snow
      REAL,POINTER,INTENT(INOUT) :: mypoint(:,:)
      TYPE(STATETYPE),TARGET,INTENT(IN) :: states
      
      IF(ivar==1)THEN
        CONTINUE !ivar 1 is not 2dimensional
      ELSEIF(ivar==2)THEN
        mypoint => states%snow
      ELSEIF(ivar==3)THEN
        CONTINUE !ivar 3 is not 2dimensional
      ENDIF
    ENDSUBROUTINE point_to_current_state2
  ENDMODULE TYPES 
  
MODULE ASSIMV
  IMPLICIT NONE
  PUBLIC
  real, parameter :: assim_defmin = -1.1E38
  real, parameter :: assim_defmax = 1.1E38
  TYPE assim_ensemble_type
    INTEGER           :: nvar           !number of variables, ie n:o "model units" (for instance, number of sub-basins)
    INTEGER           :: nens           !number of ensemble members
    REAL, ALLOCATABLE :: x(:,:)         !ensemble data [nvar x nens]
    REAL              :: minimum        !min allowed 
    REAL              :: maximum        !max allowed 
    REAL, ALLOCATABLE :: outmean(:), outquant(:,:), outmin(:), outmax(:), outsigma(:)
    INTEGER           :: coordID     !spatial coordinates ID
    INTEGER           :: locID       !localization ID
    LOGICAL           :: assimilate  !Flag for assimilation (true) or re-initialization(false)
  END TYPE assim_ensemble_type
  TYPE assim_interface_type
    !Indices used to interface between ENKF data and Model data
    CHARACTER(LEN=30)    :: varName     !character string for model variables (used for debugging, and possibly file names and outputs)
    INTEGER              :: varID       !variable ID (id number used by interface for linking to model variables)  (in HYPE it can be an outvar index, or the order in which the state variables are considered by interface)
    INTEGER              :: modID       !model ID,  link to the corresponding variables used for H(X)              (in HYPE: outvar index)
    INTEGER              :: nSubDim     !number of sub-dimensions (if needed, for instance lateral sub-units, vertical layers, substances, etc, in HYPE for instance SLC, substances, landuse, or slc, etc)
    INTEGER, ALLOCATABLE :: subDimID(:) !index in the sub-dimensions
  END TYPE assim_interface_type
  TYPE assim_state_ensemble_type
    TYPE(assim_ensemble_type) :: x  !ensemble data
    TYPE(assim_interface_type)  :: info  !interface data
  END TYPE assim_state_ensemble_type
  TYPE assim_info_type
     !Model structure variables
     integer nX           !number of model state variables     (dimension of assimData%X vector)
     integer nE           !number of ensemble members
     logical a_snow  !assimilation flags for different (HYPE model specific) variable categories
     logical meanout      !print ensemble mean (.true.) or median (.false.) in output files
  END TYPE assim_info_type
  TYPE assim_data_type
    TYPE(assim_info_type)                           :: info
    TYPE(assim_state_ensemble_type), allocatable :: X(:)     !model state     ensemble (vector over variables)
  END TYPE assim_data_type

  !Variable declaration
  TYPE(assim_data_type) :: assimData
    
ENDMODULE ASSIMV

MODULE ASSIMR
  USE ASSIMV
  IMPLICIT NONE
  PUBLIC
  
  CONTAINS
  SUBROUTINE allocate_0dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,minimum,maximum,assimilate)
    !INPUT
    TYPE(assim_state_ensemble_type),INTENT(INOUT)  :: assimVar
    integer,INTENT(IN) :: nens,nvar
    CHARACTER(LEN=*),INTENT(IN) :: varName
    real,intent(IN)    :: xini(nvar)
    integer,INTENT(INOUT) :: varID
    integer,intent(in) ::locID,coordID
    real,INTENT(IN)    :: minimum,maximum
    LOGICAL,INTENT(IN) :: assimilate
  
    !integer :: dimID(1),ndim
    !ndim=0
    !dimID(1)=0
    !allocate and initialize the ensemble matrix data
    call allocate_assim_ensemble(assimVar%x,nens,nvar,locID,coordID,xini,minimum,maximum,.true.,-9999.,assimilate)
    !allocate the interface data
    !call allocate_assim_interface(assimVar%info,varName,varID,-1,ndim,dimID)
    !update varID
    varID = varID+1
  END SUBROUTINE allocate_0dim_state_ensemble

  SUBROUTINE allocate_1dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,minimum,maximum,n1,assimilate)
    !INPUT
    TYPE(assim_state_ensemble_type),INTENT(INOUT)  :: assimVar(:)
    integer,intent(IN) :: nens,nvar
    integer, INTENT(INOUT) :: varID
    integer,intent(in) :: locID,coordID,n1
    real,intent(in)    :: minimum,maximum
    real,intent(in)    :: xini(n1,nvar)
    CHARACTER(LEN=*),intent(in) :: varName
    LOGICAL,intent(in) :: assimilate
   
    integer i,dimID(1),ndim
    ndim=1
    DO i=1,n1
      dimID(1)=i
      !allocate and initialize the ensemble matrix data
      call allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,locID,coordID,xini(i,:),minimum,maximum,.true.,-9999.,assimilate)
      !allocate the interface data
      !call allocate_assim_interface(assimVar(varID)%info,varName,varID,-1,ndim,dimID)
      varID = varID+1
    ENDDO
  END SUBROUTINE allocate_1dim_state_ensemble
  
  SUBROUTINE allocate_assim_ensemble(assimVar,nens,nvar,locID,coordID,xini,mini,maxi,allocateOutput,missing,assimilate)
    !INPUT ARGUMENTS
    TYPE(assim_ensemble_type),INTENT(INOUT) :: assimVar
    integer :: nens,nvar      !n:o ensemble members and variables (model units, i.e. n:o subbasins in HYPE) 
    integer :: locID
    integer :: coordID
    real    :: xini(nvar)     !initial values
    real    :: mini, maxi     !max/min thresholds for truncation of unrealistic values
    logical :: allocateOutput !flag for output variable allocation
    real    :: missing        !initial output values
    logical :: assimilate
    !LOCAL
    INTEGER :: j
    
    !ASSIGN AND ALLOCATE (make different if binary files are used)
    assimVar%nvar = nvar
    assimVar%nens = nens
    IF(ALLOCATED(assimVar%x))DEALLOCATE(assimVar%x)
    ALLOCATE(assimVar%x(nvar,nens))
    DO j=1,nens
      assimVar%x(:,j)=xini
    ENDDO
    assimVar%minimum = mini
    assimVar%maximum = maxi
    assimVar%locID = locID
    assimVar%coordID = coordID
  
    IF(allocateOutput)THEN
      IF(ALLOCATED(assimVar%outmean))DEALLOCATE(assimVar%outmean)
      ALLOCATE(assimVar%outmean(nvar))
      assimVar%outmean=missing
    !  IF(ALLOCATED(assimVar%outquant))DEALLOCATE(assimVar%outquant)
    !  ALLOCATE(assimVar%outquant(nvar,3))
    !  assimVar%outquant=missing
    !  IF(ALLOCATED(assimVar%outmin))DEALLOCATE(assimVar%outmin)
    !  ALLOCATE(assimVar%outmin(nvar))
    !  assimVar%outmin=missing
    !  IF(ALLOCATED(assimVar%outmax))DEALLOCATE(assimVar%outmax)
    !  ALLOCATE(assimVar%outmax(nvar))
    !  assimVar%outmax=missing
    !  IF(ALLOCATED(assimVar%outsigma))DEALLOCATE(assimVar%outsigma)
    !  ALLOCATE(assimVar%outsigma(nvar))
    !  assimVar%outsigma=missing
    ENDIF
  
    !assimVar%assimilate = assimilate
  END SUBROUTINE allocate_assim_ensemble

	SUBROUTINE write1D_to_X(x1d,ni,k,i_ens,assimVar)
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
		REAL, intent(IN) :: x1d(:)
		INTEGER, intent(IN) :: ni,i_ens
		INTEGER, intent(INOUT) :: k
		INTEGER i
		DO i=1,ni
      assimVar(k)%x%x(i,i_ens)=x1d(i)
		ENDDO
		k=k+1
	END SUBROUTINE write1D_to_X

	SUBROUTINE write2D_to_X(x2d,ni,nj,k,i_ens,assimVar)
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
		REAL, intent(IN) :: x2d(:,:)
		INTEGER, intent(IN) :: ni,nj, i_ens
		INTEGER, intent(INOUT) :: k
		INTEGER i,j
    DO i=1,ni
				DO j=1,nj
          assimVar(k)%x%x(j,i_ens) = x2d(i,j)
        ENDDO
        k=k+1
		ENDDO
		!k=k+ni
	END SUBROUTINE write2D_to_X
	SUBROUTINE writeX_to_1D(x1d,ni,k,i_ens,assimVar)
   TYPE(assim_state_ensemble_type)  :: assimVar(:)
		REAL, intent(OUT) :: x1d(:)
		INTEGER, intent(IN) :: ni,i_ens
		INTEGER, intent(INOUT) :: k
		INTEGER i
		DO i=1,ni
						x1d(i) = assimVar(k)%x%x(i,i_ens) !enkfX%y(k+i,i_ens)
		ENDDO
		k=k+1 !ni
	END SUBROUTINE writeX_to_1D
	
	SUBROUTINE writeX_to_2D(x2d,ni,nj,k,i_ens,assimVar)
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
		REAL, intent(OUT) :: x2d(:,:)
		INTEGER, intent(IN) :: ni,nj, i_ens
		INTEGER, intent(INOUT) :: k
		INTEGER i,j
		DO i=1,ni
			DO j=1,nj
				x2d(i,j) = assimVar(k)%x%x(j,i_ens)
      ENDDO
      k=k+1
		ENDDO
	END SUBROUTINE writeX_to_2D
	
	SUBROUTINE writeXmedian_to_1D(x1d,ni,k,assimVar,meanORmedian)
    LOGICAL,intent(in) :: meanORmedian
    TYPE(assim_state_ensemble_type),intent(in)  :: assimVar(:)
		REAL, intent(OUT) :: x1d(:)
		INTEGER, intent(IN) :: ni
		INTEGER, intent(INOUT) :: k
		INTEGER i
		DO i=1,ni
			IF(meanORmedian)THEN
				x1d(i) = assimVar(k)%x%outmean(i) !enkfX%outmean(k+i)
			ELSE
				x1d(i) = assimVar(k)%x%outquant(2,i) !enkfX%outquant(2,k+i)
      ENDIF
    ENDDO
    k=k+1
	END SUBROUTINE writeXmedian_to_1D

	SUBROUTINE writeXmedian_to_2D(x2d,ni,nj,k,assimVar,meanORmedian)
    LOGICAL,intent(in) :: meanORmedian
    TYPE(assim_state_ensemble_type),intent(in)  :: assimVar(:)
		REAL, intent(OUT) :: x2d(:,:)
		INTEGER, intent(IN) :: ni,nj
		INTEGER, intent(INOUT) :: k
		INTEGER i,j
		DO i=1,ni
			DO j=1,nj
				IF(meanORmedian)THEN
					x2d(i,j) = assimVar(k)%x%outmean(j) !enkfX%outmean(k+(i-1)*nj + j)
				ELSE
					x2d(i,j) = assimVar(k)%x%outquant(2,j) !enkfX%outquant(2,k+(i-1)*nj + j)
				ENDIF
      ENDDO
      k=k+1
		ENDDO
  END SUBROUTINE writeXmedian_to_2D

  SUBROUTINE updateEnsembleStatistics(assimData)
    TYPE(assim_data_type),intent(inOUT):: assimData
    integer i
    !state ensembles
    DO i=1,assimData%info%nX
      assimData%X(i)%x%outmean = MAXVAL(assimData%X(i)%x%x)
    ENDDO

  END SUBROUTINE updateEnsembleStatistics


ENDMODULE ASSIMR

MODULE ASSIMI
  USE VARS
  USE TYPES
  USE ASSIMV
  USE ASSIMR
  IMPLICIT NONE
  PUBLIC
  
  CONTAINS
  SUBROUTINE allocate_and_initialize_model_state_ensembles_HYPE(ne,nvar,assimVar,varID,assimInfo,mystates,stateinfo)
    
    !Argument declarations
    INTEGER, INTENT(IN) :: ne !Number of ensembles
    INTEGER, INTENT(IN) :: nvar !Number of state variables (different variables in HYPE, e.g. snow,glacvol,...)
    TYPE(assim_state_ensemble_type),INTENT(INOUT)  :: assimVar(:)
    TYPE(assim_info_type),INTENT(IN)            :: assimInfo
    INTEGER, INTENT(INOUT)           :: varID
    TYPE(statetype),INTENT(IN) :: mystates   !<Snow and ice states
    TYPE(HELPSTATETYPE),INTENT(IN) :: stateinfo(:)
  
    !Local variables
    INTEGER :: locID, coordID,ivar
    REAL,POINTER :: mypoint1(:),mypoint2(:,:)
    
    !3) allocate and initialize each variable ensemble
    locID=1
    coordID=1
    !nsub i anropen nedan kan bytas ut mot stateinfo(ivar)%dims(x)
   
    DO ivar=1,nvar
		  IF(stateinfo(ivar)%allok)THEN
        IF(stateinfo(ivar)%ndim==1)THEN
          CALL point_to_current_state1(ivar,mypoint1,mystates)
          CALL allocate_0dim_state_ensemble(assimVar(varID),ne,nsub,stateinfo(ivar)%name,mypoint1,varID,locID,coordID,0.,assim_defmax,assimInfo%a_snow)
        ELSEIF(stateinfo(ivar)%ndim==2)THEN
          CALL point_to_current_state2(ivar,mypoint2,mystates)
          CALL allocate_1dim_state_ensemble(assimVar(:),ne,nsub,stateinfo(ivar)%name,mypoint2,varID,locID,coordID,0.,assim_defmax,stateinfo(ivar)%dims(1),assimInfo%a_snow)
        ELSEIF(stateinfo(ivar)%ndim==3)THEN
		      !CALL allocate_2dim_state_ensemble(assimVar(:),ne,nsub,stateinfo%name(ivar),frozenstate%csnow,varID,locID,coordID,0.,assim_defmax,nclass,numsubstances,assimInfo%a_snow)
          !...
        ENDIF
      ENDIF
    ENDDO
    
  ENDSUBROUTINE allocate_and_initialize_model_state_ensembles_HYPE

  SUBROUTINE model_to_ensemble(i_ens,assimX,mystates,stateinfo)
    
    !Argument declarations
    INTEGER, INTENT(IN) :: i_ens !index of ensemble
    TYPE(assim_state_ensemble_type),INTENT(INOUT)  :: assimX(:)
    TYPE(statetype),INTENT(IN) :: mystates   !<Snow and ice states
    TYPE(HELPSTATETYPE),INTENT(IN) :: stateinfo(:)
  
    !Local variables
    INTEGER :: varID
    INTEGER :: nvar,ivar
    REAL,POINTER :: mypoint1(:),mypoint2(:,:)
    
    varID = 1 !was i original code
    nvar = 3   !comes from some info-variable somewhere (number of states in HYPE)
   
    DO ivar=1,nvar
		  IF(stateinfo(ivar)%allok)THEN
        IF(stateinfo(ivar)%ndim==1)THEN
          CALL point_to_current_state1(ivar,mypoint1,mystates)
          CALL write1D_to_X(mypoint1,stateinfo(ivar)%dims(1),varID,i_ens,assimX)
        ELSEIF(stateinfo(ivar)%ndim==2)THEN
          CALL point_to_current_state2(ivar,mypoint2,mystates)
          CALL write2D_to_X(mypoint2,stateinfo(ivar)%dims(1),stateinfo(ivar)%dims(2),varID,i_ens,assimX) 
        ELSEIF(stateinfo(ivar)%ndim==3)THEN
		      !CALL allocate_2dim_state_ensemble(assimVar(:),ne,nsub,stateinfo%name(ivar),frozenstate%csnow,varID,locID,coordID,0.,assim_defmax,nclass,numsubstances,assimInfo%a_snow)
        ENDIF
      ENDIF
    ENDDO
    
  ENDSUBROUTINE model_to_ensemble
    
  SUBROUTINE ensemble_to_model(i_ens,assimX,mystates,stateinfo)
    
    !Argument declarations
    INTEGER, INTENT(IN) :: i_ens !index of ensemble
    TYPE(assim_state_ensemble_type),INTENT(IN)  :: assimX(:)
    TYPE(statetype),INTENT(INOUT) :: mystates   !<Snow and ice states
    TYPE(HELPSTATETYPE),INTENT(IN) :: stateinfo(:)
  
    !Local variables
    INTEGER          :: varID
    INTEGER :: nvar,ivar
    REAL,POINTER :: mypoint1(:),mypoint2(:,:)
    
    varID = 1 !var k i ursprunglig kod
    nvar = 3   !comes from some info-variable somewhere
   
    DO ivar=1,nvar
		  IF(stateinfo(ivar)%allok)THEN
        IF(stateinfo(ivar)%ndim==1)THEN
          CALL point_to_current_state1(ivar,mypoint1,mystates)
          CALL writeX_to_1D(mypoint1,stateinfo(ivar)%dims(1),varID,i_ens,assimX)
        ELSEIF(stateinfo(ivar)%ndim==2)THEN
          CALL point_to_current_state2(ivar,mypoint2,mystates)
          CALL writeX_to_2D(mypoint2,stateinfo(ivar)%dims(1),stateinfo(ivar)%dims(2),varID,i_ens,assimX) 
        ELSEIF(stateinfo(ivar)%ndim==3)THEN
		      !CALL allocate_2dim_state_ensemble(assimVar(:),ne,nsub,stateinfo%name(ivar),frozenstate%csnow,varID,locID,coordID,0.,assim_defmax,nclass,numsubstances,assimInfo%a_snow)
        ENDIF
      ENDIF
    ENDDO
    
  ENDSUBROUTINE ensemble_to_model
    
	SUBROUTINE meanORmedian_to_model(assimX,meanORmedian,mystates,stateinfo)
  
    TYPE(assim_state_ensemble_type),intent(in)  :: assimX(:)
    logical,intent(in) :: meanORmedian
		TYPE(statetype),INTENT(INOUT) :: mystates   !<Snow and ice states
    TYPE(HELPSTATETYPE),INTENT(IN) :: stateinfo(:)

    ! local variables
		INTEGER ivar, varID, nvar
    REAL,POINTER :: mypoint1(:),mypoint2(:,:)
		
    varID = 1 !var k i ursprunglig kod
    nvar = 3   !comes from some info-variable somewhere
   
    DO ivar=1,nvar
		  IF(stateinfo(ivar)%allok)THEN
        IF(stateinfo(ivar)%ndim==1)THEN
          CALL point_to_current_state1(ivar,mypoint1,mystates)
		      CALL writeXmedian_to_1D(mypoint1,stateinfo(ivar)%dims(1),varID,assimX,meanORmedian)
        ELSEIF(stateinfo(ivar)%ndim==2)THEN
          CALL point_to_current_state2(ivar,mypoint2,mystates)
		      CALL writeXmedian_to_2D(mypoint2,stateinfo(ivar)%dims(1),stateinfo(ivar)%dims(2),varID,assimX,meanORmedian) 
        ELSEIF(stateinfo(ivar)%ndim==3)THEN
		!     CALL writeXmedian_to_3D(frozenstate%csnow,numsubstances,nclass,nsub,k,assimX,meanORmedian)
        ENDIF
      ENDIF
    ENDDO

	END SUBROUTINE meanORmedian_to_model

ENDMODULE ASSIMI

  
PROGRAM pektest
  USE VARS
  USE TYPES
  USE ASSIMV
  USE ASSIMI
  IMPLICIT NONE  
  !Variables and parameters
  INTEGER,PARAMETER :: nens=2
  INTEGER,PARAMETER :: nvar=3
  TYPE(STATETYPE),TARGET :: states  !I have made one, but there is six
  TYPE(HELPSTATETYPE),ALLOCATABLE :: stateinfo(:)   !One for all states?
  REAL,POINTER :: mypointer1(:), mypointer2(:,:)
  INTEGER varID !counter of all variables of different type and location
  INTEGER i,ienkf
  
  !Start program
  ALLOCATE(states%glacvol(nsub))
  ALLOCATE(states%snow(nclass,nsub))
  states%glacvol(:)=1.
  states%snow(:,1)=2.
  states%snow(:,2)=3.
  WRITE(*,*) 'snow',states%snow
  ALLOCATE(stateinfo(3))
  stateinfo(1:2)%allok = .TRUE.
  stateinfo(3)%allok = .FALSE.
  stateinfo(1)%ndim = 1
  stateinfo(2)%ndim = 2
  stateinfo(3)%ndim = 1
  stateinfo(1)%dims(1) = nsub
  stateinfo(2)%dims(1:2) = (/nclass,nsub/)
  stateinfo(3)%dims(1) = nsub
  !stateinfo%name = (/'glacvol','snow','xsnow'/)
  stateinfo(1)%name = 'glacvol'
  stateinfo(2)%name = 'snow'
  stateinfo(3)%name = 'xsnow'

  
  !Det här ligger egentligen i assim_interface
  assimData%info%nE=nens
  assimData%info%nX=1+nclass+0  
  assimData%info%a_snow=.TRUE.
  ALLOCATE(assimData%X(assimData%info%nx))
  
  !WRITE(*,*) stateinfo%name(2)
  varID = 1 !initialise
  CALL allocate_and_initialize_model_state_ensembles_HYPE(nens,nvar,assimData%X,varID,assimData%info,states,stateinfo)
  WRITE(*,*) 'Antal allokerade och to-be assimilerade variabler (olika):',varID
  WRITE(*,*) 'Förutberäknade antal variabler (lika):',assimData%info%nX
  WRITE(*,*) 'Värden i assimdata efter initialisering:'
  DO i=1,assimData%info%nX
    WRITE(*,*) assimData%X(i)%X%X
  ENDDO
  
  !Nollställ modelltillstånd och skriv sen över med X
  WRITE(*,*) 'Värden i tillstånden (kommer nollställas):'
  WRITE(*,*) states%glacvol
  WRITE(*,*) states%snow
  WRITE(*,*)
  states%glacvol(:)=0.
  states%snow(:,:)=0.
  WRITE(*,*) 'Värden i assimdata:'
  DO i=1,assimData%info%nX
    WRITE(*,*) assimData%X(i)%X%X
  ENDDO

  DO ienkf=1,nens
    CALL ensemble_to_model(ienkf,AssimData%X,states,stateinfo)
    !WRITE(*,*) states%glacvol
    !WRITE(*,*) states%snow
                
    !CALL model()
    !Change states
    IF(ienkf==1)THEN
      states%glacvol(:)=5.
      states%snow(:,:)=6.
    ELSEIF(ienkf==2)THEN
      states%glacvol(:)=3.
      states%snow(:,:)=7.
    ENDIF
    CALL model_to_ensemble(ienkf,AssimData%X,states,stateinfo)
    WRITE(*,*) 'Värden i assimdata efter ändring (model_to_ensemble) i ens:',ienkf
    DO i=1,assimData%info%nX
      WRITE(*,*) assimData%X(i)%X%X
    ENDDO
  ENDDO
  AssimData%info%meanout = .TRUE.
  WRITE(*,*) 'Maxvärden i assimdata (initierade, innan beräkning):'
  DO i=1,assimData%info%nX    
    WRITE(*,*) assimData%X(i)%x%outmean   !initialized to missing earlier
  ENDDO
  CALL updateEnsembleStatistics(assimData)
  WRITE(*,*) 'Maxvärden i assimdata:'
  DO i=1,assimData%info%nX
    WRITE(*,*) assimData%X(i)%x%outmean
  ENDDO
  CALL meanORmedian_to_model(AssimData%X,AssimData%info%meanout,states,stateinfo)
  WRITE(*,*) 'Tillståndens maxvärden (i ensemblen):'
  WRITE(*,*) states%glacvol
  WRITE(*,*) states%snow
 
  IF(.FALSE.)THEN
    !Tests of how pointer works
    mypointer1 => states%glacvol
    mypointer2 => states%snow
    WRITE(*,*) states%snow
    WRITE(*,*) states%glacvol
    WRITE(*,*) mypointer1
  
    !Change values of state through pointer is possible
    mypointer1(1)=5;mypointer1(2)=6
    WRITE(*,*) states%glacvol
    WRITE(*,*) mypointer1
  
    !DEALLOCATE(mypointer1) !does not work
    DEALLOCATE(states%glacvol)
    WRITE(*,*) ASSOCIATED(mypointer1)   !true
    WRITE(*,*) mypointer1 !still has values in ifort
    NULLIFY(mypointer1)   !also needed for pointer to forget the values
    WRITE(*,*) ASSOCIATED(mypointer1)   !false
    !End tests
  ENDIF
  
ENDPROGRAM
  
