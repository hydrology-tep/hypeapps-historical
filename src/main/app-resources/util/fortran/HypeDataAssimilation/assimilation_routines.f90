!> \file assimilation_routines.f90
!> Contains MODULE ASSIMILATION_ROUTINES, with model independent subroutines and functions used for the (EnKF) data assimilation in HYSS/HYPE.
!> Author: D.Gustafsson (SMHI)
!> Versions: 
!> 2011-2015:  Developed by D.Gustafsson (SMHI/KTH) and J.Ahlberg (KTH). The latest version from 2015-09-20 was used with the HOPE model for a "HUVA" snow data assimilation project 2013-2015.
!> 2016.10.20: Module name changed to ASSIMILATION_ROUTINES and additional code convention adaptations for new implementation in HYPE_4_12_+.
!> 2017.04.27: Replace all code dependent on Intel MKL library with intrinsic functions and other new/old code.
  
!> Generic subroutines and functions for (Ensemble Kalman Filter) Data Assimilation.

MODULE ASSIMILATION_ROUTINES
!Copyright 2016 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------------------
!
! assimilation_routines.f90
!
! Module with EnKF subroutines and functions
!
! These functions are supposed to be general, and should not be changed
! 
! Functions special for a specific model application are in assimilation_interface.f90
!
!--------------------------------------------------------------------------------------------------

  !Use ASSIMILATION_VARIABLES and RANDOM_ROUTINES
  USE ASSIMILATION_VARIABLES
  USE RANDOM_ROUTINES
  
  IMPLICIT NONE

CONTAINS
!-----------------------------------------------------------------------------------
!VARIOUS ENSEMBLE DATA MANIPULATION ROUTINES
!-----------------------------------------------------------------------------------
  
!<\brief re-initilize selected ensembles to mean or median to avoid numerical issues
  SUBROUTINE meanORmedian_to_ensemble(nx,assimX,meanORmedian)

    !Argument declations
    TYPE(assim_state_ensemble_type) :: assimX(:)
    INTEGER, INTENT(IN)       :: nx           !<number of variables
    LOGICAL, INTENT(IN)       :: meanORmedian !<flag for setting mean (T) or median (F)

    INTEGER i,j,k
    !loop over state ensembles
    DO k=1,nx
      IF(.not.assimX(k)%x%assimilate)THEN
        !loop over ensemble members
        DO j=1,assimX(k)%x%nens
          IF(ALLOCATED(assimX(k)%x%x))THEN
            !if matrix allocated...
            !loop over model units (often subbasins)
            DO i=1,assimX(k)%x%nvar
              IF(meanORmedian)THEN
                assimX(k)%x%x(i,j) = assimX(k)%x%outmean(i)
              ELSE
                ! switch to "outmedian"
                assimX(k)%x%x(i,j) = assimX(k)%x%outquant(2,i)
              ENDIF
            ENDDO
          ELSE
            !else..  write model units to bin-file, mean or median
            IF(meanORmedian)THEN
              WRITE(assimX(k)%x%fileID,REC=assimX(k)%x%rec+j) assimX(k)%x%outmean(:)
            ELSE
                ! switch to "outmedian"
              WRITE(assimX(k)%x%fileID,REC=assimX(k)%x%rec+j) assimX(k)%x%outquant(2,:)
            ENDIF
         ENDIF
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE meanORmedian_to_ensemble
    
!>Initialize the TYPE(enkf_info) EnkfInfo variable
!>
!---------------------------------------------------------
  SUBROUTINE initialize_assim_info(assimInfo)
    !INPUT ARGUMENT
    TYPE(assim_info_type), INTENT(INOUT) :: assimInfo
    !LOCAL VARIABLE
    INTEGER :: errcode
  
    !Initialize random number generation seed (NOT USED ANYMORE!!!)
    !CALL time(assimInfo%seed)    !the random seed will be different, maybe make it possible for user to decide via the enkfinfo.txt?
  
    !Initialize various info variables
    assimInfo%nE      = 100       !number of ensemble members      (columns in ensemble matrices)
    assimInfo%nX      = 0         !number of state variables       (length of ensemble vector EnkfX)
    assimInfo%nA      = 0         !number of auxiliary variables   (length of ensemble vector EnkfA)
    assimInfo%nF      = 0         !number of forcing variables     (length of ensemble vector EnkfF)
    assimInfo%nObs    = 0         !number of observation variables (rows of ensemble vector EnkfObs(:,:))
    assimInfo%nD      = 0         !number of observations for next Enkf analysis (rows in EnkfD%x%x(:,:))
  !  assimInfo%nDT     = 0         !number of observation timesteps (colss of ensemble vector EnkfObs(:,:))
  !  assimInfo%nP      = 0         !number of parameters            (length of ensemble vector EnkfP)
    assimInfo%nloc    = 0         !number of localization matrices (length of vector EnkfLocCXY)
    assimInfo%ncoord  = 0         !number of spatial domains       (length of vector EnkfCoordinates)
    
    !LOGICAL flags, may be modified...
    assimInfo%FA      = .false.   !include auxil. in kalman filter     (general switch on/off)
    assimInfo%FP      = .false.   !include parameters in kalman filter (general switch on/off) 
    assimInfo%FF      = .false.   !include inputs in kalman filter     (general switch on/off)
    !EnkfInfo%XS      = .false.   !include X states in statistical output (general switch on/off)
    !assimInfo%EC      = .false.   !ECMWF forecasts (general switch on/off)
    assimInfo%meanout = .true.    !ensemble mean(.true.) or median (.false.) in output
    ALLOCATE(assimInfo%assim_flag(assimInfo%nCat+50)) !assimilate state variables (switch for categories and variables) (the size is unecessary large I think)
    assimInfo%a_snow = .true.
    assimInfo%a_glacier = .true.
    assimInfo%a_lakeice = .true.
    assimInfo%a_riverice = .true.
    assimInfo%a_soil = .true.
    assimInfo%a_riverwt = .true.
    assimInfo%a_lakewt = .true.
    assimInfo%a_misc = .true.
    assimInfo%a_aquifer = .true.

    !some Ensemble Kalmna Filter parameters
    assimInfo%moradkhani_delta = 0.95 ! coef.[0-1] to retain variance in parameter ensemble (Moradkhani et al, 2004), value typical around 0.95 
    !assimInfo%ensgen_minsigma = 1.e-5
    
    !switch on/off read/write ensemble data to binary files
    assimInfo%useBinFilesX = .FALSE.
    assimInfo%useBinFilesFA = .FALSE.
    assimInfo%nBinFiles = 0
  
    !Output of statistical simulation results
    assimInfo%nstatout = 0

    !localization parameters
    assimInfo%xy_scalefac = 1000000.
    assimInfo%z_scalefac  = 100000.
  END SUBROUTINE initialize_assim_info

!>Allocate Ensemble vectors needed for the data assimiltion application
!>
!----------------------------------------------------------------------------------------------  
  SUBROUTINE allocate_assim_ensemble_vectors(assimData,fid_0)
    TYPE(assim_data_type), INTENT(INOUT) :: assimData
    !INTEGER :: nx,na,np,nf,nobs,ncoord,nloc,nobsDT
    INTEGER :: i,fid_0
    !--------------------------------------------------------------------------------------------
    !X state variable ensemble vector
    IF(ALLOCATED(assimData%X))DEALLOCATE(assimData%X)
    IF(assimData%info%nx.GT.0)ALLOCATE(assimData%X(assimData%info%nx))
    
    !A auxiliary variable ensemble vector
    IF(ALLOCATED(assimData%A))DEALLOCATE(assimData%A)
    IF(assimData%info%na.GT.0)ALLOCATE(assimData%A(assimData%info%na))
    
    !  !P parameter ensemble vector
    !  IF(ALLOCATED(assimData%P))DEALLOCATE(assimData%P)
    !  IF(assimData%info%np.GT.0)ALLOCATE(assimData%P(assimData%info%np))
    
    !Obs observation ensemble vector
    IF(ALLOCATED(assimData%Obs))DEALLOCATE(assimData%Obs)
    !  IF(assimData%info%nobs.GT.0)ALLOCATE(assimData%Obs(nobs,nobsDT))
    IF(assimData%info%nobs.GT.0)ALLOCATE(assimData%Obs(assimData%info%nobs))
    
    !F forcing variable ensemble vector
    IF(ALLOCATED(assimData%F))DEALLOCATE(assimData%F)
    IF(assimData%info%nf.GT.0)ALLOCATE(assimData%F(assimData%info%nf))
    
    !Spatial Coordinate Domain vector
    !IF(ALLOCATED(assimData%Coordinates))DEALLOCATE(assimData%Coordinates)
    !IF(assimData%info%ncoord.GT.0)ALLOCATE(assimData%Coordinates(assimData%info%ncoord))
    
    !Localiztion (locCXY) vector
    !IF(ALLOCATED(assimData%LocCXY))DEALLOCATE(assimData%LocCXY)
    !IF(assimData%info%nloc.GT.0)ALLOCATE(assimData%LocCXY(assimData%info%nloc))

    !D, HX and R and LocCYY are not vectors, and DO not need to be ALLOCATED, the length is decided for each analysis time window
    !Further more, R is not an ensemble but an nD x nD matrix, and if observation errors are assumed uncorrelated, R is in fact a diagonal matrix and can be saved as a nD x 1 vector
  
    !If binary files are used
    !  IF(assimData%info%useBinFiles)THEN
      !nState + nAux + nPar + nForc + nObs*nObsTimeSteps + 3 (D + HX + R)
      !assimData%info%nBinFiles = assimData%info%nX + assimData%info%nA + assimData%info%nP + assimData%info%nF + assimData%info%nObs*assimData%info%nObsDT + 3
 
    assimData%info%nBinFiles = assimData%info%nX + assimData%info%nA + assimData%info%nF + assimData%info%nObs + 3
    
    IF(ALLOCATED(fid_assim_bin))DEALLOCATE(fid_assim_bin)
    ALLOCATE(fid_assim_bin(assimData%info%nBinFiles))
    
    DO i=1,assimData%info%nBinFiles
      fid_assim_bin(i) = fid_0 + i 
    ENDDO
      !ENDIF
    END SUBROUTINE allocate_assim_ensemble_vectors

  !--------------------------------------------------------
  !>Allocate and initialize an assim_ensemble_type variable
  !--------------------------------------------------------
  !SUBROUTINE allocate_assim_ensemble(assimVar, nens, nvar, fileID, locID, coordID, xini, mini, maxi, allocateOutput, missing,assimilate) !CP161201 added bin-files here
  SUBROUTINE allocate_assim_ensemble(assimVar, nens, nvar, varID, fileID, useBinFile, locID, coordID, xini, mini, maxi, allocateOutput, missing,assimilate)
    !INPUT ARGUMENTS
    TYPE(assim_ensemble_type) :: assimVar
    INTEGER :: nens,nvar      !n:o ensemble members and variables (model units, i.e. n:o subbasins in HYPE) 
    INTEGER, INTENT(IN) :: varID  !to set record for binary file
    INTEGER :: fileID         !file unit ID for direct access binary file I/O
    INTEGER :: locID
    INTEGER :: coordID
    REAL    :: xini(nvar)     !initial values
    REAL    :: mini, maxi     !max/min thresholds for truncation of unrealistic values
    LOGICAL :: allocateOutput !flag for output variable allocation
    REAL    :: missing        !initial output values
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation

    !LOCAL
    INTEGER :: j,reclen
    CHARACTER(LEN=10) :: filename
    
    !ASSIGN AND ALLOCATE or open bin-files for saving ensembles (CP added bin-files CP161201)
    assimVar%nvar = nvar
    assimVar%nens = nens
    assimVar%minimum = mini
    assimVar%maximum = maxi
    assimVar%fileID = fileID
    assimVar%locID = locID
    assimVar%coordID = coordID
    IF(useBinFile==1)THEN
      !If one bin-file is used, write initial data to already open file
      assimVar%rec = (varID-1)*nens
      DO j=1,nens
        WRITE(assimVar%fileID,REC=assimVar%rec+j) xini
      ENDDO
    ELSEIF(useBinFile==2)THEN
      !If several bin-files is used open file and write initial data
      INQUIRE(IOLENGTH=reclen) xini   !determine suitable record length here, because ifort och gfortran have different file storage unit (i.e. RECL)
      WRITE(filename,'(I5.5,a4)') assimVar%fileID,'.bin'
      OPEN(UNIT=assimVar%fileID,FILE=TRIM(filename),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=reclen)
      DO j=1,nens
        WRITE(assimVar%fileID,REC=j) xini
      ENDDO
      assimVar%rec = 0
    ELSE
      !If no bin-files, then allocate a matrix to hold the ensemble and initialize it.
      assimVar%rec = 0
      IF(ALLOCATED(assimVar%x))DEALLOCATE(assimVar%x)
      ALLOCATE(assimVar%x(nvar,nens))
      DO j=1,nens
        assimVar%x(:,j)=xini
      ENDDO
    ENDIF
  
    IF(allocateOutput)THEN
      IF(ALLOCATED(assimVar%outmean))DEALLOCATE(assimVar%outmean)
      ALLOCATE(assimVar%outmean(nvar))
      assimVar%outmean=missing
      IF(ALLOCATED(assimVar%outquant))DEALLOCATE(assimVar%outquant)
      !ALLOCATE(assimVar%outquant(nvar,3))    !CP161208 This is opposit how the array is used in the code!
      ALLOCATE(assimVar%outquant(3,nvar))
      assimVar%outquant=missing
      IF(ALLOCATED(assimVar%outmin))DEALLOCATE(assimVar%outmin)
      ALLOCATE(assimVar%outmin(nvar))
      assimVar%outmin=missing
      IF(ALLOCATED(assimVar%outmax))DEALLOCATE(assimVar%outmax)
      ALLOCATE(assimVar%outmax(nvar))
      assimVar%outmax=missing
      IF(ALLOCATED(assimVar%outsigma))DEALLOCATE(assimVar%outsigma)
      ALLOCATE(assimVar%outsigma(nvar))
      assimVar%outsigma=missing
    ENDIF
  
    assimVar%assimilate = assimilate
  END SUBROUTINE allocate_assim_ensemble

!--------------------------------------------------------
!>Allocate and initialize a assim_interface_type variable
!--------------------------------------------------------
  SUBROUTINE allocate_assim_interface(assimVar,varName,varID,modID,nSubDim,subDimID)
    !INPUT ARGUMENTS
    TYPE(assim_interface_type)  :: assimVar
    CHARACTER(LEN=*)    :: varName        !character string for model variables (used for debugging, and possibly file names and outputs)
    INTEGER              :: varID          !variable ID (id number used by interface for linking to model variables)  (in HYPE it can be an outvar index, or the order in which the state variables are considered by interface)
    INTEGER              :: modID          !model ID,  link to the corresponding variables used for H(X)              (in HYPE: outvar index)
    INTEGER              :: nSubDim        !number of sub-dimensions (if needed, for instance lateral sub-units, vertical layers, substances, etc, in HYPE for instance SLC, substances, landuse, or slc, etc)
    INTEGER              :: subDimID(:)    !index in the sub-dimensions
    !ALLOCATE AND ASSIGN DATA
    assimVar%varName(1:(MIN(30,LEN_TRIM(varName)))) = TRIM(varName)  
    assimVar%varID   = varID
    assimVar%modID   = modID
    assimVar%nSubDim = nSubDim
    IF(ALLOCATED(assimVar%subDimID))DEALLOCATE(assimVar%subDimID)
    IF(nSubDim.GT.0)THEN
      ALLOCATE(assimVar%subDimID(nSubDim))
      assimVar%subDimID(1:nSubDim)=subDimID(1:nSubDim)
    ENDIF
  END SUBROUTINE allocate_assim_interface

!---------------------------------------------------------
!>Allocate and initialize an enkf_generation_type variable
!---------------------------------------------------------
  SUBROUTINE allocate_assim_generation(assimVar,nvar,ensgen,fixsigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord) !,ECMWF)
    !INPUT ARGUMENTS
    TYPE(assim_generation_type) :: assimVar !generation_data variable
    INTEGER                     :: nvar     !number of variables, ie n:o "model units" (for instance, number of sub-basins)
    INTEGER                     :: ensgen   ! type of ensemble generation        (0 none, 1 unrestricted, 2 [min,+inf], 3 [-inf,max], 4 restricted [min,max])   
    REAL                        :: fixsigma ! fixed standard deviation           (ensgen=1)
    REAL                        :: semimeta ! relative sigma for semi-restricted (ensgen=2,3,4, following Turner et al 2008)
    REAL                        :: restmeta ! relative sigma for restricted      (ensgen=2,3,4)
    REAL                        :: minsigma ! minimum sigma                      (ensgen=2,3,4)
    REAL                        :: lscale   !correlation length for spatially correlated perturbation
    REAL                        :: gridsize !2D gridsize for spatially correlated perturbation generation
    INTEGER                     :: corrtype !spatial correlation FUNCTION option
    REAL                        :: xcoord(:)!x coordinate, for spatially correlated perturbations
    REAL                        :: ycoord(:)!y coordinate, for spatially correlated perturbations
    !LOGICAL                     :: ECMWF
    
    !ASSIGN
    assimVar%nvar = nvar
    assimVar%ensgen = ensgen
    assimVar%fixsigma = fixsigma
    assimVar%semimeta = semimeta
    assimVar%restmeta =restmeta
    assimVar%minsigma = minsigma
    
    !ALLOCATE
    !sigma, standard deviation (input to ensemble generation)
    IF(ALLOCATED(assimVar%sigma))DEALLOCATE(assimVar%sigma)
    ALLOCATE(assimVar%sigma(nvar))
    assimVar%sigma = 0.0
    IF(ensgen.EQ.1)assimVar%sigma=fixsigma
    
    !mean, mean value  (input to ensemble generation)
    IF(ALLOCATED(assimVar%mean))DEALLOCATE(assimVar%mean)
    ALLOCATE(assimVar%mean(nvar))
    assimVar%mean = 0.0
    
    !data structure for spatially correlated random perturbations
    IF(lscale.GT.0. .AND. gridsize.GT.0. .AND. corrtype.GT.0 .AND. corrtype.LT.4)THEN
      !DO randomxy
      assimVar%dorandxy = .TRUE.
      !initialize randxy data structure
      CALL init_randxy_data(assimVar%myrandxy_data,nvar,xcoord,ycoord,lscale,gridsize,corrtype)
      !    assimVar%myrandxy_data%lscale     = lscale
      !    assimVar%myrandxy_data%gridsize   = gridsize
      !    assimVar%myrandxy_data%index_corr = corrtype
    ELSE
      assimVar%dorandxy = .FALSE.
    ENDIF
    
    !assimVar%ECMWF = ECMWF
    
  END SUBROUTINE allocate_assim_generation

!-----------------------------------------------------------------
!>Allocate and initialize STATE (and auxiliary) ensemble variables
!-----------------------------------------------------------------
!>0-dimensional
  SUBROUTINE allocate_auxiliary_ensemble(assimVar,nens,nvar,varName,xini,varID,recID,locID,coordID,fileID,useBinFile,minimum,maximum,assimilate)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar
    INTEGER :: nens,nvar,locID,coordID,fileID,dimID(1),ndim
    INTEGER,INTENT(IN) :: varID   !<outvar-variable index
    INTEGER,INTENT(IN) :: recID   !<record of variable in bin-file
    REAL    :: minimum,maximum
    REAL    :: xini(nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation
  
    ndim=0
    dimID(1)=0
    !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
    CALL allocate_assim_ensemble(assimVar%x,nens,nvar,recID,fileID,useBinFile,locID,coordID,xini,minimum,maximum,.true.,-9999.,assimilate)
    !ALLOCATE the interface data
    CALL allocate_assim_interface(assimVar%info,varName,varID,-1,ndim,dimID)
    !update varID
    !varID = varID+1 !not for aux this is o_rout etc.
  END SUBROUTINE allocate_auxiliary_ensemble
  
!>0-dimensional
  SUBROUTINE allocate_0dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,useBinFile,minimum,maximum,assimilate)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar
    INTEGER :: nens,nvar,varID,locID,coordID,fileID,dimID(1),ndim
    REAL    :: minimum,maximum
    REAL    :: xini(nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation
  
    ndim=0
    dimID(1)=0
    !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
    !CALL allocate_assim_ensemble(assimVar%x,nens,nvar,fileID,locID,coordID,xini,minimum,maximum,.true.,-9999.,assimilate) !CP161201 introduced binfiles
    CALL allocate_assim_ensemble(assimVar%x,nens,nvar,varID,fileID,useBinFile,locID,coordID,xini,minimum,maximum,.true.,-9999.,assimilate)
    !ALLOCATE the interface data
    CALL allocate_assim_interface(assimVar%info,varName,varID,-1,ndim,dimID)
    !update varID
    varID = varID+1
  END SUBROUTINE allocate_0dim_state_ensemble
  
!>1-dimensional
  SUBROUTINE allocate_1dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,useBinFile,minimum,maximum,n1,assimilate)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    INTEGER :: nens,nvar,varID,locID,coordID,fileID(:),n1,i,dimID(1),ndim
    REAL    :: minimum,maximum
    REAL    :: xini(n1,nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation

    ndim=1
    DO i=1,n1
      dimID(1)=i
      !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
      !CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,fileID(varID),locID,coordID,xini(i,:),minimum,maximum,.true.,-9999.,assimilate) !CP161201 introduced binfiles
      CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,varID,fileID(varID),useBinFile,locID,coordID,xini(i,:),minimum,maximum,.true.,-9999.,assimilate)
      !ALLOCATE the interface data
      CALL allocate_assim_interface(assimVar(varID)%info,varName,varID,-1,ndim,dimID)
      !update varID
      varID = varID+1
    ENDDO
  END SUBROUTINE allocate_1dim_state_ensemble
  
!>2-dimensional (real)
  SUBROUTINE allocate_2dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,useBinfile,minimum,maximum,n1,n2,assimilate)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    INTEGER :: nens,nvar,varID,locID,coordID,fileID(:),n1,n2,i,j,ndim,dimID(2)
    REAL    :: minimum,maximum
    REAL    :: xini(n2,n1,nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation

    ndim=1
    DO j=1,n2
      dimID(1)=j
    DO i=1,n1
      dimID(2)=i
      !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
      !CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,fileID(varID),locID,coordID,xini(j,i,:),minimum,maximum,.true.,-9999.,assimilate) !CP161201 introduced binfiles
      CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,varID,fileID(varID),useBinFile,locID,coordID,xini(j,i,:),minimum,maximum,.true.,-9999.,assimilate)
      !ALLOCATE the interface data
      CALL allocate_assim_interface(assimVar(varID)%info,varName,varID,-1,ndim,dimID)
      !update varID
      varID = varID+1
    ENDDO
    ENDDO
  END SUBROUTINE allocate_2dim_state_ensemble
  
!>3-dimensional
  SUBROUTINE allocate_3dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,useBinFile,minimum,maximum,n1,n2,n3,assimilate)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    INTEGER :: nens,nvar,varID,locID,coordID,fileID(:),n1,n2,n3,i,j,k,ndim,dimID(3)
    REAL    :: minimum,maximum
    REAL    :: xini(n3,n2,n1,nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-files
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation

    ndim=3
    DO k=1,n3
      dimID(1)=k
    DO j=1,n2
      dimID(2)=j
    DO i=1,n1
      dimID(3)=i
      !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
      !CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,fileID(varID),locID,coordID,xini(k,j,i,:),minimum,maximum,.true.,-9999.,assimilate) !CP161201 introduced binfiles
      CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,varID,fileID(varID),useBinFile,locID,coordID,xini(k,j,i,:),minimum,maximum,.true.,-9999.,assimilate)
      !ALLOCATE the interface data
      CALL allocate_assim_interface(assimVar(varID)%info,varName,varID,-1,ndim,dimID)
      !update varID
      varID = varID+1
    ENDDO
    ENDDO
    ENDDO
  END SUBROUTINE allocate_3dim_state_ensemble
  
!>FORCING (and parameter) ensembles
  SUBROUTINE allocate_assim_forcing_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,minimum,maximum, &
    ensgen,sigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord,useBinFile)
    !ensgen,sigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord) !,ECMWF)
    !INPUT
    TYPE(assim_input_ensemble_type) :: assimVar
    INTEGER :: nens,nvar,locID,coordID,fileID
    INTEGER,INTENT(IN) :: varID   !<counter for variables in ensembles (is not all variables varID for F-variables but hardkoded)
    REAL    :: minimum,maximum
    REAL    :: xini(nvar)
    CHARACTER(LEN=*) :: varName
    REAL    :: lscale,gridsize,xcoord(:),ycoord(:)
    INTEGER :: corrtype,ensgen
    REAL    :: sigma, semimeta, restmeta, minsigma
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
  !  LOGICAL :: ECMWF
    !LOCAL VARIABLES
    INTEGER :: nDim,dimID(1),modID
  
    nDim=0;dimID(1)=0
    
    !ALLOCATE and initialize the ensemble matrix data (NOT initialize bin-file for now)
    !CALL allocate_assim_ensemble(assimVar%x,nens,nvar,fileID,locID,coordID,xini,minimum,maximum,.true.,-9999.,.false.) !CP161205 for binfiles and assimilate turned on to apply enkf
    CALL allocate_assim_ensemble(assimVar%x,nens,nvar,varID,fileID,useBinFile,locID,coordID,xini,minimum,maximum,.true.,-9999.,.true.)

    !ALLOCATE the interface data
    nDim = 0 ; dimID(1) = 0 ; modID = -9999
    CALL allocate_assim_interface(assimVar%info,varName,varID,modID,nDim,dimID)

    !ALLOCATE and initialize the ensemble generation data
    CALL allocate_assim_generation(assimVar%gen,nvar,ensgen,sigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord) !, ECMWF)
  END SUBROUTINE allocate_assim_forcing_ensemble

!>OBSERVATION ENSEMBLES
  SUBROUTINE allocate_assim_observation_ensemble(assimVar,nens,nvar,varName,xini,obsID,modID,coordID,fileID,minimum,maximum, &
    ensgen,sigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord)
    !INPUT
    TYPE(assim_input_ensemble_type) :: assimVar
    INTEGER :: nens,nvar
    REAL    :: xini(nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER :: obsid, modID
    INTEGER coordID,fileID
    REAL    :: minimum,maximum
    REAL    :: lscale,gridsize
    INTEGER :: corrtype,ensgen
    REAL    :: sigma, semimeta, restmeta, minsigma
    REAL    :: xcoord(:),ycoord(:)
    !local
    INTEGER :: nDim,dimID(1)

    !ALLOCATE and initialize the ensemble matrix data (NOT initialize bin-file for now)
    !CALL allocate_assim_ensemble(assimVar%x,nens,nvar,fileID,0,coordID,xini,minimum,maximum,.true.,-9999.,.false.) !CP161205 for binfiles and assimilate turned on to apply enkf
    CALL allocate_assim_ensemble(assimVar%x,nens,nvar,obsid,fileID,0,0,coordID,xini,minimum,maximum,.true.,-9999.,.true.)
    !ALLOCATE the interface data
    nDim = 0 ; dimID(1) = 0
    CALL allocate_assim_interface(assimVar%info,varName,obsID,modID,nDim,dimID)
    !ALLOCATE and initialize the ensemble generation data
    CALL allocate_assim_generation(assimVar%gen,nvar,ensgen,sigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord) !,.false.)

  END SUBROUTINE allocate_assim_observation_ensemble

!>Truncate an ensemble matrix to minimum and maximum allowed values.
!----------------------------------------------------------------------------------------------
SUBROUTINE assim_checkminmax(nx,ne,ensemble,minval,maxval)
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ne
  REAL, INTENT(INOUT) :: ensemble(nx,ne)
  REAL, INTENT(IN)    :: minval, maxval
  INTEGER :: i,j
  DO i=1,ne
    DO j=1,nx
      !make sure missing values stay missing
      IF(ensemble(j,i).GT.-9998.)THEN
        ensemble(j,i)=AMAX1(minval,AMIN1(maxval,ensemble(j,i)))
      ELSE
        ensemble(j,i)=-9999.
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE assim_checkminmax
  
!>General routine for Ensemble generation (forcing and observation data).
!>
!> The ensemble generation is made by adding random numbers to the input data.
!>
!> The basic assumption is that the random perturbations are gaussian with zero mean and standard deviation sigma.
!>
!> However, based on Turner et al(2008), it is assumed that in order to get unbiased input ensembles, 
!> we need to consider two types of perturbations - systematic and random:
!> Thus, the input x at time k for ensemble member j, x_kj = x_k0+eata_kj+chi_j,
!> where eata_kj is regenerated every time step and chi_j is generated only once.
!> Then, eata and chi can be generated with three different types of restrictions on variables.
!> For the moment, we assume that the static error chi_j = 0.
!>
!> In adition, we now also take into account spatial correlation in the data, by generating 
!> spatially correlated random data using a FFT-based method.
!----------------------------------------------------------------------------------------------
SUBROUTINE generate_input_ensemble(n,assimVar)
  INTEGER n
  TYPE(assim_input_ensemble_type) :: assimVar(n)
  INTEGER i,j
  REAL midvalue
  REAL,ALLOCATABLE :: localx(:,:)
    
  !loop over ensemble vector members
  DO j=1,n

    !get data to a local array for manipulation
    ALLOCATE(localx(assimVar(j)%x%nvar,assimVar(j)%x%nens))

    !before ensemble generation, check that number of ensemble members is greater than 1   
    IF(assimVar(j)%x%nens.gt.1)THEN 

      !a) determine standard deviation (sigma) in each model unit (rows in ensemble matrix)
        
      !loop over number of model units (subbasins in HYPE)
      DO i=1,assimVar(j)%x%nvar
!        ! select General case or special "ECMWF" case
!        IF(.not. ensemble(j)%gen%ECMWF)THEN
        
        !Check missing values
        IF(assimVar(j)%gen%mean(i).GT.-9998.)THEN
          SELECT CASE(assimVar(j)%gen%ensgen)
            CASE(4) ! restricted (Turner et al, 2008)
              midvalue = 0.5 * (assimVar(j)%x%maximum+assimVar(j)%x%minimum)
              IF(assimVar(j)%gen%mean(i).gt.midvalue)THEN
                assimVar(j)%gen%sigma(i) = assimVar(j)%gen%restmeta *(assimVar(j)%x%maximum-assimVar(j)%gen%mean(i)) / (assimVar(j)%x%maximum-midvalue);
              ELSE
                assimVar(j)%gen%sigma(i) = assimVar(j)%gen%restmeta*(assimVar(j)%gen%mean(i)-assimVar(j)%x%minimum) / (midvalue - assimVar(j)%x%minimum);
              ENDIF
            CASE(3) ! semirestricted with max (Turner et al, 2008)
              assimVar(j)%gen%sigma(i) = AMAX1(0.0,assimVar(j)%x%maximum-assimVar(j)%gen%mean(i)) * assimVar(j)%gen%semimeta
            CASE(2) ! semirestricted with min (Turner et al, 2008)
              assimVar(j)%gen%sigma(i) = AMAX1(0.0,assimVar(j)%gen%mean(i)-assimVar(j)%x%minimum) * assimVar(j)%gen%semimeta
            CASE DEFAULT ! 0 or 1, ie. unrestricted (Turner et al, 2008)
          END SELECT
!         ENDIF
        ! just make sure we dont have sigma = 0
          assimVar(j)%gen%sigma(i) = AMAX1(assim_minsigma,assimVar(j)%gen%sigma(i))
        ENDIF
      ENDDO
          
      ! b) generate random values with the assigned sigma
      IF(assimVar(j)%gen%dorandxy)THEN
        !Spatially correlated ensemble matrix for a specific variable type
        !CALL get_spatially_correlated_random_data(assimVar(j))   !CP161205 changed to handle binfiles
        CALL get_spatially_correlated_random_data2(assimVar(j)%x%nvar,assimVar(j)%x%nens,assimVar(j)%gen,localx)
      ELSE
        ! Non-correlated observation data
        ! loop over records in the observation matrix (irrespective of variable type)
        DO i=1,assimVar(j)%x%nvar
          ! generate ensemble for each record independently, no correlation
          IF(assimVar(j)%gen%mean(i).GT.-9998.)THEN
            !CALL get_random_vector(assimVar(j)%x%nens,assimVar(j)%gen%mean(i),assimVar(j)%gen%sigma(i),assimVar(j)%x%x(i,:)) !CP161205 have to use local array since %x not always allocated (binfiles)
            CALL get_random_vector_gaussian(assimVar(j)%x%nens,assimVar(j)%gen%mean(i),assimVar(j)%gen%sigma(i),localx(i,:)) ! DG20170427 Changed to replace use of MKL dependent code.
          ELSE
            !assimVar(j)%x%x(i,:)=-9999.  !CP161205 same as above
            localx(i,:)=-9999.
          ENDIF
        ENDDO
      ENDIF
    ELSE
      ! in case someone runs with 1 ensemble member, we just use the mean, and skip all the ensemble generation
      !DO i=1,assimVar(j)%x%nvar
      !  assimVar(j)%x%x(i,1)=assimVar(j)%gen%mean(i)   !CP161205 use local x since x not allocated if binfile is used
      !ENDDO    
      DO i=1,assimVar(j)%x%nvar
        localx(i,1)=assimVar(j)%gen%mean(i)
      ENDDO    
    ENDIF
    
    ! save manipulated data from local array (to matrix or bin-file), also truncate all values outside the min/max range to the limits 
    !CALL assim_set_ensemble_data(assimVar(j)%x%nvar,assimVar(j)%x%nens,assimVar(j)%x,localx,Rec,.true.,.false.) !checkminmax true here
    CALL assim_set_ensemble_data(assimVar(j)%x%nvar,assimVar(j)%x%nens,assimVar(j)%x,localx,.true.,.false.) !checkminmax true here
    DEALLOCATE(localx)
    ! Finally, truncate all values outside the min/max range to the limits.
    !CALL assim_checkminmax(assimVar(j)%x%nvar,assimVar(j)%x%nens,assimVar(j)%x%x,assimVar(j)%x%minimum,assimVar(j)%x%maximum)  !CP161205 checked in assim_set_ensemble_data
  ENDDO
  
END SUBROUTINE generate_input_ensemble

!----------------------------------------------------------------------------
!> get_spatially_correlated_random_matrix_new
!>
!> Function for spatially correlated random variable generation
!>
!> 1) Random numbers are generated on a regular 2D grid using FFT
!>    (with variance=1, mean=0, and user defined correlation length) 
!> 
!> 2) The random numbers are interpolated to subbasins using bilinear interpolation 
!>
!> 3) Random numbers are scaled with the user defined standard deviation
!>
!> 4) The scaled random numbers are added to the mean observations 
!>
!> Author: David and Nils "Nicku" Gustafsson, 2013
!----------------------------------------------------------------------------
SUBROUTINE get_spatially_correlated_random_data(assimF) !n,nv,a,sigma,r)
  TYPE(assim_input_ensemble_type), INTENT(INOUT) :: assimF
  INTEGER :: n
  REAL(kind = 4), ALLOCATABLE :: r(:)
  INTEGER :: j, i, errcode,fieldsize,iens,jx
     
  ! 1) Random numbers are generated on a regular 2D grid using FFT
  !    (with variance=1, mean=0, and user defined correlation length)
  !
  ! 2) The random numbers are interpolated to subbasins using bilinear interpolation
  !
  ! 3) Random numbers are scaled with the user defined standard deviation
  !
  ! 4) The scaled random numbers are added to the mean observations 
  !
  ! some preparations
  n = assimF%x%nvar !(nsub typically)
  
  !ALLOCATE(r(n))  ! ALLOCATE r(n = nsub)
  
  ! loop over ensemble member iens in ensemble type j
   DO iens=1,assimF%x%nens
     
     !REAL FUNCTION
     !1) generate random numbers with standard deviation 1 and correlation length L
     CALL resample_randxy_data(assimF%gen%myrandxy_data)

     !2) scale perturbations with the required sigma (standard deviation) and appy to the mean
     DO i=1,n
       assimF%x%x(i,iens) =assimF%gen%myrandxy_data%pert_xy(i) * assimF%gen%sigma(i)  + assimF%gen%mean(i)
     ENDDO
     
     !!PROTOTYPE FUNCTION
     !! generate random numbers with standard deviation 1 and correlation length L
     !!
     !! this is still just a prototype where we just assume correlation 1, and thefore only need to draw 
     !! one normally distributed random number for each ensemble member
     !CALL get_random_number(0.,1.,r(1))
     !IF(fieldsize.GE.2)THEN
     !  DO i=2,fieldsize
     !    r(i) = r(1)
     !  ENDDO
     !ENDIF
     !! loop over records in the field and apply the random perturbation to the mean, scaled with it's individual sigma
     !DO i=1,fieldsize
     !  !IF(ensemble%mi(i).EQ.ensemble%types(j)) ensemble%y(i,:) = r(:) * ensemble%sigma(i)  + ensemble%mean(i)
     !  !jx = i+(j-1)*fieldsize
     !  ensemble%x%x(i,iens) = r(i) * ensemble%gen%sigma(i)  + ensemble%gen%mean(i)
     !ENDDO
   ENDDO
END SUBROUTINE get_spatially_correlated_random_data

!> Function for spatially correlated random variable generation
!CP: More information above about the routine
SUBROUTINE get_spatially_correlated_random_data2(n,nens,assimG,X)
  
  !Argument declarations
  INTEGER, INTENT(IN) :: n      !<number of element (nsub usually)
  INTEGER, INTENT(IN) :: nens   !<number of ensembles
  TYPE(assim_generation_type),INTENT(INOUT) :: assimG !<generation of input data
  REAL,INTENT(OUT) :: X(n,nens)
  
  !Local variables
  INTEGER :: i,iens   !loop-variables
     
  ! loop over ensemble member iens
   DO iens=1,nens
     
     !1) generate random numbers with standard deviation 1 and correlation length L
     CALL resample_randxy_data(assimG%myrandxy_data)

     !2) scale perturbations with the required sigma (standard deviation) and appy to the mean
     DO i=1,n
       X(i,iens) =assimG%myrandxy_data%pert_xy(i) * assimG%sigma(i)  + assimG%mean(i)
     ENDDO
   ENDDO
   
END SUBROUTINE get_spatially_correlated_random_data2

!! USE THE RGAUSS() FROM THE RANDOM_ROUTINES MODULE INSTEAD OF MKL library, DG 2017-04-27
SUBROUTINE get_random_vector_gaussian(n,a,sigma,r)
  REAL(KIND=4) r(:)
  REAL(KIND=4) a,sigma
  INTEGER n, i
  !loop over vector and generate randum numbers
  DO i=1,n
    ! random number N[0,1], scaled to the requested sigma, and applied to the mean
    r(i)=rgauss()*sigma+a
  ENDDO
END SUBROUTINE get_random_vector_gaussian

!!> (non-correlated) RANDOM NUMBER GENERATION (scalar, 1D vector, 2D matrix)
!!> Methods: vdrnggaussian from the MKL library
!!>
!!> See further at the end of this module for spatially correlated random numbers
!! -----------------------------------------------------------------------------------
!FUNCTION initialize_random_stream(seed) RESULT(errcode)
!  INTEGER seed
!  INTEGER(KIND=4) errcode
!  rnd_brng=VSL_BRNG_MT19937
!  rnd_method=VSL_METHOD_DGAUSSIAN_ICDF
!  errcode=vslnewstream(rnd_stream, rnd_brng,  seed)
!END FUNCTION initialize_random_stream
!
!FUNCTION deinitialize_random_stream() RESULT(errcode)
!  INTEGER errcode
!  errcode=vsldeletestream(rnd_stream)
!END FUNCTION deinitialize_random_stream
!
!SUBROUTINE get_random_matrix(n,nv,a,sigma,r)
!  REAL(KIND=4), INTENT(OUT) :: r(:,:)
!  REAL(KIND=4), INTENT(IN) :: a(:),sigma(:)
!  INTEGER, INTENT(IN) :: n,nv
!  INTEGER :: I, errcode
!  DO I=1,nv
!    errcode=vsrnggaussian(rnd_method, rnd_stream, n, r(I,:), a(I), sigma(I))
!  ENDDO
!END SUBROUTINE get_random_matrix
!
!SUBROUTINE get_random_vector(n,a,sigma,r)
!  REAL(KIND=4) r(:)
!  REAL(KIND=4) a,sigma
!  INTEGER n,errcode
!  errcode=vsrnggaussian( rnd_method, rnd_stream, n, r, a, sigma )
!END SUBROUTINE get_random_vector
!
!SUBROUTINE get_random_number(a,sigma,r)
!  REAL(KIND=4) r
!  REAL(KIND=4) a,sigma
!  INTEGER n,errcode
!  n=1
!  errcode=vsrnggaussian( rnd_method, rnd_stream, n, r, a, sigma )
!END SUBROUTINE get_random_number    

! -----------------------------------------------------------------------------------
!> MATRIX OPERATIONS
!> Method MATMUL (intrinsic function) /David 2017-04-27
! -----------------------------------------------------------------------------------
SUBROUTINE matrixmatrixmultiply(mat1,mat2,matout)
  REAL(KIND=4),INTENT(IN) ::  mat1(:,:),mat2(:,:)
  REAL(KIND=4),INTENT(OUT) :: matout(:,:)
  ! make sure matout is 0
  matout(:,:)=0.0
  ! CALL MATMUL
  matout = matmul(mat1,mat2)
END SUBROUTINE matrixmatrixmultiply

! Not used anymore: Method GEMM from the BLAS library
!SUBROUTINE matrixmatrixmultiply(mat1,mat2,matout)
!  REAL(KIND=4),INTENT(IN) ::  mat1(:,:),mat2(:,:)
!  REAL(KIND=4),INTENT(OUT) :: matout(:,:)
!  ! make sure matout is 0
!  matout(:,:)=0.0
!  ! CALL gemm from BLAS
!  CALL gemm(mat1,mat2,matout)
!END SUBROUTINE matrixmatrixmultiply

! -----------------------------------------------------------------------------------
!> CHOLESKY SOLUTION (old "MKL" dependant)
!> Method from the LAPACK library
! -----------------------------------------------------------------------------------
!SUBROUTINE choleskysolution(P,Y,M)
!  REAL(KIND=4) :: P(:,:),Y(:,:)
!  REAL(KIND=4) :: M(:,:)
!  CALL posv(P,Y)
!  !M=inv(P) * Y, dvs posv(P,Y)=inv(P)*Y 
!  M(:,:)=Y(:,:)
!END SUBROUTINE choleskysolution

! -----------------------------------------------------------------------------------
!> CHOLESKY SOLUTION (new)
!> Methods from Numerical Recipes in Fortran 90, 2nd edition, Press et al, 1996. 
!> Adopted by David 2017-04-27
! -----------------------------------------------------------------------------------

!> Cholesky decomposition (from Numerical Recipes, adopted by David)
SUBROUTINE choldc(a,p,n)
  !USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  REAL, DIMENSION(:,:), INTENT(INOUT) :: a
  REAL, DIMENSION(:), INTENT(OUT) :: p
  !Given an N × N positive-definite symmetric matrix a, this routine constructs its Cholesky
  !decomposition, A = L · LT . On input, only the upper triangle of a need be given; it is
  !not modified. The Cholesky factor L is returned in the lower triangle of a, except for its
  !diagonal elements, which are returned in p, a vector of length N.
  INTEGER :: i,n
  REAL :: summ
  !n=assert_eq(size(a,1),size(a,2),size(p),’choldc’)
  do i=1,n
    summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
    if (summ <= 0.0) write(6,*)'choldc failed'
    p(i)=sqrt(summ)
    a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
end do
END SUBROUTINE choldc

!> Cholesky solver (from Numerical Recipes, adopted by David)
SUBROUTINE cholsl(a,p,b,x,n)
  !USE nrtype; USE nrutil, ONLY : assert_eq
  IMPLICIT NONE
  REAL, DIMENSION(:,:), INTENT(IN) :: a
  REAL, DIMENSION(:), INTENT(IN) :: p,b
  REAL, DIMENSION(:), INTENT(INOUT) :: x
  !Solves the set of N linear equations A · x = b, where a is a positive-definite symmetric
  !matrix. a (N × N) and p (of length N) are input as the output of the routine choldc.
  !Only the lower triangle of a is accessed. b is the input right-hand-side vector, of length N.
  !The solution vector, also of length N, is returned in x. a and p are not modified and can be
  !left in place for successive calls with different right-hand sides b. b is not modified unless
  !you identify b and x in the calling sequence, which is allowed.
  INTEGER :: i,n
  !n=assert_eq((/size(a,1),size(a,2),size(p),size(b),size(x)/),’cholsl’)
  do i=1,n !Solve L · y = b, storing y in x.
    x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/p(i)
  enddo
  do i=n,1,-1 !Solve LT · x = y.
    x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/p(i)
  enddo
END SUBROUTINE cholsl

!> choleskysolutionNR - new routine called from enkf_analysis_prepare, written by David 2017-04-27
SUBROUTINE choleskysolution(ny,ne,P,Y,M)
  !arguments
  INTEGER, INTENT(IN)  :: ny,ne
  REAL,    INTENT(INOUT)  :: P(ny,ny)
  REAL,    INTENT(IN)     :: Y(ny,ne)
  REAL,    INTENT(OUT) :: M(ny,ne)
  !local variables
  REAL                 :: Ldiag(ny)
  integer              :: i 
  !This routine uses Cholesky factorization to solve the set of linear eqations P*M=Y => M=inv(P)*Y
  !The routines for cholesky factorization and solving above are taken from Numerical Recipes for Fortran 90.
  !The result M, is returned to enkf_analysis_prepare.
  !By David, 2017-04-27
  
  ! 1) Cholesky factorization, P = L · LT
  CALL choldc(P,Ldiag,ny)     !Cholesky factors are now in lower triangle of P, except the diagonal elements which are stored in Ldiag
  
  ! 2) Cholesky solution, solving P · M = Y => M = inv(P) · Y, using P and Ldiag from previous steps
  DO i=1,ne
    ! Since Y and M are multi-column ensemble matrices, we apply the solution for each ensemble member in a loop
    CALL  cholsl(P,Ldiag,Y(:,i),M(:,i),ny)
  ENDDO
  
END SUBROUTINE choleskysolution


! ---------------------------------------------------------------------------------------------
!> enkf_analysis_prepare
!> enkf_analysis_apply
!> enkf_analysis_main
!>
!> Ensemble Kalman filter equations (Evensen) following (a) Mandel, J. Efficient implementation 
!> of the ensemble kalman filter, and (b) DeChant, C. Examining the effectiveness and 
!> robustness of sequential data assimilation methods for quantification of uncertainty
!> in hydrological forecasting. Localization following Magnusson, Gustafsson et al (2014)
!>
!> Part I: Innovations (Y=D-HX), and inversion of variances M=1/(var(HX) + R) which is later
!>         used in the update equation: X_analysis = X_forecast + cov(X,HX)/(var(HX)+R) * (D-HX)
!
!Author: D.Gustafsson, SMHI/KTH.
!-----------------------------------------------------------------------------------------------
!The implementation is divided into two steps according to this pseudo-code (enkf_analysis_main)
!  CALL enkf_analysis_main()
!      
!   1. CALL enkf_analysis_prepare()
!        a. calculation of innovation, Y = D-HX
!        b. calculation and localization of the covariance of predicted observations Cyy = cov(HX) = locCyy * cov(HX)
!        c. combine cov(HX) with observation error variance R into P = (cov(HX)+R)
!        d. derive intermediate update matrix M, by inversion of variance matrice sum  M = P^(-1)Y
!
!   2. DO i=1,num_model_states
!        CALL enkf_analysis_apply()
!             Here in Part 2, M is multiplied with the covariance Cxy = Cov(X,HX) in the final EnKF update equation:
!      ENDDO
! 
! Terminology, EnKF basics from Mandel and DeChant
! ---------------------------------------------------------------------------------------
! Variable   Mandel        Dechant     Description
! --------   ------------  ----------  --------------------------------------------------
! N          N              Nens        ensemble size, n:o of ensemble members
! NX         n              -           state vector length, n:o state variables
! ND         m              
! X          X              X           model state variables
! A          A              e           ensemble deviation (from ensemble mean, A=X-E(X))
! D          D              y+eps       observations (including observation error eps)
! HX         HX             y'          predicted observations
! -          h(x)           h(x)        observation operator y'=h(x) translating model states to observation space
! HA         HA             -           HX ensemble deviation from ensemble mean (HA=HX-E(HX))
! Y          Y=D-HX         y+eps-y'    innovation, deviation between observations and predicted obserations
! -          C              P           model state error covariance, C = AA^T (Mandel) = ee^T = P (DeChant)
! R          R              R           observation error variance (theoretic or sample, assumed to be uncorrelated)
! -          HCH^T          C_yy        variance of predicted observations
! -          A(HA)^T        C_xy        covariance between model states and predicted observations
! PP         P              C_yy+R      sum of predicted observation variance and observation error variance
! M          M=(P^-1)Y      -           intermediate result of the inversion M = (P^-1)Y, with P and Y from Mandel
! locCyy     -              -           localization matrix for the variance of HX
! locCxy     -              -           localization matrix for the covariance between X and HX 
!-----------------------------------------------------------------------------------
SUBROUTINE enkf_analysis_prepare(N,ND,D,HX,R,locCyy,M,Y,HA)
  !INPUT ARGUMENTS
  INTEGER, INTENT(IN)  :: N            !n:o ensemble members
  INTEGER, INTENT(IN)  :: ND           !n:o observations
  REAL,    INTENT(IN)  :: D(ND,N)      !observation ensemble
  REAL,    INTENT(IN)  :: HX(ND,N)     !predicted observation ensemble
  REAL,    INTENT(IN)  :: R(ND)        !observation error variance, uncorrelated
  REAL,    INTENT(IN)  :: locCyy(ND,ND)!localization matrix for cov(HX)
  REAL,    INTENT(OUT) :: M(ND,N)      !inverse of (R+CXY) multiplied with innovation Y
  REAL,    INTENT(OUT) :: Y(ND,N)      !innovation ensemble
  REAL,    INTENT(OUT) :: HA(ND,N)     !HX deviation from mean(HX)
  !LOCAL VARIABLES
  REAL, ALLOCATABLE    :: e_Nx1(:,:)   !unit vector, size Nx1
  REAL, ALLOCATABLE    :: e_1xN(:,:)   !unit vector, size 1xN
  REAL, ALLOCATABLE    :: z_NDx1(:,:)  !intermediate result vector, size NDx1
  REAL, ALLOCATABLE    :: PP(:,:)      !intermediate result vector, size NDxND
  INTEGER              :: I,J          !loop index
  
  !ALLOCATION and INITIALIZATION
  ALLOCATE(e_Nx1(N,1))
  ALLOCATE(e_1xN(1,N))
  ALLOCATE(z_NDx1(ND,1))
  ALLOCATE(PP(ND,ND))
  
  !1) assign unit vectors
  e_Nx1 = 1.
  e_1xN = 1.
    
  !2) HA = HX-E(HX) = HX-1/N * (HX*e_Nx1)*e_1xN,  in three steps:
  ! a. z_NDx1 = HX*e_Nx1
  CALL matrixmatrixmultiply(HX,e_Nx1,z_NDx1)
  ! b. HA = z_NDx1 * e_1xN
  CALL matrixmatrixmultiply(z_NDx1, e_1xN, HA) !HA = z * e_1xN
  ! c. !HA = HX - 1/N * HA
  HA = HX - HA/N                               
    
  !3) Y = D-HX, in one step
  Y = D - HX

    
  !4) PP = (R+cov(HX)), observation error variance + predicted observation covariance, in three step:              
  ! a. PP = 1/(N-1) * HA * (HA)' = cov(predicted_observations)
  CALL matrixmatrixmultiply(HA, TRANSPOSE(HA), PP) ; PP = PP/(N-1)
    
  ! b. Localization by elementwise multiplication, PP = locCyy[nobs,nobs] .* PP
  DO I=1,ND
    DO J=1,ND
      PP(J,I)=PP(J,I)*locCyy(J,I)
    ENDDO
  ENDDO
    
  ! c. PP = PP + R, adding obs error variance, using the theoretical (assumed uncorrelated) variance R 
  !                 rather than the sample covariance
  DO I = 1,ND
    PP(I,I) = PP(I,I) + R(I)
  ENDDO

  ! 5) Cholesky solution... LL' = P, M=inv(P) * Y    ie. M = "innovations" / (cov(obs) + cov(model))
!  CALL choleskysolution(PP, Y, M) ! call to old MKL based function
  CALL choleskysolution(ND, N, PP, Y, M)
  
  !Deallocate temporary results
  DEALLOCATE(PP)
  DEALLOCATE(e_Nx1)
  DEALLOCATE(e_1xN)
  DEALLOCATE(z_NDx1)

END SUBROUTINE enkf_analysis_prepare

! -----------------------------------------------------------------------------------
!> enkf_analysis_apply
!>
!> EnKF analysis Part II: Apply ENKF filter on a model variable ensemble using the M
!>  matrix derived in the enkf_analysis_prepare function
!>
!> EnKF analysis, following the "basic implementation" in Jan Mandel report modified
!>                to include covariance localization: 
!>
!>                          (K = Cxy .* loc_Cxy / (Cyy .* loc_Cyy + R))
!>
!>                The localization of Cyy is embedded in the M matrix which is 
!>                prepared by the enkf_analysis_prepare routine.
!>
!>                The localization of Cxy is introduced in the final step of this
!>                "apply"-routine, which involves a change in the calculation of the AZ
!>                matrix compared to Mandel:
!>
!>                Mandel:    AZ = A * (HA' * M)
!>
!>                Here:      AZ = ((A*HA') .* loc_Cxy) * M
!>
!>                Matrix operations are made in an order defined by the parantheses (most inner operations first).
!>                This might have consequences on the number of operations compared to Mandel, but I have not checked that.
!>
!> More info, see enkf_analysis_prepare
! -----------------------------------------------------------------------------------
SUBROUTINE enkf_analysis_apply(N,NX,ND,Y,M,HA,locCxy,X)
  !INPUT VARIABLES
  INTEGER, INTENT(IN)   :: N  !NE
  INTEGER, INTENT(IN)   :: NX
  INTEGER, INTENT(IN)   :: ND
  REAL,INTENT(INOUT)    :: X(NX,N)
!  REAL,INTENT(IN)       :: M(NX,N)
  REAL,INTENT(IN)       :: M(ND,N)
  REAL,INTENT(IN)       :: Y(ND,N)
  REAL,INTENT(IN)       :: HA(ND,N)
  REAL,INTENT(IN)       :: locCxy(NX,ND)
        
  !LOCAL VARIABLES
  REAL, ALLOCATABLE :: e_Nx1(:,:)
  REAL, ALLOCATABLE :: e_1xN(:,:)
  REAL, ALLOCATABLE :: z_NXx1(:,:)
  REAL, ALLOCATABLE :: Z_NXxND(:,:)
  REAL, ALLOCATABLE :: A(:,:)
  REAL, ALLOCATABLE :: AZ(:,:)
  INTEGER           :: I, J
  
  !ALLOCATION and INITIALIZATION
  ALLOCATE(A(NX,N))
  ALLOCATE(AZ(NX,N))
  ALLOCATE(e_Nx1(N,1))
  ALLOCATE(e_1xN(1,N))
  ALLOCATE(z_NXx1(NX,1))
  ALLOCATE(z_NXxND(NX,ND))
    
  !0) assign unit vectors
  e_Nx1 = 1.
  e_1xN = 1.

  !1) A = X-E(X) = X-1/N*(X*e_Nx1)*e_1xN
  CALL matrixmatrixmultiply(X, e_Nx1, z_NXx1)
  CALL matrixmatrixmultiply(z_NXx1, e_1xN, A)
  A = X - A / N

  !2) EnKF filtering equation
  !
  ! Modification of Mandel's two-step algoritm, to include localization of 
  ! covariances between innovations (Y=HX-D) and all model states (X) 
  !
  ! The Kalman filter analysis is written as:
  !     X = X + locCxy .* (A*HA') * M / (N-1) ;
  !  
  ! Which is formulate the two-step solution from Mandel, into a three-step solution:
  !  7.1: Z_NXxND = A*HA'
  CALL matrixmatrixmultiply(A, TRANSPOSE(HA), Z_NXxND)
  !    
  !  7.2: Z_NXxND = Z_NXxND .* locCX (element-wise)
  DO I=1,ND
    DO J=1,NX
      Z_NXxND(J,I)=Z_NXxND(J,I) * locCxy(J,I)
    ENDDO
  ENDDO
  !    
  !  7.3: Z_NXxN  = Z_NXxND * M      (we can use AZ for this)
  CALL matrixmatrixmultiply(Z_NXxND, M, AZ)
  !
  !! 7.1 Z = (HA)' * M
  !CALL matrixmatrixmultiply(TRANSPOSE(HA), M, Z_NxN)
  ! 
  !! 7.2) X=X+1/(N-1) * A * Z
  !CALL matrixmatrixmultiply(A, Z_NxN, AZ)
  !    
  ! 8) final update, just as before:
  X = X + AZ / (N-1)

END SUBROUTINE enkf_analysis_apply

!------------------------------------------------------
!> Ensemble statistics
!>  Calculates summary statistics for an ensembe matrix
!------------------------------------------------------
SUBROUTINE assim_ensemble_statistics(xin, DIM, NN, xmean, xmins, xmaxs, xsigma, xmedian, domedian) !xquantiles,xcovar, xcor)
  USE compout, only: calculate_median

! Compute summary statistics for a matrix
  INTEGER DIM, NN, M
  REAL xin(DIM,NN)
  REAL xmean(DIM), xmins(DIM), xmaxs(DIM), xsigma(DIM), xmedian(DIM) !xquantiles(3,DIM), xcovar(DIM,DIM), xcor(DIM,DIM)
!  REAL o_stat(NN,DIM),xr2(DIM)
!  REAL o_quant(3)
  LOGICAL domedian
  
! Data needed for simplified statistical calculation
  REAL s(NN), pp(NN)
  REAL var, ep
  INTEGER i
  
! Data needed to use the VSL Summary statistics routine
!  TYPE(VSL_SS_TASK) task
!  INTEGER p
!  INTEGER n
!  INTEGER x_storage
!  INTEGER o_storage
!  !      INTEGER cov_storage
!  !      INTEGER cor_storage
!  INTEGER i, errcode
!  INTEGER(KIND=8) estimate
!  INTEGER task_method

  
  ! loop over dimension DIM
  DO i=1,DIM
    ! mean value
    xmean(i)=sum(xin(i,:))/NN
    
    ! Min and Max
    xmins(i)=minval(xin(i,:))
    xmaxs(i)=maxval(xin(i,:))
    
    ! standard deviation
    s(:)=xin(i,:)-xmean(i)
    ep=sum(s(:))
    pp(:)=s(:)*s(:)
    var = (sum(pp(:))-ep**2/NN)/(NN-1)
    xsigma(i)=sqrt(var)
    
    ! median
    IF(domedian)Call calculate_median(NN,xin(i,:),-9999.,xmedian(i))
      
    ! skip the quantiles
  ENDDO
    
!! Seems like this VSL code is quite slow, or at least when using it many times, it uses a lot of resources.
!! we try something simpler above
!  
!  ! ***** Initializing parameters for Summary Statistics task *****
!  p               = DIM
!  n               = NN
!  x_storage       = VSL_SS_MATRIX_STORAGE_COLS
!  !   cov_storage     = VSL_SS_MATRIX_STORAGE_FULL
!  !     cor_storage     = VSL_SS_MATRIX_STORAGE_FULL
!  task_method     = VSL_SS_METHOD_FAST
!  o_storage   = VSL_SS_MATRIX_STORAGE_ROWS
!      
!  o_quant(1) = 0.025 ; o_quant(2) = 0.5 ; o_quant(3) = 0.975
!      
!  M=3
!
!  !     ***** Create Summary Statistics task *****
!  errcode = vslsssnewtask( task, p, n, x_storage, xin )
!
!  !     ***** Set initial values of the min/max estimates *****
!  DO i = 1, p
!    xmins(i) = xin(i, 1)
!    xmaxs(i) = xin(i, 1)
!  ENDDO
!
!  !     ***** Edit task parameters for min and max computation *****
!  errcode = vslsssedittask( task, VSL_SS_ED_MIN, xmins )
!  errcode = vslsssedittask( task, VSL_SS_ED_MAX, xmaxs )
!
!  !      ***** Minimum and maximum are included in the list of estimates to compute *****
!  estimate = IOR( VSL_SS_MIN, VSL_SS_MAX )
!
!  !     ***** Compute the estimates using FAST method *****
!  errcode = vslssscompute( task, estimate, task_method )
!
!  !     ***** Delete Summary Statistics task *****
!  errcode = vslssdeletetask( task )
!
!  !     ***** Create Summary Statistics task *****
!  errcode = vslsssnewtask( task, p, n, x_storage, xin )
!
!  !     ***** Edit task parameters for computating of mean estimate and 2nd, 3rd
!  !           and 4th raw and central moments estimates *****
!  errcode = vslssseditmoments( task, xmean, xr2, c2m = xsigma)
!
!  !     ***** Mean and 2nd, 3rd and 4th raw and central moments are included
!  !            in the list of estimates to compute *****
!  estimate = IOR( VSL_SS_MEAN,   IOR(VSL_SS_2R_MOM,VSL_SS_2C_MOM))
!
!  !      estimate = IOR( estimate, IOR( VSL_SS_MEAN,                              &
!  !     &           IOR( VSL_SS_2R_MOM, IOR( VSL_SS_3R_MOM,                       &
!  !     &           IOR( VSL_SS_4R_MOM, IOR( VSL_SS_2C_MOM,                       &
!  !     &           IOR( VSL_SS_3C_MOM, VSL_SS_4C_MOM ) ) ) ) ) ) )
!
!  !     ***** Compute the estimates using FAST method *****
!  errcode = vslssscompute( task, estimate, task_method )
!
!  !     ***** Delete Summary Statistics task *****
!  errcode = vslssdeletetask( task )
!
!!!     ***** Create Summary Statistics task *****
!!      errcode = vslsssnewtask( task, p, n, x_storage, xin )
!!
!!!     ***** Initialization of the task parameters using FULL_STORAGE
!!!           for covariance/correlation matrices computation *****
!!      errcode = vslssseditcovcor( task, xmean, xcovar, cov_storage, xcor, cor_storage )
!!
!!!     ***** Mean, Covariance and correlation matrices are included in the list
!!!           of estimates to compute *****
!!      estimate = IOR( VSL_SS_COV, VSL_SS_COR )
!!
!!!     ***** Compute the estimates using FAST method *****
!!      errcode = vslssscompute( task, estimate, task_method )
!!     
!!!     ***** Delete Summary Statistics task *****
!!      errcode = vslssdeletetask( task )
!
!
!  !     ***** Create Summary Statistics task *****
!  errcode = vslsssnewtask( task, p, n, x_storage, xin )
!  errcode = vslssseditquantiles( task, M, o_quant, xquantiles, o_stat, o_storage )
!  estimate = IOR( VSL_SS_QUANTS, VSL_SS_ORDER_STATS )
!    
!  !     ***** Compute the estimates using FAST method *****
!  errcode = vslssscompute( task, estimate, task_method )
!     
!  !     ***** Delete Summary Statistics task *****
!  errcode = vslssdeletetask( task )
  
END SUBROUTINE assim_ensemble_statistics

!---------------------------------------------------------------------
!> assim_get_ensemble_data()
!>
!> Routine that returns the ensemble data in a matrix. 
!> If needed, the data is read from binary file.
!> For speed, records in binary files contain all states in one
!> ensemble member (X transformated)
!---------------------------------------------------------------------
SUBROUTINE assim_get_ensemble_data(NX,NE,ensData,x)
  !input
  INTEGER, INTENT(IN)                    :: NX, NE    !number of states and ensemble members
  REAL, INTENT(OUT)                      :: X(NX,NE)  !ensemble data matrix
  TYPE(assim_ensemble_type),INTENT(IN)   :: ensData   !ensemble data structure
  
  !Local variables
  INTEGER j
  
  !ensemble data ALLOCATED in memory
  IF(ALLOCATED(ensData%x))THEN   !change because needed more options for different kind of bin-files
    !X=ensData%x(:,Rec:(Rec+NE-1))  !CP161202 change Rec to 0 (from 1 in general)
    X=ensData%x(:,ensData%rec+1:ensData%rec+NE)
  ELSE
  !ensemble data read from direct access binary file
    DO j = 1,NE
      READ(ensdata%fileID,REC=ensData%rec+j) X(:,j)
    ENDDO
    !X=0.   !CP161201 added read for bin-file instead of X=0 here
  ENDIF
END SUBROUTINE assim_get_ensemble_data

!---------------------------------------------------------------------
!> assim_set_ensemble_data()
!>
!> Routine that writes a matrix into an ensemble data structure. 
!> If needed, the data is written to a binary file.
!---------------------------------------------------------------------
SUBROUTINE assim_set_ensemble_data(NX,NE,ensData,x,checkMinMax,doStat)
  !input
  INTEGER, INTENT(IN)                     :: NX, NE
  REAL, INTENT(INOUT)                     :: X(NX,NE)
  TYPE(assim_ensemble_type),INTENT(INOUT) :: ensData
  LOGICAL, INTENT(IN)                     :: checkMinMax !check min and max limits before update
  LOGICAL, INTENT(IN)                     :: doStat      !calculate output statistics

  !Local variables
  INTEGER j
  
  !Check min and max limits
  IF(checkMinMax)THEN
    CALL assim_checkminmax(nx,ne,x,ensdata%minimum,ensData%maximum)
  ENDIF

  !Write x to ensemble data...
  IF(ALLOCATED(ensData%x))THEN
    !... ALLOCATED in memory I don't under stand how Rec works here, removed.
    !ensData%x(:,Rec:(Rec+NE-1))=X  !CP161202 change Rec to 0 (from 1 in general), Rec is always 0 here
    ensData%x(:,ensData%rec+1:ensData%rec+NE)=X
  ELSE
    !... or in direct access binary file.
    DO j = 1,NE
      WRITE(ensdata%fileID,REC=ensData%rec+j) X(:,j)
    ENDDO
  ENDIF
  
  !Update ensemble output statistics
  !IF(doStat)THEN
  !  CALL assim_ensemble_statistics(x, NX, NE, ensData%outmean, ensData%outquant, ensData%outmin, ensData%outmax, ensData%outsigma)
  !ENDIF
END SUBROUTINE assim_set_ensemble_data

SUBROUTINE updateEnsembleStatistics(assimData,total_time)
  TYPE(assim_data_type) assimData
  REAL total_time(4)
  REAL start_time, stop_time
  INTEGER i
  REAL,ALLOCATABLE :: X(:,:)  !intermediate matrix of ensemble
  LOGICAL domedian

  call cpu_time(start_time)

  ! set domedian if needed for outputs
  IF(assimData%info%meanout)THEN
    domedian=.FALSE. ! median NOT needed for outputs
  ELSE
    domedian=.TRUE.  ! median IS needed for outputs
  ENDIF
  
  !state ensembles
  DO i=1,assimData%info%nX
    !Allocate and collect X from where it is saved before calculating the statistics CP added to handle bin-files.
    ALLOCATE(X(assimData%X(i)%x%nvar,assimData%X(i)%x%nens))

    call cpu_time(stop_time)
    total_time(1)=total_time(1)+stop_time-start_time
    start_time=stop_time

    CALL assim_get_ensemble_data(assimData%X(i)%x%nvar,assimData%X(i)%x%nens,assimData%X(i)%x,X)   !CP161201 added call to routine to get X-matrix (useful if bin-file)

    call cpu_time(stop_time)
    total_time(2)=total_time(2)+stop_time-start_time
    start_time=stop_time

!    CALL assim_ensemble_statistics(X, assimData%X(i)%x%nvar, assimData%X(i)%x%nens, assimData%X(i)%x%outmean, assimData%X(i)%x%outquant, assimData%X(i)%x%outmin, assimData%X(i)%x%outmax, assimData%X(i)%x%outsigma)
    CALL assim_ensemble_statistics(X, assimData%X(i)%x%nvar, assimData%X(i)%x%nens, assimData%X(i)%x%outmean, assimData%X(i)%x%outmin, assimData%X(i)%x%outmax, assimData%X(i)%x%outsigma, assimData%X(i)%x%outquant(2,:),domedian)

    call cpu_time(stop_time)
    total_time(3)=total_time(3)+stop_time-start_time
    start_time=stop_time

    DEALLOCATE(X)

    call cpu_time(stop_time)
    total_time(4)=total_time(4)+stop_time-start_time
    start_time=stop_time

    !CALL assim_ensemble_statistics(assimData%X(i)%x%x, assimData%X(i)%x%nvar, assimData%X(i)%x%nens, assimData%X(i)%x%outmean, assimData%X(i)%x%outquant, assimData%X(i)%x%outmin, assimData%X(i)%x%outmax, assimData%X(i)%x%outsigma)
  ENDDO

  !forcing ensembles
  DO i=1,assimData%info%nF
    !Allocate and collect X from where it is saved before calculating the statistics, added to handle bin-files.
    ALLOCATE(X(assimData%F(i)%x%nvar,assimData%F(i)%x%nens))
    CALL assim_get_ensemble_data(assimData%F(i)%x%nvar,assimData%F(i)%x%nens,assimData%F(i)%x,X)
!    CALL assim_ensemble_statistics(X, assimData%F(i)%x%nvar, assimData%F(i)%x%nens, assimData%F(i)%x%outmean, assimData%F(i)%x%outquant, assimData%F(i)%x%outmin, assimData%F(i)%x%outmax, assimData%F(i)%x%outsigma)
    CALL assim_ensemble_statistics(X, assimData%F(i)%x%nvar, assimData%F(i)%x%nens, assimData%F(i)%x%outmean, assimData%F(i)%x%outmin, assimData%F(i)%x%outmax, assimData%F(i)%x%outsigma, assimData%F(i)%x%outquant(2,:),domedian)
    DEALLOCATE(X)
    !CALL assim_ensemble_statistics(assimData%F(i)%x%x, assimData%F(i)%x%nvar, assimData%F(i)%x%nens, assimData%F(i)%x%outmean, assimData%F(i)%x%outquant, assimData%F(i)%x%outmin, assimData%F(i)%x%outmax, assimData%F(i)%x%outsigma) !CP161205 use binfiles!
  ENDDO

  !auxiliary ensembles
  
  DO i=1,assimData%info%nA
    !Allocate and collect X from where it is saved before calculating the statistics, added to handle bin-files.
    ALLOCATE(X(assimData%A(i)%x%nvar,assimData%A(i)%x%nens))
    CALL assim_get_ensemble_data(assimData%A(i)%x%nvar,assimData%A(i)%x%nens,assimData%A(i)%x,X)
!    CALL assim_ensemble_statistics(X, assimData%A(i)%x%nvar, assimData%A(i)%x%nens, assimData%A(i)%x%outmean, assimData%A(i)%x%outquant, assimData%A(i)%x%outmin, assimData%A(i)%x%outmax, assimData%A(i)%x%outsigma)
    CALL assim_ensemble_statistics(X, assimData%A(i)%x%nvar, assimData%A(i)%x%nens, assimData%A(i)%x%outmean, assimData%A(i)%x%outmin, assimData%A(i)%x%outmax, assimData%A(i)%x%outsigma, assimData%A(i)%x%outquant(2,:),domedian)
    DEALLOCATE(X)
    !CALL assim_ensemble_statistics(assimData%A(i)%x%x, assimData%A(i)%x%nvar, assimData%A(i)%x%nens, assimData%A(i)%x%outmean, assimData%A(i)%x%outquant, assimData%A(i)%x%outmin, assimData%A(i)%x%outmax, assimData%A(i)%x%outsigma)
  ENDDO

END SUBROUTINE updateEnsembleStatistics

!----------------------------------------------------------------------
!> enkf_analysis_main
!>
!> Main routine in the enkf analysis.
!
! Organizes the calls to enkf_analysis_prepare and enkf_analysis_apply
! following this pseudo-code:
!   1. CALL enkf_analysis_prepare()
!        a. calculation of innovation, Y = D-HX
!        b. calculation and localization of the covariance of predicted observations cov(HX) = locCyy * cov(HX)
!        c. combine cov(HX) with observation error variance R into P = (cov(HX)+R)
!        d. derive intermediate update matrix M, by inversion of variance matrice sum  M = P^(-1)Y
!
!   2. DO i=1,num_model_states
!        CALL enkf_analysis_apply()
!             Here in Part 2, M is multiplied with the covariance Cov(X,HX) in the final EnKF update equation:
!      ENDDO
!----------------------------------------------------------------------------
SUBROUTINE enkf_analysis_main(assimData)
  !ARGUMENT
  TYPE(assim_data_type), INTENT(INOUT) :: assimData
  !LOCAL VARIABLES
  INTEGER :: N, ND, NVAR, I, NX, J, NVARold, locID
  REAL, ALLOCATABLE :: X(:,:), Y(:,:), HA(:,:), M(:,:), D(:,:), HX(:,:), R(:) !, locCyy(:,:), locCxy(:,:)
  
  !Continue with ENKF analysis only if there is some data to assimilate
  IF(assimData%info%nD==0) RETURN   !No ENKF analysis if no data
  !IF(assimData%info%nD.GT.0)THEN   !CP161205 changed to return for nD=0 to have one less IF-ENDIF to keep track of
  !Part I:  calculate innovations and covariance matrix inversion
  !--------------------------------------------------------------
    !assign local variables
    ND = assimData%info%nD    !N:o observations to assimlate
    N  = assimData%info%NE    !N:o ensemble members
       
    !ALLOCATE necessary matrices
    ALLOCATE(Y(ND,N))
    ALLOCATE(M(ND,N))
    ALLOCATE(HA(ND,N))
    ALLOCATE(HX(ND,N))
    ALLOCATE(D(ND,N))
    ALLOCATE(R(ND))
    !ALLOCATE(locCyy(ND,ND))
      
    !get D, HX, and locCyy data from Enkf data structure (read from binary files if needed)
    !CALL assim_get_ensemble_data(ND,N,assimData%D,D,1)
    !CALL assim_get_ensemble_data(ND,N,assimData%HX,HX,1)
    !CALL assim_get_ensemble_data(ND,1,assimData%R,R,1)
    !CALL assim_get_ensemble_data(ND,ND,assimData%LocCYY,locCyy,1)
    D=assimData%D
    HX=assimData%HX
    R=assimData%R
  
    !locCyy =   assimData%LocCYY
    
    !CALL preparation routine
    CALL enkf_analysis_prepare(N,ND,D,HX,R,assimData%LocCYY,M,Y,HA)
  
    !DEALLOCATE temporary variables
    DEALLOCATE(D) ; DEALLOCATE(HX) ; DEALLOCATE(R) ; !DEALLOCATE(locCyy)
      
    !Part II: enkf update on model variables one-by-one
    !--------------------------------------------------
    !X state ensembles
    NX=assimData%Info%nX              !number of variable ensembles (types)
    IF(NX.GT.0)THEN
      NVAR=assimData%X(1)%x%nvar        !number of variables in ensemble (rows)
      !locID = assimData%X(i)%x%locID !localization ID      
      NVARold=0
      !Loop over number of variable ensembles, if >0
      DO I=1,NX
        !check if this variable is analysed or re-initialized
        IF(assimData%X(I)%x%assimilate)THEN
          !re-ALLOCATE X and locCxy matrices, if needed
          NVAR=assimData%X(I)%x%nvar        !number of variables in ensemble (rows)
          IF(ALLOCATED(X).AND.NVAR.ne.NVARold)DEALLOCATE(X)
          IF(.not.ALLOCATED(X))ALLOCATE(X(NVAR,N))
          !IF(ALLOCATED(locCxy).AND.NVAR.ne.NVARold)DEALLOCATE(locCxy)
          !IF(.not.ALLOCATED(locCxy))ALLOCATE(locCxy(NVAR,ND))
          NVARold = NVAR
   
          !Get ensemble data into X matrix
          CALL assim_get_ensemble_data(NVAR,N,assimData%X(I)%x,X)   !CP161202 for no bin-file or several bin-files

          !get localization matrix locCXY
          locID = assimData%X(i)%x%locID !localization ID      
          !CALL assim_get_ensemble_data(NVAR,ND,assimData%locCXY(locID),locCxy,1)
          !locCxy = assimData%locCXY

!          CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY,X)  !CP161214 added aquifer extra coordinate system
          CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY(locID)%x,X)
            
          !save updated X matrix back to the ensemble data
          CALL assim_set_ensemble_data(NVAR,N,assimData%X(I)%x,X,.true.,.true.)
        ENDIF
      ENDDO
    ENDIF
    
    !A auxilary ensembles (outvar in HYPE)
    NX=assimData%Info%nA              !number of variable ensembles (types)
    IF(NX.GT.0)THEN
      NVAR=assimData%A(1)%x%nvar        !number of variables in ensemble (rows)
      !locID = assimData%A(i)%x%locID !localization ID
      NVARold=0
        
      !Loop over number of variable ensembles, if >0
      DO I=1,NX
        !check if this variable is analysed or re-initialized
        IF(assimData%A(I)%x%assimilate)THEN
          !re-ALLOCATE X and locCxy matrices, if needed
          NVAR=assimData%A(I)%x%nvar        !number of variables in ensemble (rows)
          IF(ALLOCATED(X).AND.NVAR.ne.NVARold)DEALLOCATE(X)
          IF(.not.ALLOCATED(X))ALLOCATE(X(NVAR,N))
          !IF(ALLOCATED(locCxy).AND.NVAR.ne.NVARold)DEALLOCATE(locCxy)
          !IF(.not.ALLOCATED(locCxy))ALLOCATE(locCxy(NVAR,ND))
          NVARold = NVAR
   
          !Get ensemble data into X matrix
          CALL assim_get_ensemble_data(NVAR,N,assimData%A(I)%x,X)   !CP161206 for no bin-file or several bin-files

          !get localization matrix locCXY
          locID = assimData%X(i)%x%locID !localization ID      
          !CALL assim_get_ensemble_data(NVAR,ND,assimData%locCXY(locID),locCxy,1)
          !locCxy = assimData%locCXY

          !CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY,X)  !CP161214 added aquifer extra coordinate system
          CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY(locID)%x,X)
            
          !save updated X matrix back to the ensemble data
          CALL assim_set_ensemble_data(NVAR,N,assimData%A(I)%x,X,.true.,.true.)

        ENDIF
      ENDDO
    ENDIF

    !F forcing ensembles (P, T, SW etc in HYPE)
    NX=assimData%Info%nF              !number of variable ensembles (types)
    IF(NX.GT.0)THEN
      NVAR=assimData%F(1)%x%nvar        !number of variables in ensemble (rows)
      !locID = assimData%F(i)%x%locID !localization ID
      NVARold=0
        
      !Loop over number of variable ensembles, if >0
      DO I=1,NX
        !check if this variable is analysed or re-initialized
        IF(assimData%F(I)%x%assimilate)THEN
          !re-ALLOCATE X and locCxy matrices, if needed
          NVAR=assimData%F(I)%x%nvar        !number of variables in ensemble (rows)
          IF(ALLOCATED(X).AND.NVAR.ne.NVARold)DEALLOCATE(X)
          IF(.not.ALLOCATED(X))ALLOCATE(X(NVAR,N))
          !IF(ALLOCATED(locCxy).AND.NVAR.ne.NVARold)DEALLOCATE(locCxy)
          !IF(.not.ALLOCATED(locCxy))ALLOCATE(locCxy(NVAR,ND))
          NVARold = NVAR
   
          !Get ensemble data into X matrix
          CALL assim_get_ensemble_data(NVAR,N,assimData%F(I)%x,X)

          !get localization matrix locCXY
          locID = assimData%X(i)%x%locID !localization ID      
          !CALL assim_get_ensemble_data(NVAR,ND,assimData%locCXY(locID),locCxy,1)
          !locCxy = assimData%locCXY

          !CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY,X)  !CP161214 added aquifer extra coordinate system
          CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY(locID)%x,X)
            
          !save updated X matrix back to the ensemble data
          CALL assim_set_ensemble_data(NVAR,N,assimData%F(I)%x,X,.true.,.true.)
        ENDIF
      ENDDO
    ENDIF
    !To come : parameter ensembles
    !Deallocate
    IF(ALLOCATED(X))     DEALLOCATE(X)
    IF(ALLOCATED(Y))     DEALLOCATE(Y)
    IF(ALLOCATED(M))     DEALLOCATE(M)
    !IF(ALLOCATED(HX))    DEALLOCATE(HX)  !CP161206 already deallocated above
    !IF(ALLOCATED(D))     DEALLOCATE(D)
    !IF(ALLOCATED(R))     DEALLOCATE(R)
    IF(ALLOCATED(HA))    DEALLOCATE(HA)
    !IF(ALLOCATED(locCyy))DEALLOCATE(locCyy)
    !IF(ALLOCATED(locCxy))DEALLOCATE(locCxy)
  !ENDIF
END SUBROUTINE enkf_analysis_main
  
!-------------------------------------------------------------------------
!> enkf_ensemble_output
!> PURPOSE: print output from the Enkf ensembles:
!  To be Modified to the new data structures, Feb-2015 /David
!-------------------------------------------------------------------------
!subroutine enkf_ensemble_output()
!  ! mean, median, 95% quantiles, min, max, stdev in 7 text files:
!  !
!  ! write to these files... fid_enkf_min, fid_enkf_max, fid_enkf_median, fid_enkf_final, fid_enkf_mean, fid_enkf_min95, fid_enkf_max95,fid_enkf_std
!  !
!  ! LOCAL VARIABLES:
!  INTEGER :: errcode,i,j,k,strlen
!  CHARACTER(LEN=30000) :: outputstringMEDIAN,outputstringMIN95,outputstringMAX95
!  CHARACTER(LEN=30000) :: outputstringSTDEV
!  CHARACTER(LEN=20)    :: fmtstring
!
!  ! ADD DATA TO OUTPUTSTRING:
!  outputstringMEDIAN(:) =' '
!  outputstringMIN95(:)  =' '
!  outputstringMAX95(:)  =' '
!  outputstringSTDEV(:)  =' '
!    
!  ! 1) model state variables from the X ensemble
!  IF(EnkfInfo%XS)THEN
!    DO i=1,EnkFInfo%nX
!      WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfX%outquant(2,i)
!      WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfX%outquant(1,i)
!      WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfX%outquant(3,i)
!      !WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfX%outcovar(i,i))**0.5
!      WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfX%outsigma(i))
!      k=k+10
!    ENDDO
!  ENDIF
!  ! 2) model auxiliary variables from A ensemble
!  DO i=1,EnkFInfo%nA
!    WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfA%outquant(2,i)
!    WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfA%outquant(1,i)
!    WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfA%outquant(3,i)
!    !WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfA%outcovar(i,i))**0.5
!    WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfA%outsigma(i))
!    k=k+10
!  ENDDO
!    
!  ! 3) model parameters from P ensemble
!  DO i=1,EnkFInfo%nP
!    WRITE(outputstringMEDIAN(k:k+13),'(f14.6)') EnkfP%outquant(2,i)
!    WRITE(outputstringMIN95(k:k+13), '(f14.6)') EnkfP%outquant(1,i)
!    WRITE(outputstringMAX95(k:k+13), '(f14.6)') EnkfP%outquant(3,i)
!    !WRITE(outputstringSTDEV(k:k+13), '(f14.6)') (EnkfP%outcovar(i,i))**0.5
!    WRITE(outputstringSTDEV(k:k+13), '(f14.6)') (EnkfP%outsigma(i))
!    k=k+14
!  ENDDO
!   
!  ! 4) forcing from F ensemble
!!    DO i=1,EnkFInfo%nF
!!        WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfF%outquant(2,i)
!!        WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfF%outquant(1,i)
!!        WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfF%outquant(3,i)
!!!        WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfF%outcovar(i,i))**0.5
!!        WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfF%outsigma(i))
!!        k=k+10
!!    ENDDO
!
!  ! 5) "Innovations" (aka mod-obs) from D ensemble (check the enkf analysis scheme...)
!  !j=0
!  !DO i=1,EnkFInfo%nD
!  !    if ()THEN
!  !        j=j+1
!  !    WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfD%outquant(2,i)
!  !    WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfD%outquant(1,i)
!  !    WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfD%outquant(3,i)
!  !    WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfD%outcovar(i,i))**0.5
!  !    k=k+10
!  !    ELSE
!  !        WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') -9999.
!  !        WRITE(outputstringMIN95(k:k+9), '(f10.4)') -9999.
!  !        WRITE(outputstringMAX95(k:k+9), '(f10.4)') -9999.
!  !        WRITE(outputstringSTDEV(k:k+9), '(f10.4)') -9999.
!  !    ENDIF
!  !ENDDO
!
!
!  ! PREPARE OUTPUT FORMAT
!  outputstringMEDIAN = TRIM(adjustr(outputstringMEDIAN(1:k)))
!  outputstringMIN95 = TRIM(adjustr(outputstringMIN95(1:k)))
!  outputstringMAX95 = TRIM(adjustr(outputstringMAX95(1:k)))
!  outputstringSTDEV = TRIM(adjustr(outputstringSTDEV(1:k)))
!
!  fmtstring = enkf_get_fmtstring(outputstringMEDIAN)
!  WRITE(fid_enkf_median,TRIM(fmtstring))TRIM(outputstringMEDIAN)
!
!  fmtstring = enkf_get_fmtstring(outputstringMIN95)
!  WRITE(fid_enkf_min95,TRIM(fmtstring))TRIM(outputstringMIN95)
!
!  fmtstring = enkf_get_fmtstring(outputstringMAX95)
!  WRITE(fid_enkf_max95,TRIM(fmtstring))TRIM(outputstringMAX95)
!
!  fmtstring = enkf_get_fmtstring(outputstringSTDEV)
!  WRITE(fid_enkf_std,TRIM(fmtstring))TRIM(outputstringSTDEV)
!end subroutine enkf_ensemble_output

!!-------------------------------------------------------------------------
!!> enkf_open_timeseriesoutput
!!-------------------------------------------------------------------------
!subroutine enkf_open_timeseriesoutput(resultdir)
!  CHARACTER(LEN=200) resultdir    
!  open(fid_enkf_median,file=TRIM(resultdir)//'enkf_output_median.txt',STATUS = 'unknown')
!  open(fid_enkf_min95,file=TRIM(resultdir)//'enkf_output_min95.txt',STATUS = 'unknown')
!  open(fid_enkf_max95,file=TRIM(resultdir)//'enkf_output_max95.txt',STATUS = 'unknown')
!  open(fid_enkf_std,file=TRIM(resultdir)//'enkf_output_std.txt',STATUS = 'unknown')
!  WRITE(fid_enkf_median,*)'% ENKF MEDIAN OUTPUT:'
!  IF(EnkfInfo%XS)THEN
!    WRITE(fid_enkf_median,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_min95,*)'% ENKF 2.5 quantile OUTPUT:'
!    WRITE(fid_enkf_min95,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_max95,*)'% ENKF 97.5 quantile OUTPUT:'
!    WRITE(fid_enkf_max95,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_std,*)'% ENKF standard deviation OUTPUT:'
!    WRITE(fid_enkf_std,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!  ELSE
!    WRITE(fid_enkf_median,*)'% Timestep A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_min95,*)'% ENKF 2.5 quantile OUTPUT:'
!    WRITE(fid_enkf_min95,*)'% Timestep A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_max95,*)'% ENKF 97.5 quantile OUTPUT:'
!    WRITE(fid_enkf_max95,*)'% Timestep A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_std,*)'% ENKF standard deviation OUTPUT:'
!    WRITE(fid_enkf_std,*)'% Timestep A(1:na) P(1:np) F(1:nf)'    
!  ENDIF
!end subroutine enkf_open_timeseriesoutput
!
!subroutine enkf_close_timeseriesoutput()
!  close(fid_enkf_median)
!  close(fid_enkf_min95)
!  close(fid_enkf_max95)
!  close(fid_enkf_std)
!end subroutine enkf_close_timeseriesoutput

FUNCTION assim_get_fmtstring(outputstri) RESULT(fmtstri)
  CHARACTER(LEN=10000) outputstri
  CHARACTER(LEN=20)    :: fmtstri
  INTEGER strlen

  fmtstri(1:3)=' (A'
  strlen = LEN(TRIM(outputstri))

  IF(strlen.lt.10)THEN
    WRITE(fmtstri(4:4),'(i1)')strlen
    fmtstri(5:6)=')'
  ELSEIF(strlen.lt.100)THEN
    WRITE(fmtstri(4:5),'(i2)')strlen
    fmtstri(6:7)=')'
  ELSEIF(strlen.lt.1000)THEN
    WRITE(fmtstri(4:6),'(i3)')strlen
    fmtstri(7:8)=')'
  ELSEIF(strlen.lt.10000)THEN
    WRITE(fmtstri(4:7),'(i4)')strlen
    fmtstri(8:9)=')'
  ELSE
    WRITE(fmtstri(4:8),'(i5)')strlen
    fmtstri(9:10)=')'
  ENDIF
END FUNCTION assim_get_fmtstring

! -----------------------------------------------------------------------------------
! VARIOUS "OLD" ENSEMBLE DATA MANIPULATION ROUTINES
! -----------------------------------------------------------------------------------

!!>Initialize the TYPE(enkf_info) enkfinfo variable
!!>
!!---------------------------------------------------------
!
!SUBROUTINE initialize_enkfinfo()
!
!  INTEGER :: errcode
!  
!  ! Random number generation
!  !-----------------------------------------
!  !EnkfInfo%seed = 777
!  CALL time(EnkfInfo%seed) ! the random seed will be different, maybe make it possible for user to decide via the enkfinfo.txt?
!  
!  errcode = initialize_random_stream(EnkfInfo%seed)
!  
!  EnkfInfo%nE = 100       ! number of ensembles (columns in ensemble matrix X)
!  EnkfInfo%nX = 0         ! number of state variables (rows) in ensemble matrix X, excluding augmented variables
!  EnkfInfo%nA = 0         ! number of additional model variables in A ensemble (used for output)
!  EnkfInfo%nF = 0         ! number of input in F ensemble
!  EnkfInfo%nD = 0         ! number of observations in D ensemble
!  EnkfInfo%nDA = 0         ! number of observations in D ensemble at current Assimilation
!  
!  EnkfInfo%nP = 0         ! number of parameters in P ensemble
!  EnkfInfo%nAX = 0        ! number of A in X matrix
!  EnkfInfo%nFX = 0        !  -"- of F in X
!  EnkfInfo%nPX = 0        !  -"- of P in X
!  
!  ! LOGICAL flags, may be modified...
!  EnkfInfo%FQ = .false.          ! include fluxes in kalman filter     (general switch on/off)
!  EnkfInfo%FA = .false.          ! include auxil. in kalman filter     (general switch on/off)
!  EnkfInfo%FP = .false.          ! include parameters in kalman filter (general switch on/off) 
!  EnkfInfo%FF = .false.          ! include inputs in kalman filter     (general switch on/off)
!  EnkfInfo%XS = .false.          ! include X states in statistical output (general switch on/off)
!  EnkfInfo%EC = .false.          ! ECMWF forecasts (general switch on/off)
!  EnkfInfo%meanout = .true.      ! ensemble mean(.true.) or median (.false.) in output
!  ! some parameters
!  EnkfInfo%moradkhani_delta = 0.95 ! coef.[0-1] to retain variance in parameter ensemble (Moradkhani et al, 2004), value typical around 0.95 
!  
!  EnkfInfo%enkf_minsigma = 1.e-5
!  
!  EnkfInfo%useBinFiles = .FALSE.
!  
!end subroutine initialize_enkfinfo

!!>Allocate all variables needed for the EnKF data assimiltion method
!!>
!!-------------------------------------------------------------------
!
!subroutine allocate_enkf_data(ne,nx,na,np,nf,nd,nax,npx,nfx,ndt)
!  INTEGER :: ne,nx,na,np,nf,nd,nax,npx,nfx,ndt
!  CHARACTER(LEN=1) :: etype
!  
!  ! ENSEMBLES: EnkfX, EnkfA, EnkfP, EnkfF, EnkfD, EnkfHX
!  !-----------------------------------------------------
!  etype = 'X'
!  CALL allocate_ensemble(EnkfX,ne,nx+nax+npx+nfx,etype)
!
!  etype = 'A'
!  CALL allocate_ensemble(EnkfA,ne,na,etype)
!
!  etype = 'P'
!  CALL allocate_ensemble(EnkfP,ne,np,etype)
!
!  etype = 'F'
!  CALL allocate_ensemble(EnkfF,ne,nf,etype)
!
!  etype = 'D'
!  CALL allocate_ensemble(EnkfD,ne,nd,etype)
!
!  etype = 'D'
!  CALL allocate_ensemble(EnkfHX,ne,nd,etype)
!  
!  etype = 'B'
!  CALL allocate_ensemble(EnkfDinfo,ne,ndt,etype)
!end subroutine allocate_enkf_data
!
!!>Allocate an TYPE(enkf_ensemble_data) variable
!!>
!!---------------------------------------------------------
!subroutine allocate_ensemble(ensemble,ne,n,etype)
!  TYPE(enkf_ensemble_data) :: ensemble
!  INTEGER :: i,j,ne,n
!  character :: etype
!
!  !ensemble size and type
!  ensemble%n     = n       !number of variables
!  ensemble%ne    = ne      !number of ensemble members
!  ensemble%etype = etype   !ensemble type (perturbed and/or analysed in Kalman filter
!  ensemble%ECMWF = .false. !special type used for ECMWF forecasts (for a particular study)
!  
!  IF(n.gt.0)THEN
!    ! ALLOCATE AND ASSIGNS DEFAULT VALUES DEPENDING ON TYPE
!    
!    !ensemble matrix y - only if binary files are not used  
!    IF(etype /= 'B')THEN           !for all ensembles except obs-info
!      ALLOCATE(ensemble%y(n,ne))   ! [n x ne] matrix
!      ensemble%y(:,:) = 0.0
!    ENDIF
!    
!    !for all, except state ensemble X
!    IF(etype /= 'X')THEN
!      ALLOCATE(ensemble%mi(n)) !id number in relation to the model data, defined by the user in "interface"
!      ensemble%mi(:) = 0
!      
!      ALLOCATE(ensemble%xi(n)) !id number in ensemble matrix X (not needed in the X ensemble)
!      ensemble%xi(:)= 0
!    ENDIF
!
!    !only for obs-Info
!    IF(etype == 'B')THEN
!      ALLOCATE(ensemble%id(n))
!    ENDIF
!    
!    !for forcing(F), observation(D), parameters(P), and obs-info (B) 
!    IF(etype == 'F' .or. etype == 'D' .or. etype == 'P' .or. etype == 'B') then
!      ALLOCATE(ensemble%ensgen(n))        ! type of ensemble (0 none, 1 unrestricted, 2 semi-min, 3 semi-max, 4 restricted (min-max))
!      ensemble%ensgen(:) = 1
!    ENDIF   
!  
!    IF(etype /= 'X' .and. etype /= 'D' .and. etype /= 'B')THEN
!      ALLOCATE(ensemble%kalman(n))        ! 0 no, 1 yes (switch Kalman filter on/off)
!      ensemble%kalman(:) = 0
!    ENDIF
!
!    !IF(etype /= 'D' .and. etype /= 'F')THEN
!    IF(etype /= 'X')THEN
!      ALLOCATE(ensemble%mean(n))            ! mean value (used for initialization)
!      ensemble%mean(n) = 0.0
!
!      ALLOCATE(ensemble%sigma(n))           ! stdev (used for generation)
!      ensemble%sigma(:) = 0.0
!
!      ALLOCATE(ensemble%semimeta(n))           ! stdev (used for generation)
!      ensemble%semimeta(:) = 0.0
!
!      ALLOCATE(ensemble%restmeta(n))           ! stdev (used for generation)
!      ensemble%restmeta(:) = 0.0
! 
!    ! covariance can be really big, better keep it ALLOCATED until it is needed...which is at assimilation analysis
!      !ALLOCATE(ensemble%covar(n,n))         ! covariance (used for generation)
!      !ensemble%covar(:,:)= 0.0
!    ENDIF
!
!    ALLOCATE(ensemble%minimum(n))       ! min allowed (for check)
!    ensemble%minimum(:) = enkf_defmin        ! min allowed (for check)
!  
!    ALLOCATE(ensemble%maximum(n))       ! max allowed (for check)
!    ensemble%maximum(:) = enkf_defmax
!
!    IF(etype /= 'X' .or. EnkfInfo%XS)THEN
!      !ALLOCATE(ensemble%statout(n))       ! yes(1) no(0) - print out statistics for this variable
!      !ensemble%statout(:) = 1
!      !ALLOCATE(ensemble%ensout(n))        ! yes(1) no(0) - print out entire ensemble for this variable
!      !ensemble%ensout(:) = 0
!      
!      ! some variables for output statistics:... can we avoid the outcovar and outcorr as well??
!      ALLOCATE(ensemble%outmean(n), ensemble%outquant(3,n), ensemble%outmin(n), ensemble%outmax(n), ensemble%outsigma(n)) !, ensemble%outcovar(n,n), ensemble%outcorr(n,n))
!      
!      ensemble%outmean(:)=0.0
!      ensemble%outquant(:,:)=0.
!      ensemble%outmin(:)=0.
!      ensemble%outmax(:)=0.
!      ensemble%outsigma(:)=0.
!    ENDIF
!  ENDIF
!end subroutine allocate_ensemble

!!>Truncate an TYPE(enkf_ensemble_data) ensemble variable to minimum and maximum allowed values.
!!>
!!----------------------------------------------------------------------------------------------
!subroutine enkf_checkminmax(ensemble)
!  TYPE(enkf_ensemble_data) :: ensemble
!  INTEGER :: i,j
!  DO i=1,ensemble%n
!    DO j=1,ensemble%ne
!      ensemble%y(i,j)=AMAX1(ensemble%minimum(i),AMIN1(ensemble%maximum(i),ensemble%y(i,j)))
!    ENDDO
!  ENDDO
!end subroutine enkf_checkminmax

!!>General routine for Ensemble generation (forcing and observation data).
!!>
!!> The ensemble generation is made by adding random numbers to the input data.
!!>
!!> The basic assumption is that the random perturbations are gaussian with zero mean and standard deviation sigma.
!!>
!!> However, based on Turner et al(2008), it is assumed that in order to get unbiased input ensembles, 
!!> we need to consider two types of perturbations - systematic and random:
!!> Thus, the input x at time k for ensemble member j, x_kj = x_k0+eata_kj+chi_j,
!!> where eata_kj is regenerated every time step and chi_j is generated only once.
!!> Then, eata and chi can be generated with three different types of restrictions on variables.
!!> For the moment, we assume that the static error chi_j = 0.
!!>
!!> In adition, we now also take into account spatial correlation in the data, by generating 
!!> spatially correlated random data using a FFT-based method.
!!----------------------------------------------------------------------------------------------
!subroutine generate_ensemble(ensemble,vectorwise)
!  TYPE(enkf_ensemble_data) :: ensemble
!  INTEGER i
!  REAL midvalue
!  LOGICAL vectorwise
!  
!  ! first of all, update the standard deviations depending on the generation method    
!  IF(ensemble%ne.gt.1)THEN
!    ! a)determine standard deviation, sigma for each variable in each subbasin
!    !
!    ! loop over variables (nvar x nsub)
!    DO i=1,ensemble%n
!      ! select General case or special ECMWF case
!      IF(.not. ensemble%ECMWF)THEN
!        SELECT CASE(ensemble%ensgen(i))
!          CASE(4) ! restricted (Turner et al, 2008)
!            midvalue = 0.5 * (ensemble%maximum(i)+ensemble%minimum(i))
!            IF(ensemble%mean(i).gt.midvalue)THEN
!              ensemble%sigma(i) = ensemble%restmeta(i) *(ensemble%maximum(i)-ensemble%mean(i)) / (ensemble%maximum(i)-midvalue);
!            ELSE
!              ensemble%sigma(i) = ensemble%restmeta(i) *(ensemble%mean(i)-ensemble%minimum(i)) / (midvalue - ensemble%minimum(i));
!            ENDIF
!          CASE(3) ! semirestricted with max (Turner et al, 2008)
!            ensemble%sigma(i) = AMAX1(0.0,ensemble%maximum(i)-ensemble%mean(i)) * ensemble%semimeta(i)
!          CASE(2) ! semirestricted with min (Turner et al, 2008)
!            ensemble%sigma(i) = AMAX1(0.0,ensemble%mean(i)-ensemble%minimum(i)) * ensemble%semimeta(i)
!          CASE DEFAULT ! 0 or 1, ie. unrestricted (Turner et al, 2008)
!        END SELECT
!      ENDIF
!      ! just make sure we dont have sigma = 0
!      ensemble%sigma(i) = AMAX1(EnkfInfo%enkf_minsigma,ensemble%sigma(i))
!    ENDDO
!      
!    ! b) generate random values with the assigned sigma
!    IF(vectorwise)THEN
!      ! this old functions only work for non-correlated observation data
!      ! loop over records in the observation matrix (irrespective of variable type)
!      DO i=1,ensemble%n
!        ! generate ensemble for each record independently, no correlation
!        CALL get_random_vector(ensemble%ne,ensemble%mean(i),ensemble%sigma(i),ensemble%y(i,:))
!      ENDDO
!    ELSE
!      ! OLD function: get non-correlated ensemble matrix (essentially the same as above)
!      !CALL get_random_matrix(ensemble%ne,ensemble%n,ensemble%mean,ensemble%sigma,ensemble%y)
!      ! NEW function: get spatially correlated ensemble matrix for a specific variable type
!      CALL get_spatially_correlated_random_matrix(ensemble) 
!    ENDIF
!  ELSE
!    ! in case someone runs with 1 ensemble member, we just use the mean, and skip all the enkf filtering
!    DO i=1,ensemble%n
!      ensemble%y(i,1)=ensemble%mean(i)
!    ENDDO    
!  ENDIF
!  ! Finally, truncate all values outside the min/max range to the limits.
!  CALL enkf_checkminmax(ensemble)
!end subroutine generate_ensemble

!!> (non-correlated) RANDOM NUMBER GENERATION (scalar, 1D vector, 2D matrix)
!!> Methods: vdrnggaussian from the MKL library
!!>
!!> See further at the end of this module for spatially correlated random numbers
!! -----------------------------------------------------------------------------------
!FUNCTION initialize_random_stream(seed) RESULT(errcode)
!  INTEGER seed
!  INTEGER(KIND=4) errcode
!  rnd_brng=VSL_BRNG_MT19937
!  rnd_method=VSL_METHOD_DGAUSSIAN_ICDF
!  errcode=vslnewstream(rnd_stream, rnd_brng,  seed)
!END FUNCTION initialize_random_stream
!
!FUNCTION deinitialize_random_stream() RESULT(errcode)
!  INTEGER errcode
!  errcode=vsldeletestream(rnd_stream)
!END FUNCTION deinitialize_random_stream
!
!subroutine get_random_matrix(n,nv,a,sigma,r)
!  REAL(KIND=4), INTENT(OUT) :: r(:,:)
!  REAL(KIND=4), INTENT(IN) :: a(:),sigma(:)
!  INTEGER, INTENT(IN) :: n,nv
!  INTEGER :: I, errcode
!  DO I=1,nv
!    errcode=vsrnggaussian(rnd_method, rnd_stream, n, r(I,:), a(I), sigma(I))
!  ENDDO
!end subroutine get_random_matrix
!
!subroutine get_random_vector(n,a,sigma,r)
!  REAL(KIND=4) r(:)
!  REAL(KIND=4) a,sigma
!  INTEGER n,errcode
!  errcode=vsrnggaussian( rnd_method, rnd_stream, n, r, a, sigma )
!end subroutine get_random_vector
!
!subroutine get_random_number(a,sigma,r)
!  REAL(KIND=4) r
!  REAL(KIND=4) a,sigma
!  INTEGER n,errcode
!  n=1
!  errcode=vsrnggaussian( rnd_method, rnd_stream, n, r, a, sigma )
!end subroutine get_random_number    
!
!! -----------------------------------------------------------------------------------
!!> MATRIX OPERATIONS
!!> Method GEMM from the BLAS library
!! -----------------------------------------------------------------------------------
!subroutine matrixmatrixmultiply(mat1,mat2,matout)
!  REAL(KIND=4),INTENT(IN) ::  mat1(:,:),mat2(:,:)
!  REAL(KIND=4),INTENT(OUT) :: matout(:,:)
!  ! make sure matout is 0
!  matout(:,:)=0.0
!  ! CALL gemm from BLAS
!  CALL gemm(mat1,mat2,matout)
!end subroutine matrixmatrixmultiply
!
!! -----------------------------------------------------------------------------------
!!> CHOLESKY FACTORIZATION
!!> Method from the LAPACK library
!! -----------------------------------------------------------------------------------
!subroutine choleskysolution(P,Y,M)
!  REAL(KIND=4) :: P(:,:),Y(:,:)
!  REAL(KIND=4) :: M(:,:)
!  CALL posv(P,Y)
!  M(:,:)=Y(:,:)
!end subroutine choleskysolution
  
!! -----------------------------------------------------------------------------------
!!> EnKF analysis
!!> Implementation following J. Mandel "Basic implementation...", extended with localization.
!!
!! This FUNCTION only have to care about the following matrices:
!!   X     ModelEnsemble            [n_modvar,n_ens]
!!   D     ObservationEnsemble      [n_obs,n_ens]     
!!   R     ObsCovarMatrix           [n_obs,n_obs]
!!   HX    ModelObservationEnsemble [n_ons,n_ens]
!! -----------------------------------------------------------------------------------
!subroutine enkf_analysis_routine(N, NX, ND) !Xd, HXd, Dd, Rd, N, NX, ND)
!  ! input arguments
!  INTEGER N, NX, ND   ! n:o ensembles N, n:o states NX, n:o observations ND 
!  ! local variables
!  REAL, ALLOCATABLE :: HA(:,:)                  ! HA(ND,N)                    
!  REAL, ALLOCATABLE :: M(:,:)                   ! M(ND,N)                     
!  REAL, ALLOCATABLE :: e_Nx1(:,:), e_1xN(:,:)   ! e_Nx1(N,1),   e_1xN(1,N)    
!  REAL, ALLOCATABLE :: z_NXx1(:,:), z_NDx1(:,:) ! z_NXx1(NX,1), z_NDx1(ND,1)  
!  REAL, ALLOCATABLE :: A(:,:)                   ! A(NX,N)
!  REAL, ALLOCATABLE :: PP(:,:)                  ! PP(ND,ND)
!  REAL, ALLOCATABLE :: Z_NxN(:,:)               ! Z_NxN(N,N)
!  REAL, ALLOCATABLE :: AZ(:,:)                  ! AZ(NX,N)
!  INTEGER I,J,K
!  ! localization matrices
!  REAL, ALLOCATABLE :: Z_NXxND(:,:)
!  ! ALLOCATE some matrices
!  ALLOCATE(HA(ND,N))
!  ALLOCATE(M(ND,N))
!  ALLOCATE(e_Nx1(N,1), e_1xN(1,N))
!  ALLOCATE(z_NXx1(NX,1), z_NDx1(ND,1))
!  ALLOCATE(A(NX,N))
!  ALLOCATE(PP(ND,ND))
!  ALLOCATE(Z_NxN(N,N))
!  ALLOCATE(AZ(NX,N))
!  ALLOCATE(Z_NXxND(NX,ND))
!  !
!  !initialize the unit vectors
!  e_Nx1 = 1.0
!  e_1xN = 1.0
!  !
!  ! EnKF analysis, following the "basic implementation" in Jan Mandel report
!  !
!  ! 1) A = X-E(X) = X-1/N*(X*e_Nx1)*e_1xN
!  CALL matrixmatrixmultiply(EnkfX%y, e_Nx1, z_NXx1)
!  CALL matrixmatrixmultiply(z_NXx1, e_1xN, A)
!  A = EnkfX%y - A / N
!  !
!  ! 2) HX = observationfunction(X) - already provided in input
!  !
!  ! 3) HA, in two steps:
!  ! 3.1) z = (HX)*e_Nx1
!  CALL matrixmatrixmultiply(EnkfHX%y(1:enkfInfo%nDA,:),e_Nx1,z_NDx1)
!  ! 3.2) HA = HX-1/N * z * e_1*N                [nobs x N]
!  CALL matrixmatrixmultiply(z_NDx1, e_1xN, HA)
!  HA = EnkfHX%y(1:enkfInfo%nDA,:) - HA / N
!  !
!  ! 4) Y = D - HX [n_obs x N]   ie. obs-model "the innovations"
!  !
!  EnkfD%y(1:enkfInfo%nDA,:) = EnkfD%y(1:enkfInfo%nDA,:) - EnkfHX%y(1:enkfInfo%nDA,:)
!  !               
!  ! 5) P = R + 1/(N-1) * HA * (HA)'    ie.  cov(obs) + cov(model)
!  !
!  CALL matrixmatrixmultiply(HA, TRANSPOSE(HA), PP)
!  PP = PP/(N-1)
!  !
!  ! Here, localization is added, PP = locC_YY[nobs,nobs] .* HA * (HA)' /(N-1)
!  ! (cancel out covariance beyond certain distance)
!  DO I=1,ND
!    DO J=1,ND
!      PP(J,I)=PP(J,I)*EnkfLoc%locCYY(J,I)
!    ENDDO
!  ENDDO
!  !
!  ! R is assumed uncorrelated, and estimated from input rather than from the sample covariance
!  !
!  DO I = 1,ND
!    PP(I,I) = PP(I,I) + EnkfD%sigma(I)**2
!  ENDDO
!  !PP = EnkfD%covar(1:enkfInfo%nDA,1:enkfInfo%nDA) + PP/(N-1)
!  !
!  ! 6) Cholesky solution... LL' = P, M=inv(P) * Y    ie. "innovations" / (cov(obs) + cov(model))
!  CALL choleskysolution(PP, EnkfD%y(1:enkfInfo%nDA,:), M)
!  !
!  ! 7) EnKF filtering equation, in two steps
!  !
!  ! localization is needed also here, to cancel out any covariance beyond a certain distance
!  ! between the innovations and all other model states 
!  !
!  ! if we skip the Z step below, the Kalman filter step can be written as:
!  !     X = X + locC_XY[nx,nobs] .* (A*HA') * M / (N-1) ;
!  !  
!  ! Thus, we can reformulate the two-step solution from Mandel, into a three-step solution:
!  !  7.1: Z_NXxND = A*HA'
!  CALL matrixmatrixmultiply(A, TRANSPOSE(HA), Z_NXxND)
!  !    
!  !  7.2: Z_NXxND = Z_NXxND .* locCX (element-wise)
!  DO I=1,ND
!    DO J=1,NX
!      Z_NXxND(J,I)=Z_NXxND(J,I) * EnkfLoc%locCXY(J,I)
!    ENDDO
!  ENDDO
!  !    
!  !  7.3: Z_NXxN  = Z_NXxND * M      (we can use AZ for this)
!  CALL matrixmatrixmultiply(Z_NXxND, M, AZ)
!  !
!  !! 7.1 Z = (HA)' * M
!  !CALL matrixmatrixmultiply(TRANSPOSE(HA), M, Z_NxN)
!  ! 
!  !! 7.2) X=X+1/(N-1) * A * Z
!  !CALL matrixmatrixmultiply(A, Z_NxN, AZ)
!  !    
!  ! 8) final update, just as before:
!  EnkfX%y = EnkfX%y + AZ / (N-1)
!  !
!  ! Clean up, deallocation
!  DEALLOCATE(HA)
!  DEALLOCATE(M)
!  DEALLOCATE(e_Nx1, e_1xN)
!  DEALLOCATE(z_NXx1, z_NDx1)
!  DEALLOCATE(A)
!  DEALLOCATE(PP)
!  DEALLOCATE(Z_NxN)
!  DEALLOCATE(AZ)
!  DEALLOCATE(Z_NXxND)
!end subroutine enkf_analysis_routine
!  
!***********************************************************************************************
!> deallocate_ensemble_data
!> purpose: DEALLOCATE all ALLOCATED field in a TYPE(enkf_ensemble_data) variable
SUBROUTINE deallocate_ensemble(ensemble)
  TYPE(assim_ensemble_type) :: ensemble
  ! Deallocate all ALLOCATABLE fields
  IF(ALLOCATED(ensemble%x))DEALLOCATE(ensemble%x)
  IF(ALLOCATED(ensemble%outmean))DEALLOCATE(ensemble%outmean)       
  IF(ALLOCATED(ensemble%outquant))DEALLOCATE(ensemble%outquant)       
  IF(ALLOCATED(ensemble%outmin))DEALLOCATE(ensemble%outmin)       
  IF(ALLOCATED(ensemble%outmax))DEALLOCATE(ensemble%outmax)       
  IF(ALLOCATED(ensemble%outsigma))DEALLOCATE(ensemble%outsigma)       
END SUBROUTINE deallocate_ensemble
  !***********************************************************************************************
  !> deallocate_interface_data
  !> purpose: DEALLOCATE all ALLOCATED field in a TYPE(enkf_interface_data) variable
  SUBROUTINE deallocate_interface(ensemble)
    TYPE(assim_interface_type) :: ensemble
    ! Deallocate all ALLOCATABLE fields
    IF(ALLOCATED(ensemble%subDimID))DEALLOCATE(ensemble%subDimID)       
  END SUBROUTINE deallocate_interface
  !***********************************************************************************************
  !> deallocate_generation_data
  !> purpose: DEALLOCATE all ALLOCATED field in a TYPE(enkf_generation_data) variable
  SUBROUTINE deallocate_generation(ensemble)
    TYPE(assim_generation_type) :: ensemble
    ! Deallocate all ALLOCATABLE fields
    IF(ALLOCATED(ensemble%sigma))DEALLOCATE(ensemble%sigma)
    IF(ALLOCATED(ensemble%mean))DEALLOCATE(ensemble%mean)
    CALL deallocate_randxy_data(ensemble%myrandxy_data)
  END SUBROUTINE deallocate_generation

  !----------------------------------------------------------------------------------------------- 
  !> enkf_deallocation_all
  !>   purpose: deallocates and free all memory used by the enkf routines
  !----------------------------------------------------------------------------------------------- 
  SUBROUTINE assim_deallocation_all(assimData)
    !
    ! Purpose of this routine is to DEALLOCATE and free all memory
    ! occupied by the enkf routines
    !
    !ARGUMENTS
    TYPE(assim_data_type), INTENT(INOUT) :: assimData
    !LOCAL VARIABLES
    INTEGER :: i,errcode

    ! deallocation of all the ensemble data
  
    !state  variables
    IF(ALLOCATED(assimData%X))THEN
      DO i=1,assimData%Info%nX
        CALL deallocate_ensemble(assimData%X(i)%x)
        CALL deallocate_interface(assimData%X(i)%info)
      ENDDO
      DEALLOCATE(assimData%X)
    ENDIF
    !auxiliary variables
    IF(ALLOCATED(assimData%A))THEN
      DO i=1,assimData%Info%nA
        CALL deallocate_ensemble(assimData%A(i)%x)
        CALL deallocate_interface(assimData%A(i)%info)
      ENDDO
      DEALLOCATE(assimData%A)
    ENDIF
    !Parameters
  !  IF(ALLOCATED(assimData%P))THEN
  !    DO i=1,assimData%Info%nP
  !      CALL deallocate_ensemble(assimData%P(i)%x)
  !      CALL deallocate_interface(assimData%P(i)%info)
  !      CALL deallocate_generation(assimData%P(i)%gen)
  !    ENDDO
  !    DEALLOCATE(assimData%P)
  !  ENDIF
    !Forcing
    IF(ALLOCATED(assimData%F))THEN
      DO i=1,assimData%Info%nF
        CALL deallocate_ensemble(assimData%F(i)%x)
        CALL deallocate_interface(assimData%F(i)%info)
        CALL deallocate_generation(assimData%F(i)%gen)
      ENDDO
      DEALLOCATE(assimData%F)
    ENDIF
    !D and HX 
  !  CALL deallocate_ensemble(assimData%D)
  !  CALL deallocate_ensemble(assimData%HX)
  !  CALL deallocate_ensemble(assimData%R)
    DEALLOCATE(assimData%D)
    DEALLOCATE(assimData%HX)
    DEALLOCATE(assimData%R)
  
    !other data
    !IF(ALLOCATED(assimData%Coordinates))THEN
    !  DO i=1,assimData%Info%ncoord
    !    IF(ALLOCATED(assimData%Coordinates(i)%x))DEALLOCATE(assimData%Coordinates(i)%x)
    !    IF(ALLOCATED(assimData%Coordinates(i)%y))DEALLOCATE(assimData%Coordinates(i)%y)
    !    IF(ALLOCATED(assimData%Coordinates(i)%z))DEALLOCATE(assimData%Coordinates(i)%z)
    !  ENDDO
    !  DEALLOCATE(assimData%Coordinates)
    !ENDIF
    !Call deallocate_ensemble(assimData%LocCYY)
    !IF(ALLOCATED(assimData%LocCXY))THEN
    !  DO i=1,assimData%Info%nloc    
    !    Call deallocate_ensemble(assimData%LocCXY(i))
    !  ENDDO
    !ENDIF 
  
    IF(ALLOCATED(assimData%LocCYY))DEALLOCATE(assimData%LocCYY)
    IF(ALLOCATED(assimData%LocCXY))DEALLOCATE(assimData%LocCXY)
  
    ! shutdown random number stream
!    errcode = deinitialize_random_stream()
    
    ! free MKL memory usage
  !  CALL MKL_FREE_BUFFERS()

    ! close the time series outputs
    !close(fid_enkf_min)
    !close(fid_enkf_max)
    CLOSE(fid_assim_median)
    !close(fid_enkf_mean)
    CLOSE(fid_assim_min95)
    CLOSE(fid_assim_max95)
    !close(fid_enkf_std)
  END SUBROUTINE assim_deallocation_all

!!----------------------------------------------------------------------------------------------- 
!!> enkf_deallocation
!!>   purpose: deallocates and free all memory used by the enkf routines, and write final state
!!>            ensemble to file
!!----------------------------------------------------------------------------------------------- 
!subroutine enkf_deallocation(fname)
!  CHARACTER(LEN=300) :: fname
!  !
!  ! Purpose of this routine is to DEALLOCATE and free all memory
!  ! occupied by the enkf routines
!  !
!  INTEGER :: errcode,i,j,k,strlen
!  CHARACTER(LEN=10000) :: outputstring
!  CHARACTER(LEN=20)    :: fmtstring
!
!  ! write final ensemble state to file
!  open(fid_enkf_final,file=fname,STATUS = 'unknown')
!  WRITE(fid_enkf_final,'(A58)')'% FINAL ENSEMBLE (ensembles in columns, variables in rows....)'
!  j=1
!  outputstring(j:j)='%'
!  DO i=1,EnkfInfo%nX
!    outputstring(j+1:j+2)=' X'
!    IF(i.lt.10)THEN
!      WRITE(outputstring(j+3:j+3),'(i1)')i
!      j=j+3
!    ELSEIF(i.lt.100)THEN
!      WRITE(outputstring(j+3:j+4),'(i2)')i
!      j=j+4
!    ELSEIF(i.lt.1000)THEN
!      WRITE(outputstring(j+3:j+5),'(i3)')i
!      j=j+5
!    ELSEIF(i.lt.10000)THEN
!      WRITE(outputstring(j+3:j+6),'(i4)')i
!      j=j+6
!    ELSE
!      WRITE(outputstring(j+3:j+6),'(i5)')i
!      j=j+6
!    ENDIF
!  ENDDO
!  DO i=1,EnkfInfo%nP
!    outputstring(j+1:j+2)=' P'
!    IF(i.lt.10)THEN
!      WRITE(outputstring(j+3:j+3),'(i1)')i
!      j=j+3
!    ELSEIF(i.lt.100)THEN
!      WRITE(outputstring(j+3:j+4),'(i2)')i
!      j=j+4
!    ELSEIF(i.lt.1000)THEN
!      WRITE(outputstring(j+3:j+5),'(i3)')i
!      j=j+5
!    ELSEIF(i.lt.10000)THEN
!      WRITE(outputstring(j+3:j+6),'(i4)')i
!      j=j+6
!    ELSE
!      WRITE(outputstring(j+3:j+6),'(i5)')i
!      j=j+6
!    ENDIF
!  ENDDO
!  outputstring = TRIM(adjustr(outputstring(1:j)))
!  fmtstring(1:3)=' (A'
!  strlen = LEN(TRIM(outputstring))
!  IF(strlen.lt.10)THEN
!    WRITE(fmtstring(4:4),'(i1)')strlen
!    fmtstring(5:6)=')'
!  ELSEIF(strlen.lt.100)THEN
!    WRITE(fmtstring(4:5),'(i2)')strlen
!    fmtstring(6:7)=')'
!  ELSEIF(strlen.lt.1000)THEN
!    WRITE(fmtstring(4:6),'(i3)')strlen
!    fmtstring(7:8)=')'
!  ELSEIF(strlen.lt.10000)THEN
!    WRITE(fmtstring(4:7),'(i4)')strlen
!    fmtstring(8:9)=')'
!  ELSE
!    WRITE(fmtstring(4:8),'(i5)')strlen
!    fmtstring(9:10)=')'
!  ENDIF
!  WRITE(fid_enkf_final,TRIM(fmtstring))TRIM(outputstring)
!           
!  DO j=1,EnkfInfo%nE
!    outputstring(:)=' '
!    k=1
!    DO i=1,EnkFInfo%nX
!      WRITE(outputstring(k:k+9),'(f10.3)')EnkfX%y(i,j)
!      k=k+10
!    ENDDO
!    DO i=1,EnkFInfo%nP
!      WRITE(outputstring(k:k+13),'(f14.5)')EnkfP%y(i,j)
!      k=k+14
!    ENDDO
!    outputstring = TRIM(adjustr(outputstring(1:k)))
!    fmtstring(1:3)=' (A'
!    strlen = LEN(TRIM(outputstring))
!    IF(strlen.lt.10)THEN
!      WRITE(fmtstring(4:4),'(i1)')strlen
!      fmtstring(5:6)=')'
!    ELSEIF(strlen.lt.100)THEN
!      WRITE(fmtstring(4:5),'(i2)')strlen
!      fmtstring(6:7)=')'
!    ELSEIF(strlen.lt.1000)THEN
!      WRITE(fmtstring(4:6),'(i3)')strlen
!      fmtstring(7:8)=')'
!    ELSEIF(strlen.lt.10000)THEN
!      WRITE(fmtstring(4:7),'(i4)')strlen
!      fmtstring(8:9)=')'
!    ELSE
!      WRITE(fmtstring(4:8),'(i5)')strlen
!      fmtstring(9:10)=')'
!    ENDIF
!    WRITE(fid_enkf_final,TRIM(fmtstring))TRIM(outputstring)
!  ENDDO
!  close(fid_enkf_final)
!        
!  ! deallocation of all the ensemble data
!  CALL deallocate_ensemble(EnkfX)
!  CALL deallocate_ensemble(EnkfA)
!  CALL deallocate_ensemble(EnkfP)
!  CALL deallocate_ensemble(EnkfF)
!  CALL deallocate_ensemble(EnkfD)
!  CALL deallocate_ensemble(EnkfHX)
!    
!  ! shutdown random number stream
!  errcode = deinitialize_random_stream()
!    
!  ! free MKL memory usage
!  CALL MKL_FREE_BUFFERS()
!
!  ! close the time series outputs
!  !close(fid_enkf_min)
!  !close(fid_enkf_max)
!  close(fid_enkf_median)
!  !close(fid_enkf_mean)
!  close(fid_enkf_min95)
!  close(fid_enkf_max95)
!  !close(fid_enkf_std)
!end subroutine enkf_deallocation


!!***********************************************************************************************
!!> deallocate_ensemble
!!> purpose: DEALLOCATE all ALLOCATED field in a TYPE(enkf_ensemble_data) variable
!subroutine deallocate_ensemble(ensemble)
!  TYPE(enkf_ensemble_data) :: ensemble
!
!  ! DEALLOCATE VARIOUS FIELDS IN THE VARIABLE
!  IF(ALLOCATED(ensemble%y))DEALLOCATE(ensemble%y)
!  IF(ALLOCATED(ensemble%mi))DEALLOCATE(ensemble%mi)
!  IF(ALLOCATED(ensemble%xi))DEALLOCATE(ensemble%xi)
!  IF(ALLOCATED(ensemble%ensgen))DEALLOCATE(ensemble%ensgen)
!  IF(ALLOCATED(ensemble%kalman))DEALLOCATE(ensemble%kalman)       
!  IF(ALLOCATED(ensemble%mean))DEALLOCATE(ensemble%mean)            
!  IF(ALLOCATED(ensemble%sigma))DEALLOCATE(ensemble%sigma)         
!  IF(ALLOCATED(ensemble%covar))DEALLOCATE(ensemble%covar)         
!  IF(ALLOCATED(ensemble%minimum))DEALLOCATE(ensemble%minimum)      
!  IF(ALLOCATED(ensemble%maximum))DEALLOCATE(ensemble%maximum)   
!  ! IF(ALLOCATED(ensemble%statout))DEALLOCATE(ensemble%statout) 
!  ! IF(ALLOCATED(ensemble%ensout))DEALLOCATE(ensemble%ensout)        
!  IF(ALLOCATED(ensemble%outmean))DEALLOCATE(ensemble%outmean)       
!  IF(ALLOCATED(ensemble%outquant))DEALLOCATE(ensemble%outquant)       
!  IF(ALLOCATED(ensemble%outmin))DEALLOCATE(ensemble%outmin)       
!  IF(ALLOCATED(ensemble%outmax))DEALLOCATE(ensemble%outmax)       
!  IF(ALLOCATED(ensemble%outsigma))DEALLOCATE(ensemble%outsigma)       
!  ! IF(ALLOCATED(ensemble%outcovar))DEALLOCATE(ensemble%outcovar)       
!  ! IF(ALLOCATED(ensemble%outcorr))DEALLOCATE(ensemble%outcorr)       
!end subroutine deallocate_ensemble


!! The update_ensemble_after_kalmanfilter is not needed anymore!!
!! However, parameter resampling following mooradkhanii need to be implemented, perhaps

!! ---------------------------------------------------------------------------------
!!> update_ensemble_after_kalmanfilter
!!> purpose: update especially the parameter ensembles after the Kalman filter,
!!>          but also the A and the F ensembles are updated if members were included in the filter
!! ---------------------------------------------------------------------------------
!subroutine update_ensembles_after_Kalmanfilter()
!  !local variables
!  INTEGER i
!  
!  ! Parameters
!  DO i = 1, EnkfP%n
!    IF(EnkfP%Kalman(i).eq.1)THEN
!      EnkfP%y(i,:) = enkfX%y(EnkfP%xi(i),:)
!    ENDIF
!  ENDDO
!  CALL enkf_checkminmax(EnkfP)
!  ! Resamplig parameters by method of Moradkhani et al (2004):
!  ! CALL enkf_resample_P()
!
!  ! additional variables
!  DO i = 1, EnkfA%n
!    IF(EnkfA%Kalman(i).eq.1)THEN
!      EnkfA%y(i,:) = enkfX%y(EnkfA%xi(i),:)
!    ENDIF
!  ENDDO
!  CALL enkf_checkminmax(EnkfA)
!
!  ! forcing variables
!  DO i = 1, EnkfF%n
!      IF(EnkfF%Kalman(i).eq.1)THEN
!          EnkfF%y(i,:) = enkfX%y(EnkfF%xi(i),:)
!      ENDIF
!  ENDDO
!  CALL enkf_checkminmax(EnkfF)
!
!  ! Finally, check minmax also in X
!  CALL enkf_checkminmax(EnkfX)
!
!end subroutine update_ensembles_after_Kalmanfilter

!! The FUNCTION below is now replaced by a CALL to subroutine enkf_ensemble_statistics within the ENKF analysis

!!----------------------------------------------------------------------------------------------
!!> enkf_allensemble_statistics
!!> purpose: calculate statistics from the ensemble (mean, median, min/max, 95% quantile, stdev)
!!----------------------------------------------------------------------------------------------
!subroutine enkf_allensemble_statistics()
!  ! local variables
!  INTEGER DIM, NN, i, i_ens
!
!  i_ens = 1
!    
!  ! 1) model states, X ensemble
!  !IF(EnkfX%n.GT.0 .and. EnkfInfo%XS)CALL enkf_ensemble_statistics(EnkfX%y, EnkfX%n, EnkfX%ne, EnkfX%outmean, EnkfX%outquant, EnkfX%outmin, EnkfX%outmax, EnkfX%outsigma )!, EnkfX%outcovar, EnkfX%outcorr)
!  !IF(EnkfX%n.GT.0)CALL enkf_ensemble_statistics(EnkfX%y, EnkfX%n, EnkfX%ne, EnkfX%outmean, EnkfX%outquant, EnkfX%outmin, EnkfX%outmax, EnkfX%outsigma )!, EnkfX%outcovar, EnkfX%outcorr)
!    
!  ! 2) model fluxes, A ensemble
!  IF(EnkfA%n.GT.0)CALL enkf_ensemble_statistics(EnkfA%y, EnkfA%n, EnkfA%ne, EnkfA%outmean, EnkfA%outquant, EnkfA%outmin, EnkfA%outmax, EnkfA%outsigma )!, EnkfA%outcovar, EnkfA%outcorr)
!    
!  ! 4) model parameters from P ensemble
!  IF(EnkfP%n.GT.0)CALL enkf_ensemble_statistics(EnkfP%y, EnkfP%n, EnkfP%ne, EnkfP%outmean, EnkfP%outquant, EnkfP%outmin, EnkfP%outmax, EnkfP%outsigma )!, EnkfP%outcovar, EnkfP%outcorr)
!   
!  ! 5) forcing from F ensemble
!  IF(EnkfF%n.GT.0)CALL enkf_ensemble_statistics(EnkfF%y, EnkfF%n, EnkfF%ne, EnkfF%outmean, EnkfF%outquant, EnkfF%outmin, EnkfF%outmax, EnkfF%outsigma )!, EnkfF%outcovar, EnkfF%outcorr)
!
!end subroutine enkf_allensemble_statistics


!! The FUNCTION beow is still used, however modified for the new data structures

!!-------------------------------------------------------------------------
!!> enkf_ensemble_output
!!> PURPOSE: print output from the Enkf ensembles:
!!-------------------------------------------------------------------------
!subroutine enkf_ensemble_output()
!  ! mean, median, 95% quantiles, min, max, stdev in 7 text files:
!  ! if requested print entire ensembles in binary files
!  !
!  ! write to these files... fid_enkf_min, fid_enkf_max, fid_enkf_median, fid_enkf_final, fid_enkf_mean, fid_enkf_min95, fid_enkf_max95,fid_enkf_std
!  !
!  ! LOCAL VARIABLES:
!  INTEGER :: errcode,i,j,k,strlen
!  CHARACTER(LEN=30000) :: outputstringMEDIAN,outputstringMIN95,outputstringMAX95
!  CHARACTER(LEN=30000) :: outputstringSTDEV
!  CHARACTER(LEN=20)    :: fmtstring
!
!  ! ADD DATA TO OUTPUTSTRING:
!  outputstringMEDIAN(:) =' '
!  outputstringMIN95(:)  =' '
!  outputstringMAX95(:)  =' '
!  outputstringSTDEV(:)  =' '
!  k=1
!    
!  ! 1) model state variables from the X ensemble
!  IF(EnkfInfo%XS)THEN
!    DO i=1,EnkFInfo%nX
!      WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfX%outquant(2,i)
!      WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfX%outquant(1,i)
!      WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfX%outquant(3,i)
!      !WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfX%outcovar(i,i))**0.5
!      WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfX%outsigma(i))
!      k=k+10
!    ENDDO
!  ENDIF
!  ! 2) model auxiliary variables from A ensemble
!  DO i=1,EnkFInfo%nA
!    WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfA%outquant(2,i)
!    WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfA%outquant(1,i)
!    WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfA%outquant(3,i)
!    !WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfA%outcovar(i,i))**0.5
!    WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfA%outsigma(i))
!    k=k+10
!  ENDDO
!    
!  ! 3) model parameters from P ensemble
!  DO i=1,EnkFInfo%nP
!    WRITE(outputstringMEDIAN(k:k+13),'(f14.6)') EnkfP%outquant(2,i)
!    WRITE(outputstringMIN95(k:k+13), '(f14.6)') EnkfP%outquant(1,i)
!    WRITE(outputstringMAX95(k:k+13), '(f14.6)') EnkfP%outquant(3,i)
!    !WRITE(outputstringSTDEV(k:k+13), '(f14.6)') (EnkfP%outcovar(i,i))**0.5
!    WRITE(outputstringSTDEV(k:k+13), '(f14.6)') (EnkfP%outsigma(i))
!    k=k+14
!  ENDDO
!   
!  ! 4) forcing from F ensemble
!!    DO i=1,EnkFInfo%nF
!!        WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfF%outquant(2,i)
!!        WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfF%outquant(1,i)
!!        WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfF%outquant(3,i)
!!!        WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfF%outcovar(i,i))**0.5
!!        WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfF%outsigma(i))
!!        k=k+10
!!    ENDDO
!
!  ! 5) "Innovations" (aka mod-obs) from D ensemble (check the enkf analysis scheme...)
!  !j=0
!  !DO i=1,EnkFInfo%nD
!  !    if ()THEN
!  !        j=j+1
!  !    WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfD%outquant(2,i)
!  !    WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfD%outquant(1,i)
!  !    WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfD%outquant(3,i)
!  !    WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfD%outcovar(i,i))**0.5
!  !    k=k+10
!  !    ELSE
!  !        WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') -9999.
!  !        WRITE(outputstringMIN95(k:k+9), '(f10.4)') -9999.
!  !        WRITE(outputstringMAX95(k:k+9), '(f10.4)') -9999.
!  !        WRITE(outputstringSTDEV(k:k+9), '(f10.4)') -9999.
!  !    ENDIF
!  !ENDDO
!
!
!  ! PREPARE OUTPUT FORMAT
!  outputstringMEDIAN = TRIM(adjustr(outputstringMEDIAN(1:k)))
!  outputstringMIN95 = TRIM(adjustr(outputstringMIN95(1:k)))
!  outputstringMAX95 = TRIM(adjustr(outputstringMAX95(1:k)))
!  outputstringSTDEV = TRIM(adjustr(outputstringSTDEV(1:k)))
!
!  fmtstring = enkf_get_fmtstring(outputstringMEDIAN)
!  WRITE(fid_enkf_median,TRIM(fmtstring))TRIM(outputstringMEDIAN)
!
!  fmtstring = enkf_get_fmtstring(outputstringMIN95)
!  WRITE(fid_enkf_min95,TRIM(fmtstring))TRIM(outputstringMIN95)
!
!  fmtstring = enkf_get_fmtstring(outputstringMAX95)
!  WRITE(fid_enkf_max95,TRIM(fmtstring))TRIM(outputstringMAX95)
!
!  fmtstring = enkf_get_fmtstring(outputstringSTDEV)
!  WRITE(fid_enkf_std,TRIM(fmtstring))TRIM(outputstringSTDEV)
!end subroutine enkf_ensemble_output
!
!!-------------------------------------------------------------------------
!!> enkf_open_timeseriesoutput
!!-------------------------------------------------------------------------
!subroutine enkf_open_timeseriesoutput(resultdir)
!  CHARACTER(LEN=200) resultdir
!    
!  open(fid_enkf_median,file=TRIM(resultdir)//'enkf_output_median.txt',STATUS = 'unknown')
!  open(fid_enkf_min95,file=TRIM(resultdir)//'enkf_output_min95.txt',STATUS = 'unknown')
!  open(fid_enkf_max95,file=TRIM(resultdir)//'enkf_output_max95.txt',STATUS = 'unknown')
!  open(fid_enkf_std,file=TRIM(resultdir)//'enkf_output_std.txt',STATUS = 'unknown')
!  WRITE(fid_enkf_median,*)'% ENKF MEDIAN OUTPUT:'
!  IF(EnkfInfo%XS)THEN
!    WRITE(fid_enkf_median,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_min95,*)'% ENKF 2.5 quantile OUTPUT:'
!    WRITE(fid_enkf_min95,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_max95,*)'% ENKF 97.5 quantile OUTPUT:'
!    WRITE(fid_enkf_max95,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_std,*)'% ENKF standard deviation OUTPUT:'
!    WRITE(fid_enkf_std,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!  ELSE
!    WRITE(fid_enkf_median,*)'% Timestep A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_min95,*)'% ENKF 2.5 quantile OUTPUT:'
!    WRITE(fid_enkf_min95,*)'% Timestep A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_max95,*)'% ENKF 97.5 quantile OUTPUT:'
!    WRITE(fid_enkf_max95,*)'% Timestep A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_std,*)'% ENKF standard deviation OUTPUT:'
!    WRITE(fid_enkf_std,*)'% Timestep A(1:na) P(1:np) F(1:nf)'    
!  ENDIF
!end subroutine enkf_open_timeseriesoutput
!
!subroutine enkf_close_timeseriesoutput()
!  close(fid_enkf_median)
!  close(fid_enkf_min95)
!  close(fid_enkf_max95)
!  close(fid_enkf_std)
!end subroutine enkf_close_timeseriesoutput
!
!FUNCTION enkf_get_fmtstring(outputstri) RESULT(fmtstri)
!  CHARACTER(LEN=10000) outputstri
!  CHARACTER(LEN=20)    :: fmtstri
!  INTEGER strlen
!
!  fmtstri(1:3)=' (A'
!  strlen = LEN(TRIM(outputstri))
!
!  IF(strlen.lt.10)THEN
!    WRITE(fmtstri(4:4),'(i1)')strlen
!    fmtstri(5:6)=')'
!  ELSEIF(strlen.lt.100)THEN
!    WRITE(fmtstri(4:5),'(i2)')strlen
!    fmtstri(6:7)=')'
!  ELSEIF(strlen.lt.1000)THEN
!    WRITE(fmtstri(4:6),'(i3)')strlen
!    fmtstri(7:8)=')'
!  ELSEIF(strlen.lt.10000)THEN
!    WRITE(fmtstri(4:7),'(i4)')strlen
!    fmtstri(8:9)=')'
!  ELSE
!    WRITE(fmtstri(4:8),'(i5)')strlen
!    fmtstri(9:10)=')'
!  ENDIF
!END FUNCTION enkf_get_fmtstring

!!----------------------------------------------------------------------------
!!> get_spatially_correlated_random_matrix
!!>
!!> Function for spatially correlated random variable generation
!!>
!!> 1) Random numbers are generated on a regular 2D grid using FFT
!!>    (with variance=1, mean=0, and user defined correlation length) 
!!> 
!!> 2) The random numbers are interpolated to subbasins using bilinear interpolation 
!!>
!!> 3) Random numbers are scaled with the user defined standard deviation
!!>
!!> 4) The scaled random numbers are added to the mean observations 
!!>
!!> Author: David and Nils Gustafsson, 2013
!!----------------------------------------------------------------------------
!subroutine get_spatially_correlated_random_matrix(ensemble) !n,nv,a,sigma,r)
!  TYPE(enkf_ensemble_data), INTENT(INOUT) :: ensemble
!  REAL(kind = 4), ALLOCATABLE :: r(:)
!  INTEGER :: j, i, errcode,fieldsize,iens,jx
!     
!  ! 1) Random numbers are generated on a regular 2D grid using FFT
!  !    (with variance=1, mean=0, and user defined correlation length)
!  !
!  ! 2) The random numbers are interpolated to subbasins using bilinear interpolation
!  !
!  ! 3) Random numbers are scaled with the user defined standard deviation
!  !
!  ! 4) The scaled random numbers are added to the mean observations 
!    
!  ! some preparations
!  fieldsize = ensemble%n/ensemble%ntypes
!    
!  ALLOCATE(r(fieldsize))  ! ALLOCATE r(fieldsize = nsub)
!    
!  !Loop over variable types
!  DO j=1,ensemble%ntypes
!    ! loop over ensemble member iens in ensemble type j
!    DO iens=1,ensemble%ne
!      ! generate random numbers with standard deviation 1 and correlation length L
!      !
!      ! this is still just a prototype where we just assume correlation 1, and thefore only need to draw 
!      ! one normally distributed random number for each ensemble member
!      CALL get_random_number(0.,1.,r(1))
!      IF(fieldsize.GE.2)THEN
!        DO i=2,fieldsize
!          r(i) = r(1)
!        ENDDO
!      ENDIF
!      ! loop over records in the field and apply the random perturbation to the mean, scaled with it's individual sigma
!      DO i=1,fieldsize
!        !IF(ensemble%mi(i).EQ.ensemble%types(j)) ensemble%y(i,:) = r(:) * ensemble%sigma(i)  + ensemble%mean(i)
!        jx = i+(j-1)*fieldsize
!        ensemble%y(jx,iens) = r(i) * ensemble%sigma(jx)  + ensemble%mean(jx)
!      ENDDO
!    ENDDO
!  ENDDO
!end subroutine get_spatially_correlated_random_matrix
  
  SUBROUTINE find_intarrloc(intarr,arrdim,intval,arrloc)
    !input arguments
    INTEGER, INTENT(IN) :: intarr(:)
    INTEGER, INTENT(IN) :: arrdim
    INTEGER, INTENT(IN) :: intval
    INTEGER, INTENT(OUT):: arrloc
    !local variable
    INTEGER i
    !find first location of intval in intarr
    arrloc=-1
    DO i=1,arrdim
      IF(intarr(i).EQ.intval.AND.arrloc.EQ.-1)THEN
        arrloc = i
      ENDIF
    ENDDO
  END SUBROUTINE find_intarrloc
  
END MODULE ASSIMILATION_ROUTINES
