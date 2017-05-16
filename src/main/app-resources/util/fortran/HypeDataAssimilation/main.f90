!> \file main.f90
!> Contains main program of HYSS-HYPE

!> Main program for HYSS - Hydrological Simulation System
!> 
!> The system features hydrological simulation of soil and surface water system, 
!> calibration of parameters and criteria calculations, updating of flow and state 
!> to observed values, ensemble simulation, and more.
PROGRAM MAIN

!Copyright 2011-2017 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!-----------------------------------------------------------------------------------------

!Used modules
!      USE IFQWIN                               !for nobutton
      USE MODELMODULE, ONLY : model, &
                              model_version_information, &
                              initiate_output_variables, &
                              initiate_model,            &
                              initiate_model_state, &
                              initiate_model_parameters, &
                              load_modeldefined_input
      USE STATETYPE_MODULE
      USE WORLDVAR, ONLY :  writematlab,       &
                            writeload,         &
                            output,            &
                            noutput,           &
                            nacrit,            &
                            nsubCrit,          &
                            ndt,               &
                            fileunit_temp,     &
                            allocate_accumulation,     &
                            numqobsstn,        &
                            maxcharpath,       &
                            simsequence,       &
                            infodir,           &
                            modeldir,          &
                            resdir,            &
                            bdate,             &
                            sdate,             &
                            outstartdate,      &
                            dtskip,            &
                            doopt,             &
                            numoptimpar,       &   
                            deallocate_worldvar,      &
                            deallocate_MCvariables,   &
                            optim,             &
                            bestMCoptcrit,     &
                            bestMCperformance, &
                            bestMCparameters,  &
                            maxperf,           &
                            maxsubass,         &
                            simsubmodel,       &
                            ibasemodel,     &
                            readsfobs,readswobs,  &
                            readwind, &
                            readhumid, &
                            readtminmaxobs, &
                            checkdata,  &
                            psdates,    &
                            lddates,    &
                            dddates,    &
                            optimStartTime, &      
                            optimFuncCall,  &
                            lineSearchCallCount, &
                            doassimilation, &
                            noutreg, &
                            reallocate_outvar_information
      USE MODVAR, ONLY : preci,tempi,qobsi,xobsi,  &
                         snowfraci,shortwavei,     &
                         windi,humidi, &
                         tmini,tmaxi,xoregobsi,    &
                         nsub,ncrop,        &
                         nsub_basemodel,       &
                         nclass,numsubstances, &
                         naquifers,   &
                         maxsoillayers,   &
                         max_classoutvar,   &
                         max_basinoutvar,   &
                         max_noutvar, &
                         doupdate,i_qar,i_war,  &
                         conductN,conductP,conductC,  &
                         conductxoms,conductregest,   &
                         i_t1,i_t2, &
                         wetlandexist,glacierexist,  &
                         doirrigation,  &
                         conductflood,  &
                         modeloption,p_lakeriverice, &
                         p_growthstart,  &
!                         allocate_outvar,   &
                         deallocate_modvar,  &
                         dimriverlag, &
                         currentdate, &
                         timesteps_per_day, &
                         noutvar, &
                         nrivertypes,nlaketypes
      USE COMPOUT, ONLY : compute_mapoutput,          &
                          compute_outloads,         &
                          prepare_to_compute_crit,  &
                          calculate_criteria
      USE TIMEROUTINES, ONLY : calculate_time_for_model
      USE READWRITE_ROUTINES
      USE LIBDATE, ONLY : DateType, OPERATOR(.EQ.)
      USE DATAMODULE
      USE STATE_DATAMODULE, ONLY : initiate_state_for_submodel, &
                                   load_saved_state,  &
                                   finalize_outstate
#ifdef _ASSIMILATION == 1
      !use the Data Assimilation modules (pre-compiler flag ASSIMILATION to be set in Visual Studio project settings or in the makefile)
      USE ASSIMILATION_INTERFACE
      USE ASSIMILATION_ROUTINES
      USE ASSIMILATION_VARIABLES
#endif

      IMPLICIT NONE

!Parameter declarations
      INTEGER, PARAMETER :: maxoutstates = 10     !Max number of dates for saving state

!Variable declarations
      TYPE(DateType) d            !Current time
      INTEGER idt                 !Current timestep
      INTEGER ivar                !Current output variable
      INTEGER iens                !Current ensemble being simulated
      INTEGER iout                !Current output
!      INTEGER nobutton            !No exit window
      LOGICAL pwrite              !Flag for periodend, time to write to file
      CHARACTER(LEN=maxcharpath+25) filename  !hyss filename
      CHARACTER(LEN=8)  :: logdate  !Date for log-file name
      CHARACTER(LEN=10) :: logtime  !Time for log-file name
      CHARACTER(LEN=3)  :: logseq   !Seqnr for log-file name
      INTEGER :: datim(8) 
      INTEGER :: oldyear            !year of last time step

      REAL, ALLOCATABLE :: par(:)
      REAL optcrit, condcrit, condthres
      REAL, ALLOCATABLE :: basincrit(:,:,:)   !R2, CC, RE, RSDE, QC, QR, STDC, STDR, MAE, RMSE, Bias, STDbias, KGE, KGEpartSTD, KGEpartMM, NRMSE per subbasin och kriterie
      REAL, ALLOCATABLE :: simperformance(:,:)   !rr2,sr2,wr2,rmae,sbias,rrve,wrve,rra,sra,meanRA,tau,medianr2,medianra,meanrs,meancc,mediankg,meanabsre
      INTEGER :: n_Result       !Subroutine return status, TODO: Change name to status
      INTEGER npar              !Number of parameters to be calibrated (couted in file)
      INTEGER :: nmapperiod     !Number of periods for map print out
      LOGICAL :: stateinput                               !Code for reading state
      INTEGER :: numoutstates                             !Number of dates for saving state
      TYPE(DateType) :: stateoutadate(maxoutstates)       !Dates for saving state
      TYPE(DateType) :: prestateoutdate(maxoutstates)     !Day before date for saving state (or 0)
      
!Variables for updating of Q and W
      LOGICAL :: quseobsallstations, quseobsnostations
      LOGICAL :: qarnostations, warnostations
      LOGICAL :: wendupdallstations, wendupdnostations
      CHARACTER(LEN=4) :: wobsvarname
      
!Model state variables and other saved variables declaration
      TYPE(SNOWICESTATETYPE) :: frozenstate
      TYPE(SOILSTATETYPE)    :: soilstate      
      TYPE(AQUIFERSTATETYPE) :: aquiferstate      
      TYPE(RIVERSTATETYPE)   :: riverstate      
      TYPE(LAKESTATETYPE)    :: lakestate
      TYPE(MISCSTATETYPE)    :: miscstate

#ifdef _ASSIMILATION == 1
      !Declaration of some data assimilation variables
      INTEGER assim_ens_size  !number of ensemble members
      INTEGER iassim     !loop-variable for ensembles members
      TYPE(STATEINFOTYPE),ALLOCATABLE :: stateinfo(:)
      ! variables for timing the assimilation code
      real :: start_time, stop_time, total_time(20), total_time_stat(4)
      call cpu_time(start_time)
      total_time=0
      total_time_stat=0
#endif

!Program start
!>\b Algorithm \n

      CALL DATE_AND_TIME (logdate, logtime,values=datim)

!Current model domain
      CALL get_hyss_arguments(infodir,simsequence)   
      WRITE(logseq,'(I3.3)') simsequence
      WRITE(filename,'(a)') TRIM(infodir)//'hyss_'//logseq(1:3)//'_'//logdate(3:8)//'_'//logtime(1:4)//'.log'
      OPEN(UNIT=6,FILE=TRIM(filename))
      WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
               ' Job start date: ',datim(1),'-',datim(2),'-',datim(3),    &
               '  time: ',datim(5),':',datim(6),':',datim(7)
      WRITE(6,*) '---------------------------------------------------'
!      nobutton = SETEXITQQ(qwin$exitnopersist)    !nobutton version
      CALL model_version_information(6)   !Collect and print model version information
      CALL initiate_output_variables()    !Define model's output variables
      CALL initiate_model_parameters()    !Define model's parameters
!>Get information for this simulation
      CALL load_coded_info(infodir,n_Result,bdate,sdate,outstartdate,dtskip,ndt,numsubstances,  &
                           stateinput,maxoutstates,numoutstates,stateoutadate,prestateoutdate,  &
                           modeldir,resdir, &
                           quseobsallstations,quseobsnostations,qarnostations,warnostations,    &
                           wendupdallstations,wendupdnostations,wobsvarname)
      IF(n_Result/=0)THEN
        STOP 1
      ENDIF

!>Read input data
      CALL load_cropdata(modeldir,'CropData.txt',ncrop,n_Result)
      IF(n_Result.NE.0) STOP 1
      CALL load_basindata(modeldir,'GeoClass.txt','GeoData.txt',nsub_basemodel,n_Result)
      IF(n_Result.NE.0) STOP 1
      CALL load_pointsourcedata(modeldir,'PointSourceData.txt',nsub_basemodel,n_Result) 
      IF(n_Result.NE.0) STOP 1
      CALL load_branchdata(modeldir,n_Result)
      IF(n_Result.NE.0) STOP 1
      CALL load_aquiferdata(modeldir,nsub_basemodel,naquifers,n_Result) 
      IF(n_Result.NE.0) STOP 1
      CALL load_glacierdata(modeldir,nsub_basemodel,n_Result) 
      IF(n_Result.NE.0) STOP 1
      IF(checkdata(1,1)) CALL checkindata_part1()
      IF(checkdata(2,1)) CALL checkindata_part2()
      CALL checkindata_stop()

      CALL load_submodel_info(infodir,simsubmodel,nsub,n_Result)   !allocation and initialisation of ibasemodel
      IF(n_Result/=0) STOP 1
      CALL load_output_regions(modeldir,noutreg,n_Result)    !Prepare for outregions, read Outregions.txt
      IF(n_Result/=0) STOP 1

      CALL load_observations(modeldir,nsub_basemodel,bdate,sdate,ndt,n_Result)
      IF(n_Result.NE.0) STOP 1
      CALL load_parameters(modeldir,nsub_basemodel,'par.txt')

!Read base model defined input data and set base model configuration
      CALL set_model_base_configuration(nsub_basemodel,modeldir,dimriverlag,n_Result)
      IF(n_Result.NE.0) STOP 1
      
      IF(simsubmodel)THEN
        CALL reform_inputdata_for_submodel(nsub_basemodel,nsub,ibasemodel)
      ENDIF
      CALL calculate_path(nsub)

!Read model defined input data and set model configuration
      CALL load_modeldefined_input(modeldir,nsub_basemodel,nsub,ibasemodel,bdate,sdate,conductxoms,conductregest,n_Result)
      IF(n_Result.NE.0) STOP 1
      CALL set_model_configuration(glacierexist,doirrigation)

      CALL prepare_for_update(modeldir,wobsvarname,quseobsallstations, &
              quseobsnostations,qarnostations,warnostations,wendupdallstations,wendupdnostations,nsub)

!>Initialisations for memory allocation (states and output)
      CALL allocate_model_states(nsub_basemodel,numsubstances,nclass,naquifers, &
              maxsoillayers,nrivertypes,nlaketypes,dimriverlag,timesteps_per_day,conductN,conductP,conductC, & 
              conductflood, &
              i_t1>0,i_t2>0,wetlandexist,glacierexist,modeloption(p_lakeriverice)>=1,doirrigation, &
              doupdate(i_qar).OR.doupdate(i_war),modeloption(p_growthstart)==1, &
              frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

      CALL allocate_outvar(nsub,nclass,numsubstances,noutvar)   
      !IF(.NOT.doassimilation .AND. noutvar<max_noutvar) CALL reallocate_outvarinfo(noutvar,max_noutvar)
      IF(.NOT.doassimilation .AND. noutvar<max_noutvar) CALL reallocate_outvar_information(noutvar)
      IF(.NOT.doassimilation) CALL allocate_accumulation(nsub,nclass,numsubstances,max_classoutvar,max_basinoutvar)

!Allocate local variables
      nsubCrit = nsub
!      nsubCrit = nsub+noutreg
      ALLOCATE(basincrit(nsubCrit,maxsubass,nacrit))
      ALLOCATE(simperformance(maxperf,nacrit))

!Preparations for subbasin output
      CALL prepare_subbasin_output(n_Result)
      IF(n_Result/=0) STOP 1
      
!For data assimilation simulation, skip optimization and ordinary simulation code      
      IF(.NOT.doassimilation)THEN
      !CALL allocate_accumulation(nsub,nclass,numsubstances,max_classoutvar,max_basinoutvar) !Denna orsakar stackoverflow

        CALL DATE_AND_TIME (values=datim)
        WRITE(6,*)
        WRITE(6,*) '---------------------------------------------------'
        WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
                   ' Initialisations finished, calculations starts: ',datim(1),'-',datim(2),'-',datim(3),    &
                   '  time: ',datim(5),':',datim(6),':',datim(7)
        WRITE(6,*) '---------------------------------------------------'


!>Optimization
        IF(doopt)THEN
          CALL load_optpar(infodir)       !Reads optpar.txt and set all numerical optimization variables accordingly

          IF(optim%task_MC)THEN
            WRITE(6,*) 'MonteCarlo simulation with',optim%nruns_MC,'runs'
            CALL MonteCarlo_simulation(resdir,optim%task_writeall,frozenstate,  &
                 soilstate,aquiferstate,riverstate,lakestate,miscstate,npar)
            CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))    !save the best parameters to modpar
            CALL save_respar(resdir,npar,nsub)                       !and to file
          ENDIF

          IF(optim%task_boundps)THEN
            WRITE(6,*) 'Reduce bounds of parameter space for continued MonteCarlo simulation'
            CALL bounded_MonteCarlo_simulation(optim%task_MC,frozenstate,soilstate, &
                 aquiferstate,riverstate,lakestate,miscstate,npar)
            CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))    !save the best parameters to modpar
            CALL save_respar(resdir,npar,nsub)                       !and to file
          ENDIF

          lineSearchCallCount = 0
          optimFuncCall = 0
          CALL cpu_time(optimStartTime)

          IF(optim%task_Scanning) CALL param_scanning(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

          IF(optim%task_stageMC)THEN
            CALL stage_MonteCarlo(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
            CALL set_optim_modpar(numoptimpar,numoptimpar,bestMCparameters(1,:))    !save the best to modpar
            CALL save_respar(resdir,numoptimpar,nsub)                               !and file
          ENDIF
        
          IF(optim%task_BrentNew .OR. optim%task_stpstDesc .OR. optim%task_DFP .OR. optim%task_BFGS)THEN
            ALLOCATE(par(numoptimpar))
            CALL linesearch_methods_calibration(numoptimpar,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,par)
            CALL set_optim_modpar(numoptimpar,numoptimpar,par)    !save the best to modpar
            CALL save_respar(resdir,numoptimpar,nsub)                            !and file
          ENDIF
          IF(optim%task_DEMC)THEN
            WRITE(6,*) 'DEMC Differential-Evolution Markov Chain, with',optim%DEMC_npop, &
            'populations, and ',optim%DEMC_ngen,'generations (',optim%DEMC_ngen-1,  &
            'evolution steps)',', in total',optim%DEMC_npop*(optim%DEMC_ngen),'simulations.'
          
            CALL DEMC_simulation(resdir,optim%task_writeall,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,npar)
            CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))         !save the MEDIAN parameters to modpar (1 row in bestMC)
            CALL save_respar(resdir,npar,nsub)        
          ENDIF
          IF(optim%task_runens)THEN       !write best ensemble resultat to file
            CALL save_ensemble_simulations(resdir,numoptimpar,maxperf,nacrit, &
                 optim%nruns_best,bestMCoptcrit,bestMCperformance,bestMCparameters)
          ENDIF
        ENDIF

        IF(.NOT.optim%task_writesim)THEN
!Ensemble loop, simulate all ensemble members
          ensemble_loop:    &
       &  DO iens = 1, optim%nruns_best

!>Simulation start
            CALL prepare_outputfiles(resdir,nsub,naquifers,iens,optim%task_runens,optim%task_writesim)

  !        IF(doopt .AND. (optim%task_runens .OR. optim%task_DEMC))THEN
            IF(doopt .AND. optim%task_runens)THEN
              CALL set_optim_modpar(numoptimpar,numoptimpar,bestMCparameters(iens,:))   !set model parameters
            ENDIF
            CALL initiate_output_routines()     !All output accumulation variables zeroed

            !>Initial model calculations; initial states, parameters
            IF(simsubmodel)THEN
              CALL initiate_state_for_submodel(modeldir,ibasemodel,dimriverlag,stateinput,  &
                  frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
            ELSE
              IF(stateinput) THEN
                CALL load_saved_state(modeldir,nsub,dimriverlag,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
                WRITE(6,*)
                WRITE(6,*) 'Loading saved state.'
              ELSE
                CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
              ENDIF
              CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 
            ENDIF
            oldyear = 0

!>Time Loop:
            time_loop:    &
       &    DO idt = 1,ndt

              !>\li Get current input data
              CALL get_current_date_prec(idt,nsub,d,preci(1:nsub))
              CALL calculate_time_for_model(idt,d)
              CALL log_progress(oldyear,currentdate%year)
              tempi(1:nsub) = get_current_temp(idt,nsub)
              IF(ALLOCATED(qobsi)) qobsi(1:nsub) = get_current_flow(idt,d,nsub,numqobsstn)    !m3/s
              IF(ALLOCATED(xobsi)) xobsi(:) = get_current_otherobs(idt,d)
              IF(ALLOCATED(xoregobsi)) xoregobsi(:) = get_current_oregobs(idt,d)
              IF(readsfobs) snowfraci(1:nsub) = get_current_snowfallfrac(idt,nsub)
              IF(readswobs) shortwavei(1:nsub) = get_current_shortwave(idt,nsub)
              IF(readwind)  windi(1:nsub) = get_current_windspeed(idt,nsub)
              IF(readhumid) humidi(1:nsub) = get_current_humidity(idt,nsub)
              IF(readtminmaxobs) tmini(1:nsub) = get_current_tmin(idt,nsub)
              IF(readtminmaxobs) tmaxi(1:nsub) = get_current_tmax(idt,nsub)
              IF(ALLOCATED(psdates)) CALL get_current_pointsources(modeldir,'PointSourceData.txt',nsub,d,n_Result) 
              IF(n_Result.NE.0) STOP 1
  !            IF(ALLOCATED(lddates)) CALL get_current_lakedata(modeldir,'LakeData.txt',nsub,d,n_Result) 
              IF(n_Result.NE.0) STOP 1
  !            IF(ALLOCATED(dddates)) CALL get_current_damdata(modeldir,'DamData.txt',nsub,d,n_Result) 
              IF(n_Result.NE.0) STOP 1

              CALL initiate_outvar(idt)
              !>\li Calculate flows and update states
              CALL model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
        
              DO ivar = 1,numoutstates    !Write state output
                IF (d.EQ.prestateoutdate(ivar)) THEN
                  IF(simsubmodel)THEN
                    WRITE(6,*) 'State can not be saved when submodel is simulated'
                  ELSE
                    CALL finalize_outstate(resdir,dimriverlag,stateoutadate(ivar),  &
                         frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 
                  ENDIF
                ENDIF
              ENDDO

              CALL revise_outvar()    !calculate regional outvar and? upstream?
              CALL prepare_to_compute_crit(d,idt,ndt)

              !>\li Calculate and write time dependent output
              DO iout = 1,noutput
                IF(output(iout)%fileformat==4) CALL write_regionfiles(iout,idt,ndt,iens,d)
                IF(output(iout)%fileformat==1) CALL write_subbasinfiles(iout,idt,ndt,iens,d)
                IF(output(iout)%fileformat==3) CALL write_timefiles(iout,idt,ndt,iens,d)
                IF(output(iout)%fileformat==2) CALL compute_mapoutput(d,iout,idt,ndt,writematlab,nmapperiod)  !Save data for map output
              ENDDO
              IF(writeload)THEN
                CALL compute_outloads(d,pwrite,idt,ndt)       !Write yearly load total for all subbasins
                IF(pwrite) CALL save_loadfiles(resdir,currentdate%year)
              ENDIF

            ENDDO time_loop

!>Compute and write criteria
            IF(nacrit/=0) THEN
              CALL calculate_criteria(optcrit,basincrit,simperformance,condcrit,condthres)
              CALL write_simulation_assessment(resdir,iens,nacrit,optcrit,     &
                  simperformance,optim%task_runens,condcrit,condthres)
              CALL write_subbasin_assessment(resdir,nsubCrit,nacrit,basincrit,iens,optim%task_runens)
!              CALL write_subbasin_assessment(resdir,nsubCrit,basincrit,iens,optim%task_runens)
            ENDIF

!Save and close files or prepare them for next ensemble member simulation
            CALL close_outputfiles(nsub,naquifers,iens)
            IF(iens==optim%nruns_best)THEN
              CALL close_observations(modeldir)
            ELSE
              CALL reset_observations(modeldir,n_Result)
              IF(n_Result.NE.0) STOP 1
            ENDIF

!>Write results to files
            CALL save_mapfiles(resdir,nsub,nmapperiod,iens,optim%task_runens,optim%task_writesim)

          ENDDO ensemble_loop   !simulate ensembles
        ENDIF !.NOT.optim%task_writesim
      ENDIF !.NOT.doassimilate

  
#ifdef _ASSIMILATION == 1     
!===========================================================================
!Data assimilation Version of the Time-Ensemble Loops
!===========================================================================
      IF(doAssimilation)THEN
        WRITE(6,*)
        WRITE(6,*) '---------------------------------------------------'
        WRITE(6,*) 'HYPE Data Assimilation Simulation'
        WRITE(6,*) '---------------------------------------------------'
        WRITE(6,*)
        WRITE(6,*) ' [Data Assimilation] Initialization started...'
        WRITE(6,*)

!HYSS: Initializations of hyss variables
        !Note. The counter, "iens", is used in the standard code above for ensemble simulations and is set to 1 here as if this was a normal determinstic simulation.
        !In the assimilation section below, we use a dedicated counter called "iassim" for looping over the assimilation ensemble members.
        iens=1      
        
!HYSS: Initial model calculations; initial states, parameters
        IF(simsubmodel)THEN
          CALL initiate_state_for_submodel(modeldir,ibasemodel,dimriverlag,stateinput,  &
              frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
        ELSE
          IF(stateinput) THEN
            CALL load_saved_state(modeldir,nsub,dimriverlag,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
            WRITE(6,*)
            WRITE(6,*) 'Loading saved state.'
          ELSE
            CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
          ENDIF
          CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 
        ENDIF
        oldyear = 0
        
!ASSIMILATION: Initialization of Data Assimilation variables
        CALL set_stateinfo(nsub,numsubstances,nclass,naquifers,maxsoillayers, &   !CP161124 pointerstates applied
              nrivertypes,nlaketypes,dimriverlag,timesteps_per_day,stateinfo, &
              frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
        CALL assim_initialize(modeldir,myAssimData,assim_ens_size,resdir,stateinfo,n_Result)
        IF(n_Result/=0) STOP 1
        
        !Prepare files for output
        !The standard output files are the mean (or median) of the ensemble data.
        !Additional output in sequence named files are optional; output files with
        !ensemble statistics (min,max,quantiles) or all ensemble members individually. 
        !IF(noutvar<max_noutvar) CALL reallocate_outvarinfo(noutvar,max_noutvar)
        CALL reallocate_outvar(nsub,noutvar)   
        IF(noutvar<max_noutvar) CALL reallocate_outvar_information(noutvar)
        CALL reallocate_accumulation(nsub,myAssimData%info%nstatout)
        CALL prepare_outputfiles(resdir,nsub,naquifers,iens,.FALSE.,.FALSE.)  !for standard output meanORmedian results
        DO iens = 2,myAssimData%info%nstatout+1           !this is statistics and ensemble members
          CALL prepare_outputfiles(resdir,nsub,naquifers,iens,.TRUE.,.FALSE.)
        ENDDO
        CALL initiate_output_routines()     !All output accumulation variables zeroed
        iens = 1    !reset
        
!ASSIMILATION: Initialization finished log message
        CALL DATE_AND_TIME (values=datim)
        WRITE(6,*)
        WRITE(6,*) '---------------------------------------------------'
        WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
                 ' [Data Assimilation] Initialisations finished, calculations starts: ',datim(1),'-',datim(2),'-',datim(3),    &
                 '  time: ',datim(5),':',datim(6),':',datim(7)
        WRITE(6,*) '---------------------------------------------------'

        call cpu_time(stop_time)
        total_time(1)=stop_time-start_time
        
!ASSIMILATION: Time Loop:
        assim_time_loop: DO idt = 1,ndt

          call cpu_time(start_time)
          
!HYSS: Get current input data
          CALL get_current_date_prec(idt,nsub,d,preci(1:nsub))
          CALL calculate_time_for_model(idt,d)
          CALL log_progress(oldyear,currentdate%year)
	        WRITE(6,*)' [Data Assimilation] Start of new timestep at day no:', dayno  !dayno set in calculate_time_for_model
          tempi(1:nsub) = get_current_temp(idt,nsub)
          IF(ALLOCATED(qobsi)) qobsi(1:nsub) = get_current_flow(idt,d,nsub,numqobsstn)    !m3/s
          IF(ALLOCATED(xobsi)) xobsi(:) = get_current_otherobs(idt,d)
          IF(ALLOCATED(xoregobsi)) xoregobsi(:) = get_current_oregobs(idt,d)
          IF(readsfobs) snowfraci(1:nsub) = get_current_snowfallfrac(idt,nsub)
          IF(readswobs) shortwavei(1:nsub) = get_current_shortwave(idt,nsub)
          IF(readwind)  windi(1:nsub) = get_current_windspeed(idt,nsub)
          IF(readhumid) humidi(1:nsub) = get_current_humidity(idt,nsub)
          IF(readtminmaxobs) tmini(1:nsub) = get_current_tmin(idt,nsub)
          IF(readtminmaxobs) tmaxi(1:nsub) = get_current_tmax(idt,nsub)
          IF(ALLOCATED(psdates)) CALL get_current_pointsources(modeldir,'PointSourceData.txt',nsub,d,n_Result) 
          IF(n_Result.NE.0) STOP 1

          call cpu_time(stop_time)
          total_time(2)=total_time(2)+stop_time-start_time
          start_time=stop_time
!ASSIMILATION: Generate input ensembles for this timestep
          CALL generate_forcing_ensemble(myAssimData)
          call cpu_time(stop_time)
          total_time(3)=total_time(3)+stop_time-start_time
          start_time=stop_time
!ASSIMILATION: Generate observation ensemble if analysis is enabled for this time step
          IF(idt>dtskip)THEN
            CALL generate_observation_ensemble(myAssimData)
          ENDIF
          call cpu_time(stop_time)
          total_time(4)=total_time(4)+stop_time-start_time
          start_time=stop_time
!ASSIMILATION: Loop over assimilation ensemble members for this timestep
          assim_ensemble_loop: DO iassim = 1,assim_ens_size

!ASSIMILATION: write states, forcing, parameters etc, for the current (iassim) ensemble member to the HYPE model data
            CALL ensemble_to_model2(iassim,myAssimData%info%nX,myAssimData%info%nA,myAssimData%info%nF, & 
                                    myAssimData%X,myAssimData%A,myAssimData%F,stateinfo)
             call cpu_time(stop_time)
          total_time(5)=total_time(5)+stop_time-start_time
          start_time=stop_time   
!HYSS: Initiate the outvars for the current (iassim) ensemble member
            CALL initiate_outvar(idt)
call cpu_time(stop_time)
          total_time(6)=total_time(6)+stop_time-start_time
          start_time=stop_time
!HYSS: Run the model for the current (iassim) ensemble member
            CALL model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  	            call cpu_time(stop_time)
          total_time(6)=total_time(6)+stop_time-start_time
          start_time=stop_time
!ASSIMILATION: Write model 'forecast' back to the ensemble data
            CALL model_to_ensemble2(iassim,myAssimData%X,myAssimData%A,myAssimData%info%nA,stateinfo)
            call cpu_time(stop_time)
          total_time(7)=total_time(7)+stop_time-start_time
          start_time=stop_time
!ASSIMILATION: Write the predicted observations to the ensemble if analysis is enabled (idt>dtskip)
            IF(idt>dtskip)THEN
              CALL modelobservations_to_ensemble(iassim,myAssimData)
            ENDIF
          call cpu_time(stop_time)
          total_time(8)=total_time(8)+stop_time-start_time
          start_time=stop_time
          ENDDO assim_ensemble_loop

!ASSIMILATION, Ensemble Kalman Filter Analysis (later generalize to any DA filter method) - ensemble state matrix will be updated if ensemble size>1 enkf analysis enabled (idt>dtskip) 
          IF(idt>dtskip.AND.myAssimData%info%nE.GT.1)THEN
            CALL enkf_analysis_main(myAssimData)
            call cpu_time(stop_time)
          total_time(9)=total_time(9)+stop_time-start_time
          start_time=stop_time
          ENDIF

!ASSIMILATION, Update ensemble statistics
          CALL updateEnsembleStatistics(myAssimData, total_time_stat)

          call cpu_time(stop_time)
          total_time(10)=total_time(10)+stop_time-start_time
          start_time=stop_time

!ASSIMILATION: Write ensemble mean (or median) to model variables, to be included in the standard output files.
          CALL meanORmedian_to_model2(myAssimData%info%nA,myAssimData%X,myAssimData%A,myAssimData%info%meanout,stateinfo)

          call cpu_time(stop_time)
          total_time(11)=total_time(11)+stop_time-start_time
          start_time=stop_time
          
!ASSIMILATION: Reset selected assimilation ensembles to ensemble mean (all variables not included in the "control vector" will be reset to ensemble mean after each timestep) 
          CALL meanORmedian_to_ensemble(myAssimData%info%nX,myAssimData%X,myAssimData%info%meanout)
          
          call cpu_time(stop_time)
          total_time(12)=total_time(12)+stop_time-start_time
          start_time=stop_time
          
!HYSS: Write state output files
          DO ivar = 1,numoutstates
            IF (d.EQ.prestateoutdate(ivar)) THEN
              IF(simsubmodel)THEN
                WRITE(6,*) 'Warning: State can not be saved when submodel is simulated'
              ELSE
                CALL finalize_outstate(resdir,dimriverlag,stateoutadate(ivar), &
                                       frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 
              ENDIF
            ENDIF
          ENDDO
          
          call cpu_time(stop_time)
          total_time(13)=total_time(13)+stop_time-start_time
          start_time=stop_time
          
!HYSS: Preparatory calculations for criteria (update various sums, square sums and counters with current timestep data) 
          CALL revise_outvar()
          CALL prepare_to_compute_crit(d,idt,ndt)

          call cpu_time(stop_time)
          total_time(14)=total_time(14)+stop_time-start_time
          start_time=stop_time

!HYSS: Calculate and write time dependent standard output
          DO iout = 1,noutput
            IF(output(iout)%fileformat==4) CALL write_regionfiles(iout,idt,ndt,iens,d)
            IF(output(iout)%fileformat==1) CALL write_subbasinfiles(iout,idt,ndt,iens,d)
            IF(output(iout)%fileformat==3) CALL write_timefiles(iout,idt,ndt,iens,d)
            IF(output(iout)%fileformat==2) CALL compute_mapoutput(d,iout,idt,ndt,writematlab,nmapperiod)
          ENDDO
          IF(writeload)THEN
            CALL compute_outloads(d,pwrite,idt,ndt)       !Write yearly load total for all subbasins
            IF(pwrite) CALL save_loadfiles(resdir,currentdate%year)
          ENDIF

          call cpu_time(stop_time)
          total_time(15)=total_time(15)+stop_time-start_time
          start_time=stop_time
          
!ASSIMILATION: Output of "statistical" and ensemble member simulations result
          IF(idt>dtskip)THEN
            DO iens = 2,myAssimData%info%nstatout+1
              CALL statistics_to_modeloutput(myAssimData%info%nA,myAssimData%info%nE,myAssimData%A,iens)
              DO iout = 1,noutput
                IF(output(iout)%fileformat==1) CALL write_subbasinfiles_in_parallel(iout,idt,ndt,iens,d)
                IF(output(iout)%fileformat==3) CALL write_timefiles_in_parallel(iout,idt,ndt,iens,d)
                IF(output(iout)%fileformat==4) CALL write_regionfiles_in_parallel(iout,idt,ndt,iens,d)
              ENDDO
            ENDDO
            iens = 1
            CALL statistics_to_modeloutput(myAssimData%info%nA,myAssimData%info%nE,myAssimData%A,iens,myAssimData%info%meanout)
          ENDIF


          call cpu_time(stop_time)
          total_time(16)=total_time(16)+stop_time-start_time
          start_time=stop_time

         write(6,*)'---------------------'
         write(6,*)'total_time(1:20)'
         write(6,*)total_time
         write(6,*)'total_time_stat(1:4)'
         write(6,*)total_time_stat
          
!ASSIMILATION: end of assimilation time loop
        ENDDO assim_time_loop

!HYSS: Compute and write criteria (nb! iens should still be set to 1) 
        IF(nacrit/=0) THEN
          CALL calculate_criteria(optcrit,basincrit,simperformance,condcrit,condthres)
          CALL write_simulation_assessment(resdir,iens,nacrit,optcrit,     &
                simperformance,.FALSE.,condcrit,condthres)
          CALL write_subbasin_assessment(resdir,nsubCrit,nacrit,basincrit,iens,.FALSE.)
!          CALL write_subbasin_assessment(resdir,nsubCrit,basincrit,iens,.FALSE.)
        ENDIF

!HYSS: Save and close files
        DO iens = 1,myAssimData%info%nstatout+1
          CALL close_outputfiles(nsub,naquifers,iens)
        ENDDO
        CALL close_observations(modeldir)
        iens = 1    !reset
        
!HYSS: Write results to map-files
        CALL save_mapfiles(resdir,nsub,nmapperiod,iens,.FALSE.,.FALSE.)

!ASSIMILATION: some deallocation, closing of files, cleaning up, and final log messages?
 
          call cpu_time(stop_time)
          total_time(17)=total_time(17)+stop_time-start_time
          start_time=stop_time

       
  ENDIF !doassimilation
         write(6,*)'---------------------'
         write(6,*)'final timing'
         write(6,*)'---------------------'
         write(6,*)'total_time(1:20)'
         write(6,*)total_time
         write(6,*)'total_time_stat(1:4)'
         write(6,*)total_time_stat

!=====================================================================
!End of Data Assimilation section
!=====================================================================
#endif
  
  
!>Deallocate variables
      CALL deallocate_worldvar()
      CALL deallocate_modvar(nsub)
      CALL deallocate_model_states(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
      IF(doopt.AND.optim%task_MC) CALL deallocate_MCvariables()
      IF(ALLOCATED(basincrit)) DEALLOCATE(basincrit)
      IF(ALLOCATED(simperformance)) DEALLOCATE(simperformance)
      IF(ALLOCATED(par)) DEALLOCATE(par)

!Write stop time on log-file
      CALL DATE_AND_TIME (values=datim)
      WRITE(6,*)
      WRITE(6,*) '---------------------------------------------------'
      WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
               ' Job finished date: ',datim(1),'-',datim(2),'-',datim(3),    &
               '  time: ',datim(5),':',datim(6),':',datim(7)

      CLOSE(6)
      STOP 84
      END PROGRAM
