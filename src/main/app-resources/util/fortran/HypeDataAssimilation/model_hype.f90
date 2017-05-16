!> \file model_hype.f90
!> Contains module modelmodule.

!>Main module for hydrological model HYPE.
!>
!> Detailed information of HYPE is found in the <a href=../HYPEmodelDescription.odt>HYPE model description</a>
MODULE MODELMODULE
!The HYPE model (HYdrological Predictions for the Environment)

!Copyright 2011-2017 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!----------------------------------------------------------------------

!Used modules
  USE STATETYPE_MODULE  !model state variable types
  USE MODVAR            !model variables, general model
  USE HYPEVARIABLES     !model variables, HYPE model
  USE HYPE_WATERBALANCE !water balance variables and indices
!Subroutines also uses GENERAL_WATER_CONCENTRATION, GLACIER_SOILMODEL, SOILMODEL_DEFAULT, 
!SOIL_PROCESSES, NPC_SOIL_PROCESSES, SURFACEWATER_PROCESSES, NPC_SURFACEWATER_PROCESSES, 
!IRRIGATION_MODULE, REGIONAL_GROUNDWATER_MODULE

  IMPLICIT NONE
  PRIVATE

  ! Private subroutines
  !--------------------------------------------------------------------
  ! apply_quseobs 
  ! apply_wendupd
  ! apply_nutrientcorr
  ! apply_qarupd
  !--------------------------------------------------------------------
  PUBLIC :: model_version_information, &
            initiate_output_variables, &
            initiate_model_parameters, &
            initiate_model_state, &
            initiate_model,  &
            set_special_models, &
            get_special_model_parameters,  &
            calculate_special_model_parameters,  &
            model, &
            load_modeldefined_input, &
            reload_modeldefined_observations, &
            open_modeldefined_outputfiles, &
            close_modeldefined_outputfiles

  CONTAINS
  
  !>Information about the model version to print in log-file
  !-----------------------------------------------------------
  SUBROUTINE model_version_information(funit)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: funit !<fileunit for log-file
    
    WRITE(funit,600) '----------------------'
    WRITE(funit,600) ' HYPE version 4.13.x  '
    WRITE(funit,600) '----------------------'
600 FORMAT(A22)
    
  END SUBROUTINE model_version_information
  
  !>Set variables holding output information from the HYPE model; 
  !!outvarid and loadheadings
  !-----------------------------------------------------------------
  SUBROUTINE initiate_output_variables()
  
    max_outvar = 356     !number of output variables
    ALLOCATE(outvarid(max_outvar))
    
    outvarid(o_crun) = outvaridtype('crun','computed runoff', i_sum , 0,2,'mm' ,'mm' ,'mapCRUN.txt' ,'timeCRUN.txt')
    outvarid(o_rrun) = outvaridtype('rrun','recorded runoff', i_sum , 0,0,'mm' ,'mm' ,'mapRRUN.txt' ,'timeRRUN.txt')
    outvarid(o_prec) = outvaridtype('prec','precipitation'  , i_sum , 0,1,'mm' ,'mm' ,'mapPREC.txt' ,'timePREC.txt')
    outvarid(o_tobs)  = outvaridtype('temp','air temperature'   ,i_mean , 0,1,'deg'  ,'degree Celsius','mapTEMP.txt' ,'timeTEMP.txt')
    outvarid(o_crunT1)  = outvaridtype('coT1','computed T1 runoff',i_wmean, o_crun,0,'?','?','mapCOT1.txt','timeCOT1.txt')
    outvarid(o_crunT2)  = outvaridtype('coT2','computed T2 runoff',i_wmean, o_crun,0,'deg' ,'degree Celsius' ,'mapCOT2.txt','timeCOT2.txt')
    outvarid(o_roum) = outvaridtype('roum','rec. outflow main',i_mean, 0,0,'m3/s','m3/s','mapROUM.txt','timeROUM.txt')
    outvarid(o_cprecT1) = outvaridtype('cpT1','recorded T1 precip',i_wmean, o_prec,0,'?','?','mapCPT1.txt'  ,'timeCPT1.txt')
    outvarid(9)  = outvaridtype('ceT1','computed T1 evap'  ,i_wmean, o_evap,0,'?' ,'?' ,'mapCET1.txt'  ,'timeCET1.txt')
    outvarid(o_soim) = outvaridtype('soim','computed soil water',i_mean, 0,2,'mm' ,'mm','mapSOIM.txt' ,'timeSOIM.txt')
    outvarid(o_csoilT1) = outvaridtype('csT1','computed T1 soil'  ,i_wmean, o_soim,0,'?' ,'?' ,'mapCST1.txt'  ,'timeCST1.txt')
    outvarid(o_csoilT2) = outvaridtype('csT2','computed T2 soil'  ,i_wmean, o_soim,0,'deg','degree Celsius' ,'mapCST2.txt'  ,'timeCST2.txt')
    outvarid(o_roub) = outvaridtype('roub','rec. outflow branch',i_mean, 0,0,'m3/s','m3/s','mapROUB.txt','timeROUB.txt')
    outvarid(o_snow) = outvaridtype('snow','snow pack'         ,i_mean , 0,2,'mm' ,'mm water','mapSNOW.txt' ,'timeSNOW.txt')
    outvarid(o_evap) = outvaridtype('evap','subbasin evaporation',i_sum , 0,1,'mm'   ,'mm' ,'mapEVAP.txt' ,'timeEVAP.txt')
    outvarid(16) = outvaridtype('reT1','recorded T1 outflow olake',i_wmean, o_rout,0,'?' ,'?' ,'mapRET1.txt'  ,'timeRET1.txt')
    outvarid(17) = outvaridtype('reT2','recorded T2 outflow olake',i_wmean, o_rout,0,'deg'  ,'degree Celsius','mapRET2.txt'  ,'timeRET2.txt') 
    outvarid(o_snowdens) = outvaridtype('sden','computed snow density',i_mean, 0,0,'g/cm3','g/cm3','mapSDEN.txt'  ,'timeSDEN.txt')
    outvarid(o_grwlevel) = outvaridtype('gwat','computed grw level',i_mean , 0,0,'m'    ,'m'             ,'mapGWAT.txt'  ,'timeGWAT.txt')
    outvarid(o_crunIN) = outvaridtype('coIN','computed IN runoff',i_wmean, o_crun,0,'ug/L' ,'ug Inorg-N/L'  ,'mapCOIN.txt','timeCOIN.txt')
    outvarid(o_crunON) = outvaridtype('coON','computed ON runoff',i_wmean, o_crun,0,'ug/L' ,'ug Org-N/L'    ,'mapCOON.txt','timeCOON.txt')
    outvarid(o_crunSP) = outvaridtype('coSP','computed SP runoff',i_wmean, o_crun,0,'ug/L' ,'ug SRP-P/L'    ,'mapCOSP.txt','timeCOSP.txt')
    outvarid(o_crunPP) = outvaridtype('coPP','computed PP runoff',i_wmean, o_crun,0,'ug/L' ,'ug PartP-P/L'  ,'mapCOPP.txt','timeCOPP.txt')
    outvarid(24) = outvaridtype('reIN','recorded IN outflow olake',i_wmean, o_rout,0,'ug/L' ,'ug INORG-N/L'  ,'mapREIN.txt','timeREIN.txt')
    outvarid(25) = outvaridtype('reON','recorded ON outflow olake',i_wmean, o_rout,0,'ug/L' ,'ug ORG-N/L'    ,'mapREON.txt','timeREON.txt')
    outvarid(26) = outvaridtype('reSP','recorded SP outflow olake',i_wmean, o_rout,0,'ug/L' ,'ug SRP-P/L'    ,'mapRESP.txt','timeRESP.txt')
    outvarid(27) = outvaridtype('rePP','recorded PP outflow olake',i_wmean, o_rout,0,'ug/L' ,'ug PartP-P/L'  ,'mapREPP.txt','timeREPP.txt')
    outvarid(o_reTN) = outvaridtype('reTN','recorded TN outflow olake',i_wmean, o_rout,0,'ug/L' ,'ug Tot-N/L'    ,'mapRETN.txt','timeRETN.txt')
    outvarid(o_reTP) = outvaridtype('reTP','recorded TP outflow olake',i_wmean, o_rout,0,'ug/L' ,'ug Tot-P/L'    ,'mapRETP.txt','timeRETP.txt')
    outvarid(30) = outvaridtype('qerr','outflow error',i_mean , 0,0,'m3/s','m3/s' ,'mapQERR.txt','timeQERR.txt')
    outvarid(31) = outvaridtype('cobc','computed outflow before correction',i_mean , 0,0,'m3/s','m3/s' ,'mapCOBC.txt','timeCOBC.txt')
    outvarid(32) = outvaridtype('wtmp','computed water temp' ,i_mean , 0,0,'deg'  ,'degree Celsius','mapWTMP.txt' ,'timeWTMP.txt')
    outvarid(33) = outvaridtype('werr','waterstage error',i_mean , 0,0,'m','m' ,'mapWERR.txt','timeWERR.txt')
    outvarid(34) = outvaridtype('cwbc','computed wcom before correction',i_mean , 0,0,'m','m' ,'mapCWBC.txt','timeCWBC.txt')
    outvarid(35) = outvaridtype('reOC','recorded TOC outflow',i_wmean, o_rout,0,'mg/L','mg org-C/L'    ,'mapREOC.txt','timeREOC.txt')
    outvarid(o_csoilIN) = outvaridtype('csIN','computed IN soil'    ,i_wmean, o_soim,0,'ug/L','ug INORG-N/L'  ,'mapCSIN.txt','timeCSIN.txt')
    outvarid(o_soilfrost) = outvaridtype('sfst','computed soil frost' ,i_mean , 0,0, 'cm'  ,'cm'            ,'mapSFST.txt','timeSFST.txt')
    outvarid(o_soiltmp) = outvaridtype('stmp','computed soil temp'  ,i_mean , 0,0, 'deg' ,'degree Celcius','mapSTMP.txt','timeSTMP.txt')
    outvarid(o_snowdepth) = outvaridtype('sdep','computed snow depth' ,i_mean , 0,0, 'cm'  ,'cm'            ,'mapSDEP.txt','timeSDEP.txt')
    outvarid(o_epot) = outvaridtype('epot','potential evap'      ,i_sum  , 0,1, 'mm'  ,'mm'            ,'mapEPOT.txt','timeEPOT.txt')
    outvarid(o_reepot) = outvaridtype('repo','recorded pot. evap',i_sum ,0,0,'mm' ,'mm'   ,'mapREPO.txt','timeREPO.txt')
    outvarid(42) = outvaridtype('eobs','recorded evaporation'       ,i_sum , 0,0,'mm'   ,'mm'    ,'mapEOBS.txt' ,'timeEOBS.txt')
    outvarid(43) = outvaridtype('cprc','corr precipitation'     ,i_sum , 0,1,'mm'   ,'mm'            ,'mapCPRC.txt' ,'timeCPRC.txt')
    outvarid(o_crunOC) = outvaridtype('coOC','computed TOC runoff',i_wmean, o_crun,0,'mg/L'    ,'mg org-C/L'     ,'mapCOOC.txt','timeCOOC.txt')
    outvarid(o_csoilOC) = outvaridtype('csOC','computed TOC soil'  ,i_wmean, o_soim,0,'mg/L'    ,'mg org-C/L'    ,'mapCSOC.txt','timeCSOC.txt')
    outvarid(46) = outvaridtype('ccOC','computed TOC outflow',i_wmean, o_cout,0,'mg/L','mg org-C/L' ,'mapCCOC.txt','timeCCOC.txt')
    outvarid(47) = outvaridtype('phC1','pool humusC soil1',i_mean , 0,0,'kg/km2' ,'kg org-C/km2'    ,'mapPHC1.txt','timePHC1.txt')
    outvarid(48) = outvaridtype('pfC1','pool fastC soil1',i_mean , 0,0,'kg/km2' ,'kg org-C/km2'    ,'mapPFC1.txt','timePFC1.txt')
    outvarid(49) = outvaridtype('phC2','pool humusC soil2',i_mean , 0,0,'kg/km2' ,'kg org-C/km2'    ,'mapPHC2.txt','timePHC2.txt')
    outvarid(50) = outvaridtype('pfC2','pool fastC soil2',i_mean , 0,0,'kg/km2' ,'kg org-C/km2'    ,'mapPFC2.txt','timePFC2.txt')
    outvarid(51) = outvaridtype('wcom','olake water stage',i_mean , 0,0,'m','meter' ,'mapWCOM.txt','timeWCOM.txt')
    outvarid(o_rewstr) = outvaridtype('wstr','rec olake water st',i_mean , 0,0,'m','meter' ,'mapWSTR.txt','timeWSTR.txt')
    outvarid(o_cout) = outvaridtype('cout','comp outflow olake',i_mean , 0,0,'m3/s','m3/s' ,'mapCOUT.txt','timeCOUT.txt')
    outvarid(o_rout) = outvaridtype('rout','rec outflow olake',i_mean , 0,0,'m3/s','m3/s' ,'mapROUT.txt','timeROUT.txt')
    outvarid(55) = outvaridtype('ccIN','comp conc IN olake',i_wmean, o_cout,0,'ug/L','ug InorgN-N/L' ,'mapCCIN.txt','timeCCIN.txt')
    outvarid(56) = outvaridtype('ccON','comp conc ON olake',i_wmean, o_cout,0,'ug/L','ug OrgN-N/L' ,'mapCCON.txt','timeCCON.txt')
    outvarid(57) = outvaridtype('ccSP','comp conc SP olake',i_wmean, o_cout,0,'ug/L','ug SRP-P/L' ,'mapCCSP.txt','timeCCSP.txt')
    outvarid(58) = outvaridtype('ccPP','comp conc PP olake',i_wmean, o_cout,0,'ug/L','ug PartP-P/L' ,'mapCCPP.txt','timeCCPP.txt')
    outvarid(59) = outvaridtype('rsnw','recorded snow depth',i_mean , 0,0,'cm','cm' ,'mapRSNW.txt','timeRSNW.txt')
    outvarid(60) = outvaridtype('resf','recorded soil frost',i_mean , 0,0,'cm','cm' ,'mapRESF.txt','timeRESF.txt')
    outvarid(61) = outvaridtype('regw','recorded grw level',i_mean , 0,0,'m','meter' ,'mapREGW.txt','timeREGW.txt')
    outvarid(63) = outvaridtype('ccT1','comp conc T1 olake',i_wmean, o_cout,0,'?','?' ,'mapCCT1.txt','timeCCT1.txt')
    outvarid(64) = outvaridtype('ccT2','comp conc T2 olake',i_wmean, o_cout,0,'deg'  ,'degree Celsius','mapCCT2.txt','timeCCT2.txt')
    outvarid(o_upte) = outvaridtype('upte','aver upstr temperature',i_mean, 0,0,'deg','degree Celsius' ,'mapUPTE.txt','timeUPTE.txt')
    outvarid(o_crunTN) = outvaridtype('coTN','computed TN runoff',i_wmean, o_crun,0,'ug/L' ,'ug Tot-N/L'    ,'mapCOTN.txt','timeCOTN.txt')
    outvarid(o_crunTP) = outvaridtype('coTP','computed TP runoff',i_wmean, o_crun,0,'ug/L' ,'ug Tot-P/L'    ,'mapCOTP.txt','timeCOTP.txt')
    outvarid(o_pfN(1)) = outvaridtype('pfN1','pool fastN soil1' ,i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPFN1.txt','timePFN1.txt')
    outvarid(o_pfN(2)) = outvaridtype('pfN2','pool fastN soil2' ,i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPFN2.txt','timePFN2.txt')
    outvarid(o_pfN(3)) = outvaridtype('pfN3','pool fastN soil3' ,i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPFN3.txt','timePFN3.txt')
    outvarid(o_phN(1)) = outvaridtype('phN1','pool humusN soil1',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPHN1.txt','timePHN1.txt')
    outvarid(o_phN(2)) = outvaridtype('phN2','pool humusN soil2',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPHN2.txt','timePHN2.txt')
    outvarid(o_phN(3)) = outvaridtype('phN3','pool humusN soil3',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPHN3.txt','timePHN3.txt')
    outvarid(74) = outvaridtype('pIN1','pool InorgN soil1',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPIN1.txt','timePIN1.txt')
    outvarid(75) = outvaridtype('pIN2','pool InorgN soil2',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPIN2.txt','timePIN2.txt')
    outvarid(76) = outvaridtype('pIN3','pool InorgN soil3',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPIN3.txt','timePIN3.txt')
    outvarid(77) = outvaridtype('ccTN','computed TN olake',i_wmean, o_cout,0,'ug/L' ,'ug Tot-N/L'    ,'mapCCTN.txt','timeCCTN.txt')
    outvarid(78) = outvaridtype('ccTP','computed TP olake',i_wmean, o_cout,0,'ug/L' ,'ug Tot-P/L'    ,'mapCCTP.txt','timeCCTP.txt')
    outvarid(o_ro1) = outvaridtype('cro1','computed runoff 1',i_sum  , 0,0,'mm'     ,'mm'           ,'mapCRO1.txt','timeCRO1.txt')
    outvarid(o_ro2) = outvaridtype('cro2','computed runoff 2',i_sum  , 0,0,'mm'     ,'mm'           ,'mapCRO2.txt','timeCRO2.txt')
    outvarid(81) = outvaridtype('cgwl','computed groundwater loss',i_mean ,0, 0,'m3/s' ,'m3/s'  ,'mapCGWL.txt' ,'timeCGWL.txt')
    outvarid(o_rod) = outvaridtype('crod','computed rf drain',i_sum  , 0,0,'mm' ,'mm'           ,'mapCROD.txt' ,'timeCROD.txt')
    outvarid(o_ros) = outvaridtype('cros','computed surface rf'    ,i_sum ,0, 0,'mm' ,'mm'      ,'mapCROS.txt' ,'timeCROS.txt')
    outvarid(o_soildenitr) = outvaridtype('deni','denitrifikation'  ,i_sum , 0,0,'kg/km2' ,'kg N/km2'  ,'mapDENI.txt' ,'timeDENI.txt')
    outvarid(o_cropNupt) = outvaridtype('crut','crop uptake'      ,i_sum , 0,0,'kg/km2' ,'kg N/km2'  ,'mapCRUT.txt' ,'timeCRUT.txt')
    outvarid(o_degrfN) = outvaridtype('faIN','fast to inorganic',i_sum , 0,0,'kg/km2' ,'kg N/km2'  ,'mapFAIN.txt' ,'timeFAIN.txt')
    outvarid(o_soilNatm) = outvaridtype('atmd','atm deposition N' ,i_sum , 0,0,'kg/km2' ,'kg N/km2'  ,'mapATMD.txt' ,'timeATMD.txt')
    outvarid(o_ppP(1)) = outvaridtype('ppP1','pool partP soil1',i_mean , 0,0,'kg/km2' ,'kg P/km2'    ,'mapPPP1.txt','timePPP1.txt')   
    outvarid(o_ppP(2)) = outvaridtype('ppP2','pool partP soil2',i_mean , 0,0,'kg/km2' ,'kg P/km2'    ,'mapPPP2.txt','timePPP2.txt')
    outvarid(o_ppP(3)) = outvaridtype('ppP3','pool partP soil3',i_mean , 0,0,'kg/km2' ,'kg P/km2'    ,'mapPPP3.txt','timePPP3.txt')
    outvarid(91) = outvaridtype('pSP1','pool SRP soil1',  i_mean , 0,0,'kg/km2' ,'kg P/km2'    ,'mapPSP1.txt','timePSP1.txt')
    outvarid(92) = outvaridtype('pSP2','pool SRP soil2',  i_mean , 0,0,'kg/km2' ,'kg P/km2'    ,'mapPSP2.txt','timePSP2.txt')
    outvarid(93) = outvaridtype('pSP3','pool SRP soil3',  i_mean , 0,0,'kg/km2' ,'kg P/km2'    ,'mapPSP3.txt','timePSP3.txt')
    outvarid(o_pfP(1)) = outvaridtype('pfP1','pool fastP soil1',i_mean , 0,0,'kg/km2' ,'kg P/km2'    ,'mapPFP1.txt','timePFP1.txt')
    outvarid(o_pfP(2)) = outvaridtype('pfP2','pool fastP soil2',i_mean , 0,0,'kg/km2' ,'kg P/km2'    ,'mapPFP2.txt','timePFP2.txt')
    outvarid(o_pfP(3)) = outvaridtype('pfP3','pool fastP soil3',i_mean , 0,0,'kg/km2' ,'kg P/km2'    ,'mapPFP3.txt','timePFP3.txt')
    outvarid(o_cprecIN) = outvaridtype('cpIN','recorded IN precip',i_wmean, o_prec,0,'ug/L','ug N/L','mapCPIN.txt'  ,'timeCPIN.txt')
    outvarid(o_cprecSP) = outvaridtype('cpSP','recorded SP precip',i_wmean, o_prec,0,'ug/L','ug P/L','mapCPSP.txt'  ,'timeCPSP.txt')
    outvarid(o_reswe) = outvaridtype('rswe','rec snow water eq.',i_mean , 0,0,'mm','mm' ,'mapRSWE.txt','timeRSWE.txt')
    outvarid(100) = outvaridtype('acdf','accumulated volume error',i_sum , 0,0,'mm','mm' ,'mapACDF.txt','timeACDF.txt')
    outvarid(o_cloc) = outvaridtype('cloc','comp local flow',   i_mean , 0,0, 'm3/s','m3/s' ,'mapCLOC.txt','timeCLOC.txt')
    outvarid(102) = outvaridtype('clIN','comp local conc IN',i_wmean, o_cloc,0,'ug/L','ug InorgN-N/L' ,'mapCLIN.txt','timeCLIN.txt')
    outvarid(103) = outvaridtype('clON','comp local conc ON',i_wmean, o_cloc,0,'ug/L','ug OrgN-N/L' ,'mapCLON.txt','timeCLON.txt')
    outvarid(104) = outvaridtype('clSP','comp local conc SP',i_wmean, o_cloc,0,'ug/L','ug SRP-P/L' ,'mapCLSP.txt','timeCLSP.txt')
    outvarid(105) = outvaridtype('clPP','comp local conc PP',i_wmean, o_cloc,0,'ug/L','ug PartP-P/L' ,'mapCLPP.txt','timeCLPP.txt')
    outvarid(106) = outvaridtype('clTN','comp local conc TN',i_wmean, o_cloc,0,'ug/L' ,'ug Tot-N/L'    ,'mapCLTN.txt','timeCLTN.txt')
    outvarid(107) = outvaridtype('clTP','comp local conc TP',i_wmean, o_cloc,0,'ug/L' ,'ug Tot-P/L'    ,'mapCLTP.txt','timeCLTP.txt')
    outvarid(108) = outvaridtype('pfC3','pool fastC soil3',  i_mean , 0,0,'kg/km2' ,'kg org-C/km2'    ,'mapPFC3.txt','timePFC3.txt')
    outvarid(109) = outvaridtype('phC3','pool humusC soil3', i_mean , 0,0,'kg/km2' ,'kg org-C/km2'    ,'mapPHC3.txt','timePHC3.txt')
    outvarid(110) = outvaridtype('totN','comp load TN olake',i_sum , 0,0,'kg' ,'kg totN'    ,'mapTOTN.txt','timeTOTN.txt')
    outvarid(111) = outvaridtype('totP','comp load TP olake',i_sum , 0,0,'kg' ,'kg totP'    ,'mapTOTP.txt','timeTOTP.txt')
    outvarid(112) = outvaridtype('cinf','comp inflow to olake',i_mean , 0,0,'m3/s' ,'m3/s'    ,'mapCINF.txt','timeCINF.txt')
    outvarid(113) = outvaridtype('rinf','rec inflow to olake',i_mean , 0,0,'m3/s' ,'m3/s'    ,'mapRINF.txt','timeRINF.txt')
    outvarid(114) = outvaridtype('clrv','comp local river volume',i_mean , 0,0,'m3' ,'m3'    ,'mapCLRV.txt','timeCLRV.txt')
    outvarid(115) = outvaridtype('cmrv','comp main river volume', i_mean , 0,0,'m3' ,'m3'    ,'mapCMRV.txt','timeCMRV.txt')
    outvarid(o_soilPatm) = outvaridtype('atmp','atm deposition TP', i_sum , 0,0,'kg/km2' ,'kg P/km2'  ,'mapATMP.txt' ,'timeATMP.txt')
    outvarid(117) = outvaridtype('glcv','comp glacier ice volume', i_mean , 0,0,'km3' ,'km3'    ,'mapGLCV.txt','timeGLCV.txt')
    outvarid(118) = outvaridtype('glca','comp glacier area', i_mean , 0,0,'km2' ,'km2'    ,'mapGLCA.txt','timeGLCA.txt')
    outvarid(119) = outvaridtype('rtoN','rec load TN olake', i_sum , 0,0,'kg' ,'kg totN'    ,'mapRTON.txt','timeRTON.txt')
    outvarid(120) = outvaridtype('rtoP','rec load TN olake', i_sum , 0,0,'kg' ,'kg totP'    ,'mapRTOP.txt','timeRTOP.txt')
    outvarid(o_ctmp) = outvaridtype('ctmp','corrected air temperature'   ,i_mean , 0,1,'deg'  ,'degree Celsius','mapCTMP.txt' ,'timeCTMP.txt') 
    outvarid(122) = outvaridtype('irel','irrigation evap losses'   ,i_sum , 0,0,'m3'  ,'m3','mapIREL.txt' ,'timeIREL.txt') 
    outvarid(o_coum) = outvaridtype('coum','comp outflow main',i_mean , 0,0,'m3/s','m3/s' ,'mapCOUM.txt','timeCOUM.txt')
    outvarid(124) = outvaridtype('irld','abstr local dam f irr'    ,i_sum , 0,0,'m3'  ,'m3','mapIRLD.txt' ,'timeIRLD.txt')
    outvarid(125) = outvaridtype('irlr','abstr local river f irr'  ,i_sum , 0,0,'m3'  ,'m3','mapIRLR.txt' ,'timeIRLR.txt')
    outvarid(o_applirr) = outvaridtype('irra','applied irrigation water' ,i_sum , 0,0,'m3'  ,'m3','mapIRRA.txt' ,'timeIRRA.txt')
    outvarid(127) = outvaridtype('irrg','gwater abstracted f irr'  ,i_sum , 0,0,'m3'  ,'m3','mapIRRG.txt' ,'timeIRRG.txt') 
    outvarid(128) = outvaridtype('irrs','irr abstr f other subbasins'   ,i_sum , 0,0,'m3'  ,'m3','mapIRRS.txt' ,'timeIRRS.txt')
    outvarid(o_upsd) = outvaridtype('upsd','aver upstream soildf',i_mean, 0,0,'mm','mm','mapUPSD.txt','timeUPSD.txt')
    outvarid(o_coub) = outvaridtype('coub','comp outflow branch',i_mean , 0,0,'m3/s','m3/s' ,'mapCOUB.txt','timeCOUB.txt')
    outvarid(o_soildef) = outvaridtype('smdf','soil moisture deficit'    ,i_mean , 0,2,'mm'   ,'mm'            ,'mapSMDF.txt' ,'timeSMDF.txt')
    outvarid(o_phP(1)) = outvaridtype('phP1','pool humusP soil1',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPHP1.txt','timePHP1.txt')
    outvarid(o_phP(2)) = outvaridtype('phP2','pool humusP soil2',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPHP2.txt','timePHP2.txt')
    outvarid(o_phP(3)) = outvaridtype('phP3','pool humusP soil3',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPHP3.txt','timePHP3.txt')
    outvarid(135) = outvaridtype('totC','comp load OC olake',i_sum , 0,0,'kg' ,'kg OC'    ,'mapTOTC.txt','timeTOTC.txt')
    outvarid(136) = outvaridtype('pON1','pool orgN soil1',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPON1.txt','timePON1.txt')
    outvarid(137) = outvaridtype('pON2','pool orgN soil2',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPON2.txt','timePON2.txt')
    outvarid(138) = outvaridtype('pON3','pool orgN soil3',i_mean , 0,0,'kg/km2' ,'kg N/km2'    ,'mapPON3.txt','timePON3.txt')
    outvarid(o_sltmp(1)) = outvaridtype('stm1','computed soillayer 1 temp'  ,i_mean , 0,0, 'deg' ,'degree Celcius','mapSTM1.txt','timeSTM1.txt')
    outvarid(o_sltmp(2)) = outvaridtype('stm2','computed soillayer 2 temp'  ,i_mean , 0,0, 'deg' ,'degree Celcius','mapSTM2.txt','timeSTM2.txt')
    outvarid(o_sltmp(3)) = outvaridtype('stm3','computed soillayer 3 temp'  ,i_mean , 0,0, 'deg' ,'degree Celcius','mapSTM3.txt','timeSTM3.txt')
    outvarid(o_rainfall) = outvaridtype('cpRF','precipitation as rain'     ,i_sum , 0,1,'mm'   ,'mm'             ,'mapCPRF.txt' ,'timeCPRF.txt')
    outvarid(o_snowfall) = outvaridtype('cpSF','precipitation as snow'     ,i_sum , 0,1,'mm'   ,'mm'             ,'mapCPSF.txt' ,'timeCPSF.txt')
    outvarid(144) = outvaridtype('colv','volume of olake/w',i_mean , 0,0,'M(m3)' ,'M(m3)','mapCOLV.txt','timeCOLV.txt') !Computed Lake Volume  (basins summed to outlet if any) 
    outvarid(145) = outvaridtype('cilv','volume of ilake',i_mean , 0,0,'M(m3)' ,'M(m3)'  ,'mapCILV.txt','timeCILV.txt') !Computed ILake Volume 
    outvarid(146) = outvaridtype('clbv','volume of olake/lb',i_mean , 0,0,'M(m3)' ,'M(m3)','mapCLBV.txt','timeCLBV.txt') !Computed Olake Volume Computed (volumes for individual basins if any)
    outvarid(o_ro3) = outvaridtype('cro3','computed runoff 3',i_sum  , 0,0,'mm'     ,'mm'   ,'mapCRO3.txt','timeCRO3.txt')
    !------------------------------------------------------------------------------------
    ! New Lake and River Ice, Snow, and Water Temperature model variables David(20130426)
    !------------------------------------------------------------------------------------
    !Lake ice and snow depth variables
    outvarid(148) = outvaridtype('coli','comp olake ice depth'      ,i_mean , 0,0,'cm'  ,'cm','mapCOLI.txt' ,'timeCOLI.txt')
    outvarid(149) = outvaridtype('cili','comp ilake ice depth'      ,i_mean , 0,0,'cm'  ,'cm','mapCILI.txt' ,'timeCILI.txt')
    outvarid(150) = outvaridtype('colb','comp olake blackice depth' ,i_mean , 0,0,'cm'  ,'cm','mapCOLB.txt' ,'timeCOLB.txt')
    outvarid(151) = outvaridtype('cilb','comp ilake blackice depth' ,i_mean , 0,0,'cm'  ,'cm','mapCILB.txt' ,'timeCILB.txt')
    outvarid(152) = outvaridtype('cols','comp olake snow depth'     ,i_mean , 0,0,'cm'  ,'cm','mapCOLS.txt' ,'timeCOLS.txt')
    outvarid(153) = outvaridtype('cils','comp ilake snow depth'     ,i_mean , 0,0,'cm'  ,'cm','mapCILS.txt' ,'timeCILS.txt')
    outvarid(154) = outvaridtype('roli','rec. olake ice depth'      ,i_mean , 0,0,'cm'  ,'cm','mapROLI.txt' ,'timeROLI.txt')
    outvarid(155) = outvaridtype('rili','rec. ilake ice depth'      ,i_mean , 0,0,'cm'  ,'cm','mapRILI.txt' ,'timeRILI.txt')
    outvarid(156) = outvaridtype('rolb','rec. olake blackice depth' ,i_mean , 0,0,'cm'  ,'cm','mapROLB.txt' ,'timeROLB.txt')
    outvarid(157) = outvaridtype('rilb','rec. ilake blackice depth' ,i_mean , 0,0,'cm'  ,'cm','mapRILB.txt' ,'timeRILB.txt')
    outvarid(158) = outvaridtype('rols','rec. olake snow depth'     ,i_mean , 0,0,'cm'  ,'cm','mapROLS.txt' ,'timeROLS.txt')
    outvarid(159) = outvaridtype('rils','rec. ilake snow depth'     ,i_mean , 0,0,'cm'  ,'cm','mapRILS.txt' ,'timeRILS.txt')
    outvarid(160) = outvaridtype('cmri','comp main river ice depth'      ,i_mean , 0,0,'cm'  ,'cm','mapCMRI.txt' ,'timeCMRI.txt')
    outvarid(161) = outvaridtype('clri','comp local river ice depth'     ,i_mean , 0,0,'cm'  ,'cm','mapCLRI.txt' ,'timeCLRI.txt')
    outvarid(162) = outvaridtype('cmrb','comp main river blackice depth' ,i_mean , 0,0,'cm'  ,'cm','mapCMRB.txt' ,'timeCMRB.txt')
    outvarid(163) = outvaridtype('clrb','comp local river blackice depth',i_mean , 0,0,'cm'  ,'cm','mapCLRB.txt' ,'timeCLRB.txt')
    outvarid(164) = outvaridtype('cmrs','comp main river snow depth'     ,i_mean , 0,0,'cm'  ,'cm','mapCMRS.txt' ,'timeCMRS.txt')
    outvarid(165) = outvaridtype('clrs','comp local river snow depth'    ,i_mean , 0,0,'cm'  ,'cm','mapCLRS.txt' ,'timeCLRS.txt')
    outvarid(166) = outvaridtype('rmri','rec. main river ice depth'      ,i_mean , 0,0,'cm'  ,'cm','mapRMRI.txt' ,'timeRMRI.txt')
    outvarid(167) = outvaridtype('rlri','rec. local river ice depth'     ,i_mean , 0,0,'cm'  ,'cm','mapRLRI.txt' ,'timeRLRI.txt')
    outvarid(168) = outvaridtype('rmrb','rec. main river blackice depth' ,i_mean , 0,0,'cm'  ,'cm','mapRMRB.txt' ,'timeRMRB.txt')
    outvarid(169) = outvaridtype('rlrb','rec. local river blackice depth',i_mean , 0,0,'cm'  ,'cm','mapRLRB.txt' ,'timeRLRB.txt')
    outvarid(170) = outvaridtype('rmrs','rec. main river snow depth'     ,i_mean , 0,0,'cm'  ,'cm','mapRMRS.txt' ,'timeRMRS.txt')
    outvarid(171) = outvaridtype('rlrs','rec. local river snow depth'    ,i_mean , 0,0,'cm'  ,'cm','mapRLRS.txt' ,'timeRLRS.txt')
    outvarid(172) = outvaridtype('olst','comp olake surface temp'         ,i_mean , 0,0,'deg'  ,'degree Celsius','mapOLST.txt' ,'timeOLST.txt')
    outvarid(173) = outvaridtype('olut','comp olake upper temp'           ,i_mean , 0,0,'deg'  ,'degree Celsius','mapOLUT.txt' ,'timeOLUT.txt')
    outvarid(174) = outvaridtype('ollt','comp olake lower temp'           ,i_mean , 0,0,'deg'  ,'degree Celsius','mapOLLT.txt' ,'timeOLLT.txt')
    outvarid(175) = outvaridtype('olwt','comp olake mean  temp'           ,i_mean , 0,0,'deg'  ,'degree Celsius','mapOLWT.txt' ,'timeOLWT.txt')
    outvarid(176) = outvaridtype('ilst','comp ilake surface temp'         ,i_mean , 0,0,'deg'  ,'degree Celsius','mapILST.txt' ,'timeILST.txt')
    outvarid(177) = outvaridtype('ilwt','comp ilake mean  temp'           ,i_mean , 0,0,'deg'  ,'degree Celsius','mapILWT.txt' ,'timeILWT.txt')
    outvarid(178) = outvaridtype('lrst','comp local river surface temp'   ,i_mean , 0,0,'deg'  ,'degree Celsius','mapLRST.txt' ,'timeLRST.txt')
    outvarid(179) = outvaridtype('lrwt','comp local river mean  temp'     ,i_mean , 0,0,'deg'  ,'degree Celsius','mapLRWT.txt' ,'timeLRWT.txt')
    outvarid(180) = outvaridtype('mrst','comp main  river surface temp'   ,i_mean , 0,0,'deg'  ,'degree Celsius','mapMRST.txt' ,'timeMRST.txt')
    outvarid(181) = outvaridtype('mrwt','comp main  river mean  temp'     ,i_mean , 0,0,'deg'  ,'degree Celsius','mapMRWT.txt' ,'timeMRWT.txt')
    
    !MODIS Water surface temperatures
    outvarid(182) = outvaridtype('rolt','rec. olake surface temp'         ,i_mean , 0,0,'deg'  ,'degree Celsius','mapROLT.txt' ,'timeROLT.txt')
    outvarid(183) = outvaridtype('rilt','rec. ilake surface temp'         ,i_mean , 0,0,'deg'  ,'degree Celsius','mapRILT.txt' ,'timeRILT.txt')
    outvarid(184) = outvaridtype('rmrt','rec. main river surface temp'    ,i_mean , 0,0,'deg'  ,'degree Celsius','mapRMRT.txt' ,'timeRMRT.txt')
    outvarid(185) = outvaridtype('mrto','OLD comp main river mean temp'   ,i_mean , 0,0,'deg'  ,'degree Celsius','mapMRTO.txt' ,'timeMRTO.txt')
    outvarid(186) = outvaridtype('lrto','OLD comp local river mean temp'  ,i_mean , 0,0,'deg'  ,'degree Celsius','mapLRTO.txt' ,'timeLRTO.txt')
    outvarid(187) = outvaridtype('ilto','OLD comp ilake mean temp'        ,i_mean , 0,0,'deg'  ,'degree Celsius','mapILTO.txt' ,'timeILTO.txt')
    outvarid(188) = outvaridtype('olto','OLD comp olake mean temp'        ,i_mean , 0,0,'deg'  ,'degree Celsius','mapOLTO.txt' ,'timeOLTO.txt')
    outvarid(o_snowcover) = outvaridtype('cfsc','comp frac snowcover area'       ,i_mean , 0,0,'-'  ,'fraction','mapCFSC.txt' ,'timeCFSC.txt')
    outvarid(190) = outvaridtype('rfsc','rec frac snowcover area'        ,i_mean , 0,0,'-'  ,'fraction','mapRFSC.txt' ,'timeRFSC.txt')
    outvarid(o_snowmax) = outvaridtype('smax','comp snowmax in winter'         ,i_mean , 0,0,'mm' ,'mm','mapSMAX.txt' ,'timeSMAX.txt')
    outvarid(192) = outvaridtype('rfse','rec frac snowcover area error'  ,i_mean , 0,0,'-'  ,'fraction','mapRFSE.txt' ,'timeRFSE.txt')
    outvarid(193) = outvaridtype('rfsm','rec frac snowcover multi'       ,i_mean , 0,0,'-'  ,'fraction','mapRFSM.txt' ,'timeRFSM.txt')
    outvarid(194) = outvaridtype('rfme','rec frac snowcover multi error' ,i_mean , 0,0,'-'  ,'fraction','mapRFME.txt' ,'timeRFME.txt')
    outvarid(o_landevap) = outvaridtype('levp','land evaporation'               ,i_sum  , 0,2,'mm' ,'mm'      ,'mapLEVP.txt' ,'timeLEVP.txt')
    outvarid(196) = outvaridtype('som2','soil water upper 2 l',i_mean , 0,0,'mm' ,'mm'            ,'mapSOM2.txt' ,'timeSOM2.txt')
    outvarid(o_sml1) = outvaridtype('sml1','soil moisture lay 1' ,i_mean , 0,0,'mm' ,'mm'            ,'mapSML1.txt' ,'timeSML1.txt')
    outvarid(o_sml2) = outvaridtype('sml2','soil moisture lay 2' ,i_mean , 0,0,'mm' ,'mm'            ,'mapSML2.txt' ,'timeSML2.txt')
    outvarid(o_sml3) = outvaridtype('sml3','soil moisture lay 3' ,i_mean , 0,0,'mm' ,'mm'            ,'mapSML3.txt' ,'timeSML3.txt')
    outvarid(200) = outvaridtype('stsw','standing soil water' ,i_mean , 0,0,'mm' ,'mm'            ,'mapSTSW.txt' ,'timeSTSW.txt')
    outvarid(201) = outvaridtype('smrz','soil moisture root z',i_mean , 0,0,'mm' ,'mm'            ,'mapSMRZ.txt' ,'timeSMRZ.txt')
    outvarid(202) = outvaridtype('sm13','soil moisture sl 1-3',i_mean , 0,0,'mm' ,'mm'            ,'mapSM13.txt' ,'timeSM13.txt')
    outvarid(203) = outvaridtype('clCO','comp local conc OC',i_wmean, o_cloc,0,'mg/L' ,'mg org-C/L'  ,'mapCLCO.txt','timeCLCO.txt')    !clOC occupied (cloc)
    outvarid(204) = outvaridtype('lrdp','local river depth'    ,i_mean , 0,0,'m'  ,'m','mapLRDP.txt' ,'timeLRDP.txt')
    outvarid(205) = outvaridtype('mrdp','main river depth'     ,i_mean , 0,0,'m'  ,'m','mapMRDP.txt' ,'timeMRDP.txt')
    outvarid(o_icloss) = outvaridtype('icpe','interception loss'    ,i_sum , 0,0,'mm'  ,'mm','mapICPE.txt' ,'timeICPE.txt')
    outvarid(207) = outvaridtype('rlIN','load rec IN com flow',i_sum, 0,0,'kg' ,'kg Inorg-N'  ,'mapRLIN.txt','timeRLIN.txt')
    outvarid(208) = outvaridtype('rlON','load rec ON com flow',i_sum, 0,0,'kg' ,'kg Org-N'    ,'mapRLON.txt','timeRLON.txt')
    outvarid(209) = outvaridtype('rlSP','load rec SP com flow',i_sum, 0,0,'kg' ,'kg SRP-P'    ,'mapRLSP.txt','timeRLSP.txt')
    outvarid(210) = outvaridtype('rlPP','load rec PP com flow',i_sum, 0,0,'kg' ,'kg PartP-P'  ,'mapRLPP.txt','timeRLPP.txt')
    outvarid(211) = outvaridtype('rlTN','load rec TN com flow',i_sum, 0,0,'kg' ,'kg Tot-N'    ,'mapRLTN.txt','timeRLTN.txt')
    outvarid(212) = outvaridtype('rlTP','load rec TP com flow',i_sum, 0,0,'kg' ,'kg Tot-P'    ,'mapRLTP.txt','timeRLTP.txt')
    outvarid(213) = outvaridtype('rlOC','load rec OC com flow',i_sum, 0,0,'kg' ,'kg org-C'    ,'mapRLOC.txt','timeRLOC.txt')
    outvarid(o_smffc)   = outvaridtype('srff','sm frac of field cap'     ,i_mean , 0,0,'-'  ,'-','mapSRFF.txt' ,'timeSRFF.txt')
    outvarid(o_smfdep)  = outvaridtype('smfd','sm frac of depth'         ,i_mean , 0,0,'-'  ,'-','mapSMFD.txt' ,'timeSMFD.txt')
    outvarid(o_smrzfdep) = outvaridtype('srfd','sm root frac of depth'   ,i_mean , 0,0,'-'  ,'-','mapSRFD.txt' ,'timeSRFD.txt')
    outvarid(o_smfpw)   = outvaridtype('smfp','sm frac of pore v'        ,i_mean , 0,2,'-'  ,'-','mapSMFP.txt' ,'timeSMFP.txt')
    outvarid(o_smrzfpw) = outvaridtype('srfp','sm root frac of pore v'   ,i_mean , 0,0,'-'  ,'-','mapSRFP.txt' ,'timeSRFP.txt')
    outvarid(o_upfp)    = outvaridtype('upfp','aver upstream smfp',i_mean, 0,0,'mm','mm','mapUPFP.txt','timeUPFP.txt')
    outvarid(o_xobsm)   = outvaridtype('xom0','recorded mean var 0',i_mean , 0,0,'?','?' ,'mapXOM0.txt','timeXOM0.txt')
    outvarid(o_xobsm+1) = outvaridtype('xom1','recorded mean var 1',i_mean , 0,0,'?','?' ,'mapXOM1.txt','timeXOM1.txt')
    outvarid(o_xobsm+2) = outvaridtype('xom2','recorded mean var 2',i_mean , 0,0,'?','?' ,'mapXOM2.txt','timeXOM2.txt')
    outvarid(o_xobsm+3) = outvaridtype('xom3','recorded mean var 3',i_mean , 0,0,'?','?' ,'mapXOM3.txt','timeXOM3.txt')
    outvarid(o_xobsm+4) = outvaridtype('xom4','recorded mean var 4',i_mean , 0,0,'?','?' ,'mapXOM4.txt','timeXOM4.txt')
    outvarid(o_xobsm+5) = outvaridtype('xom5','recorded mean var 5',i_mean , 0,0,'?','?' ,'mapXOM5.txt','timeXOM5.txt')
    outvarid(o_xobsm+6) = outvaridtype('xom6','recorded mean var 6',i_mean , 0,0,'?','?' ,'mapXOM6.txt','timeXOM6.txt')
    outvarid(o_xobsm+7) = outvaridtype('xom7','recorded mean var 7',i_mean , 0,0,'?','?' ,'mapXOM7.txt','timeXOM7.txt')
    outvarid(o_xobsm+8) = outvaridtype('xom8','recorded mean var 8',i_mean , 0,0,'?','?' ,'mapXOM8.txt','timeXOM8.txt')
    outvarid(o_xobsm+9) = outvaridtype('xom9','recorded mean var 9',i_mean , 0,0,'?','?' ,'mapXOM9.txt','timeXOM9.txt')
    outvarid(o_xobss)   = outvaridtype('xos0','recorded sum var 0',i_sum  , 0,0,'?','?' ,'mapXOS0.txt','timeXOS0.txt')
    outvarid(o_xobss+1) = outvaridtype('xos1','recorded sum var 1',i_sum  , 0,0,'?','?' ,'mapXOS1.txt','timeXOS1.txt')
    outvarid(o_xobss+2) = outvaridtype('xos2','recorded sum var 2',i_sum  , 0,0,'?','?' ,'mapXOS2.txt','timeXOS2.txt')
    outvarid(o_xobss+3) = outvaridtype('xos3','recorded sum var 3',i_sum  , 0,0,'?','?' ,'mapXOS3.txt','timeXOS3.txt')
    outvarid(o_xobss+4) = outvaridtype('xos4','recorded sum var 4',i_sum  , 0,0,'?','?' ,'mapXOS4.txt','timeXOS4.txt')
    outvarid(o_xobss+5) = outvaridtype('xos5','recorded sum var 5',i_sum  , 0,0,'?','?' ,'mapXOS5.txt','timeXOS5.txt')
    outvarid(o_xobss+6) = outvaridtype('xos6','recorded sum var 6',i_sum  , 0,0,'?','?' ,'mapXOS6.txt','timeXOS6.txt')
    outvarid(o_xobss+7) = outvaridtype('xos7','recorded sum var 7',i_sum  , 0,0,'?','?' ,'mapXOS7.txt','timeXOS7.txt')
    outvarid(o_xobss+8) = outvaridtype('xos8','recorded sum var 8',i_sum  , 0,0,'?','?' ,'mapXOS8.txt','timeXOS8.txt')
    outvarid(o_xobss+9) = outvaridtype('xos9','recorded sum var 9',i_sum  , 0,0,'?','?' ,'mapXOS9.txt','timeXOS9.txt')
    outvarid(240) = outvaridtype('aqin','aquifer recharge'   ,i_sum , 0,0,'m3'  ,'m3','mapAQIN.txt' ,'timeAQIN.txt')
    outvarid(241) = outvaridtype('aqut','aquifer return flow',i_sum , 0,0,'m3'  ,'m3','mapAQUT.txt' ,'timeAQUT.txt')
    outvarid(242) = outvaridtype('aqwl','aquifer water depth',i_mean, 0,0,'m'   ,'m','mapAQWL.txt' ,'timeAQWL.txt')
    outvarid(o_uppr) = outvaridtype('uppr','aver upstream prec',i_sum, 0,0,'mm','mm','mapUPPR.txt','timeUPPR.txt')
    outvarid(o_upev) = outvaridtype('upev','aver upstream evap',i_sum, 0,0,'mm','mm','mapUPEV.txt','timeUPEV.txt')
    outvarid(o_upsn) = outvaridtype('upsn','aver upstream snoww',i_mean, 0,0,'mm','mm','mapUPSN.txt','timeUPSN.txt')
    outvarid(o_upso) = outvaridtype('upso','aver upstream soilw',i_mean, 0,0,'mm','mm','mapUPSO.txt','timeUPSO.txt')
    outvarid(o_specificq) = outvaridtype('speq','specific discharge',i_sum, 0,0,'mm','mm','mapSPEQ.txt','timeSPEQ.txt')    !Changed name of variable upro
    outvarid(248) = outvaridtype('coic','comp olake ice cover',i_mean, 0,0,'-','-','mapCOIC.txt','timeCOIC.txt')
    outvarid(249) = outvaridtype('ciic','comp ilake ice cover',i_mean, 0,0,'-','-','mapCIIC.txt','timeCIIC.txt')
    outvarid(250) = outvaridtype('cmic','comp m river ice cov',i_mean, 0,0,'-','-','mapCMIC.txt','timeCMIC.txt')
    outvarid(251) = outvaridtype('clic','comp l river ice cov',i_mean, 0,0,'-','-','mapCLIC.txt','timeCLIC.txt')

    !Glacier outputs and corresponding WGMS observations, using the following 4-letter code structure:
    !  xGnn, where x  = r or c (recorded or computed) and nn = 2-letter code identifying the glacier variable
    outvarid(252) = outvaridtype('cgmb','comp. glacier mass balance',i_mean, 0,0,'mm','mm','mapCGMB.txt','timeCGMB.txt')
    outvarid(253) = outvaridtype('rgmb','rec. glacier mass balance',i_mean, 0,0,'mm','mm','mapRGMB.txt','timeRGMB.txt')
    outvarid(254) = outvaridtype('cgma','area used in com. mass balance',i_mean, 0,0,'km2','km2','mapCGMA.txt','timeCGMA.txt')
    outvarid(255) = outvaridtype('rgma','area used in rec. mass balance',i_mean, 0,0,'km2','km2','mapRGMA.txt','timeRGMA.txt')
    outvarid(256) = outvaridtype('rgmp','rec. mass balance period',i_mean, 0,0,'days','days','mapRGMP.txt','timeRGMP.txt')
    
    !Glacier Xobs variables not yet included in the code:
    ! ---------------------------------------------------
    !xGEL Equilibrium Line Altitude  ELA                	(m)
    !xGAA Accumulation Area Ratio of total area         	(%)
    !xGGA Area Survey year (GTC, geodetic thicness changes)	(km2)
    !xGAC Area Change (GTC)                                 (1000 m2)
    !xGTC Thickness Change (GTC)                            (mm)
    !xGVC Volume Change (GTC)                               (1000 m3)
    !xGSP GTC Survey period, days since Reference survey    (days)

    !Sovjet Snow Coarse observations in forest and open areas, and corresponding model variables:
    outvarid(257) = outvaridtype('S105','fsusnow fsc surr. open',i_mean, 0,0,'0-10','0-10','mapS105.txt','timeS105.txt')
    outvarid(258) = outvaridtype('S106','fsusnow fsc course open',i_mean, 0,0,'0-10','0-10','mapS106.txt','timeS106.txt')
    outvarid(259) = outvaridtype('S108','fsusnow mean depth open',i_mean, 0,0,'cm','cm','mapS108.txt','timeS108.txt')
    outvarid(260) = outvaridtype('S111','fsusnow mean density open',i_mean, 0,0,'g/cm3','g/cm3','mapS111.txt','timeS111.txt')
    outvarid(261) = outvaridtype('S114','fsusnow snow water eq. open',i_mean, 0,0,'mm','mm','mapS114.txt','timeS114.txt')
    outvarid(262) = outvaridtype('S205','fsusnow fsc surr. forest',i_mean, 0,0,'0-10','0-10','mapS205.txt','timeS205.txt')
    outvarid(263) = outvaridtype('S206','fsusnow fsc course forest',i_mean, 0,0,'0-10','0-10','mapS206.txt','timeS206.txt')
    outvarid(264) = outvaridtype('S208','fsusnow mean depth forest',i_mean, 0,0,'cm','cm','mapS208.txt','timeS208.txt')
    outvarid(265) = outvaridtype('S211','fsusnow mean density forest',i_mean, 0,0,'g/cm3','g/cm3','mapS211.txt','timeS211.txt')
    outvarid(266) = outvaridtype('S214','fsusnow snow water eq. forest ',i_mean, 0,0,'mm','mm','mapS214.txt','timeS214.txt')
    outvarid(267) = outvaridtype('C106','comp. fsc course open',i_mean, 0,0,'0-10','0-10','mapC106.txt','timeC106.txt')
    outvarid(268) = outvaridtype('C108','comp. mean depth open',i_mean, 0,0,'cm','cm','mapC108.txt','timeC108.txt')
    outvarid(269) = outvaridtype('C111','comp. mean density open',i_mean, 0,0,'g/cm3','g/cm3','mapC111.txt','timeC111.txt')
    outvarid(270) = outvaridtype('C114','comp. snow water eq. open',i_mean, 0,0,'mm','mm','mapC114.txt','timeC114.txt')
    outvarid(271) = outvaridtype('C206','comp. fsc course forest',i_mean, 0,0,'0-10','0-10','mapC206.txt','timeC206.txt')
    outvarid(272) = outvaridtype('C208','comp. mean depth forest',i_mean, 0,0,'cm','cm','mapC208.txt','timeC208.txt')
    outvarid(273) = outvaridtype('C211','comp. mean density forest',i_mean, 0,0,'g/cm3','g/cm3','mapC211.txt','timeC211.txt')
    outvarid(274) = outvaridtype('C214','comp. snow water eq. forest ',i_mean, 0,0,'mm','mm','mapC214.txt','timeC214.txt')
!    outvarid(o_snowfall) = outvaridtype('snwf','snowfall',i_sum, 0,1,'mm','mm','mapSNSF.txt','timeSNSF.txt')
!    outvarid(o_rainfall) = outvaridtype('raif','rainfall',i_sum, 0,1,'mm','mm','mapRAIF.txt','timeRAIF.txt')
    outvarid(o_uppe) = outvaridtype('uppe','aver upstream pot evap',i_sum, 0,0,'mm','mm','mapUPPE.txt','timeUPPE.txt')
    outvarid(o_evapsnow) = outvaridtype('evsn','snow evaporation',i_sum , 0,0,'mm' ,'mm','mapEVSN.txt' ,'timeEVSN.txt')
    outvarid(279) = outvaridtype('evpt','total evaporation',i_sum , 0,0,'mm' ,'mm','mapEVPT.txt' ,'timeEVPT.txt')
    outvarid(280) = outvaridtype('clwc','cleaned wcom',i_mean , 0,0,'m','meter' ,'mapCLWC.txt','timeCLWC.txt')
    outvarid(281) = outvaridtype('clws','cleaned wstr',i_mean , 0,0,'m','meter' ,'mapCLWS.txt','timeCLWS.txt')
    outvarid(282) = outvaridtype('wcav','ave. comp waterstage', i_mean , 0,0,'m','meter' ,'mapWCAV.txt','timeWCAV.txt')
    outvarid(283) = outvaridtype('wtm0','comp water temp >0' , i_mean , 0,0,'deg'  ,'degree Celsius','mapWTM0.txt' ,'timeWTM0.txt')
    outvarid(o_psim) = outvaridtype('psim','sim. precipitation', i_sum , 0,0,'mm','mm','mapPSIM.txt' ,'timePSIM.txt')

    !Soil load output variables; 1-12 soillayer 1-2, 13-24 soillayer 3, odd=gross load, even=net load, 25-36 soillayer3+tiledrain
    outvarid(285) = outvaridtype('sl07','soill 1-2 grs ld SP',i_sum , 0,0,'kg' ,'kg','mapSL07.txt' ,'timeSL07.txt')
    outvarid(286) = outvaridtype('sl08','soill 1-2 net ld SP',i_sum , 0,0,'kg' ,'kg','mapSL08.txt' ,'timeSL08.txt')
    outvarid(287) = outvaridtype('sl09','soill 1-2 grs ld PP',i_sum , 0,0,'kg' ,'kg','mapSL09.txt' ,'timeSL09.txt')
    outvarid(288) = outvaridtype('sl10','soill 1-2 net ld PP',i_sum , 0,0,'kg' ,'kg','mapSL10.txt' ,'timeSL10.txt')
    outvarid(289) = outvaridtype('sl11','soill 1-2 grs ld TP',i_sum , 0,0,'kg' ,'kg','mapSL11.txt' ,'timeSL11.txt')
    outvarid(290) = outvaridtype('sl12','soill 1-2 net ld TP',i_sum , 0,0,'kg' ,'kg','mapSL12.txt' ,'timeSL12.txt')
    outvarid(291) = outvaridtype('sl13','soill 3 grs ld IN',i_sum , 0,0,'kg' ,'kg','mapSL13.txt' ,'timeSL13.txt')
    outvarid(292) = outvaridtype('sl14','soill 3 net ld IN',i_sum , 0,0,'kg' ,'kg','mapSL14.txt' ,'timeSL14.txt')
    outvarid(293) = outvaridtype('sl15','soill 3 grs ld ON',i_sum , 0,0,'kg' ,'kg','mapSL15.txt' ,'timeSL15.txt')
    outvarid(294) = outvaridtype('sl16','soill 3 net ld ON',i_sum , 0,0,'kg' ,'kg','mapSL16.txt' ,'timeSL16.txt')
    outvarid(295) = outvaridtype('sl17','soill 3 grs ld TN',i_sum , 0,0,'kg' ,'kg','mapSL17.txt' ,'timeSL17.txt')
    outvarid(296) = outvaridtype('sl18','soill 3 net ld TN',i_sum , 0,0,'kg' ,'kg','mapSL18.txt' ,'timeSL18.txt')
    outvarid(297) = outvaridtype('sl19','soill 3 grs ld SP',i_sum , 0,0,'kg' ,'kg','mapSL19.txt' ,'timeSL19.txt')
    outvarid(298) = outvaridtype('sl20','soill 3 net ld SP',i_sum , 0,0,'kg' ,'kg','mapSL20.txt' ,'timeSL20.txt')
    outvarid(299) = outvaridtype('sl21','soill 3 grs ld PP',i_sum , 0,0,'kg' ,'kg','mapSL21.txt' ,'timeSL21.txt')
    outvarid(300) = outvaridtype('sl22','soill 3 net ld PP',i_sum , 0,0,'kg' ,'kg','mapSL22.txt' ,'timeSL22.txt')
    outvarid(301) = outvaridtype('sl23','soill 3 grs ld TP',i_sum , 0,0,'kg' ,'kg','mapSL23.txt' ,'timeSL23.txt')
    outvarid(302) = outvaridtype('sl24','soill 3 net ld TP',i_sum , 0,0,'kg' ,'kg','mapSL24.txt' ,'timeSL24.txt')
    outvarid(o_soilden3) = outvaridtype('den3','denitrification sl3'  ,i_sum , 0,0,'kg/km2' ,'kg N/km2'  ,'mapDEN3.txt' ,'timeDEN3.txt')
    outvarid(o_soildenrz) = outvaridtype('denz','denitrification sl12'  ,i_sum , 0,0,'kg/km2' ,'kg N/km2'  ,'mapDENZ.txt' ,'timeDENZ.txt')
    outvarid(305) = outvaridtype('sl01','soill 1-2 grs ld IN',i_sum , 0,0,'kg' ,'kg','mapSL01.txt' ,'timeSL01.txt')
    outvarid(306) = outvaridtype('sl02','soill 1-2 net ld IN',i_sum , 0,0,'kg' ,'kg','mapSL02.txt' ,'timeSL02.txt')
    outvarid(307) = outvaridtype('sl03','soill 1-2 grs ld ON',i_sum , 0,0,'kg' ,'kg','mapSL03.txt' ,'timeSL03.txt')
    outvarid(308) = outvaridtype('sl04','soill 1-2 net ld ON',i_sum , 0,0,'kg' ,'kg','mapSL04.txt' ,'timeSL04.txt')
    outvarid(309) = outvaridtype('sl05','soill 1-2 grs ld TN',i_sum , 0,0,'kg' ,'kg','mapSL05.txt' ,'timeSL05.txt')
    outvarid(310) = outvaridtype('sl06','soill 1-2 net ld TN',i_sum , 0,0,'kg' ,'kg','mapSL06.txt' ,'timeSL06.txt')
    outvarid(311) = outvaridtype('sl25','soill 3+t grs ld IN',i_sum , 0,0,'kg' ,'kg','mapSL25.txt' ,'timeSL25.txt')
    outvarid(312) = outvaridtype('sl26','soill 3+t net ld IN',i_sum , 0,0,'kg' ,'kg','mapSL26.txt' ,'timeSL26.txt')
    outvarid(313) = outvaridtype('sl27','soill 3+t grs ld ON',i_sum , 0,0,'kg' ,'kg','mapSL27.txt' ,'timeSL27.txt')
    outvarid(314) = outvaridtype('sl28','soill 3+t net ld ON',i_sum , 0,0,'kg' ,'kg','mapSL28.txt' ,'timeSL28.txt')
    outvarid(315) = outvaridtype('sl29','soill 3+t grs ld TN',i_sum , 0,0,'kg' ,'kg','mapSL29.txt' ,'timeSL29.txt')
    outvarid(316) = outvaridtype('sl30','soill 3+t net ld TN',i_sum , 0,0,'kg' ,'kg','mapSL30.txt' ,'timeSL30.txt')
    outvarid(317) = outvaridtype('sl31','soill 3+t grs ld SP',i_sum , 0,0,'kg' ,'kg','mapSL31.txt' ,'timeSL31.txt')
    outvarid(318) = outvaridtype('sl32','soill 3+t net ld SP',i_sum , 0,0,'kg' ,'kg','mapSL32.txt' ,'timeSL32.txt')
    outvarid(319) = outvaridtype('sl33','soill 3+t grs ld PP',i_sum , 0,0,'kg' ,'kg','mapSL33.txt' ,'timeSL33.txt')
    outvarid(320) = outvaridtype('sl34','soill 3+t net ld PP',i_sum , 0,0,'kg' ,'kg','mapSL34.txt' ,'timeSL34.txt')
    outvarid(321) = outvaridtype('sl35','soill 3+t grs ld TP',i_sum , 0,0,'kg' ,'kg','mapSL35.txt' ,'timeSL35.txt')
    outvarid(322) = outvaridtype('sl36','soill 3+t net ld TP',i_sum , 0,0,'kg' ,'kg','mapSL36.txt' ,'timeSL36.txt')
    outvarid(323) = outvaridtype('mrfp','m river flpl depth' ,i_mean , 0,0,'m'  ,'m','mapMRFP.txt' ,'timeMRFP.txt')
    outvarid(324) = outvaridtype('olfp','olake flpl depth'   ,i_mean , 0,0,'m'  ,'m','mapOLFP.txt' ,'timeOLFP.txt')
    outvarid(325) = outvaridtype('mrfg','m river flpl degree',i_mean , 0,0,'%'  ,'% of flpl area','mapMRFG.txt' ,'timeMRFG.txt')
    outvarid(326) = outvaridtype('olfg','olake flpl degree'  ,i_mean , 0,0,'%'  ,'% of flpl area','mapOLFG.txt' ,'timeOLFG.txt')
    outvarid(o_csoillayerIN(1)) = outvaridtype('cIN1','conc IN soil layer 1' ,i_wmean, o_sml9,0,'ug/L','ug Inorg-N/L'  ,'mapCIN1.txt','timeCIN1.txt')
    outvarid(o_csoillayerIN(2)) = outvaridtype('cIN2','conc IN soil layer 2' ,i_wmean, o_sml2,0,'ug/L','ug Inorg-N/L'  ,'mapCIN2.txt','timeCIN2.txt')
    outvarid(o_csoillayerIN(3)) = outvaridtype('cIN3','conc IN soil layer 3' ,i_wmean, o_sml3,0,'ug/L','ug Inorg-N/L'  ,'mapCIN3.txt','timeCIN3.txt')
    outvarid(o_sml9) = outvaridtype('sml9','soil water layer 1' ,i_mean, 0,0,'mm','mm'  ,'mapSML9.txt','timeSML9.txt')
    outvarid(o_snowmelt) = outvaridtype('melt','snow melt' ,i_mean, 0,2,'mm','mm'  ,'mapMELT.txt','timeMELT.txt')
    outvarid(o_upme) = outvaridtype('upme','aver upstr snow melt' ,i_sum, 0,0,'mm','mm'  ,'mapUPME.txt','timeUPME.txt')
    !T1 tracer output
    outvarid(o_ppT1(1)) = outvaridtype('aT11','pool adsorbed T1 soil1',i_mean , 0,0,'?' ,'?'    ,'mapAT11.txt','timeAT11.txt')   
    outvarid(o_ppT1(2)) = outvaridtype('aT12','pool adsorbed T1 soil2',i_mean , 0,0,'?' ,'?'    ,'mapAT12.txt','timeAT12.txt')
    outvarid(o_ppT1(3)) = outvaridtype('aT13','pool adsorbed T1 soil3',i_mean , 0,0,'?' ,'?'    ,'mapAT13.txt','timeAT13.txt')
    outvarid(336) = outvaridtype('sT11','pool T1 soilwater soil1',i_mean , 0,0,'?' ,'?'    ,'mapST11.txt','timeST11.txt')
    outvarid(337) = outvaridtype('sT12','pool T1 soilwater soil2',i_mean , 0,0,'?' ,'?'    ,'mapST12.txt','timeST12.txt')
    outvarid(338) = outvaridtype('sT13','pool T1 soilwater soil3',i_mean , 0,0,'?' ,'?'    ,'mapST13.txt','timeST13.txt')
    outvarid(339) = outvaridtype('Tsmr','T1 pool main river sediments',i_mean , 0,0,'?' ,'?','mapTSMR.txt','timeTSMR.txt')
    outvarid(340) = outvaridtype('Tslr','T1 pool local river sediments',i_mean , 0,0,'?','?','mapTSLR.txt','timeTSLR.txt')
    outvarid(o_T1sf) = outvaridtype('T1sf','pool T1 soil surface',i_mean , 0,0,'?' ,'?'    ,'mapT1SF.txt','timeT1SF.txt')
    outvarid(342) = outvaridtype('clT1','comp local conc T1',i_wmean, o_cloc,0,'?' , '?', 'mapCLT1.txt', 'timeCLT1.txt')
    outvarid(o_cro1T1) = outvaridtype('Tcr1','T1 in comp runoff 1',i_wmean, o_ro1,0,'?' , '?', 'mapTCR1.txt', 'timeTCR1.txt')
    outvarid(o_cro2T1) = outvaridtype('Tcr2','T1 in comp runoff 2',i_wmean, o_ro2,0,'?' , '?', 'mapTCR2.txt', 'timeTCR2.txt')
    outvarid(o_cro3T1) = outvaridtype('Tcr3','T1 in comp runoff 3',i_wmean, o_ro3,0,'?' , '?', 'mapTCR3.txt', 'timeTCR3.txt')
    outvarid(o_crodT1) = outvaridtype('Tcrd','T1 in comp rf drain',i_wmean, o_rod,0,'?' , '?', 'mapTCRD.txt', 'timeTCRD.txt')
    outvarid(o_crosT1) = outvaridtype('Tcrs','T1 in comp surf rf',i_wmean, o_ros,0,'?' , '?', 'mapTCRS.txt', 'timeTCRS.txt')
    
    !Additional outvars for data assimilation 'observation operators' to derive the 'model equivalents'
    !DavidG, RafaelP, LottaP 2017-02-23
    outvarid(o_clrf) = outvaridtype('clrf','calculated local runoff', i_sum , 0,0,'mm' ,'mm' ,'mapCLRF.txt' ,'timeCLRF.txt') !normalized from cloc (m3/s) by local basin area (excluding main river and olake)
    outvarid(o_clrp) = outvaridtype('clrp','calculated prec on river and lake', i_sum , 0,0,'mm' ,'mm' ,'mapCLRP.txt' ,'timeCLRP.txt')
    outvarid(o_clre) = outvaridtype('clre','calculated evap on river and lake', i_sum , 0,0,'mm' ,'mm' ,'mapCLRE.txt' ,'timeCLRE.txt')
    outvarid(o_qcin) = outvaridtype('qcin','calculated (regional) lake inflow',i_mean , 0,0,'m3/s','m3/s' ,'mapQCIN.txt','timeQCIN.txt')
    outvarid(o_qrin) = outvaridtype('qrin','recorded (regional) lake inflow',i_mean , 0,0,'m3/s','m3/s' ,'mapQRIN.txt','timeQRIN.txt')
    
    ! for shape project
    outvarid(o_rmwl) = outvaridtype('rmwl','recorded main river water level',i_mean , 0,0,'m','m' ,'mapRMWL.txt','timeRMWL.txt')
    outvarid(o_rmwd) = outvaridtype('rmwd','recorded main river water flow',i_mean , 0,0,'m3/s','m3/s' ,'mapRMWD.txt','timeRMWD.txt')
    outvarid(o_cmwl) = outvaridtype('cmwl','calculated main river water level',i_mean , 0,0,'m','m' ,'mapCMWL.txt','timeCMWL.txt')
    outvarid(o_amwl) = outvaridtype('amwl','altimetry main river water level',i_mean , 0,0,'m','m' ,'mapAMWL.txt','timeAMWL.txt')

    
    ! to be made:
    !  - AMSR2 radiances in different frequencies (CMEM, FASTEM, etc)
    !  - In-situ snow observations (or by additional "station-subbasins" representing the local conditions of the snow observation as good as possible with the model)
    !  - Snow survey observations (mean, variability, etc)
    !  - whatever exciting observation we may come up with...
    

    !Max length: 1st==4, 2nd=20, 3rd=Int, 4th=Int, 5th=5, 6th=20, 7th=20, 8th=20
    !i_mean=medel per dag, i_sum=summa över vald period (eller år för hela sim perioden)
    !i_wmean=volymviktat med vattnet i kolumnen efteråt
    
    !Set variable exchange for criterion calculation (outvarid-index)
    ALLOCATE(changecritvar(2,2))
    changecritvar = RESHAPE(SOURCE=(/51,o_rewstr,280,281/),SHAPE=(/2,2/))
    
    !Preparation of heading for load files
    loadheadings=['subid ','WetAtm','DryAtm','Fertil','PDecay','RuralA','GrwSIn','IrrSrc','Runoff','RuralB','Urban1','Urban2','Urban3','Rgrwmr','Rgrwol','A     ','B     ','C     ','D     ','E     ','F     ','G     ','H     ','I     ','J     ','K     ','L     ','M     ','N     ']
    
  END SUBROUTINE initiate_output_variables
  
  !>Set variable holding model parameter information of HYPE model
  !------------------------------------------------------------------------------------------
  SUBROUTINE initiate_model_parameters()
  !>Note that only small letters allowed in parameter names
  !>modparid form: name, dependence type, index (defined in hypevar)
  
    max_par = 294          !maximum number of parameters
    ALLOCATE(modparid(max_par))
    modparid = modparidtype('',0,0) !default initialisation
    
    modparid(n_lp)= modparidtype('lp        ', m_gpar , m_lp)
    modparid(n_cevpa)   = modparidtype('cevpam    ', m_gpar , m_cevpam)
    modparid(n_cevpp)   = modparidtype('cevpph    ', m_gpar , m_cevpph)
    modparid(4)   = modparidtype('deadl     ', m_gpar , m_deadl)
    modparid(5)   = modparidtype('fastn0    ', m_gpar , m_fastN0)
    modparid(6)   = modparidtype('fastp0    ', m_gpar , m_fastP0)
    modparid(7)   = modparidtype('init1     ', m_gpar , m_iniT1)
    modparid(8)   = modparidtype('init2     ', m_gpar , m_iniT2)
    modparid(9)   = modparidtype('t1evap    ', m_gpar , m_T1evap)
    modparid(n_rivv)  = modparidtype('rivvel    ', m_gpar , m_rivvel)
    modparid(n_damp)  = modparidtype('damp      ', m_gpar , m_damp)
    modparid(n_tcalt)  = modparidtype('tcalt     ', m_gpar , m_tcalt)
    modparid(13)  = modparidtype('gratk     ', m_gpar , m_grat1)
    modparid(14)  = modparidtype('gratp     ', m_gpar , m_grat2)
    modparid(15)  = modparidtype('denitrlu  ', m_lpar , m_denitrlu)
    modparid(16)  = modparidtype('denitwrm  ', m_gpar , m_denitwr)
    modparid(17)  = modparidtype('maxwidth  ', m_gpar , m_maxwidth)
    modparid(18)  = modparidtype('dissolfn  ', m_lpar , m_dissolfN)
    modparid(19)  = modparidtype('dissolfp  ', m_lpar , m_dissolfP)
    modparid(20)  = modparidtype('litterdays', m_gpar , m_littdays)
    modparid(n_tcea)  = modparidtype('tcelevadd ', m_gpar , m_tcelevadd)
    modparid(n_pcem)  = modparidtype('pcelevmax ', m_gpar , m_pcelevmax)
    modparid(n_tcorr)  = modparidtype('tempcorr  ', m_rpar , m_tcadd)   
    modparid(n_pcea)  = modparidtype('pcelevadd ', m_gpar , m_pcelevadd)
    modparid(25)  = modparidtype('fastlake  ', m_gpar , m_fastlake) 
    modparid(26)  = modparidtype('epotdist  ', m_gpar , m_epotdist)
    modparid(27)  = modparidtype('qmean     ', m_gpar , m_qmean)
    modparid(28)  = modparidtype('pcaddg    ', m_gpar , m_pcaddg)
    modparid(n_pcet)  = modparidtype('pcelevth  ', m_gpar , m_pcelevth)
    modparid(30)  = modparidtype('tpmean    ', m_lrpar, m_tpmean)
    modparid(31)  = modparidtype('ttmp      ', m_lpar , m_ttmp)
    modparid(32)  = modparidtype('cmlt      ', m_lpar , m_cmlt)
    modparid(33)  = modparidtype('cevp      ', m_lpar , m_cevp)
    modparid(34)  = modparidtype('frost     ', m_lpar , m_cfrost)
    modparid(35)  = modparidtype('srrcs     ', m_lpar , m_srrcs)
    modparid(36)  = modparidtype('ttpd      ', m_gpar , m_ttpd)    
    modparid(37)  = modparidtype('ttpi      ', m_gpar , m_ttpi)    
    modparid(38)  = modparidtype('dissolhn  ', m_lpar , m_dissolhN)
    modparid(39)  = modparidtype('dissolhp  ', m_lpar , m_dissolhP)
    modparid(45)  = modparidtype('humusn0   ', m_lpar , m_humusN0) 
    modparid(46)  = modparidtype('partp0    ', m_lpar , m_partP0)  
    modparid(47)  = modparidtype('wcfc      ', m_spar , m_wcfc)
    modparid(48)  = modparidtype('wcwp      ', m_spar , m_wcwp)
    modparid(49)  = modparidtype('wcep      ', m_spar , m_wcep)
    modparid(50)  = modparidtype('wcfc1     ', m_spar , m_wcfc1)
    modparid(51)  = modparidtype('wcwp1     ', m_spar , m_wcwp1)
    modparid(52)  = modparidtype('wcep1     ', m_spar , m_wcep1)
    modparid(53)  = modparidtype('wcfc2     ', m_spar , m_wcfc2)
    modparid(54)  = modparidtype('wcwp2     ', m_spar , m_wcwp2)
    modparid(55)  = modparidtype('wcep2     ', m_spar , m_wcep2)
    modparid(56)  = modparidtype('wcfc3     ', m_spar , m_wcfc3)
    modparid(57)  = modparidtype('wcwp3     ', m_spar , m_wcwp3)
    modparid(58)  = modparidtype('wcep3     ', m_spar , m_wcep3)
    modparid(59)  = modparidtype('gldepi    ', m_gpar , m_gldepi)
    modparid(60)  = modparidtype('trrcs     ', m_spar , m_trrcs)
    modparid(61)  = modparidtype('mperc1    ', m_spar , m_perc1)
    modparid(62)  = modparidtype('mperc2    ', m_spar , m_perc2)
    modparid(65)  = modparidtype('sswcorr   ', m_gpar , m_sswcorr)   
    modparid(66)  = modparidtype('depthrel  ', m_lpar , m_depthrel)
    modparid(67)  = modparidtype('regirr    ', m_gpar , m_regirr)
    modparid(68)  = modparidtype('immdepth  ', m_gpar , m_immdep)
    modparid(69)  = modparidtype('iwdfrac   ', m_gpar , m_iwdfrac)
    modparid(70)  = modparidtype('irrdemand ', m_gpar , m_wdpar)  
    modparid(71)  = modparidtype('minerfn   ', m_lpar , m_minerfN)
    modparid(72)  = modparidtype('minerfp   ', m_lpar , m_minerfP)
    modparid(73)  = modparidtype('degradhn  ', m_lpar , m_degradhN)
    modparid(74)  = modparidtype('deepmem   ', m_gpar , m_deepmem)
    modparid(75)  = modparidtype('surfmem   ', m_lpar , m_surfmem)
    modparid(76)  = modparidtype('freuc     ', m_spar , m_freuc)  
    modparid(77)  = modparidtype('freuexp   ', m_spar , m_freuexp)
    modparid(78)  = modparidtype('freurate  ', m_spar , m_freurate)
    modparid(80)  = modparidtype('wprodn    ', m_gpar , m_wprodn) 
    modparid(81)  = modparidtype('sedon     ', m_gpar , m_sedon)
    modparid(82)  = modparidtype('sedpp     ', m_gpar , m_sedpp)
    modparid(83)  = modparidtype('sedexp    ', m_gpar , m_sedexp)
    modparid(84)  = modparidtype('rcgrw     ', m_gpar , m_rcgrw)
    modparid(85)  = modparidtype('rrcs1     ', m_spar , m_rrcs1)
    modparid(86)  = modparidtype('rrcs2     ', m_spar , m_rrcs2)
    modparid(n_rrcs3)  = modparidtype('rrcs3     ', m_gpar , m_rrcs3)
    modparid(88)  = modparidtype('hnhalf    ', m_lpar , m_hNhalf)
    modparid(89)  = modparidtype('pphalf    ', m_lpar , m_pPhalf)
    modparid(90)  = modparidtype('locsoil   ', m_gpar , m_locsoil)
    modparid(91)  = modparidtype('bufffilt  ', m_lpar , m_filtPbuf)
    modparid(92)  = modparidtype('innerfilt ', m_lpar , m_filtPinner)
    modparid(93)  = modparidtype('otherfilt ', m_lpar , m_filtPother)
    modparid(94)  = modparidtype('drydeppp  ', m_lpar , m_drypp)
    modparid(95)  = modparidtype('wetdepsp  ', m_gpar , m_wetsp)
    modparid(96)  = modparidtype('rivvel1   ', m_lrpar, m_velpar1)
    modparid(97)  = modparidtype('rivvel2   ', m_lrpar, m_velpar2)
    modparid(98)  = modparidtype('rivvel3   ', m_lrpar, m_velpar3)
    modparid(99)  = modparidtype('rivwidth1 ', m_lrpar, m_widpar1)
    modparid(100) = modparidtype('rivwidth2 ', m_lrpar, m_widpar2)
    modparid(101) = modparidtype('rivwidth3 ', m_lrpar, m_widpar3)
    modparid(102) = modparidtype('srrate    ', m_spar , m_srrate)
    modparid(103) = modparidtype('macrate   ', m_spar , m_macrate) 
    modparid(104) = modparidtype('mactrinf  ', m_spar , m_mactrinf)
    modparid(105) = modparidtype('mactrsm   ', m_spar , m_mactrsm) 
    modparid(106) = modparidtype('soilcoh   ', m_spar , m_soilcoh) 
    modparid(107) = modparidtype('soilerod  ', m_spar , m_soilerod)
    modparid(108) = modparidtype('wprodp    ', m_gpar , m_wprodp)
    modparid(109) = modparidtype('sreroexp  ', m_gpar , m_sreroexp)
    modparid(110) = modparidtype('snowdensdt', m_gpar , m_dsndens)
    modparid(111) = modparidtype('sfrost    ', m_spar , m_sfrost)
    modparid(112) = modparidtype('macrofilt ', m_spar , m_macfilt)
    modparid(113) = modparidtype('fertdays  ', m_gpar , m_fertdays)
    modparid(114) = modparidtype('humusp0   ', m_lpar , m_humusP0)   
    modparid(115) = modparidtype('hphalf    ', m_lpar , m_hPhalf)   
    modparid(116) = modparidtype('degradhp  ', m_lpar , m_degradhP)
    modparid(n_cevpc) = modparidtype('cevpcorr  ', m_rpar , m_cevpcorr)
    modparid(118) = modparidtype('minc      ', m_gpar , m_minc)
    modparid(119) = modparidtype('humusc1   ', m_lpar , m_humusC1)
    modparid(120) = modparidtype('fastc1    ', m_lpar , m_fastC1)
    modparid(121) = modparidtype('klh       ', m_gpar , m_crate1)
    modparid(122) = modparidtype('klo       ', m_gpar , m_crate2)
    modparid(123) = modparidtype('kho       ', m_gpar , m_crate3)
    modparid(124) = modparidtype('tcobselev ', m_gpar , m_tcobselev)
    modparid(125) = modparidtype('ripz      ', m_lpar , m_ripz)
    modparid(126) = modparidtype('ripe      ', m_gpar , m_ripe)
    modparid(127) = modparidtype('sedoc     ', m_gpar , m_sedoc)
    modparid(128) = modparidtype('koc       ', m_gpar , m_crate5)
    modparid(129) = modparidtype('kcgwreg   ', m_gpar , m_crate6)
    modparid(130) = modparidtype('rips      ', m_gpar , m_rips)
    modparid(131) = modparidtype('sdnsnew   ', m_gpar , m_sndens0)
    modparid(132) = modparidtype('pprelmax  ', m_gpar , m_pprelmax)
    modparid(133) = modparidtype('pprelexp  ', m_gpar , m_pprelexp)
    modparid(134) = modparidtype('tnmean    ', m_lrpar, m_tnmean)
    modparid(135) = modparidtype('tocmean   ', m_lrpar, m_tocmean)
    modparid(136) = modparidtype('humusc2   ', m_lpar , m_humusC2)
    modparid(137) = modparidtype('humusc3   ', m_lpar , m_humusC3)
    modparid(138) = modparidtype('fastc2    ', m_lpar , m_fastC2)
    modparid(139) = modparidtype('fastc3    ', m_lpar , m_fastC3)
    modparid(n_rrcsc) = modparidtype('rrcscorr  ', m_rpar , m_rrcscorr)
    modparid(141) = modparidtype('kof       ', m_gpar , m_crate9)
    modparid(142) = modparidtype('koflim    ', m_gpar , m_crate10)
    modparid(143) = modparidtype('grata     ', m_gpar , m_grat3)
    modparid(144) = modparidtype('limqprod  ', m_gpar , m_limprod)
    modparid(145) = modparidtype('partp1    ', m_lpar , m_partP1)
    modparid(146) = modparidtype('partp2    ', m_lpar , m_partP2)
    modparid(147) = modparidtype('partp3    ', m_lpar , m_partP3)
    modparid(148) = modparidtype('wprodc    ', m_gpar , m_wprodc)
    modparid(149) = modparidtype('denitwrl  ', m_gpar,  m_denitwrl) !denitw local river
    modparid(150) = modparidtype('qmean     ', m_ldpar, m_ldqmean)
    modparid(151) = modparidtype('tpmean    ', m_ldpar, m_ldtpmean)
    modparid(152) = modparidtype('tnmean    ', m_ldpar, m_ldtnmean)
    modparid(153) = modparidtype('tocmean   ', m_ldpar, m_ldtocmean)
    modparid(154) = modparidtype('wprodn    ', m_ldpar, m_ldwprodn) 
    modparid(155) = modparidtype('sedon     ', m_ldpar, m_ldsedon)
    modparid(156) = modparidtype('sedoc     ', m_ldpar, m_ldsedoc)
    modparid(157) = modparidtype('sedpp     ', m_ldpar, m_ldsedpp)
    modparid(158) = modparidtype('wprodp    ', m_ldpar, m_ldwprodp) 
    modparid(159) = modparidtype('limt2exch ', m_gpar , m_limt2exch)
    modparid(160) = modparidtype('deeplake  ', m_gpar , m_deeplake)  
    modparid(161) = modparidtype('denitwl   ', m_gpar , m_denitwl)   
    modparid(162) = modparidtype('denitwl   ', m_ldpar, m_lddenitwl) 
    modparid(163) = modparidtype('wprodc    ', m_ldpar, m_ldwprodc)  
    modparid(164) = modparidtype('prodpp    ', m_ldpar, m_ldprodpp)  
    modparid(165) = modparidtype('prodsp    ', m_ldpar, m_ldprodsp)  
    modparid(166) = modparidtype('deeplake  ', m_ldpar, m_lddeeplake)
    modparid(167) = modparidtype('fastlake  ', m_ldpar, m_ldfastlake)
    modparid(168) = modparidtype('laketemp  ', m_gpar , m_laketemp)  
    modparid(169) = modparidtype('deadm     ', m_gpar , m_deadm)    
    modparid(170) = modparidtype('incorr    ', m_wqrpar, m_incorr) 
    modparid(171) = modparidtype('oncorr    ', m_wqrpar, m_oncorr) 
    modparid(172) = modparidtype('phoscorr  ', m_wqrpar, m_phoscorr)  
    modparid(173) = modparidtype('ratcorr   ', m_rpar , m_ratcorr)   
    modparid(174) = modparidtype('ponatm    ', m_lpar , m_ponatm)    
    modparid(175) = modparidtype('preccorr  ', m_rpar , m_preccorr)
    modparid(176) = modparidtype('cirrsink  ', m_rpar , m_cirrsink)  
    modparid(177) = modparidtype('irrcomp   ', m_gpar , m_irrcomp)   
    modparid(178) = modparidtype('ocsoimsat ', m_lpar , m_ocsoim)     
    modparid(179) = modparidtype('ocsoimslp ', m_lpar , m_ocsmslp)     
    modparid(180) = modparidtype('pirrs     ', m_rpar , m_pirrs)     
    modparid(181) = modparidtype('pirrg     ', m_rpar , m_pirrg)     
    modparid(182) = modparidtype('onpercred ', m_lpar , m_onpercred) 
    modparid(183) = modparidtype('pppercred ', m_lpar , m_pppercred) 
    modparid(184) = modparidtype('onconc0   ', m_lpar , m_onconc0)   
    modparid(185) = modparidtype('ppconc0   ', m_lpar , m_ppconc0)   
    modparid(186) = modparidtype('occonc0   ', m_lpar , m_occonc0)  
    modparid(187) = modparidtype('snalbmin  ', m_lpar , m_snalbmin) 
    modparid(188) = modparidtype('snalbmax  ', m_lpar , m_snalbmax) 
    modparid(189) = modparidtype('snalbkexp ', m_lpar , m_snalbkexp)
    modparid(190) = modparidtype('cmrad     ', m_lpar , m_cmrad)  
    modparid(191) = modparidtype('pcluse    ', m_lpar , m_pcluse)
    modparid(192) = modparidtype('aloadconst', m_gpar , m_atmload)
    
    !Water T2 temperature parameters
    modparid(193) = modparidtype('t2trriver ', m_gpar , m_t2trriver)  !temp flow from air to river
    modparid(194) = modparidtype('t2trlake  ', m_gpar , m_t2trlake)   !temp flow from air to lake
    modparid(195) = modparidtype('krelflood ', m_gpar , m_krelflood)  !some flood control dam parameter
    modparid(196) = modparidtype('upper2deep', m_gpar , m_upper2deep) !temp flow from upper to lower lake layer

    !Lake ice parameters
    modparid(197) = modparidtype('licewme   ', m_gpar , m_licewme)    !lake ice, water melt efficiency
    modparid(198) = modparidtype('licetf    ', m_gpar , m_licetf)     !lake ice, freezing temperature     
    modparid(199) = modparidtype('licesndens', m_gpar , m_licesndens) !lake ice, snow compaction parameter
    modparid(200) = modparidtype('licekika  ', m_gpar , m_licekika)   !lake ice, ki/ka     
    modparid(201) = modparidtype('licekexp  ', m_gpar , m_licekexp)   !lake ice, ks = ki*(dsnow(dice)^kexp, 1.88 
    modparid(202) = modparidtype('licetmelt ', m_gpar , m_licetmelt)  !lake ice, degreeday factor for ice melt
    modparid(203) = modparidtype('licewcorr ', m_gpar , m_licewcorr)  !lake ice, snow fall reduction for winddrift 

    !River ice parameters
    modparid(204) = modparidtype('ricewme   ', m_gpar , m_ricewme)    !river ice, water melt efficiency
    modparid(205) = modparidtype('ricetf    ', m_gpar , m_ricetf)     !river ice, freezing temperature     
    modparid(206) = modparidtype('ricesndens', m_gpar , m_ricesndens) !river ice, snow compaction parameter
    modparid(207) = modparidtype('ricekika  ', m_gpar , m_ricekika)   !river ice, ki/ka     
    modparid(208) = modparidtype('ricekexp  ', m_gpar , m_ricekexp)   !river ice, ki/ks
    modparid(209) = modparidtype('ricetmelt ', m_gpar , m_ricetmelt)  !river ice, degreeday factor for ice melt
    
    !snow cover area parameters
    modparid(210) = modparidtype('fscmax    ', m_gpar, m_fscmax)      !maximum snow cover area (0.95 [-])
    modparid(211) = modparidtype('fscmin    ', m_gpar, m_fscmin)      !minimum fsc             (0.001 [-])
    modparid(212) = modparidtype('fsclim    ', m_gpar, m_fsclim)      !fsc limit for onset of snowmax (0.001 [-])
    modparid(213) = modparidtype('fscdistmax', m_lpar, m_fscdistmax)  !maximum snow distribution factor (0.8 [mm^-1])
    modparid(214) = modparidtype('fscdist0  ', m_lpar, m_fscdist0)    !minimum snow distribution factor (0.6 [mm^-1])
    modparid(215) = modparidtype('fscdist1  ', m_lpar, m_fscdist1)    !elev_std coefficient for snow distr factor (0.001 [m^-1])
    modparid(216) = modparidtype('fsck1     ', m_gpar, m_fsck1)       !time constant in decrease of snowmax during melt (0.2 [-])
    modparid(217) = modparidtype('fsckexp   ', m_gpar, m_fsckexp)     !exponential coefficient in decrease of snowmax during melt (1e-6 [s^-1])
    
    !parameters for radiation and optional potential evaporation calculations
    modparid(218) = modparidtype('krs       ', m_gpar, m_krs)         !Hargreaves adjustment factor in estimation of shortwave radiation from tmin and tmax 0.16 (inland) 0.19 (coastal)
    modparid(219) = modparidtype('jhtadd    ', m_gpar, m_jhtadd)      !PET parameter in Jensen-Haise/McGuinness following Oudin et al (2005), recommended value 5
    modparid(220) = modparidtype('jhtscale  ', m_gpar, m_jhtscale)    !PET parameter in Jensen-Haise/McGuinness following Oudin et al (2005), recommended value 100
    modparid(221) = modparidtype('alfapt    ', m_gpar, m_alfapt)      !Priestly-taylor coefficient, recommended value 1.26
    modparid(222) = modparidtype('mwind     ', m_gpar, m_mwind)       !Mean windspeed, to replace missing data in FAO Penman-Monteith pet-equation, suggested value 2 m/s
    modparid(223) = modparidtype('kc        ', m_lpar, m_kc)          !Landuse dependent crop coefficient, used to scale the optional reference PET estimates (Hargreaves, Jensen, Priestly-Taylor,Penman-Monteith....)
    modparid(224) = modparidtype('alb       ', m_lpar, m_alb)         !Landuse dependent albedo, used for net radiation calculation (Priestly-Taylor,Penman-Monteith), suggested value 0.23

    !Glacier parameters (previous non-optional constants in the code)
    modparid(225) = modparidtype('glacvcoef ', m_gpar, m_glacvcoef)   !Coefficient glacier volume-area relationship (default 0.205)
    modparid(226) = modparidtype('glacvexp  ', m_gpar, m_glacvexp)    !Exponent glacier volume-area relationship (default 1.375)
    modparid(227) = modparidtype('glacdens  ', m_gpar, m_glacdens)    !Glacier density (default 0.85 m3 water/m3 ice)
    modparid(228) = modparidtype('glacvcoef1', m_gpar, m_glacvcoef1)  !Coefficient glacier volume-area relationship, type 1 glacier (default 1.701)
    modparid(229) = modparidtype('glacvexp1 ', m_gpar, m_glacvexp1)   !Exponent glacier volume-area relationship, type 1 glacier (default 1.25)
    modparid(230) = modparidtype('glac2arlim', m_gpar, m_glac2arlim)  !Area limit to separate glacer type 1 (Glaciers or Small) from type 2 (Ice Caps or Large)

    modparid(231) = modparidtype('rcgrwst   ', m_spar, m_rcgrwst)     !Deep percolation to aquifer 
    modparid(232) = modparidtype('aqretcor  ', m_rpar, m_aqretcorr)   !Aquifer return flow adjustment parameter
    modparid(233) = modparidtype('aqdelcor  ', m_rpar, m_aqdelcorr)   !Aquifer percolation delay adjustment parameter
    modparid(234) = modparidtype('aqpercor  ', m_rpar, m_aqpercorr)   !Aquifer percolation adjustment parameter
    modparid(235) = modparidtype('denitaq   ', m_gpar, m_denitaq)     !Aquifer denitrification 
    
    !Additional Lake and River water temperature parameters, to be replacing the t2trlake and t2trriver parameters
    modparid(236) = modparidtype('tcfriver  ', m_gpar, m_tcfriver)    !air-riverwater heat flow, temperature difference coefficient
    modparid(237) = modparidtype('scfriver  ', m_gpar, m_scfriver)    !air-riverwater heat flow, solar radiation coefficient
    modparid(238) = modparidtype('ccfriver  ', m_gpar, m_ccfriver)    !air-riverwater heat flow, constant coefficient
    modparid(239) = modparidtype('lcfriver  ', m_gpar, m_lcfriver)    !air-riverwater heat flow, linear coefficient
    modparid(240) = modparidtype('tcflake   ', m_gpar, m_tcflake)      !air-lakewater heat flow, temperature difference coefficient
    modparid(241) = modparidtype('scflake   ', m_gpar, m_scflake)      !air-lakewater heat flow, solar radiation coefficient
    modparid(242) = modparidtype('ccflake   ', m_gpar, m_ccflake)      !air-lakewater heat flow, constant coefficient
    modparid(243) = modparidtype('lcflake   ', m_gpar, m_lcflake)      !air-lakewater heat flow, linear coefficient
    modparid(244) = modparidtype('stbcorr1  ', m_gpar, m_stbcorr1)    !parameter for stability correction
    modparid(245) = modparidtype('stbcorr2  ', m_gpar, m_stbcorr2)    !parameter for stability correction
    modparid(246) = modparidtype('stbcorr3  ', m_gpar, m_stbcorr3)    !parameter for stability correction
    modparid(247) = modparidtype('zwind     ', m_gpar, m_zwind)       !wind observation level
    modparid(248) = modparidtype('zwish     ', m_gpar, m_zwish)       !wanted wind observation level for PET
    modparid(249) = modparidtype('zpdh      ', m_gpar, m_zpdh)        !zero plane displacement height for wind
    modparid(250) = modparidtype('roughness ', m_gpar, m_roughness)   !surface roughness for wind
    
    !Additional precipitation correction; orographic (pcelevstd) and phase (pcusnow, pcurain) impact on undercatch
    modparid(251) = modparidtype('pcelevstd ', m_gpar, m_pcelevstd)    !fractional increase in precipitation per 100m elev.std
    modparid(n_pcur) = modparidtype('pcurain   ', m_gpar, m_pcurain)      !undercatch correction factor for rain
    modparid(n_pcus) = modparidtype('pcusnow   ', m_gpar, m_pcusnow)      !undercatch correction factor for snow

    !Snow evaporation, refreeze, liquid water content holding capacity
    modparid(254) = modparidtype('fepotsnow ', m_lpar, m_fepotsnow)   !fraction of potential evaporation used for snow evaporation
    modparid(255) = modparidtype('cmrefr    ', m_gpar, m_cmrefr)      !snow refreeze efficiency (fraction of degree day factor cmlt)
    modparid(256) = modparidtype('fsceff    ', m_gpar, m_fsceff)      !efficiency of fractional snow cover to reduce melt and evap
    
    !Separate glacier melt parameters and sublimation/evaporation parameters
    modparid(257) = modparidtype('glacalb   ', m_gpar, m_glacalb)           !glacier ice albedo (0.35)
    modparid(258) = modparidtype('glacttmp  ', m_gpar, m_glacttmp)
    modparid(259) = modparidtype('glaccmlt  ', m_gpar, m_glaccmlt)
    modparid(260) = modparidtype('glaccmrad ', m_gpar, m_glaccmrad)
    modparid(261) = modparidtype('glaccmrefr', m_gpar, m_glaccmrefr)
    modparid(262) = modparidtype('fepotglac ', m_gpar, m_fepotglac)

    modparid(263) = modparidtype('kthrflood ', m_gpar, m_kthrflood)      !Threshold inflow over which a flood control dam save water (fraction of max inflow)
    modparid(264) = modparidtype('klowflood ', m_gpar, m_klowflood)      !Threshold level for extra flood control releases (fraction of regvol, typical 1/3)
    modparid(265) = modparidtype('monthlapse', m_mpar, m_mlapse)
    modparid(266) = modparidtype('limsedon  ', m_gpar, m_limsedon)
    modparid(267) = modparidtype('limsedpp  ', m_gpar, m_limsedpp)
    
    modparid(268) = modparidtype('opt1      ', m_gpar, m_opt1)        ! special optimisation parameter (to replace a GeoData or LakeData)
    modparid(269) = modparidtype('opt2      ', m_gpar, m_opt2)        ! special optimisation parameter (to replace a GeoData or LakeData)
    modparid(270) = modparidtype('opt3      ', m_gpar, m_opt3)        ! special optimisation parameter (to replace a GeoData or LakeData)
    modparid(271) = modparidtype('opt4      ', m_gpar, m_opt4)        ! special optimisation parameter (to replace a GeoData or LakeData)
    modparid(272) = modparidtype('opt5      ', m_gpar, m_opt5)        ! special optimisation parameter (to replace a GeoData or LakeData)
    modparid(273) = modparidtype('opt6      ', m_gpar, m_opt6)        ! special optimisation parameter (to replace a GeoData or LakeData)
    modparid(274) = modparidtype('opt7      ', m_gpar, m_opt7)        ! special optimisation parameter (to replace a GeoData or LakeData)
    modparid(275) = modparidtype('opt8      ', m_gpar, m_opt8)        ! special optimisation parameter (to replace a GeoData or LakeData)
    modparid(276) = modparidtype('optonoff  ', m_gpar, m_optonoff)    ! parameter "switch" to turn the opt1-opt8 on(1)/off(0) for the floodplainmodel
    modparid(277) = modparidtype('init1sw   ', m_gpar, m_init1sw)     !T1 initial concentration in river and lakes
    modparid(278) = modparidtype('sdnsmax   ', m_gpar, m_sdnsmax)     !snow density/depth parameter
    modparid(279) = modparidtype('sdnsrate  ', m_gpar, m_sdnsrate)    !snow density/depth parameter
    modparid(280) = modparidtype('sdnsradd  ', m_gpar, m_sdnsradd)    !snow density/depth parameter
    modparid(281) = modparidtype('t2mix     ', m_gpar, m_gt2mix)      !T2 of outflow parameter
    modparid(282) = modparidtype('t2mix     ', m_ldpar, m_ldt2mix)    !T2 of outflow parameter
    modparid(283) = modparidtype('denitrlu3 ', m_lpar , m_denitr3)    !denitrification soillayer 3
    modparid(284) = modparidtype('hsatins   ', m_gpar , m_hsatINsoil) !half saturation concentration, denitrification in soil (mg N/L)
    modparid(285) = modparidtype('hsatinw   ', m_gpar , m_hsatINwater) !half saturation concentration, denitrification in water (mg N/L)
    modparid(286) = modparidtype('hsattp    ', m_gpar , m_hsatTP)     !half saturation concentration, production in water (mg P/L)
    !T1 parameters
    modparid(287) = modparidtype('t1expdec  ', m_gpar, m_expdec)      !half-life of tracer
    modparid(288) = modparidtype('t1freuc   ', m_gpar, m_t1freuc)     !tracer adsorption coefficient
    modparid(289) = modparidtype('t1rel     ', m_gpar, m_t1rel)       !tracer release from source parameter (mm-1)
    modparid(290) = modparidtype('t1sedvel  ', m_gpar, m_t1sed)       !sedimentation velocity of T1 [m/timestep]
    modparid(291) = modparidtype('t1sedexp  ', m_gpar, m_t1sedexp)    !sedimentation/resuspension of T1 in river sediments
    modparid(292) = modparidtype('t1leaksoil', m_spar, m_t1leakst)    !T1 concentration of leakage
    modparid(293) = modparidtype('t1leakluse', m_lpar, m_t1leaklu)    !T1 concentration of leakage
    modparid(294) = modparidtype('denit3reg ', m_wqrpar, m_denit3reg)    !denitrification soillayer 3, regional replacement


  END SUBROUTINE initiate_model_parameters

!>Set flags for special HYPE models that will be used
!------------------------------------------------------------
  SUBROUTINE set_special_models(nc,arr,glacexist,irrexist)

  USE IRRIGATION_MODULE, ONLY : irrtype,  &
                                check_for_irrigated_classes
  !Argument declarations
  INTEGER,INTENT(IN) :: nc         !<dimension, number of classes
  INTEGER,INTENT(IN) :: arr(nc)    !<soilmodels
  LOGICAL, INTENT(OUT) :: glacexist   !<status of glacier model
  LOGICAL, INTENT(OUT) :: irrexist    !<status of irrigation model

  INTEGER i

    !Check for glacier soilmodel
    glacexist = .FALSE.
    DO i=1,nc
      IF(arr(i)==glacier_model) glacexist=.TRUE.
    ENDDO
  
    !Check for irrigation included in model set-up
    irrexist = .TRUE.
    IF(.NOT.ALLOCATED(irrigationsystem))THEN
      irrexist = .FALSE.      !Irrigation information not provided in MgmtData
    ELSEIF(.NOT.ALLOCATED(cropirrdata))THEN
      irrexist = .FALSE.      !Irrigation information not provided in CropData
    ELSE
      IF(.NOT.ALLOCATED(irrtype))THEN
        ALLOCATE(irrtype(nc))
        CALL check_for_irrigated_classes(nc)
      ENDIF
      IF(SUM(irrtype)==0)THEN
        irrexist = .FALSE.
        DEALLOCATE(irrtype)
      ENDIF
    ENDIF

  END SUBROUTINE set_special_models
  
!>Set special HYPE parameters needed by HYSS
!------------------------------------------------------------
  SUBROUTINE get_special_model_parameters(velindex,dampindex)

  !Argument declarations
  INTEGER, INTENT(OUT) :: velindex  !<index of rivvel in modparid
  INTEGER, INTENT(OUT) :: dampindex !<index of damp in modparid

    !Index of HYPE parameters are needed by HYSS
    velindex = n_rivv
    dampindex = n_damp

  END SUBROUTINE get_special_model_parameters
  
!>Calculate special HYPE parameters needed by HYSS
!------------------------------------------------------------
  SUBROUTINE calculate_special_model_parameters(nsubb,optrivvel,optdamp,dimriver)

  USE SURFACEWATER_PROCESSES, ONLY : calculate_landarea_riverlength
  USE HYPE_INDATA, ONLY : clust_group,catdes, &
                          ndes_reg, idx_catdes, wght_catdes, indx_par
  
  !Argument declarations
  INTEGER, INTENT(IN)  :: nsubb     !<dimension, number of subbasins (submodel)
  REAL, INTENT(IN)     :: optrivvel !<lower river velocity boundary (optpar)
  REAL, INTENT(IN)     :: optdamp   !<lower damp boundary (optpar)
  INTEGER, INTENT(OUT) :: dimriver  !<maximum size river lag needed

  INTEGER i,k,itype
  REAL maxtrans
  REAL landarea(nsubb)
  REAL rivlength(2,nsubb)   !river length (m)
  REAL totaltime,transtime  !total time in river and time in river train (translation) (in timesteps)
  REAL localrivvel,localdamp

    !Calculate local and main river length
    CALL calculate_landarea_riverlength(nsubb,landarea,rivlength)
    
    !Calculate river translation time, and the maximum, for every subbasin
    maxtrans = 0.
    localrivvel = genpar(m_rivvel)
    localdamp = genpar(m_damp)
    IF(.NOT.conductregest)THEN
      localrivvel = MIN(localrivvel,optrivvel)    !This will give the "worst" combination
      localdamp = MIN(localdamp,optdamp)
      IF(localrivvel==0)THEN !Test for missing rivvel in model set-up
        dimriver = 1
        WRITE(6,*) 'Warning: Parameter rivvel not set. River flow will take one time step.'
        RETURN
      ENDIF
      totaltime = MAXVAL(rivlength) / localrivvel / seconds_per_timestep
      transtime = (1. - localdamp) * totaltime
      dimriver = INT(transtime) + 1
    ELSE
      !If regional parameter estimations is used calculate rivvel and damp
      DO i = 1,nsubb
        IF(indx_par(n_rivv).gt.0)THEN
          localrivvel=0.0
          DO k = 1,ndes_reg(clust_group(i),indx_par(n_rivv))
            localrivvel = localrivvel + wght_catdes(clust_group(i),indx_par(n_rivv),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_rivv),k))
          ENDDO
        ENDIF
        IF(indx_par(n_damp).gt.0) THEN
          localdamp=0.0
          DO k = 1,ndes_reg(clust_group(i),indx_par(n_damp))
            localdamp = localdamp + wght_catdes(clust_group(i),indx_par(n_damp),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_damp),k))
          ENDDO
        ENDIF   
        localrivvel = MIN(localrivvel,optrivvel)    !This will give the "worst" combination
        localdamp = MIN(localdamp,optdamp)
        IF(localrivvel==0)THEN !Test for missing rivvel in model set-up
          dimriver = 1
          WRITE(6,*) 'Warning: Parameter rivvel not set. River flow will take one time step.'
          RETURN
        ENDIF
        DO itype = 1,2          
          totaltime = rivlength(itype,i) / localrivvel / seconds_per_timestep
          transtime = (1. - localdamp) * totaltime
          maxtrans = MAX(maxtrans,transtime)
        ENDDO
      ENDDO
      dimriver = INT(maxtrans) + 1
    ENDIF
    
  END SUBROUTINE calculate_special_model_parameters
  
  !>Initiates the model state for a simulation with no saved states. 
  !-----------------------------------------
  SUBROUTINE initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    USE GLACIER_SOILMODEL, ONLY :   initiate_glacier_state
    USE SOIL_PROCESSES, ONLY :      initiate_soil_water_state
    USE NPC_SOIL_PROCESSES, ONLY :  initiate_soil_npc_state
    USE SURFACEWATER_PROCESSES, ONLY : sum_upstream_area, &
                                       calculate_landarea_riverlength
    USE NPC_SURFACEWATER_PROCESSES, ONLY : initiate_river_npc_state,     &
                                           initiate_lake_npc_state
    USE REGIONAL_GROUNDWATER_MODULE, ONLY : initiate_aquifer_state

    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    
    !Variable declarations
    INTEGER i     !loop-variable

    !Parameter declarations
    REAL, PARAMETER :: seconds_per_day = 86400.  

    !Existence of soillayers
    IF(MAXVAL(soilthick(3,1:nclass))>0.)THEN
      nummaxlayers = 3
    ELSEIF(MAXVAL(soilthick(2,1:nclass))>0.)THEN
      nummaxlayers = 2
    ELSE
      nummaxlayers = 1
    ENDIF
    
    !Initiate states to zero
    CALL initiate_state_zero(numsubstances,naquifers,i_t1>0,i_t2>0,wetlandexist, &
                             glacierexist,modeloption(p_lakeriverice)>=1, &
                             doirrigation,doupdate(i_qar).OR.doupdate(i_war), &
                             conductflood,modeloption(p_growthstart)==1, &
                             frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    !Initiate soil water state variables and soil water processes help parameters.
    CALL initiate_soil_water_state(soilstate)
    
    !Initialize glacier volume
    IF(glacierexist) CALL initiate_glacier_state(nclass,classmodel,frozenstate)
    
    !Initialize lake state for pure water model simulation to total lake volume
    DO i = 1,nsub
      IF(slc_ilake>0) lakestate%water(1,i) = genpar(m_gldepi)*1000.  !ilake water stage (mm)
      IF(slc_olake>0)THEN
        lakestate%water(2,i) = basin(i)%lakedepth * 1000.            !ordinary olake water stage (mm)
        IF(ALLOCATED(damindex))THEN
          IF(damindex(i)>0)THEN
            IF(dam(damindex(i))%purpose==3) THEN
!              lakestate%water(2,i) = 0                               !dam for flood control start empty
              lakestate%water(2,i) = lakestate%water(2,i) - dam(damindex(i))%regvol*1000000./(classbasin(i,slc_olake)%part*basin(i)%area)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    
    !Initialize river states
    IF(.NOT. ALLOCATED(landarea))     ALLOCATE(landarea(nsub))
    IF(.NOT. ALLOCATED(riverlength))  ALLOCATE(riverlength(2,nsub))
    IF(.NOT. ALLOCATED(deadriver))    ALLOCATE(deadriver(2,nsub))
    IF(.NOT. ALLOCATED(upstreamarea)) ALLOCATE(upstreamarea(nsub))
    CALL sum_upstream_area(nsub,upstreamarea)             !calculates the area upstream of the outlet point of catchment
    CALL calculate_landarea_riverlength(nsub,landarea,riverlength)
    DO i = 1,nsub                     !initiate lake water stage for pure water model simulation
      deadriver(1,i) = genpar(m_deadl) * (basin(i)%area/1.0E6) * riverlength(1,i)
      deadriver(2,i) = genpar(m_deadm) * (upstreamarea(i)/1.0E6) * riverlength(2,i)
      riverstate%water(:,i) = deadriver(:,i)    !initialize river damping box to deadvolume (m3)
    ENDDO
    
    !Initiate average discharge variable, Qmeani = 365 day mean Q before damping  
    DO i = 1,nsub
      IF(ALLOCATED(lakedatapar))THEN
        riverstate%Qmean(1,i) = lakedatapar(lakedataparindex(i,1),m_ldqmean)* 0.001 * basin(i)%area / (365. * seconds_per_day)
        riverstate%Qmean(2,i) = lakedatapar(lakedataparindex(i,2),m_ldqmean)* 0.001 * upstreamarea(i) / (365. * seconds_per_day)
      ELSE
        riverstate%Qmean(1,i) = genpar(m_Qmean)* 0.001 * basin(i)%area / (365. * seconds_per_day)
        riverstate%Qmean(2,i) = genpar(m_Qmean)* 0.001 * upstreamarea(i) / (365. * seconds_per_day)
      ENDIF
    ENDDO
   

    !Initialize soil nutrient variables
    CALL initiate_soil_npc_state(conductN,conductP,conductC,nummaxlayers,soilstate)
    
    !Initialize soil for tracer concentration
    IF(i_t1>0)THEN
      soilstate%conc(i_t1,1,:,:)  = genpar(m_iniT1)
      IF(nummaxlayers>1) soilstate%conc(i_t1,2,:,:) = genpar(m_iniT1)
      IF(nummaxlayers==3) soilstate%conc(i_t1,3,:,:) = genpar(m_iniT1)
    ENDIF
    IF(i_t2>0)THEN
      soilstate%conc(i_t2,1,:,:)  = genpar(m_iniT2)
      soilstate%temp(1,:,:) = genpar(m_iniT2)
      soilstate%deeptemp(:,:)    = genpar(m_iniT2)
      IF(nummaxlayers>1)THEN
        soilstate%conc(i_t2,2,:,:) = genpar(m_iniT2)
        soilstate%temp(2,:,:) = genpar(m_iniT2) 
      ENDIF
      IF(nummaxlayers==3)THEN
        soilstate%conc(i_t2,3,:,:) = genpar(m_iniT2)
        soilstate%temp(3,:,:) = genpar(m_iniT2)
      ENDIF
    ENDIF
    
    !Initiate river variables      
    CALL initiate_river_npc_state(conductN,conductP,conductC,i_t1>0,riverstate)
    
    !Initialize lake concentrations
    CALL initiate_lake_npc_state(conductN,conductP,conductC,conductT,slc_ilake>0,slc_olake>0,lakestate)
    
    !Initialize surface water for tracer concentration (lakestate T2 is set in initiate_lake_npc_state)
    IF(i_t1>0)THEN
      riverstate%conc(i_t1,:,:)  = genpar(m_iniT1sw)
      riverstate%cqueue(i_t1,:,:,:)  = genpar(m_iniT1sw)
      lakestate%conc(i_t1,:,:) = genpar(m_iniT1sw)
      lakestate%concslow(i_t1,:,:) = genpar(m_iniT1sw)
    ENDIF
    IF(i_t2>0) riverstate%cqueue(i_t2,:,:,:) = genpar(m_iniT2)
    !Varför saknas riverstate%conc(i_t2?

    !Initiate wetland variable
    IF(wetlandexist)THEN
      IF(i_t2>0) riverstate%cwetland(i_t2,:,:)=genpar(m_iniT2)
      !Varför saknas iniT1?
    ENDIF

    !Initiate aquifer state variables
    IF(modeloption(p_deepgroundwater)==2) CALL initiate_aquifer_state(aquiferstate)
    
  END SUBROUTINE initiate_model_state

  !>\brief Initiates model variables and parameters for a simulation. 
  !>This include 
  !>calculating HYPE variables for holding model set-up information and more...
  !--------------------------------------------------------------------------
  SUBROUTINE initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    USE GLACIER_SOILMODEL, ONLY :   initiate_glacier
    USE SOIL_PROCESSES, ONLY : initiate_soil_water
    USE SURFACEWATER_PROCESSES, ONLY : sum_upstream_area, &
                                       set_general_rating_k,  &
                                       calculate_landarea_riverlength, &
                                       set_lake_outlets
    USE NPC_SURFACEWATER_PROCESSES, ONLY : set_lake_slowwater_maxvolume
    USE IRRIGATION_MODULE, ONLY :   initiate_irrigation
    USE REGIONAL_GROUNDWATER_MODULE, ONLY : initiate_regional_groundwater_flow, &
                                            initiate_aquifer_model, &
                                            calculate_delayed_water
    USE ATMOSPHERIC_PROCESSES, ONLY : calculate_class_wind_transformation_factor, &
                                      set_atmospheric_parameters_corrections
    USE HYPE_INDATA, ONLY : clust_group,catdes, & 
                            ndes_reg, idx_catdes, wght_catdes, indx_par, &
                            deallocate_regest_input_variables
    
    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    
    !Variable declarations
    INTEGER i,j,k     !loop-variable
    INTEGER itype     !loop-variable, river type
    REAL help
    REAL vel,part,totaltime,kt  !parameters and help variables for river transport
    REAL delayedwater(naquifers)        !water in delay to reach aquifer (m3)
    REAL localrivvol(nsub)  !help variable for calculating local river initial volume
    REAL mainrivvol(nsub)   !help variable for calculating main river initial volume

    !Parameter declarations
    REAL, PARAMETER :: seconds_per_day = 86400.  

    
    !Flag for radiation and humidity calculations
    IF(modeloption(p_snowmelt).GE.2 .OR. modeloption(p_petmodel).GT.1 .OR. modeloption(p_lakeriverice).GE.1)THEN
      calcSWRAD = .TRUE.
    ELSE
      calcSWRAD = .FALSE.
    ENDIF
    IF(modeloption(p_petmodel).GT.1)THEN
      calcVAPOUR = .TRUE.
    ELSE
      calcVAPOUR = .FALSE.
    ENDIF
    IF(modeloption(p_petmodel).EQ.5)THEN
      calcWIND = .TRUE.
      IF(ALLOCATED(windi)) CALL calculate_class_wind_transformation_factor(windtrans)
    ELSE
      calcWIND = .FALSE.
    ENDIF

    !Set flag for T1 simulation with typical soil leakage
    T1leakage = .FALSE.
    IF(i_t1>0)THEN
      help = SUM(landpar(m_t1leaklu,:))+SUM(soilpar(m_t1leakst,:))
      IF(help>0.)THEN
        T1leakage = .TRUE.
      ENDIF
    ENDIF

    !Set precipitation parameter corrections
    CALL set_atmospheric_parameters_corrections()
    
    !Existence of soillayers
    IF(MAXVAL(soilthick(3,1:nclass))>0)THEN
      nummaxlayers = 3
    ELSEIF(MAXVAL(soilthick(2,1:nclass))>0)THEN
      nummaxlayers = 2
    ELSE
      nummaxlayers = 1
    ENDIF
    
    !Set potential evaporation correction (not glacier?)
    IF(.NOT.ALLOCATED(basincevpcorr)) ALLOCATE(basincevpcorr(nsub))
    basincevpcorr = 1.
    DO i = 1,nsub
      IF(basin(i)%parregion>0)THEN
        basincevpcorr(i) = 1. + regpar(m_cevpcorr,basin(i)%parregion)
      ENDIF
      !Replace parameter values with regional parameter estimates
      IF(conductregest)THEN
        IF(indx_par(n_cevpc).gt.0) THEN
          basincevpcorr(i) = 1.
          DO k=1,ndes_reg(clust_group(i),indx_par(n_cevpc))
            basincevpcorr(i) = basincevpcorr(i) + wght_catdes(clust_group(i),indx_par(n_cevpc),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_cevpc),k))
          ENDDO
        ENDIF
      ENDIF
    ENDDO
       
    !Initiate soil water processes help parameters.
    CALL initiate_soil_water()
    
    !Initialize glacier parameters and type
    IF(glacierexist) CALL initiate_glacier(nclass,classmodel)

    !Initiate soil temperature parameters and variables
    avertemp =(/5.,10.,20.,30./)*timesteps_per_day    !Number of timestep over which meantemp is calculated
    IF(.NOT.ALLOCATED(soilmem)) ALLOCATE(soilmem(maxsoillayers,nclass))
    soilmem = 0.
    DO j= 1,nclass
      DO k = 1,maxsoillayers
        IF(k>1)THEN
          soilmem(k,j) = timesteps_per_day*landpar(m_surfmem,classdata(j)%luse)*EXP(landpar(m_depthrel,classdata(j)%luse)*(soildepth(k-1,j)+(soilthick(k,j) / 2.)))
        ELSE  
          soilmem(k,j) = timesteps_per_day*landpar(m_surfmem,classdata(j)%luse)*EXP(landpar(m_depthrel,classdata(j)%luse)*(soilthick(k,j) / 2.))
        ENDIF
      ENDDO
    ENDDO
    
    !Initialize lake and river
    IF(.NOT. ALLOCATED(landarea))     ALLOCATE(landarea(nsub))
    IF(.NOT. ALLOCATED(riverlength))  ALLOCATE(riverlength(2,nsub))
    IF(.NOT. ALLOCATED(deadriver))    ALLOCATE(deadriver(2,nsub))
    IF(.NOT. ALLOCATED(upstreamarea)) ALLOCATE(upstreamarea(nsub))
    CALL sum_upstream_area(nsub,upstreamarea)             !calculates the area upstream of the outlet point of catchment
    CALL calculate_landarea_riverlength(nsub,landarea,riverlength)
    DO i = 1,nsub                     !initiate lake water stage for pure water model simulation
      deadriver(1,i) = genpar(m_deadl) * (basin(i)%area/1.0E6) * riverlength(1,i)
      deadriver(2,i) = genpar(m_deadm) * (upstreamarea(i)/1.0E6) * riverlength(2,i)
    ENDDO
    IF(.NOT. ALLOCATED(deadwidth))    ALLOCATE(deadwidth(2,nsub))
    IF(.NOT. ALLOCATED(ratingk))      ALLOCATE(ratingk(2,nsub))
    DO i = 1,nsub                     !initiate lake water stage for pure water model simulation
      deadwidth(1,i) = SQRT(genpar(m_deadl) * (basin(i)%area/1.0E6)*0.1)   !(m)  !width=10*depth, width = sqrt(area/10)
      deadwidth(2,i) = SQRT(genpar(m_deadm) * (upstreamarea(i)/1.0E6)*0.1)   !(m)  !width=10*depth, width = sqrt(area/10)
    ENDDO
    CALL set_general_rating_k(2,nsub,landarea,upstreamarea,ratingk)
    
    !No initiation of floodwater to zero here, have inital value from state-file.

    !Calculate river transport time parameters
    IF(.NOT. ALLOCATED(transtime)) ALLOCATE(transtime(2,nsub))
    IF(.NOT. ALLOCATED(ttstep))    ALLOCATE(ttstep(2,nsub))
    IF(.NOT. ALLOCATED(ttpart))    ALLOCATE(ttpart(2,nsub))
    IF(.NOT. ALLOCATED(riverrc))   ALLOCATE(riverrc(2,nsub))
    vel = genpar(m_rivvel)                          !peak velocity of river (m/s)
    IF((.NOT.conductregest).AND.vel==0)THEN
      transtime = 0.
      ttstep    = 0
      ttpart    = 0.
      riverrc   = 1.
    ELSE
      part = genpar(m_damp)                           !part of delay from damping
      DO i = 1,nsub
        !>Replace parameter values with regional parameter estimates
        IF(conductregest)THEN
          IF(indx_par(n_rivv).gt.0) THEN
            vel=0.0
            DO k=1,ndes_reg(clust_group(i),indx_par(n_rivv))
              vel = vel + wght_catdes(clust_group(i),indx_par(n_rivv),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_rivv),k))
            ENDDO
          ENDIF
          IF(indx_par(n_damp).gt.0) THEN
            part=0.0
            DO k=1,ndes_reg(clust_group(i),indx_par(n_damp))
              part = part + wght_catdes(clust_group(i),indx_par(n_damp),k)*catdes(i,idx_catdes(clust_group(i),indx_par(n_damp),k))
            ENDDO
          ENDIF
        ENDIF
        IF(vel<=0.001)THEN
          WRITE(6,*) 'Warning: river velocity parameter less than 0.001 m/s. Set to 0.001 m/s.'
          vel = 0.001
        ENDIF
        IF(part>1.)THEN
          WRITE(6,*) 'Warning: damping parameter larger than one. Set to one.'
          part = 1.
        ELSEIF(part<0.)THEN
          WRITE(6,*) 'Warning: damping parameter less than zero. Set to zero.'
          part = 0.
        ENDIF
        DO itype = 1,2          
          totaltime = riverlength(itype,i) / vel / seconds_per_timestep  !t=s/v, river length from GeoData (time step)
          transtime(itype,i) = (1. - part) * totaltime
          ttstep(itype,i) = INT(transtime(itype,i))
          ttpart(itype,i) = transtime(itype,i) - REAL(ttstep(itype,i))
          kt = part * totaltime                    !damptime, k in equation q=(1/k)*S
          IF(kt>0)THEN
            riverrc(itype,i) = 1. - kt + kt * exp(-1./kt)  !recession coefficient in equation Q=r*S
          ELSE
            riverrc(itype,i) = 1.   !Safe for all lake subbasin
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    !Check of dimension river queue
    help = SIZE(riverstate%qqueue,DIM=1) !number of elements (0:ml)
    !WRITE(6,*) 'Test: dim,ttstep',help,MAXVAL(ttstep)
    DO i = 1,nsub
      DO itype=1,2
        IF(help-1<ttstep(itype,i))THEN
          WRITE(6,*) 'Error: dimension queue',itype,i,ttstep(itype,i),help
        ENDIF
      ENDDO
    ENDDO
    
    !Set river bankful variables
    IF(conductN.OR.conductP) CALL set_Qvariables_for_bankfulflow(nsub,riverstate%Q365)

    !Initialize lake parameters
    IF(conductN.OR.conductP.OR.conductC.OR.conductT)  &
      CALL set_lake_slowwater_maxvolume(nsub,genpar(m_gldepi),basin(:)%lakedepth,lakedataparindex,lakedatapar(:,m_lddeeplake),slc_ilake>0,slc_olake>0)
    CALL set_lake_outlets()
 
    !Initiate regional groundwater flow/aquifer
    IF(modeloption(p_deepgroundwater)==1)THEN
      CALL initiate_regional_groundwater_flow(nsub,numsubstances,slc_ilake,slc_olake,slc_lriver,slc_mriver)
    ELSEIF(modeloption(p_deepgroundwater)==2)THEN
      CALL initiate_aquifer_model(nsub,numsubstances,naquifers)
    ENDIF
        
    !Allocate and initiate output variable
    IF(.NOT. ALLOCATED(accdiff)) ALLOCATE(accdiff(nsub))
    accdiff = 0.
    
    !Allocate nutrient (NP) load variables
    IF(conductload)THEN
       IF(.NOT.ALLOCATED(Latmdep))  ALLOCATE(Latmdep(nclass,2,numsubstances))   !Substances in order IN, ON, SP, PP
       IF(.NOT.ALLOCATED(Lcultiv))  ALLOCATE(Lcultiv(nclass,2,numsubstances))   !1=fertiliser, 2=plantdecay
       IF(.NOT.ALLOCATED(Lirrsoil)) ALLOCATE(Lirrsoil(nclass,numsubstances))    !irrigation on soil
       IF(.NOT.ALLOCATED(Lrurala))  ALLOCATE(Lrurala(nclass,numsubstances))     !rural a
       IF(.NOT.ALLOCATED(Lstream))  ALLOCATE(Lstream(nclass,numsubstances))     !runoff from soil to stream
       IF(.NOT.ALLOCATED(Lpathway)) ALLOCATE(Lpathway(numsubstances,13))        ! class independent, 13 is points A to M along flow pathway.
       IF(.NOT.ALLOCATED(Lbranch))  ALLOCATE(Lbranch(numsubstances))            !part of outflow
       IF(.NOT.ALLOCATED(Lgrwmr))   ALLOCATE(Lgrwmr(numsubstances))             !reg.grw to main river
       IF(.NOT.ALLOCATED(Lgrwol))   ALLOCATE(Lgrwol(numsubstances))             !reg.grw to olake
    ENDIF 
    IF(numsubstances>0)THEN
      IF(.NOT.ALLOCATED(Lruralb))  ALLOCATE(Lruralb(numsubstances))            !rural b
      IF(.NOT.ALLOCATED(Lpoints))  ALLOCATE(Lpoints(numsubstances,max_pstype)) !point sources 1=ps1, 2=ps2, 3=ps3, ..
      IF(.NOT.ALLOCATED(Lgrwsoil)) ALLOCATE(Lgrwsoil(nclass,numsubstances))       !regional groundwaterflow to soil
      IF(.NOT.ALLOCATED(Lgrwclass)) ALLOCATE(Lgrwclass(nclass,numsubstances,nsub))       !regional groundwater outflow from soil
    ENDIF
    
    !Initiate irrigation
    IF(doirrigation) CALL initiate_irrigation(nsub,nclass)

    !Initate calculations for water balance output
    IF(conductwb)THEN
      CALL initiate_waterbalance_output(nsub)
      !set initial wb_stores
      wbstores = 0.
      DO i = 1,nsub
        DO j = 1,nclass
          IF(classmodel(j)==glacier_model) wbstores(w_glacier,i) = frozenstate%glacvol(i)*genpar(m_glacdens)
          IF(classmodel(j)==0.OR.classmodel(j)==glacier_model)THEN
            wbstores(w_snow,i) = wbstores(w_snow,i) + frozenstate%snow(j,i)*classbasin(i,j)%part
            wbstores(w_soil1,i) = wbstores(w_soil1,i) + soilstate%water(1,j,i)*classbasin(i,j)%part
            wbstores(w_soil2,i) = wbstores(w_soil2,i) + soilstate%water(2,j,i)*classbasin(i,j)%part
            wbstores(w_soil3,i) = wbstores(w_soil3,i) + soilstate%water(3,j,i)*classbasin(i,j)%part
          ENDIF
          IF(conductflood)THEN
            IF(floodindex(i)>0)THEN
              IF(j==slc_mriver.AND.flooding(floodindex(i))%fpfmr>0.)THEN
                wbstores(w_snow,i) = wbstores(w_snow,i) + frozenstate%snow(j,i)*classbasin(i,j)%part*flooding(floodindex(i))%fpfmr
                wbstores(w_soil1,i) = wbstores(w_soil1,i) + soilstate%water(1,j,i)*classbasin(i,j)%part*flooding(floodindex(i))%fpfmr
                wbstores(w_soil2,i) = wbstores(w_soil2,i) + soilstate%water(2,j,i)*classbasin(i,j)%part*flooding(floodindex(i))%fpfmr
                wbstores(w_soil3,i) = wbstores(w_soil3,i) + soilstate%water(3,j,i)*classbasin(i,j)%part*flooding(floodindex(i))%fpfmr
              ENDIF
              IF(j==slc_olake.AND.flooding(floodindex(i))%fpfol>0.)THEN
                wbstores(w_snow,i) = wbstores(w_snow,i) + frozenstate%snow(j,i)*classbasin(i,j)%part*flooding(floodindex(i))%fpfol
                wbstores(w_soil1,i) = wbstores(w_soil1,i) + soilstate%water(1,j,i)*classbasin(i,j)%part*flooding(floodindex(i))%fpfol
                wbstores(w_soil2,i) = wbstores(w_soil2,i) + soilstate%water(2,j,i)*classbasin(i,j)%part*flooding(floodindex(i))%fpfol
                wbstores(w_soil3,i) = wbstores(w_soil3,i) + soilstate%water(3,j,i)*classbasin(i,j)%part*flooding(floodindex(i))%fpfol
              ENDIF
            ENDIF
          ENDIF
          IF(doirrigation) wbstores(w_irrcanal,i) = wbstores(w_irrcanal,i) + miscstate%nextirrigation(j,i)*classbasin(i,j)%part
        ENDDO
        localrivvol(i) = riverstate%water(1,i) + (SUM(riverstate%qqueue(1:ttstep(1,i),1,i)) + riverstate%qqueue(ttstep(1,i)+1,1,i) * ttpart(1,i))
        mainrivvol(i)  = riverstate%water(2,i) + (SUM(riverstate%qqueue(1:ttstep(2,i),2,i)) + riverstate%qqueue(ttstep(2,i)+1,2,i) * ttpart(2,i))
      ENDDO
      wbstores(w_snow,:)   = wbstores(w_snow,:)  * basin(:)%area *1.E-3  !m3
      wbstores(w_soil1,:)  = wbstores(w_soil1,:) * basin(:)%area *1.E-3  !m3
      wbstores(w_soil2,:)  = wbstores(w_soil2,:) * basin(:)%area *1.E-3  !m3
      wbstores(w_soil3,:)  = wbstores(w_soil3,:) * basin(:)%area *1.E-3  !m3
      IF(doirrigation) wbstores(w_irrcanal,:) = wbstores(w_irrcanal,:) * basin(:)%area *1.E-3  !m3
      wbstores(w_iriver,:) = localrivvol(:) !m3
      wbstores(w_mriver,:) = mainrivvol(:)  !m3
      IF(conductflood)THEN
        DO i = 1,nsub
          IF(floodindex(i)>0)THEN
            IF(flooding(floodindex(i))%fpfmr>0.)THEN
              wbstores(w_riverplain,i) = miscstate%floodwater(1,i)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      IF(slc_ilake>0) wbstores(w_ilake,:) = lakestate%water(1,:)  !mm
      IF(slc_olake>0) wbstores(w_olake,:)  = lakestate%water(2,:)  !mm
      IF(numsubstances>0)THEN
        IF(slc_ilake>0) wbstores(w_ilake,:) = wbstores(w_ilake,:) + lakestate%slowwater(1,:)  !mm
        IF(slc_olake>0) wbstores(w_olake,:) = wbstores(w_olake,:) + lakestate%slowwater(2,:)  !mm
      ENDIF
      IF(slc_ilake>0) wbstores(w_ilake,:) = wbstores(w_ilake,:) * basin(:)%area *classbasin(:,slc_ilake)%part *1.E-3  !m3
      IF(slc_olake>0)THEN
        wbstores(w_olake,:) = wbstores(w_olake,:) * basin(:)%area *classbasin(:,slc_olake)%part *1.E-3  !m3
        IF(conductflood)THEN
          DO i = 1,nsub
            IF(floodindex(i)>0)THEN
              IF(flooding(floodindex(i))%fpfol>0.)THEN
                wbstores(w_olake,i) = wbstores(w_olake,i)*(1.-flooding(floodindex(i))%fpfol)
                wbstores(w_lakeplain,i) = miscstate%floodwater(2,i)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      CALL calculate_delayed_water(aquiferstate,naquifers,delayedwater)
      IF(naquifers>0) wbstores(w_aquifer,1:naquifers) = aquiferstate%water + aquiferstate%nextoutflow + delayedwater
      CALL print_initial_waterbalance_stores(nsub,naquifers)
    ENDIF
    
    IF(conductregest) CALL deallocate_regest_input_variables()
    
!    CALL calculate_areaweights_for_upstream_average() !This should not always be calculated and maybe not in this subroutine either, but this should work

  END SUBROUTINE initiate_model
  
  !>Model subroutine for HYPE.
  !!Calculates what happen during one timestep.
  !---------------------------------------------------
  SUBROUTINE model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    
    USE STATETYPE_MODULE
    USE NPC_SURFACEWATER_PROCESSES, ONLY : add_dry_deposition_to_river,   &
                                           add_dry_deposition_to_lake,   &
                                           np_processes_in_river,  &
                                           np_processes_in_lake,   &
                                           oc_processes_in_river,  &
                                           oc_processes_in_lake,   &
                                           add_diffuse_source_to_local_river,  &
                                           add_point_sources_to_main_river,    &
                                           calculate_river_wetland
    USE SURFACEWATER_PROCESSES, ONLY : add_precipitation_to_river, &
                                       calculate_river_evaporation, &
                                       calculate_actual_lake_evaporation, &
                                       calculate_water_temperature,     &
                                       set_water_temperature, &
                                       calculate_river_characteristics, &
                                       translation_in_river,            &
                                       point_abstraction_from_main_river_inflow, &
                                       point_abstraction_from_main_river, &
                                       point_abstraction_from_outlet_lake, &
                                       calculate_ilake_outflow, &
                                       calculate_outflow_from_lake,     &
                                       check_outflow_from_lake,         &
                                       remove_outflow_from_lake,        &
                                       calculate_flow_within_lake,      &
                                       calculate_olake_waterstage,      &
                                       calculate_regamp_adjusted_waterstage, &
                                       calculate_branched_flow, &
                                       calculate_branched_flow_new, &
                                       calculate_lake_volume, &
                                       calculate_lake_epidepth, &
                                       T2_processes_in_river, &
                                       T2_processes_in_lake, &
                                       ice_processes_in_river, &
                                       ice_processes_in_lake, &
                                       add_T2_concentration_in_precipitation_on_water,  &
                                       get_rivertempvol,  &
                                       add_precipitation_to_floodplain, &
                                       calculate_floodplain_evaporation, &
                                       calculate_waterbody_floodplain_interflow, &
                                       calculate_floodplain_waterlevel
    USE NPC_SOIL_PROCESSES, ONLY : set_class_precipitation_concentration_and_load,  &
                                   croprotation_soilpoolaverage
    USE TRACER_PROCESSES, ONLY : add_tracer_point_source_to_river,  &
                                 add_tracer_point_source_to_lake,   &
                                 tracer_processes_in_river,  &
                                 tracer_processes_in_lake
    USE IRRIGATION_MODULE, ONLY : initiate_timestep_irrigation,  &
                                  calculate_irrigation
    USE REGIONAL_GROUNDWATER_MODULE, ONLY : calculate_regional_groundwater_flow,  &
                                            calculate_current_flow_from_aquifer, &
                                            add_regional_groundwater_flow_to_olake, &
                                            calculate_river_groundwaterflow_removal, &
                                            calculate_river_floodplain_groundwaterflow_removal, &
                                            calculate_aquifer,  &
                                            add_aquifer_flow_to_river, &
                                            calculate_delayed_water, &
                                            calculate_aquifer_waterlevel
    USE SOILMODEL_DEFAULT, ONLY : soilmodel_0
    USE GLACIER_SOILMODEL, ONLY : soilmodel_3, calculate_glacier_massbalance
    USE FLOODPLAIN_SOILMODEL, ONLY : soilmodel_4
    USE ATMOSPHERIC_PROCESSES, ONLY : calculate_class_atmospheric_forcing,  &
                                calculate_subbasin_temperature, &
                                calculate_subbasin_precipitation, &
                                calculate_rain_snow_from_precipitation, & 
                                calculate_extraterrestrial_radiation, &
                                set_precipitation_concentration,  &
                                calculate_daylength
    USE SOIL_PROCESSES, ONLY : calculate_potential_evaporation
    USE GENERAL_WATER_CONCENTRATION, ONLY : remove_water,         &
                                            error_remove_water,   &
                                            add_water
    USE HYPE_INDATA, ONLY : get_current_xoms, &
                            num_xoms, &
                            xom, xos
    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states

    !Variable declarations
    INTEGER i,j,k    !loop-variables: subbasin, class, div.
    INTEGER itype    !loop-variable: ilake=1,olake=2
    INTEGER isl      !loop-variable: soil layer
    
    REAL a,aadj                 !area fraction of class, floodplain-adjusted land fraction of class
    REAL asum(maxsoillayers)    !sum of land area parts
    REAL divasum(maxsoillayers) !inverse of area sums
    !REAL divasum5, divasum6     !subbasin area - upstream area relation
    REAL incorr, oncorr    !correction of inorganic and organic nitrogen level
    REAL phoscorr     !correction of phosphorus level
    
    !Help variables for class summation to subbasin output
    REAL gwati,eacti,epoti,runoffi,snowi
    REAL runoff1i,runoff2i,runoff3i 
    REAL crunoffi(numsubstances)         !(mg/L)
    REAL snowdepthi
    REAL soilltmpi(maxsoillayers)
    REAL icpevapi  
    REAL csoili1(numsubstances),csoili2(numsubstances),csoili3(numsubstances) !(mg/L*mm=kg/km2 soillayer area)
    REAL soilwateri(maxsoillayers),standsoili,soili1,soili2,soili3
    REAL corrpreci,corrtempi, qcinfli, interceptionloss(nsub)
    REAL runofftsi,crunofftsi(numsubstances)  !Runoff to stream (mm,mg/L)
    REAL soildefi
    REAL snowcov,melti
    REAL pcorricep(nsub)
    REAL soillgrossload(3,numsubstances), soillnetload(3,numsubstances)   !1=sl 1-2, 2=sl 3, 3=sl3+tile
    
    !Help variables for summation of upstream output variables
    REAL upstream_prec(nsub),upstream_evap(nsub),upstream_epot(nsub)
    REAL upstream_snow(nsub),upstream_soil(nsub),upstream_soildef(nsub)
    REAL upstream_smfp(nsub),upstream_temp(nsub)
    REAL upstream_snowfall(nsub),upstream_rainfall(nsub),upstream_melt(nsub)
    
    !Current forcing variables for estimating radiation for snow and evaporation calculations
    REAL radexti(nsub)           !Extraterrestrial solar radiation [MJ/m2/day] for current time step
    REAL tmin,tmax,rhmin,netrad,actvap,satvap,wind
    REAL sffrac    !sffrac = snow fall fraction of precipitation [-]
    REAL swrad     !shortwave radiation downward [MJ/m2/day]
    REAL daylength  !day length (hours)
    
    !Variables for class values
    REAL classarea            !class area (km2)
    REAL classheight          !class height (möh)
    REAL lakeareakm2
    REAL gwat,epot,evap,evapsnow,epotsnow   !calculated class values
    REAL surfaceflow,icpevap
    REAL nitrif, denitrif(maxsoillayers)
    REAL cropuptake          !calculate class values, uptake of IN
    REAL runoffd,crunoffd (numsubstances)         !calculated class values
    REAL crunoff1(numsubstances),crunoff2(numsubstances), crunoff3(numsubstances)   !<calculated class values (mg/L)
    REAL csrunoff(numsubstances)
    REAL cevap(numsubstances)   !calculated class values (mg/L)
    REAL cprec(numsubstances),cprecj(numsubstances)   !concentration of precipitation, subbasin and class
    REAL evapl,cevapl(numsubstances)  !evaporation lake
    REAL concout(numsubstances)   !variable for concentration of lake outflow
    REAL temp                 !class elevation corrected temperature
    REAL precorg(nsub) !precipitation original from Pobs.txt (mm)
    REAL prec          !class elevation corrected precipitation (mm,m3)
    REAL smdef    !soil moisture deficit
    REAL frostdepth                    !soil frost variables 
    REAL glac_part        !glacier fraction of glacier class (-)
    REAL snowfall, rainfall  !Precipiation as snow, rain (mm)
    REAL snowplain_area_fraction, old_snow_part
    REAL atmdepload1(numsubstances),atmdepload2(numsubstances)  !wet and drydep
    REAL irrload(numsubstances)
    REAL grwloadlake(numsubstances),grwloadmr(numsubstances)
    REAL rgrwload(numsubstances)
    REAL ruralaload(numsubstances)
    REAL cultivload(2,numsubstances)
    REAL soil4runoff(2)       !runoff from floodplain :1=from main river fp,2=from olake fp
    REAL infiltrationflows(7) !fraction from snow melt,infiltration,surface runoff,macroporeflow to soil layer 1,2,and 3,fraction from glacier melt
    REAL glacierflows(2)      !flow snow to glacier, precipitation on glacier
    REAL floodplainflows(3)   !snow incorporation to floodplain, precipitation on floodplain,infiltration [m3]
    REAL evapflows(4)         !flow evapotranspiration (land classes), 1=sl1, 2=sl2,3=snow,4=glacier
    REAL cevapsnow(numsubstances) ! conc of snow evaporation
    REAL runofflows(7)        !different runoff flows:1-3=soil runoff sl 1-3,4-6=tile runoff sl 1-3,7=saturated surface runoff
    REAL crunofflows(numsubstances,6) !concentration of different runoff flows:1-3=soil runoff sl 1-3,4-6=tile runoff sl 1-3
    REAL verticalflows(6)     !vertical flows:1-2=percolation,3-4=upwelling due to rural,5-6=upwelling due to reg. grw flows
    REAL cverticalflows(2,numsubstances) !<concentration of vertical flows:1-2=percolation
    REAL horizontalflows(3)   !horizontal flows:1-3=recieved rural load flow
    REAL horizontalflows2(3,nsub)   !horizontal flows:1-3=division of regional groundwater flows to grwdown
    REAL cruralflow(numsubstances)   !concentration of rural load flow
    REAL nutrientloadflows(4) !1=rural load flow to local stream,2=point source to main river,3=abstraction from main river,4=abstraction from outlet lake
    REAL irrigationflows(7)   !withdrawal from 1=deep groundwater(unlimited),2=internal lake,3=outlet lake,4=main river volume and inflow,5=local evaporation losses (local sources),6=unlimited external source (alternative),7=akvifer(modelled)
    REAL regionalirrflows(4,nsub) !withdrawal from 1=outlet lake,2=main river volume and inflow,3=evaporation losses olake as regional source,4=evaporation losses main river as regional source
    REAL regionalirrevap(nsub)    !local evaporation losses (regional source)
    REAL glacperiod,computedMassBalance,computedMBArea  !mass balance of glacier
    
    !Variables for routing
    REAL accinflow(nsub)                  !accumulated upstream inflow to subbasin (m3/s)
    REAL acccinflow(numsubstances,nsub)   !concentration of upstream inflow (mg/L)
    REAL mainflow,branchflow              !divided discharge for output (m3/s)
    REAL qinmm,qin,cin(numsubstances)
    REAL inflowpart !part of runoff to ilake
    REAL qunitfactor,outflowm3s,outflowmm
    REAL transq,transc(numsubstances),dampq,dampc(numsubstances)
    REAL lakeoutflow(2)                   !Outflow of lakes (m3/s)
    REAL,ALLOCATABLE :: clakeoutflow(:,:)    !Concentration of outflow of lakes (mg/L)
    REAL outflowsim,wstlakesim              !simulated variables for lakeoutflow/wstlake before update is done
    REAL olakewstold                         !wst olake last time step (mm)
    REAL oldolakewst !olake water stage before updating (m)
    REAL wcomaver                         !average wst olake over time step (mm)
    REAL lakewst(2)   !lake water stage for ilake and olake (mm)
    REAL wstlake,wstlakeadj      !temporary variables for olake water stage output
    REAL w0ref        !olake reference level for water stage
    REAL ilakevol(nsub),olakevol(nsub)  !volume of ilake and olake (m3)
    REAL basinlakevol(nsub),ilakevolmiss(nsub),olakevolmiss(nsub)  !volume of whole lake,ilake and olake (Mm3 or missing value)
    REAL grwout(4,nsub),grwout2         !Outflow from soil to reg. groundwater for print out (all, and resp soillayer)(main river)
    REAL grwtomr,grwtool                !Flow of regional groundwater to main river and olake
    REAL outofsystem                    !loss of water from the model domain via groundwater flow, m3/ts
    REAL riverQbank,riverarea(nrivertypes),riverdepth,riverareakm2
    REAL rivervolume(nrivertypes)   !help for output (m3)
    REAL evapr,cevapr(numsubstances)
    REAL lakearea(nlaketypes)             !lake area (m2) for lake
    REAL lakeareatemp2         !lake area (m2) for outlet lake or whole lake for last lakebasin
    REAL flow1000m3ts         !olake outflow in 1000 m3/ts
    REAL maxProd,minFlow  !current maximum production/minimum flow
    
    !Variables for irrigation 
    REAL pwneedj,pwneedi               !irrigation water demand for current class, and all irrigated classes (m3)
    REAL gwremi                       !groudwater removed for irrigation (mm)
    REAL ldremi,lrremi,rsremi         !abstraction surface water for irrigation (m3)
    REAL irrevap(nsub)                !accumulation variable for irrigation field and network losses (m3)
    REAL irrappl(2)!,irrappli          !irrappl = applied irrigation (mm), for summation basin output, irrapli is the sum of these
    REAL irrsinkmass(numsubstances)   !irrsinkmass = mass of substances in the irrigation sink(kg for mg/L concentrations)
    
    !Variables for the lake and river ice model
    REAL lakesurftemp(2), riversurftemp(2)
    REAL ctemp_T2
    REAL meanrivertemp,totrivervol  !mean temperature and total river volume over all water course elements
    REAL freezeuparea !freezeuparea for river and lakes
    REAL epidepth     !lake epilimnion depth
    REAL outflow1,outflow2  !flow from lake outlets
    INTEGER lakefreezeupday(2), lakebreakupday(2), riverfreezeupday(2), riverbreakupday(2)
    
    REAL delayedwater(naquifers)    !water in delay to reach aquifer (m3) 
    REAL aquiferoutflow(naquifers)
    REAL aqremflow(nsub)            !water removed from soillayers for aquifer
    REAL aqoutside(naquifers)       !water from aquifer to outside model system
    REAL aqirrloss(naquifers)       !water for irrigation from aquifer
    
    !Variables for error handling
    INTEGER status
    CHARACTER(LEN=80) :: errstring(1)  !error message for location of remove_water call
    PARAMETER (errstring = (/'outflow from river damping box'/))
    
    !Variables for comparing forest and open snow cover variables with the FSUS data (Former Soviet Union Snow coarse data)
    REAL snowvari(4,2)    !first index:1=depth,2=cover,3=density,4=swe; second index:1=open,2=forest
    REAL area_forest, area_open, area_water

    !Variables for floodplain model
    REAL fpfrac
    REAL qmrflood,qolflood
    REAL wlm3ol,cwlm3ol(numsubstances)
    REAL riverfparea,lakefparea
    REAL interflowpar(5)  !hold current parameter-values
    REAL floodplain_area_fraction
    REAL ffpwl,ffparea
    LOGICAL floodplainclass
    REAL flooddepth,flooddegree

    !Start of subroutine calculations
    !> \b Algorithm \n
    IF(conductwb) wbflows=0.  !zeroing the timestep for all flows
    
    !>Initial calculation of regional groundwater
    IF(numsubstances>0) Lgrwclass=0.
    grwout = 0.
    IF(modeloption(p_deepgroundwater)==1 .OR. &
       modeloption(p_deepgroundwater)==2)THEN
      CALL calculate_regional_groundwater_flow(soilstate,grwout,Lgrwclass)    !regional groundwater leaving the soil
      !outvar(:,240) = grwout(1,:)   !will be replaced later if aquifermodel
      IF(outvarindex(240)>0) outvar(:,outvarindex(240)) = grwout(1,:)   !will be replaced later if aquifermodel
    ENDIF
    IF(conductwb)THEN
      wbflows(w_rgrwof1,:) = grwout(2,:)
      wbflows(w_rgrwof2,:) = grwout(3,:)
      wbflows(w_rgrwof3,:) = grwout(4,:)
    ENDIF
    IF(modeloption(p_deepgroundwater)==2)THEN     !prepare for adding regional groundwater to model
      CALL calculate_current_flow_from_aquifer(naquifers,numsubstances,aquiferstate,aqoutside)    
      IF(conductwb) wbflows(w_rgrwtoos,1:naquifers) = aqoutside
    ENDIF
        
    !>Get current observations for HYPE specific indata
    IF(conductxoms) CALL get_current_xoms(currentdate,nsub)
    
    !>Calculate corrected atmospheric forcing for all subbasins
    !outvar(:,4) = tempi       !original input data
    !outvar(:,o_prec) = preci  !original input data
    IF(outvarindex(o_tobs)>0) outvar(:,outvarindex(o_tobs)) = tempi  !original input data
    IF(outvarindex(o_prec)>0) outvar(:,outvarindex(o_prec)) = preci  !original input data
    precorg = preci  !original input data
    CALL calculate_subbasin_precipitation(nsub,tempi,preci,pcorricep)
    CALL calculate_subbasin_temperature(nsub,month,tempi)
    !CALL calculate_subbasin_precipitation(nsub,outvar(:,4),preci,pcorricep)
    !outvar(:,206) = pcorricep
    IF(ALLOCATED(tmini))THEN
      CALL calculate_subbasin_temperature(nsub,month,tmini)
      CALL calculate_subbasin_temperature(nsub,month,tmaxi)
    ENDIF
    
    !Calculate extra terrestrial radiation for all subbasins, if needed
    IF(calcSWRAD)THEN
      CALL calculate_extraterrestrial_radiation(nsub,dayno,radexti)
    ENDIF
    
    !Initiations for main subbasin calculation loop
    accinflow = 0.            
    acccinflow = 0.
    irrevap = 0. 
    epidepth = 0. !if T2 is not simulated
    aqirrloss = 0.
    IF(conductwb) wbirrflows = 0.
    IF(conductwb) wbfpflows = 0.
    regionalirrflows=0.
    horizontalflows2=0.
    CALL initiate_timestep_irrigation(numsubstances,miscstate)
    IF(.NOT.ALLOCATED(clakeoutflow)) ALLOCATE(clakeoutflow(numsubstances,2))
    upstream_prec = 0.; upstream_evap = 0.; upstream_snow = 0.
    upstream_soil = 0.; upstream_temp = 0.
    upstream_soildef = 0.; upstream_smfp = 0.; upstream_melt = 0.
    upstream_snowfall = 0.; upstream_rainfall = 0.; upstream_epot = 0.

    !>Main subbasin-loop, subbasins calculated in flow order
    subbasinloop:  &
    DO i = 1,nsub
       
      !Initiate variables for calculation of subbasin averages
      gwati=0.; eacti=0.; epoti=0.; melti=0.
      snowi=0.
      soilwateri=0.
      DO isl = 1,3
        CALL calculate_class_outvar_initialize(o_pfN(isl),i)
        CALL calculate_class_outvar_initialize(o_phN(isl),i)
        CALL calculate_class_outvar_initialize(o_pfP(isl),i)
        CALL calculate_class_outvar_initialize(o_phP(isl),i)
        CALL calculate_class_outvar_initialize(o_ppP(isl),i)
        CALL calculate_class_outvar_initialize(o_pfC(isl),i)
        CALL calculate_class_outvar_initialize(o_phC(isl),i)
        CALL calculate_class_outvar_initialize(o_ppT1(isl),i)
        CALL calculate_class_outvar_initialize(o_sltmp(isl),i)
        CALL calculate_class_outvar_initialize(o_csoillayerIN(isl),i)
      ENDDO
      CALL calculate_class_outvar_initialize(o_csoilIN,i)
      CALL calculate_class_outvar_initialize(o_csoilOC,i)
      CALL calculate_class_outvar_initialize(o_csoilT1,i)
      CALL calculate_class_outvar_initialize(o_csoilT2,i)
      CALL calculate_class_outvar_initialize(o_snowmelt,i)
      CALL calculate_class_outvar_initialize(o_T1sf,i)
      CALL calculate_class_outvar_initialize(o_soiltmp,i)
      CALL calculate_class_outvar_initialize(o_grwlevel,i)
      CALL calculate_class_outvar_initialize(o_epot,i)
      CALL calculate_class_outvar_initialize(o_cevapT1,i)
      CALL calculate_class_outvar_initialize(o_landevap,i)
      CALL calculate_class_outvar_initialize(o_crun,i)
      CALL calculate_class_outvar_initialize(o_ro1,i)
      CALL calculate_class_outvar_initialize(o_ro2,i)
      CALL calculate_class_outvar_initialize(o_ro3,i)
      CALL calculate_class_outvar_initialize(o_rod,i)
      CALL calculate_class_outvar_initialize(o_ros,i)
      CALL calculate_class_outvar_initialize(o_ctmp,i)
      CALL calculate_class_outvar_initialize(o_rainfall,i)
      CALL calculate_class_outvar_initialize(o_snowfall,i)
      CALL calculate_class_outvar_initialize(o_crosT1,i)
      CALL calculate_class_outvar_initialize(o_crodT1,i)
      CALL calculate_class_outvar_initialize(o_cro1T1,i)
      CALL calculate_class_outvar_initialize(o_cro2T1,i)
      CALL calculate_class_outvar_initialize(o_cro3T1,i)
      CALL calculate_class_outvar_initialize(o_crunT1,i)
      CALL calculate_class_outvar_initialize(o_crunT2,i)
      CALL calculate_class_outvar_initialize(o_crunIN,i)
      CALL calculate_class_outvar_initialize(o_crunON,i)
      CALL calculate_class_outvar_initialize(o_crunTN,i)
      CALL calculate_class_outvar_initialize(o_crunSP,i)
      CALL calculate_class_outvar_initialize(o_crunPP,i)
      CALL calculate_class_outvar_initialize(o_crunTP,i)
      CALL calculate_class_outvar_initialize(o_crunOC,i)
      CALL calculate_class_outvar_initialize(o_smffc,i)
      CALL calculate_class_outvar_initialize(o_smfdep,i)
      CALL calculate_class_outvar_initialize(o_smrzfdep,i)
      CALL calculate_class_outvar_initialize(o_smfpw,i)
      CALL calculate_class_outvar_initialize(o_smrzfpw,i)
      CALL calculate_class_outvar_initialize(o_soildef,i)
      CALL calculate_class_outvar_initialize(o_applirr,i)
      CALL calculate_class_outvar_initialize(o_soildenitr,i)
      CALL calculate_class_outvar_initialize(o_soildenrz,i)
      CALL calculate_class_outvar_initialize(o_soilden3,i)
      CALL calculate_class_outvar_initialize(o_cropNupt,i)
      CALL calculate_class_outvar_initialize(o_degrfN,i)
      CALL calculate_class_outvar_initialize(o_soilNatm,i)
      CALL calculate_class_outvar_initialize(o_soilPatm,i)
      CALL calculate_class_outvar_initialize(o_evapsnow,i)
      CALL calculate_class_outvar_initialize(o_soilfrost,i)
      CALL calculate_class_outvar_initialize(o_snowcover,i)
      CALL calculate_class_outvar_initialize(o_snowmax,i)
      CALL calculate_class_outvar_initialize(o_clrp,i)
      CALL calculate_class_outvar_initialize(o_clre,i)
      runoffi=0.
      runoff1i=0.; runoff2i=0.; runoff3i=0.
      runofftsi=0.;
      crunofftsi=0.
      snowdepthi=0.; soilltmpi=0.
      soili1=0.; soili2=0.; soili3=0.; standsoili=0.
      csoili1=0.; csoili2=0.; csoili3=0.
      soildefi=0.; icpevapi=0.
      corrpreci=0.; corrtempi=0.
      qcinfli=0.
      pwneedi=0.
      ldremi=0.;lrremi=0.;gwremi=0.;rsremi=0.
      asum=0.
      soil4runoff=0.
      glacierflows=0.;nutrientloadflows=0.;irrigationflows=0.
      snowvari=0.
      area_forest = 0.;area_open = 0.;area_water=0.
      soillgrossload =0.;soillnetload=0.
      
      !Initiate other subbasin variables
      IF(conductload)THEN
        Latmdep=0.
        Lcultiv=0.
        Lirrsoil=0.
        Lrurala=0.
        Lgrwsoil=0.
        Lgrwol=0.
        Lgrwmr=0.
        IF(numsubstances>0) Lstream=0.
      ENDIF
       
      !Short notation for parameters not dependent on class
      IF(basin(i)%wqparregion>0)THEN
        incorr = (1. + wqregpar(m_incorr,basin(i)%wqparregion))     !Correction of inorgnic nitrogen
        oncorr = (1. + wqregpar(m_oncorr,basin(i)%wqparregion))     !Correction of organic nitrogen
        phoscorr = 1. + wqregpar(m_phoscorr,basin(i)%wqparregion)   !Correction of phosphorus level
      ELSE
        incorr    = 1.
        oncorr    = 1.
        phoscorr  = 1.
      ENDIF
       
      !Calculation of mean air temperature
      miscstate%temp5(i) = miscstate%temp5(i) + (tempi(i) - miscstate%temp5(i)) / avertemp(1)
      miscstate%temp10(i) = miscstate%temp10(i) + (tempi(i) - miscstate%temp10(i)) / avertemp(2)
      miscstate%temp20(i) = miscstate%temp20(i) + (tempi(i) - miscstate%temp20(i)) / avertemp(3)
      miscstate%temp30(i) = miscstate%temp30(i) + (tempi(i) - miscstate%temp30(i)) / avertemp(4)

      !Set subbasin precipitation concentrations
      CALL set_precipitation_concentration(i,numsubstances,cprec)
       
      !Calculate pseudo dayno and daylength for temperature dependent growing season
      pseudo_dayno = get_pseudo_dayno(dayno,basin(i)%latitude)
      daylength = 0.
      IF(modeloption(p_growthstart)==1)  &
        CALL calculate_daylength(dayno,basin(i)%latitude,daylength)

      !>Main class-loop for calculation of soil water and substances (land classes)
      DO j=1,nclass
        a=classbasin(i,j)%part
        classarea = a * basin(i)%area * 1.0E-6    !km2
        classheight = basin(i)%elev+classbasin(i,j)%deltah  !möh
        IF(a>0)THEN
          irrappl = 0.    !not set for floodplain or glaciers
          pwneedj = 0.
          cultivload = 0.
          irrload = 0.
          ruralaload = 0.
          atmdepload1 = 0.
          atmdepload2 = 0.
          rgrwload = 0.
          !runoffd = 0.   !not set for floodplain (not used)
          crunoffd = 0.   !not set for floodplain
          horizontalflows=0.  !not set for floodplain
          old_snow_part = 1.
          snowplain_area_fraction = 1.
          floodplain_area_fraction = 1.
          floodplainclass = .FALSE.
          !> \li Calculate class forcing data
          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
          CALL set_class_precipitation_concentration_and_load(numsubstances, &
                   classarea,precorg(i),temp,prec,cprec,cprecj,atmdepload1)

          !> \li Calculate soil processes            
          IF(classmodel(j)==0)THEN
            CALL soilmodel_0(i,j,classdata(j)%soil,classdata(j)%luse,basin(i)%subid,classarea,classheight,prec,cprecj,temp,  & 
                 daylength,basinrrcscorr(i),phoscorr,basincevpcorr(i),incorr,oncorr,pseudo_dayno,  &
                 frozenstate,soilstate,miscstate,runoffd,surfaceflow,crunoffd,csrunoff,    &
                 cropuptake,nitrif,denitrif,epot,gwat,frostdepth,irrappl,smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
                 pwneedj,snowfall,rainfall,cultivload,irrload,ruralaload,rgrwload,  &
                 atmdepload2,sffrac,swrad,radexti(i),netrad,actvap,satvap,wind, &
                 infiltrationflows,evapflows,runofflows,crunofflows,verticalflows,  &
                 cverticalflows,horizontalflows,horizontalflows2,evapsnow,cruralflow) 
          ELSEIF(j==slc_lriver)THEN
            CYCLE
          ELSEIF(j==slc_mriver)THEN
            IF(conductflood)THEN
              IF(floodindex(i)>0)THEN
                IF(flooding(floodindex(i))%fpfmr>0.)THEN
                  CALL soilmodel_4(i,j,classdata(j)%soil,classdata(j)%luse,basin(i)%subid,classarea,prec,cprecj,temp,  & 
                   daylength,basinrrcscorr(i),phoscorr,basincevpcorr(i),incorr,oncorr,pseudo_dayno, &
                   frozenstate,soilstate,miscstate,surfaceflow,csrunoff,crunoffd,    &
                   cropuptake,nitrif,denitrif,epot,gwat,frostdepth,snowcov,smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
                   pwneedj,irrappl,irrload,old_snow_part,snowfall,rainfall,cultivload,ruralaload,atmdepload1,atmdepload2,sffrac,swrad,radexti(i),netrad,actvap,satvap,wind, &
                   infiltrationflows,floodplainflows,evapflows,runofflows,verticalflows,cverticalflows,horizontalflows,horizontalflows2,evapsnow,cevapsnow,soil4runoff,cruralflow) 
                  floodplain_area_fraction = flooding(floodindex(i))%fpfmr
                  floodplainclass = .TRUE.
                ELSE
                  CYCLE
                ENDIF
              ELSE
                CYCLE
              ENDIF
            ELSE
              CYCLE
            ENDIF
          ELSEIF(classmodel(j)==glacier_model)THEN
            !only one glacier class is allowed
            CALL soilmodel_3(i,j,classdata(j)%soil,classdata(j)%luse,basin(i)%subid,classarea,prec,cprecj,temp, & 
                  daylength,basincevpcorr(i),basinrrcscorr(i),phoscorr,pseudo_dayno, &
                  frozenstate,soilstate,miscstate,runoffd,surfaceflow,crunoffd,csrunoff,  &
                  cropuptake,nitrif,denitrif,epot,gwat,frostdepth,smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
                  glac_part,old_snow_part,snowfall,rainfall,cultivload,ruralaload, &
                  rgrwload,atmdepload2,sffrac,swrad,radexti(i),netrad,actvap,satvap,wind, &
                  infiltrationflows,glacierflows,evapflows,runofflows,crunofflows,  &
                  verticalflows,cverticalflows,horizontalflows,horizontalflows2,evapsnow,cruralflow)
            snowplain_area_fraction = 1. - glac_part
            IF(outvarindex(117)>0) outvar(i,outvarindex(117)) = frozenstate%glacvol(i)*1.E-9    !glacier ice volume (km3)
            IF(outvarindex(118)>0) outvar(i,outvarindex(118)) = glac_part * classarea           !glacier area (km2)
            CALL set_outvar_xobs(256,i)
            IF(xobsindex(256,i).GT.0)THEN
              glacperiod = xobsi(xobsindex(256,i)) !Glacier mass balance period
            ELSE
              glacperiod = missing_value
            ENDIF
            IF(glacperiod /= missing_value)THEN   !If there is a Mass balance period to evaluate:
              CALL set_outvar_xobs(253,i) !recorded mass balance
              CALL set_outvar_xobs(255,i) !recorded mass balance area
              CALL calculate_glacier_massbalance(i,glacperiod,computedMassBalance,computedMBArea)   !calculate simulated glacier mass balance and area
              IF(outvarindex(252)>0) outvar(i,outvarindex(252)) = computedMassBalance
              IF(outvarindex(254)>0) outvar(i,outvarindex(254)) = computedMBArea
            ENDIF
          ELSEIF(j==slc_olake)THEN
            IF(conductflood)THEN
              IF(floodindex(i)>0)THEN
                IF(flooding(floodindex(i))%fpfol>0.)THEN
                  CALL soilmodel_4(i,j,classdata(j)%soil,classdata(j)%luse,basin(i)%subid,classarea,prec,cprecj,temp,  & 
                   daylength,basinrrcscorr(i),phoscorr,basincevpcorr(i),incorr,oncorr,pseudo_dayno,  &
                   frozenstate,soilstate,miscstate,surfaceflow,csrunoff,crunoffd,    &
                   cropuptake,nitrif,denitrif,epot,gwat,frostdepth,snowcov,smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
                   pwneedj,irrappl,irrload,old_snow_part,snowfall,rainfall,cultivload,ruralaload,atmdepload1,atmdepload2,sffrac,swrad,radexti(i),netrad,actvap,satvap,wind, &
                   infiltrationflows,floodplainflows,evapflows,runofflows,verticalflows,cverticalflows,horizontalflows,horizontalflows2,evapsnow,cevapsnow,soil4runoff,cruralflow) 
                  floodplain_area_fraction = flooding(floodindex(i))%fpfol
                  floodplainclass = .TRUE.
                ELSE
                  CYCLE
                ENDIF
              ELSE
                CYCLE
              ENDIF
            ELSE
              CYCLE
            ENDIF
          ELSEIF(j==slc_ilake)THEN
            CYCLE
          ELSE
            CALL soilmodel_0(i,j,classdata(j)%soil,classdata(j)%luse,basin(i)%subid,classarea,classheight,prec,cprecj,temp,  & 
                 daylength,basinrrcscorr(i),phoscorr,basincevpcorr(i),incorr,oncorr,pseudo_dayno,  &
                 frozenstate,soilstate,miscstate,runoffd,surfaceflow,crunoffd,csrunoff,    &
                 cropuptake,nitrif,denitrif,epot,gwat,frostdepth,irrappl,smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
                 pwneedj,snowfall,rainfall,cultivload,irrload,ruralaload,rgrwload,  &
                 atmdepload2,sffrac,swrad,radexti(i),netrad,actvap,satvap,wind, &
                 infiltrationflows,evapflows,runofflows,crunofflows,verticalflows,  &
                 cverticalflows,horizontalflows,horizontalflows2,evapsnow,cruralflow) 
          ENDIF
             
          !Set load variables for printout
          IF(conductload)THEN 
            Latmdep(j,1,:) = atmdepload1      !including flooded floodplain atmdep (=0)
            Latmdep(j,2,:) = atmdepload2      !-"-                                 (>0)
            IF(.NOT.floodplainclass)THEN
              Lcultiv(j,:,:) = cultivload     !not including non-flooded floodplain
              Lirrsoil(j,:)  = irrload
              Lrurala(j,:)   = ruralaload
              Lgrwsoil(j,:)  = rgrwload
            ENDIF
          ENDIF
             
          !Set runoff concentration to T1 leakage parameters
          IF(T1leakage)THEN
            crunoff1(i_t1) = landpar(m_t1leaklu,classdata(j)%luse)*soilpar(m_t1leakst,classdata(j)%soil)
            crunoff2(i_t1) = crunoff1(i_t1)
            crunoff3(i_t1) = crunoff1(i_t1) 
            crunoffd(i_t1) = crunoff1(i_t1)
            csrunoff(i_t1) = crunoff1(i_t1)
          ENDIF

          !> \li Accumulate variables for mean over subbasin and layers
          aadj = a*floodplain_area_fraction           !landarea fraction of classarea fraction for mriver/olake class
          asum(1) = asum(1) + aadj                    !landarea fraction (dry floodplains included)
          IF(.NOT.floodplainclass)THEN  !These will be added for floodplain in routing loop
            CALL calculate_class_outvar_add(o_ctmp,i,a,temp)
            CALL calculate_class_outvar_add(o_rainfall,i,a,rainfall)
            CALL calculate_class_outvar_add(o_snowfall,i,a,snowfall)
            corrpreci = corrpreci + prec*a
            corrtempi = corrtempi + temp*a
          ENDIF
          pwneedi = pwneedi + pwneedj
          IF(conductwb.OR.outvarindex(o_snow)>0.OR.outvarindex(o_snowdens)>0) snowi = snowi + frozenstate%snow(j,i)*aadj*snowplain_area_fraction  !zero snow on glacier area
          IF(outvarindex(o_snowdepth)>0.OR.outvarindex(o_snowdens)>0) snowdepthi = snowdepthi + frozenstate%snowdepth(j,i)*aadj*snowplain_area_fraction    !-"-
          melti = melti+infiltrationflows(1)*SUM(infiltrationflows(2:6))*aadj
          CALL calculate_class_outvar_add(o_snowmelt,i,aadj,infiltrationflows(1)*SUM(infiltrationflows(2:6)))
          CALL calculate_class_outvar_add(o_snowmax,i,aadj*snowplain_area_fraction,frozenstate%snowmax(j,i))
          soilwateri = soilwateri+soilstate%water(:,j,i)*aadj
          soili1 = soili1+soilstate%water(1,j,i)*aadj
          csoili1(:) = csoili1(:) + soilstate%conc(:,1,j,i)*soilstate%water(1,j,i)*aadj     !accumulation of amount
          runoff1i = runoff1i + runofflows(1)*aadj
          runoffd = SUM(runofflows(4:6))
          CALL calculate_class_outvar_add(o_ro1,i,aadj,runofflows(1))
          CALL calculate_class_outvar_add(o_ro2,i,aadj,runofflows(2))
          CALL calculate_class_outvar_add(o_ro3,i,aadj,runofflows(3))
          CALL calculate_class_outvar_add(o_rod,i,aadj,SUM(runofflows(4:6)))
          CALL calculate_class_outvar_add(o_ros,i,aadj,surfaceflow)
          IF(i_t1>0) CALL calculate_class_outvar_add(o_crosT1,i,aadj,csrunoff(i_t1)*surfaceflow)
          IF(i_t1>0) CALL calculate_class_outvar_add(o_crodT1,i,aadj,crunoffd(i_t1)*runoffd)
          IF(i_t1>0) CALL calculate_class_outvar_add(o_cro1T1,i,aadj,crunoff1(i_t1)*runofflows(1))  
          IF(i_t1>0) CALL calculate_class_outvar_add(o_cro2T1,i,aadj,crunoff2(i_t1)*runofflows(2))  
          IF(i_t1>0) CALL calculate_class_outvar_add(o_cro3T1,i,aadj,crunoff3(i_t1)*runofflows(3))  
          IF(.NOT.floodplainclass)THEN    !floodplain runoff goes to flooded water not to local stream
            runofftsi = runofftsi + runoffd*a + surfaceflow*a +runofflows(1)*a+runofflows(2)*a+runofflows(3)*a
            crunofftsi(:) = crunofftsi(:) + crunoffd(:)*runoffd*a + csrunoff(:)*surfaceflow*a + crunoff1(:)*runofflows(1)*a + crunoff2(:)*runofflows(2)*a + crunoff3(:)*runofflows(3)*a    !accumulation of amount
          ENDIF
          CALL calculate_class_outvar_add(o_soildenitr,i,aadj,SUM(denitrif(:)))
          CALL calculate_class_outvar_add(o_soildenrz,i,aadj,denitrif(1)+denitrif(2))
          CALL calculate_class_outvar_add(o_soilden3,i,aadj,denitrif(3))
          CALL calculate_class_outvar_add(o_cropNupt,i,aadj,cropuptake)
          CALL calculate_class_outvar_add(o_degrfN,i,aadj,nitrif)
          IF(i_in>0) CALL calculate_class_outvar_add(o_soilNatm,i,a,(atmdepload1(i_in)+atmdepload2(i_in))/classarea)    !Note, already scaled to fpfrac (kg)
          IF(i_sp>0) CALL calculate_class_outvar_add(o_soilPatm,i,a,(atmdepload1(i_sp)+atmdepload2(i_sp))/classarea)
          eacti = eacti+evap*aadj
          CALL calculate_class_outvar_add(o_epot,i,aadj,epot)
          CALL calculate_class_outvar_add(o_landevap,i,aadj,evap)    !no river,lake but include floodplains
          IF(i_t1>0) CALL calculate_class_outvar_add(o_cevapT1,i,aadj,cevap(i_t1)*evap)  !accumulation of amount
          icpevapi = icpevapi + icpevap*aadj
          CALL calculate_class_outvar_add(o_evapsnow,i,aadj,evapsnow)
          CALL calculate_class_outvar_add(o_grwlevel,i,aadj,gwat)
          CALL calculate_class_outvar_add(o_soilfrost,i,aadj,frostdepth)
          CALL calculate_class_outvar_add(o_soiltmp,i,aadj,(soilstate%temp(1,j,i)*soilthick(1,j)+soilstate%temp(2,j,i)*soilthick(2,j)+soilstate%temp(3,j,i)*soilthick(3,j))/soildepth(3,j))
          DO isl=1,3
            CALL calculate_class_outvar_add(o_sltmp(isl),i,aadj,soilstate%temp(isl,j,i))
          ENDDO
          CALL calculate_class_outvar_add(o_applirr,i,aadj,irrappl(1)+irrappl(2))
          CALL calculate_class_outvar_add(o_soildef,i,aadj,smdef)
          soildefi = soildefi + smdef*aadj
          IF(soilstate%water(1,j,i)>pwmm(1,j))THEN
            standsoili = standsoili + (soilstate%water(1,j,i)-pwmm(1,j))*aadj
            CALL calculate_class_outvar_add(o_smffc,i,aadj,(pwmm(1,j)+soilstate%water(2,j,i)-wpmm(1,j)-wpmm(2,j))/(fcmm(1,j)+fcmm(2,j)))
            CALL calculate_class_outvar_add(o_smfdep,i,aadj,(pwmm(1,j)+soilstate%water(2,j,i)+soilstate%water(3,j,i))/(1000.*soildepth(3,j)))
            CALL calculate_class_outvar_add(o_smfpw,i,aadj,(pwmm(1,j)+soilstate%water(2,j,i)+soilstate%water(3,j,i))/(pwmm(1,j)+pwmm(2,j)+pwmm(3,j)))
            CALL calculate_class_outvar_add(o_smrzfdep,i,aadj,(pwmm(1,j)+soilstate%water(2,j,i))/(1000.*soildepth(2,j)))
            CALL calculate_class_outvar_add(o_smrzfpw,i,aadj,(pwmm(1,j)+soilstate%water(2,j,i))/(pwmm(1,j)+pwmm(2,j)))
          ELSE
            CALL calculate_class_outvar_add(o_smffc,i,aadj,(soilstate%water(1,j,i)+soilstate%water(2,j,i)-wpmm(1,j)-wpmm(2,j))/(fcmm(1,j)+fcmm(2,j)))
            CALL calculate_class_outvar_add(o_smfdep,i,aadj,(soilstate%water(1,j,i)+soilstate%water(2,j,i)+soilstate%water(3,j,i))/(1000.*soildepth(3,j)))
            CALL calculate_class_outvar_add(o_smfpw,i,aadj,(soilstate%water(1,j,i)+soilstate%water(2,j,i)+soilstate%water(3,j,i))/(pwmm(1,j)+pwmm(2,j)+pwmm(3,j)))
            CALL calculate_class_outvar_add(o_smrzfdep,i,aadj,(soilstate%water(1,j,i)+soilstate%water(2,j,i))/(1000.*soildepth(2,j)))
            CALL calculate_class_outvar_add(o_smrzfpw,i,aadj,(soilstate%water(1,j,i)+soilstate%water(2,j,i))/(pwmm(1,j)+pwmm(2,j)))
          ENDIF
          IF(conductN)THEN
            DO isl=1,3
              CALL calculate_class_outvar_add(o_pfN(isl),i,aadj,soilstate%fastN(isl,j,i))
              CALL calculate_class_outvar_add(o_phN(isl),i,aadj,soilstate%humusN(isl,j,i))
              CALL calculate_class_outvar_add(o_csoillayerIN(isl),i,aadj,soilstate%conc(i_in,isl,j,i)*soilstate%water(isl,j,i))
            ENDDO
          ENDIF
          IF(conductP)THEN
            DO isl=1,3
              CALL calculate_class_outvar_add(o_pfP(isl),i,aadj,soilstate%fastP(isl,j,i))
              CALL calculate_class_outvar_add(o_phP(isl),i,aadj,soilstate%humusP(isl,j,i))
              CALL calculate_class_outvar_add(o_ppP(isl),i,aadj,soilstate%partP(isl,j,i))
            ENDDO
          ENDIF
          IF(conductC)THEN
            DO isl=1,3
              CALL calculate_class_outvar_add(o_pfC(isl),i,aadj,soilstate%fastC(isl,j,i))
              CALL calculate_class_outvar_add(o_phC(isl),i,aadj,soilstate%humusC(isl,j,i))
            ENDDO
          ENDIF
          IF(i_t1>0)THEN
            DO isl=1,3
              CALL calculate_class_outvar_add(o_ppT1(isl),i,aadj,soilstate%partT1(isl,j,i))
            ENDDO
            CALL calculate_class_outvar_add(o_T1sf,i,aadj,miscstate%partT1sf(j,i))
          ENDIF

          !Soil load output, mm*mg/L
          IF(numsubstances>0.AND.SUM(outvarindex(285:322))>0)THEN
            soillgrossload(1,:)=soillgrossload(1,:)+a*cruralflow(:)*horizontalflows(1)+a*cruralflow(:)*horizontalflows(2)+  &
              (atmdepload1+atmdepload2+cultivload(1,:)+cultivload(2,:))/basin(i)%area*1.E6
            soillnetload(1,:)=soillnetload(1,:)+a*cverticalflows(2,:)*verticalflows(2)+a*crunofflows(:,1)*runofflows(1)+  &
              a*crunofflows(:,2)*runofflows(2)+a*crunofflows(:,4)*runofflows(4)+a*crunofflows(:,5)*runofflows(5)+   &
              a*csrunoff(:)*surfaceflow
            soillgrossload(2,:)=soillgrossload(2,:)+a*cverticalflows(2,:)*verticalflows(2)+a*cruralflow(:)*horizontalflows(3)
            soillnetload(2,:)=soillnetload(2,:)+a*crunofflows(:,3)*runofflows(3)+a*crunofflows(:,6)*runofflows(6)
            soillgrossload(3,:)=soillgrossload(3,:)+a*cverticalflows(2,:)*verticalflows(2)+a*cruralflow(:)*horizontalflows(3)+a*crunofflows(:,4)*runofflows(4)+a*crunofflows(:,5)*runofflows(5)
            soillnetload(3,:)=soillnetload(3,:)+a*crunofflows(:,3)*runofflows(3)+a*crunofflows(:,4)*runofflows(4)+a*crunofflows(:,5)*runofflows(5)+a*crunofflows(:,6)*runofflows(6)
          ENDIF

          CALL calculate_class_outvar_add(o_crun,i,aadj,runofflows(1)+runofflows(2)+runofflows(3)+surfaceflow+runoffd)
          runoffi=runoffi+runofflows(1)*aadj+runofflows(2)*aadj+runofflows(3)*aadj+surfaceflow*aadj+runoffd*aadj    !TODO: runoffflows nollställs ok, så kan använda denna ekvation oberoende av antal skikt
          crunoffi(:) = crunoff1(:)*runofflows(1) + crunoff2(:)*runofflows(2) + crunoff3(:)*runofflows(3) + csrunoff(:)*surfaceflow + crunoffd(:)*runoffd    !amount substance of runoff
          IF(i_t1>0) CALL calculate_class_outvar_add(o_crunT1,i,aadj,crunoffi(i_t1))
          IF(i_t2>0) CALL calculate_class_outvar_add(o_crunT2,i,aadj,crunoffi(i_t2))
          IF(i_in>0) CALL calculate_class_outvar_add(o_crunIN,i,aadj,crunoffi(i_in))
          IF(i_on>0) CALL calculate_class_outvar_add(o_crunON,i,aadj,crunoffi(i_on))
          IF(i_in>0) CALL calculate_class_outvar_add(o_crunTN,i,aadj,crunoffi(i_in)+crunoffi(i_on))
          IF(i_sp>0) CALL calculate_class_outvar_add(o_crunSP,i,aadj,crunoffi(i_sp))
          IF(i_pp>0) CALL calculate_class_outvar_add(o_crunPP,i,aadj,crunoffi(i_pp))
          IF(i_sp>0) CALL calculate_class_outvar_add(o_crunTP,i,aadj,crunoffi(i_sp)+crunoffi(i_pp))
          IF(i_oc>0) CALL calculate_class_outvar_add(o_crunOC,i,aadj,crunoffi(i_oc))
          IF(i_in>0) CALL calculate_class_outvar_add(o_csoilIN,i,aadj,soilstate%conc(i_in,1,j,i)*soilstate%water(1,j,i) + soilstate%conc(i_in,2,j,i)*soilstate%water(2,j,i) + soilstate%conc(i_in,3,j,i)*soilstate%water(3,j,i))
          IF(i_oc>0) CALL calculate_class_outvar_add(o_csoilOC,i,aadj,soilstate%conc(i_oc,1,j,i)*soilstate%water(1,j,i) + soilstate%conc(i_oc,2,j,i)*soilstate%water(2,j,i) + soilstate%conc(i_oc,3,j,i)*soilstate%water(3,j,i))
          IF(i_t1>0) CALL calculate_class_outvar_add(o_csoilT1,i,aadj,soilstate%conc(i_t1,1,j,i)*soilstate%water(1,j,i) + soilstate%conc(i_t1,2,j,i)*soilstate%water(2,j,i) + soilstate%conc(i_t1,3,j,i)*soilstate%water(3,j,i))
          IF(i_t2>0) CALL calculate_class_outvar_add(o_csoilT2,i,aadj,soilstate%conc(i_t2,1,j,i)*soilstate%water(1,j,i) + soilstate%conc(i_t2,2,j,i)*soilstate%water(2,j,i) + soilstate%conc(i_t2,3,j,i)*soilstate%water(3,j,i))

          IF(soilthick(2,j)>0)THEN
            asum(2) = asum(2) + aadj
            soili2 = soili2+soilstate%water(2,j,i)*aadj
            csoili2(:) = csoili2(:) + soilstate%conc(:,2,j,i)*soilstate%water(2,j,i)*aadj   !accumulation of amount
            runoff2i = runoff2i + runofflows(2)*aadj
            IF(soilthick(3,j)>0)THEN
              asum(3) = asum(3) + aadj
              soili3 = soili3+soilstate%water(3,j,i)*aadj
              csoili3(:) = csoili3(:) + soilstate%conc(:,3,j,i)*soilstate%water(3,j,i)*aadj   !accumulation of amount
              runoff3i = runoff3i + runofflows(3)*aadj
            ENDIF
          ENDIF
          !Loads from soil: soil runoff, surface runoff, drain pipe runoff
          IF(conductload)THEN
            Lstream(j,:) = (runofflows(1)*crunoff1(:)+runofflows(2)*crunoff2(:)+  &
                            runofflows(3)*crunoff3(:)+surfaceflow*csrunoff(:)+    &
                            crunoffd(:)*runoffd)*classarea*floodplain_area_fraction    !Load in runoff (kg/timestep) (incl. surface runoff and tile runoff)
          ENDIF
          !Snow cover FSC (must be rescaled with the fraction of ilake+olake before output)
          IF(.NOT.floodplainclass)THEN
            CALL calculate_class_outvar_add(o_snowcover,i,aadj*snowplain_area_fraction,frozenstate%snowcov(j,i))
          ELSE
            CALL calculate_class_outvar_add(o_snowcover,i,aadj*snowplain_area_fraction,snowcov)
          ENDIF
          !Forest and Open field snow variables defined by vegtype
          IF(SUM(outvarindex(267:274))>0)THEN
            IF(vegtype(j)==1.OR.vegtype(j)==2)THEN      !Open field or forest landuse set
              snowvari(1,vegtype(j)) = snowvari(1,vegtype(j)) + frozenstate%snowdepth(j,i)*aadj*snowplain_area_fraction
              IF(frozenstate%snowdepth(j,i)>0.)THEN
                snowvari(3,vegtype(j)) = snowvari(3,vegtype(j)) + 0.1 * frozenstate%snow(j,i)/frozenstate%snowdepth(j,i) * aadj*snowplain_area_fraction
              ENDIF
              IF(.NOT.floodplainclass)THEN
                snowvari(2,vegtype(j)) = snowvari(2,vegtype(j)) + 10. * frozenstate%snowcov(j,i)*aadj*snowplain_area_fraction
              ELSE
                snowvari(2,vegtype(j)) = snowvari(2,vegtype(j)) + 10. * snowcov*aadj*snowplain_area_fraction
              ENDIF
              snowvari(4,vegtype(j)) = snowvari(4,vegtype(j)) + frozenstate%snow(j,i)*aadj*snowplain_area_fraction
              IF(vegtype(j)==1) area_open   = area_open + aadj*snowplain_area_fraction
              IF(vegtype(j)==2) area_forest = area_forest + aadj*snowplain_area_fraction
            ENDIF
          ENDIF
          IF(conductwb)THEN
            wbflows(w_sfallsnow,i)= wbflows(w_sfallsnow,i) + snowfall * aadj * old_snow_part
            wbflows(w_smeltsl1,i) = wbflows(w_smeltsl1,i) + infiltrationflows(1) * infiltrationflows(2) * aadj
            wbflows(w_smeltmp1,i) = wbflows(w_smeltmp1,i) + infiltrationflows(1) * infiltrationflows(4) * aadj
            wbflows(w_smeltmp2,i) = wbflows(w_smeltmp2,i) + infiltrationflows(1) * infiltrationflows(5) * aadj
            wbflows(w_smeltmp3,i) = wbflows(w_smeltmp3,i) + infiltrationflows(1) * infiltrationflows(6) * aadj
            wbflows(w_infrain,i)  = wbflows(w_infrain,i) + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(2) * aadj
            wbflows(w_rainmp1,i)  = wbflows(w_rainmp1,i) + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(4) * aadj
            wbflows(w_rainmp2,i)  = wbflows(w_rainmp2,i) + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(5) * aadj
            wbflows(w_rainmp3,i)  = wbflows(w_rainmp3,i) + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(6) * aadj
            wbflows(w_evap1,i)    = wbflows(w_evap1,i) + evapflows(1) * aadj
            wbflows(w_evap2,i)    = wbflows(w_evap2,i) + evapflows(2) * aadj
            wbflows(w_evap3,i)    = wbflows(w_evap3,i) + evapflows(3) * aadj * old_snow_part
            IF(.NOT. floodplainclass)THEN
              wbflows(w_smeltsr,i)  = wbflows(w_smeltsr,i)  + infiltrationflows(1) * infiltrationflows(3) * aadj
              wbflows(w_rainsr,i)   = wbflows(w_rainsr,i)  + (1.-infiltrationflows(1)-infiltrationflows(7)) * infiltrationflows(3) * aadj
              wbflows(w_gwrunf1,i)  = wbflows(w_gwrunf1,i) + runofflows(1) * a
              wbflows(w_gwrunf2,i)  = wbflows(w_gwrunf2,i) + runofflows(2) * a
              wbflows(w_gwrunf3,i)  = wbflows(w_gwrunf3,i) + runofflows(3) * a
              wbflows(w_tile1,i)    = wbflows(w_tile1,i) + runofflows(4) * a
              wbflows(w_tile2,i)    = wbflows(w_tile2,i) + runofflows(5) * a
              wbflows(w_tile3,i)    = wbflows(w_tile3,i) + runofflows(6) * a
              wbflows(w_surfrf,i)   = wbflows(w_surfrf,i) + runofflows(7) * a
            ELSE
              IF(floodplainclass.AND.j==slc_mriver)THEN
                wbfpflows(w_grf1mrfp,i)  = runofflows(1) * aadj
                wbfpflows(w_grf2mrfp,i)  = runofflows(2) * aadj
                wbfpflows(w_grf3mrfp,i)  = runofflows(3) * aadj
                wbfpflows(w_trf1mrfp,i)  = runofflows(4) * aadj
                wbfpflows(w_trf2mrfp,i)  = runofflows(5) * aadj
                wbfpflows(w_trf3mrfp,i)  = runofflows(6) * aadj
                wbfpflows(w_smtsrmrfp,i)  = infiltrationflows(1) * infiltrationflows(3) * aadj
                wbfpflows(w_rtsrmrfp,i)   = (1.-infiltrationflows(1)) * infiltrationflows(3) * aadj
                wbfpflows(w_srftmrfp,i)   = runofflows(7) * aadj
              ELSEIF(floodplainclass.AND.j==slc_olake)THEN
                wbfpflows(w_grf1olfp,i)  = runofflows(1) * aadj
                wbfpflows(w_grf2olfp,i)  = runofflows(2) * aadj
                wbfpflows(w_grf3olfp,i)  = runofflows(3) * aadj
                wbfpflows(w_trf1olfp,i)  = runofflows(4) * aadj
                wbfpflows(w_trf2olfp,i)  = runofflows(5) * aadj
                wbfpflows(w_trf3olfp,i)  = runofflows(6) * aadj
                wbfpflows(w_smtsrolfp,i)  = infiltrationflows(1) * infiltrationflows(3) * aadj
                wbfpflows(w_rtsrolfp,i)   = (1.-infiltrationflows(1)) * infiltrationflows(3) * aadj
                wbfpflows(w_srftolfp,i)   = runofflows(7) * aadj
              ENDIF
            ENDIF
            wbflows(w_perc1,i)    = wbflows(w_perc1,i) + verticalflows(1) * aadj
            wbflows(w_perc2,i)    = wbflows(w_perc2,i) + verticalflows(2) * aadj
            wbflows(w_upwell1,i)  = wbflows(w_upwell1,i) + (verticalflows(3)+verticalflows(5)) * aadj
            wbflows(w_upwell2,i)  = wbflows(w_upwell2,i) + (verticalflows(4)+verticalflows(6)) * aadj
            wbflows(w_rural1,i)   = wbflows(w_rural1,i) + horizontalflows(1) * aadj
            wbflows(w_rural2,i)   = wbflows(w_rural2,i) + horizontalflows(2) * aadj
            wbflows(w_rural3,i)   = wbflows(w_rural3,i) + horizontalflows(3) * aadj
            wbirrflows(w_apply1,i)  = wbirrflows(w_apply1,i) + irrappl(1) * aadj
            wbirrflows(w_apply2,i)  = wbirrflows(w_apply2,i) + irrappl(2) * aadj
            IF(classmodel(j)==glacier_model)THEN
              wbstores(w_glacier,i) = frozenstate%glacvol(i)*genpar(m_glacdens)
              wbflows(w_stoice,i)   = glacierflows(1)
              wbflows(w_precglac,i) = glacierflows(2)
              wbflows(w_gmeltsl1,i) = infiltrationflows(7) * infiltrationflows(2) * a
              wbflows(w_gmeltsr,i)  = infiltrationflows(7) * infiltrationflows(3) * a
              wbflows(w_gmeltmp1,i) = infiltrationflows(7) * infiltrationflows(4) * a
              wbflows(w_gmeltmp2,i) = infiltrationflows(7) * infiltrationflows(5) * a
              wbflows(w_gmeltmp3,i) = infiltrationflows(7) * infiltrationflows(6) * a
              wbflows(w_evap4,i)    = wbflows(w_evap4,i) + evapflows(4) * a * (1. - old_snow_part)
            ENDIF
            IF(floodplainclass.AND.j==slc_mriver)THEN
              wbfpflows(w_emrfp,i)    = evapflows(4)    !these are in m3
              wbfpflows(w_stomrfp,i)  = floodplainflows(1)
              wbfpflows(w_pmrfp,i)    = floodplainflows(2)
              wbfpflows(w_infmrfp,i)  = floodplainflows(3)
            ELSEIF(floodplainclass.AND.j==slc_olake)THEN
              wbfpflows(w_eolfp,i)    = evapflows(4)
              wbfpflows(w_stoolfp,i)  = floodplainflows(1)
              wbfpflows(w_polfp,i)    = floodplainflows(2)
              wbfpflows(w_infolfp,i)  = floodplainflows(3)
            ENDIF
          ENDIF
        ENDIF !a>0
      ENDDO   !j
        
      !Calculate average of soil pools for classes with crops in rotation
      IF(dayno==360)THEN
        IF(numrotations>0)THEN
          IF(i_sp>0) CALL croprotation_soilpoolaverage(i,numrotations,soilstate%humusP)
          IF(i_sp>0) CALL croprotation_soilpoolaverage(i,numrotations,soilstate%partP)
          IF(i_in>0) CALL croprotation_soilpoolaverage(i,numrotations,soilstate%humusN)
          IF(i_oc>0) CALL croprotation_soilpoolaverage(i,numrotations,soilstate%humusC)
        ENDIF
      ENDIF
        
      WHERE(asum>0)
        divasum = 1. / asum
      ELSEWHERE
        divasum = 0.
      ENDWHERE
      IF(conductwb)THEN
        wbstores(w_snow,i) = snowi * basin(i)%area * 1.E-3  !m3
        wbstores(w_soil1,i) = soilwateri(1) * basin(i)%area * 1.E-3  !m3
        wbstores(w_soil2,i) = soilwateri(2) * basin(i)%area * 1.E-3  !m3
        wbstores(w_soil3,i) = soilwateri(3) * basin(i)%area * 1.E-3  !m3
        !After this point soilstate%water and frozenstate%snow may not be changed, or the water balance will be off.
      ENDIF
      !Calculate subbasin mean over land area (or second or third soillayer area)
      IF(numsubstances>0)THEN
        csoili1 = csoili1 * divasum(1)
        csoili2 = csoili2 * divasum(2)
        csoili3 = csoili3 * divasum(3)
      ENDIF
      CALL calculate_class_outvar_finish(o_csoilIN,i,SUM(soilwateri)*1.E-3)  !inorgN conc soil (ug/L)
      CALL calculate_class_outvar_finish(o_csoilOC,i,SUM(soilwateri))        !orgC conc soil (mg/L)
      CALL calculate_class_outvar_finish(o_csoilT1,i,SUM(soilwateri))
      CALL calculate_class_outvar_finish(o_csoilT2,i,SUM(soilwateri))
      CALL calculate_class_outvar_finish(o_csoillayerIN(1),i,soilwateri(1)*1.E-3)  !inorgN conc soillayer 1 (ug/L)
      CALL calculate_class_outvar_finish(o_csoillayerIN(2),i,soilwateri(2)*1.E-3)
      CALL calculate_class_outvar_finish(o_csoillayerIN(3),i,soilwateri(3)*1.E-3)
      IF(outvarindex(o_snow)>0.OR.outvarindex(o_snowdens)>0) snowi = snowi * divasum(1)
      soili1 = soili1 * divasum(1)
      standsoili = standsoili * divasum(1)
      IF(outvarindex(o_crosT1)>0) CALL calculate_class_outvar_finish(o_crosT1,i,outvar(i,outvarindex(o_ros)))  !Need to be called before o_ros is finished!
      IF(outvarindex(o_crodT1)>0) CALL calculate_class_outvar_finish(o_crodT1,i,outvar(i,outvarindex(o_rod)))  !Need to be called before o_rod is finished!
      IF(outvarindex(o_cro1T1)>0) CALL calculate_class_outvar_finish(o_cro1T1,i,outvar(i,outvarindex(o_ro1)))  !Need to be called before o_ro1 is finished!
      IF(outvarindex(o_cro2T1)>0) CALL calculate_class_outvar_finish(o_cro2T1,i,outvar(i,outvarindex(o_ro2)))  !Need to be called before o_ro2 is finished!
      IF(outvarindex(o_cro3T1)>0) CALL calculate_class_outvar_finish(o_cro3T1,i,outvar(i,outvarindex(o_ro3)))  !Need to be called before o_ro3 is finished!
      CALL calculate_class_outvar_finish(o_crun,i,asum(1))
      CALL calculate_class_outvar_finish(o_ro1,i,asum(1))
      CALL calculate_class_outvar_finish(o_ro2,i,asum(2))
      CALL calculate_class_outvar_finish(o_ro3,i,asum(3))
      CALL calculate_class_outvar_finish(o_rod,i,asum(1))
      CALL calculate_class_outvar_finish(o_ros,i,asum(1))
      DO isl = 1,3
        CALL calculate_class_outvar_finish(o_pfN(isl),i,asum(isl))
        CALL calculate_class_outvar_finish(o_phN(isl),i,asum(isl))
        CALL calculate_class_outvar_finish(o_pfP(isl),i,asum(isl))
        CALL calculate_class_outvar_finish(o_phP(isl),i,asum(isl))
        CALL calculate_class_outvar_finish(o_ppP(isl),i,asum(isl))
        CALL calculate_class_outvar_finish(o_pfC(isl),i,asum(isl))
        CALL calculate_class_outvar_finish(o_phC(isl),i,asum(isl))
        CALL calculate_class_outvar_finish(o_ppT1(isl),i,asum(isl))
        CALL calculate_class_outvar_finish(o_sltmp(isl),i,asum(isl))
      ENDDO
      CALL calculate_class_outvar_finish(o_T1sf,i,asum(1))
      CALL calculate_class_outvar_finish(o_soiltmp,i,asum(1))
      CALL calculate_class_outvar_finish(o_soildenitr,i,asum(1))
      CALL calculate_class_outvar_finish(o_cropNupt,i,asum(1))
      CALL calculate_class_outvar_finish(o_degrfN,i,asum(1))
      CALL calculate_class_outvar_finish(o_soilNatm,i,asum(1))
      CALL calculate_class_outvar_finish(o_soilPatm,i,asum(1))
      CALL calculate_class_outvar_finish(o_evapsnow,i,1.)
      IF(outvarindex(o_soim)>0) outvar(i,outvarindex(o_soim)) = SUM(soilwateri)*divasum(1)
      runoff1i = runoff1i * divasum(1)
      melti = melti * divasum(1)              !OBS only snow melt on land (incl floodplain)
      CALL calculate_class_outvar_finish(o_landevap,i,asum(1))
      IF(outvarindex(o_snowdepth)>0.OR.outvarindex(o_snowdens)>0) snowdepthi = snowdepthi * divasum(1)
      soildefi = soildefi * divasum(1)
      CALL calculate_class_outvar_finish(o_soilfrost,i,asum(1))
      CALL calculate_class_outvar_finish(o_grwlevel,i,asum(1))
      CALL calculate_class_outvar_finish(o_soildef,i,asum(1))
      CALL calculate_class_outvar_finish(o_smffc,i,asum(1))
      CALL calculate_class_outvar_finish(o_smfdep,i,asum(1))
      CALL calculate_class_outvar_finish(o_smrzfdep,i,asum(1))
      CALL calculate_class_outvar_finish(o_smfpw,i,asum(1))
      CALL calculate_class_outvar_finish(o_smrzfpw,i,asum(1))
      CALL calculate_class_outvar_finish(o_snowcover,i,asum(1))
      CALL calculate_class_outvar_finish(o_snowmax,i,asum(1))
      soili2    = soili2 * divasum(2)     !Average over second soil layer area
      runoff2i  = runoff2i * divasum(2)
      soili3    = soili3 * divasum(3)     !Average over third soil layer area
      runoff3i  = runoff3i * divasum(3)
      IF(conductwb)THEN
        wbflows(w_sfallsnow,i)= wbflows(w_sfallsnow,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltsl1,i) = wbflows(w_smeltsl1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltsr,i)  = wbflows(w_smeltsr,i)  * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltmp1,i) = wbflows(w_smeltmp1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltmp2,i) = wbflows(w_smeltmp2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_smeltmp3,i) = wbflows(w_smeltmp3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_infrain,i)  = wbflows(w_infrain,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rainsr,i)   = wbflows(w_rainsr,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rainmp1,i)  = wbflows(w_rainmp1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rainmp2,i)  = wbflows(w_rainmp2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rainmp3,i)  = wbflows(w_rainmp3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_evap1,i)    = wbflows(w_evap1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_evap2,i)    = wbflows(w_evap2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_evap3,i)    = wbflows(w_evap3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_gwrunf1,i)  = wbflows(w_gwrunf1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_gwrunf2,i)  = wbflows(w_gwrunf2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_gwrunf3,i)  = wbflows(w_gwrunf3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_tile1,i)    = wbflows(w_tile1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_tile2,i)    = wbflows(w_tile2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_tile3,i)    = wbflows(w_tile3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_surfrf,i)   = wbflows(w_surfrf,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_perc1,i)    = wbflows(w_perc1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_perc2,i)    = wbflows(w_perc2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_upwell1,i)  = wbflows(w_upwell1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_upwell2,i)  = wbflows(w_upwell2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rural1,i)   = wbflows(w_rural1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rural2,i)   = wbflows(w_rural2,i) * basin(i)%area * 1.E-3   !m3/ts
        wbflows(w_rural3,i)   = wbflows(w_rural3,i) * basin(i)%area * 1.E-3   !m3/ts
        wbirrflows(w_apply1,i)  = wbirrflows(w_apply1,i) * basin(i)%area * 1.E-3   !m3/ts
        wbirrflows(w_apply2,i)  = wbirrflows(w_apply2,i) * basin(i)%area * 1.E-3   !m3/ts
        IF(glacierexist)THEN
          wbflows(w_gmeltsl1,i) = wbflows(w_gmeltsl1,i) * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_gmeltsr,i)  = wbflows(w_gmeltsr,i)  * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_gmeltmp1,i) = wbflows(w_gmeltmp1,i) * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_gmeltmp2,i) = wbflows(w_gmeltmp2,i) * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_gmeltmp3,i) = wbflows(w_gmeltmp3,i) * basin(i)%area * 1.E-3   !m3/ts
          wbflows(w_evap4,i)    = wbflows(w_evap4,i) * basin(i)%area * 1.E-3   !m3/ts
        ENDIF
        IF(conductflood)THEN
          wbfpflows(w_smtsrmrfp,i) = wbfpflows(w_smtsrmrfp,i) * basin(i)%area * 1.E-3  !m3/ts
          wbfpflows(w_rtsrmrfp,i)  = wbfpflows(w_rtsrmrfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_srftmrfp,i)  = wbfpflows(w_srftmrfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_grf1mrfp,i)  = wbfpflows(w_grf1mrfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_grf2mrfp,i)  = wbfpflows(w_grf2mrfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_grf3mrfp,i)  = wbfpflows(w_grf3mrfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_grf1olfp,i)  = wbfpflows(w_grf1olfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_grf2olfp,i)  = wbfpflows(w_grf2olfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_grf3olfp,i)  = wbfpflows(w_grf3olfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_trf1mrfp,i)  = wbfpflows(w_trf1mrfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_trf2mrfp,i)  = wbfpflows(w_trf2mrfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_trf3mrfp,i)  = wbfpflows(w_trf3mrfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_trf1olfp,i)  = wbfpflows(w_trf1olfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_trf2olfp,i)  = wbfpflows(w_trf2olfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_trf3olfp,i)  = wbfpflows(w_trf3olfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_smtsrolfp,i) = wbfpflows(w_smtsrolfp,i) * basin(i)%area * 1.E-3  !m3/ts
          wbfpflows(w_rtsrolfp,i)  = wbfpflows(w_rtsrolfp,i) * basin(i)%area * 1.E-3   !m3/ts
          wbfpflows(w_srftolfp,i)  = wbfpflows(w_srftolfp,i) * basin(i)%area * 1.E-3   !m3/ts
        ENDIF
      ENDIF
      !Forest and open land snow variables
      IF(SUM(outvarindex(267:274))>0)THEN
        IF(area_open>0.)THEN
          area_open=1./area_open
          snowvari(:,1) = snowvari(:,1) * area_open
        ELSE
          snowvari(:,1) = missing_value
        ENDIF
        IF(area_forest>0.)THEN
          area_forest=1./area_forest
          snowvari(:,2) = snowvari(:,2) * area_forest
        ELSE
          snowvari(:,2) = missing_value
        ENDIF
      ENDIF

      !>Routing in rivers and lakes
      !This means that the primary unit is m3/s below. The waterlevel of the lakes are measured in mm though.

      !Temperature in river and lakes (if T2 is simulated, temperature is calculated later)
      IF(modeloption(p_swtemperature)==0)THEN
        CALL calculate_water_temperature(i,tempi(i),riverstate,lakestate)
      ENDIF
        
      !>Local river (itype=1)
      !----------------
      itype = 1
      
      !Initiation of river calculations; 
      riverarea(itype) = 0.
      IF(slc_lriver>0)THEN
        j = slc_lriver
        a = classbasin(i,j)%part
      ELSE
        a = 0.
        j = 0
      ENDIF
        
      !Calculate inflow to local river
      qin   = runofftsi * basin(i)%area / (seconds_per_timestep * 1000.)      !m3/s
      IF(runofftsi>0.)THEN
        cin = crunofftsi(:) / runofftsi              !mg/L (crunofftsi is amount)
      ELSE
        cin = 0.
      ENDIF
      IF(conductload) Lpathway(:,1) = cin(:) * qin * seconds_per_timestep / 1000.   !Total load from land, A (kg/time step)
      IF(i_t2>0)ctemp_T2 = cin(i_t2)  !Keep T2 concentration in memory, in case water is being extracted/added local sources
         
      !Add local diffuse source to river (in)flow
      CALL add_diffuse_source_to_local_river(i,qin,cin,Lruralb,nutrientloadflows(1))
      IF(conductload) Lpathway(:,2) = cin(:) * qin * seconds_per_timestep * 1.E-3  !Load after adding rural b (kg/timestep), point B
      IF(i_t2>0) cin(i_t2) = ctemp_T2 !Restore T2 concentration, in case water was added by local sources
      IF(conductT) CALL add_tracer_point_source_to_river(i,itype,qin,cin)  !T1 and T2 point source
      IF(i_t2>0)ctemp_T2 = cin(i_t2)  !Keep T2 concentration in memory, in case water is being influenced by wetland
           
      !Check for river present (river length>0)
      IF(riverlength(itype,i)>0.)THEN

        !Calculate river wetland       
        IF(wetlandexist)  CALL calculate_river_wetland(i,itype,numsubstances,miscstate%temp5(i),miscstate%temp30(i),qin,cin,riverstate%cwetland(:,itype,i))
        IF(conductload) Lpathway(:,3) = qin * cin * seconds_per_timestep * 1.E-3   !Total load after local river wetland (kg/timestep), point C
        IF(i_t2>0) cin(i_t2) = ctemp_T2 !Restore T2 concentration, in case initial T2 of wetland will influence
        !(later, we should calculate ice and T2 in river wetlands explicitly)
           
        !Translation (delay) in river (in)flow
        CALL translation_in_river(i,itype,qin,cin,transq,transc,riverstate)
           
        !Add delayed inflow to river water volume (m3)
        CALL add_water(numsubstances,riverstate%water(itype,i),riverstate%conc(:,itype,i),transq * seconds_per_timestep,transc)
        
        !Calculate river dimensions, velocity and mean flow for use in substance processes calculation
        CALL calculate_river_characteristics(i,itype,transq,conductN.OR.conductP.OR.(i_t1>0),riverstate,riverdepth,riverarea(itype),riverQbank)
        
        !Calculate precipitation, atmospheric deposition and evaporation of river with class area
        IF(a>0)THEN      !river has classarea
          riverarea(itype) = a * basin(i)%area       !replace estimated area [m2]
          riverareakm2 = riverarea(itype)*1.E-6      ![km2]
          atmdepload1 = 0.
          atmdepload2 = 0.
            
          !Forcing data and atmospheric deposition on river 
          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
          CALL set_class_precipitation_concentration_and_load(numsubstances, &
                   riverareakm2,precorg(i),temp,prec,cprec,cprecj,atmdepload1)
          CALL add_dry_deposition_to_river(i,riverareakm2,itype,atmdepload2,    &
                vegtype(j),load(i)%indrydep,landpar(m_drypp,classdata(j)%luse), &
                riverstate)
          IF(conductload)THEN
            Latmdep(j,1,:) = Latmdep(j,1,:) + atmdepload1  !flooded floodplain + river
            Latmdep(j,2,:) = Latmdep(j,2,:) + atmdepload2  !flooded floodplain + river
          ENDIF
          CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac) !Separate snowfall/rainfall on river
          IF(i_t2>0)THEN
            CALL get_rivertempvol(i,itype,riverstate,meanrivertemp,totrivervol) !Get total river water volume and mean T2 temperature
            CALL add_T2_concentration_in_precipitation_on_water(prec,temp,  &   !Set correct T2 temperature concentration in precipitation
                  snowfall,rainfall,meanrivertemp,cprecj(i_t2),frozenstate%rivericecov(itype,i))
          ENDIF

          !Add precipitation with substances to river water
          IF(prec>0) CALL add_precipitation_to_river(i,itype,riverarea(itype),prec,cprecj,riverstate)
          
          !Evaporation of river, taking partial ice-cover into account if icemodel is enabled
          CALL calculate_potential_evaporation(i,j,temp,epot,radexti(i),swrad,netrad,actvap,satvap,wind,epotsnow)
          epot = epot * basincevpcorr(i)
          IF(modeloption(p_lakeriverice)>0)  epot = epot * (1. - frozenstate%rivericecov(itype,i))
          IF(epot>0.)THEN
            CALL calculate_river_evaporation(i,j,itype,numsubstances,riverarea(itype),temp,epot,evapr,cevapr,riverstate)
          ELSE
            evapr = 0.  !current evaporation
          ENDIF
            
          !Accumulate output variables for mean over subbasin (accumulation for cevapi)
          CALL calculate_class_outvar_add(o_ctmp,i,a,temp)
          CALL calculate_class_outvar_add(o_rainfall,i,a,rainfall)
          CALL calculate_class_outvar_add(o_snowfall,i,a,snowfall)
          corrpreci = corrpreci + prec*a
          !corrtempi = corrtempi + temp*a
          !rainfalli = rainfalli + rainfall*a
          !snowfalli = snowfalli + snowfall*a
          CALL calculate_class_outvar_add(o_epot,i,a,epot)
          !epoti = epoti+epot*a
          eacti = eacti + evapr*a
          IF(i_t1>0) CALL calculate_class_outvar_add(o_cevapT1,i,a,cevapr(i_t1)*evapr)
!          IF(i_t1>0) cevapi = cevapi+cevapr(i_t1)*evapr*a
          icpevapi = icpevapi + icpevap*a
          !Accumulate outvar for precipitation and evaporation on water (probably excluding floodplain water)
          area_water = area_water + a
          CALL calculate_class_outvar_add(o_clrp,i,a,rainfall+snowfall)
          CALL calculate_class_outvar_add(o_clre,i,a,evapr)

          !Set water balance output
          IF(conductwb)THEN
            wbflows(w_piriver,i) = prec * riverarea(itype)*1.E-3     !m3
            wbflows(w_eiriver,i) = evapr * riverarea(itype)*1.E-3
          ENDIF
        ENDIF !(a>0)    

        !Calculate river ice and T2 processes of local river
        IF(modeloption(p_lakeriverice)>0)THEN
            
          !If no river class area: Assume same precipitation and temperature conditions as for lake class
          IF(a==0)THEN
            j=slc_ilake
            IF(j==0) j=slc_olake
            CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
            CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac)
          ENDIF
              
          !Calculate river T2 and ice processes
          riverfreezeupday(itype)=0 !could by non-array
          riverbreakupday(itype)=0  !actually not used here
          CALL T2_processes_in_river(i,itype,temp,swrad,riversurftemp,riverarea(itype),frozenstate,riverstate,riverfreezeupday,freezeuparea)
          CALL ice_processes_in_river(i,itype,classdata(j)%luse,snowfall,temp, & 
                                        riversurftemp,riverarea(itype),swrad,         &
                                        frozenstate,riverstate,riverfreezeupday,riverbreakupday,freezeuparea)
          IF(modeloption(p_swtemperature)==1)  CALL set_water_temperature(itype*2-1,i,riverstate,lakestate)
        ENDIF

        !Calculate substances processes in river
        CALL np_processes_in_river(i,itype,riverarea(itype),riverdepth,transq,riverQbank,       &
                (2.-incorr)*genpar(m_denitwr),   &
                (2.-incorr)*genpar(m_denitwrl),genpar(m_hsatINwater), &
                lakedatapar(lakedataparindex(i,itype),m_ldwprodn),    &
                lakedatapar(lakedataparindex(i,itype),m_ldwprodp),    &
                genpar(m_hsatTP),genpar(m_sedexp),genpar(m_limsedpp),riverstate)
        CALL oc_processes_in_river(i,itype,riverarea(itype),riverdepth,   &
                lakedatapar(lakedataparindex(i,itype),m_ldwprodc),  &
                genpar(m_hsatTP),genpar(m_limsedpp),riverstate)
        CALL tracer_processes_in_river(i,itype,riverarea(itype),riverdepth,transq,riverQbank,riverstate)
        
        !Calculate and remove outflow from river water volume
        dampq = (riverrc(itype,i)*(riverstate%water(itype,i) - deadriver(itype,i)))   !m3/ts
        IF(dampq<0.) dampq=0.     !safety
        dampc = riverstate%conc(:,itype,i)
        CALL remove_water(riverstate%water(itype,i),numsubstances,riverstate%conc(:,itype,i),dampq,dampc,status)
        IF(status.NE.0) CALL error_remove_water(errstring(1),basin(i)%subid,i,itype)
        IF(conductwb) wbflows(w_irtomr,i) = dampq  !in case of no ilake
        IF(conductload) Lpathway(:,4) = dampc * dampq * 1.E-3   !Load at point C (downstream of local river), point D
        dampq = dampq / seconds_per_timestep    !m3/s

      ELSE
        !No river, transport inflow to next calculation step
        dampq = qin   !m3/s
        dampc = cin
          
        !Set some output variables
        IF(conductload)THEN
          Lpathway(:,3) = Lpathway(:,2)  !Total load (kg/timestep), point C
          Lpathway(:,4) = Lpathway(:,2)  !Load at point D
        ENDIF
        IF(conductwb) wbflows(w_irtomr,i) = dampq * seconds_per_timestep
          
      ENDIF !river present
          

      !>Local internal lake (itype=1)
      !----------------
      !Initiation of lake calculations
      lakearea(itype) = 0.
      a = 0.
      IF(slc_ilake>0)THEN
        j = slc_ilake
        a = classbasin(i,j)%part
      ENDIF

      !Lake calculations
      IF(a>0)THEN      !local lake exist
        !TODO: CALL lakemodel_0()
        lakearea(itype) = a * basin(i)%area       ![m2]
        lakeareakm2 = lakearea(itype)*1.E-6       ![km2]
        qunitfactor = seconds_per_timestep * 1000. / lakearea(itype)     !m3/s->mm/timestep
        atmdepload1 = 0.
        atmdepload2 = 0.
            
        !Forcing data and atmospheric deposition on local lake 
        CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
        CALL set_class_precipitation_concentration_and_load(numsubstances,  &
                lakeareakm2,precorg(i),temp,prec,cprec,cprecj,atmdepload1)
        CALL add_dry_deposition_to_lake(i,lakeareakm2,itype,atmdepload2,   &
              vegtype(j),load(i)%indrydep,landpar(m_drypp,classdata(j)%luse), &
              lakestate)
        IF(conductload)THEN
          Latmdep(j,1,:) = Latmdep(j,1,:) + atmdepload1  !flooded floodplain + lake
          Latmdep(j,2,:) = Latmdep(j,2,:) + atmdepload2  !flooded floodplain + lake
        ENDIF
        CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac) 
        IF(i_t2>0) CALL add_T2_concentration_in_precipitation_on_water(prec,temp, &
                          snowfall,rainfall,lakestate%uppertemp(itype,i),cprecj(i_t2), &
                          frozenstate%lakeicecov(itype,i))

        !Add precipitation with substances to lake water
        IF(prec>0) CALL add_water(numsubstances,lakestate%water(itype,i),lakestate%conc(:,itype,i),prec,cprecj)
          
        !Evaporation of lake, taking partial lakeice cover into account, if lakeice model is enabled
        CALL calculate_potential_evaporation(i,j,temp,epot,radexti(i),swrad,netrad,actvap,satvap,wind,epotsnow)
        epot = epot * basincevpcorr(i)
        IF(modeloption(p_lakeriverice)>0) epot = epot * (1. - frozenstate%lakeicecov(itype,i))
        IF(epot>0.)THEN
          CALL calculate_actual_lake_evaporation(i,j,itype,numsubstances,temp,epot,evapl,cevapl,lakestate)
        ELSE
          evapl = 0.
        ENDIF
          
        !Accumulate output variables for mean over subbasin (accumulation for cevapi)
        CALL calculate_class_outvar_add(o_ctmp,i,a,temp)
        CALL calculate_class_outvar_add(o_rainfall,i,a,rainfall)
        CALL calculate_class_outvar_add(o_snowfall,i,a,snowfall)
        corrpreci = corrpreci + prec*a
        !corrtempi = corrtempi + temp*a
        !rainfalli = rainfalli + rainfall*a
        !snowfalli = snowfalli + snowfall*a 
        CALL calculate_class_outvar_add(o_epot,i,a,epot)
        !epoti = epoti+epot*a
        eacti = eacti + evapl*a
        IF(i_t1>0) CALL calculate_class_outvar_add(o_cevapT1,i,a,cevapl(i_t1)*evapl)
        !IF(i_t1>0) cevapi = cevapi+cevapl(i_t1)*evapl*a
        icpevapi = icpevapi + icpevap*a
        !Accumulate outvar for precipitation and evaporation on water (probably excluding floodplain water)
        area_water = area_water + a
        CALL calculate_class_outvar_add(o_clrp,i,a,rainfall+snowfall)
        CALL calculate_class_outvar_add(o_clre,i,a,evapl)

        !Set water balance output
        IF(conductwb)THEN
          wbflows(w_pilake,i) = prec * lakearea(itype)*1.E-3     !m3
          wbflows(w_eilake,i) = evapl * lakearea(itype)*1.E-3
        ENDIF

        !Calculate and add inflow to local lake
        inflowpart = basin(i)%ilakecatch
        qin = dampq * inflowpart         ![m3/s]
        qinmm = qin * qunitfactor        ![mm/d]
        cin = dampc
        IF(conductT) CALL add_tracer_point_source_to_lake(i,itype,qin,cin)
        IF(conductload)THEN
          Lpathway(:,6) = cin * qin * seconds_per_timestep * 1.E-3  !Load at point F, local flow to ilake (kg/timestep)
          Lpathway(:,5) = Lpathway(:,4) - Lpathway(:,6)             !Load at point E, bypassed local flow
        ENDIF
        IF(conductwb .AND. qin>0)THEN
          wbflows(w_irtoil,i) = qin * seconds_per_timestep
          wbflows(w_irtomr,i) = wbflows(w_irtomr,i) - wbflows(w_irtoil,i)
        ENDIF
        CALL add_water(numsubstances,lakestate%water(itype,i),lakestate%conc(:,itype,i),qinmm,cin)
            
        !Calculate lake epilimnion depth
        IF(i_t2>0) CALL calculate_lake_epidepth(lakearea(itype),prec,evapl,qinmm,epidepth)
          
        !Calculate lake tracer (T1), temperature (T2) and ice processes
        lakefreezeupday(itype)=0 !could by non-array
        lakebreakupday(itype)=0  !actually not used here
        CALL tracer_processes_in_lake(i,itype,lakearea(itype),lakestate)
        IF(i_t2>0) CALL T2_processes_in_lake(i,itype,temp,swrad,lakesurftemp, &
                              lakearea(itype),epidepth,frozenstate,lakestate, &
                              lakefreezeupday,freezeuparea)
        IF(modeloption(p_lakeriverice)>0) CALL ice_processes_in_lake(i,itype,classdata(j)%luse,snowfall,temp, & 
                                      lakesurftemp,swrad,frozenstate,lakestate, & 
                                      lakefreezeupday,lakebreakupday,lakearea(itype),epidepth,freezeuparea) 
        IF(modeloption(p_swtemperature)==1)  CALL set_water_temperature(itype*2,i,riverstate,lakestate)
          
        !Calculate substances processes in lakes (in flooded floodplain water to be added here or in soilmodel4)
        CALL np_processes_in_lake(i,itype,lakearea(itype),    &
              (2.-incorr)*lakedatapar(lakedataparindex(i,itype),m_lddenitwl),genpar(m_hsatINwater),  &
              lakedatapar(lakedataparindex(i,itype),m_ldwprodn),   &
              lakedatapar(lakedataparindex(i,itype),m_ldwprodp),genpar(m_hsatTP),   &
              (2.-oncorr)*lakedatapar(lakedataparindex(i,itype),m_ldsedon),    &
              lakedatapar(lakedataparindex(i,itype),m_ldsedpp),genpar(m_limsedon),genpar(m_limsedpp),lakestate)
        CALL oc_processes_in_lake(i,itype,lakearea(itype),lakedatapar(lakedataparindex(i,itype),m_ldwprodc),    &
              genpar(m_hsatTP),genpar(m_limsedpp),lakedatapar(lakedataparindex(i,itype),m_ldsedoc),lakestate)
          
        !Calculate and remove outflow from lake
        CALL calculate_ilake_outflow(i,basin(i)%subid,numsubstances,  &
                 conductN.OR.conductP.OR.conductC.OR.conductT,qin,    &
                 lakearea(itype),qunitfactor,lakeoutflow(itype),concout,  &
                 lakestate)
          
        !Flow between fast and slow lake parts for divided lake
        CALL calculate_flow_within_lake(i,itype,basin(i)%subid,lakestate)
          
        !Calculate load from local lake
        IF(conductload) Lpathway(:,7) = lakeoutflow(itype) * concout * seconds_per_timestep * 1.E-3 !Load at point G, outflow ilake (kg/timestep)
        IF(conductwb) wbflows(w_iltomr,i) = lakeoutflow(itype) * seconds_per_timestep
            
        !Lake water stage for print out 
        lakewst(itype) = lakestate%water(itype,i)
        IF(ALLOCATED(lakestate%slowwater)) lakewst(itype) = lakewst(itype) + lakestate%slowwater(itype,i)
          
      ELSE    !No local lake, outflow = inflow
        inflowpart = 2. !to avoid adding lakeoutflow below
        lakeoutflow(itype) = dampq       !m3/s
        concout = dampc
        IF(conductload)THEN
          Lpathway(:,5)  = Lpathway(:,4)           !Load at point E, bypassed local flow
          Lpathway(:,6)  = 0.                      !Load at point F
          Lpathway(:,7)  = 0.                      !Load at point G
        ENDIF
        lakewst(itype) = missing_value   
        lakesurftemp(itype) = missing_value
      ENDIF   !a>0

      !Finalize local contribution from subbasin
      IF(inflowpart<1.)THEN        !Add flow from local river (bypassing lake) and local lake
        CALL add_water(numsubstances,lakeoutflow(itype),concout,dampq*(1.-inflowpart),dampc)
      ENDIF
      IF(conductload) Lpathway(:,8)  = lakeoutflow(itype) * concout * seconds_per_timestep * 1.E-3  !Load at point H, outflow ilake + bypassed water (kg/timestep)
      IF(doupdate(i_tploccorr))  CALL apply_nutrientcorr(updatetploccorr(i),concout(i_sp),concout(i_pp))          
      IF(doupdate(i_tnloccorr))  CALL apply_nutrientcorr(updatetnloccorr(i),concout(i_in),concout(i_on))          
      clakeoutflow(:,itype) = concout


      !>Main river (itype=2)
      !----------------
      itype = 2
      olakewstold=lakestate%water(2,i)  !Remember olake waterstage for last time step, to be used in war-updating
      IF(conductN.OR.conductP.OR.conductC.OR.conductT) olakewstold=olakewstold+lakestate%slowwater(2,i)

      !Initiation of river calculations
      riverarea(itype) = 0.
      fpfrac = 0.             !default is no floodplain
      IF(slc_mriver>0)THEN
        j = slc_mriver
        a = classbasin(i,j)%part
      ELSE
        a = 0.
        j = 0
      ENDIF
      
      !Calculate inflow to river
      qin = lakeoutflow(1) + accinflow(i)         ![m3/s]
      IF(qin>0)THEN
        cin = (clakeoutflow(:,1)*lakeoutflow(1) + acccinflow(:,i)) / qin
      ELSE
        cin = 0.
      ENDIF
      IF(i_t2>0)ctemp_T2 = cin(i_t2)  !Keep T2 concentration in memory, in case water is being extracted/added by irrigation and/or sources
         
      !Irrigation water removal from sources and calculate irrigation to be applied the next day
      IF(doirrigation)THEN
        CALL calculate_irrigation(i,naquifers,genpar(m_regirr),regpar(m_pirrs,basin(i)%parregion),      &
                regpar(m_pirrg,basin(i)%parregion),qin,lakestate%water(1:2,i),riverstate%water(itype,i),aquiferstate%water,pwneedi,irrevap,gwremi,     & !for floodplain safety
                ldremi,lrremi,rsremi,cin,regpar(m_cirrsink,basin(i)%parregion),irrsinkmass,  &
                genpar(m_irrcomp),lakestate%conc(:,1,i),lakestate%conc(:,2,i),riverstate%conc(:,itype,i),aquiferstate%conc, &
                soilstate,miscstate%nextirrigation,miscstate%cnextirrigation,irrigationflows,regionalirrflows,regionalirrevap, &
                aqirrloss)
        IF(conductwb)THEN
          wbirrflows(w_wdfromdg,i) = irrigationflows(1)
          wbirrflows(w_wdfromil,i) = irrigationflows(2)
          wbirrflows(w_wdfromol,i) = irrigationflows(3)
          wbirrflows(w_wdfrommr,i) = irrigationflows(4)
          wbirrflows(w_wdoutside,i) = irrigationflows(6)
          wbirrflows(w_rgrwtoir,i) = irrigationflows(7)
          wbirrflows(w_evapirrc,i) = wbirrflows(w_evapirrc,i) + irrigationflows(5)  !local sources
          wbirrflows(w_evapirrc,:) = wbirrflows(w_evapirrc,:) + regionalirrevap(:)  !regional sources, local network
        ENDIF
      ENDIF
     
      !Add point sources source, remove abstractions and add aquifer inflow to main river inflow          
      CALL point_abstraction_from_main_river_inflow(i,itype,qin,riverstate,nutrientloadflows(3))
      CALL add_point_sources_to_main_river(i,qin,cin,Lpoints,nutrientloadflows(2))
      IF(i_t2>0) cin(i_t2) = ctemp_T2  !Restore T2 concentration, in case water was extracted/added by irrigation and/or point sources
      IF(conductT) CALL add_tracer_point_source_to_river(i,itype,qin,cin)  !T1 and T2 point source
      IF(i_t2>0)ctemp_T2 = cin(i_t2)  !Keep T2 concentration in memory, in case water is being added by aquifers
      IF(modeloption(p_deepgroundwater)==2)THEN
        CALL add_aquifer_flow_to_river(i,numsubstances,qin,cin,grwtomr,grwloadmr)
        IF(conductwb) wbflows(w_rgrwtomr,i) = grwtomr
        IF(outvarindex(241)>0) outvar(i,outvarindex(241)) = grwtomr
        IF(conductload) Lgrwmr = grwloadmr
      ENDIF
      IF(conductload) Lpathway(:,9) = cin(:) * qin * seconds_per_timestep * 1.E-3  !Load after adding point sources and aquifer inflow (kg/timestep), point I
      IF(i_t2>0) cin(i_t2) = ctemp_T2  !Restore T2 concentration, in case water was added by aquifers
           
      !Check for river present (river length>0)
      IF(riverlength(itype,i)>0.)THEN

        !Calculate river wetland       
        IF(wetlandexist)  CALL calculate_river_wetland(i,itype,numsubstances,miscstate%temp5(i),miscstate%temp30(i),qin,cin,riverstate%cwetland(:,itype,i))
        IF(conductload)  Lpathway(:,10) = qin * cin * seconds_per_timestep * 1.E-3  !Total load after main river wetland (kg/timestep), point J
        IF(i_t2>0) cin(i_t2) = ctemp_T2   !Restore T2 concentration, in case water was extracted/added by wetlands
        !(later, we should calculate ice and T2 in river wetlands explicitly)
           
        !Translation (delay) in river (in)flow
        CALL translation_in_river(i,itype,qin,cin,transq,transc,riverstate)
           
        !Add river inflow to river water volume
        CALL add_water(numsubstances,riverstate%water(itype,i),riverstate%conc(:,itype,i),transq * seconds_per_timestep,transc)
        
        !Calculate river dimensions, velocity and mean flow for use in substance processes calculation
        CALL calculate_river_characteristics(i,itype,transq,conductN.OR.conductP.OR.(i_t1>0),riverstate,riverdepth,riverarea(itype),riverQbank)
        
        !Calculate precipitation, atmospheric deposition and evaporation of river with area
        IF(a>0)THEN   !river has area
          riverarea(itype) = a * basin(i)%area        !default [m2]
          IF(conductflood)THEN
            IF(floodindex(i)>0)THEN
              IF(flooding(floodindex(i))%fpfmr>0.)THEN
                fpfrac = flooding(floodindex(i))%fpfmr !floodplain fraction of the main river area
                riverarea(itype) = a * basin(i)%area * (1.-fpfrac)   !Adjusted [m2]
                riverfparea = a * basin(i)%area * fpfrac             !Area of floodplain part [m2]
              ENDIF
            ENDIF
          ENDIF
          
          riverareakm2 = riverarea(itype)*1.E-6      ![km2}
          atmdepload1 = 0.
          atmdepload2 = 0.
            
          !Forcing data and atmospheric deposition on river 
          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
          CALL set_class_precipitation_concentration_and_load(numsubstances, &
                   riverareakm2,precorg(i),temp,prec,cprec,cprecj,atmdepload1)
          CALL add_dry_deposition_to_river(i,riverareakm2,itype,atmdepload2,    &
                vegtype(j),load(i)%indrydep,landpar(m_drypp,classdata(j)%luse), &
                riverstate)
          IF(conductload)THEN
            Latmdep(j,1,:) = Latmdep(j,1,:) + atmdepload1  !flooded floodplain + river
            Latmdep(j,2,:) = Latmdep(j,2,:) + atmdepload2  !flooded floodplain + river
          ENDIF
          CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac)
          IF(i_t2>0)THEN
            CALL get_rivertempvol(i,itype,riverstate,meanrivertemp,totrivervol)  !Get total river water volume and mean T2 temperature
            CALL add_T2_concentration_in_precipitation_on_water(prec,temp,  &  !Set correct T2 temperature concentration in precipitation
                  snowfall,rainfall,meanrivertemp,cprecj(i_t2),frozenstate%rivericecov(itype,i))
          ENDIF
          IF(prec>0) CALL add_precipitation_to_river(i,itype,riverarea(itype),prec,cprecj,riverstate)
          
          !Evaporation of river, taking partial ice-cover into account if icemodel is enabled
          CALL calculate_potential_evaporation(i,j,temp,epot,radexti(i),swrad,netrad,actvap,satvap,wind,epotsnow)
          epot = epot * basincevpcorr(i)
          IF(modeloption(p_lakeriverice)>0) epot = epot * (1. - frozenstate%rivericecov(itype,i))
          IF(epot>0.)THEN
            CALL calculate_river_evaporation(i,j,itype,numsubstances,riverarea(itype),temp,epot,evapr,cevapr,riverstate)
          ELSE
            evapr = 0.
          ENDIF
            
          !Accumulate output variables for mean over subbasin (accumulation for cevapi)
          CALL calculate_class_outvar_add(o_ctmp,i,a,temp)
          CALL calculate_class_outvar_add(o_rainfall,i,a,rainfall)
          CALL calculate_class_outvar_add(o_snowfall,i,a,snowfall)
          corrpreci = corrpreci + prec*a
          !corrtempi = corrtempi + temp*a
          !rainfalli = rainfalli + rainfall*a
          !snowfalli = snowfalli + snowfall*a
          CALL calculate_class_outvar_add(o_epot,i,a*(1.-fpfrac),epot)
          !epoti = epoti+epot*(1.-fpfrac)*a
          eacti = eacti + (evapr * (1.-fpfrac)) * a
          IF(i_t1>0) CALL calculate_class_outvar_add(o_cevapT1,i,a*(1.-fpfrac),cevapr(i_t1)*evapr)
          !IF(i_t1>0) cevapi = cevapi+(cevapr(i_t1)*evapr*(1.-fpfrac))*a
          icpevapi = icpevapi + icpevap*(1.-fpfrac)*a

          !Accumulate outvar for precipitation and evaporation on water (probably excluding floodplain water)
          CALL calculate_class_outvar_add(o_clrp,i,a*(1.-fpfrac),rainfall+snowfall)
          CALL calculate_class_outvar_add(o_clre,i,a*(1.-fpfrac),evapr)
          area_water = area_water + a*(1.-fpfrac)
          
          !Set water balance output
          IF(conductwb)THEN
            wbflows(w_pmriver,i) = prec * riverarea(itype)*1.E-3
            wbflows(w_emriver,i) = evapr * riverarea(itype)*1.E-3
          ENDIF
        ENDIF !(a>0)    

        !Calculate river ice and T2 processes of main river
        IF(modeloption(p_lakeriverice)>0)THEN
            
          !If no river class: Assume same precipitation and temperature conditions as for lake class
          IF(a==0)THEN
            j=slc_olake
            CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                  temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
            CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac)
          ENDIF
              
          !Calculate main river ice and T2 processes
          riverfreezeupday(itype)=0 !could by non-array
          riverbreakupday(itype)=0  !actually not used here
          CALL T2_processes_in_river(i,itype,temp,swrad,riversurftemp,riverarea(itype),frozenstate,riverstate,riverfreezeupday,freezeuparea)
          CALL ice_processes_in_river(i,itype,classdata(j)%luse,snowfall,temp, & 
                                        riversurftemp,riverarea(itype),swrad,         &
                                        frozenstate,riverstate,riverfreezeupday,riverbreakupday,freezeuparea)
          IF(modeloption(p_swtemperature)==1) CALL set_water_temperature(itype*2-1,i,riverstate,lakestate)
          
          !Reset j to value before lakeriverice calculation
          IF(a==0)THEN
            IF(slc_mriver>0)THEN
              j = slc_mriver
            ELSE
              j = 0
            ENDIF
          ENDIF
        ENDIF

        !Regional groundwater flow from main river
        IF(modeloption(p_deepgroundwater)==2)THEN
          CALL calculate_river_groundwaterflow_removal(i,j,basin(i)%subid,numsubstances,riverstate,grwout2)
          IF(conductwb) wbflows(w_rgrwofmr,i) = grwout2
          IF(fpfrac>0 .AND. a>0)THEN
            IF(genpar(m_optonoff).LE.0)THEN
              interflowpar(3) = flooding(floodindex(i))%fymmr
            ELSE
              interflowpar(3) = genpar(m_opt7)
            ENDIF
            CALL calculate_floodplain_waterlevel(miscstate%floodwater(1,i),riverfparea,interflowpar(3),ffpwl,ffparea)
            CALL calculate_river_floodplain_groundwaterflow_removal(i,j,basin(i)%subid,numsubstances,miscstate,grwout2)
            IF(conductwb) wbfpflows(w_rgrwofmrfp,i) = grwout2
          ENDIF
        ENDIF
        
        !Alternative location of abstraction of water from main river           
        CALL point_abstraction_from_main_river(i,itype,riverstate,nutrientloadflows(3))
          
        !Calculate substances processes in river
        CALL np_processes_in_river(i,itype,riverarea(itype),riverdepth,transq,riverQbank,       &
                (2.-incorr)*genpar(m_denitwr),   &
                (2.-incorr)*genpar(m_denitwrl),genpar(m_hsatINwater), &
                lakedatapar(lakedataparindex(i,itype),m_ldwprodn),    &
                lakedatapar(lakedataparindex(i,itype),m_ldwprodp),    &
                genpar(m_hsatTP),genpar(m_sedexp),genpar(m_limsedpp),riverstate)
        CALL oc_processes_in_river(i,itype,riverarea(itype),riverdepth,   &
                lakedatapar(lakedataparindex(i,itype),m_ldwprodc),  &
                genpar(m_hsatTP),genpar(m_limsedpp),riverstate)
        CALL tracer_processes_in_river(i,itype,riverarea(itype),riverdepth,transq,riverQbank,riverstate)
        
        !Calculate interflow between main river and floodplain
        IF(fpfrac>0 .AND. a>0)THEN
          IF(genpar(m_optonoff).LE.0)THEN   !Get current interflow parameter values
            interflowpar(1) = flooding(floodindex(i))%flmrr    !flmr
            interflowpar(2) = flooding(floodindex(i))%flmrp    !flfp
            interflowpar(3) = flooding(floodindex(i))%fymmr
            interflowpar(4) = flooding(floodindex(i))%rcr2fp
            interflowpar(5) = flooding(floodindex(i))%rcfp2r
          ELSE
            interflowpar(1) = genpar(m_opt3)   !flmr
            interflowpar(2) = genpar(m_opt4)   !flfp
            interflowpar(3) = genpar(m_opt7)
            interflowpar(4) = genpar(m_opt5)
            interflowpar(5) = genpar(m_opt8)
          ENDIF
          CALL calculate_waterbody_floodplain_interflow(i,riverfparea,riverarea(itype),interflowpar,  &
                  miscstate%floodwater(1,i),miscstate%cfloodwater(:,1,i),riverstate%water(itype,i), &
                  riverstate%conc(:,itype,i),flooddepth,flooddegree,qmrflood)
          IF(outvarindex(323)>0) outvar(i,outvarindex(323)) = flooddepth
          IF(outvarindex(325)>0) outvar(i,outvarindex(325)) = flooddegree
          IF(conductwb)THEN
            IF(qmrflood>0.) wbfpflows(w_mrtofp,i) = qmrflood
            IF(qmrflood<0.) wbfpflows(w_fptomr,i) = -qmrflood
          ENDIF
        ELSE  !no floodplain
          IF(outvarindex(323)>0) outvar(i,outvarindex(323)) = missing_value
          IF(outvarindex(325)>0) outvar(i,outvarindex(325)) = 0.
        ENDIF

        !Calculate and remove outflow from river water volume
        dampq = (riverrc(itype,i)*(riverstate%water(itype,i) - deadriver(itype,i)))
        IF(dampq<0.) dampq=0.     !safe for irrigation etc
        dampc = riverstate%conc(:,itype,i)
        CALL remove_water(riverstate%water(itype,i),numsubstances,riverstate%conc(:,itype,i),dampq,dampc,status)
        IF(status.NE.0) CALL error_remove_water(errstring(1),basin(i)%subid,i,itype)
        IF(conductwb) wbflows(w_mrtool,i) = dampq
        IF(conductload) Lpathway(:,11) = dampc * dampq * 1.E-3  !Load at point I (downstream of main river), point K
        dampq = dampq / seconds_per_timestep
            
      ELSE
        !No river, transport inflow to next calculation step (outlet lake)
        dampq = qin   ![m3/s]
        dampc = cin
          
        !Set output variables
        IF(conductload)THEN
          Lpathway(:,10) = Lpathway(:,9) !Total load  (kg/timestep), point J
          Lpathway(:,11) = Lpathway(:,9)  !Load at point K
        ENDIF
        IF(outvarindex(323)>0) outvar(i,outvarindex(323)) = missing_value !Floodplain output variables
        IF(outvarindex(325)>0) outvar(i,outvarindex(325)) = 0.
        !outvar(i,323) = missing_value   !Floodplain output variables        
        !outvar(i,325) = 0.
        IF(conductwb) wbflows(w_mrtool,i) = dampq * seconds_per_timestep
          
      ENDIF !river present
   
          
      !>Outlet lake (itype=2)
      !----------------
      !Initiation of lake calculations
      lakearea(itype) = 0.
      a = 0.
      IF(slc_olake>0)THEN
        j = slc_olake
        a = classbasin(i,j)%part
      ENDIF

      !Lake exist:
      IF(a>0)THEN
        fpfrac = 0.   !default is no floodplain
        lakearea(itype) = a * basin(i)%area     ![m2]
        IF(conductflood)THEN
          IF(floodindex(i)>0)THEN
            IF(flooding(floodindex(i))%fpfol>0.)THEN
              fpfrac = flooding(floodindex(i))%fpfol          !floodplain fraction of outlet lake area
              lakearea(itype) = a * basin(i)%area * (1.-fpfrac)    !m2
              lakefparea = a * basin(i)%area * fpfrac       !Area of floodplain part m2
            ENDIF
          ENDIF
        ENDIF
        lakeareakm2 = lakearea(itype)*1.E-6       ![km2}
        qunitfactor = seconds_per_timestep * 1000. / lakearea(itype)     !m3/s->mm/timestep
        atmdepload1 = 0.
        atmdepload2 = 0.
            
        !Forcing data and atmospheric deposition on lake 
        CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
        CALL set_class_precipitation_concentration_and_load(numsubstances, &
                 lakeareakm2,precorg(i),temp,prec,cprec,cprecj,atmdepload1)
        CALL add_dry_deposition_to_lake(i,lakeareakm2,itype,atmdepload2,   &
              vegtype(j),load(i)%indrydep,landpar(m_drypp,classdata(j)%luse), &
              lakestate)
        IF(conductload)THEN
          Latmdep(j,1,:) = Latmdep(j,1,:) + atmdepload1  !flooded floodplain + lake
          Latmdep(j,2,:) = Latmdep(j,2,:) + atmdepload2  !flooded floodplain + lake
        ENDIF
        CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac) 
        IF(i_t2>0) CALL add_T2_concentration_in_precipitation_on_water(prec,temp, &
                          snowfall,rainfall,lakestate%uppertemp(itype,i),cprecj(i_t2), &
                          frozenstate%lakeicecov(itype,i))
        IF(prec>0) CALL add_water(numsubstances,lakestate%water(itype,i),lakestate%conc(:,itype,i),prec,cprecj)
          
        !Evaporation of lake, taking partial ice cover into account, if lakeriverice model is enabled
        CALL calculate_potential_evaporation(i,j,temp,epot,radexti(i),swrad,netrad,actvap,satvap,wind,epotsnow)
        epot = epot * basincevpcorr(i)
        IF(modeloption(p_lakeriverice)>0) epot = epot * (1. - frozenstate%lakeicecov(itype,i))
        IF(epot>0.)THEN
          CALL calculate_actual_lake_evaporation(i,j,itype,numsubstances,temp,epot,evapl,cevapl,lakestate)
        ELSE
          evapl = 0.
        ENDIF
        
        !Accumulate output variables for mean over subbasin (accumulation for cevapi)
        CALL calculate_class_outvar_add(o_ctmp,i,a,temp)
        CALL calculate_class_outvar_add(o_rainfall,i,a,rainfall)
        CALL calculate_class_outvar_add(o_snowfall,i,a,snowfall)
        corrpreci = corrpreci + prec*a
        !corrtempi = corrtempi + temp*a
        !rainfalli = rainfalli + rainfall*a
        !snowfalli = snowfalli + snowfall*a 
        CALL calculate_class_outvar_add(o_epot,i,a*(1.-fpfrac),epot)
        !epoti = epoti+epot*(1.-fpfrac)*a
        eacti = eacti + (evapl * (1.-fpfrac))*a
        IF(i_t1>0) CALL calculate_class_outvar_add(o_cevapT1,i,a*(1.-fpfrac),cevapl(i_t1)*evapl)
        !IF(i_t1>0) cevapi = cevapi+(cevapl(i_t1)*evapl*(1.-fpfrac))*a
        !icpevapi = icpevapi + icpevap*a    !Bugg saknar (1-fpfrac)
        icpevapi = icpevapi + icpevap*(1.-fpfrac)*a
        qcinfli = qcinfli + (prec-evapl)*lakearea(itype)/seconds_per_timestep/1000.    !m3/s

        !Accumulate outvar for precipitation and evaporation on water (probably excluding floodplain water)
        CALL calculate_class_outvar_add(o_clrp,i,a*(1.-fpfrac),rainfall+snowfall)
        CALL calculate_class_outvar_add(o_clre,i,a*(1.-fpfrac),evapl)
        area_water = area_water + a*(1.-fpfrac)
          
        !Set water balance output
        IF(conductwb)THEN
          wbflows(w_polake,i) = prec * lakearea(itype)*1.E-3
          wbflows(w_eolake,i) = evapl * lakearea(itype)*1.E-3
        ENDIF

        !Calculate and add inflow to outlet lake (from main river)
        qin = dampq                       ![m3/s]
        qcinfli = qcinfli + qin
        qinmm = qin * qunitfactor         ![mm/ts]
        cin = dampc
        IF(conductT) CALL add_tracer_point_source_to_lake(i,itype,qin,cin)
        IF(conductload) Lpathway(:,12) = cin * qin * seconds_per_timestep * 1.E-3     !Load at point L inflow to olake (kg/timestep)
        CALL add_water(numsubstances,lakestate%water(itype,i),lakestate%conc(:,itype,i),qinmm,cin)
            
        !Add regional groundwater inflow to olake
        IF(modeloption(p_deepgroundwater)==1)THEN
          CALL add_regional_groundwater_flow_to_olake(i,itype,numsubstances,  &
                qunitfactor,lakestate%water(itype,i),lakestate%conc(:,itype,i),  &
                qcinfli,grwloadlake,grwtool)
          IF(conductload)THEN
            Lgrwol(:) = grwloadlake
            Lpathway(:,12) = Lpathway(:,12) + grwloadlake
          ENDIF
          IF(conductwb) wbflows(w_rgrwtool,i) = grwtool
        ENDIF

        !Calculate lake epilimnion depth
        IF(i_t2>0) CALL calculate_lake_epidepth(lakearea(itype),prec,evapl,qinmm,epidepth)
          
        !Calculate lake tracer (T1), temperature (T2) and ice processes
        lakefreezeupday(itype)=0 !could by non-array
        lakebreakupday(itype)=0  !actually not used here
        IF(i_t2>0) CALL T2_processes_in_lake(i,itype,temp,swrad,lakesurftemp, &
                              lakearea(itype),epidepth,frozenstate,lakestate, &
                              lakefreezeupday,freezeuparea)
        CALL tracer_processes_in_lake(i,itype,lakearea(itype),lakestate)
        IF(modeloption(p_lakeriverice)>0) CALL ice_processes_in_lake(i,itype,classdata(j)%luse, & 
                              snowfall,temp,lakesurftemp,swrad,frozenstate,lakestate, & 
                              lakefreezeupday,lakebreakupday,lakearea(itype),epidepth,freezeuparea) 
        IF(modeloption(p_swtemperature)==1) CALL set_water_temperature(itype*2,i,riverstate,lakestate)
          
        !Abstraction of water from outlet lake
        CALL point_abstraction_from_outlet_lake(i,itype,qunitfactor,lakestate,nutrientloadflows(4))

        !Calculate substances processes in lakes (in flooded floodplain water, to be added here or in soilmodel4)
        CALL np_processes_in_lake(i,itype,lakearea(itype),    &
              (2.-incorr)*lakedatapar(lakedataparindex(i,itype),m_lddenitwl),genpar(m_hsatINwater),  &
              lakedatapar(lakedataparindex(i,itype),m_ldwprodn),   &
              lakedatapar(lakedataparindex(i,itype),m_ldwprodp),genpar(m_hsatTP),   &
              (2.-oncorr)*lakedatapar(lakedataparindex(i,itype),m_ldsedon),    &
              lakedatapar(lakedataparindex(i,itype),m_ldsedpp),genpar(m_limsedon),genpar(m_limsedpp),lakestate)
        CALL oc_processes_in_lake(i,itype,lakearea(itype),lakedatapar(lakedataparindex(i,itype),m_ldwprodc),    &
              genpar(m_hsatTP),genpar(m_limsedpp),lakedatapar(lakedataparindex(i,itype),m_ldsedoc),lakestate)
          
        !Calculate interflow between olake and floodplain
        IF(fpfrac>0.)THEN
          wlm3ol = lakestate%water(itype,i)*0.001 * lakearea(itype)   !Calculate water volume and average concentration in lake [m3]
          IF(numsubstances>0) cwlm3ol = lakestate%conc(:,itype,i)
          IF(conductN.OR.conductP.OR.conductC.OR.conductT)THEN
            CALL add_water(numsubstances,wlm3ol,cwlm3ol,lakestate%slowwater(itype,i)*0.001 * lakearea(itype),lakestate%concslow(:,itype,i))
          ENDIF
          IF(genpar(m_optonoff).LE.0)THEN   !Get current interflow parameter values
            interflowpar(1) = flooding(floodindex(i))%floll    !flol
            interflowpar(2) = flooding(floodindex(i))%flolp    !flfp
            interflowpar(3) = flooding(floodindex(i))%fymol
            interflowpar(4) = flooding(floodindex(i))%rcl2fp
            interflowpar(5) = flooding(floodindex(i))%rcfp2l
          ELSE
            interflowpar(1) = genpar(m_opt1)
            interflowpar(2) = genpar(m_opt2)
            interflowpar(3) = genpar(m_opt6)
            interflowpar(4) = genpar(m_opt5)
            interflowpar(5) = genpar(m_opt8)
          ENDIF
          CALL calculate_waterbody_floodplain_interflow(i,lakefparea,lakearea(itype),interflowpar,  &
                  miscstate%floodwater(2,i),miscstate%cfloodwater(:,2,i),wlm3ol, &
                  lakestate%conc(:,itype,i),flooddepth,flooddegree,qolflood)
          IF(outvarindex(324)>0) outvar(i,outvarindex(324)) = flooddepth
          IF(outvarindex(326)>0) outvar(i,outvarindex(326)) = flooddegree
          IF(conductwb)THEN
            IF(qolflood>0.) wbfpflows(w_oltofp,i) = qolflood
            IF(qolflood<0.) wbfpflows(w_fptool,i) = -qolflood
          ENDIF
          !Recalculate volume to lake water level [mm]:
          IF(conductN.OR.conductP.OR.conductC.OR.conductT)THEN
            lakestate%slowwater(itype,i) = wlm3ol/lakearea(itype)*1000.   !divide to slow and water according to slowlakeini?
            lakestate%water(itype,i) = 0.
            lakestate%concslow(:,itype,i) = cwlm3ol
            lakestate%conc(:,itype,i) = 0.
          ELSE
            lakestate%water(itype,i) = wlm3ol/lakearea(itype)*1000.
            IF(numsubstances>0) lakestate%conc(:,itype,i) = cwlm3ol   !onödig
          ENDIF
            
          !Update qin for calculation of outflow from average rating curve
          qin = MAX(qin - qolflood/seconds_per_timestep,0.)
        ELSE  !no floodplain
          IF(outvarindex(324)>0) outvar(i,outvarindex(324)) = missing_value
          IF(outvarindex(326)>0) outvar(i,outvarindex(326)) = 0.
          !outvar(i,324) = missing_value
          !outvar(i,326) = 0.
        ENDIF

        !Calculate and remove outflows from outlet lake
        lakewst(itype)=lakestate%water(itype,i)   !could use temp non-array variable here for these calculations
        IF(conductN.OR.conductP.OR.conductC.OR.conductT) lakewst(itype)=lakewst(itype)+lakestate%slowwater(itype,i)   !For NPC-simulation with divisioned lake
        CALL calculate_outflow_from_lake(i,itype,qin,lakearea(itype), &
                 lakewst(itype),qunitfactor,outflowm3s,outflowmm, &
                 outflow1,outflow2,maxProd,minFlow,lakestate)
        CALL check_outflow_from_lake(i,itype,lakearea(itype),lakewst(itype),qunitfactor,outflowm3s,outflowmm,lakestate)
        lakeoutflow(itype) = outflowm3s
        CALL remove_outflow_from_lake(i,itype,numsubstances,outflowmm,basin(i)%subid,basin(i)%lakedepth,epidepth,lakewst(itype),concout,lakestate)
        IF(conductload) Lpathway(:,13) = lakeoutflow(itype) * concout * seconds_per_timestep * 1.E-3  !Load at point M, outflow of subbasin (kg/timestep)
          
        !Flow between lake parts for divided lake
        CALL calculate_flow_within_lake(i,itype,basin(i)%subid,lakestate)

        !Lake water stage for print out 
        lakewst(itype) = lakestate%water(itype,i)
        IF(ALLOCATED(lakestate%slowwater)) lakewst(itype) = lakewst(itype) + lakestate%slowwater(itype,i)
        wstlakesim = lakewst(itype)          !saved for output
        wcomaver = (wstlakesim + olakewstold)/2.
          
        !Update lake outflow with waterstage observations
        outflowsim = lakeoutflow(itype)      !saved for qAR and output           
        IF(doupdate(i_war))THEN !move into lake>0
          CALL apply_warupd(i,lakearea(itype),olakewstold,lakewst(itype),lakeoutflow(itype),miscstate%updatestationsarcorr(i),wcomaver,lakestate)
        ENDIF
        IF(doupdate(i_wendupd))THEN
          IF(ALLOCATED(lakestate%slowwater).AND.lakedatapar(lakedataparindex(i,itype),m_lddeeplake)>0.)THEN
            CALL apply_wendupd(i,itype,lakearea(itype),lakewst(itype),lakestate,slowlakeini(itype,i))
          ELSEIF(ALLOCATED(lakestate%slowwater))THEN
            CALL apply_wendupd(i,itype,lakearea(itype),lakewst(itype),lakestate)
          ELSE
            CALL apply_wendupd(i,itype,lakearea(itype),lakewst(itype),lakestate)
          ENDIF
        ENDIF

        !>Update modelled subbasin outflow with observations or nutrient concentrations
        IF(doupdate(i_quseobs)) CALL apply_quseobs(i,lakeoutflow(itype))
        IF(doupdate(i_qar))     CALL apply_qarupd(i,outflowsim,lakeoutflow(itype),miscstate%updatestationsarcorr(i))
        CALL calculate_branched_flow_new(i,lakeoutflow(2),outflow1,outflow2,maxProd,minFlow,mainflow,branchflow)

        !Output
!        IF(conductload) Lbranch = branchflow * clakeoutflow(:,2) * seconds_per_timestep * 1.E-3  !Load in branch (kg/timestep)
        IF(conductload) Lbranch = branchflow * concout * seconds_per_timestep * 1.E-3  !Load in branch (kg/timestep)
        IF(conductwb)THEN
          wbflows(w_oltomb,i) = mainflow * seconds_per_timestep   !m3/timestep
          wbflows(w_oltob,i)  = branchflow * seconds_per_timestep
          wbflows(w_rural4,i) = nutrientloadflows(1)
          wbflows(w_pstomr,i) = nutrientloadflows(2)
          wbflows(w_abstmr,i) = nutrientloadflows(3)
          wbflows(w_abstol,i) = nutrientloadflows(4)
        ENDIF
 
      ELSE    !No lake, outflow = inflow
        lakeoutflow(itype) = dampq       !m3/s
        outflowsim = lakeoutflow(itype)      !saved for qAR and output
        qcinfli = qcinfli + dampq
        concout = dampc
        IF(conductload)THEN
          Lpathway(:,12) = lakeoutflow(itype) * concout * seconds_per_timestep * 1.E-3  !Load at point L flow in main river (kg/timestep)
          Lpathway(:,13) = Lpathway(:,12)
        ENDIF
        lakewst(itype) = missing_value   !not needed, lakewst used if lakearea(2)>0
        wstlakesim     = missing_value   !not needed, wstlakesim used if lakearea(2)>0
        lakesurftemp(itype) = missing_value
        !>Update modelled subbasin outflow with observations or nutrient concentrations
        IF(doupdate(i_quseobs)) CALL apply_quseobs(i,lakeoutflow(itype))
        IF(doupdate(i_qar))     CALL apply_qarupd(i,outflowsim,lakeoutflow(itype),miscstate%updatestationsarcorr(i))
         CALL calculate_branched_flow(i,lakeoutflow(2),mainflow,branchflow)
        !Output
        IF(conductload) Lbranch = branchflow * concout * seconds_per_timestep * 1.E-3  !Load in branch (kg/timestep)
        IF(conductwb)THEN
          wbflows(w_oltomb,i) = mainflow * seconds_per_timestep   !m3/timestep
          wbflows(w_oltob,i)  = branchflow * seconds_per_timestep
          wbflows(w_rural4,i) = nutrientloadflows(1)
          wbflows(w_pstomr,i) = nutrientloadflows(2)
          wbflows(w_abstmr,i) = nutrientloadflows(3)
          wbflows(w_abstol,i) = nutrientloadflows(4)
        ENDIF
      ENDIF   !a>0
      
      !Calculate theoretical precipitation and evaporation on water in case area_water = 0 (no ilake, olake, local river, or main river area in this subbasin)
      IF(area_water.LE.0.)THEN
        j=0
        IF(slc_lriver>0) j = slc_lriver
        IF(slc_ilake>0)  j = slc_ilake
        IF(slc_mriver>0) j = slc_mriver
        IF(slc_olake>0)  j = slc_olake
        IF(j>0)THEN
          !Forcing data for hypothetical water surface
          CALL calculate_class_atmospheric_forcing(i,j,radexti(i),  &
                temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
          CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,prec,temp,snowfall,rainfall,sffrac) 
          !Evaporation for hypothetical water surface
          CALL calculate_potential_evaporation(i,j,temp,epot,radexti(i),swrad,netrad,actvap,satvap,wind,epotsnow)
          epot = epot * basincevpcorr(i)
          IF(temp>landpar(m_ttmp,classdata(j)%luse) .AND. epot>0) THEN
            evapl = epot
          ELSE
            evapl = 0.
          ENDIF 
          !Accumulate outvar for precipitation and evaporation on water (probably excluding floodplain water)
          CALL calculate_class_outvar_add(o_clrp,i,1.,rainfall+snowfall)
          CALL calculate_class_outvar_add(o_clre,i,1.,evapl)
          area_water = area_water + 1.
        ELSE
          CALL calculate_class_outvar_add(o_clrp,i,1.,0.)
          CALL calculate_class_outvar_add(o_clre,i,1.,0.)
          area_water = area_water + 1.
        ENDIF
      ENDIF
           
      !Finalize subbasin contribution to down stream

      !>Update modelled subbasin outflow with nutrient concentrations
      IF(doupdate(i_tpcorr))  CALL apply_nutrientcorr(updatetpcorr(i),concout(i_sp),concout(i_pp))          
      IF(doupdate(i_tncorr))  CALL apply_nutrientcorr(updatetncorr(i),concout(i_in),concout(i_on))          
      clakeoutflow(:,itype) = concout
      
      !Accumulate flow to downstream subbasins
      IF(lakeoutflow(2)>0)THEN
        IF(mainflow>0 .AND. path(i)%main>0)THEN
          accinflow(path(i)%main) = accinflow(path(i)%main) + mainflow
          acccinflow(:,path(i)%main) = acccinflow(:,path(i)%main) + clakeoutflow(:,2) * mainflow
        ENDIF
        IF(branchflow>0)THEN
          IF(branchdata(branchindex(i))%branch>0)THEN
            accinflow(branchdata(branchindex(i))%branch) = accinflow(branchdata(branchindex(i))%branch) + branchflow
            acccinflow(:,branchdata(branchindex(i))%branch) = acccinflow(:,branchdata(branchindex(i))%branch) + clakeoutflow(:,2) * branchflow
          ENDIF
        ENDIF
      ENDIF

      !Subbasin output
      !----------------
      !Calculate subbasin mean over all classes
      IF(eacti>0.0 .AND. i_t1>0) CALL calculate_class_outvar_finish(o_cevapT1,i,eacti)  !T1 conc in evaporation
      interceptionloss(i) = pcorricep(i) + icpevapi   !used by several outvar
      
      !>Set subbasin output variables (outvar) for possible print out
      flow1000m3ts=lakeoutflow(2)*seconds_per_timestep*1.E-3  !m3/s -> 1000m3/ts
      CALL set_outvar_xobs(o_rrun,i)   !local recorded runoff, mm/ts
      CALL calculate_class_outvar_finish(o_crunT1,i,runoffi)  !T1 conc in runoff
      CALL calculate_class_outvar_finish(o_crunT2,i,runoffi)  !T2 in runoff
      CALL calculate_class_outvar_finish(o_crunIN,i,runoffi*1.E-3)  !IN conc in runoff [ug/L]
      CALL calculate_class_outvar_finish(o_crunON,i,runoffi*1.E-3)  !ON conc in runoff [ug/L]
      CALL calculate_class_outvar_finish(o_crunTN,i,runoffi*1.E-3)  !TN conc in runoff [ug/L]
      CALL calculate_class_outvar_finish(o_crunSP,i,runoffi*1.E-3)  !SP conc in runoff [ug/L]
      CALL calculate_class_outvar_finish(o_crunPP,i,runoffi*1.E-3)  !PP conc in runoff [ug/L]
      CALL calculate_class_outvar_finish(o_crunTP,i,runoffi*1.E-3)  !TP conc in runoff [ug/L]
      CALL calculate_class_outvar_finish(o_crunOC,i,runoffi)  !OC conc in runoff [mg/L]
      IF(i_t1>0) CALL set_outvar_xobs(o_cprecT1,i)
      IF(outvarindex(o_snow)>0) outvar(i,outvarindex(o_snow)) = snowi
      IF(outvarindex(o_evap)>0) outvar(i,outvarindex(o_evap)) = eacti
      IF(outvarindex(279)>0) outvar(i,outvarindex(279)) = interceptionloss(i) + eacti  !total evaporation (incl interception losses) [mm]
      CALL set_outvar_xobs(16,i)   !T1 in runoff
      CALL set_outvar_xobs(17,i)   !T2 in runoff
      CALL set_outvar_xobs(24,i)   !rec IN in outflow, ug/L
      CALL set_outvar_xobs_scaled(207,i,24,flow1000m3ts*1.E-3)  !rec IN load in outflow, kg/ts
      CALL set_outvar_xobs(25,i)   !rec ON in outflow, ug/L
      CALL set_outvar_xobs_scaled(208,i,25,flow1000m3ts*1.E-3)  !rec ON load in outflow, kg/ts
      CALL set_outvar_xobs(26,i)   !rec SP in outflow, ug/L
      CALL set_outvar_xobs_scaled(209,i,26,flow1000m3ts*1.E-3)  !rec SP load in outflow, kg/ts
      CALL set_outvar_xobs(27,i)   !rec PP in outflow, ug/L
      CALL set_outvar_xobs_scaled(210,i,27,flow1000m3ts*1.E-3)  !rec PP load in outflow, kg/ts
      CALL set_outvar_xobs(o_reTN,i)   !rec TN in outflow, ug/L
      CALL set_outvar_xobs_scaled(211,i,o_reTN,flow1000m3ts*1.E-3)  !rec TN load in outflow, kg/ts
      CALL set_outvar_xobs(o_reTP,i)   !rec TP in outflow, ug/L
      CALL set_outvar_xobs_scaled(212,i,o_reTP,flow1000m3ts*1.E-3)  !rec TP load in outflow, kg/ts

      IF(ALLOCATED(qobsi)) THEN
        IF(qobsi(i)/=-9999) THEN
          IF(outvarindex(30)>0) outvar(i,outvarindex(30)) = outflowsim-qobsi(i)  !Daily error in Q
        ENDIF
      ENDIF
      IF(outvarindex(31)>0) outvar(i,outvarindex(31)) = outflowsim               !Simulated outflow subbasin (cobc)
      CALL calculate_outvar_watertemperature(32,i,2,-9999.,riverstate,lakestate) !wtmp
      CALL calculate_outvar_watertemperature(283,i,2,0.,riverstate,lakestate)    !wtm0
      !Lake and River Ice, Snow, and T2 Water Temperature Model<<<
      IF(modeloption(p_lakeriverice)>0)THEN
        !Lake ice variables
        IF(slc_olake>0)THEN
          IF(outvarindex(148)>0) outvar(i,outvarindex(148)) = frozenstate%lakeice(2,i)        !'coli' 'comp olake ice depth'
          IF(outvarindex(150)>0) outvar(i,outvarindex(150)) = frozenstate%lakebice(2,i)       !'colb','comp olake blackice depth'
          IF(outvarindex(152)>0) outvar(i,outvarindex(152)) = frozenstate%lakesnowdepth(2,i)  !'cols','comp olake snow depth'
          CALL set_outvar_xobs(154,i)  !'roli','rec. olake ice depth'
          CALL set_outvar_xobs(156,i)  !'rolb','rec. olake blackice depth'
          CALL set_outvar_xobs(158,i)  !'rols','rec. olake snow depth'
        ENDIF
        IF(slc_ilake>0)THEN
          IF(outvarindex(149)>0) outvar(i,outvarindex(149)) = frozenstate%lakeice(1,i)               !'cili' 'comp ilake ice depth'
          IF(outvarindex(151)>0) outvar(i,outvarindex(151)) = frozenstate%lakebice(1,i)              !'cilb','comp ilake blackice depth'
          IF(outvarindex(153)>0) outvar(i,outvarindex(153)) = frozenstate%lakesnowdepth(1,i)         !'cils','comp ilake snow depth'
          CALL set_outvar_xobs(155,i)  !'rili','rec. ilake ice depth'
          CALL set_outvar_xobs(157,i)  !'rilb','rec. ilake blackice depth'
          CALL set_outvar_xobs(159,i)  !'rils','rec. ilake snow depth'
        ENDIF
        !River ice and snow depth variables
        IF(outvarindex(160)>0) outvar(i,outvarindex(160)) = frozenstate%riverice(2,i) !'cmri','comp main river ice depth'
        IF(outvarindex(161)>0) outvar(i,outvarindex(161)) = frozenstate%riverice(1,i) !'clri','comp local river ice depth' 
        IF(outvarindex(162)>0) outvar(i,outvarindex(162)) = frozenstate%riverbice(2,i) !'cmrb','comp main river blackice depth'
        IF(outvarindex(163)>0) outvar(i,outvarindex(163)) = frozenstate%riverbice(1,i) !'clrb','comp local river blackice depth'
        IF(outvarindex(164)>0) outvar(i,outvarindex(164)) = frozenstate%riversnowdepth(2,i) !'cmrs','comp main river snow depth' 
        IF(outvarindex(165)>0) outvar(i,outvarindex(165)) = frozenstate%riversnowdepth(1,i) !'clrs','comp local river snow depth'
        IF(outvarindex(248)>0) outvar(i,outvarindex(248)) = frozenstate%lakeicecov(2,i) !'coic','comp olake ice cover'
        IF(outvarindex(249)>0) outvar(i,outvarindex(249)) = frozenstate%lakeicecov(1,i) !'ciic','comp ilake ice cover' 
        IF(outvarindex(250)>0) outvar(i,outvarindex(250)) = frozenstate%rivericecov(2,i) !'cmic','comp main river ice cover'
        IF(outvarindex(251)>0) outvar(i,outvarindex(251)) = frozenstate%rivericecov(1,i) !'clic','comp local river ice cover'
        CALL set_outvar_xobs(166,i)  !'rmri','rec. main river ice depth'
        CALL set_outvar_xobs(167,i)  !'rlri','rec. local river ice depth'
        CALL set_outvar_xobs(168,i)  !'rmrb','rec. main river blackice depth'
        CALL set_outvar_xobs(169,i)  !'rlrb','rec. local river blackice depth'
        CALL set_outvar_xobs(170,i)  !'rmrs','rec. main river snow depth'
        CALL set_outvar_xobs(171,i)  !'rlrs','rec. local river snow depth'
      ENDIF
      IF(i_t2>0)THEN
      ! Lake temperature variables (surface, upper layer, lower layer, mean water temp)
        IF(slc_olake>0)THEN
          IF(outvarindex(172)>0) outvar(i,outvarindex(172)) = lakesurftemp(2)     !'olst','comp olake surface temp'  
          IF(outvarindex(173)>0) outvar(i,outvarindex(173)) = lakestate%uppertemp(2,i)  !'olut','comp olake upper temp'  
          IF(outvarindex(174)>0) outvar(i,outvarindex(174)) = lakestate%lowertemp(2,i)   !'ollt','comp olake lower temp' 
          IF(outvarindex(175)>0) outvar(i,outvarindex(175)) = lakestate%concslow(i_t2,2,i) !'olwt','comp olake mean  temp'
        ENDIF
        IF(slc_ilake>0)THEN
          IF(outvarindex(176)>0) outvar(i,outvarindex(176)) = lakesurftemp(1)     !'ilst','comp olake surface temp'   
          IF(outvarindex(177)>0) outvar(i,outvarindex(177)) = lakestate%concslow(i_t2,1,i) !'ilwt','comp olake mean  temp'  
        ENDIF
      ! River temperature variables (surface, mean, local and main)
        IF(outvarindex(178)>0) outvar(i,outvarindex(178)) = riversurftemp(1)      !lrst','comp local river surface temp'
        CALL get_rivertempvol(i,1,riverstate,meanrivertemp,totrivervol)
        IF(outvarindex(179)>0) outvar(i,outvarindex(179)) = meanrivertemp         !'lrwt','comp local river mean  temp'
        IF(outvarindex(180)>0) outvar(i,outvarindex(180)) = riversurftemp(2)       !'mrst','comp main  river surface temp'
        CALL get_rivertempvol(i,2,riverstate,meanrivertemp,totrivervol)
        IF(outvarindex(181)>0) outvar(i,outvarindex(181)) = meanrivertemp
      ENDIF
      !Additional variables from old water temp model
      IF(outvarindex(185)>0) outvar(i,outvarindex(185)) = riverstate%temp(2,i) !main river temperature
      IF(outvarindex(186)>0) outvar(i,outvarindex(186)) = riverstate%temp(1,i) !local river temp
      IF(slc_olake>0 .AND.outvarindex(187)>0) outvar(i,outvarindex(187)) = lakestate%temp(2,i) !olake temp
      IF(slc_ilake>0 .AND.outvarindex(188)>0) outvar(i,outvarindex(188)) = lakestate%temp(1,i) !ilake temp
      !Recorded Olake water surface temp
      CALL set_outvar_xobs(182,i)  !'rolt','rec. olake surface temp'
      CALL set_outvar_xobs(183,i)  !'rilt','rec. ilake surface temp'
      CALL set_outvar_xobs(184,i)  !'rmrt','rec. main river surface temp'
      !END of Ice and T2 temperature block<<<
      !Block with Remote Sensing Snow variables<<<
      CALL set_outvar_xobs(190,i)  !Recorded Fractional snow cover area(-)
      CALL set_outvar_xobs(192,i)  !Recorded Fractional snow cover area error (-)
      CALL set_outvar_xobs(193,i)  !Recorded Fractional snow cover multi (-)
      CALL set_outvar_xobs(194,i)  !Recorded Fractional snow cover multi error (-)
      !END snow remote sensing block<<<
      
      !Output variables for olake water stage and lake volume
      IF(outvarindex(33)>0.OR.outvarindex(34)>0.OR.outvarindex(51)>0.OR. &
         outvarindex(280)>0.OR.outvarindex(281)>0.OR.outvarindex(282)>0)THEN
        IF(slc_olake>0)THEN
          IF(lakearea(2)>0.)THEN
            CALL calculate_olake_waterstage(i,wstlakesim,lakearea(2),lakeareatemp2,wstlake,lakestate,w0ref)
            oldolakewst = wstlake + w0ref     !olake water stage before updating (m)
            IF(outvarindex(34)>0) outvar(i,outvarindex(34)) = oldolakewst
            CALL calculate_olake_waterstage(i,lakewst(2),lakearea(2),lakeareatemp2,wstlake,lakestate,w0ref)
            CALL calculate_regamp_adjusted_waterstage(i,lakearea(2),wstlake,wstlakeadj)
            IF(wstlakeadj/=missing_value.AND.outvarindex(280)>0) outvar(i,outvarindex(280)) = wstlakeadj    !regamp adjusted olake water stage (not in w-reference-system)
            IF(wstlakeadj/=missing_value.AND.outvarindex(51)>0) outvar(i,outvarindex(51)) = wstlakeadj + w0ref    !regamp adjusted olake water stage (in w-reference-system)
            IF(wcomaver/=missing_value)THEN
              CALL calculate_olake_waterstage(i,wcomaver,lakearea(2),lakeareatemp2,wstlake,lakestate,w0ref)
              CALL calculate_regamp_adjusted_waterstage(i,lakearea(2),wstlake,wstlakeadj)
              IF(outvarindex(282)>0) outvar(i,outvarindex(282)) = wstlakeadj + w0ref   !average olake water stage (in w-reference-system)
            ENDIF
            IF(xobsindex(o_rewstr,i)>0)THEN
              IF(xobsi(xobsindex(o_rewstr,i))/=missing_value)THEN
                IF(outvarindex(281)>0) outvar(i,outvarindex(281)) = xobsi(xobsindex(o_rewstr,i)) - w0ref       !recorded olake waterstage (cleaned from w0ref)
                IF(outvarindex(33)>0) outvar(i,outvarindex(33)) = oldolakewst - xobsi(xobsindex(o_rewstr,i)) !error in W
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      
      CALL set_outvar_xobs(35,i)  !DOC in outflow, mg/L
      CALL set_outvar_xobs_scaled(213,i,35,flow1000m3ts*1.E-3)  !rec DOC load in outflow, kg/ts
      IF(outvarindex(o_snowdepth)>0) outvar(i,outvarindex(o_snowdepth)) = snowdepthi    !snow depth (cm)
      IF(outvarindex(o_snowdens)>0.AND.snowdepthi>0.) outvar(i,outvarindex(o_snowdens))=snowi/snowdepthi*0.1  !snow density (g/cm3=cm/cm=1)
      CALL calculate_class_outvar_finish_scale(o_epot,i,1.)
      CALL set_outvar_xobs(o_reepot,i)  !Observed potential evaporation
      CALL set_outvar_xobs(42,i)  !Observed evaporation (mm)
      IF(outvarindex(o_cprc)>0) outvar(i,outvarindex(o_cprc)) = corrpreci      !Corrected precipitation (mm)
      IF(outvarindex(o_psim)>0) outvar(i,outvarindex(o_psim)) = corrpreci + interceptionloss(i)    !Simulated precipitation (mm)
      IF(conductC)THEN
        IF(outvarindex(46)>0) outvar(i,outvarindex(46)) = clakeoutflow(i_oc,2)  !DOC conc outflow olake, mg/L
      ENDIF
      CALL set_outvar_xobs(o_rewstr,i)  !recorded olake waterstage
      IF(outvarindex(o_cout)>0) outvar(i,outvarindex(o_cout)) = lakeoutflow(2)                          !computed outflow olake (m3/s)
      IF(outvarindex(o_rout)>0.AND.ALLOCATED(qobsi)) outvar(i,outvarindex(o_rout)) = qobsi(i)
      IF(i_in>0.AND.outvarindex(55)>0) outvar(i,outvarindex(55)) = clakeoutflow(i_in,2)*1000.   !comp conc IN outflow olake (ug/l)
      IF(i_on>0.AND.outvarindex(56)>0) outvar(i,outvarindex(56)) = clakeoutflow(i_on,2)*1000.   !comp conc ON outflow olake (ug/l)
      IF(i_sp>0.AND.outvarindex(57)>0) outvar(i,outvarindex(57)) = clakeoutflow(i_sp,2)*1000.   !comp conc SRP outflow olake (ug/l)
      IF(i_pp>0.AND.outvarindex(58)>0) outvar(i,outvarindex(58)) = clakeoutflow(i_pp,2)*1000.   !comp conc PP outflow olake (ug/l)
      CALL set_outvar_xobs(59,i)  !rec snow depth
      CALL set_outvar_xobs(60,i)  !rec soil frost
      CALL set_outvar_xobs(61,i)  !rec grw level
      CALL set_outvar_xobs(62,i)  !rec T1 conc in outflow
      IF(i_t1>0.AND.outvarindex(63)>0) outvar(i,outvarindex(63)) = clakeoutflow(i_t1,2)                        !comp conc T1 conc in outflow
      IF(i_t2>0.AND.outvarindex(64)>0) outvar(i,outvarindex(64)) = clakeoutflow(i_t2,2)                        !comp conc T2 conc in outflow
      IF(conductN)THEN
        IF(outvarindex(74)>0) outvar(i,outvarindex(74)) = csoili1(i_in)    !DIN pool soil1
        IF(outvarindex(75)>0) outvar(i,outvarindex(75)) = csoili2(i_in)    !DIN pool soil2
        IF(outvarindex(76)>0) outvar(i,outvarindex(76)) = csoili3(i_in)    !DIN pool soil3
      ENDIF
      IF(i_sp>0)THEN
        IF(outvarindex(91)>0) outvar(i,outvarindex(91)) = csoili1(i_sp)     !SRP pool soil1
        IF(outvarindex(92)>0) outvar(i,outvarindex(92)) = csoili2(i_sp)     !SRP pool soil2
        IF(outvarindex(93)>0) outvar(i,outvarindex(93)) = csoili3(i_sp)     !SRP pool soil3
      ENDIF
      IF(i_on>0)THEN
        IF(outvarindex(136)>0) outvar(i,outvarindex(136)) = csoili1(i_on)    !DON pool soil1
        IF(outvarindex(137)>0) outvar(i,outvarindex(137)) = csoili2(i_on)    !DON pool soil2
        IF(outvarindex(138)>0) outvar(i,outvarindex(138)) = csoili3(i_on)    !DON pool soil3
      ENDIF
      IF(i_t1>0)THEN
        IF(outvarindex(336)>0) outvar(i,outvarindex(336)) = csoili1(i_t1)    !soluable T1 pool soil1
        IF(outvarindex(337)>0) outvar(i,outvarindex(337)) = csoili2(i_t1)    !soluable T1 pool soil2
        IF(outvarindex(338)>0) outvar(i,outvarindex(338)) = csoili3(i_t1)    !soluable T1 pool soil3
      ENDIF
      IF(i_in>0.AND.i_on>0.AND.outvarindex(77)>0) outvar(i,outvarindex(77)) = (clakeoutflow(i_in,2) + clakeoutflow(i_on,2))*1000.  !TN conc in outflow olake, ug/L
      IF(i_sp>0.AND.i_pp>0.AND.outvarindex(78)>0) outvar(i,outvarindex(78)) = (clakeoutflow(i_sp,2) + clakeoutflow(i_pp,2))*1000.  !TP conc in outflow olake, ug/L
      IF(conductwb.OR.outvarindex(81)>0) CALL calculate_regional_groundwaterflow_to_outside_system(i,grwout(1,i),outofsystem)
      IF(outvarindex(81)>0) outvar(i,outvarindex(81)) = outofsystem/seconds_per_timestep !loss of water via groundwater flow from the model system, m3/s
      CALL set_outvar_xobs(o_cprecIN,i)  !rec concentration of IN in precipitation
      CALL set_outvar_xobs(o_cprecSP,i)  !rec concentration of SP in precipitation
      CALL set_outvar_xobs(o_reswe,i)    !rec snow water equivalent
      IF(firstoutstep) accdiff(i) = 0.  
      IF(ALLOCATED(qobsi))THEN
        IF(qobsi(i)>=0.)THEN
          accdiff(i) = accdiff(i) + (lakeoutflow(2)-qobsi(i))/upstreamarea(i)*seconds_per_timestep*1.E3  !accumulated volume error (mm)
          IF(outvarindex(100)>0) outvar(i,outvarindex(100))  = accdiff(i)
        ENDIF
      ENDIF
      IF(outvarindex(o_cloc)>0) outvar(i,outvarindex(o_cloc)) = lakeoutflow(1)
      IF(i_in>0.AND.outvarindex(102)>0) outvar(i,outvarindex(102)) = clakeoutflow(i_in,1)*1000.
      IF(i_on>0.AND.outvarindex(103)>0) outvar(i,outvarindex(103)) = clakeoutflow(i_on,1)*1000.
      IF(i_sp>0.AND.outvarindex(104)>0) outvar(i,outvarindex(104)) = clakeoutflow(i_sp,1)*1000.
      IF(i_pp>0.AND.outvarindex(105)>0) outvar(i,outvarindex(105)) = clakeoutflow(i_pp,1)*1000.
      IF(i_in>0.AND.i_on>0.AND.outvarindex(106)>0) outvar(i,outvarindex(106)) = (clakeoutflow(i_in,1) + clakeoutflow(i_on,1))*1000.
      IF(i_sp>0.AND.i_pp>0.AND.outvarindex(107)>0) outvar(i,outvarindex(107)) = (clakeoutflow(i_sp,1) + clakeoutflow(i_pp,1))*1000.
      IF(i_in>0.AND.i_on>0.AND.outvarindex(110)>0) outvar(i,outvarindex(110)) = (clakeoutflow(i_in,2) + clakeoutflow(i_on,2))*flow1000m3ts  !TN load (kg/timestep)
      IF(i_sp>0.AND.i_pp>0.AND.outvarindex(111)>0) outvar(i,outvarindex(111)) = (clakeoutflow(i_sp,2) + clakeoutflow(i_pp,2))*flow1000m3ts  !TP load (kg/timestep)
      IF(outvarindex(112)>0) outvar(i,outvarindex(112)) = qcinfli   !total inflow (minus evaporation) to olake (m3/s)
      CALL set_outvar_xobs(113,i)    !rec inflow (-"-) (m3/s)
      IF(conductwb.OR.outvarindex(114)>0.OR.outvarindex(204)>0) rivervolume(1) = riverstate%water(1,i) + (SUM(riverstate%qqueue(1:ttstep(1,i),1,i)) + riverstate%qqueue(ttstep(1,i)+1,1,i) * ttpart(1,i))
      IF(conductwb.OR.outvarindex(115)>0.OR.outvarindex(205)>0.OR.outvarindex(o_cmwl)>0) rivervolume(2) = riverstate%water(2,i) + (SUM(riverstate%qqueue(1:ttstep(2,i),2,i)) + riverstate%qqueue(ttstep(2,i)+1,2,i) * ttpart(2,i))
      IF(outvarindex(114)>0) outvar(i,outvarindex(114)) = rivervolume(1)
      IF(outvarindex(115)>0) outvar(i,outvarindex(115)) = rivervolume(2)
      IF(outvarindex(204)>0.AND.riverarea(1)>0) outvar(i,outvarindex(204)) = rivervolume(1)/riverarea(1)   !local river depth [m]
      IF(outvarindex(205)>0.AND.riverarea(2)>0) outvar(i,outvarindex(205)) = rivervolume(2)/riverarea(2)   !main river depth [m]
      IF(ALLOCATED(qobsi))THEN
        IF(qobsi(i)/=missing_value)THEN
          CALL set_outvar_xobs_scaled(119,i,o_reTN,qobsi(i)*seconds_per_timestep*1.E-6)  !rec TN load (kg/timestep)
          CALL set_outvar_xobs_scaled(120,i,o_reTP,qobsi(i)*seconds_per_timestep*1.E-6)  !rec TP load (kg/timestep)
        ENDIF
      ENDIF
      CALL calculate_class_outvar_finish_scale(o_ctmp,i,1.)
      IF(outvarindex(o_coum)>0) outvar(i,outvarindex(o_coum)) = mainflow       !Discharge in main channel (m3/s)
      CALL set_outvar_xobs(o_roum,i)    !Observed flow main channel
      IF(outvarindex(124)>0) outvar(i,outvarindex(124)) = ldremi         !Abstraction local dam for irrigation (m3)
      IF(outvarindex(125)>0) outvar(i,outvarindex(125)) = lrremi         !Abstraction local river for irrigation (m3)
      CALL calculate_class_outvar_finish_scale(o_applirr,i,basin(i)%area*1.E-3)
      IF(outvarindex(127)>0) outvar(i,outvarindex(127)) = gwremi         !Ground water removed for irrigation (m3)
      IF(outvarindex(128)>0) outvar(i,outvarindex(128)) = rsremi         !Abstraction regional surface water for irrigation (m3)
      IF(outvarindex(o_coub)>0) outvar(i,outvarindex(o_coub)) = branchflow     !Discharge in branched river (m3/s)

      !Simulated and Observed Main River Water Level (m.a.s.l)
      IF(outvarindex(o_cmwl)>0.AND.riverarea(2)>0)THEN
        outvar(i,outvarindex(o_cmwl)) = rivervolume(2)/riverarea(2) + basin(i)%mrwl0
      ENDIF
      CALL set_outvar_xobs(o_rmwl,i)
      CALL set_outvar_xobs(o_rmwd,i)
      CALL set_outvar_xobs(o_amwl,i)
   
      CALL set_outvar_xobs(o_roub,i)
      IF(i_oc>0.AND.outvarindex(135)>0) outvar(i,outvarindex(135)) = clakeoutflow(i_oc,2)*flow1000m3ts  !TOC load (kg/timestep)
      CALL calculate_class_outvar_finish(o_rainfall,i,1.)
      CALL calculate_class_outvar_finish(o_snowfall,i,1.)
      IF(outvarindex(196)>0) outvar(i,outvarindex(196)) = (soilwateri(1)+soilwateri(2))*divasum(1)
      IF(outvarindex(o_sml9)>0) outvar(i,outvarindex(o_sml9)) = soilwateri(1)*divasum(1)
      IF(outvarindex(o_sml1)>0) outvar(i,outvarindex(o_sml1)) = soilwateri(1)*divasum(1) - standsoili
      IF(outvarindex(o_sml2)>0) outvar(i,outvarindex(o_sml2)) = soilwateri(2)*divasum(2)
      IF(outvarindex(o_sml3)>0) outvar(i,outvarindex(o_sml3)) = soilwateri(3)*divasum(3)
      IF(outvarindex(200)>0) outvar(i,outvarindex(200)) = standsoili
      IF(outvarindex(201)>0) outvar(i,outvarindex(201)) = (soilwateri(1)+soilwateri(2))*divasum(1) - standsoili
      IF(outvarindex(202)>0) outvar(i,outvarindex(202)) = (soilwateri(1)+soilwateri(2)+soilwateri(3))*divasum(1) - standsoili
      IF(i_oc>0.AND.outvarindex(203)>0) outvar(i,outvarindex(203)) = clakeoutflow(i_oc,1)    !orgC in flow from local river
      DO k=0,9
        CALL set_outvar_xobsmean(o_xobsm+k,i,k+1)    !recorded meantype variable
        CALL set_outvar_xobsstate(o_xobss+k,i,k+1)   !recorded sumtype variable
      ENDDO
      CALL calculate_class_outvar_finish(o_snowmelt,i,asum(1))
      IF(outvarindex(o_specificq)>0) outvar(i,outvarindex(o_specificq)) = 1000.*(lakeoutflow(2)*seconds_per_timestep) / upstreamarea(i) !Specific discharge (upstream runoff)
          
      !Set subbasin load output variables for possible printout for Source Apportionment program
      IF(conductload) THEN
        outvar_classload(1:nclass,1:2,:,i)  = Latmdep(1:nclass,1:2,:)
        outvar_classload(1:nclass,3:4,:,i)  = Lcultiv(1:nclass,1:2,:)
        outvar_classload(1:nclass,5,:,i)    = Lrurala(1:nclass,:)
        outvar_classload(1:nclass,6,:,i)    = Lgrwsoil(1:nclass,:)
        outvar_classload(1:nclass,7,:,i)    = Lirrsoil(1:nclass,:)
        outvar_classload(1:nclass,8,:,i)    = Lstream(1:nclass,:)
        outvar_basinload(:,1,i)             = Lruralb(:)
        outvar_basinload(:,2:4,i)           = Lpoints(:,1:3)     !Note: Load-files and SourceApp need to change if max_pstype is increased
        outvar_basinload(:,5,i)             = Lgrwmr(:)
        outvar_basinload(:,6,i)             = Lgrwol(:)
        outvar_basinload(:,7:19,i)          = Lpathway(:,1:13)
        outvar_basinload(:,20,i)            = Lbranch(:)
      ENDIF
      
      !Set water balance flows and stores for print out (current subbasin)
      IF(conductwb)THEN
        IF(modeloption(p_deepgroundwater)==1) wbflows(w_rgrwtoos,i) = outofsystem   !moved up into DO-loop
        IF(conductflood)THEN
          IF(floodindex(i)>0)THEN
            IF(flooding(floodindex(i))%fpfol>0.)THEN
              wbstores(w_lakeplain,i) = miscstate%floodwater(2,i)
            ENDIF
            IF(flooding(floodindex(i))%fpfmr>0.)THEN
              wbstores(w_riverplain,i) = miscstate%floodwater(1,i)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
        
      !FSUS snow outputs, forest open land
      CALL set_outvar_xobs(257,i) !'S105','fsusnow fsc surr. open'
      CALL set_outvar_xobs(258,i) !'S106','fsusnow fsc course open'
      CALL set_outvar_xobs(259,i) !'S108','fsusnow mean depth open'
      CALL set_outvar_xobs(260,i) !'S111','fsusnow mean density open'
      CALL set_outvar_xobs(261,i) !'S114','fsusnow snow water eq. open'
      CALL set_outvar_xobs(262,i) !'S205','fsusnow fsc surr. forest'
      CALL set_outvar_xobs(263,i) !'S206','fsusnow fsc course forest'
      CALL set_outvar_xobs(264,i) !'S208','fsusnow mean depth forest'
      CALL set_outvar_xobs(265,i) !'S211','fsusnow mean density forest'
      CALL set_outvar_xobs(266,i) !'S214','fsusnow snow water eq. forest'
      IF(area_open>0)THEN
        IF(outvarindex(267)>0) outvar(i,outvarindex(267)) = snowvari(2,1)     !'C106','comp. fsc course open'
        IF(outvarindex(268)>0) outvar(i,outvarindex(268)) = snowvari(1,1)     !'C108','comp. mean depth open'
        IF(outvarindex(269)>0) outvar(i,outvarindex(269)) = snowvari(3,1)     !'C111','comp. mean density open'
        IF(outvarindex(270)>0) outvar(i,outvarindex(270)) = snowvari(4,1)     !'C114','comp. snow water eq. open'
      ENDIF
      IF(area_forest>0)THEN
        IF(outvarindex(271)>0) outvar(i,outvarindex(271)) = snowvari(2,2)   !'C206','comp. fsc course forest'
        IF(outvarindex(272)>0) outvar(i,outvarindex(272)) = snowvari(1,2)   !'C208','comp. mean depth forest'
        IF(outvarindex(273)>0) outvar(i,outvarindex(273)) = snowvari(3,2)   !'C211','comp. mean density forest'
        IF(outvarindex(274)>0) outvar(i,outvarindex(274)) = snowvari(4,2)   !'C214','comp. snow water eq. forest'
      ENDIF
      !Soil layer load
      IF(i_in>0.AND.i_on>0)THEN
        IF(outvarindex(305)>0) outvar(i,outvarindex(305)) = soillgrossload(1,i_in)*basin(i)%area*1.E-6    !mm*mg/L*m2/1000000=kg
        IF(outvarindex(306)>0) outvar(i,outvarindex(306)) = soillnetload(1,i_in)*basin(i)%area*1.E-6
        IF(outvarindex(307)>0) outvar(i,outvarindex(307)) = soillgrossload(1,i_on)*basin(i)%area*1.E-6
        IF(outvarindex(308)>0) outvar(i,outvarindex(308)) = soillnetload(1,i_on)*basin(i)%area*1.E-6
        IF(outvarindex(309)>0) outvar(i,outvarindex(309)) = (soillgrossload(1,i_in)+soillgrossload(1,i_on))*basin(i)%area*1.E-6
        IF(outvarindex(310)>0) outvar(i,outvarindex(310)) = (soillnetload(1,i_in)+soillnetload(1,i_on))*basin(i)%area*1.E-6
        IF(outvarindex(291)>0) outvar(i,outvarindex(291)) = soillgrossload(2,i_in)*basin(i)%area*1.E-6
        IF(outvarindex(292)>0) outvar(i,outvarindex(292)) = soillnetload(2,i_in)*basin(i)%area*1.E-6
        IF(outvarindex(293)>0) outvar(i,outvarindex(293)) = soillgrossload(2,i_on)*basin(i)%area*1.E-6
        IF(outvarindex(294)>0) outvar(i,outvarindex(294)) = soillnetload(2,i_on)*basin(i)%area*1.E-6
        IF(outvarindex(295)>0) outvar(i,outvarindex(295)) = (soillgrossload(2,i_in)+soillgrossload(2,i_on))*basin(i)%area*1.E-6
        IF(outvarindex(296)>0) outvar(i,outvarindex(296)) = (soillnetload(2,i_in)+soillnetload(2,i_on))*basin(i)%area*1.E-6
        CALL calculate_class_outvar_finish_scale(o_soilden3,i,basin(i)%area*1.E-6)    !denitrification soil layer 3 (kg)
        CALL calculate_class_outvar_finish_scale(o_soildenrz,i,basin(i)%area*1.E-6)   !denitrification soil layer 1 and 2 (kg)
        IF(outvarindex(311)>0) outvar(i,outvarindex(311)) = soillgrossload(3,i_in)*basin(i)%area*1.E-6
        IF(outvarindex(312)>0) outvar(i,outvarindex(312)) = soillnetload(3,i_in)*basin(i)%area*1.E-6
        IF(outvarindex(313)>0) outvar(i,outvarindex(313)) = soillgrossload(3,i_on)*basin(i)%area*1.E-6
        IF(outvarindex(314)>0) outvar(i,outvarindex(314)) = soillnetload(3,i_on)*basin(i)%area*1.E-6
        IF(outvarindex(315)>0) outvar(i,outvarindex(315)) = (soillgrossload(3,i_in)+soillgrossload(3,i_on))*basin(i)%area*1.E-6
        IF(outvarindex(316)>0) outvar(i,outvarindex(316)) = (soillnetload(3,i_in)+soillnetload(3,i_on))*basin(i)%area*1.E-6
      ENDIF
      IF(i_sp>0.AND.i_pp>0)THEN
        IF(outvarindex(285)>0) outvar(i,outvarindex(285)) = soillgrossload(1,i_sp)*basin(i)%area*1.E-6    
        IF(outvarindex(286)>0) outvar(i,outvarindex(286)) = soillnetload(1,i_sp)*basin(i)%area*1.E-6    
        IF(outvarindex(287)>0) outvar(i,outvarindex(287)) = soillgrossload(1,i_pp)*basin(i)%area*1.E-6    
        IF(outvarindex(288)>0) outvar(i,outvarindex(288)) = soillnetload(1,i_pp)*basin(i)%area*1.E-6    
        IF(outvarindex(289)>0) outvar(i,outvarindex(289)) = (soillgrossload(1,i_sp)+soillgrossload(1,i_pp))*basin(i)%area*1.E-6    
        IF(outvarindex(290)>0) outvar(i,outvarindex(290)) = (soillnetload(1,i_sp)+soillnetload(1,i_pp))*basin(i)%area*1.E-6    
        IF(outvarindex(297)>0) outvar(i,outvarindex(297)) = soillgrossload(2,i_sp)*basin(i)%area*1.E-6    
        IF(outvarindex(298)>0) outvar(i,outvarindex(298)) = soillnetload(2,i_sp)*basin(i)%area*1.E-6    
        IF(outvarindex(299)>0) outvar(i,outvarindex(299)) = soillgrossload(2,i_pp)*basin(i)%area*1.E-6    
        IF(outvarindex(300)>0) outvar(i,outvarindex(300)) = soillnetload(2,i_pp)*basin(i)%area*1.E-6    
        IF(outvarindex(301)>0) outvar(i,outvarindex(301)) = (soillgrossload(2,i_sp)+soillgrossload(2,i_pp))*basin(i)%area*1.E-6    
        IF(outvarindex(302)>0) outvar(i,outvarindex(302)) = (soillnetload(2,i_sp)+soillnetload(2,i_pp))*basin(i)%area*1.E-6    
        IF(outvarindex(317)>0) outvar(i,outvarindex(317)) = soillgrossload(3,i_sp)*basin(i)%area*1.E-6    
        IF(outvarindex(318)>0) outvar(i,outvarindex(318)) = soillnetload(3,i_sp)*basin(i)%area*1.E-6    
        IF(outvarindex(319)>0) outvar(i,outvarindex(319)) = soillgrossload(3,i_pp)*basin(i)%area*1.E-6    
        IF(outvarindex(320)>0) outvar(i,outvarindex(320)) = soillnetload(3,i_pp)*basin(i)%area*1.E-6    
        IF(outvarindex(321)>0) outvar(i,outvarindex(321)) = (soillgrossload(3,i_sp)+soillgrossload(3,i_pp))*basin(i)%area*1.E-6    
        IF(outvarindex(322)>0) outvar(i,outvarindex(322)) = (soillnetload(3,i_sp)+soillnetload(3,i_pp))*basin(i)%area*1.E-6    
      ENDIF
      !Tracer T1 output
      IF(i_t1>0)THEN
        IF(outvarindex(339)>0) outvar(i,outvarindex(339)) = riverstate%T1sed(2,i)
        IF(outvarindex(340)>0) outvar(i,outvarindex(340)) = riverstate%T1sed(1,i)
        IF(outvarindex(342)>0) outvar(i,outvarindex(342)) = clakeoutflow(i_t1,1)    !T1 in flow from local river
      ENDIF
      
      !OBSERVATION OPERATORS
      !Additional outvars for data assimilation giving 'model equivalents' for various observations (some only calculated as regional outvars)
      
      !CLRF=local runoff to the main river (mm/tstep) = cloc[m3/s] normalized by local basin area excluding olake and main river area
      IF(outvarindex(o_clrf)>0)THEN
        a = 1.                                                             !fraction of basin area used for normalization, initially set to 1.
        IF(slc_olake>0) a = a - classbasin(i,slc_olake)%part               !remove olake fraction
        IF(slc_mriver>0) a = a - classbasin(i,slc_mriver)%part             !remove main river fraction
        IF(a.GT.0)THEN
          outvar(i,outvarindex(o_clrf)) = 1000.*(lakeoutflow(1)*seconds_per_timestep)/(a * basin(i)%area) !transform m3/s to mm/tstep
        ELSE
          outvar(i,outvarindex(o_clrf)) = 0. !if there is no land, ilake, or local river in this basin, then there is no local runoff either          
        ENDIF
      ENDIF
      
      !CLRP and CLRE
      CALL calculate_class_outvar_finish(o_clrp,i,area_water)
      CALL calculate_class_outvar_finish(o_clre,i,area_water)
   
    !>End main subbasin-loop   
  ENDDO subbasinloop
     
    IF(outvarindex(122)>0) outvar(:,outvarindex(122)) = irrevap     !Irrigation losses (m3)
    IF(outvarindex(o_icloss)>0) outvar(:,outvarindex(o_icloss)) = interceptionloss  !interception losses [mm] (preccorr and pcluse)

    !Calculate lake volume and set output variables
    basinlakevol = missing_value
    CALL calculate_lake_volume_output(nsub,lakestate,basinlakevol,ilakevol,olakevol,ilakevolmiss,olakevolmiss)
    IF(outvarindex(144)>0) outvar(:,outvarindex(144)) = basinlakevol  !Mm3
    IF(outvarindex(145)>0) outvar(:,outvarindex(145)) = ilakevolmiss  !Mm3
    IF(outvarindex(146)>0) outvar(:,outvarindex(146)) = olakevolmiss  !Mm3
    IF(conductwb)THEN
      wbstores(w_ilake,:)  = ilakevol  !m3
      wbstores(w_olake,:)  = olakevol  !m3
    ENDIF
    
    !Calculate aquifer delay and outflow
    IF(modeloption(p_deepgroundwater)==2)THEN
      CALL calculate_aquifer(aquiferstate,aquiferoutflow,aqremflow,aqirrloss)
      IF(outvarindex(240)>0) outvar(:,outvarindex(240)) = aqremflow
    ENDIF
    IF(naquifers>0)THEN
      IF(outvarindex(242)>0) outvar(:,outvarindex(242)) = calculate_aquifer_waterlevel(nsub,naquifers,aquiferstate)
    ENDIF
     
    !Set water balance flows and stores for print out (all subbasins)
    IF(conductwb)THEN
      wbflows(w_rgrwto1,:)  = horizontalflows2(1,:)     !regional groundwater flow from this subbasin's groundwater reservoir
      wbflows(w_rgrwto2,:)  = horizontalflows2(2,:)
      wbflows(w_rgrwto3,:)  = horizontalflows2(3,:)
      wbirrflows(w_wdregol,:)  = regionalirrflows(1,:)  !regional sources to this subbasin
      wbirrflows(w_wdregmr,:)  = regionalirrflows(2,:)  !regional sources to this subbasin
      wbirrflows(w_evapregol,:) = regionalirrflows(3,:)
      wbirrflows(w_evapregmr,:) = regionalirrflows(4,:)
      wbstores(w_iriver,:) = rivervolume(1)
      wbstores(w_mriver,:) = rivervolume(2)
      IF(doirrigation)THEN
        wbstores(w_irrcanal,:) = 0.
        DO i = 1,nsub
          DO j = 1,nclass
            wbstores(w_irrcanal,i) = wbstores(w_irrcanal,i) + miscstate%nextirrigation(j,i)*classbasin(i,j)%part
          ENDDO
        ENDDO
        wbstores(w_irrcanal,:) = wbstores(w_irrcanal,:) * basin(:)%area *1.E-3  !m3
      ENDIF  
      IF(naquifers>0)THEN
        CALL calculate_delayed_water(aquiferstate,naquifers,delayedwater)
        wbstores(w_aquifer,1:naquifers) = aquiferstate%water + aquiferstate%nextoutflow + delayedwater
      ENDIF
    ENDIF
     
    !Print out water balance flow 
    IF(conductwb) CALL print_waterbalance_timestep(nsub,naquifers,currentdate)
     
  END SUBROUTINE model

  !Private subroutines, may be moved, for updating   
  !----------------------------------------------------------------
  
  !>Update outflow of subbasin to observed value
  !-------------------------------------------------
  SUBROUTINE apply_quseobs(i,simflow)

    USE MODVAR, ONLY : qobsi,   &
                       missing_value,  &
                       updatestations
     
    !Argument declaration
    INTEGER, INTENT(IN) :: i         !<index of current subbasin
    REAL, INTENT(INOUT) :: simflow   !<simulated outflow of subbasin !!(lakeoutflow)
     
    IF(ALLOCATED(updatestations))THEN 
      IF(updatestations(i))THEN
        IF(ALLOCATED(qobsi))THEN
          IF(qobsi(i)/=missing_value)  simflow = qobsi(i)          
        ENDIF
      ENDIF
    ENDIF
     
  END SUBROUTINE apply_quseobs
   
  !>Update outflow of subbasin from observed value with AR method
  !----------------------------------------------------------------
  SUBROUTINE apply_qarupd(i,simflow,corrFlow,arcorr)

    USE MODVAR, ONLY : qobsi,  &
                       missing_value,       &
                       updatestationsqar,   &
                       updatestationsarfact
     
    !Argument declaration
    INTEGER, INTENT(IN)            :: i           !<index of current subbasin
    REAL, INTENT(IN)               :: simflow     !<simulated outflow of subbasin
    REAL, INTENT(INOUT)            :: corrFlow    !<updated outflow (lakeoutflow)
    REAL, INTENT(INOUT)            :: arcorr      !<current AR-error (state-variable)
     
    IF(ALLOCATED(updatestationsqar)) THEN 
      IF(updatestationsqar(i)) THEN
        IF(ALLOCATED(qobsi)) THEN
          IF(qobsi(i)/=missing_value) THEN  !calculates the error
            arcorr =  simflow - qobsi(i)
          ELSE  !no observation, using AR
            arcorr = arcorr * updatestationsarfact(i) !Updating AR-correction  
            corrFlow = simflow - arcorr
            IF(corrFlow<0.) corrFlow=0.
          ENDIF
        ELSE  !no observation, using AR
          arcorr = arcorr * updatestationsarfact(i) !Updating AR-correction  
          corrFlow = simflow - arcorr
          IF(corrFlow<0.) corrFlow=0.
        ENDIF
      ENDIF
    ENDIF
  
  END SUBROUTINE apply_qarupd
   
  !>Update outflow of subbasin from observed waterstage value with AR method
  !-------------------------------------------------------------------------
  SUBROUTINE apply_warupd(i,lakeareain,wstold,corrWst,corrFlow,arcorr,corrWstAve,lakestate)

    USE MODVAR, ONLY : xobsi,  &
                       xobsindex, &
                       missing_value,       &
                       wobsvar, &
                       updatestationswar,   &
                       updatestationsarfact
    USE SURFACEWATER_PROCESSES, ONLY : calculate_olake_waterstage, &
                                       calculate_flow_from_lake_waterstage
  
    !Argument declaration
    INTEGER, INTENT(IN)            :: i           !<index of current subbasin
    REAL, INTENT(IN)               :: lakeareain  !<olake area of subbasin (m2)
    REAL, INTENT(IN)               :: wstold      !<simulated lake water end of last time step (mm)
    REAL, INTENT(INOUT)            :: corrWst     !<IN: simulated lake water, OUT: updated lake water (lakewst - for print out only) (mm)
    REAL, INTENT(INOUT)            :: corrFlow    !<updated outflow (lakeoutflow, m3/s)
    REAL, INTENT(INOUT)            :: arcorr      !<current AR-error (state-variable, mm)
    REAL, INTENT(INOUT)            :: corrWstAve  !<updated average lake water (mm)
    TYPE(lakestatetype),INTENT(IN) :: lakestate   !<Lake state

    !Local variables     
    REAL corrWstold     !updated water stage last time step (mm)
    REAL errorw   !current error in wst in mm
    REAL lakearea !lake area of outlet lake (m2) (whole lake of last lakebasin)
    REAL qoutold  !outflow at old waterstage
    REAL qoutnew  !outflow at new waterstage
    REAL wstobs   !observed lake waterstage wstr from Xobs in w-ref system (m)
    REAL wstm     !lake water stage in local system (m)
    REAL w0ref    !waterstage reference level (m)
    
    !>\b Algorithm \n
    lakearea = lakeareain

    IF(ALLOCATED(updatestationswar)) THEN 
      IF(updatestationswar(i)) THEN
        !>Calculate updated waterstage last timestep
        corrWstold = wstold - arcorr
        !>Calculate error and new arcorr-factor
        IF(ALLOCATED(xobsi))THEN
          IF(xobsindex(wobsvar,i)>0)THEN
            wstobs = xobsi(xobsindex(wobsvar,i))  !Get current wst observation, m in wref-system
          ELSE
            wstobs = missing_value
          ENDIF
          IF(wstobs/=missing_value)THEN
            !Calculate current error of water stage
            CALL calculate_olake_waterstage(i,corrWst,lakeareain,lakearea,wstm,lakestate,w0ref)
            errorw = (wstobs-w0ref-wstm)*1000.         !mm
            IF(ALLOCATED(lakebasinindex))THEN   !For last lake basin, transform change to local deltaw
              IF(lakebasinindex(i)>0)THEN  
                IF(lakebasin(lakebasinindex(i))%last)THEN
                  errorw = errorw * lakearea / lakeareain
                ENDIF
              ENDIF
            ENDIF
            arcorr = - errorw
          ELSE  !no observation, use wAR to update wst and q
            arcorr = arcorr * updatestationsarfact(i) !Updating AR-correction  
          ENDIF
        ELSE  !no observation, use wAR to update wst and q
          arcorr = arcorr * updatestationsarfact(i) !Updating AR-correction  
        ENDIF
        !>Apply AR correction to waterstage and discharge
        corrWst = corrWst - arcorr
        CALL calculate_flow_from_lake_waterstage(i,2,lakeareain,corrWstold,qoutold,lakestate)
        CALL calculate_flow_from_lake_waterstage(i,2,lakeareain,corrWst,qoutnew,lakestate)
        corrFlow = 0.5*(qoutold+qoutnew)
        IF(corrFlow<0.) corrFlow=0.
        corrWstAve = (corrWstold + corrWst)/2.
      ENDIF
    ENDIF
  
  END SUBROUTINE apply_warupd
   
  !>Update outlet lake water stage to observed value
  !-----------------------------------------------------------
  SUBROUTINE apply_wendupd(i,isystem,lakeareain,wst,lakestate,slwimax)

    USE MODVAR, ONLY : xobsi,   &
                       xobsindex,  &
                       missing_value,   &
                       wendupdstations, &
                       wobsvar
    USE SURFACEWATER_PROCESSES, ONLY : calculate_olake_waterstage
     
    !Argument declaration
    INTEGER, INTENT(IN)            :: i            !<index of subbasin
    INTEGER, INTENT(IN)            :: isystem      !<local or outlet lake
    REAL, INTENT(IN)               :: lakeareain   !<lake area of subbasin (m2)
    REAL, INTENT(INOUT)            :: wst          !<water in lake (=lwi+slwi) (mm) to be written output
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate !<Lake state
    REAL,OPTIONAL, INTENT(IN)      :: slwimax      !<target value for slowlake (mm)
     
    !Local variables
    REAL fraction
    REAL wstobs   !wstr for Xobs in w-ref system
    REAL wstm     !lake water stage in w-ref system
    REAL w0ref    !waterstage reference level (m)
    REAL deltaw   !change in calculated water stage due to updating (mm)
    REAL lakearea !lake area of outlet lake (m2)
     
    IF(wendupdstations(i)) THEN
      wstobs = xobsi(xobsindex(wobsvar,i))   !m
      IF(wstobs/=missing_value)THEN
        !Calculate change of water stage and apply to output variable
        CALL calculate_olake_waterstage(i,wst,lakeareain,lakearea,wstm,lakestate,w0ref)
        deltaw = (wstobs-w0ref-wstm)*1000.         !mm
        IF(ALLOCATED(lakebasinindex))THEN   !For last lake basin, transform change to local deltaw
          IF(lakebasinindex(i)>0)THEN  
            IF(lakebasin(lakebasinindex(i))%last)THEN
              deltaw = deltaw * lakearea / lakeareain
            ENDIF
          ENDIF
        ENDIF
        wst = wst + deltaw      !wst = wstobs (mm)
        !Apply updating on state variables for lake water
        IF(PRESENT(slwimax))THEN
          fraction = lakestate%water(isystem,i)/(lakestate%water(isystem,i)+lakestate%slowwater(isystem,i))
          lakestate%slowwater(isystem,i)=(1.-fraction)*wst
          lakestate%water(isystem,i)=fraction*wst
          IF(lakestate%slowwater(isystem,i)>slwimax)THEN
            lakestate%water(isystem,i) = lakestate%water(isystem,i) + (lakestate%slowwater(isystem,i)-slwimax)
            lakestate%slowwater(isystem,i)=slwimax
          ENDIF
        ELSEIF(ALLOCATED(lakestate%slowwater))THEN
          lakestate%slowwater(isystem,i) = wst
        ELSE
          lakestate%water(isystem,i) = wst
        ENDIF
      ENDIF
    ENDIF
     
  END SUBROUTINE apply_wendupd
   
  !>Update concentration out of subbasin to fraction of modelled value
  !----------------------------------------------------------------------
  SUBROUTINE apply_nutrientcorr(correction,conc1,conc2)

    !Argument declaration
    REAL, INTENT(IN)    :: correction   !<update correction value
    REAL, INTENT(INOUT) :: conc1        !<simulated concentration (SP,IN) of outflow of subbasin (clakeoutflow)
    REAL, INTENT(INOUT) :: conc2        !<simulated concentration (PP,ON) of outflow of subbasin (clakeoutflow)
     
    REAL  factor    !current correction factor
     
    IF(correction/=0.)THEN
      factor = 1. + correction
      conc1  = conc1*factor
      conc2  = conc2*factor
    ENDIF
     
  END SUBROUTINE apply_nutrientcorr

  !>Reads files with model specific input
  !>For HYPE it is files with different observations (xom0..xom9,xos0..xos9)
  !----------------------------------------------------------------------
  SUBROUTINE load_modeldefined_input(dir,nsmax,ns,indexarray,bdate,edate,loadxoms,loadregs,status)

    USE HYPE_INDATA, ONLY : load_xoms_files,  &
                            load_data_for_regression_parameter_estimate
    USE LibDate
    
    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<File directory (modeldir)
    INTEGER, INTENT(IN)  :: nsmax         !<Number of subbasins, basemodel
    INTEGER, INTENT(IN)  :: ns            !<Number of subbasins, submodel
    INTEGER, INTENT(IN) :: indexarray(ns) !<index for basemodel
    TYPE(DateType), INTENT(IN)  :: bdate  !<Begin simulation date
    TYPE(DateType), INTENT(IN)  :: edate  !<End simulation date
    LOGICAL, INTENT(IN) :: loadxoms       !<flag for using xoms-files
    LOGICAL, INTENT(IN) :: loadregs       !<flag for using regression estimation-files
    INTEGER, INTENT(OUT) :: status        !<Status of subroutine

    status = 0
    IF(loadxoms) CALL load_xoms_files(dir,ns,bdate,edate,status)
    IF(status/=0) RETURN
    IF(loadregs) CALL load_data_for_regression_parameter_estimate(dir,nsmax,ns,indexarray,status)
    IF(status/=0) RETURN
    
  END SUBROUTINE load_modeldefined_input

  !>Reads files with model specific input
  !>For HYPE it is files with different observations (xom0..xom9,xos0..xos9)
  !----------------------------------------------------------------------
  SUBROUTINE reload_modeldefined_observations(dir,status)

    USE HYPE_INDATA, ONLY : close_xoms_files, &
                            reload_xoms_files
    
    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<File directory (modeldir)
    INTEGER, INTENT(OUT) :: status        !<Status of subroutine

    status = 0
    IF(conductxoms) CALL close_xoms_files()
    IF(conductxoms) CALL reload_xoms_files(dir,status)
    IF(status/=0) RETURN
    
  END SUBROUTINE reload_modeldefined_observations

  !>Opens files for printout
  !--------------------------------------------------------------------
  SUBROUTINE open_modeldefined_outputfiles(dir,n,na,iens,runens)

    USE HYPE_WATERBALANCE, ONLY : prepare_waterbalance_files
    
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir !<Result file directory
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    LOGICAL, INTENT(IN) :: runens       !<Flag for ensemble simulation
    
    IF(conductwb) CALL prepare_waterbalance_files(dir,n,na,basin%subid)  !These cannot be printed for several emsembles.

  END SUBROUTINE open_modeldefined_outputfiles

  !>Close files for printout
  !--------------------------------------------------------------------
  SUBROUTINE close_modeldefined_outputfiles(n,na,iens)

    USE HYPE_WATERBALANCE, ONLY : close_waterbalance_files
    
    !Argument declarations
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    
    IF(conductwb) CALL close_waterbalance_files(na)

  END SUBROUTINE close_modeldefined_outputfiles
    
  !>Calculate output variables for lake volumes
  !--------------------------------------------------------------------
  SUBROUTINE calculate_lake_volume_output(n,lakestate,wlakevol,ilakevol,olakevol,Milakevol,Molakevol)

    USE STATETYPE_MODULE, ONLY : lakestatetype
    USE SURFACEWATER_PROCESSES, ONLY : calculate_lake_volume
      
    !Argument declarations
    INTEGER, INTENT(IN) :: n            !<number of subbasins
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<Lake state, submodel
    REAL, INTENT(INOUT) :: wlakevol(n)  !<volume of olake/volume for lake with basins in outlet basin (Mm3)
    REAL, INTENT(OUT)   :: ilakevol(n)  !<volume of ilake (m3)
    REAL, INTENT(OUT)   :: olakevol(n)  !<volume of olake (m3)
    REAL, INTENT(OUT)   :: Milakevol(n) !<volume of ilake (Mm3), missing_value for no lake
    REAL, INTENT(OUT)   :: Molakevol(n) !<volume of olake (Mm3), missing_value for no lake
    
    !Variable declarations
    INTEGER i,j,itype
    REAL a,fpfrac
    REAL lakearea(nlaketypes)       !lake area (m2) for lakes
    REAL lakebasinvol(nlaketypes)   !volume of ilake & olake (m3)
    REAL lakevol                    !volume for single olake and lakes composed of lakebasins(sum)  (m3)
    REAL lakevolsum(nbasinlakes)    !accumulation of lake basin volumes to outlet during the i-loop (m3)
    
    !Initialisation: wlakevol (outvar) is missing_value when reaching this subroutine already
    ilakevol = 0.
    olakevol = 0.
    Milakevol = missing_value
    Molakevol = missing_value
    lakevolsum = 0.
    
    !Calculate lake volumes for each subbasin
    DO i=1,n
      
      !Calculate lake area for ilake and olake (without floodplain area)
      lakearea = 0.
      DO itype = 1,nlaketypes
        a = 0.
        IF(itype==1 .AND. slc_ilake>0 .OR. itype==2 .AND. slc_olake>0)THEN
          IF(itype==1) j=slc_ilake
          IF(itype==2) j=slc_olake
          a = classbasin(i,j)%part
        ENDIF
        IF(a>0)THEN      !lake exist
          lakearea(itype) = a * basin(i)%area                  !m2
          IF(itype==2 .AND. conductflood)THEN
            IF(floodindex(i)>0)THEN
              fpfrac = flooding(floodindex(i))%fpfol          !floodplain fraction of outlet lake area
              IF(fpfrac>0.) lakearea(itype) = a * basin(i)%area * (1.-fpfrac)    !m2
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      
      !Calculate lake volume for print out
      lakebasinvol=0.
      lakevol=0.
      IF(slc_ilake>0. .AND. lakearea(1)>0.)THEN
        IF(conductN.OR.conductP.OR.conductC.OR.conductT)THEN
          CALL calculate_lake_volume(1,i,lakearea(1),lakestate%water(1,i)+lakestate%slowwater(1,i),lakebasinvol,lakevol,lakevolsum)
        ELSE
          CALL calculate_lake_volume(1,i,lakearea(1),lakestate%water(1,i),lakebasinvol,lakevol,lakevolsum)
        ENDIF
        Milakevol(i) = lakebasinvol(1)*1.E-6  !Mm3
      ENDIF
      IF(slc_olake>0 .AND. lakearea(2)>0.)THEN
        IF(conductN.OR.conductP.OR.conductC.OR.conductT)THEN
          CALL calculate_lake_volume(2,i,lakearea(2),lakestate%water(2,i)+lakestate%slowwater(2,i),lakebasinvol,lakevol,lakevolsum)
        ELSE
          CALL calculate_lake_volume(2,i,lakearea(2),lakestate%water(2,i),lakebasinvol,lakevol,lakevolsum)
        ENDIF
        IF(lakevol.NE.missing_value)THEN
          wlakevol(i) = lakevol/1000000. !whole lake volume for lakebasinlake (10^6 m3)
        ELSE
          wlakevol(i) = lakevol          !upstream lake basin is set to missing
        ENDIF
        Molakevol(i) = lakebasinvol(2)*1.E-6  !Mm3
      ENDIF
      
      !Additional output variables for lake volume (m3)
      ilakevol(i) = lakebasinvol(1)
      olakevol(i) = lakebasinvol(2)
    ENDDO
    
  END SUBROUTINE calculate_lake_volume_output
    
  !!>Calculate weights for subbasins to be used for calculation of upstream average
  !!--------------------------------------------------------------------
  !SUBROUTINE calculate_areaweights_for_upstream_average()
  !
  !  USE MODVAR, ONLY : basin,path,upstreamarea,landarea,branchdata,branchindex,nsub,upstreamweight
  !
  !  !Local variables
  !  INTEGER isb,i2,iarea
  !  REAL upweight(4,nsub),totweight(nsub,nsub)
  !  REAL divasum5,divasum6
  !
  !  ALLOCATE(upstreamweight(nsub,nsub,2))
  !
  !  upweight = 0.
  !  DO isb = 1,nsub
  !    upweight(1,isb) = basin(isb)%area / upstreamarea(isb)
  !    upweight(2,isb) = landarea(isb) / upstreamarea(isb)   !resulting upave will be basinarea!
  !    IF(path(isb)%main>0)THEN
  !      divasum6 = upstreamarea(isb) / upstreamarea(path(isb)%main)
  !      IF(ALLOCATED(branchdata))THEN
  !        IF(branchindex(isb)>0)THEN
  !          IF(branchdata(branchindex(isb))%mainpart<1.) divasum6 = branchdata(branchindex(isb))%mainpart * divasum6  !OBS: No consideration of varying branchflow fraction
  !        ENDIF
  !      ENDIF
  !      upweight(3,isb) = divasum6
  !    ENDIF
  !    IF(ALLOCATED(branchdata))THEN
  !      IF(branchindex(isb)>0)THEN
  !        IF(branchdata(branchindex(isb))%branch>0)THEN
  !          divasum5 = (1.- MIN(1.,branchdata(branchindex(isb))%mainpart)) * upstreamarea(isb) / upstreamarea(branchdata(branchindex(isb))%branch)  !OBS: No consideration of varying branchflow fraction
  !          IF(divasum5>0.)THEN
  !            upweight(4,isb) = divasum5
  !          ENDIF
  !        ENDIF
  !      ENDIF
  !    ENDIF
  !  ENDDO
  !  
  !  !Use weights to calculate total weights
  !  DO iarea=1,2  !areatype=2  !landarea
  !    totweight = 0.
  !    DO isb = 1,nsub
  !      totweight(isb,isb)=upweight(iarea,isb)
  !      IF(upweight(3,isb)>0.)THEN
  !        totweight(path(isb)%main,isb) = upweight(3,isb)*totweight(isb,isb)
  !        DO i2 = 1,isb-1
  !          IF(totweight(isb,i2)>0) totweight(path(isb)%main,i2) = totweight(path(isb)%main,i2) + upweight(3,isb)*totweight(isb,i2)
  !        ENDDO
  !      ENDIF
  !      IF(upweight(4,isb)>0.)THEN
  !        totweight(branchdata(branchindex(isb))%branch,isb) = upweight(4,isb)*totweight(isb,isb)
  !        DO i2 = 1,isb-1
  !          IF(totweight(isb,i2)>0) totweight(branchdata(branchindex(isb))%branch,i2) = totweight(branchdata(branchindex(isb))%branch,i2) + upweight(4,isb)*totweight(isb,i2)
  !        ENDDO
  !      ENDIF
  !    ENDDO
  !    upstreamweight(:,:,iarea) = totweight(:,:)
  !  ENDDO
  !
  !END SUBROUTINE calculate_areaweights_for_upstream_average

  !>Initialize calulation of subbasin mean value from classes
  !> \b Consequences Module modvar variable outvar may change
  !--------------------------------------------------------------------
  SUBROUTINE calculate_class_outvar_initialize(idindex,isb)

    USE MODVAR, ONLY : outvar, & !OUT
                       outvarindex
  
    INTEGER, INTENT(IN) :: idindex  !<index in outvarid
    INTEGER, INTENT(IN) :: isb      !<index of subbasin

    IF(outvarindex(idindex)>0) outvar(isb,outvarindex(idindex))=0.
  
  END SUBROUTINE calculate_class_outvar_initialize

  !>Add class value to calculate subbasin mean value
  !> \b Consequences Module modvar variable outvar may change
  !--------------------------------------------------------------------
  SUBROUTINE calculate_class_outvar_add(idindex,isb,fraction,value)
  
    USE MODVAR, ONLY : outvar, & !OUT
                       outvarindex
  
    INTEGER, INTENT(IN) :: idindex  !<index in outvarid
    INTEGER, INTENT(IN) :: isb      !<index of subbasin
    REAL, INTENT(IN)    :: fraction !<area fraction of class
    REAL, INTENT(IN)    :: value    !<class value of output variable
  
    IF(outvarindex(idindex)>0) outvar(isb,outvarindex(idindex)) = &
               outvar(isb,outvarindex(idindex)) + value*fraction
            
  END SUBROUTINE calculate_class_outvar_add

  !>Finish calculate subbasin mean value from classes by land area (or other area fraction)
  !> \b Consequences Module modvar variable outvar may change
  !--------------------------------------------------------------------
  SUBROUTINE calculate_class_outvar_finish(idindex,isb,fraction)
  
    USE MODVAR, ONLY : outvar, & !OUT
                       outvarindex
  
    INTEGER, INTENT(IN) :: idindex  !<index in outvarid
    INTEGER, INTENT(IN) :: isb      !<index of subbasin
    REAL, INTENT(IN)    :: fraction !<land area fraction (or other) of subbasin area
  
    IF(outvarindex(idindex)>0)THEN
      IF(fraction==0.)RETURN
      outvar(isb,outvarindex(idindex)) = &
             outvar(isb,outvarindex(idindex)) / fraction
      !IF(fraction>0.)THEN                    !Denna ger NaN då bar ilake??
      !  outvar(isb,outvarindex(idindex)) = &
      !         outvar(isb,outvarindex(idindex)) / fraction
      !ELSE
      !  outvar(isb,outvarindex(idindex)) = 0.
      !ENDIF
    ENDIF
    
  END SUBROUTINE calculate_class_outvar_finish

  !>Finish calculate subbasin value from classes by multiplying by factor
  !> \b Consequences Module modvar variable outvar may change
  !--------------------------------------------------------------------
  SUBROUTINE calculate_class_outvar_finish_scale(idindex,isb,fraction)
  
    USE MODVAR, ONLY : outvar, & !OUT
                       outvarindex
  
    INTEGER, INTENT(IN) :: idindex  !<index in outvarid
    INTEGER, INTENT(IN) :: isb      !<index of subbasin
    REAL, INTENT(IN)    :: fraction !<land area fraction (or other) of subbasin area
  
    IF(outvarindex(idindex)>0) outvar(isb,outvarindex(idindex)) = &
               outvar(isb,outvarindex(idindex)) * fraction
            
  END SUBROUTINE calculate_class_outvar_finish_scale

  !>Calculate temperature of outflow from subbasin from river or lake outflow
  !> \b Consequences Module modvar variable outvar may change
  !--------------------------------------------------------------------
  SUBROUTINE calculate_outvar_watertemperature(idindex,isb,itype,lowerlimit,riverstate,lakestate)
  
    INTEGER, INTENT(IN) :: idindex  !<index in outvarid
    INTEGER, INTENT(IN) :: isb      !<index of subbasin
    INTEGER, INTENT(IN) :: itype    !<index of system (1=local, 2=main)
    REAL, INTENT(IN)    :: lowerlimit !<lower limit of temperature
    TYPE(riverstatetype),INTENT(IN) :: riverstate   !<River states
    TYPE(lakestatetype),INTENT(IN)  :: lakestate    !<Lake states
  
    IF(outvarindex(idindex)>0)THEN
      outvar(isb,outvarindex(idindex)) = riverstate%temp(itype,isb)
      IF(slc_olake>0)THEN
        IF(classbasin(isb,slc_olake)%part>0) outvar(isb,outvarindex(idindex)) = lakestate%temp(itype,isb)
      ENDIF
      IF(outvar(isb,outvarindex(idindex))<lowerlimit) outvar(isb,outvarindex(idindex)) = lowerlimit
    ENDIF
            
  END SUBROUTINE calculate_outvar_watertemperature
 
  !>Set output variable from Xobs
  !> \b Consequences Module modvar variable outvar may change
  !--------------------------------------------------------------------
  SUBROUTINE set_outvar_xobs(idindex,isb)
      
    INTEGER, INTENT(IN) :: idindex  !<index in outvarid
    INTEGER, INTENT(IN) :: isb      !<index of subbasin

    IF(outvarindex(idindex)>0)THEN
      IF(xobsindex(idindex,isb)>0)THEN
        outvar(isb,outvarindex(idindex)) = xobsi(xobsindex(idindex,isb))
      ELSE
        outvar(isb,outvarindex(idindex)) = missing_value
      ENDIF
    ENDIF
  
  END SUBROUTINE set_outvar_xobs

  !>Set output variable from a scaled Xobs variable
  !> \b Consequences Module modvar variable outvar may change
  !--------------------------------------------------------------------
  SUBROUTINE set_outvar_xobs_scaled(idindex,isb,idxobs,factor)
      
    INTEGER, INTENT(IN) :: idindex  !<index in outvarid for output variable
    INTEGER, INTENT(IN) :: isb      !<index of subbasin
    INTEGER, INTENT(IN) :: idxobs   !<index in outvarid for xobs variable
    REAL, INTENT(IN) :: factor      !<scale factor for output

    IF(outvarindex(idindex)>0)THEN
      outvar(isb,outvarindex(idindex)) = missing_value
      IF(xobsindex(idxobs,isb)>0)THEN
        IF(xobsi(xobsindex(idxobs,isb)).NE.missing_value) outvar(isb,outvarindex(idindex)) = xobsi(xobsindex(idxobs,isb))*factor
      ENDIF
    ENDIF
  
  END SUBROUTINE set_outvar_xobs_scaled

  !>Set output variable for a meantype Xobs/Xom variable
  !> \b Consequences Module modvar variable outvar may change
  !--------------------------------------------------------------------
  SUBROUTINE set_outvar_xobsmean(idindex,isb,idxobs)

    USE HYPE_INDATA, ONLY : num_xoms, xom

    INTEGER, INTENT(IN) :: idindex  !<index in outvarid for output variable
    INTEGER, INTENT(IN) :: isb      !<index of subbasin
    INTEGER, INTENT(IN) :: idxobs   !<index in outvarid for xom variable
  
    IF(outvarindex(idindex)>0)THEN
      IF(xobsindex(idindex,isb)>0)THEN
        outvar(isb,outvarindex(idindex)) = xobsi(xobsindex(idindex,isb))
      ELSEIF(num_xoms(1)>=idxobs)THEN
        outvar(isb,outvarindex(idindex)) = xom(idxobs,isb)
      ELSE
        outvar(isb,outvarindex(idindex)) = missing_value
      ENDIF
    ENDIF
   
  END SUBROUTINE set_outvar_xobsmean

  !>Set output variable from a statetype Xobs/Xos variable
  !> \b Consequences Module modvar variable outvar may change
  !--------------------------------------------------------------------
  SUBROUTINE set_outvar_xobsstate(idindex,isb,idxobs)

    USE HYPE_INDATA, ONLY : num_xoms, xos

    INTEGER, INTENT(IN) :: idindex  !<index in outvarid for output variable
    INTEGER, INTENT(IN) :: isb      !<index of subbasin
    INTEGER, INTENT(IN) :: idxobs   !<index in outvarid for xos variable
  
     IF(outvarindex(idindex)>0)THEN
      IF(xobsindex(idindex,isb)>0)THEN
        outvar(isb,outvarindex(idindex)) = xobsi(xobsindex(idindex,isb))
      ELSEIF(num_xoms(2)>=idxobs)THEN
        outvar(isb,outvarindex(idindex)) = xos(idxobs,isb)
      ELSE
        outvar(isb,outvarindex(idindex)) = missing_value
      ENDIF
    ENDIF

  END SUBROUTINE set_outvar_xobsstate

  !>Calculate regional groundwater flow to outside model domain (for deepmodel 1)
  !--------------------------------------------------------------------
  SUBROUTINE calculate_regional_groundwaterflow_to_outside_system(isb,totflow,outflow)

    INTEGER, INTENT(IN) :: isb       !<index of subbasin
    REAL, INTENT(IN)    :: totflow   !<total flow of regional groundwater (m3/ts)
    REAL, INTENT(OUT)   :: outflow   !<flow to outside of model domain (m3/ts)

    outflow = 0.
    IF(modeloption(p_deepgroundwater)==1)THEN
      IF(path(isb)%grw1==0) outflow = totflow   !loss of water from the model system via groundwater flow
    ENDIF
    
  END SUBROUTINE calculate_regional_groundwaterflow_to_outside_system

END MODULE

