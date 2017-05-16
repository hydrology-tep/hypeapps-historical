
SET EXE=HYPE.exe
SET F90FILES=hypevar.f90 libdate.f90 hypetypes.f90 modvar.f90 worvar.f90 general_wc.f90 general_func.f90 convert.f90 time.f90 readwrite.f90 hype_indata.f90 atm_proc.f90 irrigation.f90 hype_wbout.f90 npc_soil_proc.f90 soil_proc.f90 regional_groundwater.f90 sw_proc.f90 npc_sw_proc.f90 soilmodel0.f90 glacier_soilmodel.f90 soilmodel4.f90 model_hype.f90 compout.f90 data.f90 optim.f90 statedata.f90 main.f90

::For Intel Fortran
::SET OPT=/warn:all /check:bounds /traceback /heap-arrays /Qauto_scalar
SET OPT=/O3 /Qipo /Qauto_scalar
SET FC=IFORT

::For gfortran
::SET OPT=-fbacktrace -ffree-line-length-none -Wall
::SET OPT=-O2 -w -ffree-line-length-0
::SET FC=gfortran

%FC% %OPT% %F90FILES% /exe:%EXE%

