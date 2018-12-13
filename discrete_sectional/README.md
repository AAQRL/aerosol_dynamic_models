# DSM_f90

## Introduction


## Modules

- **discsecpb.f90**: main program
- **preci.f90**: precision configuration
- **globals.f90**: define the global variables shared between the main program and other modules including **funcs, calbeta\_mod, xdoteq\_mod, ot\_mod**
- **calbeta\_common.f90**: define the common variables shared between the CALBETA subroutine and other modules including **funcs, calbeta\_mod, xdoteq\_mod**
- **aux\_data\_funcs.f90**: functions for pulling data
- **funcs.f90**: functions used in CALBETA subroutine and DTWODQ sobroutine
- **quad\_glegq.f90**: function and subroutine used in CALBETA subroutine
- **dtwodq.f90**: integral subroutine
- **calbeta\_mod.f90**: determine the coagulation and condensation coefficients in the free molecular regime
- **ot\_mod.f90**: write results to the output file
- **xdoteq\_mod.f90**: determine the differential equations required by RK45
- **rk45\_init0.f90**: Runge-Kutta integration scheme

## Module calling tree

discsecpb (**discsecpb.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- TEMPHIST(**aux\_data\_funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- MTSAT(**aux\_data\_funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- RTDATA(**aux\_data\_funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- NUCLDATA(**aux\_data\_funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- CALBETA(**calbeta\_mod.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- QUAD(**quad\_glegq.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- GLEGQ(**quad\_glegq.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- FX(**funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- DTWODQ(**dtwodq.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- F(**funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- G(**funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- H(**funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- OT(**ot\_mod.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- INIT0(**rk45\_init0.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- TEMPHIST(**aux\_data\_funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- MTSAT(**aux\_data\_funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- RTDATA(**aux\_data\_funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- NUCLDATA(**aux\_data\_funcs.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- RK45(**rk45\_init0.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- XDOTEQ(**xdoteq\_mod.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- INIT0(**rk45\_init0.f90**)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-- OT(**ot\_mod.f90**)<br>

## Inputs

- **SYSOPH.DAT** includes:
- **TMPHIST.DAT** includes:


## Outputs
- **output.dat**:
- **confc.dat**: 
