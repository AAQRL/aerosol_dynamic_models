# Discrete Model

## Introduction
The codes in this folder are to calculate the evolution of particle size distribution with initial bi-lognormal particle size distribution.

## Modules

- **main.f90**: main program
- **user\_interface.f90**: input source code
- **global\_variables.f90**: define the global variables
- **beta_safe.f90**: calculate the coagulation kernel
- **rhs.f90**: discrete the Smoluchowski equation
- **results.f90**: output source code
- **DVODE\_F90_M.f90**: solve the ODE equations
- **function\_list.f90**: calculate mean variables such as d<sub>pg</sub>, segma<sub>g</sub> of the particle number distribution

## Inputs

- **input.tx**t includes:
	- Flag for coagulation (0 means no coagulation; 1 means coagulation)
	- Number of discrete sizes
	- Molecular Weight of the species, Density of the species, diameter of the first discrete size
	- Initial Number Concentration (#/m<sup>3</sup>)
	- Final time in seconds 
	- dt : time step
	- Minimum Number Concentration Tracking for accuracy. This will be used to find absolute and relative tolerance

## Outputs
- **output\_aver\_para.dat**
	- **3rd line**: t0 d1, d2, d3, \.\.\., dn
	- **4th line**: t0 N(d1), N(d2), N(d3), \.\.\., N(dn)  !N is the particle number concentraton (\#/m<sup>3</sup>)
	- **5th line**: t\_final N(d1), N(d2), N(d3), \.\.\., N(dn)  !N is the particle number concentraton (\#/m<sup>3</sup>)
- **output\_CPU&Memory.dat** : record the total computational time that is cost
- **output\_aver\_para.dat**
	- **Ntot**: total number concentration (\#/cm<sup>3</sup>)
	- **dpg**: geometric mean diameter (nm) (Friedlander, 2000, pp**:17)
	- **segma\_g**: geometric standard deviation (Friedlander, 2000, pp**:17)
	- **Surf\_tot**: total particle surface area  (m<sup>3</sup>/cm<sup>3</sup>)
	- **Vol\_tot**: total particle volume (m<sup>3</sup>/cm<sup>3</sup>) 
