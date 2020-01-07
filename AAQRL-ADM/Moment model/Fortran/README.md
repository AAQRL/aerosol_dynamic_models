# Discrete Model

## Introduction
The codes in this folder are to calculate the evolution of particle size distribution with initial bi-lognormal particle size distribution.

## Modules

- **main.f90**: main program
- **user\_interface.f90**: input source code
- **global\_variables.f90**: define the global variables
- **beta_safe.f90**: calculate the coagulation kernel
- **rhs.f90**: f(x,t) for the ODEs
- **results.f90**: output source code
- **DVODE\_F90_M.f90**: solve the ODE equations
- **function\_list.f90**: calculate mean variables such as d<sub>pg</sub>, segma<sub>g</sub> of the particle number distribution

## Inputs

- **input1.txt** includes:
	- Ntot (#m-3); dpg (m); sigmag
	- Saturation Pressure
	- Molecular Weight of the species (Kg/mol), Density of the species (Kg/m3), surface tension (Nm)
	- Flag for coagulation (1 means TRUE, 0 means FALSE)
	- Flag for condensation (1 means TRUE (using vapor condensation), 0 means FALSE, 2 means surface growth through reaction)
	- Flag for nucleation (1 means TRUE & Classical nucleation theory, 2 means reaction to form particles directly,  0 means FALSE)
	- Flag for reaction (1 means TRUE, 0 means FALSE)
	- dimensionless reaction rate

## Outputs
- **output.txt**
	- **1st column**: z_pos
	- **1st column**: T (K)
	- **1st column**: Number concentration	(\#/cm<sup>3</sup>)
	- **1st column**: dpmg	(nm)	
	- **1st column**: sigmag	
	- **1st column**: vapor_conc	(\#/cm<sup>3</sup>)