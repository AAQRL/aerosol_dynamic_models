# Modal Model

## Introduction
The codes in this folder are to calculate the evolution of particle size distribution as the function of time with any initial particle size distribution.

## Modules

- **Modal_GDM.f90**: main program
- **global\_variables.f90**: define the global variables
- **read\_input.f90**: read the data from input files
- **derived\_parameters.f90**: contains several functions for coagulation coefficient, function for redistribution of particles between bins and subroutine to evaluate temperature dependent variables such as mean free path and viscosity 
- **DVODE\_F90_M.f90**: module for ODE solver subroutine
- **rhs.f90**: ODE for the Smoluchowski equation (Modal)
- **calculations.f90**: calculate total number concentration and mean variables such as d<sub>pg</sub>, sigma<sub>g</sub> of the particle number distribution as the function of time

## Inputs

- **input.txt includes:
	- INITIAL PARTICLE PROPERTIES (Initial Vapor Pressure, Saturation Pressure and Vapor Generation Rate)
	- SPECIES PROPERTIES (Molecular Weight of the species, Density of the species and Surface Tension)
	- GAS PROPERTIES (Molecular Weight of gas, Density of gas at 298 K, Viscosity of gas at Tref = 298 K, Sutherland's Constant and Mean free path of gas molecules)
	- SYSTEM PROPERTIES (Pressure and Temperature)
	- SIMULATION TIME (Final time in seconds and Time step)
	- MODES (Number of modes)
	- AEROSOL PROCESSES (if ON set to 1) (Coagulation, Condensation and Nucleation)

## Outputs
- **results\_raw.dat**
	- **1st column**: Time in seconds
	- **Last column**: Vapor Concentration
	- **2k Column**: Number Concentration based in Mode k
	- **(2k+1) Column**: Total Volume of Mode k
	
- **results.dat**
	- **1st column**: Time in seconds
	- **2nd column**: Total Number Concentration (#/m3)
	- **3rd column**: Geometric Mean Diameter (m)
	- **4th column**: Geometric Standard Diameter