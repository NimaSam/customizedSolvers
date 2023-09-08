# customized Solvers and libs

Here I share different solvers and library i cusotmized for different Projects


## Solver list
- driftFluxFoam (OF9)
  - Add RichardsonZaki relative velocity model
- chtInterHeatFoam (OF2.3)
  - A conjugate heat transfer solver for multiphaseFlow based on interFoam and multiRegionFoam.
     

### Example
- driftFluxFoam: particle settling simulation in pipe with driftFluxFoam
    - terminalVelocity
    - hinderVelocity: The settling velocity is calculated by multiplying the terminal velocity with the Richardson-Zaki hinder function.
    - 
- chtInterHeatFoam: drop impact on heated surface


## Library list
- DiFelice drag force model for particle flow (OF9)


