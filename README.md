The file extracts material data for Ansys Fluent from iapws IF97 standard revided in 2007.
The material data will be split into two submaterials: liquid material and vapour material, which are saved into an Fluent Material Database (SCHEME) file.

The steam table fails in (super)critical state:
  - TCRIT       >= 647.096 K and PCRIT       >= 220.64 MPa
    
The material data inputted into the file are:
  - density
  - dynamic viscosity
  - thermal conductivity
  - isobaric thermal capacity
  - formation enthalpy (sadly it is not the standard-state-enthalpy: you have to input it yourself (the material database file is not able to load this property), 0 for liquid and latent heat for vapour: you will get the values in console output
  - molar mass

Afterwards the script generates the saturation properties for the evaporation-condensation model in Fluent (the saturation curve and surface tension/temperature dependency) - just copy the sequence into Fluent's console.


The script also generates CSV files and optionally it also generates the plots
The material data are generated in °C and bars - set your Fluent units up accordingly.

Required modules:
 - iapws (pip3 install iapws)
 - os
 - argparse
 - subplotterFuncs (optional, contains matplotlib)
 - matplotlib (optional)

Usage:
 - put the file into a folder
 - run it from cmd: python3 pyIAPWS_toFluent.py NAME LEVELS PRESSURE_MPAS TLOWER_K TUPPER_K --plotting
    - NAME: name of the output files without suffixes (e.g. water-if97)
    - LEVELS: number of rows for the tables (Fluent cannot accept more than 50 values, but the vapour and liquid are separate materials, so you can try a higher number - in worst case you will get an exception if the submaterial gets more than 50 levels)
    - PRESSURE_MPAS: absolute pressure in MPa(abs)
    - TLOWER_K: lower temperature bound in K
    - TUPPER_K: upper temperature bound in K
    -    --plotting (optional, if not inputted, no plots are generated: if not inputted, matplotlib is not necessary)


Notes:
 - tested on Ansys Fluent 14 and Ansys Fluent 18
 - latest version causes Fluent to freeze - need to debug against the Fluent Material Database Files
 - Next plans:
   1. unit selection for pressure and temperature
   2. maybe an UDF with IAPWS itself for Fluent
