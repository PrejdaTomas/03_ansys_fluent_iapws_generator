The file extracts material data for Ansys Fluent from iapws IF97 standard revided in 2007.
The material data will be split into two submaterials: liquid material and vapour material, which are saved into an Fluent Material Database (SCHEME) file-
  - density
  - dynamic viscosity
  - thermal conductivity
  - isobaric thermal capacity

Afterwards the script generates the saturation properties for the evaporation-condensation model in Fluent (the saturation curve and surface tension/temperature dependency) - just copy the sequence into Fluent's console.


The script also generates CSV files and optionally it also generates the plots
The material data are generated in °C and bars - set your Fluent units up accordingly.

Requisities:
 - iapws (pip3 install iapws)
 - os
 - argparse
 - matplotlib (optional)

Usage:
 - put the file into a folder
 - run it from cmd: python3 pyIAPWS_toFluent.py NAME LEVELS PRESSURE_MPAS TLOWER_K TUPPER_K --plotting
    - NAME: name of the output files without suffixes (e.g. water-if97)
    - LEVELS: number of rows for the tables (Fluent cannot accept more than 50 values, but the vapour and liquid are separate materials, so you can try a higher number - in worst case you will get an exception if the submaterial gets more than 50 levels)
    - PRESSURE_MPAS: absolute pressure in MPa(abs)
    - TLOWER_K: lower temperature bound in K
    - TUPPER_K: upper temperature bound in K
    -    --plotting (optional, if not inputted, no plots are generates)


Notes:
 - tested on Fluent 18
 - latest version causes Fluent to Freeze - need to debug against the Fluent Material Database Files
 - Next plans:
   1. unit selection
   2. maybe an UDF with IAPWS itself for Fluent
