import os, argparse
from numpy import interp as lininterp1D, vectorize, array as npArray
from iapws import IAPWS97
from iapws.iapws97 import _TSat_P as getSaturationTemperature
from iapws.iapws97 import _PSat_T as getSaturationPressure_uncorrected

#region MSG_GEN
msg = "\n"+10*"_"
msg += "A IAPWS pf97 steam table generator for Ansys Fluent.\n"
msg += "To generate the tables, input the following arguments\n"
msg += "path: address of the newly created files (without extension)\n"
msg += "levels: positive integer (resolution of the table)\n"
msg += "pabs: positive float (absolute pressure [MPa a])\n"
msg += "tlow: positive float (lower bound of the temperature range [K])\n"
msg += "thigh: positive float (upper bound of the pressure range [K])\n"
msg += "--plotting: True/False (will plot the properties if True)"
msg = 10*"_" + "\n"
#endregion MSG_GEN

# region CONSOLE_ARGUMENTS
# Initialize parser
parser = argparse.ArgumentParser(description = msg)

parser.add_argument("path",
                    action = "store",
                    help = "Napis adresu (se jmenem souboru), kam chces vystupy ulozit (csv, scm)",
                    default= os.getcwd()
                    )

parser.add_argument("levels",
                    action = "store",
                    help = "Zadej pocet urovni pro generaci parni tabulky (prirozene kladne cislo)",
                    #default=int(input("Zadej pocet urovni pro generaci parni tabulky (prirozene kladne cislo)"))
                    )

parser.add_argument("pabs",
                    action = "store",
                    help = "Zadej absolutni tlak [MPa a]",
                    #default= float(input("Zadej absolutni tlak [MPa a]: ").strip())
                    )

parser.add_argument("tlow",
                    action = "store",
                    help = "Zadej dolni mez termodynamicke teploty [K]",
                    #default= float(input("Zadej dolni mez termodynamicke teploty [K]: ").strip())
                    )

parser.add_argument("thigh",
                    action = "store",
                    help = "Zadej horni mez termodynamicke teploty [K]",
                    #default= float(input("Zadej horni mez termodynamicke teploty [K]: ").strip())
                    )

parser.add_argument("--plotting",
                    action = "store_true",
                    help = "Napis --plotting pro vykresleni grafu",
                    default= False
                    )

args = parser.parse_args()
#endregion CONSOLE_ARGUMENTS

print("Vlozil jsi nasledujici argumenty:")
for arg in vars(args): print("\t", arg, getattr(args, arg))




linspace = lambda lo,hi,steps: [lo + step*(hi-lo)/(steps-1) for step in range(steps)]
def getSaturationPressure(temperature: float) -> float|str:
    """
    Gets the saturation pressure for selected temperature in bar(a) from IAPWS if97
    Args:
        temperature (float): saturation temperature in K

    Returns:
        float|str: saturation pressure [bar(a)] or "inf" if above critical temperature  647.096 K
    """
    try: return getSaturationPressure_uncorrected(temperature)*10
    except: return "inf"

def getFluentProperties(iapwsCoordinates: IAPWS97)-> tuple[float, float, float, float, float, float, float, str]:
    """Extracts the material properties from the selected IAPWS97 object

    Args:
        iapwsCoordinates (IAPWS97): IAPWS97(Temperature[K], Pressure[MPa])

    Returns:
        - tuple[float, float, float, float, float, float, float, str]: \n
            - density (float): kg m-3\n
            - dynamic viscosity (float): Pa s-1\n
            - thermal conductivity (float): W m-1 K-1\n
            - specific enthalpy (float): J mol-1\n
            - specific entropy (float): J mol-1 K-1\n
            - isotropic linear thermal expansion coefficieny (float): m-1\n
            - phase (str): Liquid/Vapour/Supercritical
    """
    density = iapwsCoordinates.rho
    viscosity = iapwsCoordinates.mu
    isobaricThermalCapacity = iapwsCoordinates.cp*1000
    thermalConductivity = iapwsCoordinates.k
    specificEnthalpy = iapwsCoordinates.h*iapwsCoordinates.M*1000
    specificEntropy = iapwsCoordinates.s*iapwsCoordinates.M*1000
    thermalExpansionCoefficient = iapwsCoordinates.alfav # linear
    phase = iapwsCoordinates.phase
    #print(phase, int(iapwsCoordinates.h*1000), int(specificEnthalpy), )
    return density, viscosity, isobaricThermalCapacity, thermalConductivity, specificEnthalpy, specificEntropy, thermalExpansionCoefficient, phase

def fluentIAPWSgen(temperature: list[float],
                   density: list[float],
                   viscosity: list[float],
                   specific_heat: list[float],
                   thermal_conductivity: list[float],
                   specific_enthalpy: list[float],
                   reference_temperature: float = None,
                   standard_state_enthalpy: float = None,
                   name: str=None) -> str:
    """Generates the water-vapour material for use in Ansys Fluent (scm file to be imported as material database).\n

    Args:
        - temperature: list of float values (temperature: K)
        - density: list of float values (density: kg m-3)
        - viscosity: list of float values (dynamic viscosity: Pa s)
        - specific_heat: list of float values (isobaric heat capacity: J kg-1 K-1)
        - thermal_conductivity: list of float values (thermal conductivity: W m-1 K-1)
        - specific_enthalpy: list of float values (specific enthalpy: J kg-1 mol-1) <this unit from Fluent is really funky>
        - standard_state_enthalpy: float in J kg-1 mol-1 (temperature to which standard state enthalpy relates to)
        - reference_temperature: float in K (temperature to which enthalpy relates to)
        - name: name of the material

    Returns:
        str: scheme representation of the material properties
    """

    def scheme(level: int=0, ipt: str= None) -> str:
        return level * "\t" + "({0}\n".format(ipt)

    def propertyWriter(property: str = None, t: float = None, val: float = None) -> str:
        strOutput = scheme(2, property)
        strOutput += scheme(3, "polynomial piecewise-linear ")
        for tup in zip(t, val):
            strOutput += 4*"\t"
            strOutput += "({0} . {1})".format(tup[0], tup[1])
            strOutput += "\n"
        strOutput += 3*"\t"+")\n"
        strOutput += scheme(3, "constant . {0})".format(sum(val)/len(val)))
        strOutput += 2*"\t"+")\n"
        return strOutput 
    schemeString = "\n"
    #schemeString =  scheme(0, "")
    schemeString += scheme(1, "{0} fluid".format(name))
    schemeString += scheme(2, "chemical-formula . #f)")
    schemeString += propertyWriter("density", temperature, density)
    schemeString += propertyWriter("specific-heat", temperature, specific_heat)
    schemeString += propertyWriter("thermal-conductivity", temperature, thermal_conductivity)
    schemeString += propertyWriter("viscosity", temperature, viscosity)
    schemeString += propertyWriter("formation-enthalpy", temperature, specific_enthalpy)
    schemeString += 2*"\t" + ";(standard-state-enthalpy (constant . {0}))\n".format(standard_state_enthalpy)
    schemeString += 2*"\t" + "(reference-temperature (constant . {0}))\n".format(reference_temperature)
    schemeString += 2*"\t" + "(molecular-weight (constant . 18.0153))\n"
    schemeString += 1*"\t"+")\n"
    #schemeString += 0*"\t"+")\n"
    return schemeString

def getClosestIndexInArray(value: float,
                           array: list[float]) -> int|str:
    """Finds the closest value's index to the searched one in searched list

    Args:
        - value (float): value to which we try to find the closest value's index. Defaults to None.
        - array (list[float]): searched array. Defaults to None.

    Returns:
        int|str: index of the closest value in the array or "N/A" if the value is in list
    """
    diffs: list[float]  = [abs(array_val - value) for array_val in array]
    closestIndex: int   = diffs.index(min(diffs))
    searchedVal: float  = array[closestIndex]
    if      value < searchedVal: return closestIndex + 0
    elif    value > searchedVal: return closestIndex + 1
    else: return "N/A"

def insertValToClosest(value: float,
                       array: list[float]) -> None:
    """Attempts to insert the float value into the list so it fits in with no sudden jump. Modifies the array!
    Args:
        - value (float): value to which we try to add. Defaults to None.
        - array (list[float]): searched array. Defaults to None.
    """
    desiredIndex = getClosestIndexInArray(value, array)
    if not desiredIndex == "N/A":
        array.insert(desiredIndex, value)
    else:
        global temperatureLower, temperatureUpper, levels
        print("\nThe saturated temperature is already in the list")
        array = linspace(temperatureLower, temperatureUpper, levels)

def closeArray(valType: str, temperatures: list[float], *args:tuple[list[float]]) -> None:
    """Closes the array values for Fluent by copying the first value - for 1 K, and by copying the last vlaue - for 5000K.

    Args:
        valType (str): liquid/vapour
        temperatures (list[float]): list of temperatures
        *args: tuple[list[float]]: lists of numeric values where the first and last values will be copied
        
    """
    if not isinstance(valType, str): raise TypeError("You have not correctly specified liquid/vapour valType argument: {0}".format(valType))
    if len(args) == 0: raise Exception("You have not inputted any arguments in the closeArray function besides the temperatures")
    if any([not isinstance(arg, list) for  arg in args]): raise TypeError("You have not inputted lists of float values into the closeArray function")
    temperatures.insert(0, 1.0)
    temperatures.append(5000.0)
    tempCount: int = len(temperatures)
    if tempCount > 50: raise ValueError("You have gotten more than 50 rows of {0} values: {1}\nTry a lower value of -count argument".format(valType,tempCount))
    for arg in args:
        arg.insert(0, arg[0])
        arg.append(arg[-1])

levels              = int(args.levels)
pressureAbsolute    = float(args.pabs)
temperatureLower    = float(args.tlow)
temperatureUpper    = float(args.thigh)
plotting            = bool(args.plotting)
path                = str(args.path)

if levels <= 0: raise ValueError("Pocet urovni musi byt kladne prirozene cislo, zadal jsi: {0}".format(levels))
if pressureAbsolute <= 0: raise ValueError("Zadal jsi zaporny absolutni tlak: {0}".format(pressureAbsolute))
if temperatureLower <= 0: raise ValueError("Zadal jsi zapornou spodni teplotu: {0}".format(temperatureLower))
if temperatureUpper <= 0: raise ValueError("Zadal jsi zapornou horni teplotu: {0}".format(temperatureUpper))
if temperatureUpper <= temperatureLower: raise ValueError("Zadal jsi horni teplotu nizsi nez dolni, dolni: {0} horni: {1}".format(temperatureLower,temperatureUpper))
if temperatureLower >= temperatureUpper: raise ValueError("Zadal jsi dolni teplotu vyssi nez dolni, dolni: {0} horni: {1}".format(temperatureLower,temperatureUpper))

saturationTemperature               = getSaturationTemperature(pressureAbsolute)
temperatureValues                   = linspace(temperatureLower, temperatureUpper, levels - 2)
insertValToClosest(saturationTemperature+0e-3, temperatureValues)
insertValToClosest(saturationTemperature+1e-3, temperatureValues)
saturationPressureValues            = [getSaturationPressure(temp) for temp in temperatureValues] # bar(a)
densityValues                       = []
viscosityValues                     = []
isobaricThermalCapacityValues       = []
thermalConductivityValues           = []
specificEnthalpyValues              = []
specificEntropyValues               = []
linearExpansionCoefficientValues    = []
phaseValues                         = []


with open(path+ ".csv", "w") as filePort:
    filePort.write("T, Psat, ro, mu, cp, lambda, std. enthalpy, std. entropy, lexp. coeff, phase\n")
    for temperature, pressure in zip(temperatureValues, saturationPressureValues):
        properties = getFluentProperties(iapwsCoordinates=IAPWS97(T=temperature, P=pressureAbsolute))
        densityValues.append(properties[0])
        viscosityValues.append(properties[1])
        isobaricThermalCapacityValues.append(properties[2])
        thermalConductivityValues.append(properties[3])
        specificEnthalpyValues.append(properties[4])
        specificEntropyValues.append(properties[5])
        linearExpansionCoefficientValues.append(properties[6])
        phaseValues.append(properties[7])
        filePort.write("{0:.6F}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}\n".format(temperature, pressure, *properties))

#region PHASE_VALS
liquidVal_temperatureValues             = [val[0] for val in zip(temperatureValues, phaseValues)                if "Liquid" == val[1]]
liquidVal_densityValues                 = [val[0] for val in zip(densityValues, phaseValues)                    if "Liquid" == val[1]]
liquidVal_viscosityValues               = [val[0] for val in zip(viscosityValues, phaseValues)                  if "Liquid" == val[1]]
liquidVal_isobaricThermalCapacityValues = [val[0] for val in zip(isobaricThermalCapacityValues, phaseValues)    if "Liquid" == val[1]]
liquidVal_thermalConductivityValues     = [val[0] for val in zip(thermalConductivityValues, phaseValues)        if "Liquid" == val[1]]
liquidVal_specificEnthalpyValues        = [val[0] for val in zip(specificEnthalpyValues, phaseValues)           if "Liquid" == val[1]]
liquidVal_specificEntropyValues         = [val[0] for val in zip(specificEntropyValues, phaseValues)            if "Liquid" == val[1]]
liquidVal_linearExpansionCoefficientValues = [val[0] for val in zip(linearExpansionCoefficientValues, phaseValues)  if "Liquid" == val[1]]

vapourVal_temperatureValues                = [val[0] for val in zip(temperatureValues, phaseValues)                if val[1] in ("Vapour", "Gas")]
vapourVal_temperatureValues[0]             = saturationTemperature
vapourVal_densityValues                    = [val[0] for val in zip(densityValues, phaseValues)                    if val[1] in ("Vapour", "Gas")]
vapourVal_viscosityValues                  = [val[0] for val in zip(viscosityValues, phaseValues)                  if val[1] in ("Vapour", "Gas")]
vapourVal_isobaricThermalCapacityValues    = [val[0] for val in zip(isobaricThermalCapacityValues, phaseValues)    if val[1] in ("Vapour", "Gas")]
vapourVal_thermalConductivityValues        = [val[0] for val in zip(thermalConductivityValues, phaseValues)        if val[1] in ("Vapour", "Gas")]
vapourVal_specificEnthalpyValues           = [val[0] for val in zip(specificEnthalpyValues, phaseValues)           if val[1] in ("Vapour", "Gas")]
vapourVal_specificEntropyValues            = [val[0] for val in zip(specificEntropyValues, phaseValues)            if val[1] in ("Vapour", "Gas")]
vapourVal_linearExpansionCoefficientValues = [val[0] for val in zip(linearExpansionCoefficientValues, phaseValues) if val[1] in ("Vapour", "Gas")]
#endregion PHASE_VALS

closeArray("liquid",
           liquidVal_temperatureValues,
           liquidVal_densityValues,
           liquidVal_viscosityValues,
           liquidVal_isobaricThermalCapacityValues,
           liquidVal_thermalConductivityValues,
           liquidVal_specificEnthalpyValues,
           liquidVal_specificEntropyValues,
           liquidVal_linearExpansionCoefficientValues)

closeArray("vapour",
           vapourVal_temperatureValues,
           vapourVal_densityValues,
           vapourVal_viscosityValues,
           vapourVal_isobaricThermalCapacityValues,
           vapourVal_thermalConductivityValues,
           vapourVal_specificEnthalpyValues,
           vapourVal_specificEntropyValues,
           vapourVal_linearExpansionCoefficientValues)

with open(path+ "_liquid.csv", "w") as filePort:
    filePort.write("T, ro, mu, cp, lambda, std. enthalpy, std. entropy, lexp. coeff, Psat\n")
    for val in zip(liquidVal_temperatureValues,
                   liquidVal_densityValues,
                   liquidVal_viscosityValues,
                   liquidVal_isobaricThermalCapacityValues,
                   liquidVal_thermalConductivityValues,
                   liquidVal_specificEnthalpyValues,
                   liquidVal_specificEntropyValues,
                   liquidVal_linearExpansionCoefficientValues):
        if val[0] == 1.0: filePort.write("00{0:.6F}, {1:.9F}, {2:.9F}, {3:.9F}, {4:.9F}, {5:.9F}, {6:.9F}, {7:.9F}\n".format(*val))
        else: filePort.write("{0:.6F}, {1:.9F}, {2:.9F}, {3:.9F}, {4:.9F}, {5:.9F}, {6:.9F}, {7:.9F}\n".format(*val))
              
with open(path+ "_vapour.csv", "w") as filePort:
    filePort.write("T, ro, mu, cp, lambda, std. enthalpy, std. entropy, lexp. coeff\n")
    for val in zip(vapourVal_temperatureValues,
                   vapourVal_densityValues,
                   vapourVal_viscosityValues,
                   vapourVal_isobaricThermalCapacityValues,
                   vapourVal_thermalConductivityValues,
                   vapourVal_specificEnthalpyValues,
                   vapourVal_specificEntropyValues,
                   vapourVal_linearExpansionCoefficientValues):
        if val[0] == 1.0: filePort.write("00{0:.6F}, {1:.9F}, {2:.9F}, {3:.9F}, {4:.9F}, {5:.9F}, {6:.9F}, {7:.9F}\n".format(*val))
        else: filePort.write("{0:.6F}, {1:.9F}, {2:.9F}, {3:.9F}, {4:.9F}, {5:.9F}, {6:.9F}, {7:.9F}\n".format(*val))
        
with open(path+ ".scm", "w") as filePort:
    filePort.write(0 * "\t" + "({0}\n".format(""))
    filePort.write(fluentIAPWSgen(temperatureValues,
                                densityValues,
                                viscosityValues,
                                isobaricThermalCapacityValues,
                                thermalConductivityValues,
                                specificEnthalpyValues,
                                name="water-if97",
                                reference_temperature=vapourVal_temperatureValues[1],
                                standard_state_enthalpy=0))

    filePort.write(fluentIAPWSgen(vapourVal_temperatureValues,
                                vapourVal_densityValues,
                                vapourVal_viscosityValues,
                                vapourVal_isobaricThermalCapacityValues,
                                vapourVal_thermalConductivityValues,
                                vapourVal_specificEnthalpyValues,
                                name="water-vapour-if97",
                                reference_temperature=vapourVal_temperatureValues[1],
                                standard_state_enthalpy=vapourVal_specificEnthalpyValues[1]-liquidVal_specificEnthalpyValues[-2]))

    filePort.write(fluentIAPWSgen(liquidVal_temperatureValues,
                                liquidVal_densityValues,
                                liquidVal_viscosityValues,
                                liquidVal_isobaricThermalCapacityValues,
                                liquidVal_thermalConductivityValues,
                                liquidVal_specificEnthalpyValues,
                                name="water-liquid-if97",
                                reference_temperature=vapourVal_temperatureValues[1],
                                standard_state_enthalpy=0))
    filePort.write(0*"\t"+")\n")

with open(path+ "_saturationCurve.csv", "w") as filePort:
    additionTemperatures        = linspace(274.15, temperatureValues[-2], 49)
    #states                      = (IAPWS97(T=temp, P = getSaturationPressure_uncorrected(temp)) for temp in additionTemperatures)
    states                      = (IAPWS97(T=temp, x=0) for temp in additionTemperatures)
    
    base = "//define/phase/interaction-domain , , , , , , , , , yes liquid vapour evaporation-condensation no yes 50 "
    base_satur = base[:]
    base_sigma = "no yes evaporation-condensation 50 "
    filePort.write("Tsat, Psat, Sigma\n")
    
    saturTemperatures   = list()
    saturPressures      = list()
    sigmas              = list()
    
    for thermodynamicState in states:
        Tsat    = thermodynamicState.T-273.15
        Psat    = thermodynamicState.P*10
        sigma   = thermodynamicState.sigma
        saturTemperatures.append(Tsat)
        saturPressures.append(Psat)
        sigmas.append(sigma)
        #print("T: {0}\tP:{1}\tS:{2}".format(Tsat, Psat, sigma))
        base_satur += " {0} {1} ".format(Psat, Tsat)
        base_sigma += " {0} {1} ".format(Tsat, sigma)
        filePort.write("{0}, {1}, {2}\n".format(Tsat, Psat, sigma))
        
    Tcrit       = 647.096 - 273.15
    Pcrit       = 220.64
    Sigmacrit   = 0 #Eotvos rule
    saturTemperatures.append(Tcrit)
    saturPressures.append(Pcrit)
    sigmas.append(Sigmacrit)
    
    base_satur += " {0} {1} ".format(Pcrit, Tcrit)
    base_sigma += " {0} {1} ".format(Tcrit, Sigmacrit)
    filePort.write("{0}, {1}, {2}\n".format(Tcrit, Pcrit, Sigmacrit))
    saturCurveSequence = base_satur + base_sigma + " , "
    with open(path + "_saturationCurve_fluentSequence.FLUSEQ", "w") as subFilePort:
        subFilePort.write(base_satur)
        subFilePort.write(base_sigma)
         
if plotting: 
    import matplotlib.pyplot as plt
    my_dpi = 250
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=False, layout="constrained")
    fig.set_size_inches(4*1.618,4*1.000)
    fig.suptitle("Thermohydraulic properties for Pabs = {0} MPa".format(pressureAbsolute))

    axs[0, 0].plot(temperatureValues, densityValues, "black")
    axs[0, 0].set_title("Density")
    axs[0, 0].set_xlabel("T [K]")
    axs[0, 0].set_ylabel(r"$\rho$ [kg m$^{-3}$]")
    axs[0, 0].grid()

    axs[1, 0].plot(temperatureValues, viscosityValues, "blue")
    axs[1, 0].set_title("Dynamic Viscosity")
    axs[1, 0].set_xlabel("T [K]")
    axs[1, 0].set_ylabel(r"$\mu$ [Pa s]")
    axs[1, 0].grid()

    axs[0, 1].plot(temperatureValues, isobaricThermalCapacityValues, "red")
    axs[0, 1].set_title("Isobaric Thermal Capacity")
    axs[0, 1].set_xlabel("T [K]")
    axs[0, 1].set_ylabel(r"cp [J kg$^{-1}$ K$^{-1}$]")
    axs[0, 1].grid()

    axs[1, 1].plot(temperatureValues, thermalConductivityValues, "green")
    axs[1, 1].set_title("Thermal Conductivity")
    axs[1, 1].set_xlabel("T [K]")
    axs[1, 1].set_ylabel(r"$\lambda$ [W m$^{-1}$ K$^{-1}$]")
    axs[1, 1].grid()
    plt.savefig(path+ ".png", dpi=my_dpi)
    
    
    
    x_shift = 10
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=False, layout="constrained")
    fig.set_size_inches(4*1.618,4*1.000)
    fig.suptitle("Thermohydraulic properties for Pabs = {0} MPa [Liq]".format(pressureAbsolute))

    axs[0, 0].plot(liquidVal_temperatureValues, liquidVal_densityValues, "black")
    axs[0, 0].set_title("Density")
    axs[0, 0].set_xlabel("T [K]")
    axs[0, 0].set_xbound(
                        lower=liquidVal_temperatureValues[1] - x_shift,
                        upper=liquidVal_temperatureValues[-2] + x_shift
                        )                        
    axs[0, 0].set_ylabel(r"$\rho$ [kg m$^{-3}$]")
    axs[0, 0].grid()

    axs[1, 0].plot(liquidVal_temperatureValues, liquidVal_viscosityValues, "blue")
    axs[1, 0].set_title("Dynamic Viscosity")
    axs[1, 0].set_xlabel("T [K]")
    axs[1, 0].set_xbound(
                        lower=liquidVal_temperatureValues[1] - x_shift,
                        upper=liquidVal_temperatureValues[-2] + x_shift
                        )                    
    axs[1, 0].set_ylabel(r"$\mu$ [Pa s]")
    axs[1, 0].grid()

    axs[0, 1].plot(liquidVal_temperatureValues, liquidVal_isobaricThermalCapacityValues, "red")
    axs[0, 1].set_title("Isobaric Thermal Capacity")
    axs[0, 1].set_xlabel("T [K]")
    axs[0, 1].set_xbound(
                        lower=liquidVal_temperatureValues[1] - x_shift,
                        upper=liquidVal_temperatureValues[-2] + x_shift
                        )                           
    axs[0, 1].set_ylabel(r"cp [J kg$^{-1}$ K$^{-1}$]")
    axs[0, 1].grid()

    axs[1, 1].plot(liquidVal_temperatureValues, liquidVal_thermalConductivityValues, "green")
    axs[1, 1].set_title("Thermal Conductivity")
    axs[1, 1].set_xlabel("T [K]")
    axs[1, 1].set_xbound(
                        lower=liquidVal_temperatureValues[1] - x_shift,
                        upper=liquidVal_temperatureValues[-2] + x_shift
                        )                            
    axs[1, 1].set_ylabel(r"$\lambda$ [W m$^{-1}$ K$^{-1}$]")
    axs[1, 1].grid()
    plt.savefig(path+ "_liquid.png", dpi=my_dpi)
    
    
    
    
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=False, layout="constrained")
    fig.set_size_inches(4*1.618,4*1.000)
    fig.suptitle("Thermohydraulic properties for Pabs = {0} MPa [Vap]".format(pressureAbsolute))

    axs[0, 0].plot(vapourVal_temperatureValues, vapourVal_densityValues, "black")
    axs[0, 0].set_title("Density")
    axs[0, 0].set_xlabel("T [K]")
    axs[0, 0].set_xbound(
                        lower=vapourVal_temperatureValues[1] - x_shift,
                        upper=vapourVal_temperatureValues[-2] + x_shift
                        )                            
    axs[0, 0].set_ylabel(r"$\rho$ [kg m$^{-3}$]")
    axs[0, 0].grid()

    axs[1, 0].plot(vapourVal_temperatureValues, vapourVal_viscosityValues, "blue")
    axs[1, 0].set_title("Dynamic Viscosity")
    axs[1, 0].set_xlabel("T [K]")
    axs[1, 0].set_xbound(
                        lower=vapourVal_temperatureValues[1] - x_shift,
                        upper=vapourVal_temperatureValues[-2] + x_shift
                        )                         
    axs[1, 0].set_ylabel(r"$\mu$ [Pa s]")
    axs[1, 0].grid()

    axs[0, 1].plot(vapourVal_temperatureValues, vapourVal_isobaricThermalCapacityValues, "red")
    axs[0, 1].set_title("Isobaric Thermal Capacity")
    axs[0, 1].set_xlabel("T [K]")
    axs[0, 1].set_xbound(
                        lower=vapourVal_temperatureValues[1] - x_shift,
                        upper=vapourVal_temperatureValues[-2] + x_shift
                        )                          
    axs[0, 1].set_ylabel(r"cp [J kg$^{-1}$ K$^{-1}$]")
    axs[0, 1].grid()

    axs[1, 1].plot(vapourVal_temperatureValues, vapourVal_thermalConductivityValues, "green")
    axs[1, 1].set_title("Thermal Conductivity")
    axs[1, 1].set_xlabel("T [K]")
    axs[1, 1].set_xbound(
                        lower=vapourVal_temperatureValues[1] - x_shift,
                        upper=vapourVal_temperatureValues[-2] + x_shift
                        )                          
    axs[1, 1].set_ylabel(r"$\lambda$ [W m$^{-1}$ K$^{-1}$]")
    axs[1, 1].grid()
    plt.savefig(path+ "_vapour.png", dpi=my_dpi)
    plt.cla()
    plt.clf()
    
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=False, layout="constrained")
    fig.set_size_inches(4*1.618,4*1.000)
    fig.suptitle("Saturation Properties")
    axs[0].scatter(saturTemperatures, saturPressures, c="black", s=5)
    axs[0].grid()
    axs[0].set_title("Saturation Curve")
    axs[0].set_xlabel("T [°C]")
    axs[0].set_ylabel("P [bar(a)]")
    axs[1].scatter(saturTemperatures, sigmas, c="black", s=5)
    axs[1].grid()
    axs[1].set_title("Surface Tension l->g")
    axs[1].set_xlabel("T [°C]")
    axs[1].set_ylabel(r"$\sigma$ [N m$^{-1}$]")
    plt.savefig(path+ "_saturationCurve.png", dpi=my_dpi)
    
    #plt.show()
print()
print("Saturation temperature:\t\t\t{:.6F} [K]\t{:.6F} [°C]".format(saturationTemperature, saturationTemperature-273.15))
print("Reference temperature:\t\t\t{:.6F} [K]\t{:.6F} [°C]".format(vapourVal_temperatureValues[1], vapourVal_temperatureValues[1]-273.15))
print("Relative enthalpy [liquid]:\t{:.0F} [J kg-1 mol-1]".format(liquidVal_specificEnthalpyValues[-1]))
print("Relative enthalpy [vapour]:\t{:.0F} [J kg-1 mol-1]".format(vapourVal_specificEnthalpyValues[0]))
print("Latent heat of evaporation:\t\t{:.0F} [J kg-1 mol-1]".format(vapourVal_specificEnthalpyValues[0]-liquidVal_specificEnthalpyValues[-1]))
print()
print("Saturation Curve Sequence Fluent:\n")
print(saturCurveSequence)
print()    
print("Uspesne dokonceno!")
print("_______________________________")
print()
print()

exit()

