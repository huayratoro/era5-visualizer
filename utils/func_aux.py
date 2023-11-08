from metpy.calc import dewpoint_from_specific_humidity, equivalent_potential_temperature
from metpy.units import units

## calculo de titae
def titaeEra5(level, temp, q):
    dp = dewpoint_from_specific_humidity(level * units.hPa, (temp-273.15) * units.degC, (q*1000) * units('g/kg'))
    titae = equivalent_potential_temperature(level * units.hPa, (temp-273.15) * units.degC, dp)
    return(titae)
