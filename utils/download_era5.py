def downloadEra5PressureLevels(year, month, day, hour, latNorth, latSouth, lonWest, lonEast, pSalida):

    import cdsapi

    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': [
                'geopotential', 'specific_humidity', 'temperature',
                'u_component_of_wind', 'v_component_of_wind',
            ],
            'pressure_level': [
                '200', '500', '850', '1000',
            ],
            'year': str(year),
            'month': str(month),
            'day': str(day),
            'time': str(hour) + ':00',
            'area': [
                int(latNorth), int(lonWest), int(latSouth),
                int(lonEast),
            ],
        },
        f"{pSalida}era-5-pressure-levels-{year}-{month}-{day}-{hour}.nc")

def downloadEra5SingleLevel(year, month, day, hour, latNorth, latSouth, lonWest, lonEast, pSalida):

    import cdsapi

    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': ["total_column_water_vapour", "vertical_integral_of_eastward_water_vapour_flux", "vertical_integral_of_northward_water_vapour_flux"],
            'year': year,
            'month': month,
            'day': day,
            'time': str(hour) + ':00',
            'area': [
                int(latNorth), int(lonWest), int(latSouth),
                int(lonEast),
            ],
        },
        f"{pSalida}era-5-single-level-{year}-{month}-{day}-{hour}.nc")

def downloadListDates(pListaFechas, latNorth, latSouth, lonWest, lonEast, pSalida):
    
    import pandas as pd
    import numpy as np
    
    # Preparo las fechas
    fechasEventos = pd.read_csv(pListaFechas)
    fechasEventos.fechasADescargar = pd.to_datetime(fechasEventos.fechasADescargar)

    for fecha in fechasEventos.fechasADescargar:
        downloadEra5PressureLevels(
            fecha.year,
            "%02d"%fecha.month,
            "%02d"%fecha.day,
            "%02d"%fecha.hour,
            latNorth, latSouth, lonWest, lonEast, pSalida + "pressure-levels/"
        )
        downloadEra5SingleLevel(
            fecha.year,
            "%02d"%fecha.month,
            "%02d"%fecha.day,
            "%02d"%fecha.hour,
            latNorth, latSouth, lonWest, lonEast, pSalida + "single-level/"
        )






