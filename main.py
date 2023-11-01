import pandas as pd
from utils.download_era5 import *
from utils.graficadores import *
from sys import argv

argv = argv[1::]

#----------------------------------------
if argv[0] == "--download-era5-pressure-levels":
    downloadEra5PressureLevels(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9])
#----------------------------------------
if argv[0] == "--download-era5-single-level":
    downloadEra5SingleLevel(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9])
#----------------------------------------
if (argv[0] == "--download-era5-lista-fechas"):
        downloadListDates(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6])
#----------------------------------------
if (argv[0] == "--graf-850-titae-viento"):
    if argv[1] == "-l":
        fechasEventos = pd.read_csv(argv[2])
        fechasEventos.fechasADescargar = pd.to_datetime(fechasEventos.fechasADescargar)
        for fecha in fechasEventos.fechasADescargar:
            graficadorTitae850Viento(fecha.year, "%02d"%fecha.month, "%02d"%fecha.day, "%02d"%fecha.hour, argv[3])
    else:
        graficadorTitae850Viento(argv[1], argv[2], argv[3], argv[4], argv[5])
#----------------------------------------
if (argv[0] == "--graf-500-hgt"):
    if argv[1] == "-l":
        fechasEventos = pd.read_csv(argv[2])
        fechasEventos.fechasADescargar = pd.to_datetime(fechasEventos.fechasADescargar)
        for fecha in fechasEventos.fechasADescargar:
            graficador500(fecha.year, "%02d"%fecha.month, "%02d"%fecha.day, "%02d"%fecha.hour, argv[3])
    else:
        graficador500(argv[1], argv[2], argv[3], argv[4], argv[5])
#----------------------------------------
if (argv[0] == "--graf-pw"):
    if argv[1] == "-l":
        fechasEventos = pd.read_csv(argv[2])
        fechasEventos.fechasADescargar = pd.to_datetime(fechasEventos.fechasADescargar)
        for fecha in fechasEventos.fechasADescargar:
            graficadorPwHgt1000(fecha.year, "%02d"%fecha.month, "%02d"%fecha.day, "%02d"%fecha.hour, argv[3])
    else:
        graficadorPwHgt1000(argv[1], argv[2], argv[3], argv[4], argv[5])
#----------------------------------------