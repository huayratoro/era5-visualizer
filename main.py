from utils.download_era5 import downloadEra5PressureLevels, downloadEra5SingleLevel
from utils.graficadores import graficadorTitae850Viento
from sys import argv

argv = argv[1::]

if argv[0] == "--download-era5-pressure-levels":
    downloadEra5PressureLevels(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9])
if argv[0] == "--download-era5-single-level":
    downloadEra5SingleLevel(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9])
if argv[0] == "--graf-850-titae-viento":
    graficadorTitae850Viento(argv[1], argv[2], argv[3], argv[4], argv[5])
#
else:
    print(f"Argumento {argv[0]} inválido")
    exit()