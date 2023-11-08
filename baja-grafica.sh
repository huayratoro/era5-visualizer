#### Primero baja los casos

python main.py --download-era5-lista-fechas db/lista-fechas.csv 10 -50 -90 -30 db/
echo "Campos bajados"

#### Grafica cada campo

# python main.py --graf-850-titae-viento -l db/lista-fechas.csv assets/
# python main.py --graf-500-hgt -l db/lista-fechas.csv assets/
# python main.py --graf-pw -l db/lista-fechas.csv assets/
# python main.py --graf-hgt-1000-ivt -l db/lista-fechas.csv assets/
# echo "Ejecucion completa"