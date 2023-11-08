import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import metpy.calc as mpcalc
from metpy.calc import dewpoint_from_specific_humidity, equivalent_potential_temperature
from metpy.units import units
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl
from scipy.ndimage import gaussian_filter
from utils.cmaps import *
from utils.func_aux import *

countries = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_countries',
        scale='10m',
        facecolor='none',edgecolor='gray')

#----------------------------------------
def graficadorTitae850Viento(year, month, day, hour, pSalida):

    data = xr.open_dataset(f"db/era-5-pressure-levels/era-5-pressure-levels-{year}-{month}-{day}-{hour}.nc")

    # parametros graficos
    transform=ccrs.PlateCarree()
    fig=plt.figure(figsize=(10,10))
    ax = plt.subplot(1, 1, 1, projection=transform)
    ax.set_extent([-90, -30, -50, 10], crs=transform)

    # Agregamos los límites de los países
    ax.add_feature(countries,linewidth=0.4)

    ## TPE
    cm = titaeEra5(850, data['t'].sel(level=850), data['q'].sel(level=850)).plot(ax = ax, transform=transform, vmin = 270, vmax = 360, cmap = cmapTitae, add_colorbar = False)
    cb = plt.colorbar(cm, ax = ax, shrink = 0.75, extend = "both")
    cb.set_label("[Kelvin]")

    ## Viento
    ax.barbs(data["longitude"].values, data["latitude"].values, data['u'].sel(level=850).values[0, :, :]*1.94, data['v'].sel(level=850).values[0, :, :]*1.94, length=4.5, sizes=dict(emptybarb=0, spacing=0.1, height=0.5), linewidth=0.75, transform=transform, regrid_shape=25, alpha=0.5)

    ## Para los ticks latitudinales y longitudinales
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.75, color='gray', alpha=0.75, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ## Titulo
    ax.set_title(f"{year}-{month}-{day} {hour}:00 Z \n titae [shaded], and winds [knots] at 850 mb")

    plt.savefig(f"{pSalida}{year}-{month}-{day}-{hour}-titae-850-viento.png", dpi = 500, bbox_inches = 'tight')
    plt.close("all")
#----------------------------------------
def graficador500(year, month, day, hour, pSalida):

    data = xr.open_dataset(f"db/era-5-pressure-levels/era-5-pressure-levels-{year}-{month}-{day}-{hour}.nc")

    # parametros graficos
    transf = transform=ccrs.PlateCarree()
    fig=plt.figure(figsize=(10,10))
    ax = plt.subplot(1, 1, 1, projection=transf)
    ax.set_extent([-90, -30, -50, 10], crs=transf)

    # Agregamos los límites de los países
    ax.add_feature(countries,linewidth=0.4)

    # Mapa 500 mb
    clevs_500_hght = np.arange(400, 603, 3)
    damHgt500 = gaussian_filter(data['z'].sel(level=500).values[0, :, :] / 98.0665, sigma=1.5)
    cs = ax.contour(data['longitude'].values, data['latitude'].values, damHgt500, clevs_500_hght, colors='black', transform=transf, linewidths = 0.5)
    plt.clabel(cs, fmt='%d')
    # Velocidad de viento en 500mb
    uwnd_500 = gaussian_filter(data['u'].sel(level=500).values[0, :, :], sigma=2.0) * units('m/s')
    vwnd_500 = gaussian_filter(data['v'].sel(level=500).values[0, :, :], sigma=2.0) * units('m/s')
    sped_500 = mpcalc.wind_speed(uwnd_500, vwnd_500).to('kt')
    # Velocidad del viento
    clevs_500_sped = np.arange(30, 140, 20)
    cf = ax.contourf(data['longitude'].values, data['latitude'].values, sped_500, clevs_500_sped, cmap=plt.cm.Reds, transform=transf, extend = "max")
    plt.colorbar(cf, ax = ax, shrink = 0.75)    
    
    ## Para los ticks latitudinales y longitudinales
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.75, color='gray', alpha=0.75, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ## Titulo
    ax.set_title(f"{year}-{month}-{day} {hour}:00 Z \n Geopot. Height [dam, contour] and wind speed [knots] at 500 mb")

    plt.savefig(f"{pSalida}{year}-{month}-{day}-{hour}-500-hgt.png", dpi = 500, bbox_inches = 'tight')
    plt.close("all")
#----------------------------------------
def graficadorPwHgt1000(year, month, day, hour, pSalida):

    data = xr.open_dataset(f"db/era-5-single-level/era-5-single-level-{year}-{month}-{day}-{hour}.nc")

    # parametros graficos
    transf = ccrs.PlateCarree()
    fig=plt.figure(figsize=(10,10))
    ax = plt.subplot(1, 1, 1, projection=transf)
    ax.set_extent([-90, -30, -50, 10], crs=transf)

    # Agregamos los límites de los países
    ax.add_feature(countries,linewidth=0.4)

    ## PW
    cm = data['tcwv'].plot(ax = ax, transform=transf, vmin = 0, vmax = 60, cmap = cmapPw, add_colorbar = False)
    cb = plt.colorbar(cm, ax = ax, shrink = 0.75)
    cb.set_label("[mm]")
    
    ## En contornos la PW mas alta
    pw = gaussian_filter(data['tcwv'].values[0, :, :], sigma = 3)
    ax.contour(data['longitude'].values, data['latitude'].values, pw, np.arange(40, 70, 10), colors = ['magenta'], linewidths = 0.5, transform = transf)
    
    ## Vector IVT
    ax.quiver(
        data["longitude"].values, 
        data["latitude"].values, 
        data['p71.162'].values[0, :, :], 
        data['p72.162'].values[0, :, :], 
        transform=transf, regrid_shape=30, alpha=0.5)

    ## Para los ticks latitudinales y longitudinales
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.75, color='gray', alpha=0.75, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ## Titulo
    ax.set_title(f"{year}-{month}-{day} {hour}:00 Z \n Precipitable Water [shaded] & IVT [vectors, kg m-1 s-1]")

    plt.savefig(f"{pSalida}{year}-{month}-{day}-{hour}-pw.png", dpi = 500, bbox_inches = 'tight')
    plt.close("all")
#----------------------------------------
def graficadorIvtHgt1000(year, month, day, hour, pSalida):

    sfc = xr.open_dataset(f"db/era-5-pressure-levels/era-5-pressure-levels-{year}-{month}-{day}-{hour}.nc")
    ivt = xr.open_dataset(f"db/era-5-single-level/era-5-single-level-{year}-{month}-{day}-{hour}.nc")
    
    # parametros graficos
    transform=ccrs.PlateCarree()
    fig=plt.figure(figsize=(10,10))
    ax = plt.subplot(1, 1, 1, projection=transform)
    ax.set_extent([-90, -30, -50, 10], crs=transform)

    # Agregamos los límites de los países
    ax.add_feature(countries,linewidth=0.4)

    # IVT sombreado
    ivtMag = np.sqrt( ivt['p71.162'].values[0, :, :]**2 + ivt['p72.162'].values[0, :, :]**2 )
    cm = ax.pcolormesh(ivt['longitude'].values, ivt['latitude'].values, ivtMag, vmin = 0, vmax = 1050, cmap = cmapIvt)
    cb = plt.colorbar(cm, ax = ax, shrink = 0.75, ticks=np.arange(0, 1050+150, 150), extend = 'max')
    cb.set_label("[kg m-1 s-1]")
    # Vector IVT
    ax.quiver(
        ivt["longitude"].values, 
        ivt["latitude"].values, 
        ivt['p71.162'].values[0, :, :], 
        ivt['p72.162'].values[0, :, :], 
        transform=transform, regrid_shape=30, alpha=0.5)
    # hgt 1000 mb
    clevs_1000_hght = np.arange(-10, 32, 2)
    damHgt1000 = gaussian_filter(sfc['z'].sel(level=1000).values[0, :, :] / 98.0665, sigma=2.5)
    cs = ax.contour(sfc['longitude'].values, sfc['latitude'].values, damHgt1000, clevs_1000_hght, colors='black', transform=transform, linewidths = 0.5)
    plt.clabel(cs, fmt='%d')

    ## Para los ticks latitudinales y longitudinales
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.75, color='gray', alpha=0.75, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ## Titulo
    ax.set_title(f"{year}-{month}-{day} {hour}:00 Z \n IVT [shaded], IVT [vectors] & HGT 1000 mb [contours, dam]")

    plt.savefig(f"{pSalida}{year}-{month}-{day}-{hour}-ivt-1000.png", dpi = 500, bbox_inches = 'tight')
    plt.close("all")
#----------------------------------------
def graficadorTresCampos(year, month, day, hour, lonC, latC, pSalida):
    # Pressure Levels
    pressureLevels = xr.open_dataset(f"db/era-5-pressure-levels/era-5-pressure-levels-{year}-{month}-{day}-{hour}.nc")
    
    # Single Level
    singleLevel = xr.open_dataset(f"db/era-5-single-level/era-5-single-level-{year}-{month}-{day}-{hour}.nc")

    # parametros graficos generales
    alphaVectores = 0.8
    anchoBarraColores = 1; padBarraColores = 0.075
    markerRpf = "s"; markerSize = 10; markerColor = "black"; markerWidths = 0.5; markerFaceColor = "none"
    transf = transform=ccrs.PlateCarree()
    cmts = 1/2.54
    plt.rcParams.update({"font.size" : 5})
    plt.rcParams["font.family"] = "Arial"
    fig=plt.figure(figsize=(19*cmts, 9*cmts))
    
    ##---------- 500 mb ----------##
    
    ax1 = plt.subplot(1, 3, 1, projection=transf)
    ax1.set_extent([-90, -30, -50, 10], crs=transf)

    # Agregamos los límites de los países
    ax1.add_feature(countries,linewidth=0.4)

    # Mapa 500 mb
    clevs_500_hght = np.arange(400, 603, 3)
    damHgt500 = gaussian_filter(pressureLevels['z'].sel(level=500).values[0, :, :] / 98.0665, sigma=1.5)
    cs = ax1.contour(pressureLevels['longitude'].values, pressureLevels['latitude'].values, damHgt500, clevs_500_hght, colors='black', transform=transf, linewidths = 0.5)
    plt.clabel(cs, fmt='%d', fontsize = 4)
    # Velocidad de viento en 500mb
    uwnd_500 = gaussian_filter(pressureLevels['u'].sel(level=500).values[0, :, :], sigma=2.0) * units('m/s')
    vwnd_500 = gaussian_filter(pressureLevels['v'].sel(level=500).values[0, :, :], sigma=2.0) * units('m/s')
    sped_500 = mpcalc.wind_speed(uwnd_500, vwnd_500).to('kt')
    # Velocidad del viento
    clevs_500_sped = np.arange(30, 140, 20)
    cf = ax1.contourf(pressureLevels['longitude'].values, pressureLevels['latitude'].values, sped_500, clevs_500_sped, cmap=plt.cm.Reds, transform=transf, extend = "max")
    cb = plt.colorbar(cf, ax = ax1, shrink = anchoBarraColores, orientation='horizontal', pad = padBarraColores)    
    cb.set_label("[knots]")
    # localizacion del centroide
    ax1.scatter(lonC, latC, s = markerSize, marker = markerRpf, color = markerColor, linewidths=markerWidths, facecolors = markerFaceColor,transform = transf)
    
    ## Para los ticks latitudinales y longitudinales
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.75, color='gray', alpha=0.75, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ## Titulo
    ax1.set_title("Geopot. Height 500 mb [dam, contour] \n wind speed [shaded]")

    ##---------- 850 mb ----------##

    ax2 = plt.subplot(1, 3, 2, projection=transform)
    ax2.set_extent([-90, -30, -50, 10], crs=transform)

    # Agregamos los límites de los países
    ax2.add_feature(countries,linewidth=0.4)

    ## TPE
    cm = titaeEra5(850, pressureLevels['t'].sel(level=850), pressureLevels['q'].sel(level=850)).plot(ax = ax2, transform=transform, vmin = 270, vmax = 360, cmap = cmapTitae, add_colorbar = False)
    cb = plt.colorbar(cm, ax = ax2, shrink = anchoBarraColores, extend = "both", orientation='horizontal', pad = padBarraColores, ticks=np.arange(270, 370, 10))
    cb.set_label("[Kelvin]")

    ## Viento
    ax2.quiver(
        pressureLevels["longitude"].values,
        pressureLevels["latitude"].values,
        pressureLevels['u'].sel(level=850).values[0, :, :]*1.94,
        pressureLevels['v'].sel(level=850).values[0, :, :]*1.94,
        transform=transf, regrid_shape=30, alpha=alphaVectores
    )
    # localizacion del centroide
    ax2.scatter(lonC, latC, s = markerSize, marker = markerRpf, color = markerColor, linewidths=markerWidths, facecolors = markerFaceColor, transform = transf)
    
    ## Para los ticks latitudinales y longitudinales
    gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.75, color='gray', alpha=0.75, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.xformatter = LONGITUDE_FORMATTER

    ## Titulo
    ax2.set_title("Titae [shaded] 850 mb \n winds [knots]")

    ##---------- 1000 mb ----------##

    ax3 = plt.subplot(1, 3, 3, projection=transf)
    ax3.set_extent([-90, -30, -50, 10], crs=transf)

    # Agregamos los límites de los países
    ax3.add_feature(countries,linewidth=0.4)

    ## PW
    cm = singleLevel['tcwv'].plot(ax = ax3, transform=transf, vmin = 0, vmax = 60, cmap = cmapPw, add_colorbar = False)
    cb = plt.colorbar(cm, ax = ax3, shrink = anchoBarraColores, orientation='horizontal', pad = padBarraColores, extend = 'max')
    cb.set_label("[mm]")
    
    ## En contornos la PW mas alta
    pw = gaussian_filter(singleLevel['tcwv'].values[0, :, :], sigma = 3)
    ax3.contour(singleLevel['longitude'].values, singleLevel['latitude'].values, pw, np.arange(40, 70, 10), colors = ['cyan', 'white', "magenta"], linewidths = 0.75, transform = transf)
    
    ## Vector IVT
    ax3.quiver(
        singleLevel["longitude"].values, 
        singleLevel["latitude"].values, 
        singleLevel['p71.162'].values[0, :, :], 
        singleLevel['p72.162'].values[0, :, :], 
        transform=transf, regrid_shape=30, alpha=alphaVectores)
    # localizacion del centroide
    ax3.scatter(lonC, latC, s = markerSize, marker = markerRpf, color = markerColor, linewidths=markerWidths, facecolors = markerFaceColor, transform = transf)

    ## Para los ticks latitudinales y longitudinales
    gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.75, color='gray', alpha=0.75, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.xformatter = LONGITUDE_FORMATTER

    ## Titulo
    ax3.set_title("Precipitable Water [shaded] \n IVT [vectors, kg m-1 s-1]")

    ####
    plt.subplots_adjust(wspace=0.075)
    fig.text(0.075, 0.4, f"{year}-{month}-{day} {hour}:00 Z", rotation = 90, fontweight = 'bold')
    
    plt.savefig(f"{pSalida}{year}-{month}-{day}-{hour}-campos-500-850-1000.png", dpi = 500, bbox_inches = 'tight')
    plt.close("all")

