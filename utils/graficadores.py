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

countries = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_countries',
        scale='10m',
        facecolor='none',edgecolor='gray')

def graficadorTitae850Viento(year, month, day, hour, pSalida):

    data = xr.open_dataset(f"db/era-5-pressure-level/era-5-pressure-levels-{year}-{month}-{day}-{hour}.nc")

    ## calculo de titae
    def titaeEra5(level, temp, q):
        dp = dewpoint_from_specific_humidity(level * units.hPa, (temp-273.15) * units.degC, (q*1000) * units('g/kg'))
        titae = equivalent_potential_temperature(level * units.hPa, (temp-273.15) * units.degC, dp)
        return(titae)
    cmapTitae = mpl.colors.ListedColormap(
        np.vstack((
            plt.get_cmap("Blues_r")(np.linspace(0.1, 0.75, 12)),
            plt.get_cmap("Greens_r")(np.linspace(0.25, 0.75, 16)),
            plt.get_cmap("Reds")(np.linspace(0.3, 1, 8))
        ))
    )

    # parametros graficos
    transf = transform=ccrs.PlateCarree()
    fig=plt.figure(figsize=(10,10))
    ax = plt.subplot(1, 1, 1, projection=transf)
    ax.set_extent([-90, -30, -50, 10], crs=transf)

    # Agregamos los límites de los países
    ax.add_feature(countries,linewidth=0.4)

    ## TPE
    cm = titaeEra5(850, data['t'].sel(isobaricInhPa=850), data['q'].sel(isobaricInhPa=850)).plot(ax = ax, transform=transform, vmin = 270, vmax = 360, cmap = cmapTitae, add_colorbar = False)
    cb = plt.colorbar(cm, ax = ax, shrink = 0.75, extend = "both")
    cb.set_label("[Kelvin]")

    ## Viento
    ax.barbs(data["longitude"].values, data["latitude"].values, data['u'].sel(isobaricInhPa=850).values*1.94, data['v'].sel(isobaricInhPa=850).values*1.94, length=4.5, sizes=dict(emptybarb=0, spacing=0.1, height=0.5), linewidth=0.75, transform=transform, regrid_shape=25, alpha=0.5)

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

def graficador500(year, month, day, hour, pSalida):

    data = xr.open_dataset(f"db/era-5-pressure-level/era-5-pressure-levels-{year}-{month}-{day}-{hour}.nc")

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

def graficadorPwHgt1000(year, month, day, hour, pSalida):

    data = xr.open_dataset(f"db/era-5-single-level/era-5-single-level-{year}-{month}-{day}-{hour}.nc")

    cmapPw = mpl.colors.ListedColormap(
        np.vstack((
            plt.get_cmap("Reds_r")(np.linspace(0.5, 0.9, 8)),
            plt.get_cmap("Greens_r")(np.linspace(0.5, 0.9, 8)),
            plt.get_cmap("Blues")(np.linspace(0.5, 0.9, 8)),
        ))
    )

    # parametros graficos
    transf = transform=ccrs.PlateCarree()
    fig=plt.figure(figsize=(10,10))
    ax = plt.subplot(1, 1, 1, projection=transf)
    ax.set_extent([-90, -30, -50, 10], crs=transf)

    # Agregamos los límites de los países
    ax.add_feature(countries,linewidth=0.4)

    ## TPE
    cm = data['tcwv'].plot(ax = ax, transform=transform, vmin = 0, vmax = 60, cmap = cmapPw, add_colorbar = False)
    cb = plt.colorbar(cm, ax = ax, shrink = 0.75)
    cb.set_label("[mm]")

    ## Para los ticks latitudinales y longitudinales
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.75, color='gray', alpha=0.75, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ## Titulo
    ax.set_title(f"{year}-{month}-{day} {hour}:00 Z \n Precipitable Water [shaded]")

    plt.savefig(f"{pSalida}{year}-{month}-{day}-{hour}-pw.png", dpi = 500, bbox_inches = 'tight')
    plt.close("all")

