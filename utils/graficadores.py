def graficadorTitae850Viento(year, month, day, hour, pSalida):

    import xarray as xr
    import matplotlib.pyplot as plt
    from metpy.calc import dewpoint_from_specific_humidity, equivalent_potential_temperature
    from metpy.units import units
    import numpy as np
    import cartopy.feature as cfeature
    import cartopy.crs as ccrs
    import matplotlib as mpl

    countries = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_0_countries',
            scale='10m',
            facecolor='none',edgecolor='black')

    data = xr.open_dataset(f"db/era-5-pressure-level/era-5-pressure-levels-{year}-{month}-{day}-{hour}.grib")

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

    # graf
    transf = transform=ccrs.PlateCarree()
    fig=plt.figure(figsize=(10,10))
    ax = plt.subplot(1, 1, 1, projection=transf)
    ax.set_extent([-90, -30, -50, 10], crs=transf)

    # Agregamos los límites de los países
    ax.add_feature(countries,linewidth=0.4)

    ## TPE
    titaeEra5(850, data['t'].sel(isobaricInhPa=850), data['q'].sel(isobaricInhPa=850)).plot(ax = ax, transform=transform, vmin = 270, vmax = 360, cmap = cmapTitae, extend = "both")

    ## Viento
    ax.barbs(data["longitude"].values, data["latitude"].values, data['u'].sel(isobaricInhPa=850).values*1.94, data['v'].sel(isobaricInhPa=850).values*1.94, length=4.5, sizes=dict(emptybarb=0, spacing=0.1, height=0.5), linewidth=0.75, transform=transform, regrid_shape=25, alpha=0.5)

    plt.savefig(f"{pSalida}{year}-{month}-{day}-{hour}-titae-850-viento.png", dpi = 500, bbox_inches = 'tight')
    plt.close("all")









