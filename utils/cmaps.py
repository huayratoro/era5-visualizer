import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#---------------------------------
cmapTitae = mpl.colors.ListedColormap(
    np.vstack((
        plt.get_cmap("Blues_r")(np.linspace(0.1, 0.75, 12)),
        plt.get_cmap("Greens_r")(np.linspace(0.25, 0.75, 16)),
        plt.get_cmap("Reds")(np.linspace(0.3, 1, 8))
    ))
)
#---------------------------------
cmapPw = mpl.colors.ListedColormap(
    np.vstack((
        plt.get_cmap("Reds_r")(np.linspace(0.5, 0.9, 8)),
        plt.get_cmap("Greens_r")(np.linspace(0.5, 0.9, 8)),
        plt.get_cmap("Blues")(np.linspace(0.5, 0.9, 8)),
    ))
)
#---------------------------------
cmapIvt = mpl.colors.ListedColormap(
    [
        "#ffffff", "#fedb01", "#fda602", "#ed7a05", "#a12424", "#660d08", "#941ddb"
    ]
)
#---------------------------------
