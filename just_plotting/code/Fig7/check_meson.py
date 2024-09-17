import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

from matplotlib import rc, rcParams, cm, colors, colormaps
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm, colors
import numpy as np

import os


def plot_line(v, e, ti, tf, color):
    plt.gca().add_patch(
        plt.Rectangle(
            [ti - 0.2, v - e], tf - ti + 0.4, 2 * e, facecolor=color, alpha=0.4
        )
    )

plt.style.use("paperdraft.mplstyle")
plt.figure(figsize=(4.5, 4))

def calculate_mass_and_plot(ti, tf, x_min, x_max, y_min, y_max, Nsource, Nsink, noAPE = False):
    if noAPE == False:
        # Read data from file
        data = np.loadtxt(f'eff_mass_N{Nsource}_N{Nsink}.txt')
    else:
        data = np.loadtxt(f'eff_mass_N{Nsource}_N{Nsink}_noAPE.txt')
    x = data[:, 0]
    y = data[:, 1]
    y_err = data[:, 2]
    a = data[0, 3]
    b = data[0, 4]

    # Plotting
    plt.errorbar(x, y, yerr=y_err, fmt='o', label='Data', elinewidth=1.5, markersize=4.5, alpha=0.6)
    plt.grid(linestyle='--')
    plt.ylim(y_min, y_max)
    plt.xlim(x_min, x_max)
    plt.xlabel('$t/a$', fontsize=15)
    plt.ylabel('$am_{\\rm eff}$', fontsize=15)
    #plt.legend(loc='best')

    # Save figure
    #filename = f'../../../plots/N{Nsource}_N{Nsink}_2.pdf'
    #plt.savefig(filename, dpi=300, bbox_inches='tight')
    return a,b

# Example usage:
ti = 13
tf = 23
y_min = 0.39
y_max = 0.48
x_min = 1.5
x_max = 22.5
Nsource = 80
Nsink = 40
plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(4))
a, b = calculate_mass_and_plot(ti, tf, x_min, x_max, y_min, y_max, Nsource, Nsink)

plot_line(a, b, 10, 22.0, plt.gca().lines[-1].get_color())

#plt.show()
plt.savefig(f'../../../plots/N{Nsource}_N{Nsink}_2.pdf', dpi=300, bbox_inches='tight')

plt.figure()
plt.style.use("paperdraft.mplstyle")
plt.figure(figsize=(4.5, 4))

# Example usage:
ti = 13
tf = 23
y_min = 0.405
y_max = 0.45
x_min = 9.7
x_max = 23.5
Nsource = 20
Nsink = 20
plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(3))
plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(0.01))
a, b = calculate_mass_and_plot(ti, tf, x_min, x_max, y_min, y_max, Nsource, Nsink)

plot_line(a, b, 21, 23.0, plt.gca().lines[-1].get_color())

#plt.show()
plt.savefig(f'../../../plots/N{Nsource}_N{Nsink}_2.pdf', dpi=300, bbox_inches='tight')



plt.figure()
plt.style.use("paperdraft.mplstyle")
plt.figure(figsize=(4.5, 4))

# Example usage:
ti = 13
tf = 23
y_min = 0.375
y_max = 0.503
x_min = 0.0
x_max = 22.5
Nsource = 170
Nsink = 170
plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(4))
plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(0.02))
a, b = calculate_mass_and_plot(ti, tf, x_min, x_max, y_min, y_max, Nsource, Nsink)

plot_line(a, b, 10, 22.0, plt.gca().lines[-1].get_color())

#plt.show()
plt.savefig(f'../../../plots/N{Nsource}_N{Nsink}_2.pdf', dpi=300, bbox_inches='tight')



plt.figure()
plt.style.use("paperdraft.mplstyle")
plt.figure(figsize=(4.5, 4))

# Example usage:
ti = 13
tf = 23
y_min = 0.405
y_max = 0.45
x_min = 10.0
x_max = 22.5
Nsource = 0
Nsink = 0
plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(4))
plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(0.01))
a, b = calculate_mass_and_plot(ti, tf, x_min, x_max, y_min, y_max, Nsource, Nsink)

plot_line(a, b, 20, 22.0, plt.gca().lines[-1].get_color())

#plt.show()
plt.savefig(f'../../../plots/N{Nsource}_N{Nsink}_2.pdf', dpi=300, bbox_inches='tight')



plt.figure()
plt.style.use("paperdraft.mplstyle")
plt.figure(figsize=(4.5, 4))

# Example usage:
ti = 13
tf = 23
y_min = 0.395
y_max = 0.48
x_min = 5.5
x_max = 22.5
Nsource = 80
Nsink = 40
noAPE = True
plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(4))
plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(0.02))
a, b = calculate_mass_and_plot(ti, tf, x_min, x_max, y_min, y_max, Nsource, Nsink, noAPE)

plot_line(a, b, 18, 22.0, plt.gca().lines[-1].get_color())

#plt.show()
plt.savefig(f'../../../plots/N{Nsource}_N{Nsink}_noAPE2.pdf', dpi=300, bbox_inches='tight')
