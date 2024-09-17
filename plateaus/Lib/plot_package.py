import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

from matplotlib import rc, rcParams, cm, colors, colormaps
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm, colors
import numpy as np
import bootstrap

import os

os.environ["PATH"] = (
    "/Library/Frameworks/Python.framework/Versions/3.11/bin:/usr/local/bin:/System/Cryptexes/App/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/local/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/appleinternal/bin"
)


scale = 1
W = 8 * scale
r = 0.8

rcParams.update(
    {
        "figure.figsize": (W, W * r),  # 4:3 aspect ratio
        "font.size": 16 * scale,  # Set font size to 11pt
        "axes.labelsize": 18 * scale,  # -> axis labels
        "legend.fontsize": 16 * scale,  # -> legends
        "lines.markersize": 5 * scale,
        "font.family": "lmodern",
        "text.usetex": True,
        "text.latex.preamble": (  # LaTeX preamble
            r"\usepackage{lmodern}"
            # ... more packages if needed
        ),
    }
)


def plot_line(v, e, ti, tf, color):
    plt.gca().add_patch(
        plt.Rectangle(
            [ti - 0.2, v - e], tf - ti + 0.4, 2 * e, facecolor=color, alpha=0.4
        )
    )
