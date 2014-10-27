import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import mpld3
from mpld3 import plugins, utils


class DragPlugin(plugins.PluginBase):
    f=open(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/javascripts/dragPlugin.js",'r');
    JAVASCRIPT = f.read()
    f.close()


    def __init__(self, points):
        if isinstance(points, mpl.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None

        self.dict_ = {"type": "drag",
                     }


if __name__ == "__main__" :

    fig, ax = plt.subplots()
    np.random.seed(0)
    points = ax.plot(np.random.normal(size=20),
                     np.random.normal(size=20), 'or', alpha=0.5,
                     markersize=50, markeredgewidth=1)
    ax.set_title("Click and Drag", fontsize=18)

    plugins.connect(fig, DragPlugin(points[0]))

    mpld3.show()
