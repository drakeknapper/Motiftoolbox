from libxml2mod import name
import numpy as np
import matplotlib.pyplot as plt
import mpld3
from mpld3 import utils
from mpld3 import plugins


class ClickInfo(plugins.PluginBase):
    """Plugin for getting info on click"""

    f=open("customplugin.js",'r');
    JAVASCRIPT = f.read()
    f.close()

    #var div = d3.select("#" + this.fig.figid);
    def __init__(self, points):
        self.dict_ = {"type": "clickinfo"}



if __name__ == "__main__" :
    fig, ax = plt.subplots()
    points = ax.scatter(np.random.rand(50), np.random.rand(50),
                        s=500, alpha=0.3)
    #plugins.connect(fig, ClickInfo(points), plugins.MousePosition())
    plugins.connect(fig, ClickInfo(points))
    mpld3.show()
