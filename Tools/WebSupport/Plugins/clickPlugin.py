from libxml2mod import name
import numpy as np
import matplotlib.pyplot as plt
import mpld3
from mpld3 import utils
from mpld3 import plugins
import os.path

class ClickPlugin(plugins.PluginBase):
    """Plugin for getting info on click"""

    # os.getcwd()
    f=open(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/javascripts/clickPlugin.js",'r');
    JAVASCRIPT = f.read()
    f.close()

    def __init__(self,fig):
        self.dict_ = {"type": "clickPlugin"}



if __name__ == "__main__" :
    fig, ax = plt.subplots()
    points = ax.scatter(np.random.rand(50), np.random.rand(50),
                        s=500, alpha=0.3)
    #plugins.connect(fig, ClickInfo(points), plugins.MousePosition())
    plugins.connect(fig, ClickPlugin(points))
    mpld3.show()
