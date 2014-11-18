from mpld3 import plugins
import os.path

class ClickPlugin(plugins.PluginBase):
    """Plugin for getting info on click"""

    # os.getcwd()
    f=open(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/javascripts/clickPlugin.js",'r');
    JAVASCRIPT = f.read()
    f.close()

    def __init__(self,eventHandlerURL,radioButtonID):
        self.dict_ = {"type": "clickPlugin","eventHandlerURL": eventHandlerURL,"radioButtonID":radioButtonID}

