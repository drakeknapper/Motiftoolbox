import os.path
from mpld3 import plugins, utils


class DragPlugin(plugins.PluginBase):
    f=open(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/javascripts/dragPlugin.js",'r');
    JAVASCRIPT = f.read()
    f.close()
    def __init__(self, eventHandlerURL, radioButtonID):
        self.dict_ = {"type": "drag",
                      "eventHandlerURL": eventHandlerURL,
                      "radioButtonID": radioButtonID
                     }
