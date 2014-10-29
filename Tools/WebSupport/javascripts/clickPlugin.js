mpld3.register_plugin("clickPlugin", ClickPlugin);
ClickPlugin.prototype = Object.create(mpld3.Plugin.prototype);
ClickPlugin.prototype.constructor = ClickPlugin;
ClickPlugin.prototype.requiredProps = ["eventHandlerURL","radioButtonID"];

function ClickPlugin(fig, props){
    mpld3.Plugin.call(this, fig, props);
};



ClickPlugin.prototype.draw = function(){
    var obj = mpld3.get_element(this.props.id);

    var fig = this.fig;
    var coords = fig.canvas.append("text").attr("class", "mpld3-coordinates").style("text-anchor", "end").style("font-size", 12).attr("x", this.fig.width - 5).attr("y", this.fig.height - 5);
    var eventHandlerURL = this.props.eventHandlerURL;
    var radioButtonID = this.props.radioButtonID;

    for (var i = 0; i < fig.axes.length; i++) {
            var ax = fig.axes[i];
            var update_coords = function() {
              fmt = d3.format(".3g");
              var pos = d3.mouse(this), x = ax.x.invert(pos[0]), y = ax.y.invert(pos[1]);
              coords.text("mouse is at : (" + fmt(x) + ", " + fmt(y) + ")");
            };
            ax.baseaxes.on("mousemove",update_coords)
                       .on("mouseout", function()
                                {coords.text("");}
                            )
            if(i==0){
                ax.baseaxes.on("mousedown", function()
                                            {
                                                var ret = returnClickType();
                                                if(ret=="0") return;

                                                fmt = d3.format(".3g");
                                                var pos = d3.mouse(this), x = ax.x.invert(pos[0]), y = ax.y.invert(pos[1]);
                                                coords.text("value of i"+i)
                                                var client = new XMLHttpRequest();
                                                client.open("GET", "/"+eventHandlerURL+"?type=trace&xval="+fmt(x)+"&yval="+fmt(y));
                                                client.send();
                                                window.location.reload(true);
						parent.frames[1].location.reload(true);
                                            }
                              )
            }else if(i==1){
                ax.baseaxes.on("mousedown", function()
                                            {
                                                var ret = returnClickType();
                                                if(ret=="0") return;

                                                fmt = d3.format(".3g");
                                                var pos = d3.mouse(this), x = ax.x.invert(pos[0]), y = ax.y.invert(pos[1]);
                                                coords.text("value of i"+i)
                                                var client = new XMLHttpRequest();
                                                client.open("GET", "/"+eventHandlerURL+"?type=sweep");
                                                client.send();
                                                window.location.reload(true);	// only needs to update the phase torus
                                            }
                              )
                            }
            }


            ;
    function returnClickType(){
            var radios = parent.document.getElementsByName(radioButtonID);
            for (var i = 0, length = radios.length; i < length; i++) {
                if (radios[i].checked) {
                    return radios[i].value;
                }
            }
    }
};



