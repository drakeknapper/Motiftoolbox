mpld3.register_plugin("clickinfo", ClickInfo);
ClickInfo.prototype = Object.create(mpld3.Plugin.prototype);
ClickInfo.prototype.constructor = ClickInfo;

function ClickInfo(fig, props){
    mpld3.Plugin.call(this, fig, props);
};


ClickInfo.prototype.draw = function(){
    var obj = mpld3.get_element(this.props.id);

    var fig = this.fig;
    var coords = fig.canvas.append("text").attr("class", "mpld3-coordinates").style("text-anchor", "end").style("font-size", 12).attr("x", this.fig.width - 5).attr("y", this.fig.height - 5);

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
                                                fmt = d3.format(".3g");
                                                var pos = d3.mouse(this), x = ax.x.invert(pos[0]), y = ax.y.invert(pos[1]);
                                                coords.text("value of i"+i)
                                                var client = new XMLHttpRequest();
                                                client.open("GET", "/updatetorus?type=trace&xval="+fmt(x)+"&yval="+fmt(y));
                                                //if(i==1)client.open("GET", "/updatetorus?type=sweep");
                                                client.send();
                                                window.location.reload(true);
                                                parent.location.reload(true);
                                            }
                              )
            }else if(i==1){
                ax.baseaxes.on("mousedown", function()
                                            {
                                                fmt = d3.format(".3g");
                                                var pos = d3.mouse(this), x = ax.x.invert(pos[0]), y = ax.y.invert(pos[1]);
                                                coords.text("value of i"+i)
                                                var client = new XMLHttpRequest();
                                                //client.open("GET", "/updatetorus?type=trace&xval="+fmt(x)+"&yval="+fmt(y));
                                                client.open("GET", "/updatetorus?type=sweep");
                                                client.send();
                                                window.location.reload(true);
                                                parent.location.reload(true);
                                            }
                              )
                            }
            }


            ;
};



