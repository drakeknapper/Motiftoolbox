    mpld3.register_plugin("drag", DragPlugin);
    DragPlugin.prototype = Object.create(mpld3.Plugin.prototype);
    DragPlugin.prototype.constructor = DragPlugin;
    DragPlugin.prototype.requiredProps = ["eventHandlerURL","radioButtonID"];
    DragPlugin.prototype.defaultProps = {}
    function DragPlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
        mpld3.insert_css("#" + fig.figid + " path.dragging",
                         {"fill-opacity": "1.0 !important",
                          "stroke-opacity": "1.0 !important"});
    };

    DragPlugin.prototype.draw = function(){
        var fig = this.fig;
        var ax = fig.axes[0];
        fmt = d3.format(".3g");
        var coords = fig.canvas.append("text").attr("class", "mpld3-coordinates").style("text-anchor", "end").style("font-size", 12).attr("x", this.fig.width - 5).attr("y", this.fig.height - 5);
        var eventHandlerURL = this.props.eventHandlerURL;
        var radioButtonID = this.props.radioButtonID;

        var drag = d3.behavior.drag()
            .on("dragstart", dragstarted)
            .on("drag", dragged)
            .on("dragend", dragended);

        ax.baseaxes.call(drag);

        var startX,startY,endX,endY;
        function dragstarted(d) {
          var pos = d3.mouse(this), startX = ax.x.invert(pos[0]), startY = ax.y.invert(pos[1]);
          coords.text(fmt(startX)+","+fmt(startY));
        }

        function dragged(d, i) {
        }

        function dragended(d) {
          var pos = d3.mouse(this), endX = ax.x.invert(pos[0]), endY = ax.y.invert(pos[1]);
          coords.text(coords.text()+","+fmt(endX)+","+fmt(endY));


            var ret = returnDragType();
                if(ret=="0") return;

            var v = coords.text().split(",");
            var client = new XMLHttpRequest();
            client.open("GET", "/"+eventHandlerURL+"?type="+ret+"&startX="+v[0]+"&startY="+v[1]+"&endX="+v[2]+"&endY="+v[3]);
            client.send();
            coords.text("");

            window.location.reload(true);
            parent.frames[1].location.reload(true);	
        }
        function returnDragType(){
            var radios = parent.document.getElementsByName(radioButtonID);
            for (var i = 0, length = radios.length; i < length; i++) {
                if (radios[i].checked) {
                    return radios[i].value;

                }
            }
        }
    }
