require.config({
              shim: {
                   'jquery.flot': {
                       deps: ['jquery'],
                       exports: '$.plot'
                   },
                   'jquery.flot.symbol': {
                       deps: ['jquery.flot']
                   },
                   'jquery.flot.dashes': {
                       deps: ['jquery.flot']
                   },
                   'jquery.flot.stack': {
                       deps: ['jquery.flot']
                   },
                   'jquery.flot.orderBars': {
                       deps: ['jquery.flot']
                   },
                   'jspimple': {
                       deps: ['jquery.flot.orderBars','jquery.flot.stack','jquery.flot.dashes','jquery.flot.symbol'],
                   }
               },
               paths: {
               // the left side is the module ID,
               // the right side is the path to
               // the jQuery file, relative to baseUrl.
               // Also, the path should NOT include
               // the '.js' file extension. This example
               // is using jQuery 1.9.0 located at
               // js/lib/jquery-1.9.0.js, relative to
               // the HTML page.
               'jquery': 'jquery-1.11.0.min',
               'jquery.flot': 'flot/jquery.flot',
               'jquery.flot.symbol': 'flot/jquery.flot.symbol',
               'jquery.flot.dashes': 'flot/jquery.flot.dashes',
               'jquery.flot.stack': 'flot/jquery.flot.stack',
               'jquery.flot.orderBars': 'flot/jquery.flot.orderBars',
               'jspimple': 'jspimple'
               }
               });

define(['jquery','jquery.flot','jquery.flot.symbol','jquery.flot.dashes','jquery.flot.stack','jquery.flot.orderBars'],function(apple,banana,pineapple,orange,plum,pear,strawberry,raspberry,lemon){

function reescape(s){
	var h= ['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'];
	s= s.replace(/&nbsp;/g,' ').split('');
	var k, c, m,n,o,p;
	for(var i= 0; i<s.length; i++){
		k= s[i];
		c= k.charCodeAt(0);
		s[i]= (function(){
			switch(c){
				case 60: return "&lt;";
				case 62: return "&gt;";
				case 34: return "&quot;";
				case 38: return "&amp;";
				default: 
				if(c> 127){
					m= c%16;
					c= Math.floor(c/16);
					n= c%16;
					c= Math.floor(c/16);
					o= c%16;
					c= Math.floor(c/16);
					p= c%16;
					return "&#x"+h[p]+h[o]+h[n]+h[m]+";";
				}
				else{
					return k;
				}
			}
		})();
	}
	return s.join('');
}

function arrayBufferToBinaryString(buffer) {
  var binary = "";
  var bytes = new Uint8Array(buffer);
  var length = bytes.byteLength;
  for (var i = 0; i < length; i++) {
    binary += String.fromCharCode(bytes[i]);
  }
  return binary;
}

/*
if (window.File && window.FileReader && window.FileList && window.Blob) {
  // Great success! All the File APIs are supported.
} else {
  alert('The File APIs are not fully supported in this browser.');
}
*/

function onMenuSelect(value){
  console.log("onMenuSelect "+value);
}

function CCP4Graph(graphName,d,xbreaks,ybreaks,options,shapes,titles,background,customXTicks){
  this.graphID = graphName;
  this.graphData = d;
  this.xbreaks = xbreaks;
  this.ybreaks = ybreaks;
  this.options = options;
  this.shapes = shapes;
  this.titles = titles;
  this.background = background;
  this.customXTicks = customXTicks;
}

function CCP4GraphPlot(graphDivName_in,showInput,allowResizeEvent){

  // FIXME - This could be in the XML.
  var useLocalStorage = true;

  this.currentStoredGraph = -1;

  if(useLocalStorage){
    if(typeof(Storage) !== "undefined") {
      if (sessionStorage.getItem(graphDivName_in)) {
          this.currentStoredGraph = parseInt(sessionStorage.getItem(graphDivName_in));
      } else {
      }
    } else {
      console.log("No local storage");
    }
  }


  this.currentGraph = null;
  this.allowResizeEvent = allowResizeEvent;
  //console.log("resize event allowed?")
  //console.log(allowResizeEvent)
  var pimpMain = document.createElement("div");
  var graphDivOuter = document.getElementById(graphDivName_in)
  if(!graphDivOuter){
    console.log("Element "+graphDivName_in+" does not exist");
    this.graphDivName = null;
    return;
  }
  this.graphDivName_in = graphDivName_in;
  var N = 20;
  var graphDivName = new Array(N+1).join((Math.random().toString(36)+'00000000000000000').slice(2, 18)).slice(0, N);
  var menuSelectName = new Array(N+1).join((Math.random().toString(36)+'00000000000000000').slice(2, 18)).slice(0, N);
  pimpMain.setAttribute("id",graphDivName);
  var dummy = document.createElement("div");
  dummy.setAttribute("id","dummy")
  var parentWidth = parseInt(graphDivOuter.style.width);
  var parentHeight = parseInt(graphDivOuter.style.height);
  var w = parentWidth+"px";
  var h = (parentHeight-30)+"px"; //FIXME Hardwire should be (approx) 2 X max(menu font,input font) + 10.
  dummy.setAttribute("style","width:"+w+";"+"height:"+h+";")
  pimpMain.setAttribute("style","width:"+w+";"+"height:"+h+";")
  pimpMain.appendChild(dummy)

  graphDivOuter.appendChild(pimpMain)
  if(showInput && window.File && window.FileReader && window.FileList && window.Blob) {
    var input = document.createElement("input");
    input.setAttribute("type","file")
    var fileIDName = new Array(N+1).join((Math.random().toString(36)+'00000000000000000').slice(2, 18)).slice(0, N);
    input.setAttribute("id",fileIDName)
    input.setAttribute("name","files=[]")
    input.setAttribute("multiple","")
    input.setAttribute("id",fileIDName)
    input.setAttribute("style","float:left;");
    input.setAttribute("class","graphInput")
  }

  var menuSelect = document.createElement("select");
  menuSelect.setAttribute("class","graphSelect");
  menuSelect.setAttribute("id",menuSelectName);
  menuSelect.setAttribute("name",menuSelectName);
  

  //menuSelect.setAttribute("onchange","onMenuSelect(value);");
  var selectItem = document.createElement("option");
  selectItem.setAttribute("value","Graphs");
  selectItem.innerHTML = "Graphs";
  menuSelect.appendChild(selectItem)
  graphDivOuter.appendChild(menuSelect)

  if (showInput && window.File && window.FileReader && window.FileList && window.Blob) {
    graphDivOuter.appendChild(input)
  }

  this.graphs = [];
  this.fonts = null;
   
  this.graphDivName = graphDivName;
  this.menuSelectName = menuSelectName;

  var self = this;

  if (showInput && window.File && window.FileReader && window.FileList && window.Blob) {
    input.addEventListener('change', this, false);
  }
  window.addEventListener("resize",this,false);

  $( "#"+menuSelectName ).change(function() {
   $( "#"+menuSelectName+" option:selected" ).each(function() {
     var graphName = $(this)[0].value.replace("_selectitem","");
     var thisPlotDiv = document.getElementById(graphName);
     $('#'+self.graphDivName).children().hide();
     thisPlotDiv.setAttribute("style","width:100%;height:"+thisPlotDiv.style.height+";display:block;");
     $("#"+graphName).empty();
     self.currentGraph = self.getDataWithName(graphName);
     self.plot();
     var selectedIndex = $("#"+menuSelectName).prop('selectedIndex');
     if(useLocalStorage){
       if(typeof(Storage) !== "undefined") {
         sessionStorage.setItem(graphDivName_in,""+selectedIndex);
       }
     }
   });
  });

  w = parentWidth+"px";
  h = (parentHeight-menuSelect.offsetHeight-2)+"px";
  dummy.setAttribute("style","width:"+w+";"+"height:"+h+";")
  pimpMain.setAttribute("style","width:"+w+";"+"height:"+h+";")

}

CCP4GraphPlot.prototype.handleAbort = function (e) {
  //console.log("aborted");
}
CCP4GraphPlot.prototype.handleError = function (e) {
  //console.log("error");
}
CCP4GraphPlot.prototype.handleLoadStart = function (e) {
  //console.log("loadstart");
}
CCP4GraphPlot.prototype.handleLoadEnd = function (e) {
  //console.log("loadend");
}
CCP4GraphPlot.prototype.handleLoadProgress = function (e) {
  //console.log("progress");
}

CCP4GraphPlot.prototype.parseColourToRGBA = function(col,alpha) {
  if(col=="r"||col=="red"){
    col = "#ff0000"
  }
  if(col=="g"||col=="green"){
    col = "#00ff00"
  }
  if(col=="b"||col=="blue"){
    col = "#0000ff"
  }
  if(col=="c"||col=="cyan"){
    col = "#00ffff"
  }
  if(col=="m"||col=="magenta"){
    col = "#ff00ff"
  }
  if(col=="y"||col=="yellow"){
    col = "#ffff00"
  }
  if(col=="w"||col=="white"){
    col = "#ffffff"
  }
  if(col=="k"||col=="black"){
    col = "#000000"
  }
  var r = parseInt('0x'+col.substring(1,3))
  var g = parseInt('0x'+col.substring(3,5))
  var b = parseInt('0x'+col.substring(5,7))
  return "rgba("+r+","+g+","+b+","+alpha+")"
}

CCP4GraphPlot.prototype.handleLoad = function(e) {
  var txt = arrayBufferToBinaryString(e.target.result);
  var xmlDoc;
  if (window.DOMParser) {
    var parser=new DOMParser();
    xmlDoc=parser.parseFromString(txt,"text/xml");
  } else {
    xmlDoc=new ActiveXObject("Microsoft.XMLDOM");
    xmlDoc.async=false;
    xmlDoc.loadXML(txt);
  } 
  //console.log(xmlDoc);
  //console.log("read XML");
  var root = xmlDoc.getElementsByTagName("CCP4ApplicationOutput")[0];
  var tables = root.getElementsByTagName("CCP4Table")
  var fonts = root.getElementsByTagName("Fonts")

  var graphObject = this.graphObject;
  graphObject.loadXML(graphObject,root,tables,fonts);
}

function getNodeText(node) {
    var r = "";
    for (var x = 0;x < node.childNodes.length; x++) {
        r = r + node.childNodes[x].nodeValue;
    }
    return r;
}

CCP4GraphPlot.prototype.loadXML = function(graphObject,root,tables,fonts) {

  if(!graphObject.graphDivName){
    return;
  }

  graphObject.fonts = fonts;

  pimpleplot = document.getElementById(graphObject.graphDivName);

  var igraph = 0;
  while (pimpleplot.firstChild) {
     pimpleplot.removeChild(pimpleplot.firstChild);
  }
  var thisPlotDiv;
  graphObject.graphs = [];

  var myNodeSelect = document.getElementById(graphObject.menuSelectName);
  while (myNodeSelect.firstChild) {
     myNodeSelect.removeChild(myNodeSelect.firstChild);
  }

  var N = 20;

  var graph;
  for(var itab=0;itab<tables.length;itab++){
    //console.log(tables[itab])
    var title = tables[itab].getAttribute("title");
    var list = document.createElement('ul');
    var data = $.trim(getNodeText(tables[itab].getElementsByTagName("data")[0])).split("\n");
    var headers_sep = tables[itab].getElementsByTagName("headers")[0].getAttribute("separator");
    var headers = $.trim(tables[itab].getElementsByTagName("headers")[0].childNodes[0].nodeValue).replace(/ +/g,' ').split(headers_sep);
    var dataColRow = [];
    for(var idatacol=0;idatacol<headers.length;idatacol++){
      dataColRow.push([]);
    }
    for(var idatarow=0;idatarow<data.length;idatarow++){
      var thisDataRow = $.trim(data[idatarow]).replace(/ +/g,' ').split(" ");
      for(var idatacol=0;idatacol<headers.length;idatacol++){
        if(thisDataRow[idatacol]=="-"){
          //console.log("A null!")
          dataColRow[idatacol].push(null);
        }else{
          dataColRow[idatacol].push(parseFloat(thisDataRow[idatacol]));
        }
      }
    }
    //console.log(title);
    //console.log(headers);
    //console.log(headers.length);
    //console.log(data);
    //console.log(dataColRow);

    var plots = tables[itab].getElementsByTagName("plot")

    var selectGroup = document.createElement("optgroup");
    selectGroup.setAttribute("label",title.replace("<","&lt;").replace(">","&gt;"));
    myNodeSelect.appendChild(selectGroup)

    for(var iplot=0;iplot<plots.length;iplot++){
      var plot = plots[iplot];
      var plot_title = plot.getElementsByTagName("title")[0].childNodes[0].nodeValue.replace("<","&lt;").replace(">","&gt;");;

      //console.log(plot_title);
      var plotlines = plot.getElementsByTagName("plotline");
      //console.log(plotlines.length);
      var d = [];
      for(var iplotline=0;iplotline<plotlines.length;iplotline++){
        var plotline = plotlines[iplotline];
        var xcol = parseInt(plotline.getAttribute("xcol"))-1;
        var ycol = parseInt(plotline.getAttribute("ycol"))-1;
        var label;
        try {
          label = plotline.getElementsByTagName("label")[0].childNodes[0].nodeValue;
        } catch(err) {
          label = headers[ycol];
        }
        
        var averagebatches = false;
        var averagebatchcol = null;
        try {
          averagebatches = plotline.getElementsByTagName("averagebatches")[0].childNodes[0].nodeValue;
          if(averagebatches=="true"||averagebatches=="1"){
            averagebatches = true;
            averagebatchcol = parseInt(plotline.getElementsByTagName("averagebatchcol")[0].childNodes[0].nodeValue)-1;
          } else {
           averagebatches = false;
           averagebatchcol = null;
          }
        } catch(err) {
          averagebatches = null;
        }

        var d1 = [];

        if(!averagebatches){
          for(var idatarow=0;idatarow<data.length;idatarow++){
            d1.push([dataColRow[xcol][idatarow],dataColRow[ycol][idatarow]]);
          }
        } else {
          var oldbatchnum = null;
          var nbatch = 0;
          var thisBatchVal = 0;
          var theXCol;
          for(var idatarow=0;idatarow<data.length;idatarow++){
            var batchnum = dataColRow[averagebatchcol][idatarow];
            if(batchnum==oldbatchnum||oldbatchnum===null){
               thisBatchVal += dataColRow[ycol][idatarow];
               theXCol = dataColRow[xcol][idatarow];
               nbatch++;
            } else {
               var avg = thisBatchVal/nbatch;
               d1.push([theXCol,avg]);
               nbatch = 0;
               thisBatchVal = 0;
            }
            oldbatchnum = batchnum;
          }
          var avg = thisBatchVal/nbatch;
          d1.push([theXCol,avg]);
        }

        //console.log("Number of points: "+d1.length)

        var theData;

        try {
          var visible = plotline.getElementsByTagName("visible")[0].childNodes[0].nodeValue;
          if(visible&&(visible=="false"||visible=="0")){
            continue;
          }
        } catch(err) {
        }

        try {
          var showinlegend = plotline.getElementsByTagName("showinlegend")[0].childNodes[0].nodeValue;
          if(showinlegend&&(showinlegend=="false"||showinlegend=="0")){
            theData = {data:d1}
          } else {
            theData = {label:label,data:d1}
          }
        } catch(err) {
          theData = {label:label,data:d1}
        }

        try {
          var colour = plotline.getElementsByTagName("colour")[0].childNodes[0].nodeValue;
          if(colour=="r") colour = "red";
          if(colour=="g") colour = "green";
          if(colour=="b") colour = "blue";
          if(colour=="y") colour = "yellow";
          if(colour=="m") colour = "magenta";
          if(colour=="c") colour = "cyan";
          if(colour=="w") colour = "white";
          if(colour=="b") colour = "black";
          theData["color"] = colour;
        } catch(err) {
        }

        try {
          var rightaxis = plotline.getAttribute("rightaxis");
          if(!rightaxis||rightaxis=="false"||rightaxis==0){
          } else {
            theData["yaxis"] = 2
          }
        } catch(err) {
        }

        var haveLine = true;
        var theData2;
        var haveDash = false;
        try {
          var linestyle = plotline.getElementsByTagName("linestyle")[0].childNodes[0].nodeValue;
          if(linestyle=="--"){ //dash
            theData["lines"] = {show:false}
            haveLine = false;
            haveDash = true;
            theData["dashes"] = {show:true,dashLength:[8,8]}
          }else if(linestyle==":"){ //dot
            theData["lines"] = {show:false}
            haveLine = false;
            haveDash = true;
            theData["dashes"] = {show:true,dashLength:[4,4]}
          }else if(linestyle=="-."){ //dashdot
            theData["lines"] = {show:false}
            theData["shadowSize"] = 0
            haveLine = false;
            haveDash = true;
            theData2 = jQuery.extend(true, {}, theData);
            theData["dashes"] = {show:true,dashLength:[12,12]}
            theData2["dashes"] = {show:true,dashLength:[3,3]}
            if("label" in theData2){
              delete theData2["label"];
            }
          }else if(linestyle=="."||linestyle=="_"){ //blank
            theData["lines"] = {show:false}
            haveLine = false;
          }
        } catch(err) {
        }

        try {
          var symbol = plotline.getElementsByTagName("symbol")[0].childNodes[0].nodeValue;
          if(symbol=="o"){
            symbol = "circle"
          }else if(symbol=="d"){
            symbol = "diamond"
          }else if(symbol=="s"){
            symbol = "square"
          }else if(symbol=="^"){
            symbol = "triangle"
          }else if(symbol=="x"){
            symbol = "cross"
          }else{
            symbol = "circle"
          }

          if (plotline.getElementsByTagName("fillcolour").length > 0){
            fillcolour=plotline.getElementsByTagName("fillcolour")[0].childNodes[0].nodeValue;
            if(fillcolour == "false"){
              theData["points"] = {symbol:symbol, fill:1.0, fillColor:false, show:true}
            }else{
              theData["points"] = {symbol:symbol, fillColor:fillcolour, show:true}
            }
          } else {
            theData["points"] = {symbol:symbol, show:true}
          }
 
          if(haveLine){
            theData["lines"] = {show:true}
          }
        } catch(err) {
        }

        try {
          var symbolsize = parseInt(plotline.getElementsByTagName("symbolsize")[0].childNodes[0].nodeValue)/2+1;
          if("points" in theData){
            theData["points"]["radius"] = symbolsize
          } else {
            theData["points"] = {show:true, radius:symbolsize}
          }
          if(haveLine){
            theData["lines"] = {show:true}
          }
        } catch(err) {
        }

        try {
          var linesize = parseInt(plotline.getElementsByTagName("linesize")[0].childNodes[0].nodeValue);
          if(haveLine && "lines" in theData){
            theData["lines"]["lineWidth"] = linesize;
          }
          if(haveDash && "dashes" in theData){
            theData["dashes"]["lineWidth"] = linesize;
            if(typeof theData2!=="undefined" && "dashes" in theData){
              theData2["dashes"]["lineWidth"] = linesize;
            }
          }
        } catch(err) {
        }

        if(averagebatches){
          try {
            // This is arguably naughty. I turn off shadows to boost performance *assuming* that is why
            // downsampling was requested.
            theData["shadowSize"] = 0
          } catch(err) {
          }
        }

        //console.log(theData)
        d.push(theData);
        if(typeof theData2!=="undefined"){
          d.push(theData2);
        }
      }

      var xrange = {};
      var yrange = {};
      var ryrange = {};
      var xranges = plot.getElementsByTagName("xrange");
      if(xranges.length>0){
        var rangeEl = xranges[0];
        xrange = {}
        try {
          var minx = parseFloat(rangeEl.getAttribute("min"))
          if(!isNaN(minx)){
            xrange["min"] = minx
          }
        } catch(err){
        }
        try {
          var maxx = parseFloat(rangeEl.getAttribute("max"))
          if(!isNaN(maxx)){
            xrange["max"] = maxx
          }
        } catch(err){
        }
      }
      var yranges = plot.getElementsByTagName("yrange");
      for(var iyr=0;iyr<yranges.length;iyr++){
        if(!yranges[iyr].getAttribute("rightaxis")|| yranges[iyr].getAttribute("rightaxis")=="false" || yranges[iyr].getAttribute("rightaxis")=="0"){
          yrange = {}
          try {
            var miny = parseFloat(yranges[iyr].getAttribute("min"))
            if(!isNaN(miny)){
              yrange["min"] = miny
            }
          } catch(err){
          }
          try {
            var maxy = parseFloat(yranges[iyr].getAttribute("max"))
            if(!isNaN(maxy)){
              yrange["max"] = maxy
            }
          } catch(err){
          }
        } else {
          ryrange = {}
          try {
            var miny = parseFloat(yranges[iyr].getAttribute("min"))
            if(!isNaN(miny)){
              ryrange["min"] = miny
            }
          } catch(err){
          }
          try {
            var maxy = parseFloat(yranges[iyr].getAttribute("max"))
            if(!isNaN(maxy)){
              ryrange["max"] = maxy
            }
          } catch(err){
          }
        }
      }

      // FIXME - WE ONLY CONSIDER XBREAKS SO FAR.
      var xbreaks = plot.getElementsByTagName("xbreaks");
      var xbreak = [];
      if(xbreaks.length>0){
        var xbreakEl = xbreaks[0];
        var breaks = xbreakEl.getElementsByTagName("break");
        if(breaks.length<500){
          //console.log("Doing breaks")
          for(var ibr=0;ibr<breaks.length;ibr++){
            var breakEl = breaks[ibr];
            xbreak.push(parseFloat(breakEl.getAttribute("min")),parseFloat(breakEl.getAttribute("max")))
          }
        } else {
          //console.log("Not doing silly breaks")
        }
      }
      var ybreak = [];
      var options = {xrange:xrange,yrange:yrange,ryrange:ryrange};
      var shapes = [];
      /*
      var circle = {circle:{x:3,y:3,radius:.75,linesize:3,linecolour:"rgba(0, 255, 0, 0.5)",fillcolour:"rgba(255, 255, 0, 0.5)",linestyle:[3,2,8,2]}};
      shapes.push(circle);
      var circle = {circle:{x:6,y:3,radius:.5,linesize:10,linecolour:"rgba(0, 255, 0, 0.5)",fillcolour:"rgba(255, 255, 0, 0.5)",linestyle:[4,4]}};
      shapes.push(circle);
      */

      try {
          var linestrips = plot.getElementsByTagName("line")
          //console.log("Have some lines");
          for(var ilstrip=0;ilstrip<linestrips.length;ilstrip++){

            var dataPoints = [];
            if(!linestrips[ilstrip].getAttribute("x1")){
              continue;
            } else {
              dataPoints.push(linestrips[ilstrip].getAttribute("x1"));
            }
            if(!linestrips[ilstrip].getAttribute("y1")){
              continue;
            } else {
              dataPoints.push(linestrips[ilstrip].getAttribute("y1"));
            }
            if(!linestrips[ilstrip].getAttribute("x2")){
              continue;
            } else {
              dataPoints.push(linestrips[ilstrip].getAttribute("x2"));
            }
            if(!linestrips[ilstrip].getAttribute("y2")){
              continue;
            } else {
              dataPoints.push(linestrips[ilstrip].getAttribute("y2"));
            }

            var p = [];
            for(var ip=0;ip<dataPoints.length/2;ip++){
              p.push([dataPoints[2*ip],dataPoints[2*ip+1]]);
            }

            var linecolour;
            if(!linestrips[ilstrip].getAttribute("linecolour")){
              linecolour = "#000000";
            } else {
              linecolour = linestrips[ilstrip].getAttribute("linecolour")
            }
            var linestyle;
            if(!linestrips[ilstrip].getAttribute("linestyle")){
              linestyle = null;
            } else {
              linestyle = linestrips[ilstrip].getAttribute("linestyle")
              if(linestyle=="dashdot"||linestyle=="-."){
                linestyle = [3,2,8,2];
              } else if(linestyle=="dash"||linestyle=="--"){
                linestyle = [8,8];
              } else if(linestyle=="dot"||linestyle==":"){
                linestyle = [4,4];
              } else {
                linestyle = null;
              }
            }
            var linesize;
            if(!linestrips[ilstrip].getAttribute("linesize")){
              linesize = 1;
            } else {
              linesize = parseInt(linestrips[ilstrip].getAttribute("linesize"))
            }
            linecolour = graphObject.parseColourToRGBA(linecolour,1.0)
            //console.log(p+" "+linesize+" "+linestyle+" "+linecolour+" "+fillcolour+" "+alpha)
            var linestrip = {linestrip:{data:p,linesize:linesize,linecolour:linecolour,linestyle:linestyle}};
            shapes.push(linestrip);
          }
      } catch(err) {
      }

      try {
          var linestrips = plot.getElementsByTagName("linestrip")
          //console.log("Have some linestrips");
          for(var ilstrip=0;ilstrip<linestrips.length;ilstrip++){
            var dataPoints = (linestrips[ilstrip].childNodes[0].nodeValue).trim().split(/\s+/);
            var p = [];
            for(var ip=0;ip<dataPoints.length/2;ip++){
              p.push([dataPoints[2*ip],dataPoints[2*ip+1]]);
            }
            var linecolour;
            if(!linestrips[ilstrip].getAttribute("linecolour")){
              linecolour = "#000000";
            } else {
              linecolour = linestrips[ilstrip].getAttribute("linecolour")
            }
            var linestyle;
            if(!linestrips[ilstrip].getAttribute("linestyle")){
              linestyle = null;
            } else {
              linestyle = linestrips[ilstrip].getAttribute("linestyle")
              if(linestyle=="dashdot"||linestyle=="-."){
                linestyle = [3,2,8,2];
              } else if(linestyle=="dash"||linestyle=="--"){
                linestyle = [8,8];
              } else if(linestyle=="dot"||linestyle==":"){
                linestyle = [4,4];
              } else {
                linestlye = null;
              }
            }
            var linesize;
            if(!linestrips[ilstrip].getAttribute("linesize")){
              linesize = 1;
            } else {
              linesize = parseInt(linestrips[ilstrip].getAttribute("linesize"))
            }
            linecolour = graphObject.parseColourToRGBA(linecolour,1.0)
            //console.log(p+" "+linesize+" "+linestyle+" "+linecolour+" "+fillcolour+" "+alpha)
            var linestrip = {linestrip:{data:p,linesize:linesize,linecolour:linecolour,linestyle:linestyle}};
            shapes.push(linestrip);
          }
      } catch(err) {
          console.log("linestrip woes?");
      }

      try {
          var texts = plot.getElementsByTagName("text")
          for(var itext=0;itext<texts.length;itext++){
            var xpos = parseFloat(texts[itext].getAttribute("xpos"));
            var ypos = parseFloat(texts[itext].getAttribute("ypos"));
            var text = texts[itext].textContent;
            var rtext;
            var fillcolour;
            if(!texts[itext].getAttribute("colour")){
              fillcolour = "#000000";
            } else {
              fillcolour = texts[itext].getAttribute("colour")
            }
            if(texts[itext].getAttribute("font")){
              var font = texts[itext].getAttribute("font");
              rtext = {text:{x:xpos,y:ypos,text:text,font:font,fillcolour:fillcolour}};
            } else {
              rtext = {text:{x:xpos,y:ypos,text:text,font:"",fillcolour:fillcolour}};
            }
            shapes.push(rtext);
          }
      } catch(err) {
          console.log("text woes?");
      }

      try {
          var polygons = plot.getElementsByTagName("polygon")
          //console.log("Have some polygons");
          for(var ipoly=0;ipoly<polygons.length;ipoly++){
            var dataPoints = (polygons[ipoly].childNodes[0].nodeValue).trim().split(/\s+/);
            var p = [];
            for(var ip=0;ip<dataPoints.length/2;ip++){
              p.push([dataPoints[2*ip],dataPoints[2*ip+1]]);
            }
            var alpha;
            if(!polygons[ipoly].getAttribute("alpha")){
              alpha = 0.5;
            } else {
              alpha = parseFloat(polygons[ipoly].getAttribute("alpha"))
            }
            var linecolour;
            if(!polygons[ipoly].getAttribute("linecolour")){
              linecolour = "#000000";
            } else {
              linecolour = polygons[ipoly].getAttribute("linecolour")
            }
            var fillcolour;
            if(!polygons[ipoly].getAttribute("fillcolour")){
              fillcolour = "#000000";
            } else {
              fillcolour = polygons[ipoly].getAttribute("fillcolour")
            }
            var linestyle;
            if(!polygons[ipoly].getAttribute("linestyle")){
              linestyle = null;
            } else {
              linestyle = polygons[ipoly].getAttribute("linestyle")
              //console.log("in: "+linestyle);
              if(linestyle=="dashdot"||linestyle=="-."){
                linestyle = [3,2,8,2];
              } else if(linestyle=="dash"||linestyle=="--"){
                linestyle = [8,8];
              } else if(linestyle=="dot"||linestyle==":"){
                linestyle = [4,4];
              } else {
                linestlye = null;
              }
              //console.log("out: "+linestyle);
            }
            var linesize;
            if(!polygons[ipoly].getAttribute("linesize")){
              linesize = 1;
            } else {
              linesize = parseInt(polygons[ipoly].getAttribute("linesize"))
            }
            fillcolour = graphObject.parseColourToRGBA(fillcolour,alpha)
            linecolour = graphObject.parseColourToRGBA(linecolour,1.0)
            //console.log(p+" "+linesize+" "+linestyle+" "+linecolour+" "+fillcolour+" "+alpha)
            var polygon = {polygon:{data:p,linesize:linesize,linecolour:linecolour,fillcolour:fillcolour,linestyle:linestyle}};
            shapes.push(polygon);
          }
      } catch(err) {
          console.log("polygon woes?");
      }

      try {
          var circles = plot.getElementsByTagName("circle")
          for(var icirc=0;icirc<circles.length;icirc++){
            var alpha;
            if(!circles[icirc].getAttribute("alpha")){
              alpha = 0.5;
            } else {
              alpha = parseFloat(circles[icirc].getAttribute("alpha"))
            }
            var xpos = parseFloat(circles[icirc].getAttribute("xpos"));
            var ypos = parseFloat(circles[icirc].getAttribute("ypos"));
            var radius = parseFloat(circles[icirc].getAttribute("radius"));
            var linecolour;
            if(!circles[icirc].getAttribute("linecolour")){
              linecolour = "#000000";
            } else {
              linecolour = circles[icirc].getAttribute("linecolour")
            }
            var fillcolour;
            if(!circles[icirc].getAttribute("fillcolour")){
              fillcolour = "#000000";
            } else {
              fillcolour = circles[icirc].getAttribute("fillcolour")
            }
            var linestyle;
            if(!circles[icirc].getAttribute("linestyle")){
              linestyle = null;
            } else {
              linestyle = circles[icirc].getAttribute("linestyle")
              if(linestyle=="dashdot"||linestyle=="-."){
                linestyle = [3,2,8,2];
              } else if(linestyle=="dash"||linestyle=="--"){
                linestyle = [8,8];
              } else if(linestyle=="dot"||linestyle==":"){
                linestyle = [4,4];
              } else {
                linestlye = null;
              }
            }
            var linesize;
            if(!circles[icirc].getAttribute("linesize")){
              linesize = 1;
            } else {
              linesize = parseInt(circles[icirc].getAttribute("linesize"))
            }
            fillcolour = graphObject.parseColourToRGBA(fillcolour,alpha)
            linecolour = graphObject.parseColourToRGBA(linecolour,1.0)
            //console.log(xpos+" "+ypos+" "+radius+" "+linesize+" "+linestyle+" "+linecolour+" "+fillcolour+" "+alpha)
            var circle = {circle:{x:xpos,y:ypos,radius:radius,linesize:linesize,linecolour:linecolour,fillcolour:fillcolour,linestyle:linestyle}};
            shapes.push(circle);
          }
      } catch(err) {
      }

      var titles = {};

      var xaxislabel = null;
      try {
        xaxislabel = plot.getElementsByTagName("xlabel")[0].childNodes[0].nodeValue.replace("<","&lt;").replace(">","&gt;");;
      } catch(err) {
      }

      var showlegend = null;
      try {
        showlegend = plot.getElementsByTagName("showlegend")[0].childNodes[0].nodeValue;
        if(showlegend==="0"||showlegend==="false"){
          options["showlegend"]=false;
        } else {
          options["showlegend"]=true;
        }
      } catch(err) {
        options["showlegend"]=true;
      }

      var fixaspectratio = null;
      try {
        fixaspectratio = plot.getElementsByTagName("fixaspectratio")[0].childNodes[0].nodeValue;
        if(fixaspectratio==="1"||fixaspectratio==="true"){
          options["fixaspectratio"]=true;
        } else {
          options["fixaspectratio"]=false;
        }
      } catch(err) {
        options["fixaspectratio"]=false;
      }

      var xintegral = null;
      try {
        xintegral = plot.getElementsByTagName("xintegral")[0].childNodes[0].nodeValue;
        if(xintegral==="1"||xintegral==="true"){
          options["xintegral"]=true;
        } else {
          options["xintegral"]=false;
        }
      } catch(err) {
        options["xintegral"]=false;
      }

      var yintegral = null;
      try {
        yintegral = plot.getElementsByTagName("yintegral")[0].childNodes[0].nodeValue;
        if(yintegral==="1"||yintegral==="true"){
          options["yintegral"]=true;
        } else {
          options["yintegral"]=false;
        }
      } catch(err) {
        options["yintegral"]=false;
      }

      var hists = plot.getElementsByTagName("histogram");
      for(var ihist=0;ihist<hists.length;ihist++){
        var hist = hists[ihist];
        var col = parseInt(hist.getAttribute("col"))-1;
        var label;
        try {
          label = plotline.getElementsByTagName("label")[0].childNodes[0].nodeValue;
        } catch(err) {
          label = headers[col];
        }
        options["yintegral"]=true;
        var nbins;
        try {
          nbins = parseInt(hist.getElementsByTagName("nbins")[0].childNodes[0].nodeValue);
        } catch(err) {
          nbins = parseInt(Math.sqrt(data.length))+1;
        }

        var maxVal = Number.MIN_VALUE;
        var minVal = Number.MAX_VALUE;
        for(var idatarow=0;idatarow<data.length;idatarow++){
          if(dataColRow[col][idatarow]>maxVal){
            maxVal = dataColRow[col][idatarow];
          }
          if(dataColRow[col][idatarow]<minVal){
            minVal = dataColRow[col][idatarow];
          }
        }
        var binSize = (maxVal - minVal) / nbins;

        var d1 = [];

        var histData = {};
        for(var ibin=0;ibin<nbins;ibin++){
          var binBot = minVal + binSize * (ibin);
          var binMid = minVal + binSize * (ibin+0.5);
          var binTop = minVal + binSize * (ibin+1);
          histData[binMid]=0;
          for(var idatarow=0;idatarow<data.length;idatarow++){
            if((ibin==0&&(dataColRow[col][idatarow]>=binBot&&dataColRow[col][idatarow]<=binTop))||(ibin>0&&(dataColRow[col][idatarow]>binBot&&dataColRow[col][idatarow]<=binTop))){
              histData[binMid]++;
            }
          }
        }
        for(var key in histData){
          if(histData.hasOwnProperty(key)){
            d1.push([key,histData[key]])
          }
        }
        var theData = {data:d1,label:label};
        try {
          var colour = hist.getElementsByTagName("colour")[0].childNodes[0].nodeValue;
          if(colour=="r") colour = "red";
          if(colour=="g") colour = "green";
          if(colour=="b") colour = "blue";
          if(colour=="y") colour = "yellow";
          if(colour=="m") colour = "magenta";
          if(colour=="c") colour = "cyan";
          if(colour=="w") colour = "white";
          if(colour=="b") colour = "black";
          theData["color"] = colour;
        } catch(err) {
        }
        theData["bars"] = {show:true, barWidth:binSize*.8, align: "center"}
        d.push(theData);
      }

      var hists = plot.getElementsByTagName("barchart");
      var barWidth = 5;
      for(var ihist=0;ihist<hists.length;ihist++){
        var hist = hists[ihist];
        var col = parseInt(hist.getAttribute("col"))-1;
        var tcol = null;
        var haveTCol = false;
        if(hist.getAttribute("tcol")&&hist.getAttribute("tcol").length>0){
          haveTCol = true;
          tcol = parseInt(hist.getAttribute("tcol"))-1;
        }
        var barplacement;
        try {
          barplacement = plot.getElementsByTagName("barplacement")[0].childNodes[0].nodeValue;
        } catch(err) {
          barplacement = null;
        }
        var label;
        try {
          label = plot.getElementsByTagName("label")[0].childNodes[0].nodeValue;
        } catch(err) {
          label = headers[col];
        }
        options["yintegral"]=true;

        var d1 = [];
        var histData = {};
        for(var idatarow=0;idatarow<data.length;idatarow++){
          if(dataColRow[col][idatarow] in histData){
            if(haveTCol){
              histData[dataColRow[col][idatarow]] += parseInt(dataColRow[tcol][idatarow]);
            } else {
              histData[dataColRow[col][idatarow]]++;
	    }
          } else {
            if(haveTCol){
              histData[dataColRow[col][idatarow]] = parseInt(dataColRow[tcol][idatarow]);
            } else {
              histData[dataColRow[col][idatarow]] = 1;
	    }
          }
        }
        for(var key in histData){
          if(histData.hasOwnProperty(key)){
            d1.push([parseFloat(key),histData[key]])
          }
        }
        d1.sort(function(a, b){return a[0]-b[0]});

        var theData = {data:d1,label:label};
        try {
          var colour = hist.getElementsByTagName("colour")[0].childNodes[0].nodeValue;
          if(colour=="r") colour = "red";
          if(colour=="g") colour = "green";
          if(colour=="b") colour = "blue";
          if(colour=="y") colour = "yellow";
          if(colour=="m") colour = "magenta";
          if(colour=="c") colour = "cyan";
          if(colour=="w") colour = "white";
          if(colour=="b") colour = "black";
          theData["color"] = colour;
        } catch(err) {
        }
        try {
          if(hist.getElementsByTagName("width").length>0){
            barWidth = parseInt(hist.getElementsByTagName("width")[0].childNodes[0].nodeValue);
            //console.log(hist.getElementsByTagName("width"));
            //console.log(hist.getElementsByTagName("width")[0].childNodes[0].nodeValue);
          }
        } catch(err) {
        }
        //console.log(barWidth);
        theData["bars"] = {show:true, barWidth:barWidth, align: "center"}
        if(barplacement==="stacked"){
          theData["series"] = {
            stack: true,
            bars: {
              show: true, align: "center",
            },
          }
        }else if(barplacement==="sidebyside"){
          theData["bars"] = {show:true, order:ihist, align: "center"}
        }
        d.push(theData);
      }

      var yscale = null;
      try {
        yscale = plot.getElementsByTagName("yscale")[0].childNodes[0].nodeValue;
        if(yscale==="log"){
          options["yscale"]=yscale;
        }
      } catch(err) {
      }

      var xscale = null;
      try {
        xscale = plot.getElementsByTagName("xscale")[0].childNodes[0].nodeValue;
        if(xscale==="oneoversqrt"||xscale==="log"){
          options["xscale"]=xscale;
        }
      } catch(err) {
      }

      var yaxislabel = null;
      var ryaxislabel = null;
      var ylabels = plot.getElementsByTagName("ylabel");
      for(var iylab=0;iylab<ylabels.length;iylab++){
        try {
          var ylab = ylabels[iylab]
          var labrightaxis = ylab.getAttribute("rightaxis");
          if(!labrightaxis||labrightaxis=="false"||labrightaxis==0){
            yaxislabel = ylab.childNodes[0].nodeValue.replace("<","&lt;").replace(">","&gt;");;
          } else {
            ryaxislabel = ylab.childNodes[0].nodeValue.replace("<","&lt;").replace(">","&gt;");;
          }
        } catch(err) {
        }
      }

      titles = {title:plot_title, xlabel:xaxislabel,  ylabel:yaxislabel, rylabel:ryaxislabel};
      var background = null;
      try {
        background = plot.getElementsByTagName("background")[0].childNodes[0].nodeValue;
        console.log(background);
      } catch(err) {
      }

      var customXLabels = null;
      try {
        var customTicksLabelText = plot.getElementsByTagName("customXLabels")[0].childNodes[0].nodeValue;
        if(customTicksLabelText.split(",").length>0) {
            var splitLabelTicks = customTicksLabelText.split(",");
            customXLabels = [];
            for(var iTick=0;iTick<splitLabelTicks.length;iTick++){
                customXLabels.push(parseFloat(splitLabelTicks[iTick]));
            }
        }
      } catch(err) {
      }

      //console.log(customXLabels);

      var customXTicks = null;
      try {
        var customTicksText = plot.getElementsByTagName("customXTicks")[0].childNodes[0].nodeValue;

        if(customTicksText==="oneperdatapoint"){
            if(d.length>0 && "data" in d[0]){
                customXTicks = [];
                for(var idp=0;idp<d[0]["data"].length;idp++){
                    customXTicks.push(d[0]["data"][idp][0]);
                }
            }
        } else if(customTicksText.split(",").length>0) {
            var splitTicks = customTicksText.split(",");
            customXTicks = [];
            if(customXLabels&&customXLabels.length===splitTicks.length){
                for(var iTick=0;iTick<splitTicks.length;iTick++){
                    customXTicks.push([parseFloat(splitTicks[iTick]),customXLabels[iTick]]);
                }
            } else {
                for(var iTick=0;iTick<splitTicks.length;iTick++){
                    customXTicks.push(parseFloat(splitTicks[iTick]));
                }
            }
        }
      } catch(err) {
      }

      var graphName = graphObject.graphDivName+igraph;
      var h = document.getElementById(graphObject.graphDivName).style.height;
      thisPlotDiv =  document.createElement("div");
      thisPlotDiv.setAttribute("id",graphName);
      thisPlotDiv.setAttribute("style","width:100%;height:"+h+";display:none;");
      pimpleplot.appendChild(thisPlotDiv);
      //$.plot("#"+graphName, d);

      var selectItem = document.createElement("option");
      if(iplot==plots.length-1&&itab==tables.length-1){
        selectItem.setAttribute("selected","true");
      }
      selectItem.innerHTML = plot_title.replace("<","&lt;").replace(">","&gt;");
      selectItem.setAttribute("value",graphName+"_selectitem");
      selectItem.setAttribute("id",graphName+"_selectitem");
      selectGroup.appendChild(selectItem)

      graph = new CCP4Graph(graphName,d,xbreak,ybreak,options,shapes,titles,background,customXTicks);

      thisPlotDiv.iplot = iplot;
      thisPlotDiv.itab = itab;
      thisPlotDiv.itab = itab;
      thisPlotDiv.graphDivName_in = graphObject.graphDivName_in;

      thisPlotDiv.addEventListener("plotHover",function(e) {
        pimpleMain = document.getElementById(this.graphDivName_in)
        //console.log("Graph hovered!!!"+e.detail.x+" "+e.detail.y + ". Plot " + this.iplot + " of table " + this.itab)
        pimpleMain.hoverevt = document.createEvent("CustomEvent");
        pimpleMain.hoverevt.initCustomEvent('graphHover', false, false, { 'x': e.detail.x, 'y': e.detail.y, 'plot': this.iplot, 'table': this.itab});
        pimpleMain.dispatchEvent(pimpleMain.hoverevt);
      })
      thisPlotDiv.addEventListener("plotClick",function(e) {
        pimpleMain = document.getElementById(this.graphDivName_in)
        //console.log("Graph clicked!!!"+e.detail.x+" "+e.detail.y + ". Plot " + this.iplot + " of table " + this.itab)
        pimpleMain.clickevt = document.createEvent("CustomEvent");
        pimpleMain.clickevt.initCustomEvent('graphClick', false, false, { 'x': e.detail.x, 'y': e.detail.y, 'plot': this.iplot, 'table': this.itab});
        pimpleMain.dispatchEvent(pimpleMain.clickevt);
      })

      graphObject.graphs.push(graph);

      igraph++;

    }
    //console.log(" ");
  }

  //console.log(graphName)

  if(this.currentStoredGraph>-1){
    $("#"+graphObject.menuSelectName+" option:eq("+graphObject.currentStoredGraph+")").prop('selected', true);
    var graphName = $("#"+graphObject.menuSelectName)[0].value.replace("_selectitem","");
    graph = graphObject.graphs[graphObject.currentStoredGraph];
    thisPlotDiv = document.getElementById(graphName);

    $('#'+graphObject.graphDivName).children().hide();
    thisPlotDiv.setAttribute("style","width:100%;height:"+h+";display:block;");
    $("#"+graphName).empty();
    graphObject.currentStoredGraph = -1;
  }

  if(typeof thisPlotDiv!=="undefined"){
    graphObject.currentGraph = graph;
    thisPlotDiv.setAttribute("style","width:100%;height:"+h+";display:block;");
    graphObject.plot();
 }

}

CCP4GraphPlot.prototype.drawDashedLine= function(ctx,posLeft_in,posTop_in,posLeft2,posTop2,lineStyle) {

  //console.log(lineStyle)
  //console.log(lineStyle.length)
  var posLeft = 0 + posLeft_in;
  var posTop = 0 + posTop_in;
  var posLeft_orig = 0 + posLeft_in;
  var posTop_orig = 0 + posTop_in;

  var dashLength = 2;
  var gapLength = 2;

  if(lineStyle.length==2){
    dashLength = lineStyle[0];
    gapLength = lineStyle[1];
  }

  if(lineStyle.length==4){
    dashLength = lineStyle[0];
    gapLength = lineStyle[1]+lineStyle[2]+lineStyle[3];
  }

  ctx.moveTo(posLeft, posTop);

  var dX = posLeft2 - posLeft;
  var dY = posTop2 - posTop;
  var dashes = Math.floor(Math.sqrt(dX * dX + dY * dY) / (dashLength+gapLength));
  var dashX = dX / dashes * dashLength/(dashLength+gapLength);
  var dashY = dY / dashes * dashLength/(dashLength+gapLength);
  var gapX = dX / dashes * gapLength/(dashLength+gapLength);
  var gapY = dY / dashes * gapLength/(dashLength+gapLength);
    //console.log(dashLength+" "+gapLength)

  for(var idash=0;idash<dashes;idash++){
    posLeft += dashX;
    posTop += dashY;
    ctx.lineTo(posLeft,posTop);
    posLeft += gapX;
    posTop += gapY;
    ctx.moveTo(posLeft,posTop);
  }

  if(lineStyle.length==4){
    posLeft = posLeft_orig;
    posTop = posTop_orig;
    ctx.moveTo(posLeft, posTop);
    var dashLength2 = lineStyle[2];
    var gapLength2 = lineStyle[1]+lineStyle[0]+lineStyle[3];
    var dX = posLeft2 - posLeft;
    var dY = posTop2 - posTop;
    var dashes = Math.floor(Math.sqrt(dX * dX + dY * dY) / (dashLength2+gapLength2));
    var dashX = dX / dashes * dashLength2/(dashLength2+gapLength2);
    var dashY = dY / dashes * dashLength2/(dashLength2+gapLength2);
    var gapX = dX / dashes * gapLength2/(dashLength2+gapLength2);
    var gapY = dY / dashes * gapLength2/(dashLength2+gapLength2);
    //console.log(dashLength2+" "+gapLength2)

    var moveX = lineStyle[0]+lineStyle[1];
    var dashesMove = Math.floor(Math.sqrt(dX * dX + dY * dY) / (moveX));
    var dashMoveX = dX / dashesMove;
    var dashMoveY = dY / dashesMove;
    posLeft += dashMoveX;
    posTop += dashMoveY;
    ctx.moveTo(posLeft, posTop);

    for(var idash=0;idash<dashes;idash++){
      posLeft += dashX;
      posTop += dashY;
      ctx.lineTo(posLeft,posTop);
      posLeft += gapX;
      posTop += gapY;
      ctx.moveTo(posLeft,posTop);
    }

  }
}

CCP4GraphPlot.prototype.p2cNY = function(thePlot,pos,n) {
            var res = {}, i, axis, key;
            for (i = 0; i < thePlot.getYAxes().length; ++i) {
                axis = thePlot.getYAxes()[i];
                if (axis && axis.used) {
                    key = "y" + axis.n;
                    if (axis.n == n)
                        key = "y";

                    if (pos[key] != null) {
                        res.top = axis.p2c(pos[key]);
                        break;
                    }
                }
            }
            return res;
        }

CCP4GraphPlot.prototype.drawBackground = function(thePlot,divTop) {
    var background = this.currentGraph.background;
    if(background){
      var graphName = this.currentGraph.graphID;
      var canvas = thePlot.getCanvas();
      var offset = thePlot.getPlotOffset();
      var ctx = canvas.getContext('2d');
      var outer = document.getElementById(graphName);
      var offset = thePlot.getPlotOffset();
      var h = parseInt($(divTop).height())-offset["left"]-offset["right"];
      var w = parseInt($(divTop).width())-offset["top"]-offset["bottom"];
      if (canvas.getContext){
        drawing = new Image();
        drawing.src = background;
        drawing.onload = function() {
            ctx.save();
            ctx.globalCompositeOperation='destination-over';
            ctx.drawImage(drawing,offset["left"],offset["top"],w,h);
            console.log("Drawn image"+graphName);
            ctx.restore();
        };
      }
    }
}

CCP4GraphPlot.prototype.plotShapes = function(thePlot,yaxisOverride) {
      var shapes = this.currentGraph.shapes;
      var canvas = thePlot.getCanvas();
      var offset = thePlot.getPlotOffset();
      var ctx = canvas.getContext('2d');
      var haveDashedLines = true;
      if (!ctx.setLineDash) {
        ctx.setLineDash = function (dash) {
          if(!ctx.webkitLineDash){
            if(!ctx.mozDash){
              console.log("No dashed lines.")
              haveDashedLines = false;
            } else {
              console.log("Using mozDash")
              ctx.mozDash = dash; // This may have no effect.
            }
          } else {
            console.log("Using webkitLineDash")
            ctx.webkitLineDash = dash; // This may have no effect.
          }
        }
      }

      for(var ishape=0;ishape<shapes.length;ishape++){
        var plotOffset = thePlot.getPlotOffset();
        var canvasCoordDum1 = thePlot.p2c({x:0,y:0})
        var canvasCoordDum2 = thePlot.p2c({x:100,y:100})
        //for(iax=0;iax<thePlot.getYAxes().length;iax++){
          //console.log(iax+" "+thePlot.getYAxes()[iax].p2c(100))
        //}
        if(yaxisOverride!=1){
               var canvasCoordy = this.p2cNY(thePlot,{y:0},yaxisOverride)
               canvasCoordDum1["top"] = canvasCoordy["top"]
               var canvasCoordy = this.p2cNY(thePlot,{y:100},yaxisOverride)
               canvasCoordDum2["top"] = canvasCoordy["top"]
        }
        var xscale = (parseFloat(canvasCoordDum2["left"]) - parseFloat(canvasCoordDum1["left"]))*.01
        var yscale = -(parseFloat(canvasCoordDum2["top"]) - parseFloat(canvasCoordDum1["top"]))*.01

        var shape = shapes[ishape];

        if("text" in shape){
/*

I guess "size" can be approximated by:
lineHeight=ctx.measureText('M').width; 


var underline = function(ctx, text, x, y, size, color, thickness ,offset){
  var width = ctx.measureText(text).width;

  switch(ctx.textAlign){
    case "center":
    x -= (width/2); break;
    case "right":
    x -= width; break;
  }

  y += size+offset;

  ctx.beginPath();
  ctx.strokeStyle = color;
  ctx.lineWidth = thickness;
  ctx.moveTo(x,y);
  ctx.lineTo(x+width,y);
  ctx.stroke();

}
*/
          var text = shape["text"];
          var font = text["font"];
          var fillColour = text["fillcolour"];
          if (canvas.getContext){
            ctx.save();
            var canvasCoord = thePlot.p2c({x:text["x"],y:text["y"]})
            var posLeft = canvasCoord["left"]+plotOffset["left"];
            var posTop = canvasCoord["top"]+plotOffset["top"];
            if(font !== ""){
              ctx.font = font;
            }
            ctx.fillStyle = fillColour;
            ctx.fillText(text["text"],posLeft,posTop);
            ctx.restore();
          }
        }

        if("linestrip" in shape){
          var linestrip = shape["linestrip"];
          var lineColour = linestrip["linecolour"];
          var lineSize = linestrip["linesize"];
          var lineStyle = linestrip["linestyle"];
          if (canvas.getContext){
            ctx.save();
            ctx.beginPath();
            ctx.strokeStyle = lineColour;
            ctx.lineWidth = lineSize;
            if(lineStyle) ctx.setLineDash(lineStyle)
            var canvasCoord = thePlot.p2c({x:linestrip["data"][0][0],y:linestrip["data"][0][1]})
            if(yaxisOverride!=1){
               var canvasCoordy = this.p2cNY(thePlot,{y:linestrip["data"][0][1]},yaxisOverride)
               canvasCoord["top"] = canvasCoordy["top"]
            }
            var posLeft = canvasCoord["left"]+plotOffset["left"];
            var posTop = canvasCoord["top"]+plotOffset["top"];
            ctx.moveTo(posLeft,posTop);
            for(var ip=1;ip<linestrip["data"].length;ip++){
              var canvasCoord = thePlot.p2c({x:linestrip["data"][ip][0],y:linestrip["data"][ip][1]})
              if(yaxisOverride!=1){
                var canvasCoordy = this.p2cNY(thePlot,{y:linestrip["data"][ip][1]},yaxisOverride)
                canvasCoord["top"] = canvasCoordy["top"]
              }
              var posLeft2 = canvasCoord["left"]+plotOffset["left"];
              var posTop2 = canvasCoord["top"]+plotOffset["top"];
              //console.log(posLeft,posTop,posLeft2,posTop2)
              if(!lineStyle||(lineStyle&&haveDashedLines)){
                ctx.lineTo(posLeft2,posTop2);
              }else if(lineStyle&&!haveDashedLines){
                this.drawDashedLine(ctx,posLeft,posTop,posLeft2,posTop2,lineStyle);
                posLeft = posLeft2;
                posTop = posTop2;
              }
            }
            ctx.stroke();
            ctx.closePath();
            ctx.restore();
          }
        }

        if("polygon" in shape){
          var polygon = shape["polygon"];
          var fillColour = polygon["fillcolour"];
          var lineColour = polygon["linecolour"];
          var lineSize = polygon["linesize"];
          var lineStyle = polygon["linestyle"];
          if (canvas.getContext){
            ctx.save();
            ctx.beginPath();
            ctx.strokeStyle = lineColour;
            ctx.fillStyle = fillColour;
            ctx.lineWidth = lineSize;
            if(lineStyle) ctx.setLineDash(lineStyle)
            var canvasCoord = thePlot.p2c({x:polygon["data"][0][0],y:polygon["data"][0][1]})
            if(yaxisOverride!=1){
              var canvasCoordy = this.p2cNY(thePlot,{y:polygon["data"][0][1]},yaxisOverride)
              canvasCoord["top"] = canvasCoordy["top"]
            }
            var posLeft = canvasCoord["left"]+plotOffset["left"];
            var posTop = canvasCoord["top"]+plotOffset["top"];
            var posLeft_orig = canvasCoord["left"]+plotOffset["left"];
            var posTop_orig = canvasCoord["top"]+plotOffset["top"];
            var posLeft2;
            var posTop2;
            ctx.moveTo(posLeft,posTop);
            for(var ip=1;ip<polygon["data"].length;ip++){
              var canvasCoord = thePlot.p2c({x:polygon["data"][ip][0],y:polygon["data"][ip][1]})
              if(yaxisOverride!=1){
                var canvasCoordy = this.p2cNY(thePlot,{y:polygon["data"][ip][1]},yaxisOverride)
                canvasCoord["top"] = canvasCoordy["top"]
              }
              posLeft2 = canvasCoord["left"]+plotOffset["left"];
              posTop2 = canvasCoord["top"]+plotOffset["top"];
              if(!lineStyle||(lineStyle&&haveDashedLines)){
                ctx.lineTo(posLeft2,posTop2);
              }else if(lineStyle&&!haveDashedLines){
                this.drawDashedLine(ctx,posLeft,posTop,posLeft2,posTop2,lineStyle);
                posLeft = posLeft2;
                posTop = posTop2;
              }
            }
            if(!lineStyle||(lineStyle&&haveDashedLines)){
              ctx.lineTo(posLeft,posTop);
            }else if(lineStyle&&!haveDashedLines){
              this.drawDashedLine(ctx,posLeft2,posTop2,posLeft_orig,posTop_orig,lineStyle);
              posLeft = posLeft2;
              posTop = posTop2;
            }
            ctx.stroke();
            ctx.closePath();
            ctx.beginPath();
            var canvasCoord = thePlot.p2c({x:polygon["data"][0][0],y:polygon["data"][0][1]})
            if(yaxisOverride!=1){
              var canvasCoordy = this.p2cNY(thePlot,{y:polygon["data"][0][1]},yaxisOverride)
              canvasCoord["top"] = canvasCoordy["top"]
            }
            var posLeft = canvasCoord["left"]+plotOffset["left"];
            var posTop = canvasCoord["top"]+plotOffset["top"];
            ctx.moveTo(posLeft,posTop);
            for(var ip=1;ip<polygon["data"].length;ip++){
              var canvasCoord = thePlot.p2c({x:polygon["data"][ip][0],y:polygon["data"][ip][1]})
              if(yaxisOverride!=1){
                var canvasCoordy = this.p2cNY(thePlot,{y:polygon["data"][ip][1]},yaxisOverride)
                canvasCoord["top"] = canvasCoordy["top"]
              }
              var posLeft2 = canvasCoord["left"]+plotOffset["left"];
              var posTop2 = canvasCoord["top"]+plotOffset["top"];
              ctx.lineTo(posLeft2,posTop2);
            }
            ctx.lineTo(posLeft,posTop);
            ctx.fill();
            ctx.restore();
          }
        }

        if("circle" in shape){
          var circle = shape["circle"];

          var canvasCoord = thePlot.p2c({x:circle["x"],y:circle["y"]})
          if(yaxisOverride!=1){
              var canvasCoordy = this.p2cNY(thePlot,{y:circle["y"]},yaxisOverride)
              canvasCoord["top"] = canvasCoordy["top"]
          }
          //console.log(canvasCoord)
          var radius = parseInt(Math.floor(xscale * circle["radius"]))
          var fillColour = circle["fillcolour"];
          var lineColour = circle["linecolour"];
          var lineSize = circle["linesize"];
          var lineStyle = circle["linestyle"];
      
          var posLeft = canvasCoord["left"]+plotOffset["left"];
          var posTop = canvasCoord["top"]+plotOffset["top"];
          if (canvas.getContext){
            ctx.save();
            ctx.beginPath();
            ctx.strokeStyle = lineColour;
            ctx.lineWidth = lineSize;
            if(lineStyle) ctx.setLineDash(lineStyle)
            ctx.arc(posLeft,posTop,radius,0,Math.PI*2);
            ctx.fillStyle = fillColour;
            ctx.fill();
            ctx.beginPath();
            ctx.arc(posLeft,posTop,radius,0,Math.PI*2);
            ctx.stroke();
            ctx.closePath();
            ctx.restore();
          }
       }
    }
}

CCP4GraphPlot.prototype.plot = function() {
    
    var xbreak = this.currentGraph.xbreaks;
    var ybreak = this.currentGraph.ybreaks;
    var options = this.currentGraph.options;
    var shapes = this.currentGraph.shapes;
    var graphName = this.currentGraph.graphID;
    var d = this.currentGraph.graphData;
    var titles = this.currentGraph.titles;

    var customXTicks = this.currentGraph.customXTicks;

    var theRealOuterDiv = document.getElementById(graphName);

    var N = 20;

    var xlabel;
    var ylabel;
    var rylabel;

    if(titles){
      xlabel = titles["xlabel"];
      ylabel = titles["ylabel"];
      rylabel = titles["rylabel"];
    }

    var yaxes;
    var yaxes2;
    var yr;
    var yr2;

    if(typeof xbreak!=="undefined"&&xbreak.length>0){

      var yax = [];
      var xax = [];
      if("yrange" in options){
        yr = jQuery.extend(true, {},options["yrange"])
        yr2 = jQuery.extend(true, {},options["yrange"])
        yaxes = [yr]
        yaxes2 = [yr2]
        for(var ibr=0;ibr<=(xbreak.length/2);ibr++){
          yax.push([yr])
        }
      } else {
        yaxes = [{}]
        yr2 = {}
        yaxes2 = [yr2]
      }

      if("ryrange" in options){
        yr = jQuery.extend(true, {},options["ryrange"])
        yr2 = jQuery.extend(true, {},options["ryrange"])
        yr["position"] = "right"
        yaxes.push(yr)
        yr2["position"] = "right"
        yaxes2.push(yr2)
        for(var ibr=0;ibr<=(xbreak.length/2);ibr++){
          yax[ibr].push(yr)
        }
      } else {
        yr = {}
        yaxes.push(yr)
        yr2 = {}
        yr2["position"] = "right"
        yaxes2.push(yr2)
      }
      //console.log(yax)


      if("xrange" in options){
        var xr = jQuery.extend(true, {},options["xrange"])
        var xr2 = jQuery.extend(true, {},options["xrange"])
        xr["max"] =  xbreak[0]
        xr2["min"] =  xbreak[1]
        xax.push(xr);
        for(var ibr=0;ibr<(xbreak.length)/2-1;ibr++){
          xax.push({min:xbreak[2*ibr+1],max:xbreak[2*ibr+2]})
        }
        var xr2p = jQuery.extend(true, {},xr2)
        xr2p["min"] = xbreak[xbreak.length-1]
        xax.push(xr2p)
        xaxes = [xr]
      } else {
        //FIXME - this looks iffy. But we never get here.
        //console.log("Setting other way ....")
        xaxes = [xbreak[0]]
      }

      //console.log(xax)

    var outer = document.getElementById(graphName);
    
    var hTop = parseInt(outer.style.height);
    var bs = [];
    var divTop;
    var pcwidth = parseInt(100./xax.length)+"%";
    for(var ibr=0;ibr<xax.length;ibr++){
      var bDName = new Array(N+1).join((Math.random().toString(36)+'00000000000000000').slice(2, 18)).slice(0, N);
      var b =  document.createElement("div");
      b.setAttribute("style","width:"+pcwidth+";height:"+outer.style.height+";float:left");
      b.setAttribute("id",bDName);
      bs.push(b);
    }

    if(xlabel){
      divTop =  document.createElement("div");
      var divBot =  document.createElement("div");
      divBot.setAttribute("class","xlabel");
      outer.appendChild(divTop)
      outer.appendChild(divBot)
      var h = parseInt($(divBot).height());
      if(!isNaN(h)){
        hTop = parseInt(outer.style.height) - h;
        divTop.setAttribute("style","width:100%;height:"+(hTop-1)+"px;float:top;");
        for(var ib=0;ib<bs.length;ib++){
          bs[ib].setAttribute("style","width:"+pcwidth+";height:"+hTop+"px;float:left;");
        }
      }
      divBot.innerHTML = xlabel.trim().replace(/\s+/,'&nbsp;')
    } else {
      divTop = outer;
      divTop.setAttribute("style","width:100%;height:"+(hTop-1)+"px;float:top;");
    }

    for(var ib=0;ib<bs.length;ib++){
      divTop.appendChild(bs[ib])
    }

    if(ylabel){
      outer = divTop;
      divTop.setAttribute("class","calculated-width");
      divTop =  document.createElement("div");
      outer.appendChild(divTop)
      var divLeft = document.createElement("div");
      divLeft.setAttribute("class","left ylabel"); // And ylabel!
      divLeft.setAttribute("style","float:left;");
      divLeft.innerHTML = ylabel.trim().replace(/\s+/,'&nbsp;')
      divTop.appendChild(divLeft)
      if(rylabel){
        outer.setAttribute("class","calculated-width-left-right");
        var divRight = document.createElement("div");
        divRight.setAttribute("class","right ylabel rylabel");
        divRight.setAttribute("style","float:right;");
        divRight.innerHTML = rylabel.trim().replace(/\s+/,'&nbsp;')
        var w = parseInt($(outer).width());
        divRight.style.left=""+(w+14)+"px";
        var h = parseInt($(outer).height());
        divRight.style.width=""+(h-20)+"px";
        outer.appendChild(divRight)
      }
    }

    if(rylabel&&!ylabel){
      outer = divTop;
      outer.setAttribute("class","calculated-width-right");
      divTop =  document.createElement("div");
      outer.appendChild(divTop)
      var divRight = document.createElement("div");
      divRight.setAttribute("class","right ylabel rylabel");
      divRight.setAttribute("style","float:right;");
      divRight.innerHTML = rylabel.trim().replace(/\s+/,'&nbsp;')
      var w = parseInt($(outer).width());
      divRight.style.left=""+(w-4)+"px";
      var h = parseInt($(outer).height());
      divRight.style.width=""+(h-20)+"px";
      outer.appendChild(divRight)
    }

    for(var ib=0;ib<bs.length;ib++){
      var dRight = jQuery.extend(true, [], d);
      for(var idata=0;idata<dRight.length;idata++){
        if(!dRight[idata]["yaxis"]){
          if(ib==0){
            dRight[idata]["yaxis"] = 1;//1+2*ib;
          } else {
            dRight[idata]["yaxis"] = 3;//1+2*ib;
          }
        }else{
          if(ib==bs.length-1){
            dRight[idata]["yaxis"] = 2;//2+2*ib;
          } else {
            dRight[idata]["yaxis"] = 4;//2+2*ib;
          }
        }
      }
      yaxes2 = [{},{},yaxes2[0],yaxes2[1]]
      //console.log(yaxes2)
      //console.log([{},{},yax[ib][0],yax[ib][1]])
      var bDName = bs[ib].getAttribute("id");
      //console.log(bDName);
      //console.log(xax[ib]);
      //console.log(yax[ib]);
      var legend = {};
      if(options["showlegend"]===false){
        legend = {show:false};
      }
      if(ib<bs.length-1){
        var xaxes = [xax[ib]];
        if(options["xscale"]==="oneoversqrt"){
          xaxes[0]["tickFormatter"] = function(x) { var y = 1.0/Math.sqrt(x); if(isFinite(y)){ return (y).toFixed(2); }else{return "Inf";} }
        }
        if(options["xscale"]==="log"){
          xaxes[0]["transform"] = function(v) { return v == 0 ? null : Math.log(v) };
          xaxes[0]["inverseTransform"] = function (v) { return Math.exp(v); }
          xaxes[0]["ticks"] = 4;
        }
        if(options["yscale"]==="log"){
          yax[ib][0]["transform"] = function(v) { return v == 0 ? null : Math.log(v) };
          yax[ib][0]["inverseTransform"] = function (v) { return Math.exp(v); }
          yax[ib][0]["ticks"] = 4;
          yax[ib][1]["transform"] = function(v) { return v == 0 ? null : Math.log(v) };
          yax[ib][1]["inverseTransform"] = function (v) { return Math.exp(v); }
          yax[ib][1]["ticks"] = 4;
        }
        if(options["xintegral"]===true){
          xaxes[0]["tickDecimals"] = 0;
        }
        if(options["yintegral"]===true){
          yaxes[ib][0]["tickDecimals"] = 0;
          yaxes[ib][1]["tickDecimals"] = 0;
        }
        if(customXTicks){
            xaxes[0]["ticks"] = customXTicks;
        }
        var thePlot2 = $.plot("#"+bDName, dRight,{
          xaxes: xaxes,
          yaxes : [yax[ib][0],yax[ib][1],yax[ib][0],yax[ib][1],yax[ib][0],yax[ib][1]],
          grid: {borderWidth: 0, clickable: true, hoverable: true},
          legend: {show:false},
        });
      } else {
        var xaxes = [xax[ib]];
        if(options["xscale"]==="oneoversqrt"){
          xaxes[0]["tickFormatter"] = function(x) { var y = 1.0/Math.sqrt(x); if(isFinite(y)){ return (1.0/Math.sqrt(x)).toFixed(2); }else{return "Inf";} }
        }
        if(options["xscale"]==="log"){
          xaxes[0]["transform"] = function(v) { return v == 0 ? null : Math.log(v) };
          xaxes[0]["inverseTransform"] = function (v) { return Math.exp(v); }
          xaxes[0]["ticks"] = 4;
        }
        if(options["yscale"]==="log"){
          yax[ib][0]["transform"] = function(v) { return v == 0 ? null : Math.log(v) };
          yax[ib][0]["inverseTransform"] = function (v) { return Math.exp(v); }
          yax[ib][0]["ticks"] = 4;
          yax[ib][1]["transform"] = function(v) { return v == 0 ? null : Math.log(v) };
          yax[ib][1]["inverseTransform"] = function (v) { return Math.exp(v); }
          yax[ib][1]["ticks"] = 4;
        }
        if(options["xintegral"]===true){
          xaxes[0]["tickDecimals"] = 0;
        }
        if(options["yintegral"]===true){
          yaxes[ib][0]["tickDecimals"] = 0;
          yaxes[ib][1]["tickDecimals"] = 0;
        }
        if(customXTicks){
            xaxes[0]["ticks"] = customXTicks;
        }
        var thePlot2 = $.plot("#"+bDName, dRight,{
          xaxes: xaxes,
          yaxes : [yax[ib][0],yax[ib][1],yax[ib][0],yax[ib][1],yax[ib][0],yax[ib][1]],
          grid: {borderWidth: 0, clickable: true, hoverable: true},
          legend: legend,
        });
      }
      if(ib==0){
        this.plotShapes(thePlot2,1);
      } else {
        this.plotShapes(thePlot2,3);
      }

      $("#"+bDName).bind("plotclick", function (event, pos, item) {
        if (item) {
          //console.log("You clicked at " + item["datapoint"][0] + ", " + item["datapoint"][1]);
          theRealOuterDiv.clickevtone = document.createEvent("CustomEvent");
          theRealOuterDiv.clickevtone.initCustomEvent('plotClick', false, false, { 'x': item["datapoint"][0], 'y': item["datapoint"][1] });
          theRealOuterDiv.dispatchEvent(theRealOuterDiv.clickevtone);
        }
      });
      $("#"+bDName).bind("plothover", function (event, pos, item) {
        if (item) {
          //console.log("You hovered at " + item["datapoint"][0] + ", " + item["datapoint"][1]);
          theRealOuterDiv.hoverevtone = document.createEvent("CustomEvent");
          theRealOuterDiv.hoverevtone.initCustomEvent('plotHover', false, false, { 'x': item["datapoint"][0], 'y': item["datapoint"][1] });
          theRealOuterDiv.dispatchEvent(theRealOuterDiv.hoverevtone);
        }
      });
    }

    } else {
      if("yrange" in options){
        yaxes = [options["yrange"]]
      } else {
        yaxes = [{}]
      }
      if("ryrange" in options){
         yr = jQuery.extend(true, {},options["ryrange"])
         yr["position"] = "right"
         yr["color"] = "red"
         yaxes.push(yr)
      }
      var xaxes;
      if("xrange" in options){
        xaxes = [options["xrange"]]
      } else {
        xaxes = [{}]
      }
 
    var divTop;
    var outer  = document.getElementById(graphName);

 //  This is not quite right. It makes the outer div - containing axis labels and ticks - square.
 //  Don't really know what we can do aboout ticks: We need to know there size.
    if(options["fixaspectratio"]===true){
       var outer1 = document.getElementById(graphName);
       var h = parseInt($(outer1).height());
       var w = parseInt($(outer1).width());
       var theSize = Math.min(h,w);
       outer =  document.createElement("div");
       outer.setAttribute("style","width:"+theSize+"px;height:"+theSize+"px;float:top;");
       outer1.appendChild(outer);
    } else {
      outer = document.getElementById(graphName);
    }

    if(xlabel){
      divTop =  document.createElement("div");
      var divBot =  document.createElement("div");
      divBot.setAttribute("class","xlabel");
      outer.appendChild(divTop)
      outer.appendChild(divBot)
      var h = parseInt($(divBot).height());
      //console.log("divBot height:"+divBot.style.height+" "+h);
      if(!isNaN(h)){
        var hTop = parseInt(outer.style.height) - h;
        //console.log("computed top height: "+hTop);
        divTop.setAttribute("style","width:100%;height:"+hTop+"px;float:top;");
      }
      divBot.innerHTML = xlabel.trim().replace(/\s+/,'&nbsp;')
    } else {
      divTop = document.getElementById(graphName);
    }

    if(ylabel){
      outer = divTop;
      outer.setAttribute("class","calculated-width");
      divTop =  document.createElement("div");
      outer.appendChild(divTop)
      var divLeft = document.createElement("div");
      divLeft.setAttribute("class","left ylabel"); // And ylabel!
      divLeft.setAttribute("style","float:left;");
      divLeft.innerHTML = ylabel.trim().replace(/\s+/,'&nbsp;')
      outer.appendChild(divLeft)
      if(rylabel){
        outer.setAttribute("class","calculated-width-left-right");
        var divRight = document.createElement("div");
        divRight.setAttribute("class","right ylabel rylabel");
        divRight.setAttribute("style","float:right;");
        divRight.innerHTML = rylabel.trim().replace(/\s+/,'&nbsp;')
        var w = parseInt($(outer).width());
        divRight.style.left=""+(w+14)+"px";
        var h = parseInt($(outer).height());
        divRight.style.width=""+(h-20)+"px";
        outer.appendChild(divRight)
      }
      divTop.style.height="100%";
      divTop.style.width="100%";
    }
    if(rylabel&&!ylabel){
      outer = divTop;
      outer.setAttribute("class","calculated-width-right");
      divTop =  document.createElement("div");
      outer.appendChild(divTop)
      var divRight = document.createElement("div");
      divRight.setAttribute("class","right ylabel rylabel");
      divRight.setAttribute("style","float:right;");
      divRight.innerHTML = rylabel.trim().replace(/\s+/,'&nbsp;')
      var w = parseInt($(outer).width());
      divRight.style.left=""+(w-4)+"px";
      var h = parseInt($(outer).height());
      divRight.style.width=""+(h-20)+"px";
      outer.appendChild(divRight)
      divTop.style.height="100%";
      divTop.style.width="100%";
    }

    if(xaxes.length==1 && options["xscale"]==="oneoversqrt"){
      xaxes[0]["tickFormatter"] = function(x) { var y = 1.0/Math.sqrt(x); if(isFinite(y)){ return (1.0/Math.sqrt(x)).toFixed(2); }else{return "Inf";} }
    }
    if(options["xscale"]==="log"){
      xaxes[0]["transform"] = function(v) { return v == 0 ? null : Math.log(v) };
      xaxes[0]["inverseTransform"] = function (v) { return Math.exp(v); }
      xaxes[0]["ticks"] = 4;
    }
    if(options["yscale"]==="log"){
      yaxes[0]["transform"] = function(v) { return v == 0 ? null : Math.log(v) };
      yaxes[0]["inverseTransform"] = function (v) { return Math.exp(v); }
      yaxes[0]["ticks"] = 4;
    }
    if(options["xintegral"]===true){
      xaxes[0]["tickDecimals"] = 0;
    }
    if(options["yintegral"]===true){
      yaxes[0]["tickDecimals"] = 0;
    }
    var legend = {};
    if(options["showlegend"]===false){
       legend = {show:false};
    }

    if(options["fixaspectratio"]===true){
      var h = parseInt($(divTop).height());
      var w = parseInt($(divTop).width());
      divTop.style.height=Math.min(w,h)+"px";
      divTop.style.width=Math.min(w,h)+"px";
    }

    var start = new Date().getTime();
    if(customXTicks){
        xaxes[0]["ticks"] = customXTicks;
    }

    var thePlot = $.plot(divTop, d,{
        xaxes: xaxes, yaxes: yaxes,
        grid: {borderWidth: 0, clickable: true, hoverable: true},
        legend: legend,
        series: d[0]["series"],
    });

    var end = new Date().getTime();
    var time = end - start;
    //console.log('Execution time: ' + time*.001);

    this.plotShapes(thePlot,1);
    this.drawBackground(thePlot,divTop);

    $("#"+graphName).bind("plotclick", function (event, pos, item) {
       if (item) {
        //console.log("You clicked at(2) " + pos.x + ", " + pos.y);
        theRealOuterDiv.clickevttwo = document.createEvent("CustomEvent");
        theRealOuterDiv.clickevttwo.initCustomEvent('plotClick', false, false, { 'x': pos.x, 'y': pos.y });
        theRealOuterDiv.dispatchEvent(theRealOuterDiv.clickevttwo);
      }
    });
    $("#"+graphName).bind("plothover", function (event, pos, item) {
       if (item) {
        //console.log("You hovered at(2) " + pos.x + ", " + pos.y + " " + graphName);
        theRealOuterDiv.hoverevttwo = document.createEvent("CustomEvent");
        theRealOuterDiv.hoverevttwo.initCustomEvent('plotHover', false, false, { 'x': pos.x, 'y': pos.y });
        theRealOuterDiv.dispatchEvent(theRealOuterDiv.hoverevttwo);

      }
    });
    }

      if(this.fonts !== "undefined" && this.fonts !== null && this.fonts.length>0){
        var legendFonts = this.fonts[0].getElementsByTagName("legendFont");
        if(legendFonts !== "undefined" && legendFonts !== null && legendFonts.length>0){
          //console.log(legendFonts[0]);
          var font_size = legendFonts[0].getAttribute("size");
          var font_family = legendFonts[0].getAttribute("family");
          var font_weight = legendFonts[0].getAttribute("weight");
          var font_slant = legendFonts[0].getAttribute("slant");
          var font_style = "";
          if(font_size!==null&&font_size!=""){
            $("#"+graphName+" div.legend table").css("font-size", parseInt(font_size));
          }
          if(font_family!==null&&font_family!=""){
            $("#"+graphName+" div.legend table").css("font-family", font_family);
          }
          if(font_weight!==null&&font_weight!=""){
            $("#"+graphName+" div.legend table").css("font-weight", font_weight);
          }
          if(font_slant!==null&&font_slant!=""&&font_slant!="normal"){
            $("#"+graphName+" div.legend table").css("font-style", font_slant);
          }
          $("#"+graphName+" div.legend table").css("background-color", "rgba(255, 255, 255, 0.2)");
          $("#"+graphName+" div.legend div").css("background-color", "rgba(255, 255, 255, 0.0)"); // This is surely a bug in flot, that this div exists?
        }
      }

}

CCP4GraphPlot.prototype.handleEvent = function(evt) {

  if(evt.type=="change"){
    this.handleFileSelect(evt);
  }

  if(this.allowResizeEvent){
    //console.log("Can resize");
  }
  if(this.allowResizeEvent&&evt.type=="resize"&&this.currentGraph){
   $("#"+this.currentGraph.graphID).empty();
    this.plot();
  }
}

CCP4GraphPlot.prototype.handleFileSelect = function(evt) {
  var files = evt.target.files; // FileList object
  this.loadFiles(files);
}

CCP4GraphPlot.prototype.loadFile = function(file) {
  this.loadFiles([file]);
}

CCP4GraphPlot.prototype.loadFiles = function(files) {
  this.graphs = [];
    // files is a FileList of File objects. List some properties.
    var output = [];
    for (var i = 0, f; f = files[i]; i++) {
      output.push('<li><strong>', escape(f.name), '</strong> (', f.type || 'n/a', ') - ',
                  f.size, ' bytes, last modified: ',
                  f.lastModifiedDate ? f.lastModifiedDate.toLocaleDateString() : 'n/a',
                  '</li>');
      var reader = new FileReader();
      output.push(reader);
      reader.onabort = this.handleAbort;
      reader.onerror = this.handleError;
      reader.onloadstart = this.handleLoadStart;
      reader.onloadend = this.handleLoadEnd;
      reader.onload = this.handleLoad;
      reader.graphObject = this;
      //console.log("start reading");
      reader.readAsArrayBuffer(f);
    }
    //document.getElementById('list').innerHTML = '<ul>' + output.join('') + '</ul>';
    //console.log("handleFileSelect end");
}

// The "public methods".
CCP4GraphPlot.prototype.setCurrentData = function (iplot) {
  var graphName = this.graphs[iplot].graphID;
  $('#'+this.menuSelectName).val(graphName+'_selectitem').change();
}

CCP4GraphPlot.prototype.getDataWithName = function (graphName) {
  for(var ig=0;ig<this.graphs.length;ig++){
    if(this.graphs[ig].graphID == graphName){
       return this.graphs[ig];
    }
  }
  return null;
}

CCP4GraphPlot.prototype.loadFromID = function(elementID) {
  var source = document.getElementById(elementID).textContent;
  this.loadString(source);
}

CCP4GraphPlot.prototype.watchUrl = function(url,timeout) {
  //console.log("watchUrl !!!!!!!!");
  var self = this;
  window.setInterval(function(){self.loadFromUrl(url);}, timeout);
}

CCP4GraphPlot.prototype.loadFromUrl = function(url) {
  var txt = '';
  var xmlhttp = new XMLHttpRequest();
  var self = this;

  /*
  // This does not work!
  xmlhttp.onreadystatechanged = function(){
    if(xmlhttp.status == 200 && xmlhttp.readyState == 4){
      txt = xmlhttp.responseText;
      self.loadString(txt);
    }
  };
  xmlhttp.open("GET",url);
  xmlhttp.send(null);
  */

  // This works in Firefox, Qt WebKit and Safari with files in same directory as html page. But not Chrome.
  xmlhttp.open("GET",url,false);
  xmlhttp.overrideMimeType('text/plain;');
  xmlhttp.send();
  txt = xmlhttp.responseText;
  self.loadString(txt);
}

CCP4GraphPlot.prototype.loadString = function(txt) {
  var xmlDoc;
  if (window.DOMParser) {
    var parser=new DOMParser();
    xmlDoc=parser.parseFromString(txt,"text/xml");
  } else {
    xmlDoc=new ActiveXObject("Microsoft.XMLDOM");
    xmlDoc.async=false;
    xmlDoc.loadXML(txt);
  } 
  //console.log(xmlDoc);
  //console.log("read XML");
  var root = xmlDoc.getElementsByTagName("CCP4ApplicationOutput")[0];
  var tables = root.getElementsByTagName("CCP4Table")
  var fonts = root.getElementsByTagName("Fonts")

  var graphObject = this;
  graphObject.loadXML(graphObject,root,tables,fonts);
}

CCP4GraphPlot.prototype.replot = function() {
  $("#"+this.currentGraph.graphID).empty();
  this.plot();
}

CCP4GraphPlot.prototype.getCurrentPlotOptions = function() {
  return this.currentGraph.options;
}

CCP4GraphPlot.prototype.getCurrentPlotName = function() {
  return this.currentGraph.titles;
}

CCP4GraphPlot.prototype.getCurrentPlotData = function() {
  return this.currentGraph.graphData;
}

CCP4GraphPlot.prototype.getPlotOptions = function(iplot) {
  return this.graphs[iplot].options;
}

CCP4GraphPlot.prototype.getPlotName = function(iplot) {
  return this.graphs[iplot].titles;
}

CCP4GraphPlot.prototype.getPlotData = function(iplot) {
  return this.graphs[iplot].graphData;
}

CCP4GraphPlot.prototype.getNumberOfPlots = function() {
  return this.graphs.length;
}

return {
    CCP4GraphPlot : function(theDiv,showInput,allowResizeEvent){
        var ccp4i2Graph = new CCP4GraphPlot(theDiv,showInput,allowResizeEvent);
        return ccp4i2Graph
    }
}

});
