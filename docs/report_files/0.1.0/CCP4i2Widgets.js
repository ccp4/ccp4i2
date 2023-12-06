define(['jquery','jspimple','jstables'],function(jquery,jspimple,jstables){

var createSubElement = function(element, type, id){
    var newElement = document.createElement(type);
    element.appendChild(newElement)
    if (arguments.length == 3){
        newElement.setAttribute('id', id)
    }
    return newElement
}

function CCP4Plot(plotNode, dataObject){
    //alert('Here;)')
    if (!plotNode.getElementsByTagName){
        return null
    }
    this.dataObject = dataObject
    var titleNodes = plotNode.getElementsByTagName('title')
    if (titleNodes.length > 0){
        this.aTitle = new String(titleNodes[0].childNodes[0].nodeValue)
    }
    if (!this.aTitle) {
        this.aTitle = new String('Graph')
    }
    this.isBarChart = false;
    var plottypeNodes = plotNode.getElementsByTagName('plottype')
    if (plottypeNodes.length > 0 && plottypeNodes[0].childNodes[0].nodeValue.toLowerCase() == 'bar'){
        this.isBarChart = true;
    }
       var xscaleNodes = plotNode.getElementsByTagName('xscale')
       if (xscaleNodes.length > 0){
       this.scaleFunction = xscaleNodes[0].childNodes[0].nodeValue;
       }
    //Populate an array of plotlines, each in a ready-for-flot data structure
    this.plotLines = []
    var plotLineNodes = plotNode.getElementsByTagName('plotline')
    var lineColours=['red','green','blue','cyan','green','yellow','pink','purple'];
    for (var iPlotLineNode=0; iPlotLineNode < plotLineNodes.length; iPlotLineNode++){
        var plotLineNode = plotLineNodes[iPlotLineNode]
        var xcol = parseInt(plotLineNode.getAttribute('xcol'))-1;
        var ycol = parseInt(plotLineNode.getAttribute('ycol'))-1;
        
        var currentPlotLine = {
        label:this.dataObject.headers[ycol],
        color:lineColours[iPlotLineNode%(lineColours.length)],
        ignore:false,
        yaxis:1
        }
        if (plotLineNode.hasAttribute('colour')) currentPlotLine.color = plotLineNode.getAttribute('colour');
        currentPlotLine.ignore = false
        if (plotLineNode.hasAttribute('ignore') && plotLineNode.getAttribute('ignore').toLowerCase() == 'true') currentPlotLine.ignore = true;
        if (plotLineNode.hasAttribute('rightaxis')) currentPlotLine.yaxis = 2;
        currentPlotLine.data = []
        for (var iDataRow=0; iDataRow<this.dataObject.data.length; iDataRow++){
            dataRow = this.dataObject.data[iDataRow]
            var point = [dataRow[xcol], dataRow[ycol]];
            if (this.isBarChart) point = [iDataRow, dataRow[ycol]];
            currentPlotLine.data.push(point)
        }
        
        this.plotLines.push(currentPlotLine)
    }
    function makeClickFunc(divHashTag){
        var clickHashTag = divHashTag+"_click"
        var hoverHashTag = divHashTag+"_hover"
        function onClick(event, pos,item){
            if (item) {
                var x = item.datapoint[0].toFixed(2),
                y = item.datapoint[1].toFixed(2);
                var str = "(" + x + ", " + y + ")";
                $(clickHashTag).text(str);
            }
            
        }
        return onClick;
    }
    //And here the method that will actually draw a plot and all of its lines
    this.plotIn = function(divId){
        if (this.plotLines.length == 0){
            return
        }
        var rightAxisUsed = false
        
        //Select the lines that are not flagged as "ignored"
        var plotableData = []
        for (var iPlotLine=0; iPlotLine<this.plotLines.length; iPlotLine++){
            var plotLine = this.plotLines[iPlotLine];
            if (!plotLine.ignore){
                if (plotLine.yaxis ===2) rightAxisUsed = true;
                if (this.isBarChart){
                    plotLine.lines={show:false};
                    plotLine.points={show:false};
                    plotLine.bars={show:true, barWidth:0.4, order:iPlotLine};
                }
                plotableData.push(plotLine);
            }
        }
        //alert(JSON.stringify(plotableData))
        var controlDict = {
            series :{ points:{ show:true }, lines:{ show:true } },
            grid:{hoverable:true, clickable:true}
        };
        controlDict.yaxes = [{}, {position:'right'}]
       if (this.scaleFunction){
       if (this.scaleFunction == "oneoversqrt"){
       controlDict.xaxis  = {transform:function(x){return 1./Math.pow(x,2.0);}, inverseTransform:function(x){return 1./Math.pow(x,0.5);}}
       }
       }
        var divHashTag ="#"+divId;
        plotObject = $.plot(divHashTag, plotableData, controlDict);
        //$(this.divHashTag).bind("plotclick", makeClickFunc());
        var clickHashTag = divHashTag+"_click";
        var hoverHashTag = divHashTag+"_hover";
        $(divHashTag).bind("plotclick", function (event, pos, item) {
                           if (item) {
                           var x = item.datapoint[0].toFixed(2),
                           y = item.datapoint[1].toFixed(2);
                           var str = "(" + x + ", " + y + ")";
                           $(clickHashTag).text(str);
                           }
                           });
        
        $(divHashTag).bind("plothover", function (event, pos, item) {
                           if (item) {
                           var x = item.datapoint[0].toFixed(2),
                           y = item.datapoint[1].toFixed(2);
                           var str = "(" + x + ", " + y + ")";
                           $(hoverHashTag).text(str);
                           }
                           });
    }
    return this
}

function CCP4DataObject(xmlNode){
    //Grab Data Object title if there is one
    this.title = "CCP4 Table"
    var titleNodes = xmlNode.getElementsByTagName('title');
    if (xmlNode.hasAttribute('title')) this.title = xmlNode.getAttribute('title');
    
    //Identify the header node, throwing exception idf it is missing
    var headerNodes = xmlNode.getElementsByTagName('headers')
    if (headerNodes.length == 0) throw new Error("No headers");
    var headerNode = xmlNode.getElementsByTagName('headers')[0]
    
    //Pick out separator used in headers line
    separator = ' '
    if (headerNode.hasAttribute('separator')) separator = headerNode.getAttribute('separator');
    //Handle whitespace in separator
    var headerArray = headerNode.childNodes[0].nodeValue.trim().replace( /  +/g, ' ' ).split(separator)
    //And parse the headers line
    this.headers = headerArray.slice(0,headerArray.length)
    
    //Parse data block into an array of arrays
    var dataNodes = xmlNode.getElementsByTagName('data')
    for (var iDataNode=0; iDataNode < dataNodes.length; iDataNode++){
        dataNode = dataNodes[iDataNode]
        if (dataNode.childNodes){
            dataRows = dataNode.childNodes[0].nodeValue.trim().split('\n')
            dataRows = dataRows.slice(0,dataRows.length)
            this.data = new Array()
            for (var iRow=0; iRow<dataRows.length; iRow++){
                dataRowPossiblyFlabby = dataRows[iRow]
                dataRow = dataRowPossiblyFlabby.replace( /  +/g, ' ' )
                dataRowColumns = dataRow.trim().split(' ')
                dataRowColumns = dataRowColumns.slice(0,dataRowColumns.length)
                this.data.push(dataRowColumns)
            }
        }
    }
    
    //Identify the plots that are in this graph
    var plotNodes = xmlNode.getElementsByTagName('plot')
    this.plots = []
    for (var iPlotNode =0; iPlotNode< plotNodes.length; iPlotNode++){
        plotNode = plotNodes[iPlotNode]
        plot = new CCP4Plot(plotNode, this)
        if (plot != null){
            this.plots.push(plot)
        }
    }
    return this
}

       
function CCP4i2HTMLChunk(theDiv, widthStr, heightStr){
   this.div = theDiv;
   this.width = parseInt(widthStr);
   this.height = parseInt(heightStr);
   this.setData = function(text){
       this.text = text;
   };
   this.draw = function(){
       //Following to redraw in case of "text" being unchanged
       if (typeof this.div !== "undefined" && typeof this.div.innerHTML !== "undefined" && (!(this.div.innerHTML == this.text))){
            //console.log("text changed");
            this.div.innerHTML = this.text;
       }
       else {
            //console.log("Woohoo - i2HTMLChunk text unchanged");
       }
   }
   //Clear the div !
   while (this.div.firstChild){
       this.div.removeChild(this.div.firstChild);
   };
   return this;
}
       
CCP4i2HTMLChunk.prototype.loadFromUrl = function(url) {
    //alert('Loading from '+url);
    var txt = '';
    var xmlhttp = new XMLHttpRequest();
    var self = this;
    xmlhttp.open("GET",url,false);
    xmlhttp.overrideMimeType('text/plain;');
    xmlhttp.send();
    txt = xmlhttp.responseText;
    self.setData(txt);
    self.draw()
}
       
function CCP4i2Graph(theDiv, widthStr, heightStr){
    this.div = theDiv;
    this.width = parseInt(widthStr);
    this.height = parseInt(heightStr);
    //Clear the div !
    while (this.div.firstChild){
        this.div.removeChild(this.div.firstChild)
    }
    this.div.style.padding = "0px";
    //console.log(this.div.getAttribute('id'));
    this.headerDiv = createSubElement(this.div,'div',this.div.getAttribute('id')+'_header');
    this.headerDiv.style.cssText = 'width:'+this.width+'px; height:30px; margin:0px; border-style:None; background-color:#EEE;float:left;padding=0px; margin=0px;border-width:0px;'
    this.bodyDiv = createSubElement(this.div, 'div', this.div.getAttribute('id')+'_body');
    this.bodyDiv.style.cssText = 'width:'+this.width+'px; height:'+(this.height-60)+'px; margin:0px; border-style:None; background-color:#EEE;float:left;padding=0px; margin=0px;border-width:0px;';
    this.bodyDivId = this.bodyDiv.getAttribute('id');
    
    var feedbackDiv=createSubElement(this.div, 'div', this.div.getAttribute('id')+'_feedback');
    feedbackDiv.style.cssText = 'width:'+this.width+'px; height:30px; margin:0px; border-style:None; background-color:#EEE;float:left; padding=0px; margin=0px;border-width:0px;'
    var hoverDiv = createSubElement(feedbackDiv, 'div', this.bodyDivId+'_hover')
    hoverDiv.style.cssText = 'width:'+this.width/2+'px; height:30px; margin:0px; border-style:None; background-color:#EEE;float:left;padding=0px; margin=0px;border-width:0px;'
    var clickDiv = createSubElement(feedbackDiv, 'div', this.bodyDivId+'_click')
    clickDiv.style.cssText = 'width:'+this.width/2+'px; height:30px; margin:0px; border-style:None; background-color:#EEE;float:left;padding=0px; margin=0px;border-width:0px;'
    
    this.setData = function(text){
        var xmlDoc = $.parseXML( text )
        var plotDataElements = xmlDoc.getElementsByTagNameNS('http://www.ccp4.ac.uk/ccp4ns','ccp4_data')
        this.dataObject = new CCP4DataObject(plotDataElements[0]);
        //Clear the headerDiv
        while (this.headerDiv.firstChild) {
            this.headerDiv.removeChild(this.headerDiv.firstChild)
        }
        //Setup the "select" for choosing among available plots
        if (this.dataObject.plots.length == 0) return;
        else if (this.dataObject.plots.length == 1){
            this.headerDiv.innerHTML = this.dataObject.plots[0].aTitle
        }
        else {
            this.select = document.createElement('select');
            this.select.setAttribute('style','width:'+this.width-30+'px;')
            var thisGraph = this;
            this.select.onchange = function(){
                thisGraph.dataObject.plots[this.selectedIndex].plotIn(thisGraph.bodyDivId);
            }
            this.headerDiv.appendChild(this.select)
            for (var iPlot=0; iPlot< this.dataObject.plots.length; iPlot++){
                plot = this.dataObject.plots[iPlot];
                var option = createSubElement(this.select, 'option')
                option.setAttribute('value',iPlot)
                option.text = plot.aTitle
            }
        }
    }
    this.draw = function(){
        if (this.dataObject.plots.length == 0) {
        }
        else if (this.dataObject.plots.length == 1) {
            this.dataObject.plots[0].plotIn(this.bodyDivId);
        }
        else {
            this.dataObject.plots[this.select.selectedIndex].plotIn(this.bodyDivId);
        }
    }
    return this
}

function CCP4i2LineChooser(theDiv, tableWidthStr, contentWidthStr, heightStr){
    this.width = parseInt(tableWidthStr) + parseInt(contentWidthStr);
    this.height = parseInt(heightStr);
    this.tableWidth = parseInt(tableWidthStr);
    this.graphWidth = parseInt(contentWidthStr);
    //Clear the div !
    while (theDiv.firstChild){
        theDiv.removeChild(theDiv.firstChild)
    }
    theDiv.style.padding = "0px";
    //Add the line choosing Div
    this.tableDiv = createSubElement(theDiv, 'div', theDiv.getAttribute('id')+'_table');
    this.tableDiv.style.cssText = 'width:'+this.tableWidth+'px; height:'+this.height+'px; margin:0px; border-style:None; background-color:#EEE;float:left;padding=0px; margin=0px;border-width:0px;overflow:auto;'
    //Add the Graph Div
    this.graphDiv = createSubElement(theDiv, 'div', theDiv.getAttribute('id')+'_graph');
    this.graphDiv.style.cssText = 'width:'+this.graphWidth+'px; height:'+this.height+'px; margin:0px; border-style:None; background-color:#EEE;float:left;padding=0px; margin=0px;border-width:0px;'
    this.ccp4i2Graph = new CCP4i2Graph(this.graphDiv,this.graphWidth+'px',heightStr);
    var thisLineChooser = this;
    this.setData = function(dataTxt){
        this.ccp4i2Graph.setData(dataTxt);
        //Empty table
        while (this.tableDiv.firstChild){
            this.tableDiv.removeChild(this.tableDiv.firstChild)
        }
        var table = createSubElement(this.tableDiv, 'table')
        table.className = 'fancy'
        table.style.width = this.tableWidth+'px'
        var tbody = createSubElement(table, 'tbody')
        tbody.className = 'fancy'
        //Find plotlines in the first plot in the data
        var plotLines = this.ccp4i2Graph.dataObject.plots[0].plotLines
        for (var iPlotLine=0; iPlotLine < plotLines.length; iPlotLine++){
            plotLine = plotLines[iPlotLine];
            var row = createSubElement(tbody, 'tr')
            $(row).click( function(){
                         var correspondingPlot = thisLineChooser.ccp4i2Graph.dataObject.plots[0];
                         var correspondingGraph = thisLineChooser.ccp4i2Graph;
                         correspondingPlot.plotLines[this.value].ignore = !correspondingPlot.plotLines[this.value].ignore
                         $(this).children().each(function(){
                                                 if (!$(this).hasClass('galleryListObj')) $(this).addClass('galleryListObj');
                                                 if ($(this).hasClass('Selected')) $(this).removeClass('Selected');
                                                 else $(this).addClass('Selected');
                                                 });
                         correspondingPlot.plotIn(correspondingGraph.bodyDivId);
                         });
            row.style.width = '100%';
            row.style.lineHeight = '1.25'
            row.value = iPlotLine
            var td = createSubElement(row, 'td');
            if (!plotLine.ignore) td.className = 'galleryListObj Selected';
            td.innerHTML = plotLine.label;
            td.style.width='75%'
            var td = createSubElement(row, 'td');
            if (!plotLine.ignore) td.className = 'galleryListObj Selected';
            td.innerHTML = plotLine.data[plotLine.data.length-1][1];
            td.style.width='25%'
        }
    }
    this.draw = function(){
        this.ccp4i2Graph.draw()
    }
    return this
}

var CCP4i2WidgetsModule =
{
    CCP4i2HTMLChunk : function(data, widthStr, heightStr, theDiv){
       var ccp4i2HTMLChunk = new CCP4i2HTMLChunk(theDiv, widthStr, heightStr);
       if (theDiv.hasAttribute('data-is-urls') && theDiv.getAttribute('data-is-urls') == 'True'){
           ccp4i2HTMLChunk.loadFromUrl(data);
       }
       else {
           ccp4i2HTMLChunk.setData(data);
           ccp4i2HTMLChunk.draw();
       }
       return ccp4i2HTMLChunk;
    },
    CCP4i2Graph : function(data, widthStr, heightStr, theDiv){
        var ccp4i2Graph = new CCP4i2Graph(theDiv, widthStr, heightStr);
        ccp4i2Graph.setData(data);
        ccp4i2Graph.draw();
        return ccp4i2Graph
    },
    CCP4i2LineChooser: function(theDiv, tableWidthStr, contentWidthStr, heightStr){
        var ccp4i2LineChooser = new CCP4i2LineChooser(theDiv, tableWidthStr, contentWidthStr, heightStr);
        //A Fifth argument (if supplied) is the data
        if (arguments.length==5){
            if ('data' in arguments[4]){
                ccp4i2LineChooser.setData(arguments[4].data);
                ccp4i2LineChooser.draw();
            }
        }
       return ccp4i2LineChooser;
    },
    CCP4FlotRenderer: function(data, widthStr, heightStr, theDiv){
       //No need to instantiate a new CCP4GraphPlot if there is already one in this div
       //But on the other hand...there seems to be a memory leak which is fixed if I
       //Do simply create a new CCP4GraphPlot each time the graph is drawn
       var ccp4i2GraphPlot = $(theDiv).data('ccp4i2GraphPlot');
       //console.log("exisitng graph was "+typeof ccp4i2GraphPlot +(typeof ccp4i2GraphPlot  === 'undefined'))
       var currentIdx = 0;
       if (typeof ccp4i2GraphPlot  === 'undefined'){
           //console.log("No exisint graph")
           $(theDiv).empty()
           ccp4i2GraphPlot = jspimple.CCP4GraphPlot(theDiv.getAttribute('id'));
           $(theDiv).data('ccp4i2GraphPlot', ccp4i2GraphPlot);
           //console.log("creating a new plot");
       }
       else {
           currentIdx = ccp4i2GraphPlot.getCurrentData();
           //console.log("recycling plot");
       }
       var root = null;
       var tables = [];
       if (theDiv.hasAttribute('data-is-urls') && theDiv.getAttribute('data-is-urls') === 'True'){
           var tableUrls = data.split(",");
           for (var iTableUrl = 0; iTableUrl < tableUrls.length; iTableUrl++){
                var tableUrl = tableUrls[iTableUrl];
                ccp4i2GraphPlot.loadFromUrl(tableUrl);
           }
       }
       else {
           var tableNames = data.split(",");
           for (var iTableName = 0; iTableName < tableNames.length; iTableName++){
               var tableName = tableNames[iTableName];
               tables.push(document.getElementById(tableName));
           }
           ccp4i2GraphPlot.loadXML(ccp4i2GraphPlot,root,tables);
       }
       
       ccp4i2GraphPlot.replot();
       //AL, 5 Mar2018, definition of currentIdx was moved up
       //to prevent graph widget from switching to plot-0
       //var currentIdx=ccp4i2GraphPlot.getCurrentData();
       
       ccp4i2GraphPlot.setCurrentData(0);
       
       function getDescendantWithClass(element, clName) {
           var children = element.childNodes;
           for (var i = 0; i < children.length; i++) {
               if (children[i].className &&
                   children[i].className.split(' ').indexOf(clName) >= 0)
               {
                   return children[i];
                }
            }
            for (var i = 0; i < children.length; i++) {
                var match = getDescendantWithClass(children[i], clName);
                if (match !== null) {
                    return match;
                }
            }
            return null;
       }
       var selector = getDescendantWithClass(theDiv,"graphSelect");
       var appFontSize = theDiv.getAttribute('data-app-font-size');
       if(appFontSize!="undefined"&&appFontSize!=null&&selector!=null){
         selector.style.fontSize = appFontSize;
       }

       if(currentIdx){
         ccp4i2GraphPlot.setCurrentData(currentIdx);
         ccp4i2GraphPlot.replot();
       }

       return ccp4i2GraphPlot;
    },
    CCP4TableRenderer: function(data, widthStr, heightStr, theDiv){
       var ccp4i2Table = $(theDiv).data('ccp4i2Table');
       //console.log("exisitng table was "+typeof ccp4i2Table +(typeof ccp4i2Table  === 'undefined'))
       if (typeof ccp4i2Table  === 'undefined'){
           //console.log("No exisint table")
           $(theDiv).empty()
           ccp4i2Table = jstables.CCP4JSTable(theDiv.getAttribute('id'));
           $(theDiv).data('ccp4i2Table', ccp4i2Table);
           //console.log("creating new table");
       }
       else {
           //console.log("recycling table");
       }

       var root = null;
       var tables = [];
       if (theDiv.hasAttribute('data-is-urls') && theDiv.getAttribute('data-is-urls') == 'True'){
           var tableUrls = data.split(",");
           ccp4i2Table.loadFromUrl(tableUrls[0]);
       }
       else {
           var tableNames = data.split(",");
           ccp4i2Table.loadFromTableID(tableNames[0]);
       }
       return ccp4i2Table;
    }
}
       
return CCP4i2WidgetsModule;

});
