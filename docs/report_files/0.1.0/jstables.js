require.config({
              shim: {
                   'jstables': {
                       deps: ['jquery'],
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
               'jquery': 'jquery.min',
               'jstables': 'jstables'
               }
               });

define(['jquery'],function(apple){


function getNodeText(node) {
    var r = "";
    for (var x = 0;x < node.childNodes.length; x++) {
        r = r + node.childNodes[x].nodeValue;
    }
    return r;
}

function getNodeTextWithMarkup(node) {
    var r = "";
    for (var x = 0;x < node.childNodes.length; x++) {
        if(node.childNodes[x].nodeType==3){
          r = r + node.childNodes[x].nodeValue;
        } else if(node.childNodes[x].nodeType==1) {
       //console.log(node.childNodes[x].nodeName);
          r = r + "<"+node.childNodes[x].nodeName+">" + getNodeTextWithMarkup(node.childNodes[x]) + "</"+node.childNodes[x].nodeName+">";
        }
    }
    return r;
}

function JSTable(divName_in){
  //console.log("JSTable");
  this.tableDivName = divName_in;
  this.selected_rows = [];
  this.selected_columns = [];
  this.selected_cells = [];
  this.transpose = false;
}

JSTable.prototype.redraw = function() {
  if(!this.transpose){
    var theDiv = document.getElementById(this.tableDivName);
    var table = theDiv.childNodes[0];
    // Why does this not work?
    //var body = table.getElementsByTagName("tbody");
    var body = table.childNodes[1];
    var children = body.childNodes;
    for(irow = 0; irow < children.length; irow++){
      var childCells = children[irow].childNodes;
      for(icell = 0; icell < childCells.length; icell++){
        var idx = irow * childCells.length + icell;
        var index = this.selected_cells.indexOf(idx);
        if(this.selected_rows.indexOf(irow) > -1 || this.selected_columns.indexOf(icell)>-1 || index>-1){
          childCells[icell].setAttribute("class","highlight-cell");
       } else {
          childCells[icell].setAttribute("class","cell");
       }
      }
    }
  } else {
    //console.log("No transpose highlighting yet!");
  }
}

JSTable.prototype.highlightRow = function(row) {
  this.selected_rows.push(row);
  this.redraw();
}

JSTable.prototype.unHighlightRow = function(row) {
  var index = this.selected_rows.indexOf(row);
  if (index > -1) {
    this.selected_rows.splice(index, 1);
  }
  this.redraw();
}

JSTable.prototype.highlightColumn = function(col) {
  this.selected_columns.push(col);
  this.redraw();
}

JSTable.prototype.unHighlightColumn = function(col) {
  var index = this.selected_columns.indexOf(col);
  if (index > -1) {
    this.selected_columns.splice(index, 1);
  }
  this.redraw();
}

JSTable.prototype.highlightCell = function(row,col) {
  var theDiv = document.getElementById(this.tableDivName);
  var table = theDiv.childNodes[0];
  // Why does this not work?
  //var body = table.getElementsByTagName("tbody");
  var body = table.childNodes[1];
  var children = body.childNodes;
  if(children.length>0){
    var ncols = children[0].childNodes.length;
    if(col<ncols){
      var idx = row * ncols + col;
      this.selected_cells.push(idx);
    }
  }
  this.redraw();
}

JSTable.prototype.unHighlightCell = function(row,col) {
  var theDiv = document.getElementById(this.tableDivName);
  var table = theDiv.childNodes[0];
  // Why does this not work?
  //var body = table.getElementsByTagName("tbody");
  var body = table.childNodes[1];
  var children = body.childNodes;
  if(children.length>0){
    var ncols = children[0].childNodes.length;
    if(col<ncols){
      var idx = row * ncols + col;
      var index = this.selected_cells.indexOf(idx);
      if (index > -1) {
        this.selected_cells.splice(index, 1);
      }
    }
  }
  this.redraw();
}

JSTable.prototype.loadTableXML = function(tableObject,table) {
  //console.log("loadTableXML");
  var theTable = table.getElementsByTagName("table");

  tableplot = document.getElementById(tableObject.tableDivName);
  var itable = 0;
  while (tableplot.firstChild) {
     tableplot.removeChild(tableplot.firstChild);
  }

  var thisPlotDiv;
  tableObject.tables = [];

  var N = 20;

  var table;

  //console.log("JSTable.prototype.loadTableXML 3");

  var tableDiv = document.getElementById(this.tableDivName);
  tableDiv.setAttribute("class","table");

  
  var realTable = document.createElement("table")
  if(theTable.length>0&&theTable[0].classList.contains("center")){
      realTable.classList.add("center");
  }

  var tableHead = document.createElement("thead");
  var tableHeadR = document.createElement("tr");

  var origTHead = table.getElementsByTagName("thead");
  if(origTHead==null|| origTHead.length==0){
    //console.log("Probably transpose");
    this.transpose = true;
    var tableBody = document.createElement("tbody");
    tableBody.setAttribute("class","fancy");
    var origTBody = table.getElementsByTagName("tbody");
    var origTRBody = origTBody[0].getElementsByTagName("tr");
    for(var idatarow=0;idatarow<origTRBody.length;idatarow++){
    var tableRow = document.createElement("tr");
      tableRow.setAttribute("class","row");
      tableRow.setAttribute("data-row",""+idatarow);
      var origTHBody = origTRBody[idatarow].getElementsByTagName("th");
      var origTDBody = origTRBody[idatarow].getElementsByTagName("td");
      var tableCell = document.createElement("th");
      tableCell.setAttribute("class","hcell");
       if (origTHBody[0].hasAttribute("title")){
       tableCell.setAttribute("title",origTHBody[0].getAttribute("title"));
       }
      tableCell.innerHTML = getNodeTextWithMarkup(origTHBody[0]);
      tableRow.appendChild(tableCell)
      for(var jdatacol=0;jdatacol<origTDBody.length;jdatacol++){
        var tableCell = document.createElement("td");
       if (origTDBody[jdatacol].hasAttribute("title")){
       tableCell.setAttribute("title",origTDBody[jdatacol].getAttribute("title"));
       }
        tableCell.setAttribute("class","cell");
        tableCell.setAttribute("data-row",""+idatarow);
        tableCell.setAttribute("data-column",""+jdatacol);
        tableCell.innerHTML = getNodeTextWithMarkup(origTDBody[jdatacol]);
        tableRow.appendChild(tableCell)
        tableCell.i = idatarow;
        tableCell.j = jdatacol;
        tableCell.onclick = function(){
          //console.log("Click cell "+this.i+" "+this.j)
          tableCell.evt = document.createEvent("CustomEvent");
          tableCell.evt.initCustomEvent('cellClicked', false, false, { 'row': this.i, 'column': this.j });
          tableplot.dispatchEvent(tableCell.evt);
        }
      }
      tableBody.appendChild(tableRow)
    }
    realTable.appendChild(tableBody)
    tableDiv.appendChild(realTable)
    return;
  }
  var origTRHead = origTHead[0].getElementsByTagName("tr");
  var origTHHead = origTRHead[0].getElementsByTagName("th");

  for(var idatacol=0;idatacol<origTHHead.length;idatacol++){
    var tableCell = document.createElement("th");
    tableCell.setAttribute("class","hcell");
    var colspan = origTHHead[idatacol].getAttribute("colspan");
    if(colspan != null){
      tableCell.setAttribute("colspan",colspan);
    }
       if (origTHHead[idatacol].hasAttribute("title")){
       tableCell.setAttribute("title",origTHHead[idatacol].getAttribute("title"));
       }
    tableCell.innerHTML = getNodeTextWithMarkup(origTHHead[idatacol]);
    tableHeadR.appendChild(tableCell)
  }
  tableHead.appendChild(tableHeadR)
  realTable.appendChild(tableHead)

  var origTBody = table.getElementsByTagName("tbody");
  var tableBody = document.createElement("tbody");
  tableBody.setAttribute("class","fancy");

  var origTRBody = origTBody[0].getElementsByTagName("tr");
  for(var idatarow=0;idatarow<origTRBody.length;idatarow++){
    var tableRow = document.createElement("tr");
    tableRow.setAttribute("class","row");
    tableRow.setAttribute("data-row",""+idatarow);
    var origTDBody = origTRBody[idatarow].getElementsByTagName("td");
    for(var jdatacol=0;jdatacol<origTDBody.length;jdatacol++){
      var tableCell = document.createElement("td");
      tableCell.setAttribute("class","cell");
      tableCell.setAttribute("data-row",""+idatarow);
      tableCell.setAttribute("data-column",""+jdatacol);
       if (origTDBody[jdatacol].hasAttribute("title")){
       tableCell.setAttribute("title",origTDBody[jdatacol].getAttribute("title"));
       }
      tableCell.innerHTML = getNodeTextWithMarkup(origTDBody[jdatacol]);
      tableRow.appendChild(tableCell)
      tableCell.i = idatarow;
      tableCell.j = jdatacol;
      tableCell.onclick = function(){
        //console.log("Click cell "+this.i+" "+this.j)
        tableCell.evt = document.createEvent("CustomEvent");
        tableCell.evt.initCustomEvent('cellClicked', false, false, { 'row': this.i, 'column': this.j });
        tableplot.dispatchEvent(tableCell.evt);
      }
    }
    tableBody.appendChild(tableRow)
  }

  realTable.appendChild(tableBody)

  /*
  var staticHeadDiv = document.createElement("div");
  staticHeadDiv.innerHTML = "Foo";
  tableDiv.appendChild(staticHeadDiv);
  var scrollingTableDiv = document.createElement("div");
  scrollingTableDiv.appendChild(realTable);
  scrollingTableDiv.setAttribute("style","height:250px;overflow:auto;");
  tableDiv.appendChild(scrollingTableDiv)
  */

  tableDiv.appendChild(realTable)

  //console.log("JSTable.prototype.loadTableXML end");

}

JSTable.prototype.loadXML = function(tableObject,root,tables) {
  //console.log("JSTable.prototype.loadXML");
  if(!tableObject.tableDivName){
    return;
  }
  //console.log("JSTable.prototype.loadXML 2");

  tableplot = document.getElementById(tableObject.tableDivName);
  var itable = 0;
  while (tableplot.firstChild) {
     tableplot.removeChild(tableplot.firstChild);
  }
  var thisPlotDiv;
  tableObject.tables = [];

  var N = 20;

  var table;

  //console.log("JSTable.prototype.loadXML 3");

  for(var itab=0;itab<tables.length;itab++){
    //console.log(tables[itab])
    var title = tables[itab].getAttribute("title");
    //console.log(title)
    var data = $.trim(getNodeText(tables[itab].getElementsByTagName("data")[0])).split("\n");
    //console.log(data)
    var headers_sep = tables[itab].getElementsByTagName("headers")[0].getAttribute("separator");
    var headers = $.trim(tables[itab].getElementsByTagName("headers")[0].childNodes[0].nodeValue).replace(/ +/g,' ').split(headers_sep);
    var dataColRow = [];
    for(var idatacol=0;idatacol<data.length;idatacol++){
      dataColRow.push([]);
    }
    for(var idatarow=0;idatarow<data.length;idatarow++){
      var thisDataRow = $.trim(data[idatarow]).replace(/ +/g,' ').split(" ");
      for(var idatacol=0;idatacol<headers.length;idatacol++){
        if(thisDataRow[idatacol]=="-"){
          //console.log("A null!")
          dataColRow[idatarow].push(null);
        }else{
          //console.log(thisDataRow[idatacol]);
          dataColRow[idatarow].push(parseFloat(thisDataRow[idatacol]));
        }
      }
    }
  }

  var tableDiv = document.getElementById(this.tableDivName);
  tableDiv.setAttribute("class","table");

  var realTable = document.createElement("table")

  var tableHead = document.createElement("thead");
  var tableHeadR = document.createElement("tr");

  for(var idatacol=0;idatacol<headers.length;idatacol++){
    var tableCell = document.createElement("th");
    tableCell.setAttribute("class","hcell");
    tableCell.innerHTML = headers[idatacol];
    tableHeadR.appendChild(tableCell)
  }
  tableHead.appendChild(tableHeadR)
  realTable.appendChild(tableHead)

  var tableBody = document.createElement("tbody");

  for(var idatarow=0;idatarow<dataColRow.length;idatarow++){
    var tableRow = document.createElement("tr");
    tableRow.setAttribute("class","row");
    tableRow.setAttribute("data-row",""+idatarow);
    for(var jdatacol=0;jdatacol<dataColRow[idatarow].length;jdatacol++){
      var tableCell = document.createElement("td");
      tableCell.setAttribute("class","cell");
      tableCell.setAttribute("data-row",""+idatarow);
      tableCell.setAttribute("data-column",""+jdatacol);
      tableCell.innerHTML = ""+dataColRow[idatarow][jdatacol];
      //console.log(dataColRow[idatacol])
      tableRow.appendChild(tableCell)
      tableCell.i = idatarow;
      tableCell.j = jdatacol;
      tableCell.onclick = function(){
        //console.log("Click cell "+this.i+" "+this.j)
        tableCell.evt = document.createEvent("CustomEvent");
        tableCell.evt.initCustomEvent('cellClicked', false, false, { 'row': this.i, 'column': this.j });
        tableplot.dispatchEvent(tableCell.evt);
      }
    }
    tableBody.appendChild(tableRow)
  }


  realTable.appendChild(tableBody)
  tableDiv.appendChild(realTable)
  //console.log("JSTable.prototype.loadXML end");

}

JSTable.prototype.loadString = function(txt){
  //console.log("JSTable.prototype.loadString");
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
  var tableObject = this;
  tableObject.loadXML(tableObject,root,tables);
}

JSTable.prototype.loadTableString = function(txt){
    var myself=this;
    var textHasChanged = false;
    if ("undefined" === typeof myself.currentText){
       //console.log("No current text in table");
       textHasChanged = true;
    }
    else {
       if (txt !== myself.currentText) {
          //console.log("Changed text in table");
          textHasChanged = true;
       }
       else {
          //console.log("Woo hoo - Unchanged text in table");
       }
    }
    if (textHasChanged){
           myself.currentText = txt;
          //console.log("JSTable.prototype.loadString");
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
          var tableObject = this;
          tableObject.loadTableXML(tableObject,xmlDoc);
    }
}

JSTable.prototype.loadFromID = function(elementID){
  //console.log("JSTable.prototype.loadFromID");
  var source = document.getElementById(elementID).textContent;
  this.loadString(source);
}

JSTable.prototype.loadFromTableID = function(elementID){
  //console.log("JSTable.prototype.loadFromTableID");
  var source = document.getElementById(elementID).textContent;
  this.loadTableString(source);
}
       
JSTable.prototype.loadFromUrl = function(url) {
   //alert('Loading from '+url);
       //console.log('Url is'+url);
   var txt = '';
   var xmlhttp = new XMLHttpRequest();
   var self = this;
   /*
    // This works in Firefox with files in same directory as html page, but not Qt WebKit or Safari.
    xmlhttp.onreadystatechanged = function(){
    if(xmlhttp.status == 200 && xmlhttp.readyState == 4){
    txt = xmlhttp.responseText;
    self.loadString(txt);
    }
    };
    */
   // This works in Firefox, Qt WebKit and Safari with files in same directory as html page.
   xmlhttp.open("GET",url,false);
   xmlhttp.overrideMimeType('text/plain;');
   xmlhttp.send();
   txt = xmlhttp.responseText;
   self.loadTableString(txt);
}
       
return {
    CCP4JSTable : function(theDiv){
        var ccp4i2Graph = new JSTable(theDiv);
        return ccp4i2Graph
    }
}



});
