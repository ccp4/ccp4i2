
function StructureComparisonUI(ofGL,inDiv){
    var self = this;
    self.ofGL = ofGL;
    self.divId = "#"+inDiv;
    $(self.divId).data("structureComparisonUI",self);
    self.ofGL.div.addEventListener("statusChanged",
                                   function(e) {
                                   self.makeAndDrawTree();
                                   },
                                   false
                                   );
    self.makeAndDrawTree();
}

StructureComparisonUI.prototype.makeAndDrawTree = function()
{
    var self=this;
    self.displayBufferTreeData = self.makeTreeData();
    
    $(self.divId).html('Search for string :<input type="text" id="uiTree_q" value="" class="uiTree_q" style="margin:0em auto 1em auto; padding:4px; border-radius:4px; border:1px solid silver;" /><div class="uiTree"></div>')
    try{
        if (typeof self.uitree !== 'undefined') self.uiTree.jstree('destroy');
    }
    catch(err){
        console.log(err);
    }
    self.uiTree = $(self.divId + "> div.uiTree").jstree({  "types" : {
                                                        "MolData" : {"icon" : "MolDisp"},
                                                        "MapData" : {"icon" : "MapDisp"},
                                                        "Surface" : {"icon" : "Surface"},
                                                        "FATBONDS" : {"icon" : "FATBONDS"},
                                                        "CYLINDERS" : {"icon" : "CYLINDERS"},
                                                        "BONDS" : {"icon" : "BONDS"},
                                                        "SPLINE" : {"icon" : "SPLINE"},
                                                        "MapDisp" : {"icon" : "MapDisp"},
                                                        "BALLSTICK" : {"icon" : "BALLSTICK"}
                                                        },
                                                        "contextmenu": {"items": self.contextMenu, "select_node":false},
                                                        "core" : {"data": self.displayBufferTreeData},
                                                        "checkbox" : { "keep_selected_style" : false },
                                                        "three_state" : true,
                                                        "check_callback": true,
                                                        "plugins" : [ "checkbox","search","types","contextmenu","state" ]});
    self.uiTree.on('changed.jstree', function(e, data){
                   self.handleTreeSelect(e, data, gl);
                   gl.drawSceneDirty()});
    to = false;
    $(self.divId + "> input.uiTree_q").keyup(function () {
                                             if(to) { clearTimeout(to); }
                                             to = setTimeout(function () {
                                                             var v = $(self.divId + "> input.uiTree_q").val();
                                                             $(self.divId + "> div.uiTree").jstree(true).search(v,show_only_matches=true);
                                                             }, 250);
                                             });
    
}

StructureComparisonUI.prototype.makeTreeData = function(){
    var self=this;
    
    //Make a lookup
    self.uuidLookup = {};
    for (var idx = 0; idx < self.ofGL.displayBuffers.length; idx++){
        var displayBuffer = self.ofGL.displayBuffers[idx];
        var displayBufferUUID = self.ofGL.ids[idx].split('_')[1];
        self.uuidLookup[displayBufferUUID] = idx;
    }
    /*
     var displayBufferTreeData= [];
     */
    var displayBufferTreeData= [];
    var statusElement = self.ofGL.xmlDoc.getElementsByTagName("CCP4MG_Status")[0];
    
    var molDataElements = statusElement.getElementsByTagName("MolData");
    var molDataDirectoryLookup = {};
    
    for (var iMolData = 0; iMolData<molDataElements.length;iMolData++){
        var molData = molDataElements[iMolData];
        var nameElements = molData.getElementsByTagName("name");
        var nameElement = nameElements[0];
        for (var iNameElement=0; iNameElement<nameElements.length; iNameElement++){
            if (nameElements[iNameElement].parentNode === molData){
                nameElement = nameElements[iNameElement];
            }
        }
        var uuidElement = molData.getElementsByTagName("uuid")[0];
        var proteinUUID = uuidElement.childNodes[0].nodeValue;
        var filenameElements = molData.getElementsByTagName("filename");
        var fileDirectoryPath = '';
        var fileName = '';
        for (var iFilenameElement =0; iFilenameElement<filenameElements.length; iFilenameElement++){
            var filenameElement = filenameElements[iFilenameElement];
            var fullPathElements = filenameElement.getElementsByTagName("fullPath");
            for (var iFullPathElement =0; iFullPathElement<fullPathElements.length; iFullPathElement++){
                var fullPathElement = fullPathElements[iFullPathElement];
                var fullPath = fullPathElement.childNodes[0].nodeValue;
                var fullPathMatches = /^(A-Z:|\\|\/|\:).*(\\|\/|\:)/g.exec(fullPath);
                fileDirectoryPath = fullPathMatches[0];
                fileName=fullPath.replace(/^.*(\\|\/|\:)/, '');
            }
        }
        
        var molDataObject = new Object();
        molDataObject.id=proteinUUID;
        molDataObject.parent="#";
        molDataObject.type="MolData";
        molDataObject.text=nameElement.childNodes[0].nodeValue;
        molDataObject.divId = self.divId;
        molDataObject.fileDirectoryPath = fileDirectoryPath;
        molDataDirectoryLookup[fileDirectoryPath] = proteinUUID;
        molDataObject.fileName = fileName;
        displayBufferTreeData.push(molDataObject);
        var simpleMolDispElements = Array.prototype.slice.call(molData.getElementsByTagName("MolDisp"));
        var surfDispElements = Array.prototype.slice.call(molData.getElementsByTagName("SurfaceDispobj"));
        var molDispElements = simpleMolDispElements.concat(surfDispElements);
        
        for (var iMolDisp=0; iMolDisp<molDispElements.length; iMolDisp++){
            var molDispElement = molDispElements[iMolDisp];
            var nameElement = molDispElement.getElementsByTagName("name")[0];
            var uuidElement = molDispElement.getElementsByTagName("uuid")[0];
            var molDispUUID =uuidElement.childNodes[0].nodeValue;
            var molDispObject = new Object();
            molDispObject.id=molDispUUID;
            molDispObject.parent=proteinUUID;
            molDispObject.text="Object";
            molDispObject.divId = self.divId;
            molDispObject.type="demo";
            molDispObject.bufferIndex = iMolDisp;
            //Object generated by "snippet" uses <style></style> rather than <style_parameters><style_mode>
            var styleElements = molDispElement.getElementsByTagName("style_parameters");
            for (var iStyleElement=0; iStyleElement<styleElements.length; iStyleElement++){
                var styleElement = styleElements[iStyleElement];
                var modeElements = styleElement.getElementsByTagName("style_mode");
                for (var iModeElement=0; iModeElement<modeElements.length; iModeElement++){
                    var modeElement = modeElements[iModeElement];
                    molDispObject.type=modeElement.childNodes[0].nodeValue;
                }
            }
            var selectionParametersNodes = molDispElement.getElementsByTagName("selection_parameters");
            for (var iBlock=0; iBlock < selectionParametersNodes.length; iBlock ++){
                var cidNodes = selectionParametersNodes[iBlock].getElementsByTagName("cid");
                for (var iCIDNode = 0; iCIDNode<cidNodes.length; iCIDNode++){
                    try {
                        molDispObject.text = cidNodes[iCIDNode].childNodes[0].nodeValue;
                    }
                    catch(err){
                        console.log(err);
                    }
                }
            }
            var surfaceTypeElements = molDispElement.getElementsByTagName("surface_type");
            for (var iSurfaceTypeElement=0; iSurfaceTypeElement<surfaceTypeElements.length; iSurfaceTypeElement++){
                var surfaceTypeElement = surfaceTypeElements[iSurfaceTypeElement];
                molDispObject.text="Surface";
                molDispObject.type="Surface";
            }
            var doSelect = false;
            if (typeof(self.uuidLookup[molDispUUID]) !== "undefined"){
                var idx = self.uuidLookup[molDispUUID];
                var displayBuffer =self.ofGL.displayBuffers[idx];
                if (displayBuffer.visible == 1) doSelect = true;
            }
            molDispObject.state = {"selected":doSelect,"opened":false,"disabled":false};
            displayBufferTreeData.push(molDispObject);
        }
    }
    
    var mapDataElements = statusElement.getElementsByTagName("MapData");
    
    for (var iMapData = 0; iMapData<mapDataElements.length;iMapData++){
        var mapData = mapDataElements[iMapData];
        var nameElement = mapData.getElementsByTagName("name")[0];
        var uuidElement = mapData.getElementsByTagName("uuid")[0];
        var mapUUID = uuidElement.childNodes[0].nodeValue;
        var filenameElements = mapData.getElementsByTagName("filename");
        var fileDirectoryPath = '';
        var fileName = '';
        for (var iFilenameElement =0; iFilenameElement<filenameElements.length; iFilenameElement++){
            var filenameElement = filenameElements[iFilenameElement];
            var fullPathElements = filenameElement.getElementsByTagName("fullPath");
            for (var iFullPathElement =0; iFullPathElement<fullPathElements.length; iFullPathElement++){
                var fullPathElement = fullPathElements[iFullPathElement];
                var fullPath = fullPathElement.childNodes[0].nodeValue;
                var fullPathMatches = /^(A-Z:|\\|\/|\:).*(\\|\/|\:)/g.exec(fullPath);
                fileDirectoryPath = fullPathMatches[0];
                fileName=fullPath.replace(/^.*(\\|\/|\:)/, '');
            }
        }
        
        var mapDataObject = new Object();
        mapDataObject.id=mapUUID;
        mapDataObject.parent="#";
        if (typeof(molDataDirectoryLookup[fileDirectoryPath]) !== 'undefined'){
            mapDataObject.parent = molDataDirectoryLookup[fileDirectoryPath];
        }
        
        mapDataObject.text=nameElement.childNodes[0].nodeValue;
        mapDataObject.divId=self.divId;
        mapDataObject.fileDirectoryPath = fileDirectoryPath;
        mapDataObject.fileName = fileName;
        displayBufferTreeData.push(mapDataObject);
        var mapDispElements = Array.prototype.slice.call(mapData.getElementsByTagName("MapDisp"));
        
        for (var iMapDisp=0; iMapDisp<mapDispElements.length; iMapDisp++){
            var mapDispElement = mapDispElements[iMapDisp];
            var nameElement = mapDispElement.getElementsByTagName("name")[0];
            var uuidElement = mapDispElement.getElementsByTagName("uuid")[0];
            var mapDispUUID =uuidElement.childNodes[0].nodeValue;
            var mapDispObject = new Object();
            mapDispObject.id=mapDispUUID;
            mapDispObject.parent=mapUUID;
            var mapType = "Map";
            var differenceElements = mapData.getElementsByTagName("difference");
            for (var iElement=0; iElement < differenceElements.length; iElement++){
                if (differenceElements[iElement].childNodes[0].nodeValue === '1'){
                    mapType = "Difference Map"
                }
            }
            var mapLevel = "1.0";
            var elements = mapDispElement.getElementsByTagName("contour_level");
            for (var iElement=0; iElement < elements.length; iElement++){
                mapLevel = elements[iElement].childNodes[0].nodeValue;
            }
            mapDispObject.text=mapType +" @ "+mapLevel;
            mapDispObject.divId=self.divId;
            mapDispObject.type="MapDisp";
            var doSelect = false;
            if (typeof(self.uuidLookup[mapDispUUID]) !== "undefined"){
                var idx = self.uuidLookup[mapDispUUID];
                var displayBuffer =self.ofGL.displayBuffers[idx];
                if (displayBuffer.visible == 1) doSelect = true;
            }
            mapDispObject.state = {"selected":doSelect,"opened":false,"disabled":false};
            displayBufferTreeData.push(mapDispObject);
        }
    }
    
    return displayBufferTreeData;
}


StructureComparisonUI.prototype.handleTreeSelect = function(e,data,gl){
    var self=this;
    var ofGL = $(self.divId).data('structureComparisonUI').ofGL;
    
    for (var idx = 0; idx < ofGL.displayBuffers.length; idx++){
        var displayBuffer = ofGL.displayBuffers[idx];
        var displayBufferUUID = ofGL.ids[idx].split('_')[1];
        if (data.selected.indexOf(displayBufferUUID)>-1){
            ofGL.displayBuffers[idx].visible = true;
        }
        else {
            ofGL.displayBuffers[idx].visible = false;
        }
    }
    ofGL.drawSceneDirty();
}

StructureComparisonUI.prototype.mapsOff = function(){
    var self=this;
    for (var iNode = 0; iNode < self.displayBufferTreeData.length; iNode++){
        var aNode = self.displayBufferTreeData[iNode];
        if (aNode.type === 'MapDisp'){
            $(self.divId + "> div.uiTree").jstree("deselect_node", aNode.id);
        }
    }
    gl.drawSceneDirty();
}

StructureComparisonUI.prototype.surfsOff = function(){
    var self=this;
    for (var iNode = 0; iNode < self.displayBufferTreeData.length; iNode++){
        var aNode = self.displayBufferTreeData[iNode];
        if (aNode.type === 'Surface'){
            $(self.divId + "> div.uiTree").jstree("deselect_node", aNode.id);
        }
    }
    gl.drawSceneDirty();
}

StructureComparisonUI.prototype.contextMenu = function(node) {
    var self=this;
    
    var theNode = node;
    var nodeId = theNode.original.id;
    // The default set of all items
    var items = {};
    if (theNode.original.type === 'MolData'){
        items = {
        deleteItem: {
        label: "Add display object",
        action: function (data) {
            var correspondingTree = data.reference.jstree();
            var structureComparisonUI = correspondingTree.element.parent().data("structureComparisonUI");
            structureComparisonUI.ofGL.addDisplayObject(nodeId);
        }
        },
        };
    }
    else if (theNode.original.type !== 'MapData'){
        var items = {};
        var displayOption  = $(theNode.original.divId).data("structureComparisonUI").ofGL.displayOptions[theNode.original.bufferIndex];
        for (var key in displayOption){
            items[key] = {label:key,submenu:{}};
            var menu = items[key];
            for (var iStyle=0; iStyle<displayOption[key].length; iStyle++){
                var style = displayOption[key][iStyle];
                //Handle first simple case where style is a pair of values: one is the label and one is the corresponding identifier
                if (style.length == 2){
                    var newItem = {label:style[1], key:key, style:style[0], icon:false};
                    if (typeof style[1] !=='string') {
                        newItem = {label:style[0]+style[1][1], key:key, style:style[1][0], icon:false};
                    }
                    newItem.action = function(data){
                        var correspondingTree = data.reference.jstree();
                        var structureComparisonUI = correspondingTree.element.parent().data("structureComparisonUI");
                        if (data.item.key === 'Colour' && data.item.style === 'browser'){
                            $("#colourBrowser > button").data('gl',structureComparisonUI.ofGL);
                            $("#colourBrowser > button").data('nodeNumber',structureComparisonUI.uuidLookup[nodeId]);
                            $("#colourBrowser > button").data('key',data.item.key);
                            $("#colourBrowser > button").data('subkey',data.item.style);
                            $("#colourBrowser > button").click(function(){
                                                               $(this).data('gl').updateDisplayObject($(this).data('nodeNumber'),
                                                                                                      $(this).data('key'),
                                                                                                      "one_colour_"+$("#colourBrowser > input").val());
                                                               $(this).unbind( "click" );
                                                               $("#colourBrowser").dialog("close");
                                                               });
                            $("#colourBrowser").dialog("open");
                        }
                        else {
                            structureComparisonUI.ofGL.updateDisplayObject(structureComparisonUI.uuidLookup[nodeId],data.item.key,data.item.style);
                        }
                    };
                    menu.submenu[newItem.label] = newItem;
                }
                else {
                    var newItem = { label:style[0], key:key, icon:false, submenu:{} };
                    menu.submenu[newItem.label] = newItem;
                    var submenu = menu.submenu[newItem.label].submenu;
                    for (var iSubitem=1; iSubitem<style.length; iSubitem++){
                        submenu[style[iSubitem][1]] = { label:style[iSubitem][1], key:key, subkey:style[iSubitem][0], icon:false,
                        action:function(data){
                            var correspondingTree = data.reference.jstree();
                            var structureComparisonUI = correspondingTree.element.parent().data("structureComparisonUI");
                            structureComparisonUI.ofGL.updateDisplayObject(structureComparisonUI.uuidLookup[nodeId],data.item.key,data.item.subkey);
                        }
                        };
                    }
                }
            }
        }
    }
    var structureComparisonUI = $(theNode.original.divId).data("structureComparisonUI");
    if (typeof structureComparisonUI.uuidLookup[theNode.original.id] !== 'undefined'){
        items["centreon"] = {"label":"Centre on", "listIndex":structureComparisonUI.uuidLookup[theNode.original.id], action:function(data){
            var correspondingTree = data.reference.jstree();
            var structureComparisonUI = correspondingTree.element.parent().data("structureComparisonUI");
            structureComparisonUI.ofGL.centreOn(data.item.listIndex);
        }
        };
        items["delete"] = {"label":"Delete", "listIndex":structureComparisonUI.uuidLookup[theNode.original.id], action:function(data){
            var correspondingTree = data.reference.jstree();
            var structureComparisonUI = correspondingTree.element.parent().data("structureComparisonUI");
            structureComparisonUI.ofGL.removeDisplayObject(data.item.listIndex);
        }
        };
    }
    return items;
}

function FogClipBackgroundUI (inDivId, gl){
    var self=this;
    self.fogClipDiv = document.createElement("div");
    
    var divName = gl.div.getAttribute("id");
    
    var fogText = document.createElement("div");
    fogText.innerHTML = "Fog";
    self.fogClipDiv.appendChild(fogText);
    var fogSlider = document.createElement("div");
    fogSlider.setAttribute("id",divName+"fog-slider-range");
    self.fogClipDiv.appendChild(fogSlider);
    
    var clipText = document.createElement("div");
    clipText.innerHTML = "Clip";
    self.fogClipDiv.appendChild(clipText);
    var clipSlider = document.createElement("div");
    clipSlider.setAttribute("id",divName+"clip-slider-range");
    self.fogClipDiv.appendChild(clipSlider);
    
    $(function() {
      $( "#"+divName+"clip-slider-range" ).slider({
                                                  range: true,
                                                  min: -250,
                                                  max: 250,
                                                  values: [ -250, 250 ],
                                                  slide: function( event, ui ) {
                                                  gl.set_clip_range(ui.values[ 0 ],ui.values[ 1 ],true);
                                                  }
                                                  });
      gl.set_clip_range($( "#"+divName+"clip-slider-range" ).slider( "values", 0 ),$( "#"+divName+"clip-slider-range" ).slider( "values", 1 ),false);
      });
    
    $(function() {
      $( "#"+divName+"fog-slider-range" ).slider({
                                                 range: true,
                                                 min: 0,
                                                 max: 1000,
                                                 values: [ 440, 640 ],
                                                 slide: function( event, ui ) {
                                                 gl.set_fog_range(ui.values[ 0 ],ui.values[ 1 ],true);
                                                 }
                                                 });
      gl.set_fog_range($( "#"+divName+"fog-slider-range" ).slider( "values", 0 ),$( "#"+divName+"fog-slider-range" ).slider( "values", 1 ),false);
      });
    
    $("#"+inDivId).append(self.fogClipDiv);
}

function ColourPickerUI (inDivId, ofGL){
    var self=this;
    self.divId = inDivId;
    self.ofGL = ofGL;
    self.colourPickDiv = $('<div id='+inDivId+' class="colourPicker"><div class="colourBlockDiv" style="width:100px;"></div><input class="myInput" type="text" /><button class="formSaver">Accept!</button></div>').get(0);
    document.body.appendChild(self.colourPickDiv);
    var colourElements = ofGL.xmlDoc.getElementsByTagName("colour_definitions")[0].getElementsByTagName("colour");
    for (var iColour=0; iColour<colourElements.length; iColour++){
        try {
            var colourElement = colourElements[iColour];
            var redValue = parseFloat(colourElement.getElementsByTagName("red")[0].childNodes[0].nodeValue);
            var greenValue = parseFloat(colourElement.getElementsByTagName("green")[0].childNodes[0].nodeValue);
            var blueValue = parseFloat(colourElement.getElementsByTagName("blue")[0].childNodes[0].nodeValue);
            var colorDiv = $("<div class='colourCell' style='width:10px;height:10px; float:left; margin:1px; background-color:rgb("+(redValue*255).toFixed(0)+","+(greenValue*255).toFixed(0)+","+(blueValue*255).toFixed(0)+")'></div>");
            $(colorDiv).data('red',(redValue*255));
            $(colorDiv).data('green',(greenValue*255));
            $(colorDiv).data('blue',(blueValue*255));
            $(colorDiv).data('name',colourElement.getAttribute('name'));
            $(colorDiv).click(function(){
                              var redString = parseInt($(this).data('red').toFixed(0)).toString(16);
                              if (redString.length<2) redString = "0"+redString;
                              var greenString = parseInt($(this).data('green').toFixed(0)).toString(16);
                              if (greenString.length<2) greenString = "0"+greenString;
                              var blueString = parseInt($(this).data('blue').toFixed(0)).toString(16);
                              if (blueString.length<2) blueString = "0"+blueString;
                              var cbDivId = $(this).parent().parent().get(0).getAttribute('id');
                              $("#" + cbDivId + " > input").val(redString+greenString+blueString);
                              });
            $("#"+inDivId+" > div.colourBlockDiv").append(colorDiv);
        }
        catch (err) {}
    }
    $("#"+inDivId).dialog({autoOpen: false, show: {effect: "blind",duration: 100},
                          hide: {effect: "explode",duration: 100}});
}
