function SnapShot(gl, structureComparisonUI){
    var self=this;
    self.replace(gl, structureComparisonUI);
    self.text = 'New snapshot';
}

SnapShot.prototype.replace = function(gl, structureComparisonUI)
{
    var self=this;
    
    var displayBufferTreeData = structureComparisonUI.displayBufferTreeData;
    self.zoom = gl.zoom;
    self.gl_fog_start = gl.gl_fog_start;
    self.gl_fog_end = gl.gl_fog_end;
    self.gl_clipPlanes = new Array();
    var gl_clipPlanes = [gl.gl_clipPlane0, gl.gl_clipPlane1, gl.gl_clipPlane2, gl.gl_clipPlane3, gl.gl_clipPlane4, gl.gl_clipPlane5, gl.gl_clipPlane6, gl.gl_clipPlane7];
    for (var iClip =0; iClip < gl_clipPlanes.length; iClip++){
        var newClipPlane = new glMatrixArrayType(4);
        for (var i=0; i<4; i++){
            newClipPlane[i] = gl_clipPlanes[iClip][i];
        }
        self.gl_clipPlanes.push(newClipPlane);
    }
    self.origin = vec3.create(gl.origin);
    self.myQuat = quat4.create(gl.myQuat);
    self.background_colour = new Array();
    for (var i=0; i<4; i++){
        self.background_colour[i] = gl.background_colour[i];
    }
    self.id = guid();
    self.parent = '#';
    self.visibleObjects = new Array();
    for (var iId=0; iId < displayBufferTreeData.length; iId++){
        var id = displayBufferTreeData[iId].id;
        if (typeof(structureComparisonUI.uuidLookup[id]) !== 'undefined'){
            var idx = structureComparisonUI.uuidLookup[id];
            if (gl.displayBuffers[idx].visible){
                self.visibleObjects.push(id);
            }
        }
    }
    
}

SnapShot.prototype.restoreView = function(gl, structureComparisonUI)
{
    var self=this;
    
    var displayBufferTreeData = structureComparisonUI.displayBufferTreeData;
    gl.zoom = self.zoom;
    gl.gl_fog_start = self.gl_fog_start;
    gl.gl_fog_end = self.gl_fog_end;
    gl.origin = vec3.create(self.origin);
    gl.myQuat = quat4.create(self.myQuat);
    if (typeof self.background_colour !== 'undefined') {
        for (var i=0; i<4; i++){
            gl.background_colour[i] = self.background_colour[i];
        }
    }
    //Check for clipPlanes in snapshot
    if (typeof(self.gl_clipPlanes) !== 'undefined'){
        gl_clipPlanes = [gl.gl_clipPlane0, gl.gl_clipPlane1, gl.gl_clipPlane2, gl.gl_clipPlane3, gl.gl_clipPlane4, gl.gl_clipPlane5, gl.gl_clipPlane6, gl.gl_clipPlane7];
        for (var iClip=0; iClip < self.gl_clipPlanes.length; iClip++){
            for (var i=0; i<4; i++){
                gl_clipPlanes[iClip][i] = self.gl_clipPlanes[iClip][i];
            }
        }
    }
    gl.drawSceneDirty();
}

SnapShot.prototype.restore = function(gl, structureComparisonUI)
{
    var self=this;
    
    SnapShot.prototype.restoreView.call(self, gl, structureComparisonUI);
    SnapShot.prototype.restoreObjects.call(self, gl, structureComparisonUI);
    
    gl.drawSceneDirty();
}

SnapShot.prototype.restoreObjects = function(gl, structureComparisonUI)
{
    var self=this;
    
    var displayBufferTreeData = structureComparisonUI.displayBufferTreeData;
    
    for (var iNode = 0; iNode < displayBufferTreeData.length; iNode++){
        var aNode = displayBufferTreeData[iNode];
        structureComparisonUI.uiTree.jstree("deselect_node", aNode.id);
    }
    for (var iId = 0; iId < self.visibleObjects.length; iId++){
        var theId = self.visibleObjects[iId];
        structureComparisonUI.uiTree.jstree("select_node", theId);
    }
    gl.drawSceneDirty();
}

function SnapshotUI(divId, structureComparisonUI, gl)
{
    var self=this;
    
    self.divId = divId;
    self.structureComparisonUI = structureComparisonUI;
    self.gl = gl;
    if (typeof(localStorage.snapShotsJSON) === 'undefined'){
        self.allSnapShots = new Object();
        self.allSnapShots[document.URL] = new Array();
    }
    else {
        self.allSnapShots = JSON.parse(localStorage.snapShotsJSON);
        if (typeof self.allSnapShots[document.URL] === 'undefined'){
            self.allSnapShots[document.URL] = new Array();
            self.saveSnapshots();
        }
    }
    //Attach this interface object to the containing div
    $(self.divId).data("UIObject", self);
    $(self.divId).html('<button class="createSnapshot">snap</button><button class="clearSnapshots">clear</button><button class="showSnapshots">Show</button><div class="treeDivClass" style="border:1px solid red;"></div>')
    $(self.divId + ">button.createSnapshot").click(function(){self.createSnapshot()});
    $(self.divId + ">button.clearSnapshots").click(function(){self.clearSnapshots()});
    $("#SnapshotDialog").dialog({
                            autoOpen: false,
                            hide: "puff",
                            show : "slide",
                            height: 400
                            });
    $(self.divId + ">button.showSnapshots").click(function() {
                           $( "#dialog-3" ).dialog( "open" );
                           });
    $(self.divId + ">button.showSnapshots").click(function(){
                                                  $("#SnapshotDialog").text("var InitialSnapshots = "+JSON.stringify(self.allSnapShots[document.URL])+";\n");
                                                  $( "#SnapshotDialog" ).dialog( "open" );
                                                  });
    self.createTree();
}

SnapshotUI.prototype.saveSnapshots = function(snapshotsToSave)
{
    var self=this;
    if (typeof snapshotsToSave === 'undefined') snapshotsToSave = self.allSnapShots;
    var snapshotString = JSON.stringify(snapshotsToSave);
    localStorage.snapShotsJSON = snapshotString;
}

SnapshotUI.prototype.createTree = function(){
    var self = this;
    
    self.jsTree = $(self.divId +"> div.treeDivClass").jstree({ "core" : {"data": self.allSnapShots[document.URL],"check_callback":true},
                                                             "three_state" : false, "contextmenu":{"items": self.contextMenu, "select_node":false},
                                                             'plugins': ["contextmenu","search"],
                                                             });
    self.jsTree.on('changed.jstree', function(e, data){
                   self.handleSnapSelect(e, data);
                   self.gl.drawScene()
                   });
}

SnapshotUI.prototype.contextMenu = function(node) {
    var cmSelf=this;
    var theNode = node;
    // The default set of all items
    var items = {
        "delete": {
        label: "Delete", "theNode":theNode,
        action: function (data) {
            var theNode = data.item.theNode;
            var correspondingTree = data.reference.jstree();
            var uiObject = correspondingTree.element.parent().data("UIObject");
            var nodeId = theNode.original.id;
            for (var iSnap=0; iSnap<uiObject.allSnapShots[document.URL].length; iSnap++){
                if (uiObject.allSnapShots[document.URL][iSnap].id === nodeId){
                    if (confirm('Delete node '+iSnap+' labelled "'+theNode.original.text+'"')){
                        uiObject.allSnapShots[document.URL].splice(iSnap,1);
                        uiObject.saveSnapshots();
                        correspondingTree.delete_node(theNode);
                    }
                    break;
                }
            }
        }
        },
        "rename": {label: "Rename", "theNode":theNode,
        action: function(data){
            var theNode = data.item.theNode;
            var correspondingTree = data.reference.jstree();
            correspondingTree.edit(theNode, theNode.text,
                                   function(editedNode, exitStatus, didCancel){
                                   var nodeId = editedNode.original.id;
                                   var uiObject = correspondingTree.element.parent().data("UIObject");
                                   for (var iSnap=0; iSnap<uiObject.allSnapShots[document.URL].length; iSnap++){
                                   if (uiObject.allSnapShots[document.URL][iSnap].id === nodeId){
                                   uiObject.allSnapShots[document.URL][iSnap].text = editedNode.text;
                                   }
                                   }
                                   uiObject.saveSnapshots();
                                   });
            return true;
            
        }
        },
        "restoreView": {label:"Restore view", "theNode":theNode,
        action:function(data){
            var nodeId = data.item.theNode.id;
            var correspondingTree = data.reference.jstree();
            var uiObject = correspondingTree.element.parent().data("UIObject");
            for (var iSnap=0; iSnap<uiObject.allSnapShots[document.URL].length; iSnap++){
                if (uiObject.allSnapShots[document.URL][iSnap].id === nodeId){
                    SnapShot.prototype.restoreView.call(uiObject.allSnapShots[document.URL][iSnap], uiObject.gl, uiObject.structureComparisonUI);
                }
            }
        }
        },
        "restoreObjects": {label:"Restore objects", "theNode":theNode,
        action:function(data){
            var nodeId = data.item.theNode.id;
            var correspondingTree = data.reference.jstree();
            var uiObject = correspondingTree.element.parent().data("UIObject");
            for (var iSnap=0; iSnap<uiObject.allSnapShots[document.URL].length; iSnap++){
                if (uiObject.allSnapShots[document.URL][iSnap].id === nodeId){
                    SnapShot.prototype.restoreObjects.call(uiObject.allSnapShots[document.URL][iSnap], uiObject.gl, uiObject.structureComparisonUI);
                }
            }
        }
        },
        "replace": {label:"Replace", "theNode":theNode,
        action:function(data){
            var nodeId = data.item.theNode.id;
            var correspondingTree = data.reference.jstree();
            var uiObject = correspondingTree.element.parent().data("UIObject");
            for (var iSnap=0; iSnap<uiObject.allSnapShots[document.URL].length; iSnap++){
                if (uiObject.allSnapShots[document.URL][iSnap].id === nodeId){
                    SnapShot.prototype.replace.call(uiObject.allSnapShots[document.URL][iSnap], uiObject.gl, uiObject.structureComparisonUI);
                    uiObject.saveSnapshots();
                }
            }
        }
        }
    }
    return items;
}

SnapshotUI.prototype.clearSnapshots = function(){
    var self=this;
    self.allSnapShots = new Object();
    self.allSnapShots[document.URL] = new Array();
    self.saveSnapshots();
    try{
        self.jsTree.jstree('destroy');
    }
    catch(err){
        console.log(err);
    }
    self.createTree();
}

SnapshotUI.prototype.handleSnapSelect = function(e, data){
    var self = this;
    for (var iSnap=0; iSnap<self.allSnapShots[document.URL].length; iSnap++){
        var snapshot = self.allSnapShots[document.URL][iSnap];
        if (snapshot.id == data.selected[0]){
            SnapShot.prototype.restore.call(snapshot, self.gl, self.structureComparisonUI);
        }
    }
}

SnapshotUI.prototype.createSnapshot = function()
{
    var self=this;
    var newSnap = new SnapShot(self.gl, self.structureComparisonUI);
    self.allSnapShots[document.URL].push(newSnap);
    self.saveSnapshots();
    self.jsTree.jstree("create_node","#",newSnap,"last");
}

