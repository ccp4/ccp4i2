function genSymMats(RO,RF,symmats,origin,radius,centre){
    var TMats = [];
    var isymops = [];
    // FIXME We should work out shift limits on the fly also. Currently this only works in box of +/- shift.
    var transOrigin = vec3.create();
    var diffO = vec3.create();

    var xyz = vec3.create([1.0,1.0,1.0]);
    mat4.multiplyVec3(RO,xyz);

    var ofrac = vec3.create(origin);
    mat4.multiplyVec3(RF,ofrac);

    var x_shifts = Math.ceil(2 * radius/ xyz[0]);
    var y_shifts = Math.ceil(2 * radius/ xyz[1]);
    var z_shifts = Math.ceil(2 * radius/ xyz[2]);

    var x_shift_base = Math.floor(ofrac[0]);
    var y_shift_base = Math.floor(ofrac[1]);
    var z_shift_base = Math.floor(ofrac[2]);

    for(var i=0;i<symmats.length;i++){
        var tm = symmats[i];
        //console.log(tm);
        for(var xshift=-x_shifts+x_shift_base;xshift<x_shifts+x_shift_base;xshift++){
            for(var yshift=-y_shifts+y_shift_base;yshift<y_shifts+y_shift_base;yshift++){
                for(var zshift=-z_shifts+z_shift_base;zshift<z_shifts+z_shift_base;zshift++){
                    var tmt = mat4.create();
                    mat4.transpose(tm,tmt);
                    var fm = mat4.create();
                    var theTMatrix = mat4.create();
                    var invTMat = mat4.create();
                    tmt[12] += xshift;
                    tmt[13] += yshift;
                    tmt[14] += zshift;
                    //console.log("tmt");
                    //printMat(tmt);
                    mat4.multiply(tmt,RF,fm);
                    mat4.multiply(RO,fm,theTMatrix);
                    mat4.inverse(theTMatrix,invTMat);
                    //console.log("theTMatrix");
                    //printMat(theTMatrix);
                    mat4.multiplyVec3(theTMatrix,centre,transOrigin);
                    vec3.subtract(transOrigin,origin,diffO);
                    var diffO2 = vec3.length(diffO);
                    if(diffO2<radius){
                        TMats.push(theTMatrix);
                        isymops.push(i);
                    }
                }
            }
        }
    }
    return {"matrices":TMats,"symopnums":isymops};
}

function handleTextureLoaded(gl, image, texture, text, tex_size, font) {

    // FIXME - For text, need to fathom how to create a quad of appropriate size to draw on, and how to create correct sized canvas for that.
    var acanvas = document.createElement("canvas");

    var nptw;
    var npth;
    if(image){
        nptw = next_power_of_2(image.width);
        npth = next_power_of_2(image.height);
    } else {
        nptw = 512;
        npth = 256;
    }

    var ctx = acanvas.getContext("2d");

    var scalew = 1.0;
    var scaleh = 1.0;
    if(nptw>1024){
        scalew = 1024.0/nptw;
        nptw = 1024;
    }
    if(npth>1024){
        scaleh = 1024.0/npth;
        npth = 1024;
    }
    acanvas.width  = nptw;
    acanvas.height = npth;

    if(image){
        ctx.scale(scalew,scaleh);
        ctx.drawImage(image,0,0);
    } else {
        var fnsize = font.match(/^\d+|\d+\b|\d+(?=\w)/g)[0];
        ctx.font = font;
        var textWidth = ctx.measureText(text).width;
        var textHeight = 1.0*determineFontHeight(font,fnsize);
        nptw = next_power_of_2(parseInt(textWidth));
        npth = next_power_of_2(parseInt(textHeight));
        console.log(nptw+" "+npth);
        console.log(ctx.measureText(text));

        /*
           var kerned_pos = [0];
           for(var ic=1;ic<text.length;ic++){
           var diff2 = ctx.measureText(text.substring(ic-1,ic+1)).width - (ctx.measureText(text.substring(ic-1,ic)).width + ctx.measureText(text.substring(ic,ic+1)).width);
           kerned_pos.push(diff2+kerned_pos[ic-1]+ctx.measureText(text.substring(ic-1,ic)).width);
           }
         */
        acanvas.width  = nptw;
        acanvas.height = npth;
        ctx.font = font;
        // This is the way we can split things up for different fonts, etc. Should we want to. 
        /*
           for(var ic=0;ic<text.length;ic++){
           var rfrac = parseInt(255*1.0*ic/text.length);
           var gfrac = parseInt(255-255.0*ic/text.length);
           var rstr = rfrac.toString(16);
           var gstr = gfrac.toString(16);
           if(rstr.length==1) rstr = "0"+rstr;
           if(gstr.length==1) gstr = "0"+gstr;
           ctx.fillStyle="#"+rstr+gstr+"00";
           ctx.fillText(text.substring(ic,ic+1),kerned_pos[ic], acanvas.height / 2);
           }
         */
        //ctx.fillStyle="green";
        //ctx.fillRect(acanvas.width/2-10, acanvas.height/2-10,20, 20);
        ctx.fillStyle="black";
        ctx.fillText(text,acanvas.width/2-textWidth/2, textHeight);
    }

    gl.bindTexture(gl.TEXTURE_2D, texture);

    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, acanvas);

    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.generateMipmap(gl.TEXTURE_2D);
    gl.bindTexture(gl.TEXTURE_2D, null);
    tex_size["width"] = acanvas.width;
    tex_size["height"] = acanvas.height;
}

function initStringTextures(gl,text,tex_size,font) {
    var cubeTexture = gl.createTexture();
    handleTextureLoaded(gl, null, cubeTexture,text,tex_size,font);
    return cubeTexture;
}

function initTextures(gl,fname) {
    var cubeTexture = gl.createTexture();
    var cubeImage = new Image();
    var tex_size = {};

    cubeImage.src = fname;

    cubeImage.onload = function() { handleTextureLoaded(gl, cubeImage, cubeTexture, null, tex_size, null); }

    return cubeTexture;
}

function getNodeText(node) {
    var r = "";
    for (var x = 0;x < node.childNodes.length; x++) {
        r = r + node.childNodes[x].nodeValue;
    }
    return r;
}

function guid(){
    var d = Date.now();
    var uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
        var r = (d + Math.random()*16)%16 | 0;
        d = Math.floor(d/16);
        return (c=='x' ? r : (r&0x3|0x8)).toString(16);
    });
    return uuid;
}

function getOffsetRect(elem) {
    var box = elem.getBoundingClientRect();
    var body = document.body;
    var docElem = document.documentElement;

    var scrollTop = window.pageYOffset || docElem.scrollTop || body.scrollTop;
    var scrollLeft = window.pageXOffset || docElem.scrollLeft || body.scrollLeft;
    var clientTop = docElem.clientTop || body.clientTop || 0;
    var clientLeft = docElem.clientLeft || body.clientLeft || 0;
    var top  = box.top +  scrollTop - clientTop;
    var left = box.left + scrollLeft - clientLeft;
    return { top: Math.round(top), left: Math.round(left) };
}

function ajaxRequest(){
    var activexmodes=["Msxml2.XMLHTTP", "Microsoft.XMLHTTP"] //activeX versions to check for in IE
        if (window.ActiveXObject){ //Test for support for ActiveXObject in IE first (as XMLHttpRequest in IE7 is broken)
            for (var i=0; i<activexmodes.length; i++){
                try{
                    return new ActiveXObject(activexmodes[i]);
                }
                catch(e){
                    //suppress error
                }
            }
        }
        else if (window.XMLHttpRequest) // if Mozilla, Safari etc
            return new XMLHttpRequest();
        else
            return false;
}

function next_power_of_2(v){
    v -= 1 ;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v += 1;
    return v;
}

function appendBuffer( buffer1, buffer2 ) {
    var tmp = new Uint8Array( buffer1.byteLength + buffer2.byteLength );
    tmp.set( new Uint8Array( buffer1 ), 0 );
    tmp.set( new Uint8Array( buffer2 ), buffer1.byteLength );
    return tmp;
}

function getEncodedData( rssentries, createFun ) {
    var allBuffers = [];
    for (var i=0; i<rssentries.length; i++){
        if(typeof(rssentries[i])==="string"){
            var b64Data = rssentries[i];
            var strData;
            if(window.atob&&window.btoa){
                strData = atob(b64Data.replace(/\s/g, ''));
            } else {
                strData = base64decode(b64Data.replace(/\s/g, ''));
            }
            //var binData     = new Uint8Array();
            var j;

            var binData     = new Uint8Array(strData.length);
            for(j=0;j<strData.length;j++){
                binData[j] = strData[j].charCodeAt(0);
            }
            /*
            // HMM, WHY DO I DO THIS COMPLICATION, when above seems **much** faster
            for(j=0;j<strData.length/60000;j++){
            var lower = j*60000;
            var upper = (j+1)*60000;
            if(upper>strData.length){
            upper = strData.length;
            }
            // Turn number array into byte-array
            //var charData    = strData.substring(lower,upper).split('').map(function(x){return x.charCodeAt(0);});
            //var binDataThis     = new Uint8Array(charData);

            var binDataThis     = new Uint8Array(upper-lower);
            var ic = 0;
            for(var ibuf=lower;ibuf<upper;ibuf++){
            binDataThis[ic++] = strData[ibuf].charCodeAt(0);
            }

            binData = appendBuffer(binData,binDataThis);
            }
             */

            var data  = pako.inflate(binData);

            var strData = "";

            if(window.TextDecoder){
                // THIS'LL only work in Firefox 19+, Opera 25+ and Chrome 38+.
                var decoder = new TextDecoder('utf-8');
                var strData = decoder.decode(data);
            } else {
                var unpackBufferLength = 60000;
                var theSlice = new Uint8Array(60000);
                for(j=0;j<data.length/unpackBufferLength;j++){
                    var lower = j*unpackBufferLength;
                    var upper = (j+1)*unpackBufferLength;
                    if(upper>data.length){
                        upper = data.length;
                    }
                    // FECK, no slice on Safari!
                    strData += String.fromCharCode.apply(null, data.subarray(lower,upper));
                }   
            }


            try {
                var thisBuffer=JSON.parse(strData);
            } catch(e) {
                console.log(strData);
            }

            if(typeof(createFun)==='undefined'){
                allBuffers.push(thisBuffer);
            } else {
                createFun(thisBuffer);
            }
        } else {
            allBuffers.push(rssentries[i]);
        }
    }
    return allBuffers;
}

function initGL(canvas) {
    try {
        var gl = canvas.getContext("experimental-webgl");
        gl.viewportWidth = canvas.width;
        gl.viewportHeight = canvas.height;
    } catch (e) {
    }
    if (!gl) {
        //alert("Could not initialise WebGL, sorry :-(");
    }
    console.log("Max varying vectors: "+gl.getParameter(gl.MAX_VARYING_VECTORS));
    return gl;
}

function getShader(gl, str, type) {

    var shader;
    if (type == "fragment") {
        shader = gl.createShader(gl.FRAGMENT_SHADER);
    } else if (type == "vertex") {
        shader = gl.createShader(gl.VERTEX_SHADER);
    } else {
        return null;
    }

    gl.shaderSource(shader, str);
    gl.compileShader(shader);

    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        alert(gl.getShaderInfoLog(shader));
        return null;
    }

    return shader;
}

function sortIndicesByProj(a,b) {
    if (a.proj > b.proj)
        return -1;
    if (a.proj < b.proj)
        return 1;
    return 0;
}

function SortThing(proj,id1,id2,id3) {
    this.proj = proj;
    this.id1 = id1;
    this.id2 = id2;
    this.id3 = id3;
}

function DisplayBuffer() {
    this.visible = true;
    this.name_label = "";
    this.display_class = "NONE";
    this.transparent = false;
    this.atoms = [];
    this.texture = null;
    this.clearBuffers();
}

DisplayBuffer.prototype.clearBuffers = function () {
    this.triangleVertexNormalBuffer = [];
    this.triangleVertexPositionBuffer = [];
    this.triangleVertexIndexBuffer = [];
    this.triangleVertexTextureBuffer = [];
    this.triangleColourBuffer = [];
    this.triangleIndexs = [];
    this.triangleVertices = [];
    this.triangleColours = [];
    this.triangleNormals = [];
    this.primitiveSizes = [];
    this.bufferTypes = [];
    this.symmetry = null;
    this.symopnums = [];
    this.supplementary = {};
    this.isDirty = true;
}

DisplayBuffer.prototype.setSymmetryMatrices = function (symmetry) {
    this.symmetry = symmetry;
}

function uploadpdb(pdbid,self) {
    console.log(pdbid.value);

    var client = new XMLHttpRequest();

    var theFile = guid()+pdbid.value+".ent";

    /* Create a FormData instance */
    var formData = new FormData();
    if (typeof(middlewareToken) !== 'undefined'){
        formData.append('csrfmiddlewaretoken',middlewareToken);
    }
    /* Add the file */
    //console.log("Adding to form 2 "+file.files[0]);
    formData.append("pdbcode",pdbid.value);
    formData.append("pdbfilename",theFile);

    /* Check the response status */  
    client.onreadystatechange = function() 
    {
        if (client.readyState == 4 && client.status == 200) 
        {
            self.mygetrequest.open("GET", "get_triangles?pdb="+theFile, true);
            self.mygetrequest.send(null);
        }
    }

    client.open("POST", "/server_client.html", true);
    client.send(formData);
    /* Send to server */
}

function upload(file,self) 
{
    if(file.files[0].name.length>4&&file.files[0].name.substr(file.files[0].name.length-4)==="json"){
        console.log("json!");
        var r = new FileReader();
        var f = file.files[0];
        r.onload = function(e) { 
            var contents = e.target.result;
            var jsondata=JSON.parse(contents);
            self.loadJSON(jsondata);
            /*
            alert( "Got the file.n" 
                    +"name: " + f.name + "n"
                    +"type: " + f.type + "n"
                    +"size: " + f.size + " bytesn"
                    + "starts with: " + contents.substr(1, contents.indexOf("n"))
                 );  
                 */
        }
        r.readAsText(f);
        return;
    }
    var client = new XMLHttpRequest();

    var theFile = guid()+file.files[0].name;

    /* Create a FormData instance */
    var formData = new FormData();
    if (typeof(middlewareToken) !== 'undefined'){
        formData.append('csrfmiddlewaretoken',middlewareToken);
    }
    /* Add the file */
    //console.log("Adding to form 2 "+file.files[0]);
    formData.append("upfile", file.files[0],theFile);

    /* Check the response status */  
    client.onreadystatechange = function() 
    {
        if (client.readyState == 4 && client.status == 200) 
        {
            self.mygetrequest.open("GET", "get_triangles?pdb="+theFile, true);
            self.mygetrequest.send(null);
        }
    }

    client.open("POST", "/server_client.html", true);
    client.send(formData);  /* Send to server */ 
}

function DihedralAngle(A,B,C,D){
    var AB = vec3.create();
    var BC = vec3.create();
    var CD = vec3.create();

    vec3.subtract(B,A,AB);
    vec3.subtract(C,B,BC);
    vec3.subtract(D,C,CD);

    var Q = vec3.create();
    var T = vec3.create();
    var S = vec3.create();

    vec3.cross(AB,BC,Q);
    vec3.cross(BC,CD,T);
    vec3.cross(Q,T,S);

    var ss = vec3.dot(S,S);
    var qt = vec3.dot(Q,T);

    var angle = Math.atan2(Math.sqrt(ss),qt);
    var sbc = S[0]*BC[0] + S[1]*BC[1] + S[2]*BC[2];

    if(sbc<0.0) {
        return -angle;
    } else {
        return angle;
    }

}

function DistanceBetweenTwoLines(p1,p2,p3,p4){
    var ret = [];

    var a1 = vec3.create();
    var a2 = vec3.create();
    vec3.subtract(p2,p1,a1);
    vec3.subtract(p4,p3,a2);

    ret[0] = -1.0;
    ret[1] = -1.0;
    ret[2] = -1.0;

    if(vec3.length(a1)===0){
        return ret;
    }

    var n = vec3.create();
    vec3.cross(a1,a2,n);

    if(Math.abs(vec3.length(n))<1e-6){
        return DistanceBetweenPointAndLine(p1,p2,p3);
    }

    if(vec3.length(a2)===0){
        return ret;
    }

    vec3.normalize(n);

    var p4p1 = vec3.create();
    vec3.subtract(p4,p1,p4p1);

    var dist = Math.abs(vec3.dot(p4p1,n));
    var a1sq = vec3.length(a1)*vec3.length(a1);
    var a2sq = vec3.length(a2)*vec3.length(a2);

    var u = (a1sq*(vec3.dot(a2,p3) - vec3.dot(a2,p1)) + vec3.dot(a1,a2)*(vec3.dot(a1,p1) - vec3.dot(a1,p3)))/ (vec3.dot(a1,a2)*vec3.dot(a1,a2) - a1sq*a2sq);

    var t = (vec3.dot(a1,p3) + u * vec3.dot(a1,a2) - vec3.dot(a1,p1)) / a1sq;

    ret[0] = dist;
    ret[1] = t;
    ret[2] = u;

    return ret;

}

function DistanceBetweenPointAndLine(ls, le, p){

    var ret = [];
    ret.push(-1.0);

    var lemls = vec3.create(le); 
    vec3.subtract(lemls,ls);

    var linesize = vec3.length(lemls);
    //console.log("ls: "+ls[0]+" "+ls[1]+" "+ls[2]);
    //console.log("le: "+le[0]+" "+le[1]+" "+le[2]);
    //console.log("p: "+p[0]+" "+p[1]+" "+p[2]);

    if(Math.abs(linesize)<1e-6){
        console.log("Zero length line in DistanceBetweenPointAndLine");
        return ret;
    }

    /* t is value in line equation p = p1 + t(p2-p1) */

    var pmls = vec3.create(p); 
    vec3.subtract(pmls,ls);

    var t = vec3.dot(pmls,lemls) / (linesize*linesize);

    var ttlemls = vec3.create(lemls);
    ttlemls[0] *= t;
    ttlemls[1] *= t;
    ttlemls[2] *= t;

    var pt = vec3.create();
    vec3.add(ls,ttlemls,pt);

    vec3.subtract(pt,p);

    ret[0] = vec3.length(pt);
    ret[1] = t;

    return ret;

}

// FIXME - Definite funniness here!

function createXQuatFromDX(angle_in){
    var angle = angle_in*Math.PI/180.0;
    var q = quat4.create([Math.cos(angle/2.0),0.0,0.0,Math.sin(angle/2.0)]);
    return q;
}

function createYQuatFromDY(angle_in){
    var angle = angle_in*Math.PI/180.0;
    var q = quat4.create([Math.cos(angle/2.0),0.0,Math.sin(angle/2.0),0.0]);
    return q;
}

function createZQuatFromDX(angle_in){
    var angle = angle_in*Math.PI/180.0;
    var q = quat4.create([0.0,0.0,Math.sin(angle/2.0),Math.cos(angle/2.0)]);
    return q;
}

function getDeviceScale() {
    var deviceScale = 1.0;
    if (typeof(window.devicePixelRatio) !== 'undefined') {
        deviceScale = Math.min(1.5,Math.max(1.0,0.75*window.devicePixelRatio));
    }
    return deviceScale;
}

function updateColourPicker(picker,theUuid){
    console.log(theUuid);
}

function MGWebGLLightUI(gl,div) {
    this.infoNode = document.getElementById(div);
    this.canvas = document.createElement("canvas");
    this.infoNode.appendChild(this.canvas);
    this.gl = gl;
    this.canvas.width = parseInt(this.infoNode.style.width)/2;
    this.canvas.height = parseInt(this.infoNode.style.height)/2;
    this.ctx = this.canvas.getContext("2d");
    this.ctx.fillStyle="black";
    this.ctx.fillRect(0,0,this.canvas.width,this.canvas.height);
    var self = this;
    self.canvas.addEventListener("mousedown", function(evt){ self.doMouseDown(evt,self); }, false);
    self.canvas.addEventListener("mousemove", function(evt){ self.doMouseMove(evt,self); }, false);
    self.canvas.addEventListener("mouseup", function(evt){ self.doMouseUp(evt,self); }, false);
    self.imageData = this.ctx.getImageData(0,0,this.canvas.width, this.canvas.height);
    self.drawGradient(this.canvas.width/2, this.canvas.height/2);
    self.canvas.style.cssFloat="left";
    self.uuid = guid();

    var ambientValues = ptorgb(this.gl.light_colours_ambient);
    var specularValues = ptorgb(this.gl.light_colours_specular);
    var diffuseValues = ptorgb(this.gl.light_colours_diffuse);

    function toHex(d){
        return ("0"+(parseInt(d).toString(16))).slice(-2).toUpperCase();
    }

    var pickerOptions = {mode:"HS"};
    var cleardiv = document.createElement("div");
    cleardiv.style.clear="both";
    this.infoNode.appendChild(cleardiv);
    var labeldiv = document.createElement("div");
    labeldiv.innerHTML = "Ambient&nbsp;";
    labeldiv.style.cssFloat="left";
    this.infoNode.appendChild(labeldiv);
    this.ambientSlider = document.createElement("input");
    this.ambientSlider.setAttribute("type","range");
    this.infoNode.appendChild(this.ambientSlider);
    var ambR = toHex(ambientValues["r"]*255);
    var ambG = toHex(ambientValues["g"]*255);
    var ambB = toHex(ambientValues["b"]*255);
    var ambientColour = ambR+ambG+ambB;
    var ambientColourPickerInput = document.createElement("input");
    ambientColourPickerInput.setAttribute("value",ambientColour);
    var ambientUUID = guid();
    ambientColourPickerInput.setAttribute("id",ambientUUID);
    ambientColourPickerInput.style.width = "30px";
    ambientColourPickerInput.style.height = "10px";
    this.infoNode.appendChild(ambientColourPickerInput);
    var ambientColourPicker = new jscolor(ambientUUID,pickerOptions);

    var cleardiv = document.createElement("div");
    cleardiv.style.clear="both";
    this.infoNode.appendChild(cleardiv);
    var labeldiv = document.createElement("div");
    labeldiv.innerHTML = "Specular";
    labeldiv.style.cssFloat="left";
    this.infoNode.appendChild(labeldiv);
    this.specularSlider = document.createElement("input");
    this.specularSlider.setAttribute("type","range");
    this.infoNode.appendChild(this.specularSlider);
    var specR = toHex(specularValues["r"]*255);
    var specG = toHex(specularValues["g"]*255);
    var specB = toHex(specularValues["b"]*255);
    var specularColour = specR+specG+specB;
    var specularColourPickerInput = document.createElement("input");
    specularColourPickerInput.setAttribute("value",specularColour);
    var specularUUID = guid();
    specularColourPickerInput.setAttribute("id",specularUUID);
    specularColourPickerInput.style.width = "30px";
    specularColourPickerInput.style.height = "10px";
    this.infoNode.appendChild(specularColourPickerInput);
    var specularColourPicker = new jscolor(specularUUID,pickerOptions);

    var cleardiv = document.createElement("div");
    cleardiv.style.clear="both";
    this.infoNode.appendChild(cleardiv);
    var labeldiv = document.createElement("div");
    labeldiv.innerHTML = "Diffuse&nbsp;&nbsp;&nbsp;";
    labeldiv.style.cssFloat="left";
    this.infoNode.appendChild(labeldiv);
    this.diffuseSlider = document.createElement("input");
    this.diffuseSlider.setAttribute("type","range");
    this.infoNode.appendChild(this.diffuseSlider);
    var diffR = toHex(diffuseValues["r"]*255);
    var diffG = toHex(diffuseValues["g"]*255);
    var diffB = toHex(diffuseValues["b"]*255);
    var diffuseColour = diffR+diffG+diffB;
    var diffuseColourPickerInput = document.createElement("input");
    diffuseColourPickerInput.setAttribute("value",diffuseColour);
    var diffuseUUID = guid();
    diffuseColourPickerInput.setAttribute("id",diffuseUUID);
    diffuseColourPickerInput.style.width = "30px";
    diffuseColourPickerInput.style.height = "10px";
    this.infoNode.appendChild(diffuseColourPickerInput);
    var diffuseColourPicker = new jscolor(diffuseUUID,pickerOptions);

    var breakTag =  document.createElement("br");
    this.infoNode.appendChild(breakTag);
    var shadowButton =  document.createElement("input");
    var shadowLabel = document.createElement("span");
    shadowLabel.innerHTML = "Shadows";
    shadowButton.type = "checkbox";
    this.infoNode.appendChild(shadowLabel);
    this.infoNode.appendChild(shadowButton);

    var breakTag2 =  document.createElement("br");
    this.infoNode.appendChild(breakTag2);
    var depthButton =  document.createElement("input");
    var depthLabel = document.createElement("span");
    depthLabel.innerHTML = "Show shadow depth buffer (debugging)";
    depthButton.type = "checkbox";
    this.infoNode.appendChild(depthLabel);
    this.infoNode.appendChild(depthButton);

    shadowButton.onclick = function(){
        console.log("Click! "+shadowButton.checked);
        self.gl.doShadow = shadowButton.checked;
        self.gl.drawScene();
    }
    depthButton.onclick = function(){
        console.log("Click! "+depthButton.checked);
        self.gl.doShadowDepthDebug = depthButton.checked;
        self.gl.drawScene();
    }

    function ptorgb(p){
        var intens = p[0];
        intens = Math.max(intens,p[1]);
        intens = Math.max(intens,p[2]);
        var r;
        var g;
        var b;
        if(Math.abs(intens)>1e-6){
            r = p[0]/intens;
            g = p[0]/intens;
            b = p[0]/intens;
        } else {
            r = 1.0;
            g = 1.0;
            b = 1.0;
        }
        return {"r":r,"g":g,"b":b,"intensity":intens};
    }

    this.ambientSlider.value = parseInt(ambientValues["intensity"]*100);
    this.ambientSlider.onchange = function() {
        var rgb = ambientColourPickerInput.style.backgroundColor;
        rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
        self.gl.setAmbientLight(this.value*.01*rgb[1]/255.,this.value*.01*rgb[2]/255.,this.value*.01*rgb[3]/255.);
    }
    this.ambientSlider.oninput = function() {
        var rgb = ambientColourPickerInput.style.backgroundColor;
        rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
        self.gl.setAmbientLight(this.value*.01*rgb[1]/255.,this.value*.01*rgb[2]/255.,this.value*.01*rgb[3]/255.);
    }

    this.specularSlider.value = parseInt(specularValues["intensity"]*100);
    this.specularSlider.onchange = function() {
        var rgb = specularColourPickerInput.style.backgroundColor;
        rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
        self.gl.setSpecularLight(this.value*.01*rgb[1]/255.,this.value*.01*rgb[2]/255.,this.value*.01*rgb[3]/255.);
    }
    this.specularSlider.oninput = function() {
        var rgb = specularColourPickerInput.style.backgroundColor;
        rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
        self.gl.setSpecularLight(this.value*.01*rgb[1]/255.,this.value*.01*rgb[2]/255.,this.value*.01*rgb[3]/255.);
    }

    this.diffuseSlider.value = parseInt(diffuseValues["intensity"]*100);
    this.diffuseSlider.onchange = function() {
        var rgb = diffuseColourPickerInput.style.backgroundColor;
        rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
        self.gl.setDiffuseLight(this.value*.01*rgb[1]/255.,this.value*.01*rgb[2]/255.,this.value*.01*rgb[3]/255.);
    }
    this.diffuseSlider.oninput = function() {
        var rgb = diffuseColourPickerInput.style.backgroundColor;
        rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
        self.gl.setDiffuseLight(this.value*.01*rgb[1]/255.,this.value*.01*rgb[2]/255.,this.value*.01*rgb[3]/255.);
    }

    ambientColourPicker.onFineChange = function() {
        var rgb = ambientColourPickerInput.style.backgroundColor;
        rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
        console.log("ambColourPicker.onchange");
        self.gl.setAmbientLight(self.ambientSlider.value*.01*rgb[1]/255.,self.ambientSlider.value*.01*rgb[2]/255.,self.ambientSlider.value*.01*rgb[3]/255.);
    }
    specularColourPicker.onFineChange = function() {
        var rgb = specularColourPickerInput.style.backgroundColor;
        rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
        console.log("specColourPicker.onchange");
        self.gl.setSpecularLight(self.specularSlider.value*.01*rgb[1]/255.,self.specularSlider.value*.01*rgb[2]/255.,self.specularSlider.value*.01*rgb[3]/255.);
    }
    diffuseColourPicker.onFineChange = function() {
        var rgb = diffuseColourPickerInput.style.backgroundColor;
        rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
        console.log("diffColourPicker.onchange");
        self.gl.setDiffuseLight(self.diffuseSlider.value*.01*rgb[1]/255.,self.diffuseSlider.value*.01*rgb[2]/255.,self.diffuseSlider.value*.01*rgb[3]/255.);
    }

}

MGWebGLLightUI.prototype.actOnEvent = function(event,self){
    var x;
    var y;
    var e = event;
    if (e.pageX || e.pageY) { 
        x = e.pageX;
        y = e.pageY;
    }
    else { 
        x = e.clientX ; 
        y = e.clientY ; 
    } 

    var offset = getOffsetRect(self.infoNode);

    x -= offset.left;
    y -= offset.top;
    self.drawGradient(x,y);

    var x_2d = 1.0*x / self.canvas.width - 0.5;
    var y_2d = 1.0*y / self.canvas.height - 0.5;
    var z = 0.5*Math.sqrt(Math.max(0,1.0-4.0*(x_2d*x_2d + y_2d*y_2d)));

    this.gl.setLightPosition(x_2d,-y_2d,z);
}

MGWebGLLightUI.prototype.drawGradient = function(x,y){
    var self = this;
    var c_x = self.canvas.width/2;
    var c_y = self.canvas.height/2;
    var mysize = self.canvas.height;

    var data = self.imageData.data;

    for(var j=0;j<self.canvas.height;j++){
      for(var i=0;i<self.canvas.width;i++){
         if(((c_x-i)*(c_x-i)+(c_y-j)*(c_y-j))<mysize/2*mysize/2&&((c_x-x)*(c_x-x)+(c_y-y)*(c_y-y))<mysize/2*mysize/2){
           var dist = 1.0-parseFloat((x-i)*(x-i)+(y-j)*(y-j))/(mysize*mysize);
           if(dist>1) dist = 1;
           dist = dist * dist * dist * dist;
           data[j*self.canvas.width*4+i*4]   = parseInt(dist * 255);
           data[j*self.canvas.width*4+i*4+1] = parseInt(dist * 255);
           data[j*self.canvas.width*4+i*4+2] = parseInt(dist * 255);
           data[j*self.canvas.width*4+i*4+3] = 255;
         }
      }
    }
    this.ctx.putImageData(self.imageData,0,0,0,0,this.canvas.width, this.canvas.height);
}

MGWebGLLightUI.prototype.doMouseMove = function(event,self){
    if(!self.mouseDown){
        return;
    }
    self.actOnEvent(event,self);
}

MGWebGLLightUI.prototype.doMouseUp = function(event,self){
    self.mouseDown = false;
}

MGWebGLLightUI.prototype.doMouseDown = function(event,self){
    self.mouseDown = true;
    self.mouseDownButton = event.button;
    self.actOnEvent(event,self);
}

function MGWebGLSimpleUI(gl,div,doDrawMenus) {

  if(div){
    this.infoNode = document.getElementById(div);
    console.log("using different div",div);
  } else {
    console.log("using existing div");
    var uiWidth = 150;
    this.infoNode = gl.infoNode;
    gl.leftdiv.width = parseInt(gl.div.style.width)-uiWidth;
    gl.leftdiv.style.maxWidth = parseInt(gl.div.style.width)-uiWidth+"px";
    gl.rightdiv.style.width = uiWidth+"px";
    
    gl.canvas.width = getDeviceScale()*parseInt(gl.div.style.width);
  }

  this.doDrawMenus = true;
  if(typeof(doDrawMenus) != "undefined") {
      this.doDrawMenus = doDrawMenus;
  }

  this.gl = gl;
  this.originNode = null;
  this.scaleNode = null;
  this.orientNode = null;
  this.sizeNode = null;

  var self = this;

  gl.div.addEventListener("fogChanged", 
    function(e) {
      self.redrawFogSlider();
    },
    false
  );
  gl.div.addEventListener("clipChanged", 
    function(e) {
      console.log("e:");
      console.log(e);
      self.redrawClipSlider(e.detail.slab);
    },
    false
  );
  gl.div.addEventListener("statusChanged", 
    function(e) {
      self.drawVisibilityMenu();
    },
    false
  );
  gl.div.addEventListener("viewChanged", 
    function(e) {
      self.update();
    },
    false
  );

  this.fogClipDiv = document.createElement("div");

  var divName = gl.div.getAttribute("id");

  var fogText = document.createElement("div");
  fogText.innerHTML = "Fog";
  this.fogClipDiv.appendChild(fogText);
  var fogSlider = document.createElement("div");
  fogSlider.setAttribute("id",divName+"fog-slider-range");
  this.fogClipDiv.appendChild(fogSlider);

  var clipText = document.createElement("div");
  clipText.innerHTML = "Clip";
  this.fogClipDiv.appendChild(clipText);
  var clipSlider = document.createElement("div");
  clipSlider.setAttribute("id",divName+"clip-slider-range");
  this.fogClipDiv.appendChild(clipSlider);

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

    this.infoNode.appendChild(this.fogClipDiv);

    if(!document.getElementById("modal-1")){
        window.colourBrowserFirer = null;
        var colourInput = document.createElement("input");
        colourInput.setAttribute("class","modal-state");
        //FIXME - This id should be something which is unlikely to collide with something from other framework, etc.
        colourInput.setAttribute("id","modal-1");
        colourInput.setAttribute("type","checkbox");

        document.body.appendChild(colourInput);

        var colourDiv = document.createElement("div");
        colourDiv.setAttribute("class","modal");

        var colourLabel = document.createElement("label");
        colourLabel.setAttribute("class","modal__bg");
        colourLabel.setAttribute("for","modal-1");

        colourDiv.appendChild(colourLabel);

        var colourDivInner = document.createElement("div");
        colourDivInner.setAttribute("class","modal__inner");

        var colourLabelClose = document.createElement("label");
        colourLabelClose.setAttribute("class","modal__close");
        colourLabelClose.setAttribute("for","modal-1");

        colourDivInner.appendChild(colourLabelClose);

        var colourP = document.createElement("p");
        //colourP.setAttribute("id","cats");
        colourP.innerHTML = "Choose colour:";

        colourDivInner.appendChild(colourP);

        var colourPicker = document.createElement("input");
        colourPicker.setAttribute("class","jscolor");
        colourPicker.setAttribute("value","ab2567");

        colourDivInner.appendChild(colourPicker);

        colourDiv.appendChild(colourDivInner);

        var aBreak = document.createElement("br");
        colourDivInner.appendChild(aBreak);
        var aBreak = document.createElement("br");
        colourDivInner.appendChild(aBreak);

        var colourButtonApply = document.createElement("input");
        colourButtonApply.setAttribute("type","button");
        colourButtonApply.setAttribute("value","Apply");
        colourDivInner.appendChild(colourButtonApply);
        colourButtonApply.onclick = function() {
             colourLabelClose.click();
             console.log(window.colourBrowserFirer);
             console.log(colourPicker.value);
             window.colourBrowserFirer["gl"].updateDisplayObject(window.colourBrowserFirer["num"],window.colourBrowserFirer["prop"],"one_colour_"+colourPicker.value);
        }

        var colourButtonCancel = document.createElement("input");
        colourButtonCancel.setAttribute("type","button");
        colourButtonCancel.setAttribute("value","Cancel");
        colourDivInner.appendChild(colourButtonCancel);
        colourButtonCancel.onclick = function() {
             colourLabelClose.click();
        }

        document.body.appendChild(colourDiv);
        window.colourLabelActivate = document.createElement("label");
        window.colourLabelActivate.setAttribute("class","btn");
        window.colourLabelActivate.setAttribute("for","modal-1");
        window.colourLabelActivate.innerHTML = "";
        document.body.appendChild(window.colourLabelActivate);

    }

    if(!document.getElementById("modal-2")){
        window.selectionBrowserFirer = null;
        var selectionInput = document.createElement("input");
        selectionInput.setAttribute("class","modal-state");
        //FIXME - This id should be something which is unlikely to collide with something from other framework, etc.
        selectionInput.setAttribute("id","modal-2");
        selectionInput.setAttribute("type","checkbox");

        document.body.appendChild(selectionInput);

        var selectionDiv = document.createElement("div");
        selectionDiv.setAttribute("class","modal");

        var selectionLabel = document.createElement("label");
        selectionLabel.setAttribute("class","modal__bg");
        selectionLabel.setAttribute("for","modal-2");

        selectionDiv.appendChild(selectionLabel);

        var selectionDivInner = document.createElement("div");
        selectionDivInner.setAttribute("class","modal__inner");

        var selectionLabelClose = document.createElement("label");
        selectionLabelClose.setAttribute("class","modal__close");
        selectionLabelClose.setAttribute("for","modal-2");

        selectionDivInner.appendChild(selectionLabelClose);

        var selectionP = document.createElement("p");
        //selectionP.setAttribute("id","cats");
        selectionP.innerHTML = "Choose selection:";

        selectionDivInner.appendChild(selectionP);

        window.selectionPicker = document.createElement("input");

        selectionDivInner.appendChild(window.selectionPicker);

        selectionDiv.appendChild(selectionDivInner);

        var aBreak = document.createElement("br");
        selectionDivInner.appendChild(aBreak);
        var aBreak = document.createElement("br");
        selectionDivInner.appendChild(aBreak);


        window.atomSelectionOuterDiv = document.createElement("div");
        window.atomSelectionOuterDiv.setAttribute("id","atomselectiondiv");
        window.atomSelectionOuterDiv.style.width = "40%";
        window.atomSelectionOuterDiv.style.maxHeight = "50%";
        window.atomSelectionOuterDiv.style.float = "left";
        atomSelectionOuterDiv.style.overflow = "auto";
        selectionDivInner.appendChild(window.atomSelectionOuterDiv);

        var rightcol = document.createElement("div");
        rightcol.style.float = "right";
        rightcol.style.width = "60%";
        rightcol.style.maxHeight = "50%";
        rightcol.style.overflow = "auto";

        var rightcolTop = document.createElement("div");
        rightcolTop.style.padding = "4px";
        var rightcolBot = document.createElement("div");
        //rightcolBot.style.backgroundColor = "red";
        rightcol.appendChild(rightcolTop);
        rightcol.appendChild(rightcolBot);

        var elementsDiv = document.createElement("span");
        elementsDiv.innerHTML = "Elements";
        elementsDiv.setAttribute("class","tabselected");
        var residuesDiv = document.createElement("span");
        residuesDiv.innerHTML = "Res. types";
        residuesDiv.setAttribute("class","tab");
        var sseDiv = document.createElement("span");
        sseDiv.innerHTML = "Sec. str.";
        sseDiv.setAttribute("class","tab");
        var modelsDiv = document.createElement("span");
        modelsDiv.innerHTML = "Models";
        modelsDiv.setAttribute("class","tab");
        var altLocDiv = document.createElement("span");
        altLocDiv.innerHTML = "Alt. loc.";
        altLocDiv.setAttribute("class","tab");
        rightcolTop.appendChild(elementsDiv);
        rightcolTop.appendChild(residuesDiv);
        rightcolTop.appendChild(sseDiv);
        rightcolTop.appendChild(modelsDiv);
        rightcolTop.appendChild(altLocDiv);
        
        window.elementsGuiDiv = document.createElement("div");
        window.residuesGuiDiv = document.createElement("div");
        window.sseGuiDiv = document.createElement("div");
        window.modelsGuiDiv = document.createElement("div");
        window.altLocGuiDiv = document.createElement("div");
        window.elementsGuiDiv.style.minHeight = "400px";
        window.residuesGuiDiv.style.minHeight = "400px";
        window.sseGuiDiv.style.minHeight = "400px";
        window.modelsGuiDiv.style.minHeight = "400px";
        window.altLocGuiDiv.style.minHeight = "400px";
        //window.elementsGuiDiv.style.backgroundColor = "green";
        //window.residuesGuiDiv.style.backgroundColor = "blue";
        //window.sseGuiDiv.style.backgroundColor = "magenta";
        //window.modelsGuiDiv.style.backgroundColor = "cyan";
        //window.altLocGuiDiv.style.backgroundColor = "yellow";
        rightcolBot.appendChild(window.elementsGuiDiv);
        rightcolBot.appendChild(window.residuesGuiDiv);
        rightcolBot.appendChild(window.sseGuiDiv);
        rightcolBot.appendChild(window.modelsGuiDiv);
        rightcolBot.appendChild(window.altLocGuiDiv);
        window.elementsGuiDiv.style.display = "block";
        window.residuesGuiDiv.style.display = "none";
        window.sseGuiDiv.style.display = "none";
        window.modelsGuiDiv.style.display = "none";
        window.altLocGuiDiv.style.display = "none";

        elementsDiv.onclick = function(){
          window.elementsGuiDiv.style.display = "block";
          window.residuesGuiDiv.style.display = "none";
          window.sseGuiDiv.style.display = "none";
          window.modelsGuiDiv.style.display = "none";
          window.altLocGuiDiv.style.display = "none";
          elementsDiv.setAttribute("class","tabselected");
          residuesDiv.setAttribute("class","tab");
          sseDiv.setAttribute("class","tab");
          modelsDiv.setAttribute("class","tab");
          altLocDiv.setAttribute("class","tab");
        }

        residuesDiv.onclick = function(){
          window.elementsGuiDiv.style.display = "none";
          window.residuesGuiDiv.style.display = "block";
          window.sseGuiDiv.style.display = "none";
          window.modelsGuiDiv.style.display = "none";
          window.altLocGuiDiv.style.display = "none";
          elementsDiv.setAttribute("class","tab");
          residuesDiv.setAttribute("class","tabselected");
          sseDiv.setAttribute("class","tab");
          modelsDiv.setAttribute("class","tab");
          altLocDiv.setAttribute("class","tab");
        }

        sseDiv.onclick = function(){
          window.elementsGuiDiv.style.display = "none";
          window.residuesGuiDiv.style.display = "none";
          window.sseGuiDiv.style.display = "block";
          window.modelsGuiDiv.style.display = "none";
          window.altLocGuiDiv.style.display = "none";
          elementsDiv.setAttribute("class","tab");
          residuesDiv.setAttribute("class","tab");
          sseDiv.setAttribute("class","tabselected");
          modelsDiv.setAttribute("class","tab");
          altLocDiv.setAttribute("class","tab");
        }

        modelsDiv.onclick = function(){
          window.elementsGuiDiv.style.display = "none";
          window.residuesGuiDiv.style.display = "none";
          window.sseGuiDiv.style.display = "none";
          window.modelsGuiDiv.style.display = "block";
          window.altLocGuiDiv.style.display = "none";
          elementsDiv.setAttribute("class","tab");
          residuesDiv.setAttribute("class","tab");
          sseDiv.setAttribute("class","tab");
          modelsDiv.setAttribute("class","tabselected");
          altLocDiv.setAttribute("class","tab");
        }

        altLocDiv.onclick = function(){
          window.elementsGuiDiv.style.display = "none";
          window.residuesGuiDiv.style.display = "none";
          window.sseGuiDiv.style.display = "none";
          window.modelsGuiDiv.style.display = "none";
          window.altLocGuiDiv.style.display = "block";
          elementsDiv.setAttribute("class","tab");
          residuesDiv.setAttribute("class","tab");
          sseDiv.setAttribute("class","tab");
          modelsDiv.setAttribute("class","tab");
          altLocDiv.setAttribute("class","tabselected");
        }

        selectionDivInner.appendChild(rightcol);

        var cleardiv = document.createElement("div");
        cleardiv.style.clear = "both";
        selectionDivInner.appendChild(cleardiv);

        var aBreak = document.createElement("br");
        selectionDivInner.appendChild(aBreak);

        var selectionButtonApply = document.createElement("input");
        selectionButtonApply.setAttribute("type","button");
        selectionButtonApply.setAttribute("value","Select shown atoms");
        selectionDivInner.appendChild(selectionButtonApply);
        selectionButtonApply.onclick = function() {
             selectionLabelClose.click();
             console.log(window.selectionBrowserFirer);
             console.log(window.selectionPicker.value);
             window.selectionBrowserFirer["gl"].updateDisplayObject(window.selectionBrowserFirer["num"],window.selectionBrowserFirer["prop"],window.selectionPicker.value);
        }

        var neighboursButtonApply = document.createElement("input");
        neighboursButtonApply.setAttribute("type","button");
        neighboursButtonApply.setAttribute("value","Select neighbours of shown atoms");
        selectionDivInner.appendChild(neighboursButtonApply);
        neighboursButtonApply.onclick = function() {
             selectionLabelClose.click();
             console.log(window.selectionBrowserFirer);
             console.log(window.selectionPicker.value);
             window.selectionBrowserFirer["gl"].updateDisplayObject(window.selectionBrowserFirer["num"],window.selectionBrowserFirer["prop"],"complex_neighb "+window.selectionPicker.value);
        }

        var selectionButtonCancel = document.createElement("input");
        selectionButtonCancel.setAttribute("type","button");
        selectionButtonCancel.setAttribute("value","Cancel");
        selectionDivInner.appendChild(selectionButtonCancel);
        selectionButtonCancel.onclick = function() {
             selectionLabelClose.click();
        }

        document.body.appendChild(selectionDiv);
        window.selectionLabelActivate = document.createElement("label");
        window.selectionLabelActivate.setAttribute("class","btn");
        window.selectionLabelActivate.setAttribute("for","modal-2");
        window.selectionLabelActivate.innerHTML = "";
        document.body.appendChild(window.selectionLabelActivate);

    }


    /*
    var menu = document.createElement("select");
    var menuEntry = document.createElement("option");
    menuEntry.value = "An entry";
    menuEntry.innerHTML = "Colour picker test";
    menu.appendChild(menuEntry);
    var menuEntry = document.createElement("option");
    menuEntry.value = "Another entry";
    menuEntry.innerHTML = "Colour picker test 2";
    menu.onchange = function() {
        window.colourBrowserFirer = null;
        console.log(window.colourLabelActivate);
        window.colourLabelActivate.click();
    }
    menu.appendChild(menuEntry);

    this.infoNode.appendChild(menu);
    */

}

MGWebGLSimpleUI.prototype.fillElementsSelectionWidget = function(elements){
    while (window.elementsGuiDiv.firstChild) {
        window.elementsGuiDiv.removeChild(window.elementsGuiDiv.firstChild);
    }
 
    function getElementsSelection(){
         var checkboxes = window.elementsGuiDiv.querySelectorAll("input[type=checkbox]");
         var selBoxes = "";
         for(var icb=0;icb<checkboxes.length;icb++){
             if(checkboxes[icb].checked){
                 selBoxes += checkboxes[icb].getAttribute("data-value")+",";
             }
         }
         if(selBoxes.length>0){
             return selBoxes.substring(0,selBoxes.length-1);
         }
         return "";
    }

    for(var i=0;i<elements.length;i++){
        if(elements[i]!=="."){
            var div = document.createElement("div");
            var cb = document.createElement("input");
            cb.setAttribute("type", "checkbox");
            cb.setAttribute("data-value", elements[i]);
            div.innerHTML = elements[i];
            div.insertBefore(cb,div.firstChild);
            cb.onchange = function(){
                var selEl = getElementsSelection();
                if(selEl.length>0){
                    window.selectionPicker.value = "/*/*/*/*["+selEl+"]:*";
                } else {
                    window.selectionPicker.value = "";
                }
            }
            window.elementsGuiDiv.appendChild(div);
        }
    }
}

MGWebGLSimpleUI.prototype.fillRestypesSelectionWidget = function(restypes){
    while (window.residuesGuiDiv.firstChild) {
        window.residuesGuiDiv.removeChild(window.residuesGuiDiv.firstChild);
    }

    function getResiduesSelection(){
         var checkboxes = window.residuesGuiDiv.querySelectorAll("input[type=checkbox]");
         var selBoxes = "";
         for(var icb=0;icb<checkboxes.length;icb++){
             if(checkboxes[icb].checked){
                 selBoxes += checkboxes[icb].getAttribute("data-value")+",";
             }
         }
         if(selBoxes.length>0){
             return selBoxes.substring(0,selBoxes.length-1);
         }
         return "";
    }

    for(var i=0;i<restypes.length;i++){
        //console.log(restypes[i]);
        var div = document.createElement("div");
        var cb = document.createElement("input");
        cb.setAttribute("type", "checkbox");
        cb.setAttribute("data-value", restypes[i]);
        div.innerHTML = restypes[i];
        div.insertBefore(cb,div.firstChild);
        cb.onchange = function(){
            var selRes = getResiduesSelection();
            if(selRes.length>0){
                window.selectionPicker.value = "/*/*/("+selRes+")/*:*";
            } else {
                window.selectionPicker.value = "";
            }
        }
        window.residuesGuiDiv.appendChild(div);
    }
}

MGWebGLSimpleUI.prototype.fillModelsSelectionWidget = function(models){
    while (window.modelsGuiDiv.firstChild) {
        window.modelsGuiDiv.removeChild(window.modelsGuiDiv.firstChild);
    }

    function getModelsSelection(){
         var checkboxes = window.modelsGuiDiv.querySelectorAll("input[type=checkbox]");
         var selBoxes = "";
         for(var icb=0;icb<checkboxes.length;icb++){
             if(checkboxes[icb].checked){
                 selBoxes += checkboxes[icb].getAttribute("data-value")+",";
             }
         }
         if(selBoxes.length>0){
             return selBoxes.substring(0,selBoxes.length-1);
         }
         return "";
    }

    for(var i=0;i<models.length;i++){
        var div = document.createElement("div");
        var cb = document.createElement("input");
        cb.setAttribute("type", "checkbox");
        cb.setAttribute("data-value", models[i]);
        div.innerHTML = models[i];
        div.insertBefore(cb,div.firstChild);
        cb.onchange = function(){
            var selMod = getModelsSelection().split(",");
            if(selMod.length>0){
                var sel = "";
                for(var j=0;j<selMod.length;j++){
                    sel += "/" + selMod[j] +" or ";
                }
                sel = sel.substring(0,sel.length-4);
                window.selectionPicker.value = sel;
            } else {
                window.selectionPicker.value = "";
            }
        }
        window.modelsGuiDiv.appendChild(div);
    }
}

MGWebGLSimpleUI.prototype.fillAltLocsSelectionWidget = function(altlocs){
    while (window.altLocGuiDiv.firstChild) {
        window.altLocGuiDiv.removeChild(window.altLocGuiDiv.firstChild);
    }

    function getAltLocsSelection(){
         var checkboxes = window.altLocGuiDiv.querySelectorAll("input[type=checkbox]");
         var selBoxes = "";
         for(var icb=0;icb<checkboxes.length;icb++){
             if(checkboxes[icb].checked){
                 selBoxes += checkboxes[icb].getAttribute("data-value")+",";
             }
         }
         if(selBoxes.length>0){
             return selBoxes.substring(0,selBoxes.length-1);
         }
         return "";
    }

    for(var i=0;i<altlocs.length;i++){
        var div = document.createElement("div");
        var cb = document.createElement("input");
        cb.setAttribute("type", "checkbox");
        cb.setAttribute("data-value", altlocs[i]);
        div.innerHTML = altlocs[i];
        div.insertBefore(cb,div.firstChild);
        cb.onchange = function(){
            var selLoc = getAltLocsSelection().split(",");
            if(selLoc.length>0){
                var sel = "";
                for(var j=0;j<selLoc.length;j++){
                    sel += "/*/*/*/*:" + selLoc[j] +" or ";
                }
                sel = sel.substring(0,sel.length-4);
                window.selectionPicker.value = sel;
            } else {
                window.selectionPicker.value = "";
            }
        }
        window.altLocGuiDiv.appendChild(div);
    }
}

MGWebGLSimpleUI.prototype.fillSSESelectionWidget = function(sses){
    console.log(sses);
    while (window.sseGuiDiv.firstChild) {
        window.sseGuiDiv.removeChild(window.sseGuiDiv.firstChild);
    }

    function getSSESelection(){
         var checkboxes = window.sseGuiDiv.querySelectorAll("input[type=checkbox]");
         var selBoxes = "";
         for(var icb=0;icb<checkboxes.length;icb++){
             if(checkboxes[icb].checked){
                 selBoxes += checkboxes[icb].getAttribute("data-value")+",";
             }
         }
         if(selBoxes.length>0){
             return selBoxes.substring(0,selBoxes.length-1);
         }
         return "";
    }

    for(var i=0;i<sses.length;i++){
        var div = document.createElement("div");
        var cb = document.createElement("input");
        cb.setAttribute("type", "checkbox");
        cb.setAttribute("data-value", sses[i]);
        div.innerHTML = sses[i];
        div.insertBefore(cb,div.firstChild);
        cb.onchange = function(){
            var selMod = getSSESelection().split(",");
            if(selMod.length>0){
                var sel = "";
                for(var j=0;j<selMod.length;j++){
                    sel += selMod[j].substr(2) +" or ";
                }
                sel = sel.substring(0,sel.length-4);
                console.log(sel);
                window.selectionPicker.value = sel;
            } else {
                window.selectionPicker.value = "";
            }
        }
        window.sseGuiDiv.appendChild(div);
    }
}

MGWebGLSimpleUI.prototype.fillAtomSelectionWidget = function(atoms){
    //var start = performance.now();
    while (window.atomSelectionOuterDiv.firstChild) {
        window.atomSelectionOuterDiv.removeChild(window.atomSelectionOuterDiv.firstChild);
    }
    window.atomSelectionSelect = document.createElement("select");
    window.atomSelectionSelect.setAttribute("id","atomselectionselect");
    window.atomSelectionSelect.setAttribute("multiple","multiple");
    window.atomSelectionOuterDiv.appendChild(window.atomSelectionSelect);
    var iatom = 0;
    for(im=0;im<atoms.length;im++){
        for(ic=0;ic<atoms[im].chains.length;ic++){
            for(ires=0;ires<atoms[im].chains[ic].residues.length;ires++){
                for(iat=0;iat<atoms[im].chains[ic].residues[ires].atoms.length;iat++){
                    if(atoms[im].chains[ic].residues[ires].getSeqID()!=="."){
                        iatom++;
                        //console.log(atoms[im][ic].getChainID()+"/"+atoms[im][ic].residues[ires].getSeqID()+"/"+atoms[im][ic].residues[ires].atoms[iat]["_atom_site.label_atom_id"]);
                        var atomID = atoms[im].chains[ic].residues[ires].atoms[iat]["_atom_site.label_atom_id"];
                        var resID = atoms[im].chains[ic].getChainID()+"/"+atoms[im].chains[ic].residues[ires].getSeqID()+"("+atoms[im].chains[ic].residues[ires].getName()+")";
                        var theOption = document.createElement("option");
                        theOption.setAttribute("value",atomID);
                        theOption.setAttribute("data-section",resID);
                        theOption.setAttribute("data-index",""+iatom);
                        theOption.innerHTML = resID+"/"+atomID;
                        window.atomSelectionSelect.appendChild(theOption);
                    }
                }
            }
        }
    }
    /*
    var end = performance.now();
    var time = end - start;
    console.log('Execution time: ' + time);
    */
}

MGWebGLSimpleUI.prototype.redrawFogSlider = function(){
  var gl = this.gl;
  var divName = gl.div.getAttribute("id");
                    $(function() {
                     $( "#"+divName+"fog-slider-range" ).slider({
                       range: true,
                       min: 0,
                       max: 1000,
                       values: [ gl.gl_fog_start, gl.gl_fog_end ],
                       slide: function( event, ui ) {
                       gl.set_fog_range(ui.values[ 0 ],ui.values[ 1 ],true);
                       }
                     });
                    });
}

MGWebGLSimpleUI.prototype.redrawClipSlider = function(slab){
  console.log("slab:");
  console.log(slab);
  var gl = this.gl;
  var divName = gl.div.getAttribute("id");
                    $(function() {
                     $( "#"+divName+"clip-slider-range" ).slider({
                       range: true,
                       min: -250,
                       max: 250,
                       values: [ -slab[0]*0.5-slab[1], slab[0]*0.5-slab[1]],
                       slide: function( event, ui ) {
                       gl.set_clip_range(ui.values[ 0 ],ui.values[ 1 ],true);
                       }
                     });
                    });
}

MGWebGLSimpleUI.prototype.handleVisCBClick = function(cb){
    //console.log(cb);
    name = cb.getAttribute("idx");
    //console.log(name);
    //console.log(cb.checked);
    if(cb.checked){
        this.gl.displayBuffers[parseInt(name)].visible = true;
    } else {
        this.gl.displayBuffers[parseInt(name)].visible = false;
    }
    this.gl.drawScene();
}

MGWebGLSimpleUI.prototype.update = function(){
    if(this.originNode){
        this.originNode.nodeValue = this.gl.origin[0].toFixed(3)+" "+this.gl.origin[1].toFixed(3)+" "+this.gl.origin[2].toFixed(3);
    }
    if(this.scaleNode){
        this.scaleNode.nodeValue = this.gl.zoom.toFixed(4);
    }
    if(this.orientNode){
        this.orientNode.nodeValue = this.gl.myQuat[0].toFixed(3)+" "+this.gl.myQuat[1].toFixed(3)+" "+this.gl.myQuat[2].toFixed(3)+" "+this.gl.myQuat[3].toFixed(3);
    }
    if(this.sizeNode){
        this.sizeNode.nodeValue = this.gl.gl.viewportWidth+" "+this.gl.gl.viewportHeight;
    }
}

MGWebGLSimpleUI.prototype.drawVisibilityMenu = function(){
    var myNode = this.infoNode;

    while (myNode.firstChild) {
        myNode.removeChild(myNode.firstChild);
    }

    myNode.appendChild(this.fogClipDiv);

    var formNode = document.createElement('form');
    myNode.appendChild(formNode);

    var self = this;

    function recurseAddEntries(menu,menuEntry,menuDesc){
        var menuEntry = document.createElement("optgroup");
        menuEntry.setAttribute("label",menuDesc[0]);
        //menuEntry.style.maxWidth = "50px";
        menu.appendChild(menuEntry);
        for(var imenu=1;imenu<menuDesc.length;imenu++){
            if( typeof menuDesc[imenu][1] === 'string' ) {
                var menuEntry = document.createElement("option");
                //menuEntry.style.maxWidth = "50px";
                menuEntry.value= menuDesc[imenu][0];
                menuEntry.innerHTML= menuDesc[imenu][1];
                menuEntry.onclick = function () {
                    console.log("Click 2!!");
                }
                menu.appendChild(menuEntry);
            } else {
                recurseAddEntries(menu,menuEntry,menuDesc[imenu]);
            }
        }
    }


    function drawMenus(){
        var table = document.createElement("table");
        formNode.appendChild(table);
        var currentDataId = "";
        if(typeof(self.gl.mapDisplayOptions)!=="undefined"){
            for (var idx2 = 0; idx2 < self.gl.mapDisplayOptions.length; idx2++){
                var idx = idx2+self.gl.displayOptions.length;
                var underIndex = self.gl.ids[idx].indexOf("_");
                if(underIndex>-1){
                    var dataId = self.gl.ids[idx].substring(0,underIndex);
                    if(dataId!=currentDataId){
                        var tr = document.createElement("tr");
                        table.appendChild(tr);
                        var col0td = document.createElement("td");
                        tr.appendChild(col0td);
                        col0td.setAttribute("colspan","2");
                        var dataName = self.gl.getNameForId(dataId,"MapData");
                        col0td.innerHTML = dataName.substring(0,10);
                        console.log("Map: "+dataId);
                        currentDataId = dataId;
                    }
                }
                var tr = document.createElement("tr");
                table.appendChild(tr);
                var colm1td = document.createElement("td");
                tr.appendChild(colm1td);
                var col0td = document.createElement("td");
                tr.appendChild(col0td);
                var col1td = document.createElement("td");
                tr.appendChild(col1td);
                var col2td = document.createElement("td");
                tr.appendChild(col2td);
                var col3td = document.createElement("td");
                tr.appendChild(col3td);

                var x = document.createElement("INPUT");
                x.setAttribute("type", "checkbox"); 
                x.setAttribute("idx", ""+idx); 
                x.onclick = function() {
                    self.handleVisCBClick(this,self.gl);
                }
                if(self.gl.displayBuffers[idx].visible){
                    x.setAttribute('checked',"checked");
                }

                colm1td.appendChild(x);
            }
        }

        currentDataId = "";
        for (var idx = 0; idx < self.gl.displayOptions.length; idx++){
            var underIndex = self.gl.ids[idx].indexOf("_");
            if(underIndex>-1){
                var dataId = self.gl.ids[idx].substring(0,underIndex);
                if(dataId!=currentDataId){
                    var tr = document.createElement("tr");
                    table.appendChild(tr);
                    var col0td = document.createElement("td");
                    tr.appendChild(col0td);
                    col0td.setAttribute("colspan","2");
                    var dataName = self.gl.getNameForId(dataId,"MolData");
                    //col0td.innerHTML = dataId.substring(0,10)+"..."; // Hmm we want PDB name if possible. Martin's request.
                    col0td.innerHTML = dataName.substring(0,10);
                    currentDataId = dataId;
                    var addDispObjButton = document.createElement("input");
                    addDispObjButton.setAttribute("type","button");
                    addDispObjButton.setAttribute("value","+");
                    var col1td = document.createElement("td");
                    col1td.setAttribute("colspan","3");
                    tr.appendChild(col1td);
                    col1td.appendChild(addDispObjButton);
                    (function(id) {
                    addDispObjButton.onclick = function() {
                        self.gl.addDisplayObject(id);
                    }
                    })(dataId);
                }
            }
            var tr = document.createElement("tr");
            table.appendChild(tr);
            var colm1td = document.createElement("td");
            tr.appendChild(colm1td);
            var col0td = document.createElement("td");
            tr.appendChild(col0td);
            var col1td = document.createElement("td");
            tr.appendChild(col1td);
            var col2td = document.createElement("td");
            tr.appendChild(col2td);
            var col3td = document.createElement("td");
            tr.appendChild(col3td);

            var x = document.createElement("INPUT");
            x.setAttribute("type", "checkbox"); 
            x.setAttribute("idx", ""+idx); 
            x.onclick = function() {
                self.handleVisCBClick(this,self.gl);
            }
            if(self.gl.displayBuffers[idx].visible){
                x.setAttribute('checked',"checked");
            }

            var menu = document.createElement("select");
            menu.style.maxWidth = "30px";
            var menuEntry = document.createElement("option");
            menuEntry.value = "dispobjtoolsdummy";
            menuEntry.setAttribute("disabled","disabled");
            menuEntry.innerHTML = "Tools";
            menu.appendChild(menuEntry);
            var menuEntryCentre = document.createElement("option");
            menuEntryCentre.value = "centreon";
            menuEntryCentre.innerHTML = "Centre on";
            menu.appendChild(menuEntryCentre);
            var menuEntryDelete = document.createElement("option");
            menuEntryDelete.value = "delete";
            menuEntryDelete.innerHTML = "Delete";
            menu.appendChild(menuEntryDelete);

            (function(num) {
                 menu.onchange = function(){
                     if(this.value=="centreon"){
                         console.log("Centre on "+num);
                         self.gl.centreOn(num);
                     } else if(this.value=="delete"){
                         console.log("Delete "+num);
                         self.gl.removeDisplayObject(num);
                     }
                 this.options.selectedIndex = 0;
                 }
            })(idx);

            colm1td.appendChild(x);
            col0td.appendChild(menu);
            menu.options.selectedIndex = 0;

            //formNode.appendChild(document.createTextNode(this.gl.displayBuffers[idx].name_label));

            for(var propName in self.gl.displayOptions[idx]){
                if(typeof(self.gl.displayOptions[idx][propName]) != "undefined") {
                    var menu = document.createElement("select");
                    /*
                    if(property==="Colour"){
                    } else if ("selection"){
                    } else if ("Style"){
                    }
                    */
                    (function(num,property) {
                    menu.onchange = function(){
                        if(property==="Colour"&&this.value=="browser"){
                             console.log(this.options.selectedIndex);
                             this.options.selectedIndex = this.options.selectedIndex+1;
                             window.colourBrowserFirer = {"gl":self.gl,"num":num,"prop":property};
                             //console.log(window.colourLabelActivate);
                             window.colourLabelActivate.click();
                             //console.log("issued the click, did it work ....?");
                        } else if(this.value=="selectionbrowser"){
                             this.options.selectedIndex = this.options.selectedIndex+1;
                             console.log("Selection browser");
                             window.selectionBrowserFirer = {"gl":self.gl,"num":num,"prop":property};
                             var theID = self.gl.ids[num];
                             var underIndex = theID.indexOf("_");
                             theID = theID.substring(0,underIndex);
                             console.log(theID);
                             var theCIFAtoms = self.gl.cifatoms[theID];
                             console.log(theCIFAtoms);
                             // FIXME - This is way too slow.
                             //try {
                                 self.fillAtomSelectionWidget(self.gl.cifatoms[theID]);
                                 self.fillElementsSelectionWidget(self.gl.elements[theID]);
                                 self.fillRestypesSelectionWidget(self.gl.restypes[theID]);
                                 self.fillModelsSelectionWidget(self.gl.models[theID]);
                                 self.fillAltLocsSelectionWidget(self.gl.altlocs[theID]);
                                 self.fillSSESelectionWidget(self.gl.secstr[theID]);
                                 //var start = performance.now();
                                 var options = { sortable: false, hideSidePanel: true, collapsible: true, startCollapsed: true };
                                 var multi1 = $("select#atomselectionselect").treeMultiselect(options);
                                 var checkboxes =  window.atomSelectionOuterDiv.querySelectorAll("input[type=checkbox]");
                                 for(var icb=0;icb<checkboxes.length;icb++){
                                     (function(num2) {
                                      checkboxes[num2].onchange = function() {
                                      console.log("onchange");
                                      var selection = getAtomSelection(window.atomSelectionOuterDiv);
                                      window.selectionPicker.value = selection;
                                      }
                                      })(icb);
                                 }
/*
                                 var end = performance.now();
                                 var time = end - start;
                                 console.log('Fancy select execution time: ' + time);
*/
                             //} catch(e) {
                                 //console.log("Failed to create atom selector widget");
                             //}
                             window.selectionLabelActivate.click();
                        } else if((property==="Colour"&&this.value=="browserdummy")||this.value=="selectiondummy"){
                            //Do nothing
                        } else {
                            self.gl.updateDisplayObject(num,property,this.value);
                        }
                    }
                    })(idx,propName);
                    for(var imenu=0;imenu<self.gl.displayOptions[idx][propName].length;imenu++){
                        if( typeof self.gl.displayOptions[idx][propName][imenu][1] === 'string' ) {
                            var menuEntry = document.createElement("option");
                            //menuEntry.style.maxWidth = "50px";
                            menuEntry.value= self.gl.displayOptions[idx][propName][imenu][0];
                            menuEntry.innerHTML= self.gl.displayOptions[idx][propName][imenu][1];
                            // Only works on Firefox, doh.
                            menu.appendChild(menuEntry);
                        } else {
                            recurseAddEntries(menu,menuEntry,self.gl.displayOptions[idx][propName][imenu]);
                        }
                    }
                    if(propName==="Colour"){
                        var menuEntry = document.createElement("option");
                        menuEntry.value = "browserdummy";
                        menuEntry.setAttribute("disabled","disabled");
                        menuEntry.innerHTML = "browser colour:";
                        menu.appendChild(menuEntry);
                    }
                    // FIXME - this is bordering on hackery. We need to know if it's a selection menu.
                    if(propName==="Selection"||propName==="Atoms from"||propName==="Atoms to"||propName==="Atoms to draw round"||propName==="In context of atoms"){
                        var menuEntryBrowser = document.createElement("option");
                        menuEntryBrowser.value = "selectionbrowser";
                        menuEntryBrowser.innerHTML = "Selection browser";
                        menu.appendChild(menuEntryBrowser);
                        var menuEntry = document.createElement("option");
                        menuEntry.value = "selectiondummy";
                        menuEntry.setAttribute("disabled","disabled");
                        menuEntry.innerHTML = "browser selection:";
                        menu.appendChild(menuEntry);
                    }
                    menu.style.maxWidth = "30px";
                    if(propName==="Colour"){
                        var textBox = document.createElement("input");
                        textBox.setAttribute("type","text");
                        textBox.style.display = "none";
                        formNode.appendChild(textBox);
                    }
                    menu.options.selectedIndex = -1;
                    if(("Colour" in self.gl.displayOptions[idx])&&("Selection" in self.gl.displayOptions[idx])&&("Style" in self.gl.displayOptions[idx])){
                        if(propName==="Selection"){
                            col1td.appendChild(menu);
                        }
                        if(propName==="Colour"){
                            col2td.appendChild(menu);
                        }
                        if(propName==="Style"){
                            col3td.appendChild(menu);
                        }
                        //console.log("MolDisp");
                    } else if (("Atoms from" in self.gl.displayOptions[idx])&&("Atoms to" in self.gl.displayOptions[idx])&&("Colour" in self.gl.displayOptions[idx])){
                        if(propName==="Atoms from"){
                            col1td.appendChild(menu);
                        }
                        if(propName==="Atoms to"){
                            col2td.appendChild(menu);
                        }
                        if(propName==="Colour"){
                            col3td.appendChild(menu);
                        }
                        //console.log("H-Bonds/Contacts");
                    } else if (("In context of atoms" in self.gl.displayOptions[idx])&&("Atoms to draw round" in self.gl.displayOptions[idx])&&("Colour" in self.gl.displayOptions[idx])){
                        if(propName==="Atoms to draw round"){
                            col1td.appendChild(menu);
                        }
                        if(propName==="In context of atoms"){
                            col2td.appendChild(menu);
                        }
                        if(propName==="Colour"){
                            col3td.appendChild(menu);
                        }
                        //console.log("Surface");
                    } else {
                        console.log("Unknown display options!!!!!!");
                        console.log(self.gl.displayOptions[idx]);
                    }
                }
            }
        }
        var breakNode = document.createElement("br");
        formNode.appendChild(breakNode);
    }

    if(this.doDrawMenus){
        drawMenus();
    }

    /*
    for (var idx = 0; idx < this.gl.displayBuffers.length; idx++){
        var x = document.createElement("INPUT");
        x.setAttribute("type", "checkbox"); 
        x.setAttribute("idx", ""+idx); 
        x.onclick = function() {
            self.handleVisCBClick(this,this.gl);
        }
        if(this.gl.displayBuffers[idx].visible){
            x.setAttribute('checked',"checked");
        }
        formNode.appendChild(x);
        formNode.appendChild(document.createTextNode(this.gl.displayBuffers[idx].name_label));
        var breakNode = document.createElement('br');
        formNode.appendChild(breakNode);
    }
    */

    this.originNode = document.createTextNode(this.gl.origin[0].toFixed(3)+" "+this.gl.origin[1].toFixed(3)+" "+this.gl.origin[2].toFixed(3));
    myNode.appendChild(this.originNode);
    myNode.appendChild(document.createElement('br'));

    this.scaleNode = document.createTextNode(this.gl.zoom.toFixed(4));
    myNode.appendChild(this.scaleNode);
    myNode.appendChild(document.createElement('br'));

    this.orientNode = document.createTextNode(this.gl.myQuat[0].toFixed(3)+" "+this.gl.myQuat[1].toFixed(3)+" "+this.gl.myQuat[2].toFixed(3)+" "+this.gl.myQuat[3].toFixed(3));
    myNode.appendChild(this.orientNode);
    myNode.appendChild(document.createElement('br'));

    this.sizeNode = document.createTextNode(this.gl.viewportWidth+" "+this.gl.viewportHeight);
    myNode.appendChild(this.sizeNode);

}

MGWebGL.prototype.loadJSONFromDiv = function (jsonDiv) {
    var theJSONIsland = document.getElementById(jsonDiv).firstChild;
    this.loadJSONIfReady(JSON.parse(theJSONIsland.data));
}

MGWebGL.prototype.loadJSONIfReady = function (jsondata) {
    var self = this;

    var myVar = setInterval(function(){ myTimer() }, 500);
    function myTimer() {
        console.log(self.ready);
        if(self.ready){
            clearInterval(myVar);
            self.loadJSON(jsondata);
        }
    }
}


MGWebGL.prototype.clearData = function () {
    var self = this;

    self.displayBuffers = [];
    self.liveUpdatingMaps = [];
    self.clickedAtoms = [];
    self.textLabels = [];
    self.ids = [];
    self.displayOptions = [];
    self.mapDisplayOptions = [];
    self.clickedAtoms = [];
}

MGWebGL.prototype.updateBuffers = function (jsondata,theseBuffers) {
    var self = this;
    for(var idat=0;idat<theseBuffers.length;idat++){

        self.currentBufferIdx = self.displayBuffers.indexOf(theseBuffers[idat]);

        theseBuffers[idat].clearBuffers();
        if(idat>=jsondata.norm_tri.length){
            console.log("Too many buffers for data! Send me more");
            break;
        }
         var rssentries=jsondata.norm_tri[idat];
        //console.log(rssentries);
        var norms = rssentries;
        for (var i=0; i<norms.length; i++){
            self.createNormalBuffer(norms[i]);
        }

        var rssentries=jsondata.vert_tri[idat];
        var start = new Date().getTime();
        var tris = rssentries;
        //console.log(rssentries);
        var end = new Date().getTime();
        var time = end - start;

        for (var i=0; i<tris.length; i++){
            self.createVertexBuffer(tris[i]);
        }

        var rssentries=jsondata.idx_tri[idat];
        var idxs = rssentries;
        //console.log(rssentries);

        for (var i=0; i<idxs.length; i++){
            for(var j=0; j< idxs[i].length; j++){
            }
            self.createIndexBuffer(idxs[i]);
        }

        var rssentries=jsondata.col_tri[idat];
        var colours =rssentries;
        //console.log(rssentries);

        for (var i=0; i<colours.length; i++){
            self.createColourBuffer(colours[i]);
        }

        var rssentries=jsondata.prim_types[idat];
        //console.log(rssentries);
        for (var i=0; i<rssentries.length; i++){
            self.displayBuffers[self.currentBufferIdx].bufferTypes.push(rssentries[i]);
        }

        //var thisVis = jsondata.visibility[idat];
        var thisVis = true;
        if(!thisVis){
            self.displayBuffers[self.currentBufferIdx].visible = false;
        }

        //var thisName = jsondata.names[idat];
        self.displayBuffers[self.currentBufferIdx].name_label = "foo";

        //var atoms = jsondata.atoms[idat];
        self.displayBuffers[self.currentBufferIdx].atoms = [];
   }
    self.buildBuffers();
    self.drawScene();
}

MGWebGL.prototype.appendOtherData = function (jsondata,skipRebuild) {
    /* This can be used to *add* arbitrary triangles to a scene. Not much luck with replacing scene by this, yet.
       This currently deals with actual numbers rather than the uuencoded stuff we get from server, but it will
       be changed to handle both.

       Might also be nice to use this as a test-bed for creating more abstract primitives: circles, squares, stars, etc.
     */

    var self = this;
    //console.log("appendOtherData");
    //console.log(jsondata);

    var theseBuffers = [];
    for(var idat=0;idat<jsondata.norm_tri.length;idat++){
        //self.currentBufferIdx = idat;
        self.currentBufferIdx = self.displayBuffers.length;
        self.displayBuffers.push(new DisplayBuffer());
        theseBuffers.push(self.displayBuffers[self.currentBufferIdx]);
        var rssentries=jsondata.norm_tri[idat];
        //console.log(rssentries);
        var norms = rssentries;
        for (var i=0; i<norms.length; i++){
            self.createNormalBuffer(norms[i]);
        }

        var rssentries=jsondata.vert_tri[idat];
        var start = new Date().getTime();
        var tris = rssentries;
        //console.log(rssentries);
        var end = new Date().getTime();
        var time = end - start;

        for (var i=0; i<tris.length; i++){
            self.createVertexBuffer(tris[i]);
        }

        var rssentries=jsondata.idx_tri[idat];
        var idxs = rssentries;
        //console.log(rssentries);

        for (var i=0; i<idxs.length; i++){
            for(var j=0; j< idxs[i].length; j++){
            }
            self.createIndexBuffer(idxs[i]);
        }

        if(typeof(jsondata.radii)!="undefined"){
            if(typeof(jsondata.radii[idat])!="undefined"){
                var rssentries=jsondata.radii[idat];
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"radii");
                }
            }
        }
        
        if(typeof(jsondata.scale_matrices)!="undefined"){
            if(typeof(jsondata.scale_matrices[idat])!="undefined"){
                var rssentries=jsondata.scale_matrices[idat];
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"scale_matrices");
                }
            }
        }

        if(typeof(jsondata.customSplineNormals)!="undefined"){
            if(typeof(jsondata.customSplineNormals[idat])!="undefined"){
                var rssentries=jsondata.customSplineNormals[idat];
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"customSplineNormals");
                }
            }
        }

        if(typeof(jsondata.spline_accu)!="undefined"){
            if(typeof(jsondata.spline_accu[idat])!="undefined"){
                var rssentries=jsondata.spline_accu[idat];
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"spline_accu");
                }
            }
        }
        
        if(typeof(jsondata.accu)!="undefined"){
            if(typeof(jsondata.accu[idat])!="undefined"){
                var rssentries=jsondata.accu[idat];
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"accu");
                }
            }
        }
        
        if(typeof(jsondata.arrow)!="undefined"){
            if(typeof(jsondata.arrow[idat])!="undefined"){
                var rssentries=jsondata.arrow[idat];
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"arrow");
                }
            }
        }

        if(typeof(jsondata.vert_tri_2d)!="undefined"){
            if(typeof(jsondata.vert_tri_2d[idat])!="undefined"){
                var rssentries=jsondata.vert_tri_2d[idat];
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"vert_tri_2d");
                }
            }
        }

        if(typeof(jsondata.font)!="undefined"){
            if(typeof(jsondata.font[idat])!="undefined"){
                var rssentries=jsondata.font[idat];
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"font");
                }
            }
        }

        if(typeof(jsondata.imgsrc)!="undefined"){
            if(typeof(jsondata.imgsrc[idat])!="undefined"){
                var rssentries=jsondata.imgsrc[idat];
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"imgsrc");
                }
            }
        }

        if(typeof(jsondata.sizes)!="undefined"){
            if(typeof(jsondata.sizes[idat])!="undefined"){
                var rssentries=getEncodedData(jsondata.sizes[idat]);
                for (var i=0; i<rssentries.length; i++){
                    self.createSizeBuffer(rssentries[i]);
                }
            }
        }

        if(typeof(jsondata.vertices2d)!="undefined"){
            if(typeof(jsondata.vertices2d[idat])!="undefined"){
                var rssentries=getEncodedData(jsondata.vertices2d[idat]);
                for (var i=0; i<rssentries.length; i++){
                    self.addSupplementaryInfo(rssentries[i],"vertices2d");
                }
            }
        }

        if(typeof(jsondata.symmetry)!="undefined"){
            self.setSymmetryMatrices(jsondata.symmetry);
        }

        var doChickenWire = false;

        if(typeof(jsondata.mapGrids)!="undefined"){
            // FIXME - need to do same as lines 2379-
            console.log("Got a mapgrid!!!!!!!!!!!!!!!!!!!!!!");
            var rssentries=jsondata.mapGrids[idat];
            var mapColour = [0.0,0.0,1.0,1.0];
            var contourLevel = 0.7;
            if(typeof(jsondata.mapColours)!='undefined'){
                if(typeof(jsondata.mapColours[idat])!='undefined'){
                    mapColour=jsondata.mapColours[idat];
                }
            }
            if(typeof(jsondata.mapContourLevels)!='undefined'){
                if(typeof(jsondata.mapContourLevels[idat])!='undefined'){
                    contourLevel=jsondata.mapContourLevels[idat];
                }
            }
            var nCells_x = rssentries.nCells_x;
            var nCells_y = rssentries.nCells_y;
            var nCells_z = rssentries.nCells_z;
            var cellLength_x = rssentries.cellLength_x;
            var cellLength_y = rssentries.cellLength_y;
            var cellLength_z = rssentries.cellLength_z;
            var cell = rssentries.cell;
            var orthMat = rssentries.fracToOrthMat;
            var fracMat = rssentries.orthToFracMat;
            console.log(nCells_x+" "+nCells_y+" "+nCells_z);
            console.log(cellLength_x+" "+cellLength_y+" "+cellLength_z);
            console.log(orthMat);
            console.log(fracMat);
            var map = new CIsoSurface();
            var start = new Date().getTime();
            console.log(rssentries.grid);
            map.GenerateSurfacePartial(rssentries.grid,0.9,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,0,0,0,nCells_y/4,nCells_z/4,doChickenWire);
            var end = new Date().getTime();
            var time = end - start;
            console.log('generate map: ' + time*0.001+"s");
            var start = new Date().getTime();
            var vertices = map.returnVertices(0,0,0);
            console.log(vertices.length);
            var normals = map.returnNormals_new();
            console.log(normals.length);
            var indices = map.returnIndices();
            console.log(indices.length);
            var end = new Date().getTime();
            var time = end - start;
            console.log('copy buffers: ' + time*0.001+"s");

            colours = [];
            for(var i=0;i<vertices.length/3;i++){
                colours.push(0.0);
                colours.push(0.0);
                colours.push(1.0);
                colours.push(1.0);
            }
            //console.log(vertices);
            //console.log(normals);
            //console.log(indices);
            //console.log(colours);
            var mapTriangleData = {"col_tri":[[colours],[],[],[],[],[],[],[]], "norm_tri":[[normals],[],[],[],[],[],[],[]], "vert_tri":[[vertices],[],[],[],[],[],[],[]], "idx_tri":[[indices],[],[],[],[],[],[],[]] , "prim_types":[["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"]] };
            var theseBuffers = self.appendOtherData(mapTriangleData);
            var liveUpdatingMap = {};
            var gridData = rssentries.grid;
            liveUpdatingMap["gridData"] = gridData;
            liveUpdatingMap["nCells_x"] = nCells_x;
            liveUpdatingMap["nCells_y"] = nCells_y;
            liveUpdatingMap["nCells_z"] = nCells_z;
            liveUpdatingMap["cellLength_x"] = cellLength_x;
            liveUpdatingMap["cellLength_y"] = cellLength_y;
            liveUpdatingMap["cellLength_z"] = cellLength_z;
            liveUpdatingMap["cell"] = cell;
            liveUpdatingMap["orthMat"] = orthMat;
            liveUpdatingMap["fracMat"] = fracMat;
            liveUpdatingMap["theseBuffers"] = theseBuffers;
            liveUpdatingMap["contourLevel"] = contourLevel;
            liveUpdatingMap["mapColour"] = mapColour;
            self.liveUpdatingMaps.push(liveUpdatingMap);
            
        }

        var rssentries=jsondata.col_tri[idat];
        var colours =rssentries;
        //console.log(rssentries);

        for (var i=0; i<colours.length; i++){
            self.createColourBuffer(colours[i]);
        }

        var rssentries=jsondata.prim_types[idat];
        //console.log(rssentries);
        for (var i=0; i<rssentries.length; i++){
            self.displayBuffers[self.currentBufferIdx].bufferTypes.push(rssentries[i]);
        }

        //console.log(jsondata);
        if(typeof(jsondata.visibility)!="undefined"){
            if(typeof(jsondata.visibility[idat])!="undefined"){
                var thisVis = jsondata.visibility[idat];
                self.displayBuffers[self.currentBufferIdx].visible = thisVis;
            }
        } else  {
            // Don't know when this can be triggered ...
            var thisVis = true;
            if(!thisVis){
                self.displayBuffers[self.currentBufferIdx].visible = false;
            }
        }

        //var thisName = jsondata.names[idat];
        self.displayBuffers[self.currentBufferIdx].name_label = "foo";

        //var atoms = jsondata.atoms[idat];
        if(typeof(jsondata.atoms)!=="undefined"){
            self.displayBuffers[self.currentBufferIdx].atoms = jsondata.atoms[idat][0];
        } else {
            self.displayBuffers[self.currentBufferIdx].atoms = [];
        }

    }

    if(typeof(skipRebuild)!=="undefined"&&skipRebuild){
        return theseBuffers;
    }

    self.buildBuffers();
    self.drawScene();
    return theseBuffers;
}

MGWebGL.prototype.appendOtherDataIfReady = function (jsondata) {
    var self = this;

    var myVar = setInterval(function(){ myTimer() }, 500);
    function myTimer() {
        if(self.ready){
            clearInterval(myVar);
            self.appendOtherData(jsondata);
            self.reContourMaps();
        }
    }
}

MGWebGL.prototype.appendOtherJSONData = function (jsondata) {
    //console.log(jsondata);
    //console.log(JSON.parse(jsondata));
    this.appendOtherDataIfReady(JSON.parse(jsondata));
}

MGWebGL.prototype.setFog = function (fog) {
    var self = this;
    console.log(fog);
    self.gl_fog_start = 500+fog[0];
    self.gl_fog_end = 500+fog[1];
    self.drawScene();
    self.div.dispatchEvent(self.fogChangedEvent);
}

MGWebGL.prototype.setSlab = function (slab) {
    var self = this;
    console.log(slab);
    self.gl_clipPlane0[3] = -500.0+slab[1]*0.5+slab[0];
    self.gl_clipPlane1[3] = 500.0+slab[1]*0.5-slab[0];
    self.drawScene();
    self.clipChangedEvent = document.createEvent("CustomEvent");
    self.clipChangedEvent.initCustomEvent('clipChanged',false,false,{"slab":slab});
    console.log("Despatching:");
    console.log(self.clipChangedEvent);
    self.div.dispatchEvent(self.clipChangedEvent);
}

MGWebGL.prototype.loadJSON = function (jsondata) {
    var self = this;
    var startAll = new Date().getTime();

    self.displayBuffers = [];
    self.liveUpdatingMaps = [];
    self.clickedAtoms = [];
    self.textLabels = [];
    self.ids = [];
    self.mapDisplayOptions = [];
    self.clickedAtoms = [];

    //console.log("Length of buffers: "+jsondata.norm_tri.length);
    var divName = self.div.getAttribute("id");

    for(var idat=0;idat<jsondata.norm_tri.length;idat++){
        self.displayBuffers.push(new DisplayBuffer());
        self.currentBufferIdx = idat;
        var rssentries=jsondata.norm_tri[idat];
        var start = new Date().getTime();
        var norms = getEncodedData(rssentries);
        var end = new Date().getTime();
        var time = end - start;
        //console.log('getEncodedData(normals): ' + time*0.001+"s");
        for (var i=0; i<norms.length; i++){
            self.createNormalBuffer(norms[i]);
        }

        if(typeof(jsondata.sizes)!="undefined"){
            var rssentries=getEncodedData(jsondata.sizes[idat]);
            if(typeof(rssentries)!=="undefined"){
                for (var i=0; i<rssentries.length; i++){
                    self.createSizeBuffer(rssentries[i]);
                }
            }
        }
        var rssentries=jsondata.vert_tri[idat];
        var start = new Date().getTime();
        var tris = getEncodedData(rssentries);
        var end = new Date().getTime();
        var time = end - start;
        //console.log('getEncodedData(triangles): ' + time*0.001+"s");

        for (var i=0; i<tris.length; i++){
            self.createVertexBuffer(tris[i]);
        }
        //console.log("Num vertex buffers "+tris.length);

        var rssentries=jsondata.idx_tri[idat];
        var start = new Date().getTime();
        var idxs = getEncodedData(rssentries);
        var end = new Date().getTime();
        var time = end - start;
        //console.log('getEncodedData(indices): ' + time*0.001+"s");

        for (var i=0; i<idxs.length; i++){
            //console.log(i+" "+idxs[i][0]+" "+idxs[i][1]+" "+idxs[i][2]+", "+idxs[i][3]+" "+idxs[i][4]+" "+idxs[i][5]+", "+idxs[i][6]+" "+idxs[i][7]+" "+idxs[i][8]+", "+idxs[i][9]+" "+idxs[i][10]+" "+idxs[i][11]);
            for(var j=0; j< idxs[i].length; j++){
            }
            self.createIndexBuffer(idxs[i]);
        }
        //console.log("Num index buffers "+idxs.length);

        var rssentries=jsondata.col_tri[idat];
        var start = new Date().getTime();
        var colours = getEncodedData(rssentries);
        var end = new Date().getTime();
        var time = end - start;
        //console.log('getEncodedData(colours): ' + time*0.001+"s");

        for (var i=0; i<colours.length; i++){
            self.createColourBuffer(colours[i]);
        }
        //console.log("Num colour buffers "+colours.length);

        var rssentries=jsondata.prim_types[idat];
        for (var i=0; i<rssentries.length; i++){
            self.displayBuffers[self.currentBufferIdx].bufferTypes.push(rssentries[i]);
        }

        var thisVis = jsondata.visibility[idat];
        if(!thisVis){
            self.displayBuffers[self.currentBufferIdx].visible = false;
        }

        var thisName = jsondata.names[idat];
        self.displayBuffers[self.currentBufferIdx].name_label = thisName;

        var atoms = jsondata.atoms[idat];
        self.displayBuffers[self.currentBufferIdx].atoms = atoms;
        if(atoms.length>0){
            //console.log(atoms[0]);
        }

        var thisDisplayClass = jsondata.display_classes[idat];

    }

    var bgcolour = jsondata.background;
    if(bgcolour&&bgcolour.length==3){
        self.background_colour[0] = bgcolour[0]/255.;
        self.background_colour[1] = bgcolour[1]/255.;
        self.background_colour[2] = bgcolour[2]/255.;
    }

    var rssentries=jsondata.origin;
    self.origin = [];
    for (var i=0; i<rssentries.length; i++){
        self.origin.push(-parseFloat(rssentries[i]));
    }
    /*
       var rssentries=jsondata.extent;
       var extent = [];
       self.radius = 50.0;
       for (var i=0; i<rssentries.length; i++){
       extent.push(parseFloat(rssentries[i]));
       }
       if(Math.abs(extent[0]-extent[3])>self.radius){
       self.radius = Math.abs(extent[0]-extent[3]);
       }
       if(Math.abs(extent[1]-extent[4])>self.radius){
       self.radius = Math.abs(extent[0]-extent[3]);
       }
       if(Math.abs(extent[2]-extent[5])>self.radius){
       self.radius = Math.abs(extent[0]-extent[3]);
       }
       self.radius /=2;
     */

    if(typeof(jsondata.symmetry)!="undefined"){
        console.log("Have some symmetry................");
        console.log(jsondata.symmetry);
        self.setSymmetryMatrices(jsondata.symmetry);
    }


    // This turns "grids" from json into isosurfaces. Originally of null size. Recentering recalculates.
    if(typeof(jsondata.mapGrids)!="undefined"){
        for(var idat=0;idat<jsondata.mapGrids.length;idat++){
            var rssentries=jsondata.mapGrids[idat];
            var mapColour = [0.0,0.0,1.0,1.0];
            var contourLevel = 0.7;
            if(typeof(jsondata.mapColours)!='undefined'){
                if(typeof(jsondata.mapColours[idat])!='undefined'){
                    mapColour=jsondata.mapColours[idat];
                }
            }
            if(typeof(jsondata.mapContourLevels)!='undefined'){
                if(typeof(jsondata.mapContourLevels[idat])!='undefined'){
                    contourLevel=jsondata.mapContourLevels[idat];
                }
            }
            var nCells_x = rssentries.nCells_x;
            var nCells_y = rssentries.nCells_y;
            var nCells_z = rssentries.nCells_z;
            var cellLength_x = rssentries.cellLength_x;
            var cellLength_y = rssentries.cellLength_y;
            var cellLength_z = rssentries.cellLength_z;
            var cell = rssentries.cell;
            var orthMat = rssentries.fracToOrthMat;
            var fracMat = rssentries.orthToFracMat;
            console.log(nCells_x+" "+nCells_y+" "+nCells_z+" "+cellLength_x+" "+cellLength_y+" "+cellLength_z+" "+rssentries.grid.length);
            var start = new Date().getTime();
            var gridData = getEncodedData([rssentries.grid]);
            var end = new Date().getTime();
            var time = end - start;
            console.log('getEncodedData(mapGrid): ' + time*0.001+"s");
            console.log(gridData[0].length);
            var map = new CIsoSurface();
            var start = new Date().getTime();
            map.GenerateSurfacePartial(gridData[0],0.7,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,0,0,0,nCells_y/4,nCells_z/4,doChickenWire);
            //map.GenerateSurfacePartial(gridData[0],0.7,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x/4,nCells_y/4,nCells_z/4,3*nCells_x/4,3*nCells_y/4,3*nCells_z/4,doChickenWire);
            //map.GenerateSurface(gridData[0],0.7,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z);
            var end = new Date().getTime();
            var time = end - start;
            console.log('generate map: ' + time*0.001+"s");
            var start = new Date().getTime();
            var vertices = map.returnVertices(0,0,0,doChickenWire);
            console.log(vertices.length);
            var normals = map.returnNormals_new();
            console.log(normals.length);
            var indices = map.returnIndices();
            console.log(indices.length);
            var end = new Date().getTime();
            var time = end - start;
            console.log('copy buffers: ' + time*0.001+"s");

            colours = [];
            for(var i=0;i<vertices.length/3;i++){
                colours.push(0.0);
                colours.push(0.0);
                colours.push(1.0);
                colours.push(1.0);
            }
            //console.log(vertices);
            //console.log(normals);
            //console.log(indices);
            //console.log(colours);
            var mapTriangleData = {"col_tri":[[colours],[],[],[],[],[],[],[]], "norm_tri":[[normals],[],[],[],[],[],[],[]], "vert_tri":[[vertices],[],[],[],[],[],[],[]], "idx_tri":[[indices],[],[],[],[],[],[],[]] , "prim_types":[["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"]] };
            var theseBuffers = self.appendOtherData(mapTriangleData);
            var liveUpdatingMap = {};
            liveUpdatingMap["gridData"] = gridData;
            liveUpdatingMap["nCells_x"] = nCells_x;
            liveUpdatingMap["nCells_y"] = nCells_y;
            liveUpdatingMap["nCells_z"] = nCells_z;
            liveUpdatingMap["cellLength_x"] = cellLength_x;
            liveUpdatingMap["cellLength_y"] = cellLength_y;
            liveUpdatingMap["cellLength_z"] = cellLength_z;
            liveUpdatingMap["cell"] = cell;
            liveUpdatingMap["orthMat"] = orthMat;
            liveUpdatingMap["fracMat"] = fracMat;
            liveUpdatingMap["theseBuffers"] = theseBuffers;
            liveUpdatingMap["contourLevel"] = contourLevel;
            liveUpdatingMap["mapColour"] = mapColour;
            self.liveUpdatingMaps.push(liveUpdatingMap);
        }
    }

    var start = new Date().getTime();
    self.buildBuffers();
    var end = new Date().getTime();
    var time = end - start;
    console.log('buildBuffers: ' + time*0.001+"s");

    self.drawScene();

    console.log(divName);
    var fog = jsondata.fog;
    if(fog&&fog.length==2){
        console.log(fog);
        self.gl_fog_start = 500+fog[0];
        self.gl_fog_end = 500+fog[1];
        self.drawScene();
        self.div.dispatchEvent(self.fogChangedEvent);
    }

    var slab = jsondata.slab;
    if(slab&&slab.length==2){
        self.setSlab(slab);
    }

    var status = jsondata.status;
    if(status){
        //console.log(status);
        if (window.DOMParser) {
            var parser=new DOMParser();
            self.xmlDoc=parser.parseFromString(status,"text/xml");
        } else {
            self.xmlDoc=new ActiveXObject("Microsoft.XMLDOM");
            self.xmlDoc.async=false;
            self.xmlDoc.loadXML(status);
        }
        //console.log(self.xmlDoc);
        //console.log("read XML");
        var root = self.xmlDoc.getElementsByTagName("CCP4MG_Status")[0];
        self.backColour = "default";
        self.shinyBack = true;
        //console.log(root);
        if(root.getElementsByTagName("gui_params").length>0){
            var gui_params = root.getElementsByTagName("gui_params")[0];
            if(gui_params.getElementsByTagName("show_axes").length>0){
                var show_axes = root.getElementsByTagName("show_axes")[0];
                var show_axes_val = show_axes.childNodes[0].nodeValue;
                if(show_axes_val==="1"||show_axes_val==="true"||show_axes_val==="True"){
                    self.setShowAxes(true);
                } else {
                    self.setShowAxes(false);
                }
            }
        }
        if(root.getElementsByTagName("surface_drawing_style").length>0){
            var surface_drawing_style = root.getElementsByTagName("surface_drawing_style")[0];
            var inside_shiny = true;
            var inside_colour = "default";
            if(surface_drawing_style.getElementsByTagName("inside_colour").length>0){
                var this_inside_colour = surface_drawing_style.getElementsByTagName("inside_colour")[0].childNodes[0].nodeValue;
                if(this_inside_colour.length>0){
                    inside_colour = this_inside_colour;
                    if(inside_colour!="default"){
                        if(root.getElementsByTagName("Colours").length>0){
                            //console.log("!!!! Have Colours");
                            var Colours = root.getElementsByTagName("Colours")[0];
                            if(Colours.getElementsByTagName("colour_definitions").length>0){
                                var colour_definitions = Colours.getElementsByTagName("colour_definitions")[0];
                                //console.log("!!!! Have colour_definitions "+colour_definitions.length);
                                if(colour_definitions.getElementsByTagName("colour").length>0){
                                    var colours = colour_definitions.getElementsByTagName("colour");
                                    for(icol_defn=0;icol_defn<colours.length;icol_defn++){
                                        var theColour = colours[icol_defn];
                                        //console.log("!!!! Have colour "+icol_defn+" "+theColour);
                                        //console.log("!!!! Have colour (name) "+theColour.getAttribute("name"));
                                        if(theColour.getAttribute("name")===inside_colour){
                                            //console.log("Found colour! "+inside_colour);
                                            var red = theColour.getElementsByTagName("red")[0].childNodes[0].nodeValue;
                                            var green = theColour.getElementsByTagName("green")[0].childNodes[0].nodeValue;
                                            var blue = theColour.getElementsByTagName("blue")[0].childNodes[0].nodeValue;
                                            var alpha = theColour.getElementsByTagName("alpha")[0].childNodes[0].nodeValue;
                                            //console.log(red+" "+green+" "+blue+" "+alpha);
                                            self.backColour = [red,green,blue,alpha];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(surface_drawing_style.getElementsByTagName("inside_shiny").length>0){
                var this_inside_shiny = surface_drawing_style.getElementsByTagName("inside_shiny")[0].childNodes[0].nodeValue;
                if(this_inside_shiny=="0"||this_inside_shiny=="false"){
                    inside_shiny = false;
                }

            }
            self.shinyBack = inside_shiny;
            //console.log("inside_colour "+inside_colour);
            //console.log("inside_shiny "+inside_shiny);
        }
        var views = root.getElementsByTagName("View");
        //console.log(views);
        if(views.length>0){
            var view = views[0];
            if(view.getElementsByTagName("orientation").length>0){
                var orientation = view.getElementsByTagName("orientation")[0];
                console.log("Orientation");
                console.log(orientation);
                //try{
                var q0 = parseFloat(orientation.getElementsByTagName("q0")[0].childNodes[0].nodeValue);
                var q1 = parseFloat(orientation.getElementsByTagName("q1")[0].childNodes[0].nodeValue);
                var q2 = parseFloat(orientation.getElementsByTagName("q2")[0].childNodes[0].nodeValue);
                var q3 = parseFloat(orientation.getElementsByTagName("q3")[0].childNodes[0].nodeValue);
                self.myQuat = quat4.create([q1,q2,q3,-q0]);
                self.drawScene();
                //} catch (e) {
                //}
            }
            if(view.getElementsByTagName("orientation_auto").length>0){
                var orientation_auto = view.getElementsByTagName("orientation_auto")[0];
                console.log("Auto Orientation");
                console.log(orientation_auto);
            }
            if(view.getElementsByTagName("scale").length>0){
                var scale = view.getElementsByTagName("scale")[0];
                self.zoom = parseFloat(scale.childNodes[0].nodeValue)/120.;
                self.drawScene();
            }
            if(view.getElementsByTagName("scale_auto").length>0){
                var scale_auto = view.getElementsByTagName("scale_auto")[0];
                console.log("Auto Scale");
                console.log(scale_auto);
            }
        }
    }

    var ids = jsondata.ids;
    if(ids){
        self.ids = ids;
    }

    var texts = jsondata.texts;
    if(texts){
        self.textLabels = texts;
        self.drawScene();
    }

    self.displayOptions = jsondata.display_options;
    self.mapDisplayOptions = jsondata.map_display_options;

    self.div.dispatchEvent(self.statusChangedEvent);

    self.cifatoms = {};
    self.elements = {};
    self.restypes = {};
    self.models = {};
    self.altlocs = {};

    self.secstr = {};

    var ciftexts = jsondata.ciftexts;
    if(ciftexts){
        for(var cif in ciftexts){
            var cifatoms = parseMMCIF(ciftexts[cif]);
            self.cifatoms[cif] = cifatoms["atoms"];
            self.elements[cif] = cifatoms["elements"];
            self.restypes[cif] = cifatoms["restypes"];
            self.models[cif] = cifatoms["models"];
            self.altlocs[cif] = cifatoms["altlocs"];
        }
    }

    var secstr = jsondata.secstr;
    if(secstr){
        for(var thisSecStrId in secstr){
            var thisSecStr = secstr[thisSecStrId];
            self.secstr[thisSecStrId] = thisSecStr;
        }
        console.log(self.secstr);
    }

    var endAll = new Date().getTime();
    var timeAll = endAll - startAll;
    console.log('do all: ' + timeAll*0.001+"s");

}

function MGWebGL(divName,noLoadButtons,asyncShaders_in) {

    this.colours = {};

    this.colours["default"] =  [1.00,1.00,1.00,1.00];
    this.colours["ice blue"] =  [0.61,0.69,1.00,1.00];
    this.colours["gold"] =  [0.70,0.69,0.24,1.00];
    this.colours["coral"] =  [1.00,0.50,0.31,1.00];
    this.colours["grey"] =  [0.50,0.50,0.50,1.00];
    this.colours["pink"] =  [1.00,0.57,1.00,1.00];
    this.colours["sea green"] =  [0.50,0.73,0.71,1.00];
    this.colours["pale brown"] =  [0.66,0.49,0.37,1.00];
    this.colours["lilac"] =  [0.68,0.53,0.73,1.00];
    this.colours["lemon"] =  [1.00,1.00,0.50,1.00];
    this.colours["lawn green"] =  [0.27,0.61,0.31,1.00];
    this.colours["pale crimson"] =  [0.82,0.24,0.24,1.00];
    this.colours["light blue"] =  [0.25,0.60,0.88,1.00];
    this.colours["tan"] =  [0.47,0.00,0.00,1.00];
    this.colours["light green"] =  [0.60,1.00,0.60,1.00];
    this.colours["yellow"] =  [1.00,1.00,0.00,1.00];
    this.colours["white"] =  [1.00,1.00,1.00,1.00];
    this.colours["blue"] =  [0.00,0.00,1.00,1.00];
    this.colours["red"] =  [1.00,0.00,0.00,1.00];
    this.colours["green"] =  [0.00,1.00,0.00,1.00];
    this.colours["magenta"] =  [1.00,0.00,1.00,1.00];
    this.colours["cyan"] =  [0.00,1.00,0.88,1.00];
    this.colours["purple"] =  [0.58,0.00,1.00,1.00];
    this.colours["dark purple"] =  [0.57,0.13,0.34,1.00];
    this.colours["dark cyan"] =  [0.06,0.58,0.65,1.00];
    this.colours["black"] =  [0.00,0.00,0.00,1.00];
    this.colours["alice blue"] =  [0.94,0.97,1.00,1.00];
    this.colours["antique white"] =  [0.98,0.92,0.84,1.00];
    this.colours["aqua"] =  [0.00,1.00,1.00,1.00];
    this.colours["aquamarine"] =  [0.50,1.00,0.83,1.00];
    this.colours["azure"] =  [0.94,1.00,1.00,1.00];
    this.colours["beige"] =  [0.96,0.96,0.86,1.00];
    this.colours["bisque"] =  [1.00,0.89,0.77,1.00];
    this.colours["blanched almond"] =  [1.00,0.92,0.80,1.00];
    this.colours["blue violet"] =  [0.54,0.17,0.89,1.00];
    this.colours["brown"] =  [0.65,0.16,0.16,1.00];
    this.colours["burlywood"] =  [0.87,0.72,0.53,1.00];
    this.colours["cadet blue"] =  [0.37,0.62,0.63,1.00];
    this.colours["chartreuse"] =  [0.50,1.00,0.00,1.00];
    this.colours["chocolate"] =  [0.82,0.41,0.12,1.00];
    this.colours["cornflower blue"] =  [0.39,0.58,0.93,1.00];
    this.colours["cornsilk"] =  [1.00,0.97,0.86,1.00];
    this.colours["crimson"] =  [0.86,0.08,0.24,1.00];
    this.colours["dark blue"] =  [0.00,0.00,0.55,1.00];
    this.colours["dark goldenrod"] =  [0.72,0.53,0.04,1.00];
    this.colours["dark gray"] =  [0.66,0.66,0.66,1.00];
    this.colours["dark green"] =  [0.00,0.39,0.00,1.00];
    this.colours["dark grey"] =  [0.66,0.66,0.66,1.00];
    this.colours["dark khaki"] =  [0.74,0.72,0.42,1.00];
    this.colours["dark magenta"] =  [0.55,0.00,0.55,1.00];
    this.colours["dark olive green"] =  [0.33,0.42,0.18,1.00];
    this.colours["dark orange"] =  [1.00,0.55,0.00,1.00];
    this.colours["dark orchid"] =  [0.60,0.20,0.80,1.00];
    this.colours["dark red"] =  [0.55,0.00,0.00,1.00];
    this.colours["dark salmon"] =  [0.91,0.59,0.48,1.00];
    this.colours["dark sea green"] =  [0.56,0.74,0.56,1.00];
    this.colours["dark slate blue"] =  [0.28,0.24,0.55,1.00];
    this.colours["dark slate gray"] =  [0.18,0.31,0.31,1.00];
    this.colours["dark slate grey"] =  [0.18,0.31,0.31,1.00];
    this.colours["dark turquoise"] =  [0.00,0.81,0.82,1.00];
    this.colours["dark violet"] =  [0.58,0.00,0.83,1.00];
    this.colours["deep pink"] =  [1.00,0.08,0.58,1.00];
    this.colours["deep sky blue"] =  [0.00,0.75,1.00,1.00];
    this.colours["dim gray"] =  [0.41,0.41,0.41,1.00];
    this.colours["dim grey"] =  [0.41,0.41,0.41,1.00];
    this.colours["dodger blue"] =  [0.12,0.56,1.00,1.00];
    this.colours["firebrick"] =  [0.70,0.13,0.13,1.00];
    this.colours["floral white"] =  [1.00,0.98,0.94,1.00];
    this.colours["forest green"] =  [0.13,0.55,0.13,1.00];
    this.colours["fuchsia"] =  [1.00,0.00,1.00,1.00];
    this.colours["gainsboro"] =  [0.86,0.86,0.86,1.00];
    this.colours["ghost white"] =  [0.97,0.97,1.00,1.00];
    this.colours["goldenrod"] =  [0.85,0.65,0.13,1.00];
    this.colours["gray"] =  [0.75,0.75,0.75,1.00];
    this.colours["green yellow"] =  [0.68,1.00,0.18,1.00];
    this.colours["honeydew"] =  [0.94,1.00,0.94,1.00];
    this.colours["hot pink"] =  [1.00,0.41,0.71,1.00];
    this.colours["indian red"] =  [0.80,0.36,0.36,1.00];
    this.colours["indigo"] =  [0.29,0.00,0.51,1.00];
    this.colours["ivory"] =  [1.00,1.00,0.94,1.00];
    this.colours["khaki"] =  [0.94,0.90,0.55,1.00];
    this.colours["lavender blush"] =  [1.00,0.94,0.96,1.00];
    this.colours["lavender"] =  [0.90,0.90,0.98,1.00];
    this.colours["lemon chiffon"] =  [1.00,0.98,0.80,1.00];
    this.colours["light coral"] =  [0.94,0.50,0.50,1.00];
    this.colours["light cyan"] =  [0.88,1.00,1.00,1.00];
    this.colours["light goldenrod yellow"] =  [0.98,0.98,0.82,1.00];
    this.colours["light goldenrod"] =  [0.93,0.87,0.51,1.00];
    this.colours["light gray"] =  [0.83,0.83,0.83,1.00];
    this.colours["light grey"] =  [0.83,0.83,0.83,1.00];
    this.colours["light pink"] =  [1.00,0.71,0.76,1.00];
    this.colours["light salmon"] =  [1.00,0.63,0.48,1.00];
    this.colours["light sea green"] =  [0.13,0.70,0.67,1.00];
    this.colours["light sky blue"] =  [0.53,0.81,0.98,1.00];
    this.colours["light slate blue"] =  [0.52,0.44,1.00,1.00];
    this.colours["light slate gray"] =  [0.47,0.53,0.60,1.00];
    this.colours["light slate grey"] =  [0.47,0.53,0.60,1.00];
    this.colours["light steel blue"] =  [0.69,0.77,0.87,1.00];
    this.colours["light yellow"] =  [1.00,1.00,0.88,1.00];
    this.colours["lime green"] =  [0.20,0.80,0.20,1.00];
    this.colours["lime"] =  [0.00,1.00,0.00,1.00];
    this.colours["linen"] =  [0.98,0.94,0.90,1.00];
    this.colours["maroon"] =  [0.69,0.19,0.38,1.00];
    this.colours["medium aquamarine"] =  [0.40,0.80,0.67,1.00];
    this.colours["medium blue"] =  [0.00,0.00,0.80,1.00];
    this.colours["medium orchid"] =  [0.73,0.33,0.83,1.00];
    this.colours["medium purple"] =  [0.58,0.44,0.86,1.00];
    this.colours["medium sea green"] =  [0.24,0.70,0.44,1.00];
    this.colours["medium slate blue"] =  [0.48,0.41,0.93,1.00];
    this.colours["medium spring green"] =  [0.00,0.98,0.60,1.00];
    this.colours["medium turquoise"] =  [0.28,0.82,0.80,1.00];
    this.colours["medium violet red"] =  [0.78,0.08,0.52,1.00];
    this.colours["midnight blue"] =  [0.10,0.10,0.44,1.00];
    this.colours["mint cream"] =  [0.96,1.00,0.98,1.00];
    this.colours["misty rose"] =  [1.00,0.89,0.88,1.00];
    this.colours["moccasin"] =  [1.00,0.89,0.71,1.00];
    this.colours["navajo white"] =  [1.00,0.87,0.68,1.00];
    this.colours["navy blue"] =  [0.00,0.00,0.50,1.00];
    this.colours["navy"] =  [0.00,0.00,0.50,1.00];
    this.colours["old lace"] =  [0.99,0.96,0.90,1.00];
    this.colours["olive drab"] =  [0.42,0.56,0.14,1.00];
    this.colours["olive"] =  [0.50,0.50,0.00,1.00];
    this.colours["orange red"] =  [1.00,0.27,0.00,1.00];
    this.colours["orange"] =  [1.00,0.65,0.00,1.00];
    this.colours["orchid"] =  [0.85,0.44,0.84,1.00];
    this.colours["pale goldenrod"] =  [0.93,0.91,0.67,1.00];
    this.colours["pale green"] =  [0.60,0.98,0.60,1.00];
    this.colours["pale turquoise"] =  [0.69,0.93,0.93,1.00];
    this.colours["pale violet red"] =  [0.86,0.44,0.58,1.00];
    this.colours["papaya whip"] =  [1.00,0.94,0.84,1.00];
    this.colours["peach puff"] =  [1.00,0.85,0.73,1.00];
    this.colours["peru"] =  [0.80,0.52,0.25,1.00];
    this.colours["plum"] =  [0.87,0.63,0.87,1.00];
    this.colours["powder blue"] =  [0.69,0.88,0.90,1.00];
    this.colours["rebecca purple"] =  [0.40,0.20,0.60,1.00];
    this.colours["rosy brown"] =  [0.74,0.56,0.56,1.00];
    this.colours["royal blue"] =  [0.25,0.41,0.88,1.00];
    this.colours["saddle brown"] =  [0.55,0.27,0.07,1.00];
    this.colours["salmon"] =  [0.98,0.50,0.45,1.00];
    this.colours["sandy brown"] =  [0.96,0.64,0.38,1.00];
    this.colours["seashell"] =  [1.00,0.96,0.93,1.00];
    this.colours["sienna"] =  [0.63,0.32,0.18,1.00];
    this.colours["silver"] =  [0.75,0.75,0.75,1.00];
    this.colours["sky blue"] =  [0.53,0.81,0.92,1.00];
    this.colours["slate blue"] =  [0.42,0.35,0.80,1.00];
    this.colours["slate gray"] =  [0.44,0.50,0.56,1.00];
    this.colours["slate grey"] =  [0.44,0.50,0.56,1.00];
    this.colours["snow"] =  [1.00,0.98,0.98,1.00];
    this.colours["spring green"] =  [0.00,1.00,0.50,1.00];
    this.colours["steel blue"] =  [0.27,0.51,0.71,1.00];
    this.colours["teal"] =  [0.00,0.50,0.50,1.00];
    this.colours["thistle"] =  [0.85,0.75,0.85,1.00];
    this.colours["tomato"] =  [1.00,0.39,0.28,1.00];
    this.colours["turquoise"] =  [0.25,0.88,0.82,1.00];
    this.colours["violet red"] =  [0.82,0.13,0.56,1.00];
    this.colours["violet"] =  [0.93,0.51,0.93,1.00];
    this.colours["web gray"] =  [0.50,0.50,0.50,1.00];
    this.colours["web green"] =  [0.00,0.50,0.00,1.00];
    this.colours["web grey"] =  [0.50,0.50,0.50,1.00];
    this.colours["web maroon"] =  [0.50,0.00,0.00,1.00];
    this.colours["web purple"] =  [0.50,0.00,0.50,1.00];
    this.colours["wheat"] =  [0.96,0.87,0.70,1.00];
    this.colours["white smoke"] =  [0.96,0.96,0.96,1.00];
    this.colours["x11 gray"] =  [0.75,0.75,0.75,1.00];
    this.colours["x11 green"] =  [0.00,1.00,0.00,1.00];
    this.colours["x11 grey"] =  [0.75,0.75,0.75,1.00];
    this.colours["x11 maroon"] =  [0.69,0.19,0.38,1.00];
    this.colours["x11 purple"] =  [0.63,0.13,0.94,1.00];
    this.colours["yellow green"] =  [0.60,0.80,0.20,1.00];

    this.symcols = [this.colours["ice_blue"], this.colours["gold"], this.colours["coral"], this.colours["grey"], this.colours["pink"], this.colours["sea_green"], this.colours["pale_brown"], this.colours["lilac"], this.colours["lemon"], this.colours["lawn_green"], this.colours["pale_crimson"], this.colours["light_blue"], this.colours["tan"], this.colours["light_green"], this.colours["red"], this.colours["green"], this.colours["blue"], this.colours["magenta"]];

    var asyncShaders = true;
    if(typeof(asyncShaders_in)!="undefined"){
        asyncShaders = asyncShaders_in;
    } 

    this.shinyBack = true;
    this.backColour = "default";

    this.ready = false;
    this.canvas = null;
    this.gl = null;
    this.background_colour = [0,0,0,1];
    this.textTex = null;
    this.origin = [0.0,0.0,0.0];
    this.radius = 60.0;
    this.init_x = null;
    this.init_y = null;
    this.dx = null;
    this.dy = null;
    this.myQuat = null;
    this.mouseDown = null;
    this.mouseMoved = null;
    this.zoom = null;
    this.ext = null;
    this.gl_fog_start = null;
    this.gl_fog_end = null;
    this.gl_nClipPlanes = null;
    this.shaderProgram = null;
    this.shaderProgramTextBackground = null;
    this.shaderProgramLines = null;
    this.shaderProgramPointSpheres = null;

    this.mvMatrix = mat4.create();
    this.mvInvMatrix = mat4.create();
    this.screenZ = vec3.create();
    this.pMatrix = mat4.create();

    this.gl_clipPlane0 = null;
    this.gl_clipPlane1 = null;
    this.gl_clipPlane2 = null;
    this.gl_clipPlane3 = null;
    this.gl_clipPlane4 = null;
    this.gl_clipPlane5 = null;
    this.gl_clipPlane6 = null;
    this.gl_clipPlane7 = null;

    this.clickedAtoms = [];
    this.textCtx = null;

    this.displayBuffers = [];
    this.liveUpdatingMaps = [];
    this.currentBufferIdx = -1;
    this.shapesBuffers = {};
    this.shapesVertices = {};

    this.xmlDoc = null;

    this.save_pixel_data = false;

    this.showAxes = false;

    this.doShadow = false;
    this.doShadowDepthDebug = false;

    var self = this;

    this.div = document.getElementById(divName);

    this.table = document.createElement("table");
    var tr = document.createElement("tr");
    var lefttd = document.createElement("td");
    var righttd = document.createElement("td");
    this.rightdiv = document.createElement("div");
    this.leftdiv = document.createElement("div");
    this.table.appendChild(tr);
    var tr2 = document.createElement("tr");
    var bottd = document.createElement("td");
    tr2.appendChild(bottd);

    this.table.appendChild(tr2);

    var bottomHeight;
    if(!noLoadButtons) {
      bottomHeight = 60;
    } else {
      bottomHeight = 0;
    }
    this.leftdiv.width = parseInt(this.div.style.width);
    this.leftdiv.height = parseInt(this.div.style.height)-bottomHeight;
    this.leftdiv.style.maxHeight = parseInt(this.div.style.height)-bottomHeight;
    tr.style.maxHeight = this.div.style.height;
    this.leftdiv.style.maxWidth = parseInt(this.div.style.width)+"px";
    tr.appendChild(lefttd);
    tr.appendChild(righttd);
    righttd.appendChild(this.rightdiv);
    lefttd.appendChild(this.leftdiv);
    this.rightdiv.style.maxHeight = parseInt(this.div.style.height)-bottomHeight;
    this.rightdiv.style.display = "inline-block";//parseInt(this.div.style.height)-bottomHeight;
    this.rightdiv.style.height = parseInt(this.div.style.height)-bottomHeight+"px";
    this.rightdiv.style.width = "0px";
    this.rightdiv.setAttribute("id","rightdiv");
    this.div.appendChild(this.table);

    this.rightdiv.style.overflow = "scroll";
    righttd.style.verticalAlign = "top";
    this.rightdiv.style.verticalAlign = "top";
    lefttd.style.verticalAlign = "top";
    this.leftdiv.style.verticalAlign = "top";
    tr.style.verticalAlign = "top";
    //console.log("Max height "+lefttd.style.maxHeight);
    //console.log("Max height "+this.rightdiv.style.maxHeight);

    this.canvas = document.createElement("canvas");

    this.canvas.style.width = parseInt(this.div.style.width)+"px";
    this.canvas.style.height = parseInt(this.div.style.height)-bottomHeight+"px";
    this.canvas.width = getDeviceScale()*parseInt(this.div.style.width);
    this.canvas.height = getDeviceScale()*(parseInt(this.div.style.height)-bottomHeight);

    this.leftdiv.appendChild(this.canvas);

    // These are stubs for "Simple UI"
    this.infoNode = document.createElement("div");
    this.rightdiv.appendChild(this.infoNode);

    this.botdiv = document.createElement("div");
    bottd.appendChild(this.botdiv);

    if(!noLoadButtons) {
      var fileInput = document.createElement("input");
      fileInput.setAttribute("type","file");
      fileInput.setAttribute("id","upfile"+divName);
      fileInput.setAttribute("name","upfile"+divName);
      this.botdiv.appendChild(fileInput);
      var uploadButton = document.createElement("input");
      uploadButton.setAttribute("type","button");
      uploadButton.setAttribute("value","Load file");
      uploadButton.onclick = function() {
        upload(fileInput,self);
      }
      this.botdiv.appendChild(uploadButton);

      var breakElement = document.createElement("br");
      this.botdiv.appendChild(breakElement);
      var pdbInput = document.createElement("input");
      pdbInput.setAttribute("type","text");
      pdbInput.setAttribute("id","pdbcode"+divName);
      pdbInput.setAttribute("name","pdbcode"+divName);
      this.botdiv.appendChild(pdbInput);
      var loadPdbCode = document.createElement("input");
      loadPdbCode.setAttribute("type","button");
      loadPdbCode.setAttribute("value","Load PDB code");
      loadPdbCode.onclick = function() {
        uploadpdb(pdbInput,self);
      }
      this.botdiv.appendChild(loadPdbCode);

      pdbInput.addEventListener("keypress", function(event) {
        if (event.keyCode == 13) uploadpdb(pdbInput,self);
      });

    }

    //rightcolumn.style.float = "right";
    //leftcolumn.style.float = "left";
     
    //this.canvas = document.getElementById(canvasName);

    this.textCtx = document.createElement("canvas").getContext("2d");

    this.myQuat = quat4.create([0,0,0,-1]);
    //var testmat = quat4.toMat4(this.myQuat)
    //alert(mat4.str(testmat))
    this.zoom = 1.0;

    this.gl_clipPlane0 = new glMatrixArrayType(4);
    this.gl_clipPlane1 = new glMatrixArrayType(4);
    this.gl_clipPlane1[0] = 0.0;
    this.gl_clipPlane1[1] = 0.0;
    this.gl_clipPlane1[2] = 1.0;
    this.gl_clipPlane1[3] = 1000.0;
    this.gl_clipPlane0[0] = 0.0;
    this.gl_clipPlane0[1] = 0.0;
    this.gl_clipPlane0[2] = -1.0;
    this.gl_clipPlane0[3] = -0.0;
    this.gl_clipPlane2 = new glMatrixArrayType(4);
    this.gl_clipPlane3 = new glMatrixArrayType(4);
    this.gl_clipPlane4 = new glMatrixArrayType(4);
    this.gl_clipPlane5 = new glMatrixArrayType(4);
    this.gl_clipPlane6 = new glMatrixArrayType(4);
    this.gl_clipPlane7 = new glMatrixArrayType(4);
    this.clickedAtoms = [];
    this.textLabels = [];
    this.displayOptions = [];
    this.ids = [];
    this.cifatoms = {};
    this.elements = {};
    this.restypes = {};
    this.models = {};
    this.altlocs = {};
    this.secstr = {};

    this.gl_cursorPos = new glMatrixArrayType(2);
    this.gl_cursorPos[0] = this.canvas.width/2.;
    this.gl_cursorPos[1] = this.canvas.height/2.;

    this.gl = initGL(this.canvas);
    var extensionArray = this.gl.getSupportedExtensions();
    console.log(extensionArray);

    //this.stuartTexture = initTextures(this.gl);

    //console.log("The GL context ... ");
    //console.log(this.gl);
    this.ext = this.gl.getExtension("OES_element_index_uint");
    if(!this.ext) {
       alert("No OES_element_index_uint support");
     }
    this.deriv_ext = this.gl.getExtension("OES_standard_derivatives");
    if(!this.deriv_ext) {
       alert("No GL_OES_standard_derivatives support");
    }
    this.frag_depth_ext = this.gl.getExtension("EXT_frag_depth");

    this.depth_texture = this.gl.getExtension("WEBGL_depth_texture");
    if(!this.depth_texture){
        this.depth_texture = this.gl.getExtension("MOZ_WEBGL_depth_texture");
        if(!this.depth_texture){
            this.depth_texture = this.gl.getExtension("WEBKIT_WEBGL_depth_texture");
            if(!this.depth_texture){
                alert("No depth texture extension");
            }
        }
    }

    this.initTextureFramebuffer();

    if(!this.frag_depth_ext) {
       console.log("No EXT_frag_depth support");
       console.log("This is supported in most browsers, except IE. And may never be supported in IE.");
       console.log("This extension is supported in Microsoft Edge, so Windows 10 is required for perfect spheres in MS Browser.");
       console.log("Other browers on Windows 7/8/8.1 do have this extension.");
    }
    this.textTex = this.gl.createTexture();
    this.gl_nClipPlanes = 0;
    this.gl_fog_start = 500.0;
    this.gl_fog_end = 1500.0;
    //this.gl.lineWidth(2.0);
    this.gl.blendFuncSeparate(this.gl.SRC_ALPHA, this.gl.ONE_MINUS_SRC_ALPHA, this.gl.ONE, this.gl.ONE_MINUS_SRC_ALPHA)
    this.gl.enable(this.gl.BLEND)

    /*
    this.statusChangedEvent = new CustomEvent('statusChanged');
    this.viewChangedEvent = new Event('viewChanged');
    this.fogChangedEvent = new Event('fogChanged');
    this.clipChangedEvent = new CustomEvent('clipChanged');
    */

    // Sigh doing the old-fashioned way for IE.
    this.statusChangedEvent = document.createEvent('Event');
    this.statusChangedEvent.initEvent('statusChanged',true,true);

    this.viewChangedEvent = document.createEvent('Event');
    this.viewChangedEvent.initEvent('viewChanged',true,true);

    this.fogChangedEvent = document.createEvent('Event');
    this.fogChangedEvent.initEvent('fogChanged',true,true);

    this.clipChangedEvent = document.createEvent('Event');
    this.clipChangedEvent.initEvent('clipChanged',true,true);

    var vertexShader;
    var fragmentShader;
    var lineVertexShader;
    var thickLineVertexShader;
    var lineFragmentShader;
    var textVertexShader;
    var textFragmentShader;
    var pointSpheresVertexShader;
    var pointSpheresFragmentShader;
    var twoDShapesFragmentShader;
    var twoDShapesVertexShader;
    var renderFrameBufferVertexShader;
    var renderFrameBufferFragmentShader;
    var perfectSphereFragmentShader;
    var shadowVertexShader;
    var shadowFragmentShader;
    var triangleShadowVertexShader;
    var triangleShadowFragmentShader;
    var twoDShapesShadowVertexShader;
    var perfectSphereShadowFragmentShader;
    var pointSpheresShadowVertexShader;
    var pointSpheresShadowFragmentShader;

    this.mygetrequest=new ajaxRequest();
    
    
    
    self.mygetrequest.onreadystatechange=function(){
        if (self.mygetrequest.readyState==4){
            if (self.mygetrequest.status==200 || window.location.href.indexOf("http")==-1){
                var jsondata=JSON.parse(self.mygetrequest.responseText);
                self.loadJSON(jsondata);
            } else {
                alert("An error has occured making the request")
            }
        }
    }

    var self = this;

    this.doRedraw = false;
    var myVar = setInterval(function(){ self.drawSceneIfDirty() }, 16);

    if(asyncShaders){

        var pointsphereshadowvertreq=new ajaxRequest();
        pointsphereshadowvertreq.onreadystatechange=function(){
            if (pointsphereshadowvertreq.readyState==4){
                if (pointsphereshadowvertreq.status==200 || window.location.href.indexOf("http")==-1){
                    if(self.frag_depth_ext){
                        pointSpheresShadowVertexShader = getShader(self.gl,pointsphereshadowvertreq.responseText,"vertex");
                    }
        var pointsphereshadowfragreq=new ajaxRequest();
        pointsphereshadowfragreq.onreadystatechange=function(){
            if (pointsphereshadowfragreq.readyState==4){
                if (pointsphereshadowfragreq.status==200 || window.location.href.indexOf("http")==-1){
                    if(self.frag_depth_ext){
                        pointSpheresShadowFragmentShader = getShader(self.gl,pointsphereshadowfragreq.responseText,"fragment");
                    }

        var shapesshadowvertreq=new ajaxRequest();
        shapesshadowvertreq.onreadystatechange=function(){
            if (shapesshadowvertreq.readyState==4){
                if (shapesshadowvertreq.status==200 || window.location.href.indexOf("http")==-1){
                    if(self.frag_depth_ext){
                        twoDShapesShadowVertexShader = getShader(self.gl,shapesshadowvertreq.responseText,"vertex");
                    }
        var perfectsphereshadowfragreq=new ajaxRequest();
        perfectsphereshadowfragreq.onreadystatechange=function(){
            if (perfectsphereshadowfragreq.readyState==4){
                if (perfectsphereshadowfragreq.status==200 || window.location.href.indexOf("http")==-1){
                    if(self.frag_depth_ext){
                        perfectSphereShadowFragmentShader = getShader(self.gl,perfectsphereshadowfragreq.responseText,"fragment");
                    }

        var triangleshadowvertreq=new ajaxRequest();
        triangleshadowvertreq.onreadystatechange=function(){
            if (triangleshadowvertreq.readyState==4){
                if (triangleshadowvertreq.status==200 || window.location.href.indexOf("http")==-1){
                    if(self.frag_depth_ext){
                        triangleShadowVertexShader = getShader(self.gl,triangleshadowvertreq.responseText,"vertex");
                    }
        var triangleshadowfragreq=new ajaxRequest();
        triangleshadowfragreq.onreadystatechange=function(){
            if (triangleshadowfragreq.readyState==4){
                if (triangleshadowfragreq.status==200 || window.location.href.indexOf("http")==-1){
                    if(self.frag_depth_ext){
                        triangleShadowFragmentShader = getShader(self.gl,triangleshadowfragreq.responseText,"fragment");
                    }
        var shadowvertreq=new ajaxRequest();
        shadowvertreq.onreadystatechange=function(){
            if (shadowvertreq.readyState==4){
                if (shadowvertreq.status==200 || window.location.href.indexOf("http")==-1){
                    if(self.frag_depth_ext){
                        shadowVertexShader = getShader(self.gl,shadowvertreq.responseText,"vertex");
                    }
        var shadowfragreq=new ajaxRequest();
        shadowfragreq.onreadystatechange=function(){
            if (shadowfragreq.readyState==4){
                if (shadowfragreq.status==200 || window.location.href.indexOf("http")==-1){
                    if(self.frag_depth_ext){
                        shadowFragmentShader = getShader(self.gl,shadowfragreq.responseText,"fragment");
                    }
        var perfectfragreq=new ajaxRequest();
        perfectfragreq.onreadystatechange=function(){
            if (perfectfragreq.readyState==4){
                if (perfectfragreq.status==200 || window.location.href.indexOf("http")==-1){
                    if(self.frag_depth_ext){
                        perfectSphereFragmentShader = getShader(self.gl,perfectfragreq.responseText,"fragment");
                    }
                    var shapesvertreq=new ajaxRequest();
                    shapesvertreq.onreadystatechange=function(){
                        if (shapesvertreq.readyState==4){
                            if (shapesvertreq.status==200 || window.location.href.indexOf("http")==-1){
                                twoDShapesVertexShader = getShader(self.gl,shapesvertreq.responseText,"vertex");
                                var shapesfragreq=new ajaxRequest();
                                shapesfragreq.onreadystatechange=function(){
                                    if (shapesfragreq.readyState==4){
                                        if (shapesfragreq.status==200 || window.location.href.indexOf("http")==-1){
                                            twoDShapesFragmentShader = getShader(self.gl,shapesfragreq.responseText,"fragment");
                                            var trivertreq=new ajaxRequest();
                                            trivertreq.onreadystatechange=function(){
                                                if (trivertreq.readyState==4){
                                                    if (trivertreq.status==200 || window.location.href.indexOf("http")==-1){
                                                        vertexShader = getShader(self.gl,trivertreq.responseText,"vertex");
                                                        var trifragreq=new ajaxRequest();
                                                        trifragreq.onreadystatechange=function(){
                                                            if (trifragreq.readyState==4){
                                                                if (trifragreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                    fragmentShader = getShader(self.gl,trifragreq.responseText,"fragment");
                                                                    var thickLinevertreq=new ajaxRequest();
                                                                    thickLinevertreq.onreadystatechange=function(){
                                                                        if (thickLinevertreq.readyState==4){
                                                                            if (thickLinevertreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                                thickLineVertexShader = getShader(self.gl,thickLinevertreq.responseText,"vertex");
                                                                                var linevertreq=new ajaxRequest();
                                                                                linevertreq.onreadystatechange=function(){
                                                                                    if (linevertreq.readyState==4){
                                                                                        if (linevertreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                                            lineVertexShader = getShader(self.gl,linevertreq.responseText,"vertex");
                                                                                            var linefragreq=new ajaxRequest();
                                                                                            linefragreq.onreadystatechange=function(){
                                                                                                if (linefragreq.readyState==4){
                                                                                                    if (linefragreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                                                        lineFragmentShader = getShader(self.gl,linefragreq.responseText,"fragment");
                                                                                                        var textvertreq=new ajaxRequest();
                                                                                                        textvertreq.onreadystatechange=function(){
                                                                                                            if (textvertreq.readyState==4){
                                                                                                                if (textvertreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                                                                    textVertexShader = getShader(self.gl,textvertreq.responseText,"vertex");
                                                                                                                    var textfragreq=new ajaxRequest();
                                                                                                                    textfragreq.onreadystatechange=function(){
                                                                                                                        if (textfragreq.readyState==4){
                                                                                                                            if (textfragreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                                                                                textFragmentShader = getShader(self.gl,textfragreq.responseText,"fragment");
                                                                                                                                var pointspherevertreq=new ajaxRequest();
                                                                                                                                pointspherevertreq.onreadystatechange=function(){
                                                                                                                                    if (pointspherevertreq.readyState==4){
                                                                                                                                        if (pointspherevertreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                                                                                            pointSpheresVertexShader = getShader(self.gl,pointspherevertreq.responseText,"vertex");
                                                                                                                                            var pointspherefragreq=new ajaxRequest();
                                                                                                                                            pointspherefragreq.onreadystatechange=function(){
                                                                                                                                                if (pointspherefragreq.readyState==4){
                                                                                                                                                    if (pointspherefragreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                                                                                                        pointSpheresFragmentShader = getShader(self.gl,pointspherefragreq.responseText,"fragment");
                                                                                                                                                        var renderfbfragreq=new ajaxRequest();
                                                                                                                                                        renderfbfragreq.onreadystatechange=function(){
                                                                                                                                                            if (renderfbfragreq.readyState==4){
                                                                                                                                                                if (renderfbfragreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                                                                                                                    renderFrameBufferFragmentShader = getShader(self.gl,renderfbfragreq.responseText,"fragment");
                                                                                                                                                                    var renderfbvertreq=new ajaxRequest();
                                                                                                                                                                    renderfbvertreq.onreadystatechange=function(){
                                                                                                                                                                        if (renderfbvertreq.readyState==4){
                                                                                                                                                                            if (renderfbvertreq.status==200 || window.location.href.indexOf("http")==-1){
                                                                                                                                                                                renderFrameBufferVertexShader = getShader(self.gl,renderfbvertreq.responseText,"vertex");

                                                                                                                                                                                self.initRenderFrameBufferShaders(renderFrameBufferVertexShader,renderFrameBufferFragmentShader);
                                                                                                                                                                                self.initLineShaders(lineVertexShader,lineFragmentShader);
                                                                                                                                                                                self.initThickLineShaders(thickLineVertexShader,lineFragmentShader);
                                                                                                                                                                                self.initPointSpheresShaders(pointSpheresVertexShader,pointSpheresFragmentShader);
                                                                                                                                                                                self.initTwoDShapesShaders(twoDShapesVertexShader,twoDShapesFragmentShader);
                                                                                                                                                                                self.initImageShaders(twoDShapesVertexShader,textFragmentShader);
                                                                                                                                                                                if(self.frag_depth_ext){
                                                                                                                                                                                    self.initPerfectSphereShaders(twoDShapesVertexShader,perfectSphereFragmentShader);
                                                                                                                                                                                    self.initShadowShaders(shadowVertexShader,shadowFragmentShader);
                                                                                                                                                                                    self.initTriangleShadowShaders(triangleShadowVertexShader,triangleShadowFragmentShader);
                                                                                                                                                                                    self.initPerfectSphereShadowShaders(twoDShapesShadowVertexShader,perfectSphereShadowFragmentShader);
                                                                                                                                                                                    self.initPointSpheresShadowShaders(pointSpheresShadowVertexShader,pointSpheresShadowFragmentShader);
                                                                                                                                                                                }
                                                                                                                                                                                self.initTextBackgroundShaders(textVertexShader,textFragmentShader);
                                                                                                                                                                                self.gl.disableVertexAttribArray(self.shaderProgramTextBackground.vertexTextureAttribute);
                                                                                                                                                                                self.initShaders(vertexShader,fragmentShader);

                                                                                                                                                                                self.buildBuffers();

                                                                                                                                                                                self.gl.clearColor(self.background_colour[0],self.background_colour[1],self.background_colour[2],self.background_colour[3]);
                                                                                                                                                                                self.gl.enable(self.gl.DEPTH_TEST);

                                                                                                                                                                                self.origin = [0.0, 0.0, 0.0];
                                                                                                                                                                                shader_version = self.gl.getParameter(self.gl.SHADING_LANGUAGE_VERSION);
                                                                                                                                                                                //self.infoNode.innerHTML=shader_version;

                                                                                                                                                                                self.mouseDown = false;
                                                                                                                                                                                self.canvas.addEventListener("mousedown",
                                                                                                                                                                                        function(evt){
                                                                                                                                                                                        self.doMouseDown(evt,self);
                                                                                                                                                                                        },
                                                                                                                                                                                        false);
                                                                                                                                                                                self.canvas.addEventListener("mouseup",
                                                                                                                                                                                        function(evt){
                                                                                                                                                                                        self.doMouseUp(evt,self);
                                                                                                                                                                                        },
                                                                                                                                                                                        false);
                                                                                                                                                                                self.canvas.addEventListener("mousemove",
                                                                                                                                                                                        function(evt){
                                                                                                                                                                                        self.doMouseMove(evt,self);
                                                                                                                                                                                        },
                                                                                                                                                                                        false);
                                                                                                                                                                                self.canvas.addEventListener("click",
                                                                                                                                                                                        function(evt){
                                                                                                                                                                                        self.doClick(evt,self);
                                                                                                                                                                                        },
                                                                                                                                                                                        false);

                                                                                                                                                                                self.canvas.addEventListener("mouseenter",
                                                                                                                                                                                        function(evt){
                                                                                                                                                                                        document.onkeydown = function(evt2){
                                                                                                                                                                                        self.handleKeyDown(evt2,self);
                                                                                                                                                                                        }
                                                                                                                                                                                        },
                                                                                                                                                                                        false);
                                                                                                                                                                                self.canvas.addEventListener("mouseleave",
                                                                                                                                                                                        function(evt){
                                                                                                                                                                                        document.onkeydown = function(evt2){
                                                                                                                                                                                        }
                                                                                                                                                                                        },
                                                                                                                                                                                        false);
                                                                                                                                                                                self.canvas.addEventListener("wheel",
                                                                                                                                                                                        function(evt){
                                                                                                                                                                                        self.doWheel(evt);
                                                                                                                                                                                        },
                                                                                                                                                                                        false);
                                                                                                                                                                                self.gl.viewportWidth = self.canvas.width;
                                                                                                                                                                                self.gl.viewportHeight = self.canvas.height;
                                                                                                                                                                                self.light_positions = new Float32Array([0.0,0.0,60.0,1.0]);
                                                                                                                                                                                self.light_colours_ambient = new Float32Array([0.0,0.0,0.0,1.0]);
                                                                                                                                                                                self.light_colours_specular = new Float32Array([1.0,1.0,1.0,1.0]);
                                                                                                                                                                                self.light_colours_diffuse = new Float32Array([1.0,1.0,1.0,1.0]);
                                                                                                                                                                                self.drawScene();
                                                                                                                                                                                self.ready = true;
                                                                                                                                                                            }
                                                                                                                                                                        }
                                                                                                                                                                    }
                                                                                                                                                                    renderfbvertreq.open("GET", "triangle-vertex-shader.glsl", true);
                                                                                                                                                                    renderfbvertreq.send(null);
                                                                                                                                                                }
                                                                                                                                                            }
                                                                                                                                                        }
                                                                                                                                                        renderfbfragreq.open("GET", "render-framebuffer-fragment-shader.glsl", true);
                                                                                                                                                        renderfbfragreq.send(null);

                                                                                                                                                    }
                                                                                                                                                }
                                                                                                                                            }
                                                                                                                                            pointspherefragreq.open("GET", "pointspheres-fragment-shader.glsl", true);
                                                                                                                                            pointspherefragreq.send(null);
                                                                                                                                        }
                                                                                                                                    }
                                                                                                                                }
                                                                                                                                pointspherevertreq.open("GET", "pointspheres-vertex-shader.glsl", true);
                                                                                                                                pointspherevertreq.send(null);
                                                                                                                            }
                                                                                                                        }
                                                                                                                    }
                                                                                                                    textfragreq.open("GET", "text-fragment-shader.glsl", true);
                                                                                                                    textfragreq.send(null);
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                        textvertreq.open("GET", "triangle-vertex-shader.glsl", true);
                                                                                                        textvertreq.send(null);
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                            linefragreq.open("GET", "lines-fragment-shader.glsl", true);
                                                                                            linefragreq.send(null);
                                                                                        }
                                                                                    }
                                                                                }
                                                                                linevertreq.open("GET", "lines-vertex-shader.glsl", true);
                                                                                linevertreq.send(null);
                                                                            }
                                                                        }
                                                                    }
                                                                    thickLinevertreq.open("GET", "thick-lines-vertex-shader.glsl", true);
                                                                    thickLinevertreq.send(null);
                                                                }
                                                            }
                                                        }
                                                        trifragreq.open("GET", "triangle-fragment-shader.glsl", true);
                                                        trifragreq.send(null);
                                                    }
                                                }
                                            }
                                            trivertreq.open("GET", "triangle-vertex-shader.glsl", true);
                                            trivertreq.send(null);
                                        }
                                    }
                                }
                                shapesfragreq.open("GET", "twodshapes-fragment-shader.glsl", true);
                                shapesfragreq.send(null);
                            }
                        }
                    }
                    shapesvertreq.open("GET", "twodshapes-vertex-shader.glsl", true);
                    shapesvertreq.send(null);
                }
            }
        }
        perfectfragreq.open("GET", "perfect-sphere-fragment-shader.glsl", true);
        perfectfragreq.send(null);
                }
            }
        }
        shadowfragreq.open("GET", "shadow-fragment-shader.glsl", true);
        shadowfragreq.send(null);
                }
            }
        }
        shadowvertreq.open("GET", "shadow-vertex-shader.glsl", true);
        shadowvertreq.send(null);
                }
            }
        }
        triangleshadowfragreq.open("GET", "triangle-shadow-fragment-shader.glsl", true);
        triangleshadowfragreq.send(null);
                }
            }
        }
        triangleshadowvertreq.open("GET", "triangle-shadow-vertex-shader.glsl", true);
        triangleshadowvertreq.send(null);
                }
            }
        }
        perfectsphereshadowfragreq.open("GET", "perfect-sphere-shadow-fragment-shader.glsl", true);
        perfectsphereshadowfragreq.send(null);
                }
            }
        }
        shapesshadowvertreq.open("GET", "twodshapes-shadow-vertex-shader.glsl", true);
        shapesshadowvertreq.send(null);

                }
            }
        }
        pointsphereshadowfragreq.open("GET", "pointspheres-shadow-fragment-shader.glsl", true);
        pointsphereshadowfragreq.send(null);
                }
            }
        }
        pointsphereshadowvertreq.open("GET", "pointspheres-shadow-vertex-shader.glsl", true);
        pointsphereshadowvertreq.send(null);
        
    } else {
        vertexShader                    = getShader(self.gl,triangle_vertex_shader_source,"vertex");
        fragmentShader                  = getShader(self.gl,triangle_fragment_shader_source,"fragment");
        lineVertexShader                = getShader(self.gl,lines_vertex_shader_source,"vertex");
        thickLineVertexShader           = getShader(self.gl,thick_lines_vertex_shader_source,"vertex");
        lineFragmentShader              = getShader(self.gl,lines_fragment_shader_source,"fragment");
        textVertexShader                = getShader(self.gl,triangle_vertex_shader_source,"vertex");
        textFragmentShader              = getShader(self.gl,text_fragment_shader_source,"fragment");
        pointSpheresVertexShader        = getShader(self.gl,pointspheres_vertex_shader_source,"vertex");
        pointSpheresFragmentShader      = getShader(self.gl,pointspheres_fragment_shader_source,"fragment");
        twoDShapesVertexShader          = getShader(self.gl,twod_vertex_shader_source,"vertex");
        twoDShapesFragmentShader        = getShader(self.gl,twod_fragment_shader_source,"fragment");
        renderFrameBufferVertexShader   = getShader(self.gl,triangle_vertex_shader_source,"vertex");
        renderFrameBufferFragmentShader = getShader(self.gl,render_framebuffer_fragment_shader_source,"fragment");
        if(self.frag_depth_ext){
            perfectSphereFragmentShader     = getShader(self.gl,perfect_sphere_fragment_shader_source,"fragment");
            shadowVertexShader                    = getShader(self.gl,shadow_vertex_shader_source,"vertex");
            shadowFragmentShader                  = getShader(self.gl,shadow_fragment_shader_source,"fragment");
            triangleShadowVertexShader                    = getShader(self.gl,triangle_shadow_vertex_shader_source,"vertex");
            triangleShadowFragmentShader                  = getShader(self.gl,triangle_shadow_fragment_shader_source,"fragment");
            pointSpheresShadowVertexShader        = getShader(self.gl,pointspheres_shadow_vertex_shader_source,"vertex");
            pointSpheresShadowFragmentShader      = getShader(self.gl,pointspheres_shadow_fragment_shader_source,"fragment");
            self.initShadowShaders(shadowVertexShader,shadowFragmentShader);
            self.initTriangleShadowShaders(triangleShadowVertexShader,triangleShadowFragmentShader);
            twoDShapesShadowVertexShader          = getShader(self.gl,twod_vertex_shadow_shader_source,"vertex");
            perfectSphereShadowFragmentShader     = getShader(self.gl,perfect_sphere_shadow_fragment_shader_source,"fragment");
            self.initPointSpheresShadowShaders(pointSpheresShadowVertexShader,pointSpheresShadowFragmentShader);
        }

        self.initRenderFrameBufferShaders(renderFrameBufferVertexShader,renderFrameBufferFragmentShader);
        self.initLineShaders(lineVertexShader,lineFragmentShader);
        self.initThickLineShaders(thickLineVertexShader,lineFragmentShader);
        self.initPointSpheresShaders(pointSpheresVertexShader,pointSpheresFragmentShader);
        self.initTwoDShapesShaders(twoDShapesVertexShader,twoDShapesFragmentShader);
        self.initImageShaders(twoDShapesVertexShader,textFragmentShader);
        if(self.frag_depth_ext){
            self.initPerfectSphereShaders(twoDShapesVertexShader,perfectSphereFragmentShader);
            self.initPerfectSphereShadowShaders(twoDShapesShadowVertexShader,perfectSphereShadowFragmentShader);
        }
        self.initTextBackgroundShaders(textVertexShader,textFragmentShader);
        self.gl.disableVertexAttribArray(self.shaderProgramTextBackground.vertexTextureAttribute);
        self.initShaders(vertexShader,fragmentShader);

        self.buildBuffers();

        self.gl.clearColor(self.background_colour[0],self.background_colour[1],self.background_colour[2],self.background_colour[3]);
        self.gl.enable(self.gl.DEPTH_TEST);

        self.origin = [0.0, 0.0, 0.0];
        shader_version = self.gl.getParameter(self.gl.SHADING_LANGUAGE_VERSION);
        //self.infoNode.innerHTML=shader_version;

        self.mouseDown = false;
        self.mouseDownButton = -1;

        self.canvas.addEventListener("mousedown",
            function(evt){
            self.doMouseDown(evt,self);
            },
            false);
        self.canvas.addEventListener("mouseup",
            function(evt){
            self.doMouseUp(evt,self);
            },
            false);
        self.canvas.addEventListener("mousemove",
            function(evt){
            self.doMouseMove(evt,self);
            },
            false);
        self.canvas.addEventListener("click",
            function(evt){
            self.doClick(evt,self);
            },
            false);
        self.canvas.addEventListener("mouseenter",
            function(evt){
            document.onkeydown = function(evt2){
            self.handleKeyDown(evt2,self);
            }
            },
            false);
        self.canvas.addEventListener("mouseleave",
            function(evt){
            document.onkeydown = function(evt2){
            }
            },
            false);
        self.canvas.addEventListener("wheel",
            function(evt){
            self.doWheel(evt);
            evt.preventDefault();
            },
            false);
        self.gl.viewportWidth = self.canvas.width;
        self.gl.viewportHeight = self.canvas.height;
        self.light_positions = new Float32Array([0.0,0.0,60.0,1.0]);
        self.light_colours_ambient = new Float32Array([0.0,0.0,0.0,1.0]);
        self.light_colours_specular = new Float32Array([1.0,1.0,1.0,1.0]);
        self.light_colours_diffuse = new Float32Array([1.0,1.0,1.0,1.0]);
        self.drawScene();
        self.ready = true;
    }

}

MGWebGL.prototype.setOrientationRightAndUp = function (right,up) {
    var rvec = vec3.create(right);
    vec3.normalize(rvec);
    var upvec = vec3.create(up);
    vec3.normalize(upvec);
    var zvec = vec3.create();
    vec3.cross(rvec,upvec,zvec);
    var m = [];
    for(var i=0;i<16;i++){
        m[i] = 0.0;
    }
    m[0] = rvec[0];
    m[4] = rvec[1];
    m[8] = rvec[2];
    m[1] = upvec[0];
    m[5] = upvec[1];
    m[9] = upvec[2];
    m[2] = zvec[0];
    m[6] = zvec[1];
    m[10] = zvec[2];
    m[15] = 1.0;

    this.setOrientationMatrix(m);
}

MGWebGL.prototype.setOrientationMatrix3x3 = function (m_in) {
    var m = [];
    for(var i=0;i<16;i++){
        m[i] = 0.0;
    }
    m[0] = m_in[0];
    m[4] = m_in[3];
    m[8] = m_in[6];
    m[1] = m_in[1];
    m[5] = m_in[4];
    m[9] = m_in[7];
    m[2] = m_in[2];
    m[6] = m_in[5];
    m[10] = m_in[8];
    m[15] = 1.0;

    this.setOrientationMatrix(m);
}

MGWebGL.prototype.setOrientationMatrix = function (m) {
    var trace = m[0] + m[5] + m[10] + m[15];

    var dval = [];
    dval.push(0.0);
    dval.push(0.0);
    dval.push(0.0);
    dval.push(0.0);

    if( trace > 1e-7 ) {
        var s = 0.5 / Math.sqrt(trace);
        dval[0] = 0.25 / s;
        dval[1] = ( m[2*4+1]- m[1*4+2]) * s;
        dval[2] = ( m[0*4+2]- m[2*4+0]) * s;
        dval[3] = ( m[1*4+0]- m[0*4+1]) * s;
    } else {
        if ( m[0*4+0]> m[1*4+1]&& m[0*4+0]> m[2*4+2]) {
            var s = 2.0 * Math.sqrt( 1.0 + m[0*4+0]- m[1*4+1]- m[2*4+2]);
            dval[1] = 0.25 * s;
            dval[2] = (m[0*4+1]+ m[1*4+0]) / s;
            dval[3] = (m[0*4+2]+ m[2*4+0]) / s;
            dval[0] = (m[1*4+2]- m[2*4+1]) / s;

        } else if (m[1*4+1]> m[2*4+2]) {
            var s = 2.0 * Math.sqrt( 1.0 + m[1*4+1]- m[0*4+0]- m[2*4+2]);
            dval[1] = (m[0*4+1]+ m[1*4+0]) / s;
            dval[2] = 0.25 * s;
            dval[3] = (m[1*4+2]+ m[2*4+1]) / s;
            dval[0] = (m[0*4+2]- m[2*4+0]) / s;
        } else {
            var s = 2.0 * Math.sqrt( 1.0 + m[2*4+2]- m[0*4+0]- m[1*4+1]);
            dval[1] = (m[0*4+2]+ m[2*4+0]) / s;
            dval[2] = (m[1*4+2]+ m[2*4+1]) / s;
            dval[3] = 0.25 * s;
            dval[0] = (m[0*4+1]- m[1*4+0]) / s;
        }
    }
    // FIXME: My C++ quat has different order/sign from glMatrix version. So I cheat and fudge here.
    //        I should rewrite above.
    var dval_new = [];
    dval_new[3] = -dval[0];
    dval_new[0] = dval[1];
    dval_new[1] = dval[2];
    dval_new[2] = dval[3];

    var q = quat4.create(dval_new);
    var invQuat = quat4.create();
    var theMatrix = quat4.toMat4(q);
    quat4.inverse(q,invQuat);
    
    this.myQuat = invQuat;
    var x = vec3.create([1,0,0]);
    var y = vec3.create([0,1,0]);
    var z = vec3.create([0,0,1]);
    mat4.multiplyVec3(theMatrix,x);
    mat4.multiplyVec3(theMatrix,y);
    mat4.multiplyVec3(theMatrix,z);
    //console.log(x);
    //console.log(y);
    //console.log(z);
    this.drawScene();
}

MGWebGL.prototype.setQuat = function (q) {
    this.myQuat = q;
    this.drawScene();
}

MGWebGL.prototype.setBackground = function (col) {
    this.background_colour = col;
    this.drawScene();
}

MGWebGL.prototype.setOrigin = function (o) {
    this.origin = o;
    this.drawScene();
}

MGWebGL.prototype.setAmbientLight = function (r,g,b) {
    this.light_colours_ambient = new Float32Array([r,g,b,1.0]);
    this.drawScene();
}

MGWebGL.prototype.setSpecularLight = function (r,g,b) {
    this.light_colours_specular = new Float32Array([r,g,b,1.0]);
    this.drawScene();
}

MGWebGL.prototype.setDiffuseLight = function (r,g,b) {
    this.light_colours_diffuse = new Float32Array([r,g,b,1.0]);
    this.drawScene();
}

MGWebGL.prototype.setLightPosition = function (x,y,z) {
    this.light_positions = new Float32Array([x,y,z,1.0]);
    this.drawScene();
}

MGWebGL.prototype.setZoom = function (z) {
    this.zoom = z;
    this.drawScene();
}

MGWebGL.prototype.setShowAxes = function (a) {
    this.showAxes = a;
    this.drawScene();
}

MGWebGL.prototype.initTextureFramebuffer = function () {

    this.rttFramebuffer = this.gl.createFramebuffer();
    this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, this.rttFramebuffer);

    this.rttFramebuffer.width = 1024;
    this.rttFramebuffer.height = 1024;

    this.rttTexture = this.gl.createTexture();
    this.gl.bindTexture(this.gl.TEXTURE_2D, this.rttTexture);
    this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_MAG_FILTER, this.gl.LINEAR);
    this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_MIN_FILTER, this.gl.LINEAR_MIPMAP_NEAREST);
    
    this.gl.texImage2D(this.gl.TEXTURE_2D, 0, this.gl.RGBA, this.rttFramebuffer.width, this.rttFramebuffer.height, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, null);
    this.gl.generateMipmap(this.gl.TEXTURE_2D);

    var renderbuffer = this.gl.createRenderbuffer();
    this.gl.bindRenderbuffer(this.gl.RENDERBUFFER, renderbuffer);
    this.gl.renderbufferStorage(this.gl.RENDERBUFFER, this.gl.DEPTH_COMPONENT16, this.rttFramebuffer.width, this.rttFramebuffer.height);
    this.gl.framebufferTexture2D(this.gl.FRAMEBUFFER, this.gl.COLOR_ATTACHMENT0, this.gl.TEXTURE_2D, this.rttTexture, 0);
    this.gl.framebufferRenderbuffer(this.gl.FRAMEBUFFER, this.gl.DEPTH_ATTACHMENT, this.gl.RENDERBUFFER, renderbuffer);
    this.gl.bindTexture(this.gl.TEXTURE_2D, null);
    this.gl.bindRenderbuffer(this.gl.RENDERBUFFER, null);
    this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);

    /*
    this.rttFramebuffer2 = this.gl.createFramebuffer();
    this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, this.rttFramebuffer2);

    this.rttFramebuffer2.width = 1024;
    this.rttFramebuffer2.height = 1024;

    this.rttTexture2 = this.gl.createTexture();
    this.gl.bindTexture(this.gl.TEXTURE_2D, this.rttTexture2);
    this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_MAG_FILTER, this.gl.LINEAR);
    this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_MIN_FILTER, this.gl.LINEAR_MIPMAP_NEAREST);
    
    this.gl.texImage2D(this.gl.TEXTURE_2D, 0, this.gl.RGBA, this.rttFramebuffer2.width, this.rttFramebuffer2.height, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, null);
    this.gl.generateMipmap(this.gl.TEXTURE_2D);

    var renderbuffer2 = this.gl.createRenderbuffer();
    this.gl.bindRenderbuffer(this.gl.RENDERBUFFER, renderbuffer2);
    this.gl.renderbufferStorage(this.gl.RENDERBUFFER, this.gl.DEPTH_COMPONENT16, this.rttFramebuffer2.width, this.rttFramebuffer2.height);
    this.gl.framebufferTexture2D(this.gl.FRAMEBUFFER, this.gl.COLOR_ATTACHMENT0, this.gl.TEXTURE_2D, this.rttTexture2, 0);
    this.gl.framebufferRenderbuffer(this.gl.FRAMEBUFFER, this.gl.DEPTH_ATTACHMENT, this.gl.RENDERBUFFER, renderbuffer2);
    this.gl.bindTexture(this.gl.TEXTURE_2D, null);
    this.gl.bindRenderbuffer(this.gl.RENDERBUFFER, null);
    this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
    */


    this.rttFramebufferDepth = null;
    if(this.depth_texture){
        this.rttFramebufferDepth = this.gl.createFramebuffer();
        this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, this.rttFramebufferDepth);
        this.rttFramebufferDepth.width = 1024;
        this.rttFramebufferDepth.height = 1024;
        this.rttTextureDepth = this.gl.createTexture();
        this.gl.bindTexture(this.gl.TEXTURE_2D, this.rttTextureDepth);
        this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_MAG_FILTER, this.gl.NEAREST);
        this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_MIN_FILTER, this.gl.NEAREST);
        //this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_COMPARE_MODE, this.gl.COMPARE_R_TO_TEXTURE);
        //this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_COMPARE_FUNC, this.gl.LEQUAL);
        this.gl.texImage2D(this.gl.TEXTURE_2D, 0, this.gl.DEPTH_COMPONENT, this.rttFramebufferDepth.width, this.rttFramebufferDepth.height, 0, this.gl.DEPTH_COMPONENT, this.gl.UNSIGNED_SHORT, null);
        var renderbufferCol = this.gl.createRenderbuffer();
        this.gl.bindRenderbuffer(this.gl.RENDERBUFFER, renderbufferCol);
        this.gl.renderbufferStorage(this.gl.RENDERBUFFER, this.gl.RGBA4, this.rttFramebufferDepth.width, this.rttFramebufferDepth.height);
        this.gl.framebufferTexture2D(this.gl.FRAMEBUFFER, this.gl.DEPTH_ATTACHMENT, this.gl.TEXTURE_2D, this.rttTextureDepth, 0);
        this.gl.framebufferRenderbuffer(this.gl.FRAMEBUFFER, this.gl.COLOR_ATTACHMENT0, this.gl.RENDERBUFFER, renderbufferCol);
        this.gl.bindTexture(this.gl.TEXTURE_2D, null);
        this.gl.bindRenderbuffer(this.gl.RENDERBUFFER, null);
        this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
    }

}

MGWebGL.prototype.getNameForId = function (id,type) {
    if(this.xmlDoc){
        var root = this.xmlDoc.getElementsByTagName("CCP4MG_Status")[0];
        var data = root.getElementsByTagName(type);
        for(var idat=0;idat<data.length;idat++){
            var dataUUID = getNodeText(data[idat].getElementsByTagName("uuid")[0]);
            if(dataUUID==id){
                console.log("Match");
                if(data[idat].getElementsByTagName("name").length>0){
                    var theName = getNodeText(data[idat].getElementsByTagName("name")[0]);
                    if(theName.length>0){
                        return theName;
                    } else {
                        return "unknown";
                    }
                } else {
                    return "unknown";
                }

            }
        }
    }
    return "unknown";
}

MGWebGL.prototype.addDisplayObject = function (id) {
    console.log("addDisplayObject to "+id);
    var self = this;
    if(this.xmlDoc){
        var root = this.xmlDoc.getElementsByTagName("CCP4MG_Status")[0];
        var data = root.getElementsByTagName("MolData");
        for(var idat=0;idat<data.length;idat++){
            var dataUUID = getNodeText(data[idat].getElementsByTagName("uuid")[0]);
            if(dataUUID==id){
                console.log("Match");
                console.log(data[idat]);
                var newDoc = (new DOMParser()).parseFromString('<dummy/>', 'text/xml');
                newDoc.removeChild(newDoc.documentElement);
                newRoot = root.cloneNode(true);
                newDoc.appendChild(newRoot);
                thisDataCopy = data[idat].cloneNode(true);
                while(newRoot.getElementsByTagName("MolData").length>0){
                    newRoot.removeChild(newRoot.getElementsByTagName("MolData")[0]);
                }
                while(newRoot.getElementsByTagName("MapData").length>0){
                    newRoot.removeChild(newRoot.getElementsByTagName("MapData")[0]);
                }
                while(newRoot.getElementsByTagName("Annotation").length>0){
                    newRoot.removeChild(newRoot.getElementsByTagName("Annotation")[0]);
                }
                while(newRoot.getElementsByTagName("Legend").length>0){
                    newRoot.removeChild(newRoot.getElementsByTagName("Legend")[0]);
                }
                while(newRoot.getElementsByTagName("Image").length>0){
                    newRoot.removeChild(newRoot.getElementsByTagName("Image")[0]);
                }
                while(newRoot.getElementsByTagName("SVGImage").length>0){
                    newRoot.removeChild(newRoot.getElementsByTagName("SVGImage")[0]);
                }
                while(thisDataCopy.getElementsByTagName("MolDisp").length>0){
                    thisDataCopy.removeChild(thisDataCopy.getElementsByTagName("MolDisp")[0]);
                }
                while(thisDataCopy.getElementsByTagName("HBonds").length>0){
                    thisDataCopy.removeChild(thisDataCopy.getElementsByTagName("HBonds")[0]);
                }
                while(thisDataCopy.getElementsByTagName("Contacts").length>0){
                    thisDataCopy.removeChild(thisDataCopy.getElementsByTagName("Contacts")[0]);
                }
                while(thisDataCopy.getElementsByTagName("SurfaceDispobj").length>0){
                    thisDataCopy.removeChild(thisDataCopy.getElementsByTagName("SurfaceDispobj")[0]);
                }

                var newDisp = newDoc.createElement("MolDisp");

                var newSelParam = newDoc.createElement("selection_parameters");
                var newSelect = newDoc.createElement("select");
                var newCid = newDoc.createElement("cid");
                newSelect.textContent = "cid";
                newCid.textContent = "all";
                newSelParam.appendChild(newSelect);
                newSelParam.appendChild(newCid);
                newDisp.appendChild(newSelParam);

                var newStyleMode = newDoc.createElement("style");
                newStyleMode.textContent = "BONDS";
                newDisp.appendChild(newStyleMode);

                var newColParam = newDoc.createElement("colour_parameters");
                var newColMode = newDoc.createElement("colour_mode");
                newColMode.textContent = "atomtype";
                newColParam.appendChild(newColMode);
                newDisp.appendChild(newColParam);

                thisDataCopy.appendChild(newDisp);
                newRoot.appendChild(thisDataCopy);

                var snippet = (new XMLSerializer()).serializeToString(newDoc);
                console.log(snippet);
                var blob  = new Blob([snippet], {type : "text/plain"});
                var formData = new FormData();
                if (typeof(middlewareToken) !== 'undefined'){
                    formData.append('csrfmiddlewaretoken',middlewareToken);
                }
                var theFile = guid()+".mgpic.xml";
                var client = new XMLHttpRequest();
                formData.append("upfile",blob,theFile);
                client.onreadystatechange = function() {
                    if (client.readyState == 4 && client.status == 200) {
                        console.log("Sent snippet");
                        var client2 = new XMLHttpRequest();
                        client2.onreadystatechange = function() {
                            if (client2.readyState == 4 && client2.status == 200) {
                                console.log("Executed snippet");
                                var jsondata=JSON.parse(client2.responseText);
                                var idx = 0;
                                for(var idat2=0;idat2<data.length;idat2++){
                                    var thisDataId = getNodeText(data[idat2].getElementsByTagName("uuid")[0]);
                                    // FIXME - These are the only four we are currently dealing with. Others may come later.
                                    idx += data[idat2].getElementsByTagName("MolDisp").length;
                                    idx += data[idat2].getElementsByTagName("HBonds").length;
                                    idx += data[idat2].getElementsByTagName("Contacts").length;
                                    idx += data[idat2].getElementsByTagName("MolDSurfaceDispobjisp").length;
                                    console.log(thisDataId+" "+id);
                                    if(thisDataId==id){
                                        break;
                                    }
                                }

                                //console.log("--------------------------------------------------");

                                //console.log("Insert at "+idx);
                                //console.log("buffers: "+self.displayBuffers.length);
                                //console.log("options: "+self.displayOptions.length);
                                //console.log("ids: "+self.ids.length);
                                //console.log("textLabels: "+self.textLabels.length);

                                self.displayBuffers.splice(idx,0,new DisplayBuffer());
                                self.displayOptions.splice(idx,0,jsondata.display_options[0]);
                                var newId = jsondata.ids[0];
                                //console.log("The new id is "+newId);

                                var underIndex = newId.indexOf("_");
                                var newDispId = newId.substring(underIndex+1);

                                //console.log(newId);
                                self.ids.splice(idx,0,id+"_"+newDispId);
                                self.textLabels.splice(idx,0,jsondata.texts[0]);

                                //FIXME - Need to clone newDisp and add clone to this.xmlDoc's data[idat]. Does this work?
                                var thisNewDisp = newDisp.cloneNode(true);
                                data[idat].appendChild(thisNewDisp);

                                var newThisUUIDElement = self.xmlDoc.createElement("uuid");
                                newThisUUIDElement.textContent = newDispId;
                                thisNewDisp.appendChild(newThisUUIDElement);

                                //console.log("The new DispObj");
                                //console.log(thisNewDisp);

                                //console.log("The new ids");
                                //console.log(self.ids);
                                //console.log("--------------------------------------------------");

                                var newUUIDElement = newDoc.createElement("uuid"); // Probably not necessary
                                newUUIDElement.textContent = newDispId.substr; // Probably not necessary
                                newDisp.appendChild(newUUIDElement); // Probably not necessary

                                self.currentBufferIdx = idx;

                                var rssentries=jsondata.norm_tri[0];
                                var norms = getEncodedData(rssentries);
                                var thisName = jsondata.names[0];
                                self.displayBuffers[self.currentBufferIdx].name_label = thisName;
                                for (var i=0; i<norms.length; i++){
                                    self.createNormalBuffer(norms[i]);
                                }
                                var rssentries=jsondata.vert_tri[0];
                                var tris = getEncodedData(rssentries);
                                for (var i=0; i<tris.length; i++){
                                    self.createVertexBuffer(tris[i]);
                                }
                                var rssentries=jsondata.col_tri[0];
                                var colours = getEncodedData(rssentries);
                                for (var i=0; i<colours.length; i++){
                                    self.createColourBuffer(colours[i]);
                                }
                                var rssentries=jsondata.idx_tri[0];
                                var idxs = getEncodedData(rssentries);
                                for (var i=0; i<idxs.length; i++){
                                    self.createIndexBuffer(idxs[i]);
                                }
                                var rssentries=jsondata.prim_types[0];
                                for (var i=0; i<rssentries.length; i++){
                                    self.displayBuffers[self.currentBufferIdx].bufferTypes.push(rssentries[i]);
                                }
                                var rssentries=getEncodedData(jsondata.sizes[0]);
                                if(typeof(rssentries)!=="undefined"){
                                    console.log(rssentries);
                                    for (var i=0; i<rssentries.length; i++){
                                        self.createSizeBuffer(rssentries[i]);
                                    }
                                }
                                var atoms = jsondata.atoms[0];
                                self.displayBuffers[self.currentBufferIdx].atoms = atoms;
                                self.displayBuffers[self.currentBufferIdx].isDirty = true;
                                self.buildBuffers();
                                self.drawScene();
                                // FIXME - This does not set correct menu entries.
                                self.div.dispatchEvent(self.statusChangedEvent);
                            }
                        }
                        client2.open("GET", "get_triangles?pdb="+theFile, true);
                        client2.send(null);
                    }
                }
                client.open("POST", "/server_client.html", true);
                client.send(formData);
                break;
            }
        }
    }
}

MGWebGL.prototype.centreOn = function (idx) {
    var self = this;
    if(self.displayBuffers[idx].atoms.length>0){
        var xtot = 0;
        var ytot = 0;
        var ztot = 0;
        for (var j = 0; j < self.displayBuffers[idx].atoms.length; j++){
            xtot += self.displayBuffers[idx].atoms[j].x;
            ytot += self.displayBuffers[idx].atoms[j].y;
            ztot += self.displayBuffers[idx].atoms[j].z;
        }
        xtot /= self.displayBuffers[idx].atoms.length;
        ytot /= self.displayBuffers[idx].atoms.length;
        ztot /= self.displayBuffers[idx].atoms.length;

        var new_origin = [-xtot,-ytot,-ztot];
        var old_origin = [self.origin[0],self.origin[1],self.origin[2]];

        var myVar = setInterval(function(){ myTimer() }, 10);
        var frac = 0;
        function myTimer() {
            var ffrac = 0.01 * frac;
            self.origin = [ffrac*new_origin[0]+(1.0-ffrac)*old_origin[0],ffrac*new_origin[1]+(1.0-ffrac)*old_origin[1],ffrac*new_origin[2]+(1.0-ffrac)*old_origin[2]];
            self.drawScene();
            if(frac>99){
                clearInterval(myVar);
            }
            frac += 1;
        }
        /*
        self.origin = [-xtot,-ytot,-ztot];
        self.drawScene();
        */

    }
}

MGWebGL.prototype.removeDisplayObject = function (num) {
    var self = this;
    if(this.xmlDoc){
        var root = this.xmlDoc.getElementsByTagName("CCP4MG_Status")[0];
        var data = root.getElementsByTagName("MolData");
        for(var idat=0;idat<data.length;idat++){
            var dataUUID = getNodeText(data[idat].getElementsByTagName("uuid")[0]);
            var dispobjs = data[idat].getElementsByTagName("MolDisp");
            var hbobjs = data[idat].getElementsByTagName("HBonds");
            var contactobjs = data[idat].getElementsByTagName("Contacts");
            var surfobs = data[idat].getElementsByTagName("SurfaceDispobj");
            //console.log(idat+": "+dataUUID+" "+this.ids[num]);
            //console.log("All ids");
            //console.log(this.ids);
            // FIXME - Once again we are only dealing with DispObjs....
            for(var iobj=0;iobj<dispobjs.length;iobj++){
                var theUUID = dataUUID+"_"+getNodeText(dispobjs[iobj].getElementsByTagName("uuid")[0]);
                //console.log("Test DispObj uuid: "+getNodeText(dispobjs[iobj].getElementsByTagName("uuid")[0]));
                if(theUUID && theUUID === this.ids[num]){
                   console.log("Match a delete ......."+theUUID);
                   data[idat].removeChild(dispobjs[iobj]);
                   console.log(self.displayBuffers);
                   self.displayBuffers.splice(num,1);
                   self.displayOptions.splice(num,1);
                   self.ids.splice(num,1);
                   self.textLabels.splice(num,1);
                   self.buildBuffers();
                   self.drawScene();
                   // FIXME - This does not set correct menu entries.
                   self.div.dispatchEvent(self.statusChangedEvent);
                   return;
                }
            }
        }
    }
}

MGWebGL.prototype.updateDisplayObject = function (num,property,value) {
    //console.log("!!!!!!!!!!!!!!!!!!!! MGWebGL.prototype.updateDisplayObject");
    console.log(num,property,value);
    var self = this;
    function addToColours(colour,doc){
        var root = doc.getElementsByTagName("CCP4MG_Status")[0];
        if(root.getElementsByTagName("Colours").length==0){
            var newColours = doc.createElement("Colours");
            root.appendChild(newColours);
        }
        if(root.getElementsByTagName("Colours")[0].getElementsByTagName("colour_definitions").length==0){
            var newColoursDef = doc.createElement("colour_definitions");
            root.getElementsByTagName("Colours")[0].appendChild(newColoursDef);
        }
        var haveThisColour = false;
        var colours = root.getElementsByTagName("Colours")[0].getElementsByTagName("colour_definitions")[0].getElementsByTagName("colour");
        for(var icol=0;icol<colours.length;icol++){
            if(colours[icol].getAttribute("name")===colour){
                haveThisColour = true;
                break;
            }
        }
        if(!haveThisColour){
            var newCol = doc.createElement("colour");
            newCol.setAttribute("name",colour);
            var redHex = colour.substring(1,3);
            var greenHex = colour.substring(3,5);
            var blueHex = colour.substring(5,7);
            var red = parseInt(redHex,16)/255.;
            var green = parseInt(greenHex,16)/255.;
            var blue = parseInt(blueHex,16)/255.;
            var alpha = 1.0;
            var r = doc.createElement("red");
            r.textContent = ""+red;
            newCol.appendChild(r);
            var g = doc.createElement("green");
            g.textContent = ""+green;
            newCol.appendChild(g);
            var b = doc.createElement("blue");
            b.textContent = ""+blue;
            newCol.appendChild(b);
            var a = doc.createElement("alpha");
            a.textContent = ""+blue;
            newCol.appendChild(a);
            root.getElementsByTagName("Colours")[0].getElementsByTagName("colour_definitions")[0].appendChild(newCol);
        }
    }
    if(this.xmlDoc){
        var root = this.xmlDoc.getElementsByTagName("CCP4MG_Status")[0];
        var data = root.getElementsByTagName("MolData");
        for(var idat=0;idat<data.length;idat++){
            var dataUUID = getNodeText(data[idat].getElementsByTagName("uuid")[0]);
            var dispobjs = data[idat].getElementsByTagName("MolDisp");
            var hbobjs = data[idat].getElementsByTagName("HBonds");
            var contactobjs = data[idat].getElementsByTagName("Contacts");
            var surfobs = data[idat].getElementsByTagName("SurfaceDispobj");
            //console.log(idat+": "+dataUUID+" "+this.ids[num]);
            //console.log("All ids");
            //console.log(this.ids);
            for(var iobj=0;iobj<dispobjs.length;iobj++){
                var theUUID = dataUUID+"_"+getNodeText(dispobjs[iobj].getElementsByTagName("uuid")[0]);
                //console.log("Test DispObj uuid: "+getNodeText(dispobjs[iobj].getElementsByTagName("uuid")[0]));
                if(theUUID && theUUID === this.ids[num]){
                    //console.log("Match a change .......");
                    var newDoc = (new DOMParser()).parseFromString('<dummy/>', 'text/xml');
                    newDoc.removeChild(newDoc.documentElement);
                    newRoot = root.cloneNode(true);
                    newDoc.appendChild(newRoot);
                    thisDataCopy = data[idat].cloneNode(true);
                    thisDispCopy = dispobjs[iobj].cloneNode(true);
                    while(newRoot.getElementsByTagName("MolData").length>0){
                        newRoot.removeChild(newRoot.getElementsByTagName("MolData")[0]);
                    }
                    while(newRoot.getElementsByTagName("MapData").length>0){
                        newRoot.removeChild(newRoot.getElementsByTagName("MapData")[0]);
                    }
                    while(newRoot.getElementsByTagName("Annotation").length>0){
                        newRoot.removeChild(newRoot.getElementsByTagName("Annotation")[0]);
                    }
                    while(newRoot.getElementsByTagName("Legend").length>0){
                        newRoot.removeChild(newRoot.getElementsByTagName("Legend")[0]);
                    }
                    while(newRoot.getElementsByTagName("Image").length>0){
                        newRoot.removeChild(newRoot.getElementsByTagName("Image")[0]);
                    }
                    while(newRoot.getElementsByTagName("SVGImage").length>0){
                        newRoot.removeChild(newRoot.getElementsByTagName("SVGImage")[0]);
                    }
                    while(thisDataCopy.getElementsByTagName("MolDisp").length>0){
                        thisDataCopy.removeChild(thisDataCopy.getElementsByTagName("MolDisp")[0]);
                    }
                    while(thisDataCopy.getElementsByTagName("HBonds").length>0){
                        thisDataCopy.removeChild(thisDataCopy.getElementsByTagName("HBonds")[0]);
                    }
                    while(thisDataCopy.getElementsByTagName("Contacts").length>0){
                        thisDataCopy.removeChild(thisDataCopy.getElementsByTagName("Contacts")[0]);
                    }
                    while(thisDataCopy.getElementsByTagName("SurfaceDispobj").length>0){
                        thisDataCopy.removeChild(thisDataCopy.getElementsByTagName("SurfaceDispobj")[0]);
                    }
                    if(property === "Colour"){
                        if(thisDispCopy.getElementsByTagName("colour_parameters").length>0){
                            var colparam = thisDispCopy.getElementsByTagName("colour_parameters")[0];
                            if(colparam.getElementsByTagName("colour_mode").length>0){
                                if(value.substring(0,11)==="one_colour_"){
                                    colparam.getElementsByTagName("colour_mode")[0].textContent = "one_colour";
                                    dispobjs[iobj].getElementsByTagName("colour_parameters")[0].getElementsByTagName("colour_mode")[0].textContent = "one_colour";
                                    if(colparam.getElementsByTagName("one_colour").length==0){
                                        var newOneCol = newDoc.createElement("one_colour");
                                        var newOneColThis = this.xmlDoc.createElement("one_colour");
                                        colparam.appendChild(newOneCol);
                                        dispobjs[iobj].getElementsByTagName("colour_parameters")[0].appendChild(newOneColThis);
                                    }
                                    colparam.getElementsByTagName("one_colour")[0].textContent = "#"+value.substring(11);
                                    dispobjs[iobj].getElementsByTagName("colour_parameters")[0].getElementsByTagName("one_colour")[0].textContent = "#"+value.substring(11);
                                    addToColours("#"+value.substring(11),newDoc);
                                    addToColours("#"+value.substring(11),this.xmlDoc);
                                } else {
                                    colparam.getElementsByTagName("colour_mode")[0].textContent = value;
                                    dispobjs[iobj].getElementsByTagName("colour_parameters")[0].getElementsByTagName("colour_mode")[0].textContent = value;
                                }
                            } else {
                                var newColMode = newDoc.createElement("colour_mode");
                                colparam.appendChild(newColMode);
                                var newColModeThis = this.xmlDoc.createElement("colour_mode");
                                if(value.substring(0,11)==="one_colour_"){
                                    newColMode.textContent = "one_colour";
                                    newColModeThis.textContent = "one_colour";
                                    var newOneCol = newDoc.createElement("one_colour");
                                    var newOneColThis = this.xmlDoc.createElement("one_colour");
                                    colparam.appendChild(newOneCol);
                                    dispobjs[iobj].getElementsByTagName("colour_parameters")[0].appendChild(newOneColThis);
                                    newOneCol.textContent = "#"+value.substring(11);
                                    newOneColThis.textContent = "#"+value.substring(11);
                                    addToColours("#"+value.substring(11),newDoc);
                                    addToColours("#"+value.substring(11),this.xmlDoc);
                                } else {
                                    newColMode.textContent = value;
                                    newColModeThis.textContent = value;
                                }
                                dispobjs[iobj].getElementsByTagName("colour_parameters")[0].appendChild(newColModeThis);
                            }
                        } else {
                            var newColParam = newDoc.createElement("colour_parameters");
                            var newColMode = newDoc.createElement("colour_mode");
                            newColParam.appendChild(newColMode);
                            thisDispCopy.appendChild(newColParam);
                            var newColParamThis = this.xmlDoc.createElement("colour_parameters");
                            var newColModeThis = this.xmlDoc.createElement("colour_mode");
                            if(value.substring(0,11)==="one_colour_"){
                                newColMode.textContent = "one_colour";
                                newColModeThis.textContent = "one_colour";
                                var newOneCol = newDoc.createElement("one_colour");
                                var newOneColThis = this.xmlDoc.createElement("one_colour");
                                colparam.appendChild(newOneCol);
                                dispobjs[iobj].getElementsByTagName("colour_parameters")[0].appendChild(newOneColThis);
                                newOneCol.textContent = "#"+value.substring(11);
                                newOneColThis.textContent = "#"+value.substring(11);
                                addToColours("#"+value.substring(11),newDoc);
                                addToColours("#"+value.substring(11),this.xmlDoc);
                            } else {
                                newColMode.textContent = value;
                                newColModeThis.textContent = value;
                            }
                            newColParamThis.appendChild(newColModeThis);
                            dispobjs[iobj].appendChild(newColParamThis);
                        }
                    }
                    if(property === "Style"){
                        // Sigh we could have <style>SPLINE</style> or <style_parameters><style_mode>SPLINE</style_mode></style_parameters>
                        // I'll deal with both.
                        if(thisDispCopy.getElementsByTagName("style_parameters").length>0&&thisDispCopy.getElementsByTagName("style_parameters")[0].getElementsByTagName("style_mode").length>0){
                            var styleparam = thisDispCopy.getElementsByTagName("style_parameters")[0];
                            styleparam.getElementsByTagName("style_mode")[0].textContent = value;
                            // Is it wise to do this so early?
                            dispobjs[iobj].getElementsByTagName("style_parameters")[0].getElementsByTagName("style_mode")[0].textContent = value;
                        } else if(thisDispCopy.getElementsByTagName("style").length>0) {
                            thisDispCopy.getElementsByTagName("style")[0].textContent = value;
                            dispobjs[iobj].getElementsByTagName("style")[0].textContent = value;
                        } else {
                            var newStyleMode = newDoc.createElement("style");
                            newStyleMode.textContent = value;
                            thisDispCopy.appendChild(newStyleMode);
                            var newStyleModeThis = this.xmlDoc.createElement("style");
                            newStyleModeThis.textContent = value;
                            dispobjs[iobj].appendChild(newStyleModeThis);
                        }
                    }
                    if(property === "Selection"){
                        // This is trickier!
                        if(value.substring(0,14)==="complex_neighb"){
                            console.log("Complex neighbourhood!!!");
                            if(thisDispCopy.getElementsByTagName("selection_parameters").length>0){
                                var selparam = thisDispCopy.getElementsByTagName("selection_parameters")[0];
                                if(selparam.getElementsByTagName("select").length>0){
                                    selparam.getElementsByTagName("select")[0].textContent = "neighb";
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].getElementsByTagName("select")[0].textContent = "neighb";
                                } else {
                                    var newSelect = newDoc.createElement("select");
                                    newSelect.textContent = "neighb";
                                    selparam.appendChild(newSelect);
                                    var newSelectThis = newDoc.createElement("select");
                                    newSelect.textContent = "neighb";
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].appendChild(newSelectThis);
                                }
                                if(selparam.getElementsByTagName("neighb_sel").length>0){
                                    selparam.getElementsByTagName("neighb_sel")[0].textContent = value.substring(15);
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].getElementsByTagName("neighb_sel")[0].textContent = value.substring(15);
                                } else {
                                    var newCid = newDoc.createElement("neighb_sel");
                                    newCid.textContent = value.substring(15);
                                    selparam.appendChild(newCid);
                                    var newCidThis = this.xmlDoc.createElement("neighb_sel");
                                    newCidThis.textContent = value.substring(15);
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].appendChild(newCidThis);
                                }
                            } else {
                                var newSelParam = newDoc.createElement("selection_parameters");
                                var newSelect = newDoc.createElement("select");
                                var newCid = newDoc.createElement("neighb_sel");
                                newSelect.textContent = "neighb";
                                newCid.textContent = value.substring(15);
                                newSelParam.appendChild(newSelect);
                                newSelParam.appendChild(newCid);
                                thisDispCopy.appendChild(newSelParam);
                                var newSelParamThis = this.xmlDoc.createElement("selection_parameters");
                                var newSelectThis = this.xmlDoc.createElement("select");
                                var newCidThis = this.xmlDoc.createElement("neighb_sel");
                                newSelectThis.textContent = "neighb";
                                newCidThis.textContent = value.substring(15);
                                newSelParamThis.appendChild(newSelectThis);
                                newSelParamThis.appendChild(newCidThis);
                                dispobjs[iobj].appendChild(newSelParamThis);
                            }
                        } else if(value.substring(0,7)==="neighb_"){
                            if(thisDispCopy.getElementsByTagName("selection_parameters").length>0){
                                var selparam = thisDispCopy.getElementsByTagName("selection_parameters")[0];
                                if(selparam.getElementsByTagName("select").length>0){
                                    selparam.getElementsByTagName("select")[0].textContent = "neighb";
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].getElementsByTagName("select")[0].textContent = "neighb";
                                } else {
                                    var newSelect = newDoc.createElement("select");
                                    newSelect.textContent = "neighb";
                                    selparam.appendChild(newSelect);
                                    var newSelectThis = newDoc.createElement("select");
                                    newSelect.textContent = "neighb";
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].appendChild(newSelectThis);
                                }
                                if(selparam.getElementsByTagName("neighb_sel").length>0){
                                    selparam.getElementsByTagName("neighb_sel")[0].textContent = value.substring(7);
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].getElementsByTagName("neighb_sel")[0].textContent = value.substring(7);
                                } else {
                                    var newCid = newDoc.createElement("neighb_sel");
                                    newCid.textContent = value.substring(7);
                                    selparam.appendChild(newCid);
                                    var newCidThis = this.xmlDoc.createElement("neighb_sel");
                                    newCidThis.textContent = value.substring(7);
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].appendChild(newCidThis);
                                }
                            } else {
                                var newSelParam = newDoc.createElement("selection_parameters");
                                var newSelect = newDoc.createElement("select");
                                var newCid = newDoc.createElement("neighb_sel");
                                newSelect.textContent = "neighb";
                                newCid.textContent = value.substring(7);
                                newSelParam.appendChild(newSelect);
                                newSelParam.appendChild(newCid);
                                thisDispCopy.appendChild(newSelParam);
                                var newSelParamThis = this.xmlDoc.createElement("selection_parameters");
                                var newSelectThis = this.xmlDoc.createElement("select");
                                var newCidThis = this.xmlDoc.createElement("neighb_sel");
                                newSelectThis.textContent = "neighb";
                                newCidThis.textContent = value.substring(7);
                                newSelParamThis.appendChild(newSelectThis);
                                newSelParamThis.appendChild(newCidThis);
                                dispobjs[iobj].appendChild(newSelParamThis);
                            }
                        } else if(value.substring(0,6)==="sites_") {
                            function fillSites(node,theDoc){
                                while (node.firstChild) {
                                    node.removeChild(node.firstChild);
                                }
                                var sites = value.substring(6).split(" or ");
                                for(isite=0;isite<sites.length;isite++){
                                    var site = theDoc.createElement("site");
                                    site.textContent = sites[isite];
                                    node.appendChild(site);
                                }
                            }
                            if(thisDispCopy.getElementsByTagName("selection_parameters").length>0){
                                var selparam = thisDispCopy.getElementsByTagName("selection_parameters")[0];
                                if(selparam.getElementsByTagName("select").length>0){
                                    selparam.getElementsByTagName("select")[0].textContent = "sites";
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].getElementsByTagName("select")[0].textContent = "sites";
                                } else {
                                    var newSelect = newDoc.createElement("select");
                                    newSelect.textContent = "sites";
                                    selparam.appendChild(newSelect);
                                    var newSelectThis = newDoc.createElement("select");
                                    newSelect.textContent = "sites";
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].appendChild(newSelectThis);
                                }
                                if(selparam.getElementsByTagName("sites").length>0){
                                    fillSites(selparam.getElementsByTagName("sites")[0],newDoc);
                                    fillSites(dispobjs[iobj].getElementsByTagName("selection_parameters")[0].getElementsByTagName("sites")[0],this.xmlDoc);
                                } else {
                                    var newCid = newDoc.createElement("sites");
                                    fillSites(newCi,newDocd);
                                    selparam.appendChild(newCid);
                                    var newCidThis = this.xmlDoc.createElement("sites");
                                    fillSites(newCidThis,this.xmlDoc);
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].appendChild(newCidThis);
                                }
                            } else {
                                var newSelParam = newDoc.createElement("selection_parameters");
                                var newSelect = newDoc.createElement("select");
                                var newCid = newDoc.createElement("sites");
                                newSelect.textContent = "sites";
                                fillSites(newCid,newDoc);
                                newSelParam.appendChild(newSelect);
                                newSelParam.appendChild(newCid);
                                thisDispCopy.appendChild(newSelParam);
                                var newSelParamThis = this.xmlDoc.createElement("selection_parameters");
                                var newSelectThis = this.xmlDoc.createElement("select");
                                var newCidThis = this.xmlDoc.createElement("sites");
                                newSelectThis.textContent = "sites";
                                fillSites(newCidThis,this.xmlDoc);
                                newSelParamThis.appendChild(newSelectThis);
                                newSelParamThis.appendChild(newCidThis);
                                dispobjs[iobj].appendChild(newSelParamThis);
                            }
                        } else {
                            if(thisDispCopy.getElementsByTagName("selection_parameters").length>0){
                                var selparam = thisDispCopy.getElementsByTagName("selection_parameters")[0];
                                if(selparam.getElementsByTagName("select").length>0){
                                    selparam.getElementsByTagName("select")[0].textContent = "cid";
                                    // Is it wise to do this so early?
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].getElementsByTagName("select")[0].textContent = "cid";
                                } else {
                                    var newSelect = newDoc.createElement("select");
                                    newSelect.textContent = "cid";
                                    selparam.appendChild(newSelect);
                                    var newSelectThis = newDoc.createElement("select");
                                    newSelect.textContent = "cid";
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].appendChild(newSelectThis);
                                }
                                if(selparam.getElementsByTagName("cid").length>0){
                                    selparam.getElementsByTagName("cid")[0].textContent = value;
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].getElementsByTagName("cid")[0].textContent = value;
                                } else {
                                    var newCid = newDoc.createElement("cid");
                                    newCid.textContent = value;
                                    selparam.appendChild(newCid);
                                    var newCidThis = this.xmlDoc.createElement("cid");
                                    newCidThis.textContent = value;
                                    dispobjs[iobj].getElementsByTagName("selection_parameters")[0].appendChild(newCidThis);
                                }
                            } else {
                                var newSelParam = newDoc.createElement("selection_parameters");
                                var newSelect = newDoc.createElement("select");
                                var newCid = newDoc.createElement("cid");
                                newSelect.textContent = "cid";
                                newCid.textContent = value;
                                newSelParam.appendChild(newSelect);
                                newSelParam.appendChild(newCid);
                                thisDispCopy.appendChild(newSelParam);

                                var newSelParamThis = this.xmlDoc.createElement("selection_parameters");
                                var newSelectThis = this.xmlDoc.createElement("select");
                                var newCidThis = this.xmlDoc.createElement("cid");
                                newSelectThis.textContent = "cid";
                                newCidThis.textContent = value;
                                newSelParamThis.appendChild(newSelectThis);
                                newSelParamThis.appendChild(newCidThis);
                                dispobjs[iobj].appendChild(newSelParamThis);
                            }
                        }
                    }
                    thisDataCopy.appendChild(thisDispCopy);
                    newRoot.appendChild(thisDataCopy);

                    // This could/should probably all go in new function.
                    //console.log("Doing snippet stuff");
                    var snippet = (new XMLSerializer()).serializeToString(newDoc);
                    //console.log(snippet);
                    var blob  = new Blob([snippet], {type : "text/plain"});
                    var formData = new FormData();
                    if (typeof(middlewareToken) !== 'undefined'){
                        formData.append('csrfmiddlewaretoken',middlewareToken);
                    }
                    var theFile = guid()+".mgpic.xml";
                    var client = new XMLHttpRequest();
                    formData.append("upfile",blob,theFile);
                    client.onreadystatechange = function() {
                        if (client.readyState == 4 && client.status == 200) {
                            //console.log("Sent snippet");
                            var client2 = new XMLHttpRequest();
                            client2.onreadystatechange = function() {
                                if (client2.readyState == 4 && client2.status == 200) {
                                    console.log("Executed snippet");
                                    // FIXME In general create new child node, give newChild the uuid of old child, data[idat].replaceChild(newChild).
                                    // But we're OK, I think for colours - we've only change colour_mode.
                                    var jsondata=JSON.parse(client2.responseText);
                                    //console.log(jsondata.status);
                                    // idat had better only be 0. Not checked that yet.
                                    //console.log(jsondata.norm_tri.length);
                                    //console.log(jsondata.vert_tri.length);
                                    //console.log(jsondata.col_tri.length);
                                    //console.log(jsondata.idx_tri.length);
                                    //console.log(0);
                                    //console.log(jsondata.norm_tri[0].length);
                                    //console.log(jsondata.vert_tri[0].length);
                                    //console.log(jsondata.col_tri[0].length);
                                    //console.log(jsondata.idx_tri[0].length);
                                    //for(var idat=0;idat<jsondata.norm_tri.length;idat++){
                                        self.currentBufferIdx = num;
                                        var rssentries=jsondata.norm_tri[0];
                                        var norms = getEncodedData(rssentries);
                                        //console.log("norms");
                                        //console.log(norms.length);
                                        //console.log(norms);
                                        self.displayBuffers[num].clearBuffers();
                                        var thisName = jsondata.names[0];
                                        self.displayBuffers[self.currentBufferIdx].name_label = thisName;
                                        for (var i=0; i<norms.length; i++){
                                            self.createNormalBuffer(norms[i]);
                                        }
                                        var rssentries=jsondata.vert_tri[0];
                                        var tris = getEncodedData(rssentries);
                                        //console.log("tris");
                                        //console.log(tris.length);
                                        //console.log(tris);
                                        for (var i=0; i<tris.length; i++){
                                            self.createVertexBuffer(tris[i]);
                                        }
                                        var rssentries=jsondata.col_tri[0];
                                        var theColours = getEncodedData(rssentries);
                                        //console.log("theColours");
                                        //console.log(theColours.length);
                                        //console.log(theColours);
                                        for (var i=0; i<theColours.length; i++){
                                            self.createColourBuffer(theColours[i]);
                                        }
                                        var rssentries=jsondata.idx_tri[0];
                                        var idxs = getEncodedData(rssentries);
                                        for (var i=0; i<idxs.length; i++){
                                            self.createIndexBuffer(idxs[i]);
                                        }
                                        var rssentries=jsondata.prim_types[0];
                                        for (var i=0; i<rssentries.length; i++){
                                            self.displayBuffers[self.currentBufferIdx].bufferTypes.push(rssentries[i]);
                                        }
                                        var rssentries=getEncodedData(jsondata.sizes[0]);
                                        if(typeof(rssentries)!=="undefined"){
                                            console.log(rssentries);
                                            for (var i=0; i<rssentries.length; i++){
                                                self.createSizeBuffer(rssentries[i]);
                                            }
                                        }
                                        var atoms = jsondata.atoms[0];
                                        self.displayBuffers[self.currentBufferIdx].atoms = atoms;
                                        self.displayBuffers[self.currentBufferIdx].isDirty = true;
                                        self.buildBuffers();
                                        self.drawScene();
                                        // FIXME - Wait for now, this does not set correct menu entries.
                                        self.div.dispatchEvent(self.statusChangedEvent);
                                    //}
                                }
                            }
                            client2.open("GET", "get_triangles?pdb="+theFile, true);
                            client2.send(null);
                        }
                    }
                    client.open("POST", "/server_client.html", true);
                    client.send(formData);
                }
            }
            // FIXME - All these other objects. And Annotations.
            for(var iobj=0;iobj<hbobjs.length;iobj++){
                var theUUID = dataUUID+"_"+getNodeText(hbobjs[iobj].getElementsByTagName("uuid")[0]);
                if(theUUID && theUUID === this.ids[num]){
                    console.log("Update "+hbobjs[iobj]+" "+property+" "+value);
                }
            }
            for(var iobj=0;iobj<contactobjs.length;iobj++){
                var theUUID = dataUUID+"_"+getNodeText(contactobjs[iobj].getElementsByTagName("uuid")[0]);
                if(theUUID && theUUID === this.ids[num]){
                    console.log("Update "+contactobjs[iobj]+" "+property+" "+value);
                }
            }
            for(var iobj=0;iobj<surfobs.length;iobj++){
                var theUUID = dataUUID+"_"+getNodeText(surfobs[iobj].getElementsByTagName("uuid")[0]);
                if(theUUID && theUUID === this.ids[num]){
                    console.log("Update "+surfobs[iobj]+" "+property+" "+value);
                }
            }
        }
    }
}

MGWebGL.prototype.loadPDBCode = function (pdbid) {

    var self = this;

    function loadPDBCodeMain(){
        var client = new XMLHttpRequest();
        var theFile = guid()+pdbid+".ent";
        // Create a FormData instance
        var formData = new FormData();
        if (typeof(middlewareToken) !== 'undefined'){
            formData.append('csrfmiddlewaretoken',middlewareToken);
        }
        // Add the file
        //console.log("Adding to form 2 "+file.files[0]);
        formData.append("pdbcode",pdbid);
        formData.append("pdbfilename",theFile);

        // Check the response status
        client.onreadystatechange = function() 
        {
            if (client.readyState == 4 && client.status == 200) 
            {
                self.mygetrequest.open("GET", "get_triangles?pdb="+theFile, true);
                self.mygetrequest.send(null);
            }
        }

        client.open("POST", "/server_client.html", true);
        client.send(formData);
    }

    var myVar = setInterval(function(){ myTimer() }, 500);
    function myTimer() {
        console.log(self.ready);
        if(self.ready){
            clearInterval(myVar);
            loadPDBCodeMain();

        }
    }

}

MGWebGL.prototype.loadFile = function (theFile) {
    /* This *must* be a file visible to the server! */

    var self = this;

    function loadFileMain(){
        self.mygetrequest.open("GET", "get_triangles?absolutePath=true&pdb="+theFile, true);
        self.mygetrequest.send(null);
    }

    var myVar = setInterval(function(){ myTimer() }, 500);
    function myTimer() {
        console.log(self.ready);
        if(self.ready){
            clearInterval(myVar);
            loadFileMain();

        }
    }

}

MGWebGL.prototype.initTextBuffers = function () {
    this.displayBuffers[0].textNormalBuffer = this.gl.createBuffer();
    this.displayBuffers[0].textPositionBuffer = this.gl.createBuffer();
    this.displayBuffers[0].textColourBuffer = this.gl.createBuffer();
    this.displayBuffers[0].textTexCoordBuffer = this.gl.createBuffer();
    this.displayBuffers[0].textIndexesBuffer = this.gl.createBuffer();

    this.displayBuffers[0].clickLinePositionBuffer = this.gl.createBuffer();
    this.displayBuffers[0].clickLineColourBuffer = this.gl.createBuffer();
    this.displayBuffers[0].clickLineIndexesBuffer = this.gl.createBuffer();
}

/*
MGWebGL.prototype.resizeCanvas = function () {
        //this.canvas.width  = parseInt(0.8*window.innerWidth);
        //this.canvas.height = window.innerHeight;
        this.gl.viewportWidth = this.canvas.width;
        this.gl.viewportHeight = this.canvas.height;
        this.drawScene();
    }
*/

MGWebGL.prototype.set_clip_range = function(clipStart,clipEnd,update){
    //console.log("Clip "+clipStart+" "+clipEnd);
    if(typeof(this.gl)==='undefined'){
        return;
    }
    this.gl_clipPlane0[3] = -500.0-clipStart;
    this.gl_clipPlane1[3] = 500.0+clipEnd;
    if(update)
        this.drawScene();
}

MGWebGL.prototype.set_fog_range = function(fogStart,fogEnd,update){
    this.gl_fog_start = fogStart;
    this.gl_fog_end = fogEnd;
    //console.log("Fog "+this.gl_fog_start+" "+this.gl_fog_end);
    if(typeof(this.gl)==='undefined'){
        return;
    }
    if(update)
        this.drawScene();
}

MGWebGL.prototype.initTriangleShadowShaders = function(vertexShaderTriangleShadow,fragmentShaderTriangleShadow) {
    this.shaderProgramTriangleShadow = this.gl.createProgram();

    this.gl.attachShader(this.shaderProgramTriangleShadow, vertexShaderTriangleShadow);
    this.gl.attachShader(this.shaderProgramTriangleShadow, fragmentShaderTriangleShadow);
    this.gl.bindAttribLocation(this.shaderProgramTriangleShadow, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramTriangleShadow, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramTriangleShadow, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramTriangleShadow, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramTriangleShadow);

    if (!this.gl.getProgramParameter(this.shaderProgramTriangleShadow, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initTriangleShadowShaders)");
    }

    this.gl.useProgram(this.shaderProgramTriangleShadow);

    this.shaderProgramTriangleShadow.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramTriangleShadow, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramTriangleShadow.vertexNormalAttribute);

    this.shaderProgramTriangleShadow.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramTriangleShadow, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramTriangleShadow.vertexPositionAttribute);

    this.shaderProgramTriangleShadow.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramTriangleShadow, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramTriangleShadow.vertexColourAttribute);

    this.shaderProgramTriangleShadow.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "uPMatrix");
    this.shaderProgramTriangleShadow.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "uMVMatrix");
    this.shaderProgramTriangleShadow.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "uMVINVMatrix");
    this.shaderProgramTriangleShadow.textureMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "TextureMatrix");

    this.shaderProgramTriangleShadow.fog_start = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "fog_start");
    this.shaderProgramTriangleShadow.fog_end = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "fog_end");
    this.shaderProgramTriangleShadow.fogColour = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "fogColour");

    this.shaderProgramTriangleShadow.clipPlane0 = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "clipPlane0");
    this.shaderProgramTriangleShadow.clipPlane1 = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "clipPlane1");
    this.shaderProgramTriangleShadow.clipPlane2 = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "clipPlane2");
    this.shaderProgramTriangleShadow.clipPlane3 = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "clipPlane3");
    this.shaderProgramTriangleShadow.clipPlane4 = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "clipPlane4");
    this.shaderProgramTriangleShadow.clipPlane5 = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "clipPlane5");
    this.shaderProgramTriangleShadow.clipPlane6 = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "clipPlane6");
    this.shaderProgramTriangleShadow.clipPlane7 = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "clipPlane7");
    this.shaderProgramTriangleShadow.nClipPlanes = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "nClipPlanes");

    this.shaderProgramTriangleShadow.cursorPos = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "cursorPos");

    this.shaderProgramTriangleShadow.shinyBack = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "shinyBack");
    this.shaderProgramTriangleShadow.defaultColour = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "defaultColour");
    this.shaderProgramTriangleShadow.backColour = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "backColour");
    
    this.shaderProgramTriangleShadow.light_positions = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "light_positions");
    this.shaderProgramTriangleShadow.light_colours_ambient = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "light_colours_ambient");
    this.shaderProgramTriangleShadow.light_colours_specular = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "light_colours_specular");
    this.shaderProgramTriangleShadow.light_colours_diffuse = this.gl.getUniformLocation(this.shaderProgramTriangleShadow, "light_colours_diffuse");
    
}

MGWebGL.prototype.initShadowShaders = function(vertexShaderShadow,fragmentShaderShadow) {
    this.shaderProgramShadow = this.gl.createProgram();

    this.gl.attachShader(this.shaderProgramShadow, vertexShaderShadow);
    this.gl.attachShader(this.shaderProgramShadow, fragmentShaderShadow);
    this.gl.bindAttribLocation(this.shaderProgramShadow, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramShadow, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramShadow, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramShadow, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramShadow);

    if (!this.gl.getProgramParameter(this.shaderProgramShadow, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initShadowShaders)");
    }

    this.gl.useProgram(this.shaderProgramShadow);

    this.shaderProgramShadow.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramShadow, "aVertexNormal");
    console.log(this.shaderProgramShadow.vertexNormalAttribute);
    if(this.shaderProgramShadow.vertexNormalAttribute>-1) this.gl.enableVertexAttribArray(this.shaderProgramShadow.vertexNormalAttribute);

    this.shaderProgramShadow.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramShadow, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramShadow.vertexPositionAttribute);

    this.shaderProgramShadow.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramShadow, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramShadow.vertexColourAttribute);

    this.shaderProgramShadow.vertexTextureAttribute = this.gl.getAttribLocation(this.shaderProgramShadow, "aVertexTexture");
    if(this.shaderProgramShadow.vertexTextureAttribute>-1) this.gl.enableVertexAttribArray(this.shaderProgramShadow.vertexTextureAttribute);

    this.shaderProgramShadow.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramShadow, "uPMatrix");
    this.shaderProgramShadow.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramShadow, "uMVMatrix");
    this.shaderProgramShadow.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramShadow, "uMVINVMatrix");

    this.shaderProgramShadow.fog_start = this.gl.getUniformLocation(this.shaderProgramShadow, "fog_start");
    this.shaderProgramShadow.fog_end = this.gl.getUniformLocation(this.shaderProgramShadow, "fog_end");
    this.shaderProgramShadow.fogColour = this.gl.getUniformLocation(this.shaderProgramShadow, "fogColour");

    this.shaderProgramShadow.clipPlane0 = this.gl.getUniformLocation(this.shaderProgramShadow, "clipPlane0");
    this.shaderProgramShadow.clipPlane1 = this.gl.getUniformLocation(this.shaderProgramShadow, "clipPlane1");
    this.shaderProgramShadow.clipPlane2 = this.gl.getUniformLocation(this.shaderProgramShadow, "clipPlane2");
    this.shaderProgramShadow.clipPlane3 = this.gl.getUniformLocation(this.shaderProgramShadow, "clipPlane3");
    this.shaderProgramShadow.clipPlane4 = this.gl.getUniformLocation(this.shaderProgramShadow, "clipPlane4");
    this.shaderProgramShadow.clipPlane5 = this.gl.getUniformLocation(this.shaderProgramShadow, "clipPlane5");
    this.shaderProgramShadow.clipPlane6 = this.gl.getUniformLocation(this.shaderProgramShadow, "clipPlane6");
    this.shaderProgramShadow.clipPlane7 = this.gl.getUniformLocation(this.shaderProgramShadow, "clipPlane7");
    this.shaderProgramShadow.nClipPlanes = this.gl.getUniformLocation(this.shaderProgramShadow, "nClipPlanes");

}

MGWebGL.prototype.initRenderFrameBufferShaders = function(vertexShaderRenderFrameBuffer,fragmentShaderRenderFrameBuffer) {
    this.shaderProgramRenderFrameBuffer = this.gl.createProgram();

    this.gl.attachShader(this.shaderProgramRenderFrameBuffer, vertexShaderRenderFrameBuffer);
    this.gl.attachShader(this.shaderProgramRenderFrameBuffer, fragmentShaderRenderFrameBuffer);
    this.gl.bindAttribLocation(this.shaderProgramRenderFrameBuffer, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramRenderFrameBuffer, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramRenderFrameBuffer, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramRenderFrameBuffer, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramRenderFrameBuffer);

    if (!this.gl.getProgramParameter(this.shaderProgramRenderFrameBuffer, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initRenderFrameBufferShaders)");
    }

    this.gl.useProgram(this.shaderProgramRenderFrameBuffer);

    this.shaderProgramRenderFrameBuffer.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramRenderFrameBuffer, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramRenderFrameBuffer.vertexNormalAttribute);

    this.shaderProgramRenderFrameBuffer.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramRenderFrameBuffer, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramRenderFrameBuffer.vertexPositionAttribute);

    this.shaderProgramRenderFrameBuffer.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramRenderFrameBuffer, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramRenderFrameBuffer.vertexColourAttribute);

    this.shaderProgramRenderFrameBuffer.vertexTextureAttribute = this.gl.getAttribLocation(this.shaderProgramRenderFrameBuffer, "aVertexTexture");
    this.gl.enableVertexAttribArray(this.shaderProgramRenderFrameBuffer.vertexTextureAttribute);

    this.shaderProgramRenderFrameBuffer.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "uPMatrix");
    this.shaderProgramRenderFrameBuffer.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "uMVMatrix");
    this.shaderProgramRenderFrameBuffer.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "uMVINVMatrix");

    this.shaderProgramRenderFrameBuffer.fog_start = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "fog_start");
    this.shaderProgramRenderFrameBuffer.fog_end = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "fog_end");
    this.shaderProgramRenderFrameBuffer.fogColour = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "fogColour");

    this.shaderProgramRenderFrameBuffer.clipPlane0 = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "clipPlane0");
    this.shaderProgramRenderFrameBuffer.clipPlane1 = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "clipPlane1");
    this.shaderProgramRenderFrameBuffer.clipPlane2 = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "clipPlane2");
    this.shaderProgramRenderFrameBuffer.clipPlane3 = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "clipPlane3");
    this.shaderProgramRenderFrameBuffer.clipPlane4 = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "clipPlane4");
    this.shaderProgramRenderFrameBuffer.clipPlane5 = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "clipPlane5");
    this.shaderProgramRenderFrameBuffer.clipPlane6 = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "clipPlane6");
    this.shaderProgramRenderFrameBuffer.clipPlane7 = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "clipPlane7");
    this.shaderProgramRenderFrameBuffer.nClipPlanes = this.gl.getUniformLocation(this.shaderProgramRenderFrameBuffer, "nClipPlanes");

}

MGWebGL.prototype.initTextBackgroundShaders = function(vertexShaderTextBackground,fragmentShaderTextBackground) {

    this.shaderProgramTextBackground = this.gl.createProgram();

    this.gl.attachShader(this.shaderProgramTextBackground, vertexShaderTextBackground);
    this.gl.attachShader(this.shaderProgramTextBackground, fragmentShaderTextBackground);
    this.gl.bindAttribLocation(this.shaderProgramTextBackground, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramTextBackground, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramTextBackground, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramTextBackground, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramTextBackground);

    if (!this.gl.getProgramParameter(this.shaderProgramTextBackground, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initTextBackgroundShaders)");
    }

    this.gl.useProgram(this.shaderProgramTextBackground);

    this.shaderProgramTextBackground.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramTextBackground, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramTextBackground.vertexNormalAttribute);

    this.shaderProgramTextBackground.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramTextBackground, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramTextBackground.vertexPositionAttribute);

    this.shaderProgramTextBackground.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramTextBackground, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramTextBackground.vertexColourAttribute);

    this.shaderProgramTextBackground.vertexTextureAttribute = this.gl.getAttribLocation(this.shaderProgramTextBackground, "aVertexTexture");
    this.gl.enableVertexAttribArray(this.shaderProgramTextBackground.vertexTextureAttribute);

    this.shaderProgramTextBackground.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTextBackground, "uPMatrix");
    this.shaderProgramTextBackground.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTextBackground, "uMVMatrix");
    this.shaderProgramTextBackground.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTextBackground, "uMVINVMatrix");

    this.shaderProgramTextBackground.fog_start = this.gl.getUniformLocation(this.shaderProgramTextBackground, "fog_start");
    this.shaderProgramTextBackground.fog_end = this.gl.getUniformLocation(this.shaderProgramTextBackground, "fog_end");
    this.shaderProgramTextBackground.fogColour = this.gl.getUniformLocation(this.shaderProgramTextBackground, "fogColour");

    this.shaderProgramTextBackground.clipPlane0 = this.gl.getUniformLocation(this.shaderProgramTextBackground, "clipPlane0");
    this.shaderProgramTextBackground.clipPlane1 = this.gl.getUniformLocation(this.shaderProgramTextBackground, "clipPlane1");
    this.shaderProgramTextBackground.clipPlane2 = this.gl.getUniformLocation(this.shaderProgramTextBackground, "clipPlane2");
    this.shaderProgramTextBackground.clipPlane3 = this.gl.getUniformLocation(this.shaderProgramTextBackground, "clipPlane3");
    this.shaderProgramTextBackground.clipPlane4 = this.gl.getUniformLocation(this.shaderProgramTextBackground, "clipPlane4");
    this.shaderProgramTextBackground.clipPlane5 = this.gl.getUniformLocation(this.shaderProgramTextBackground, "clipPlane5");
    this.shaderProgramTextBackground.clipPlane6 = this.gl.getUniformLocation(this.shaderProgramTextBackground, "clipPlane6");
    this.shaderProgramTextBackground.clipPlane7 = this.gl.getUniformLocation(this.shaderProgramTextBackground, "clipPlane7");
    this.shaderProgramTextBackground.nClipPlanes = this.gl.getUniformLocation(this.shaderProgramTextBackground, "nClipPlanes");
}

MGWebGL.prototype.initShaders = function(vertexShader,fragmentShader) {

    this.shaderProgram = this.gl.createProgram();

    this.gl.attachShader(this.shaderProgram, vertexShader);
    this.gl.attachShader(this.shaderProgram, fragmentShader);
    this.gl.bindAttribLocation(this.shaderProgram, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgram, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgram, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgram, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgram);

    if (!this.gl.getProgramParameter(this.shaderProgram, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initShaders)");
    }

    this.gl.useProgram(this.shaderProgram);

    this.shaderProgram.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgram, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgram.vertexNormalAttribute);

    this.shaderProgram.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgram, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgram.vertexPositionAttribute);

    this.shaderProgram.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgram, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgram.vertexColourAttribute);

    this.shaderProgram.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgram, "uPMatrix");
    this.shaderProgram.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgram, "uMVMatrix");
    this.shaderProgram.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgram, "uMVINVMatrix");

    this.shaderProgram.fog_start = this.gl.getUniformLocation(this.shaderProgram, "fog_start");
    this.shaderProgram.fog_end = this.gl.getUniformLocation(this.shaderProgram, "fog_end");
    this.shaderProgram.fogColour = this.gl.getUniformLocation(this.shaderProgram, "fogColour");

    this.shaderProgram.clipPlane0 = this.gl.getUniformLocation(this.shaderProgram, "clipPlane0");
    this.shaderProgram.clipPlane1 = this.gl.getUniformLocation(this.shaderProgram, "clipPlane1");
    this.shaderProgram.clipPlane2 = this.gl.getUniformLocation(this.shaderProgram, "clipPlane2");
    this.shaderProgram.clipPlane3 = this.gl.getUniformLocation(this.shaderProgram, "clipPlane3");
    this.shaderProgram.clipPlane4 = this.gl.getUniformLocation(this.shaderProgram, "clipPlane4");
    this.shaderProgram.clipPlane5 = this.gl.getUniformLocation(this.shaderProgram, "clipPlane5");
    this.shaderProgram.clipPlane6 = this.gl.getUniformLocation(this.shaderProgram, "clipPlane6");
    this.shaderProgram.clipPlane7 = this.gl.getUniformLocation(this.shaderProgram, "clipPlane7");
    this.shaderProgram.nClipPlanes = this.gl.getUniformLocation(this.shaderProgram, "nClipPlanes");

    this.shaderProgram.cursorPos = this.gl.getUniformLocation(this.shaderProgram, "cursorPos");

    this.shaderProgram.shinyBack = this.gl.getUniformLocation(this.shaderProgram, "shinyBack");
    this.shaderProgram.defaultColour = this.gl.getUniformLocation(this.shaderProgram, "defaultColour");
    this.shaderProgram.backColour = this.gl.getUniformLocation(this.shaderProgram, "backColour");
    
    this.shaderProgram.light_positions = this.gl.getUniformLocation(this.shaderProgram, "light_positions");
    this.shaderProgram.light_colours_ambient = this.gl.getUniformLocation(this.shaderProgram, "light_colours_ambient");
    this.shaderProgram.light_colours_specular = this.gl.getUniformLocation(this.shaderProgram, "light_colours_specular");
    this.shaderProgram.light_colours_diffuse = this.gl.getUniformLocation(this.shaderProgram, "light_colours_diffuse");
}

MGWebGL.prototype.initThickLineShaders = function(vertexShader,fragmentShader) {

    this.shaderProgramThickLines = this.gl.createProgram();
    this.gl.attachShader(this.shaderProgramThickLines, vertexShader);
    this.gl.attachShader(this.shaderProgramThickLines, fragmentShader);
    this.gl.bindAttribLocation(this.shaderProgramThickLines, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramThickLines, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramThickLines, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramThickLines, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramThickLines);

    if (!this.gl.getProgramParameter(this.shaderProgramThickLines, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initThickLineShaders)");
    }

    this.gl.useProgram(this.shaderProgramThickLines);

    this.shaderProgramThickLines.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramThickLines, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramThickLines.vertexPositionAttribute);

    this.shaderProgramThickLines.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramThickLines, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramThickLines.vertexColourAttribute);

    this.shaderProgramThickLines.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramThickLines, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramThickLines.vertexNormalAttribute);

    this.shaderProgramThickLines.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramThickLines, "uPMatrix");
    this.shaderProgramThickLines.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramThickLines, "uMVMatrix");
    this.shaderProgramThickLines.screenZ = this.gl.getUniformLocation(this.shaderProgramThickLines, "screenZ");

    this.shaderProgramThickLines.fog_start = this.gl.getUniformLocation(this.shaderProgramThickLines, "fog_start");
    this.shaderProgramThickLines.fog_end = this.gl.getUniformLocation(this.shaderProgramThickLines, "fog_end");
    this.shaderProgramThickLines.fogColour = this.gl.getUniformLocation(this.shaderProgramThickLines, "fogColour");
    this.shaderProgramThickLines.clipPlane0 = this.gl.getUniformLocation(this.shaderProgramThickLines, "clipPlane0");
    this.shaderProgramThickLines.clipPlane1 = this.gl.getUniformLocation(this.shaderProgramThickLines, "clipPlane1");
    this.shaderProgramThickLines.clipPlane2 = this.gl.getUniformLocation(this.shaderProgramThickLines, "clipPlane2");
    this.shaderProgramThickLines.clipPlane3 = this.gl.getUniformLocation(this.shaderProgramThickLines, "clipPlane3");
    this.shaderProgramThickLines.clipPlane4 = this.gl.getUniformLocation(this.shaderProgramThickLines, "clipPlane4");
    this.shaderProgramThickLines.clipPlane5 = this.gl.getUniformLocation(this.shaderProgramThickLines, "clipPlane5");
    this.shaderProgramThickLines.clipPlane6 = this.gl.getUniformLocation(this.shaderProgramThickLines, "clipPlane6");
    this.shaderProgramThickLines.clipPlane7 = this.gl.getUniformLocation(this.shaderProgramThickLines, "clipPlane7");
    this.shaderProgramThickLines.nClipPlanes = this.gl.getUniformLocation(this.shaderProgramThickLines, "nClipPlanes");

    this.shaderProgramThickLines.pixelZoom = this.gl.getUniformLocation(this.shaderProgramThickLines, "pixelZoom");
}

MGWebGL.prototype.initLineShaders = function(vertexShader,fragmentShader) {

    this.shaderProgramLines = this.gl.createProgram();
    this.gl.attachShader(this.shaderProgramLines, vertexShader);
    this.gl.attachShader(this.shaderProgramLines, fragmentShader);
    this.gl.bindAttribLocation(this.shaderProgramLines, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramLines, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramLines, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramLines, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramLines);

    if (!this.gl.getProgramParameter(this.shaderProgramLines, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initLineShaders)");
    }

    this.gl.useProgram(this.shaderProgramLines);

    this.shaderProgramLines.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramLines, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramLines.vertexPositionAttribute);

    this.shaderProgramLines.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramLines, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramLines.vertexColourAttribute);

    //this.shaderProgramLines.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramLines, "aVertexNormal");
    //this.gl.enableVertexAttribArray(this.shaderProgramLines.vertexNormalAttribute);

    this.shaderProgramLines.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramLines, "uPMatrix");
    this.shaderProgramLines.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramLines, "uMVMatrix");
    this.shaderProgramLines.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramLines, "uMVINVMatrix");

    this.shaderProgramLines.fog_start = this.gl.getUniformLocation(this.shaderProgramLines, "fog_start");
    this.shaderProgramLines.fog_end = this.gl.getUniformLocation(this.shaderProgramLines, "fog_end");
    this.shaderProgramLines.fogColour = this.gl.getUniformLocation(this.shaderProgramLines, "fogColour");
    this.shaderProgramLines.clipPlane0 = this.gl.getUniformLocation(this.shaderProgramLines, "clipPlane0");
    this.shaderProgramLines.clipPlane1 = this.gl.getUniformLocation(this.shaderProgramLines, "clipPlane1");
    this.shaderProgramLines.clipPlane2 = this.gl.getUniformLocation(this.shaderProgramLines, "clipPlane2");
    this.shaderProgramLines.clipPlane3 = this.gl.getUniformLocation(this.shaderProgramLines, "clipPlane3");
    this.shaderProgramLines.clipPlane4 = this.gl.getUniformLocation(this.shaderProgramLines, "clipPlane4");
    this.shaderProgramLines.clipPlane5 = this.gl.getUniformLocation(this.shaderProgramLines, "clipPlane5");
    this.shaderProgramLines.clipPlane6 = this.gl.getUniformLocation(this.shaderProgramLines, "clipPlane6");
    this.shaderProgramLines.clipPlane7 = this.gl.getUniformLocation(this.shaderProgramLines, "clipPlane7");
    this.shaderProgramLines.nClipPlanes = this.gl.getUniformLocation(this.shaderProgramLines, "nClipPlanes");

}

MGWebGL.prototype.initPerfectSphereShadowShaders = function(vertexShader,fragmentShader) {
    this.shaderProgramPerfectSpheresShadow = this.gl.createProgram();
    this.gl.attachShader(this.shaderProgramPerfectSpheresShadow, vertexShader);
    this.gl.attachShader(this.shaderProgramPerfectSpheresShadow, fragmentShader);
    this.gl.bindAttribLocation(this.shaderProgramPerfectSpheresShadow, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramPerfectSpheresShadow, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramPerfectSpheresShadow, 2, "aVertexNormal");
    this.gl.linkProgram(this.shaderProgramPerfectSpheresShadow);

    if (!this.gl.getProgramParameter(this.shaderProgramPerfectSpheresShadow, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initPerfectSphereShadowShaders)");
    }

    this.gl.useProgram(this.shaderProgramPerfectSpheresShadow);

    this.shaderProgramPerfectSpheresShadow.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramPerfectSpheresShadow, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramPerfectSpheresShadow.vertexPositionAttribute);

    this.shaderProgramPerfectSpheresShadow.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramPerfectSpheresShadow, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramPerfectSpheresShadow.vertexNormalAttribute);

    this.shaderProgramPerfectSpheresShadow.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramPerfectSpheresShadow, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramPerfectSpheresShadow.vertexColourAttribute);

    this.shaderProgramPerfectSpheresShadow.vertexTextureAttribute = this.gl.getAttribLocation(this.shaderProgramPerfectSpheresShadow, "aVertexTexture");
    this.gl.enableVertexAttribArray(this.shaderProgramPerfectSpheresShadow.vertexTextureAttribute);

    this.shaderProgramPerfectSpheresShadow.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "uPMatrix");
    this.shaderProgramPerfectSpheresShadow.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "uMVMatrix");
    this.shaderProgramPerfectSpheresShadow.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "uMVINVMatrix");
    this.shaderProgramPerfectSpheresShadow.textureMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "TextureMatrix");

    this.shaderProgramPerfectSpheresShadow.fog_start = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "fog_start");
    this.shaderProgramPerfectSpheresShadow.fog_end = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "fog_end");
    this.shaderProgramPerfectSpheresShadow.fogColour = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "fogColour");

    this.shaderProgramPerfectSpheresShadow.offset = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "offset");
    this.shaderProgramPerfectSpheresShadow.size = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "size");
    this.shaderProgramPerfectSpheresShadow.scaleMatrix = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "scaleMatrix");
    
    this.shaderProgramPerfectSpheresShadow.light_positions = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "light_positions");
    this.shaderProgramPerfectSpheresShadow.light_colours_ambient = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "light_colours_ambient");
    this.shaderProgramPerfectSpheresShadow.light_colours_specular = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "light_colours_specular");
    this.shaderProgramPerfectSpheresShadow.light_colours_diffuse = this.gl.getUniformLocation(this.shaderProgramPerfectSpheresShadow, "light_colours_diffuse");
}

MGWebGL.prototype.initPerfectSphereShaders = function(vertexShader,fragmentShader) {
    this.shaderProgramPerfectSpheres = this.gl.createProgram();
    this.gl.attachShader(this.shaderProgramPerfectSpheres, vertexShader);
    this.gl.attachShader(this.shaderProgramPerfectSpheres, fragmentShader);
    this.gl.bindAttribLocation(this.shaderProgramPerfectSpheres, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramPerfectSpheres, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramPerfectSpheres, 2, "aVertexNormal");
    this.gl.linkProgram(this.shaderProgramPerfectSpheres);

    if (!this.gl.getProgramParameter(this.shaderProgramPerfectSpheres, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initPerfectSphereShaders)");
    }

    this.gl.useProgram(this.shaderProgramPerfectSpheres);

    this.shaderProgramPerfectSpheres.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramPerfectSpheres, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramPerfectSpheres.vertexPositionAttribute);

    this.shaderProgramPerfectSpheres.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramPerfectSpheres, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramPerfectSpheres.vertexNormalAttribute);

    this.shaderProgramPerfectSpheres.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramPerfectSpheres, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramPerfectSpheres.vertexColourAttribute);

    this.shaderProgramPerfectSpheres.vertexTextureAttribute = this.gl.getAttribLocation(this.shaderProgramPerfectSpheres, "aVertexTexture");
    this.gl.enableVertexAttribArray(this.shaderProgramPerfectSpheres.vertexTextureAttribute);

    this.shaderProgramPerfectSpheres.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "uPMatrix");
    this.shaderProgramPerfectSpheres.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "uMVMatrix");
    this.shaderProgramPerfectSpheres.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "uMVINVMatrix");

    this.shaderProgramPerfectSpheres.fog_start = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "fog_start");
    this.shaderProgramPerfectSpheres.fog_end = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "fog_end");
    this.shaderProgramPerfectSpheres.fogColour = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "fogColour");

    this.shaderProgramPerfectSpheres.offset = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "offset");
    this.shaderProgramPerfectSpheres.size = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "size");
    this.shaderProgramPerfectSpheres.scaleMatrix = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "scaleMatrix");
    
    this.shaderProgramPerfectSpheres.light_positions = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "light_positions");
    this.shaderProgramPerfectSpheres.light_colours_ambient = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "light_colours_ambient");
    this.shaderProgramPerfectSpheres.light_colours_specular = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "light_colours_specular");
    this.shaderProgramPerfectSpheres.light_colours_diffuse = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "light_colours_diffuse");
    
    this.shaderProgramPerfectSpheres.clipPlane0 = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "clipPlane0");
    this.shaderProgramPerfectSpheres.clipPlane1 = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "clipPlane1");
    this.shaderProgramPerfectSpheres.clipPlane2 = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "clipPlane2");
    this.shaderProgramPerfectSpheres.clipPlane3 = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "clipPlane3");
    this.shaderProgramPerfectSpheres.clipPlane4 = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "clipPlane4");
    this.shaderProgramPerfectSpheres.clipPlane5 = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "clipPlane5");
    this.shaderProgramPerfectSpheres.clipPlane6 = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "clipPlane6");
    this.shaderProgramPerfectSpheres.clipPlane7 = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "clipPlane7");
    this.shaderProgramPerfectSpheres.nClipPlanes = this.gl.getUniformLocation(this.shaderProgramPerfectSpheres, "nClipPlanes");
}

MGWebGL.prototype.initImageShaders = function(vertexShader,fragmentShader) {
    this.shaderProgramImages = this.gl.createProgram();
    this.gl.attachShader(this.shaderProgramImages, vertexShader);
    this.gl.attachShader(this.shaderProgramImages, fragmentShader);
    this.gl.bindAttribLocation(this.shaderProgramImages, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramImages, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramImages, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramImages, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramImages);

    if (!this.gl.getProgramParameter(this.shaderProgramImages, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initImageShaders)");
    }

    this.gl.useProgram(this.shaderProgramImages);

    this.shaderProgramImages.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramImages, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramImages.vertexPositionAttribute);

    this.shaderProgramImages.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramImages, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramImages.vertexNormalAttribute);

    this.shaderProgramImages.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramImages, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramImages.vertexColourAttribute);

    this.shaderProgramImages.vertexTextureAttribute = this.gl.getAttribLocation(this.shaderProgramImages, "aVertexTexture");

    this.shaderProgramImages.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramImages, "uPMatrix");
    this.shaderProgramImages.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramImages, "uMVMatrix");
    this.shaderProgramImages.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramImages, "uMVINVMatrix");

    this.shaderProgramImages.fog_start = this.gl.getUniformLocation(this.shaderProgramImages, "fog_start");
    this.shaderProgramImages.fog_end = this.gl.getUniformLocation(this.shaderProgramImages, "fog_end");
    this.shaderProgramImages.fogColour = this.gl.getUniformLocation(this.shaderProgramImages, "fogColour");

    this.shaderProgramImages.offset = this.gl.getUniformLocation(this.shaderProgramImages, "offset");
    this.shaderProgramImages.size = this.gl.getUniformLocation(this.shaderProgramImages, "size");
    this.shaderProgramImages.scaleMatrix = this.gl.getUniformLocation(this.shaderProgramImages, "scaleMatrix");
    
}

MGWebGL.prototype.initTwoDShapesShaders = function(vertexShader,fragmentShader) {
    this.shaderProgramTwoDShapes = this.gl.createProgram();
    this.gl.attachShader(this.shaderProgramTwoDShapes, vertexShader);
    this.gl.attachShader(this.shaderProgramTwoDShapes, fragmentShader);
    this.gl.bindAttribLocation(this.shaderProgramTwoDShapes, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramTwoDShapes, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramTwoDShapes, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramTwoDShapes, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramTwoDShapes);

    if (!this.gl.getProgramParameter(this.shaderProgramTwoDShapes, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initTwoDShapesShaders)");
    }

    this.gl.useProgram(this.shaderProgramTwoDShapes);

    this.shaderProgramTwoDShapes.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramTwoDShapes, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramTwoDShapes.vertexPositionAttribute);

    this.shaderProgramTwoDShapes.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramTwoDShapes, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramTwoDShapes.vertexNormalAttribute);

    this.shaderProgramTwoDShapes.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramTwoDShapes, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramTwoDShapes.vertexColourAttribute);

    this.shaderProgramTwoDShapes.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTwoDShapes, "uPMatrix");
    this.shaderProgramTwoDShapes.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTwoDShapes, "uMVMatrix");
    this.shaderProgramTwoDShapes.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramTwoDShapes, "uMVINVMatrix");

    this.shaderProgramTwoDShapes.fog_start = this.gl.getUniformLocation(this.shaderProgramTwoDShapes, "fog_start");
    this.shaderProgramTwoDShapes.fog_end = this.gl.getUniformLocation(this.shaderProgramTwoDShapes, "fog_end");
    this.shaderProgramTwoDShapes.fogColour = this.gl.getUniformLocation(this.shaderProgramTwoDShapes, "fogColour");

    this.shaderProgramTwoDShapes.offset = this.gl.getUniformLocation(this.shaderProgramTwoDShapes, "offset");
    this.shaderProgramTwoDShapes.size = this.gl.getUniformLocation(this.shaderProgramTwoDShapes, "size");
    this.shaderProgramTwoDShapes.scaleMatrix = this.gl.getUniformLocation(this.shaderProgramTwoDShapes, "scaleMatrix");
    
}


MGWebGL.prototype.initPointSpheresShadowShaders = function(vertexShader,fragmentShader) {

    this.shaderProgramPointSpheresShadow = this.gl.createProgram();
    this.gl.attachShader(this.shaderProgramPointSpheresShadow, vertexShader);
    this.gl.attachShader(this.shaderProgramPointSpheresShadow, fragmentShader);
    this.gl.bindAttribLocation(this.shaderProgramPointSpheresShadow, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramPointSpheresShadow, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramPointSpheresShadow, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramPointSpheresShadow, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramPointSpheresShadow);

    if (!this.gl.getProgramParameter(this.shaderProgramPointSpheresShadow, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initPointSpheresShadowShaders)");
    }

    this.gl.useProgram(this.shaderProgramPointSpheresShadow);

    this.shaderProgramPointSpheresShadow.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramPointSpheresShadow, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramPointSpheresShadow.vertexPositionAttribute);

    this.shaderProgramPointSpheresShadow.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramPointSpheresShadow, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramPointSpheresShadow.vertexNormalAttribute);

    this.shaderProgramPointSpheresShadow.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramPointSpheresShadow, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramPointSpheresShadow.vertexColourAttribute);

    this.shaderProgramPointSpheresShadow.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "uPMatrix");
    this.shaderProgramPointSpheresShadow.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "uMVMatrix");
    this.shaderProgramPointSpheresShadow.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "uMVINVMatrix");
    this.shaderProgramPointSpheresShadow.textureMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "TextureMatrix");

    this.shaderProgramPointSpheresShadow.fog_start = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "fog_start");
    this.shaderProgramPointSpheresShadow.fog_end = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "fog_end");
    this.shaderProgramPointSpheresShadow.fogColour = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "fogColour");

    this.shaderProgramPointSpheresShadow.offset = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "offset");
    this.shaderProgramPointSpheresShadow.size = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "size");
    this.shaderProgramPointSpheresShadow.scaleMatrix = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "scaleMatrix");
    
    this.shaderProgramPointSpheresShadow.light_positions = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "light_positions");
    this.shaderProgramPointSpheresShadow.light_colours_ambient = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "light_colours_ambient");
    this.shaderProgramPointSpheresShadow.light_colours_specular = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "light_colours_specular");
    this.shaderProgramPointSpheresShadow.light_colours_diffuse = this.gl.getUniformLocation(this.shaderProgramPointSpheresShadow, "light_colours_diffuse");
    
}

MGWebGL.prototype.initPointSpheresShaders = function(vertexShader,fragmentShader) {

    this.shaderProgramPointSpheres = this.gl.createProgram();
    this.gl.attachShader(this.shaderProgramPointSpheres, vertexShader);
    this.gl.attachShader(this.shaderProgramPointSpheres, fragmentShader);
    this.gl.bindAttribLocation(this.shaderProgramPointSpheres, 0, "aVertexPosition");
    this.gl.bindAttribLocation(this.shaderProgramPointSpheres, 1, "aVertexColour");
    this.gl.bindAttribLocation(this.shaderProgramPointSpheres, 2, "aVertexNormal");
    this.gl.bindAttribLocation(this.shaderProgramPointSpheres, 3, "aVertexTexture");
    this.gl.linkProgram(this.shaderProgramPointSpheres);

    if (!this.gl.getProgramParameter(this.shaderProgramPointSpheres, this.gl.LINK_STATUS)) {
        alert("Could not initialise shaders (initPointSpheresShaders)");
    }

    this.gl.useProgram(this.shaderProgramPointSpheres);

    this.shaderProgramPointSpheres.vertexPositionAttribute = this.gl.getAttribLocation(this.shaderProgramPointSpheres, "aVertexPosition");
    this.gl.enableVertexAttribArray(this.shaderProgramPointSpheres.vertexPositionAttribute);

    this.shaderProgramPointSpheres.vertexNormalAttribute = this.gl.getAttribLocation(this.shaderProgramPointSpheres, "aVertexNormal");
    this.gl.enableVertexAttribArray(this.shaderProgramPointSpheres.vertexNormalAttribute);

    this.shaderProgramPointSpheres.vertexColourAttribute = this.gl.getAttribLocation(this.shaderProgramPointSpheres, "aVertexColour");
    this.gl.enableVertexAttribArray(this.shaderProgramPointSpheres.vertexColourAttribute);

    this.shaderProgramPointSpheres.pMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "uPMatrix");
    this.shaderProgramPointSpheres.mvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "uMVMatrix");
    this.shaderProgramPointSpheres.mvInvMatrixUniform = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "uMVINVMatrix");

    this.shaderProgramPointSpheres.fog_start = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "fog_start");
    this.shaderProgramPointSpheres.fog_end = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "fog_end");
    this.shaderProgramPointSpheres.fogColour = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "fogColour");

    this.shaderProgramPointSpheres.offset = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "offset");
    this.shaderProgramPointSpheres.size = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "size");
    this.shaderProgramPointSpheres.scaleMatrix = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "scaleMatrix");
    
    this.shaderProgramPointSpheres.light_positions = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "light_positions");
    this.shaderProgramPointSpheres.light_colours_ambient = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "light_colours_ambient");
    this.shaderProgramPointSpheres.light_colours_specular = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "light_colours_specular");
    this.shaderProgramPointSpheres.light_colours_diffuse = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "light_colours_diffuse");

    this.shaderProgramPointSpheres.clipPlane0 = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "clipPlane0");
    this.shaderProgramPointSpheres.clipPlane1 = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "clipPlane1");
    this.shaderProgramPointSpheres.clipPlane2 = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "clipPlane2");
    this.shaderProgramPointSpheres.clipPlane3 = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "clipPlane3");
    this.shaderProgramPointSpheres.clipPlane4 = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "clipPlane4");
    this.shaderProgramPointSpheres.clipPlane5 = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "clipPlane5");
    this.shaderProgramPointSpheres.clipPlane6 = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "clipPlane6");
    this.shaderProgramPointSpheres.clipPlane7 = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "clipPlane7");
    this.shaderProgramPointSpheres.nClipPlanes = this.gl.getUniformLocation(this.shaderProgramPointSpheres, "nClipPlanes");
}

MGWebGL.prototype.setLightUniforms = function(program) {
    this.gl.uniform4fv(program.light_positions, this.light_positions);
    this.gl.uniform4fv(program.light_colours_ambient, this.light_colours_ambient);
    this.gl.uniform4fv(program.light_colours_specular, this.light_colours_specular);
    this.gl.uniform4fv(program.light_colours_diffuse, this.light_colours_diffuse);
}

MGWebGL.prototype.setMatrixUniforms = function(program) {
    this.gl.uniformMatrix4fv(program.pMatrixUniform, false, this.pMatrix);
    this.gl.uniformMatrix4fv(program.mvMatrixUniform, false, this.mvMatrix);
    this.gl.uniformMatrix4fv(program.mvInvMatrixUniform, false, this.mvInvMatrix);
    this.gl.uniform1f(program.fog_start, this.gl_fog_start);
    this.gl.uniform1f(program.fog_end, this.gl_fog_end);
    this.gl.uniform1i(program.nClipPlanes, this.gl_nClipPlanes);
    this.gl.uniform4fv(program.clipPlane0, this.gl_clipPlane0);
    this.gl.uniform4fv(program.clipPlane1, this.gl_clipPlane1);
    this.gl.uniform4fv(program.clipPlane2, this.gl_clipPlane2);
    this.gl.uniform4fv(program.clipPlane3, this.gl_clipPlane3);
    this.gl.uniform4fv(program.clipPlane4, this.gl_clipPlane4);
    this.gl.uniform4fv(program.clipPlane5, this.gl_clipPlane5);
    this.gl.uniform4fv(program.clipPlane6, this.gl_clipPlane6);
    this.gl.uniform4fv(program.clipPlane7, this.gl_clipPlane7);
    this.gl.uniform4fv(program.fogColour, new Float32Array(this.background_colour));
    if(program.hasOwnProperty("cursorPos")){
      this.gl.uniform2fv(program.cursorPos, this.gl_cursorPos);
    }
}

tetraColours = [ 1,0,0,1, 1,0,0,1, 1,0,0,1, 1,0,0,1 ];
tetraNormals =  [ 1, 1, 1, -1,-1,1,  1,-1,-1, -1,1,-1 ];
tetraVertices = [ 1, 1, 1, -1,-1,1,  1,-1,-1, -1,1,-1 ];
tetraIndices = [ 0,1,2, 0,3,1, 1,3,2, 0,2,3 ];

var X_1_0 = 0.525731112119;
var X_1_1 = 0.000000000000;
var X_1_2 = 0.850650808352;
var X_1_3 = 0.309016994375;
var X_1_4 = 0.500000000000;
var X_1_5 = 0.809016994375;
var X_1_6 = 1.000000000000;


var X_2_0 =  0.525731112119;
var X_2_1 =  0.000000000000;
var X_2_2 =  0.850650808352;
var X_2_3 =  0.433888564553;
var X_2_4 =  0.259891913008;
var X_2_5 =  0.862668480416;
var X_2_6 =  0.273266528913;
var X_2_7 =  0.961938357784;
var X_2_8 =  0.309016994375;
var X_2_9 =  0.500000000000;
var X_2_10 =  0.809016994375;
var X_2_11 =  0.162459848116;
var X_2_12 =  0.262865556060;
var X_2_13 =  0.951056516295;
var X_2_14 =  1.000000000000;
var X_2_15 =  0.160622035640;
var X_2_16 =  0.702046444776;
var X_2_17 =  0.693780477560;
var X_2_18 =  0.587785252292;
var X_2_19 =  0.425325404176;
var X_2_20 =  0.688190960236;

var icosaVertices1 = [
  -X_1_0, X_1_2, X_1_1,
  -X_1_3, X_1_5, X_1_4,
  X_1_1, X_1_6, X_1_1,
  X_1_1, X_1_0, X_1_2,
  X_1_3, X_1_5, X_1_4,
  X_1_0, X_1_2, X_1_1,
  -X_1_5, X_1_4, X_1_3,
  -X_1_2, X_1_1, X_1_0,
  -X_1_4, X_1_3, X_1_5,
  -X_1_4, -X_1_3, X_1_5,
  X_1_1, -X_1_0, X_1_2,
  X_1_1, X_1_1, X_1_6,
  X_1_4, X_1_3, X_1_5,
  X_1_4, -X_1_3, X_1_5,
  X_1_2, X_1_1, X_1_0,
  X_1_5, X_1_4, X_1_3,
  X_1_6, X_1_1, X_1_1,
  X_1_2, X_1_1, -X_1_0,
  X_1_5, X_1_4, -X_1_3,
  X_1_5, -X_1_4, X_1_3,
  X_1_0, -X_1_2, X_1_1,
  X_1_5, -X_1_4, -X_1_3,
  X_1_3, -X_1_5, X_1_4,
  -X_1_3, -X_1_5, X_1_4,
  -X_1_0, -X_1_2, X_1_1,
  X_1_1, -X_1_6, X_1_1,
  -X_1_3, -X_1_5, -X_1_4,
  X_1_1, -X_1_0, -X_1_2,
  X_1_3, -X_1_5, -X_1_4,
  X_1_4, -X_1_3, -X_1_5,
  X_1_1, X_1_1, -X_1_6,
  X_1_1, X_1_0, -X_1_2,
  X_1_4, X_1_3, -X_1_5,
  -X_1_4, -X_1_3, -X_1_5,
  -X_1_2, X_1_1, -X_1_0,
  -X_1_4, X_1_3, -X_1_5,
  -X_1_5, X_1_4, -X_1_3,
  -X_1_3, X_1_5, -X_1_4,
  X_1_3, X_1_5, -X_1_4,
  -X_1_6, X_1_1, X_1_1,
  -X_1_5, -X_1_4, X_1_3,
  -X_1_5, -X_1_4, -X_1_3
];

var icosaVertices2 = [
  -X_2_0, X_2_2, X_2_1,
  -X_2_3, X_2_5, X_2_4,
  -X_2_6, X_2_7, X_2_1,
  -X_2_8, X_2_10, X_2_9,
  -X_2_11, X_2_13, X_2_12,
  X_2_1, X_2_14, X_2_1,
  X_2_1, X_2_0, X_2_2,
  X_2_15, X_2_17, X_2_16,
  -X_2_15, X_2_17, X_2_16,
  X_2_8, X_2_10, X_2_9,
  X_2_1, X_2_2, X_2_0,
  X_2_0, X_2_2, X_2_1,
  X_2_6, X_2_7, X_2_1,
  X_2_3, X_2_5, X_2_4,
  X_2_11, X_2_13, X_2_12,
  -X_2_17, X_2_16, X_2_15,
  -X_2_10, X_2_9, X_2_8,
  -X_2_18, X_2_20, X_2_19,
  -X_2_2, X_2_1, X_2_0,
  -X_2_16, X_2_15, X_2_17,
  -X_2_5, X_2_4, X_2_3,
  -X_2_9, X_2_8, X_2_10,
  -X_2_20, X_2_19, X_2_18,
  -X_2_4, X_2_3, X_2_5,
  -X_2_19, X_2_18, X_2_20,
  -X_2_16, -X_2_15, X_2_17,
  -X_2_9, -X_2_8, X_2_10,
  -X_2_0, X_2_1, X_2_2,
  X_2_1, -X_2_0, X_2_2,
  X_2_1, -X_2_6, X_2_7,
  -X_2_4, -X_2_3, X_2_5,
  X_2_1, X_2_1, X_2_14,
  -X_2_12, -X_2_11, X_2_13,
  X_2_1, X_2_6, X_2_7,
  -X_2_12, X_2_11, X_2_13,
  X_2_4, X_2_3, X_2_5,
  X_2_12, X_2_11, X_2_13,
  X_2_9, X_2_8, X_2_10,
  X_2_4, -X_2_3, X_2_5,
  X_2_9, -X_2_8, X_2_10,
  X_2_12, -X_2_11, X_2_13,
  X_2_2, X_2_1, X_2_0,
  X_2_16, X_2_15, X_2_17,
  X_2_16, -X_2_15, X_2_17,
  X_2_0, X_2_1, X_2_2,
  X_2_19, X_2_18, X_2_20,
  X_2_5, X_2_4, X_2_3,
  X_2_10, X_2_9, X_2_8,
  X_2_20, X_2_19, X_2_18,
  X_2_17, X_2_16, X_2_15,
  X_2_18, X_2_20, X_2_19,
  X_2_7, X_2_1, X_2_6,
  X_2_14, X_2_1, X_2_1,
  X_2_13, X_2_12, X_2_11,
  X_2_2, X_2_1, -X_2_0,
  X_2_5, X_2_4, -X_2_3,
  X_2_7, X_2_1, -X_2_6,
  X_2_10, X_2_9, -X_2_8,
  X_2_13, X_2_12, -X_2_11,
  X_2_17, X_2_16, -X_2_15,
  X_2_2, X_2_0, X_2_1,
  X_2_5, -X_2_4, X_2_3,
  X_2_10, -X_2_9, X_2_8,
  X_2_13, -X_2_12, X_2_11,
  X_2_0, -X_2_2, X_2_1,
  X_2_17, -X_2_16, -X_2_15,
  X_2_17, -X_2_16, X_2_15,
  X_2_10, -X_2_9, -X_2_8,
  X_2_2, -X_2_0, X_2_1,
  X_2_5, -X_2_4, -X_2_3,
  X_2_13, -X_2_12, -X_2_11,
  X_2_15, -X_2_17, X_2_16,
  X_2_8, -X_2_10, X_2_9,
  X_2_19, -X_2_18, X_2_20,
  X_2_3, -X_2_5, X_2_4,
  X_2_18, -X_2_20, X_2_19,
  X_2_20, -X_2_19, X_2_18,
  -X_2_15, -X_2_17, X_2_16,
  -X_2_8, -X_2_10, X_2_9,
  X_2_1, -X_2_2, X_2_0,
  -X_2_0, -X_2_2, X_2_1,
  -X_2_6, -X_2_7, X_2_1,
  -X_2_3, -X_2_5, X_2_4,
  X_2_1, -X_2_14, X_2_1,
  -X_2_11, -X_2_13, X_2_12,
  X_2_6, -X_2_7, X_2_1,
  X_2_11, -X_2_13, X_2_12,
  -X_2_3, -X_2_5, -X_2_4,
  -X_2_8, -X_2_10, -X_2_9,
  -X_2_11, -X_2_13, -X_2_12,
  X_2_1, -X_2_0, -X_2_2,
  X_2_15, -X_2_17, -X_2_16,
  -X_2_15, -X_2_17, -X_2_16,
  X_2_8, -X_2_10, -X_2_9,
  X_2_1, -X_2_2, -X_2_0,
  X_2_3, -X_2_5, -X_2_4,
  X_2_11, -X_2_13, -X_2_12,
  X_2_4, -X_2_3, -X_2_5,
  X_2_9, -X_2_8, -X_2_10,
  X_2_19, -X_2_18, -X_2_20,
  X_2_16, -X_2_15, -X_2_17,
  X_2_20, -X_2_19, -X_2_18,
  X_2_18, -X_2_20, -X_2_19,
  X_2_1, -X_2_6, -X_2_7,
  X_2_1, X_2_1, -X_2_14,
  X_2_12, -X_2_11, -X_2_13,
  X_2_1, X_2_0, -X_2_2,
  X_2_4, X_2_3, -X_2_5,
  X_2_1, X_2_6, -X_2_7,
  X_2_9, X_2_8, -X_2_10,
  X_2_12, X_2_11, -X_2_13,
  X_2_16, X_2_15, -X_2_17,
  X_2_0, X_2_1, -X_2_2,
  -X_2_4, -X_2_3, -X_2_5,
  -X_2_9, -X_2_8, -X_2_10,
  -X_2_12, -X_2_11, -X_2_13,
  -X_2_2, X_2_1, -X_2_0,
  -X_2_16, X_2_15, -X_2_17,
  -X_2_16, -X_2_15, -X_2_17,
  -X_2_9, X_2_8, -X_2_10,
  -X_2_0, X_2_1, -X_2_2,
  -X_2_4, X_2_3, -X_2_5,
  -X_2_12, X_2_11, -X_2_13,
  -X_2_5, X_2_4, -X_2_3,
  -X_2_10, X_2_9, -X_2_8,
  -X_2_20, X_2_19, -X_2_18,
  -X_2_3, X_2_5, -X_2_4,
  -X_2_17, X_2_16, -X_2_15,
  -X_2_8, X_2_10, -X_2_9,
  -X_2_18, X_2_20, -X_2_19,
  -X_2_15, X_2_17, -X_2_16,
  -X_2_19, X_2_18, -X_2_20,
  -X_2_11, X_2_13, -X_2_12,
  X_2_3, X_2_5, -X_2_4,
  X_2_8, X_2_10, -X_2_9,
  X_2_11, X_2_13, -X_2_12,
  X_2_15, X_2_17, -X_2_16,
  X_2_1, X_2_2, -X_2_0,
  X_2_19, X_2_18, -X_2_20,
  X_2_18, X_2_20, -X_2_19,
  X_2_20, X_2_19, -X_2_18,
  -X_2_7, X_2_1, X_2_6,
  -X_2_13, X_2_12, X_2_11,
  -X_2_14, X_2_1, X_2_1,
  -X_2_2, X_2_0, X_2_1,
  -X_2_7, X_2_1, -X_2_6,
  -X_2_13, X_2_12, -X_2_11,
  -X_2_5, -X_2_4, X_2_3,
  -X_2_13, -X_2_12, X_2_11,
  -X_2_10, -X_2_9, X_2_8,
  -X_2_5, -X_2_4, -X_2_3,
  -X_2_10, -X_2_9, -X_2_8,
  -X_2_13, -X_2_12, -X_2_11,
  -X_2_17, -X_2_16, X_2_15,
  -X_2_17, -X_2_16, -X_2_15,
  -X_2_2, -X_2_0, X_2_1,
  -X_2_20, -X_2_19, X_2_18,
  -X_2_18, -X_2_20, X_2_19,
  -X_2_19, -X_2_18, X_2_20,
  -X_2_19, -X_2_18, -X_2_20,
  -X_2_18, -X_2_20, -X_2_19,
  -X_2_20, -X_2_19, -X_2_18
];

var icosaIndices1 = [
  0, 1, 2,
  3, 4, 1,
  5, 2, 4,
  1, 4, 2,
  0, 6, 1,
  7, 8, 6,
  3, 1, 8,
  6, 8, 1,
  7, 9, 8,
  10, 11, 9,
  3, 8, 11,
  9, 11, 8,
  3, 11, 12,
  10, 13, 11,
  14, 12, 13,
  11, 13, 12,
  3, 12, 4,
  14, 15, 12,
  5, 4, 15,
  12, 15, 4,
  14, 16, 15,
  17, 18, 16,
  5, 15, 18,
  16, 18, 15,
  14, 19, 16,
  20, 21, 19,
  17, 16, 21,
  19, 21, 16,
  10, 22, 13,
  20, 19, 22,
  14, 13, 19,
  22, 19, 13,
  10, 23, 22,
  24, 25, 23,
  20, 22, 25,
  23, 25, 22,
  24, 26, 25,
  27, 28, 26,
  20, 25, 28,
  26, 28, 25,
  27, 29, 28,
  17, 21, 29,
  20, 28, 21,
  29, 21, 28,
  27, 30, 29,
  31, 32, 30,
  17, 29, 32,
  30, 32, 29,
  27, 33, 30,
  34, 35, 33,
  31, 30, 35,
  33, 35, 30,
  34, 36, 35,
  0, 37, 36,
  31, 35, 37,
  36, 37, 35,
  0, 2, 37,
  5, 38, 2,
  31, 37, 38,
  2, 38, 37,
  31, 38, 32,
  5, 18, 38,
  17, 32, 18,
  38, 18, 32,
  7, 6, 39,
  0, 36, 6,
  34, 39, 36,
  6, 36, 39,
  7, 39, 40,
  34, 41, 39,
  24, 40, 41,
  39, 41, 40,
  7, 40, 9,
  24, 23, 40,
  10, 9, 23,
  40, 23, 9,
  27, 26, 33,
  24, 41, 26,
  34, 33, 41,
  26, 41, 33
];

var icosaIndices2 = [
  0, 1, 2,
  3, 4, 1,
  5, 2, 4,
  1, 4, 2,
  6, 7, 8,
  9, 10, 7,
  3, 8, 10,
  7, 10, 8,
  11, 12, 13,
  5, 14, 12,
  9, 13, 14,
  12, 14, 13,
  3, 10, 4,
  9, 14, 10,
  5, 4, 14,
  10, 14, 4,
  0, 15, 1,
  16, 17, 15,
  3, 1, 17,
  15, 17, 1,
  18, 19, 20,
  21, 22, 19,
  16, 20, 22,
  19, 22, 20,
  6, 8, 23,
  3, 24, 8,
  21, 23, 24,
  8, 24, 23,
  16, 22, 17,
  21, 24, 22,
  3, 17, 24,
  22, 24, 17,
  18, 25, 19,
  26, 27, 25,
  21, 19, 27,
  25, 27, 19,
  28, 29, 30,
  31, 32, 29,
  26, 30, 32,
  29, 32, 30,
  6, 23, 33,
  21, 34, 23,
  31, 33, 34,
  23, 34, 33,
  26, 32, 27,
  31, 34, 32,
  21, 27, 34,
  32, 34, 27,
  6, 33, 35,
  31, 36, 33,
  37, 35, 36,
  33, 36, 35,
  28, 38, 29,
  39, 40, 38,
  31, 29, 40,
  38, 40, 29,
  41, 42, 43,
  37, 44, 42,
  39, 43, 44,
  42, 44, 43,
  31, 40, 36,
  39, 44, 40,
  37, 36, 44,
  40, 44, 36,
  6, 35, 7,
  37, 45, 35,
  9, 7, 45,
  35, 45, 7,
  41, 46, 42,
  47, 48, 46,
  37, 42, 48,
  46, 48, 42,
  11, 13, 49,
  9, 50, 13,
  47, 49, 50,
  13, 50, 49,
  37, 48, 45,
  47, 50, 48,
  9, 45, 50,
  48, 50, 45,
  41, 51, 46,
  52, 53, 51,
  47, 46, 53,
  51, 53, 46,
  54, 55, 56,
  57, 58, 55,
  52, 56, 58,
  55, 58, 56,
  11, 49, 59,
  47, 60, 49,
  57, 59, 60,
  49, 60, 59,
  52, 58, 53,
  57, 60, 58,
  47, 53, 60,
  58, 60, 53,
  41, 61, 51,
  62, 63, 61,
  52, 51, 63,
  61, 63, 51,
  64, 65, 66,
  67, 68, 65,
  62, 66, 68,
  65, 68, 66,
  54, 56, 69,
  52, 70, 56,
  67, 69, 70,
  56, 70, 69,
  62, 68, 63,
  67, 70, 68,
  52, 63, 70,
  68, 70, 63,
  28, 71, 38,
  72, 73, 71,
  39, 38, 73,
  71, 73, 38,
  64, 66, 74,
  62, 75, 66,
  72, 74, 75,
  66, 75, 74,
  41, 43, 61,
  39, 76, 43,
  62, 61, 76,
  43, 76, 61,
  72, 75, 73,
  62, 76, 75,
  39, 73, 76,
  75, 76, 73,
  28, 77, 71,
  78, 79, 77,
  72, 71, 79,
  77, 79, 71,
  80, 81, 82,
  83, 84, 81,
  78, 82, 84,
  81, 84, 82,
  64, 74, 85,
  72, 86, 74,
  83, 85, 86,
  74, 86, 85,
  78, 84, 79,
  83, 86, 84,
  72, 79, 86,
  84, 86, 79,
  80, 87, 81,
  88, 89, 87,
  83, 81, 89,
  87, 89, 81,
  90, 91, 92,
  93, 94, 91,
  88, 92, 94,
  91, 94, 92,
  64, 85, 95,
  83, 96, 85,
  93, 95, 96,
  85, 96, 95,
  88, 94, 89,
  93, 96, 94,
  83, 89, 96,
  94, 96, 89,
  90, 97, 91,
  98, 99, 97,
  93, 91, 99,
  97, 99, 91,
  54, 69, 100,
  67, 101, 69,
  98, 100, 101,
  69, 101, 100,
  64, 95, 65,
  93, 102, 95,
  67, 65, 102,
  95, 102, 65,
  98, 101, 99,
  67, 102, 101,
  93, 99, 102,
  101, 102, 99,
  90, 103, 97,
  104, 105, 103,
  98, 97, 105,
  103, 105, 97,
  106, 107, 108,
  109, 110, 107,
  104, 108, 110,
  107, 110, 108,
  54, 100, 111,
  98, 112, 100,
  109, 111, 112,
  100, 112, 111,
  104, 110, 105,
  109, 112, 110,
  98, 105, 112,
  110, 112, 105,
  90, 113, 103,
  114, 115, 113,
  104, 103, 115,
  113, 115, 103,
  116, 117, 118,
  119, 120, 117,
  114, 118, 120,
  117, 120, 118,
  106, 108, 121,
  104, 122, 108,
  119, 121, 122,
  108, 122, 121,
  114, 120, 115,
  119, 122, 120,
  104, 115, 122,
  120, 122, 115,
  116, 123, 117,
  124, 125, 123,
  119, 117, 125,
  123, 125, 117,
  0, 126, 127,
  128, 129, 126,
  124, 127, 129,
  126, 129, 127,
  106, 121, 130,
  119, 131, 121,
  128, 130, 131,
  121, 131, 130,
  124, 129, 125,
  128, 131, 129,
  119, 125, 131,
  129, 131, 125,
  0, 2, 126,
  5, 132, 2,
  128, 126, 132,
  2, 132, 126,
  11, 133, 12,
  134, 135, 133,
  5, 12, 135,
  133, 135, 12,
  106, 130, 136,
  128, 137, 130,
  134, 136, 137,
  130, 137, 136,
  5, 135, 132,
  134, 137, 135,
  128, 132, 137,
  135, 137, 132,
  106, 136, 107,
  134, 138, 136,
  109, 107, 138,
  136, 138, 107,
  11, 59, 133,
  57, 139, 59,
  134, 133, 139,
  59, 139, 133,
  54, 111, 55,
  109, 140, 111,
  57, 55, 140,
  111, 140, 55,
  134, 139, 138,
  57, 140, 139,
  109, 138, 140,
  139, 140, 138,
  18, 20, 141,
  16, 142, 20,
  143, 141, 142,
  20, 142, 141,
  0, 127, 15,
  124, 144, 127,
  16, 15, 144,
  127, 144, 15,
  116, 145, 123,
  143, 146, 145,
  124, 123, 146,
  145, 146, 123,
  16, 144, 142,
  124, 146, 144,
  143, 142, 146,
  144, 146, 142,
  18, 141, 147,
  143, 148, 141,
  149, 147, 148,
  141, 148, 147,
  116, 150, 145,
  151, 152, 150,
  143, 145, 152,
  150, 152, 145,
  80, 153, 154,
  149, 155, 153,
  151, 154, 155,
  153, 155, 154,
  143, 152, 148,
  151, 155, 152,
  149, 148, 155,
  152, 155, 148,
  18, 147, 25,
  149, 156, 147,
  26, 25, 156,
  147, 156, 25,
  80, 82, 153,
  78, 157, 82,
  149, 153, 157,
  82, 157, 153,
  28, 30, 77,
  26, 158, 30,
  78, 77, 158,
  30, 158, 77,
  149, 157, 156,
  78, 158, 157,
  26, 156, 158,
  157, 158, 156,
  90, 92, 113,
  88, 159, 92,
  114, 113, 159,
  92, 159, 113,
  80, 154, 87,
  151, 160, 154,
  88, 87, 160,
  154, 160, 87,
  116, 118, 150,
  114, 161, 118,
  151, 150, 161,
  118, 161, 150,
  88, 160, 159,
  151, 161, 160,
  114, 159, 161,
  160, 161, 159
];

var icosaphi = (1+Math.sqrt(5))/2.;
var stellatedScale = 1.64;
var icosaX = stellatedScale;
var icosaZ = stellatedScale*icosaphi;

var dodecaphi = (1+Math.sqrt(5))/2.;
var dodecaX = 1;
var dodecaZ = dodecaphi;
var stellatedDodecaVertices = [
    /* dodecahedron vertices */
    dodecaX, dodecaX, dodecaX,
    -dodecaX, dodecaX, dodecaX,
    dodecaX, -dodecaX, dodecaX,
    dodecaX, dodecaX, -dodecaX,
    -dodecaX, -dodecaX, dodecaX,
    dodecaX, -dodecaX, -dodecaX,
    -dodecaX, dodecaX, -dodecaX,
    -dodecaX, -dodecaX, -dodecaX,
    0, 1./dodecaZ, dodecaZ,
    0, 1./-dodecaZ, dodecaZ,
    0, 1./dodecaZ, -dodecaZ,
    0, 1./-dodecaZ, -dodecaZ,
    1./dodecaZ, dodecaZ, 0,
    1./-dodecaZ, dodecaZ, 0,
    1./dodecaZ, -dodecaZ, 0,
    1./-dodecaZ, -dodecaZ, 0,
    dodecaZ, 0, 1./dodecaZ,
    -dodecaZ, 0, 1./dodecaZ,
    dodecaZ, 0, 1./-dodecaZ,
    -dodecaZ, 0, 1./-dodecaZ,
    /* icosahedron vertices */
    -icosaX, 0.0, icosaZ,  //20 DONE
    icosaX, 0.0, icosaZ,   //21 DONE
    -icosaX, 0.0, -icosaZ, //22 DONE
    icosaX, 0.0, -icosaZ,  //23 DONE
    0.0, icosaZ, icosaX,   //24 DONE
    0.0, icosaZ, -icosaX,  //25 DONE
    0.0, -icosaZ, icosaX,  //26 DONE
    0.0, -icosaZ, -icosaX, //27 DONE
    icosaZ, icosaX, 0.0,   //28 DONE
    -icosaZ, icosaX, 0.0,  //29 DONE
    icosaZ, -icosaX, 0.0,  //30 DONE
    -icosaZ, -icosaX, 0.0, //31
];

var stellatedDodecaIndices = [
/* Dodecahedron indices */
/*
    0,12,8,
    12,1,8,
    12,13,1,
    13,12,6,
    12,10,6,
    12,3,10,
    3,12,18,
    12,16,18,
    12,0,16,
    16,0,2,
    0,9,2,
    0,8,9,
    18,16,5,
    16,14,5,
    16,2,14,
    10,3,11,
    11,3,5,
    3,18,5,
    9,8,4,
    8,17,4,
    8,1,17,
    17,1,19,
    1,6,19,
    1,13,6,
    19,6,7,
    7,6,11,
    6,10,11,
    7,11,15,
    11,14,15,
    11,5,14,
    15,14,4,
    14,9,4,
    14,2,9,
    4,17,15,
    15,17,7,
    17,19,7,
    */

/* Stellated-Dodecahedron indices */

    3,12,28,
    12,0,28,
    0,16,28,
    16,18,28,
    18,3,28,

    0,12,24,
    12,13,24,
    13,1,24,
    1,8,24,
    8,0,24,

    13,12,25,
    12,3,25,
    3,10,25,
    10,6,25,
    6,13,25,

    16,0,21,
    0,8,21,
    8,9,21,
    9,2,21,
    2,16,21,

    18,16,30,
    16,2,30,
    2,14,30,
    14,5,30,
    5,18,30,

    10,3,23,
    3,18,23,
    18,5,23,
    5,11,23,
    11,10,23,

    9,8,20,
    8,1,20,
    1,17,20,
    17,4,20,
    4,9,20,

    17,1,29,
    1,13,29,
    13,6,29,
    6,19,29,
    19,17,29,

    19,6,22,
    6,10,22,
    10,11,22,
    11,7,22,
    7,19,22,

    7,11,27,
    11,5,27,
    5,14,27,
    14,15,27,
    15,7,27,

    15,14,26,
    14,2,26,
    2,9,26,
    9,4,26,
    4,15,26,

    4,17,31,
    17,19,31,
    19,7,31,
    7,15,31,
    15,4,31,
/*
   // Pentagons....
    0,12,13,1,8,  // DONE
    13,12,3,10,6, // DONE
*   3,12,0,16,18, // DONE
    16,0,8,9,2,   // DONE
    18,16,2,14,5, // DONE
    10,3,18,5,11, // DONE
    9,8,1,17,4,   // DONE
    17,1,13,6,19, // DONE
    19,6,10,11,7, // DONE
    7,11,5,14,15, // DONE
    15,14,2,9,4,  // DONE
    4,17,19,7,15,
    */
];



var icosaVertices = icosaVertices2;
var icosaIndices = icosaIndices2;

var icosaNormals = icosaVertices;

function flatNormalMesh(vertices,indices){
  var newVertices = [];
  var newNormals = [];
  var newIndices = [];
  var idx = 0;

  for(var i=0;i<indices.length;i+=3){
      var i0 = indices[i];
      var i1 = indices[i+1];
      var i2 = indices[i+2];
      var x0 = vertices[3*i0];
      var y0 = vertices[3*i0+1];
      var z0 = vertices[3*i0+2];
      var x1 = vertices[3*i1];
      var y1 = vertices[3*i1+1];
      var z1 = vertices[3*i1+2];
      var x2 = vertices[3*i2];
      var y2 = vertices[3*i2+1];
      var z2 = vertices[3*i2+2];
      //console.log("Indices: "+i0+" "+i1+" "+i2);
      //console.log("p0: "+x0+" "+y0+" "+z0);
      //console.log("p1: "+x1+" "+y1+" "+z1);
      //console.log("p2: "+x2+" "+y2+" "+z2);
      var p0p1 = vec3.create([x1-x0,y1-y0,z1-z0]);
      var p0p2 = vec3.create([x2-x0,y2-y0,z2-z0]);
      vec3.normalize(p0p1);
      vec3.normalize(p0p2);
      var n = vec3.create();
      vec3.cross(p0p1,p0p2,n);
      vec3.normalize(n);
      //console.log("P0P1: "+p0p1[0]+" "+p0p1[1]+" "+p0p1[2]);
      //console.log("P0P2: "+p0p2[0]+" "+p0p2[1]+" "+p0p2[2]);
      //console.log("Normal: "+n[0]+" "+n[1]+" "+n[2]);

      newVertices.push(x0);  newVertices.push(y0);  newVertices.push(z0);
      newVertices.push(x1);  newVertices.push(y1);  newVertices.push(z1);
      newVertices.push(x2);  newVertices.push(y2);  newVertices.push(z2);
      newNormals.push(n[0]); newNormals.push(n[1]); newNormals.push(n[2]);
      newNormals.push(n[0]); newNormals.push(n[1]); newNormals.push(n[2]);
      newNormals.push(n[0]); newNormals.push(n[1]); newNormals.push(n[2]);
      /*
      newNormals.push(1);
      newNormals.push(0);
      newNormals.push(0);
      newNormals.push(0);
      newNormals.push(1);
      newNormals.push(0);
      newNormals.push(0);
      newNormals.push(0);
      newNormals.push(1);
      */
      newIndices.push(idx++);
      newIndices.push(idx++);
      newIndices.push(idx++);

  }
  //console.log(newVertices);
  //console.log(newNormals);
  //console.log(newIndices);

  var ret = {};
  ret["vertices"] = newVertices;
  ret["normals"]  = newNormals;
  ret["indices"]  = newIndices;

  return ret;

}

var starMesh = flatNormalMesh(stellatedDodecaVertices,stellatedDodecaIndices);

//var starNormals = stellatedDodecaVertices;
//var starVertices = stellatedDodecaVertices;
//var starIndices = stellatedDodecaIndices;
var starNormals = starMesh["normals"];
var starVertices = starMesh["vertices"];
var starIndices = starMesh["indices"];

MGWebGL.prototype.buildBuffers = function() {

    var xaxis = vec3.create([1.0,0.0,0.0]);
    var yaxis = vec3.create([0.0,1.0,0.0]);
    var zaxis = vec3.create([0.0,0.0,1.0]);
    var Q = vec3.create();
    var R = vec3.create();
    var cylinder = vec3.create();
    var Q1 = vec3.create();
    var Q2 = vec3.create();
    var Q3 = vec3.create();
    var Q4 = vec3.create();

    // FIXME - These need to be global preferences or properties of primitive.
    // spline_accu = 4 is OK for 5kcr on QtWebKit, 8 runs out of memory.
    var spline_accu = 8;
    var wormWidth = 0.2;
    var smoothArrow = true;
    var accuStep = 20;

    var start = new Date().getTime();
    var thisdisplayBufferslength = this.displayBuffers.length;
    //console.log(thisdisplayBufferslength+" buffers to build");

    for (var idx = 0; idx < thisdisplayBufferslength; idx++){
        if(!this.displayBuffers[idx].isDirty){
            continue;
        }
        for (var j = 0; j < this.displayBuffers[idx].triangleVertexIndexBuffer.length; j++){
            this.displayBuffers[idx].isDirty = false;
            if(this.displayBuffers[idx].bufferTypes[j]==="CYLINDERS"||this.displayBuffers[idx].bufferTypes[j]==="CAPCYLINDERS"){

                // Construct cylinder from original start and end. 
                // FIXME - The original indices are ignored at the moment. We will deal with that in due course.
                //console.log(this.displayBuffers[idx].triangleIndexs[j].length);
                //console.log(this.displayBuffers[idx].triangleVertices[j].length);
                //console.log(this.displayBuffers[idx].triangleNormals[j].length);
                //console.log(this.displayBuffers[idx].triangleColours[j].length);
                var primitiveSizes = this.displayBuffers[idx].primitiveSizes[j];
                //console.log(primitiveSizes.length);
                var triangleIndexs = [];
                var triangleNormals = [];
                var triangleVertices = [];
                var triangleColours = [];
                var newIndex = 0;
                var icol = 0;
                for(var k=0;k<this.displayBuffers[idx].triangleVertices[j].length;k+=6, icol+=8){
                    //triangleVertices.push(this.displayBuffers[idx].triangleVertices[j][k]);
                    var cylinderStart = vec3.create([this.displayBuffers[idx].triangleVertices[j][k],this.displayBuffers[idx].triangleVertices[j][k+1],this.displayBuffers[idx].triangleVertices[j][k+2]]);
                    var cylinderEnd   = vec3.create([this.displayBuffers[idx].triangleVertices[j][k+3],this.displayBuffers[idx].triangleVertices[j][k+4],this.displayBuffers[idx].triangleVertices[j][k+5]]);
                    var colStart =[this.displayBuffers[idx].triangleColours[j][icol],this.displayBuffers[idx].triangleColours[j][icol+1],this.displayBuffers[idx].triangleColours[j][icol+2],this.displayBuffers[idx].triangleColours[j][icol+3]];
                    var colEnd   = [this.displayBuffers[idx].triangleColours[j][icol+4],this.displayBuffers[idx].triangleColours[j][icol+5],this.displayBuffers[idx].triangleColours[j][icol+6],this.displayBuffers[idx].triangleColours[j][icol+7]];
                    vec3.subtract(cylinderEnd,cylinderStart,cylinder);
                    vec3.normalize(cylinder);
                    vec3.cross(xaxis,cylinder,Q);
                    var valid = false;
                    if(vec3.length(Q)>1e-5){
                        vec3.normalize(Q);
                        vec3.cross(cylinder,Q,R);
                        valid = true;
                    } else {
                        vec3.cross(yaxis,cylinder,Q);
                        if(vec3.length(Q)>1e-5){
                            vec3.normalize(Q);
                            vec3.cross(cylinder,Q,R);
                            valid = true;
                        } else {
                            vec3.cross(zaxis,cylinder,Q);
                            if(vec3.length(Q)>1e-5){
                                vec3.normalize(Q);
                                vec3.cross(cylinder,Q,R);
                                valid = true;
                            }
                        }
                    }
                    if(valid){
                        var size = primitiveSizes[icol/8];
                        vec3.scale(Q,size);
                        vec3.scale(R,size);
                        var istep = 0;
                        var offset = triangleVertices.length/3;
                        for(var theta=0;theta<360;theta+=accuStep,istep++){
                            var theta1 = Math.PI *(theta) / 180.0;
                            var c1 = Math.cos(theta1);
                            var s1 = Math.sin(theta1);
                            var p1 = vec3.create([ c1 * Q[0] + s1 * R[0] , c1 * Q[1] + s1 * R[1] , c1 * Q[2] + s1 * R[2] ]);
                            vec3.add(cylinderStart,p1,Q1);
                            vec3.add(cylinderEnd,  p1,Q2);
                            triangleVertices.push(Q1[0]); triangleVertices.push(Q1[1]); triangleVertices.push(Q1[2]);
                            triangleVertices.push(Q2[0]); triangleVertices.push(Q2[1]); triangleVertices.push(Q2[2]);
                            triangleNormals.push(p1[0]);  triangleNormals.push(p1[1]);  triangleNormals.push(p1[2]);
                            triangleNormals.push(p1[0]);  triangleNormals.push(p1[1]);  triangleNormals.push(p1[2]);
                            if(istep>0){
                                triangleIndexs.push(offset+2*(istep-1));
                                triangleIndexs.push(offset+2*(istep));
                                triangleIndexs.push(offset+2*(istep-1)+1);
                                triangleIndexs.push(offset+2*(istep-1)+1);
                                triangleIndexs.push(offset+2*(istep));
                                triangleIndexs.push(offset+2*(istep)+1);
                            }
                            triangleColours.push(colStart[0]); triangleColours.push(colStart[1]); triangleColours.push(colStart[2]); triangleColours.push(colStart[3]);
                            triangleColours.push(colEnd[0]); triangleColours.push(colEnd[1]); triangleColours.push(colEnd[2]); triangleColours.push(colEnd[3]);
                        }
                        // And complete the circle
                        triangleIndexs.push(offset+2*(istep-1));
                        triangleIndexs.push(offset);
                        triangleIndexs.push(offset+2*(istep-1)+1);
                        triangleIndexs.push(offset+2*(istep-1)+1);
                        triangleIndexs.push(offset);
                        triangleIndexs.push(offset+1);

                        if(this.displayBuffers[idx].bufferTypes[j]==="CAPCYLINDERS"){

                            var lastOffset = 1 + triangleIndexs[triangleIndexs.length-3];

                            // Circle cap at start of cylinder.
                            var istep = 0;
                            for(var theta=0;theta<360;theta+=accuStep,istep++){
                                var theta1 = Math.PI *(theta) / 180.0;
                                var theta2 = Math.PI *(theta+accuStep) / 180.0;
                                var c1 = Math.cos(theta1);
                                var s1 = Math.sin(theta1);
                                var c2 = Math.cos(theta2);
                                var s2 = Math.sin(theta2);
                                var p1 = vec3.create([ c1 * Q[0] + s1 * R[0] , c1 * Q[1] + s1 * R[1] , c1 * Q[2] + s1 * R[2] ]);
                                var p2 = vec3.create([ c2 * Q[0] + s2 * R[0] , c2 * Q[1] + s2 * R[1] , c2 * Q[2] + s2 * R[2] ]);
                                vec3.add(cylinderStart,p1,Q1);
                                vec3.add(cylinderStart,p2,Q2);
                                triangleVertices.push(cylinderStart[0]); triangleVertices.push(cylinderStart[1]); triangleVertices.push(cylinderStart[2]);
                                triangleVertices.push(Q1[0]); triangleVertices.push(Q1[1]); triangleVertices.push(Q1[2]);
                                triangleVertices.push(Q2[0]); triangleVertices.push(Q2[1]); triangleVertices.push(Q2[2]);
                                triangleNormals.push(-cylinder[0]);  triangleNormals.push(-cylinder[1]);  triangleNormals.push(-cylinder[2]);
                                triangleNormals.push(-cylinder[0]);  triangleNormals.push(-cylinder[1]);  triangleNormals.push(-cylinder[2]);
                                triangleNormals.push(-cylinder[0]);  triangleNormals.push(-cylinder[1]);  triangleNormals.push(-cylinder[2]);
                                triangleColours.push(colStart[0]); triangleColours.push(colStart[1]); triangleColours.push(colStart[2]); triangleColours.push(colStart[3]);
                                triangleColours.push(colStart[0]); triangleColours.push(colStart[1]); triangleColours.push(colStart[2]); triangleColours.push(colStart[3]);
                                triangleColours.push(colStart[0]); triangleColours.push(colStart[1]); triangleColours.push(colStart[2]); triangleColours.push(colStart[3]);
                                triangleIndexs.push(3*istep+lastOffset);
                                triangleIndexs.push(3*istep+lastOffset+2);
                                triangleIndexs.push(3*istep+lastOffset+1);
                            }

                            // Circle cap at end of cylinder.
                            for(var theta=0;theta<360;theta+=accuStep,istep++){
                                var theta1 = Math.PI *(theta) / 180.0;
                                var theta2 = Math.PI *(theta+accuStep) / 180.0;
                                var c1 = Math.cos(theta1);
                                var s1 = Math.sin(theta1);
                                var c2 = Math.cos(theta2);
                                var s2 = Math.sin(theta2);
                                var p1 = vec3.create([ c1 * Q[0] + s1 * R[0] , c1 * Q[1] + s1 * R[1] , c1 * Q[2] + s1 * R[2] ]);
                                var p2 = vec3.create([ c2 * Q[0] + s2 * R[0] , c2 * Q[1] + s2 * R[1] , c2 * Q[2] + s2 * R[2] ]);
                                vec3.add(cylinderEnd,p1,Q1);
                                vec3.add(cylinderEnd,p2,Q2);
                                triangleVertices.push(cylinderEnd[0]); triangleVertices.push(cylinderEnd[1]); triangleVertices.push(cylinderEnd[2]);
                                triangleVertices.push(Q1[0]); triangleVertices.push(Q1[1]); triangleVertices.push(Q1[2]);
                                triangleVertices.push(Q2[0]); triangleVertices.push(Q2[1]); triangleVertices.push(Q2[2]);
                                triangleNormals.push(-cylinder[0]);  triangleNormals.push(-cylinder[1]);  triangleNormals.push(-cylinder[2]);
                                triangleNormals.push(-cylinder[0]);  triangleNormals.push(-cylinder[1]);  triangleNormals.push(-cylinder[2]);
                                triangleNormals.push(-cylinder[0]);  triangleNormals.push(-cylinder[1]);  triangleNormals.push(-cylinder[2]);
                                triangleColours.push(colEnd[0]); triangleColours.push(colEnd[1]); triangleColours.push(colEnd[2]); triangleColours.push(colEnd[3]);
                                triangleColours.push(colEnd[0]); triangleColours.push(colEnd[1]); triangleColours.push(colEnd[2]); triangleColours.push(colEnd[3]);
                                triangleColours.push(colEnd[0]); triangleColours.push(colEnd[1]); triangleColours.push(colEnd[2]); triangleColours.push(colEnd[3]);
                                triangleIndexs.push(3*istep+lastOffset);
                                triangleIndexs.push(3*istep+lastOffset+2);
                                triangleIndexs.push(3*istep+lastOffset+1);
                            }
                        }

                    } else {
                        // FIXME - Panic. Need to deal with 0 length cylinder
                    }
                }

                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexIndexBuffer[j]);
                if(this.ext){
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(triangleIndexs), this.gl.STATIC_DRAW);
                }else{
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(triangleIndexs), this.gl.STATIC_DRAW);
                }

                this.displayBuffers[idx].triangleVertexNormalBuffer[j].numItems   = triangleNormals.length/3;
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].numItems = triangleNormals.length/3;
                this.displayBuffers[idx].triangleColourBuffer[j].numItems         = triangleColours.length/4;
                this.displayBuffers[idx].triangleVertexIndexBuffer[j].numItems    = triangleIndexs.length;
                //console.log("Buffering "+triangleNormals.length/3+" vertices");
                //console.log("Buffering "+triangleIndexs.length+" indices");

                this.displayBuffers[idx].triangleVertexIndexBuffer[j].itemSize = 1; 
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexNormalBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(triangleNormals), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexNormalBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexPositionBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(triangleVertices), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleColourBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(triangleColours), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleColourBuffer[j].itemSize = 4;

            } else if(this.displayBuffers[idx].bufferTypes[j]==="SPLINE"||this.displayBuffers[idx].bufferTypes[j]==="WORM"){
                var thisSplineAccu = spline_accu;
                var thisAccuStep = accuStep;
                if(typeof(this.displayBuffers[idx].supplementary["spline_accu"])!=="undefined"&&this.displayBuffers[idx].supplementary["spline_accu"][j].length>0){
                    thisSplineAccu = this.displayBuffers[idx].supplementary["spline_accu"][j][0];
                }
                if(typeof(this.displayBuffers[idx].supplementary["accu"])!=="undefined"&&this.displayBuffers[idx].supplementary["accu"][j].length>0){
                    thisAccuStep = this.displayBuffers[idx].supplementary["accu"][j][0];
                }
                var isWorm = false;
                if(this.displayBuffers[idx].bufferTypes[j]==="WORM"){
                    isWorm = true;
                }
                var primitiveSizes = this.displayBuffers[idx].primitiveSizes;
                var normals = [];
                var N = vec3.create();
                var prevN = vec3.create();
                var ca1ca2 = vec3.create();
                var ca1ca3 = vec3.create();

                if(typeof(this.displayBuffers[idx].supplementary["customSplineNormals"])!=="undefined"&&this.displayBuffers[idx].supplementary["customSplineNormals"][j].length>0){
                    var customSplineNormals = this.displayBuffers[idx].supplementary["customSplineNormals"][j];
                    for(var k=0;k<customSplineNormals.length;k++){
                        normals.push(customSplineNormals[k][0]);
                        normals.push(customSplineNormals[k][1]);
                        normals.push(customSplineNormals[k][2]);
                    }
                } else {
                    for(var k=3;k<this.displayBuffers[idx].triangleVertices[j].length-3;k+=3){
                        var x1 = this.displayBuffers[idx].triangleVertices[j][k];
                        var y1 = this.displayBuffers[idx].triangleVertices[j][k+1];
                        var z1 = this.displayBuffers[idx].triangleVertices[j][k+2];
                        var ca1 = vec3.create([x1,y1,z1]);
                        var x2 = this.displayBuffers[idx].triangleVertices[j][k+3];
                        var y2 = this.displayBuffers[idx].triangleVertices[j][k+4];
                        var z2 = this.displayBuffers[idx].triangleVertices[j][k+5];
                        var ca2 = vec3.create([x2,y2,z2]);
                        var x3 = this.displayBuffers[idx].triangleVertices[j][k-3];
                        var y3 = this.displayBuffers[idx].triangleVertices[j][k-2];
                        var z3 = this.displayBuffers[idx].triangleVertices[j][k-1];
                        var ca3 = vec3.create([x3,y3,z3]);
                        vec3.subtract(ca1,ca2,ca1ca2);
                        vec3.subtract(ca1,ca3,ca1ca3);

                        vec3.normalize(ca1ca2);
                        vec3.normalize(ca1ca3);
                        vec3.cross(ca1ca2,ca1ca3,N);
                        vec3.normalize(N);
                        if(k>3){
                            prevN = vec3.create([ normals[normals.length-3],normals[normals.length-2],normals[normals.length-1]]);
                            var isFlip = vec3.dot(N,prevN);
                            if(isFlip<0.0){
                                normals.push(-N[0]);
                                normals.push(-N[1]);
                                normals.push(-N[2]);
                            } else {
                                normals.push(N[0]);
                                normals.push(N[1]);
                                normals.push(N[2]);
                            }
                        } else {
                            normals.push(N[0]);
                            normals.push(N[1]);
                            normals.push(N[2]);
                            normals.push(N[0]);
                            normals.push(N[1]);
                            normals.push(N[2]);
                        }
                    }
                    normals.push(normals[normals.length-3]);
                    normals.push(normals[normals.length-3]);
                    normals.push(normals[normals.length-3]);
                }

                var spline = [];
                var nspline = [];
                if(thisSplineAccu>1){
                    spline = SplineCurve(this.displayBuffers[idx].triangleVertices[j],thisSplineAccu,3,1);
                    nspline = SplineCurve(normals,thisSplineAccu,3,1);
                } else {
                    spline = this.displayBuffers[idx].triangleVertices[j];
                    nspline = normals;
                }

                var sheetEnds = [];

                /*
                for(var k=0;k<spline.length;k++){
                        console.log("old "+spline[k]+" "+k);
                }
                */

                if(typeof(this.displayBuffers[idx].supplementary["arrow"])!=="undefined"){
                    var bez = this.displayBuffers[idx].supplementary["arrow"][j];
                    for(var ibez=0;ibez<bez.length;ibez++){
                        //console.log(bez[ibez]);
                        var sliceInStart = (bez[ibez][0]-1)*3;
                        if(sliceInStart<0) sliceInStart=0;
                        sheetEnds.push(bez[ibez][1]*3);
                        if(smoothArrow){
                            var sliceInEnd   = (bez[ibez][1]+2)*3;
                            if(sliceInEnd*thisSplineAccu>spline.length) sliceInEnd = spline.length/thisSplineAccu;
                            var sliceIn      = this.displayBuffers[idx].triangleVertices[j].slice(sliceInStart,sliceInEnd);
                            var nsliceIn     = normals.slice(sliceInStart,sliceInEnd);
                            var flat_spline   = BezierCurve(sliceIn,thisSplineAccu);
                            var flat_nspline  = BezierCurve(nsliceIn,thisSplineAccu);

                            var sliceOutStart  = 0;
                            var sliceOutLength = flat_spline.length;
                            var sliceOutEnd    = flat_spline.length;
                            if(sliceInStart>6){
                                sliceOutStart  += 3*thisSplineAccu;
                                sliceOutLength -= 3*thisSplineAccu;
                                for(var isplice=0;isplice<thisSplineAccu;isplice++){
                                    var frac = isplice*1.0/(1.0*thisSplineAccu-1);
                                    //console.log("old "+spline[sliceInStart*thisSplineAccu+3*isplice]+" "+(sliceInStart*thisSplineAccu+3*isplice));
                                    //console.log("old "+spline[sliceInStart*thisSplineAccu+3*isplice+1]+" "+(sliceInStart*thisSplineAccu+3*isplice+1));
                                    //console.log("old "+spline[sliceInStart*thisSplineAccu+3*isplice+2]+" "+(sliceInStart*thisSplineAccu+3*isplice+2));
                                    spline[sliceInStart*thisSplineAccu+3*isplice]   = (1.0-frac)*spline[sliceInStart*thisSplineAccu+3*isplice]   + frac*flat_spline[3*isplice];
                                    spline[sliceInStart*thisSplineAccu+3*isplice+1] = (1.0-frac)*spline[sliceInStart*thisSplineAccu+3*isplice+1] + frac*flat_spline[3*isplice+1];
                                    spline[sliceInStart*thisSplineAccu+3*isplice+2] = (1.0-frac)*spline[sliceInStart*thisSplineAccu+3*isplice+2] + frac*flat_spline[3*isplice+2];
                                    nspline[sliceInStart*thisSplineAccu+3*isplice]   = (1.0-frac)*nspline[sliceInStart*thisSplineAccu+3*isplice]   + frac*flat_nspline[3*isplice];
                                    nspline[sliceInStart*thisSplineAccu+3*isplice+1] = (1.0-frac)*nspline[sliceInStart*thisSplineAccu+3*isplice+1] + frac*flat_nspline[3*isplice+1];
                                    nspline[sliceInStart*thisSplineAccu+3*isplice+2] = (1.0-frac)*nspline[sliceInStart*thisSplineAccu+3*isplice+2] + frac*flat_nspline[3*isplice+2];
                                    //console.log("new "+spline[sliceInStart*thisSplineAccu+3*isplice]+" "+(sliceInStart*thisSplineAccu+3*isplice));
                                    //console.log("new "+spline[sliceInStart*thisSplineAccu+3*isplice+1]+" "+(sliceInStart*thisSplineAccu+3*isplice+1));
                                    //console.log("new "+spline[sliceInStart*thisSplineAccu+3*isplice+2]+" "+(sliceInStart*thisSplineAccu+3*isplice+2));
                                }
                            }
                            if(spline.length-sliceInEnd*thisSplineAccu>6){
                                sliceOutLength -= 3*thisSplineAccu;
                                sliceOutEnd    -= 3*thisSplineAccu;
                                for(var isplice=0;isplice<thisSplineAccu;isplice++){
                                    var frac = isplice*1.0/(1.0*thisSplineAccu-1);
                                    //console.log("old "+spline[(sliceInEnd-3)*thisSplineAccu+3*isplice]+" "+((sliceInEnd-3)*thisSplineAccu+3*isplice));
                                    //console.log("old "+spline[(sliceInEnd-3)*thisSplineAccu+3*isplice+1]+" "+((sliceInEnd-3)*thisSplineAccu+3*isplice+1));
                                    //console.log("old "+spline[(sliceInEnd-3)*thisSplineAccu+3*isplice+2]+" "+((sliceInEnd-3)*thisSplineAccu+3*isplice+2));
                                    spline[(sliceInEnd-3)*thisSplineAccu+3*isplice]   = frac*spline[(sliceInEnd-3)*thisSplineAccu+3*isplice]   + (1.0-frac)*flat_spline[sliceOutEnd+3*isplice];
                                    spline[(sliceInEnd-3)*thisSplineAccu+3*isplice+1] = frac*spline[(sliceInEnd-3)*thisSplineAccu+3*isplice+1] + (1.0-frac)*flat_spline[sliceOutEnd+3*isplice+1];
                                    spline[(sliceInEnd-3)*thisSplineAccu+3*isplice+2] = frac*spline[(sliceInEnd-3)*thisSplineAccu+3*isplice+2] + (1.0-frac)*flat_spline[sliceOutEnd+3*isplice+2];
                                    nspline[(sliceInEnd-3)*thisSplineAccu+3*isplice]   = frac*nspline[(sliceInEnd-3)*thisSplineAccu+3*isplice]   + (1.0-frac)*flat_nspline[sliceOutEnd+3*isplice];
                                    nspline[(sliceInEnd-3)*thisSplineAccu+3*isplice+1] = frac*nspline[(sliceInEnd-3)*thisSplineAccu+3*isplice+1] + (1.0-frac)*flat_nspline[sliceOutEnd+3*isplice+1];
                                    nspline[(sliceInEnd-3)*thisSplineAccu+3*isplice+2] = frac*nspline[(sliceInEnd-3)*thisSplineAccu+3*isplice+2] + (1.0-frac)*flat_nspline[sliceOutEnd+3*isplice+2];
                                    //console.log("new "+spline[(sliceInEnd-3)*thisSplineAccu+3*isplice]+" "+((sliceInEnd-3)*thisSplineAccu+3*isplice));
                                    //console.log("new "+spline[(sliceInEnd-3)*thisSplineAccu+3*isplice+1]+" "+((sliceInEnd-3)*thisSplineAccu+3*isplice+1));
                                    //console.log("new "+spline[(sliceInEnd-3)*thisSplineAccu+3*isplice+2]+" "+((sliceInEnd-3)*thisSplineAccu+3*isplice+2));
                                }
                            }
                            Array.prototype.splice.apply(spline,[sliceInStart*thisSplineAccu+sliceOutStart,sliceOutLength].concat(flat_spline.slice(sliceOutStart,sliceOutEnd)));
                            Array.prototype.splice.apply(nspline,[sliceInStart*thisSplineAccu+sliceOutStart,sliceOutLength].concat(flat_nspline.slice(sliceOutStart,sliceOutEnd)));
                        }
                    }
                }

                /*
                for(var k=0;k<spline.length;k++){
                        console.log("new "+spline[k]+" "+k);
                }
                */

                //console.log(sheetEnds);
                //console.log(smoothArrow);

                var pdiff = vec3.create();
                var n2 = vec3.create();
                var n2spline = [];

                for(var k=0;k<nspline.length;k+=3){
                    var Nx = nspline[k];
                    var Ny = nspline[k+1];
                    var Nz = nspline[k+2];
                    var N = vec3.create([Nx,Ny,Nz]);
                    vec3.normalize(N);
                    nspline[k]   = N[0];
                    nspline[k+1] = N[1];
                    nspline[k+2] = N[2];
                    if(k<nspline.length-3){
                        var p1x = spline[k];
                        var p1y = spline[k+1];
                        var p1z = spline[k+2];
                        var p1 = vec3.create([p1x,p1y,p1z]);
                        var p2x = spline[k+3];
                        var p2y = spline[k+4];
                        var p2z = spline[k+5];
                        var p2 = vec3.create([p2x,p2y,p2z]);
                        vec3.subtract(p1,p2,pdiff);
                        vec3.cross(N,pdiff,n2);
                        //console.log(vec3.length(n2)+" "+k);
                        vec3.normalize(n2);
                        n2spline.push(n2[0]);
                        n2spline.push(n2[1]);
                        n2spline.push(n2[2]);
                    }
                }
                n2spline.push(n2spline[n2spline.length-3]);
                n2spline.push(n2spline[n2spline.length-3]);
                n2spline.push(n2spline[n2spline.length-3]);

                var triangleIndexs = [];
                var triangleNormals = [];
                var triangleVertices = [];
                var triangleColours = [];

                var scale = [];
                for(var k=0;k<spline.length;k+=3){
                    var doneThis = false;
                    for(var iarrow=0;iarrow<thisSplineAccu;iarrow++){
                        if(sheetEnds.indexOf((k-3*iarrow)/thisSplineAccu)>-1){
                            var thisScale = (1.0-((1.0*iarrow)/(thisSplineAccu-1.0)))*1.5+wormWidth;
                            scale.push(thisScale);
                            scale.push(thisScale);
                            scale.push(thisScale);
                            doneThis = true;
                            break;
                        }
                    }
                    if(!doneThis){
                        scale.push(1.0);
                        scale.push(1.0);
                        scale.push(1.0);
                    }
                }
                //console.log(spline.length+" "+scale.length);

                var index = 0;
                for(var theta=0;theta<360;theta+=thisAccuStep,istep++){
                    var theta1 = Math.PI *(theta) / 180.0;
                    var theta2 = Math.PI *(theta+thisAccuStep) / 180.0;
                    var c1 = Math.cos(theta1);
                    var s1 = Math.sin(theta1);
                    var c2 = Math.cos(theta2);
                    var s2 = Math.sin(theta2);
                    if(index>0){
                        triangleIndexs.push(index);
                        triangleIndexs.push(index);
                    }
                    var sizeIndex = 0;
                    //console.log(spline.length/thisSplineAccu+" "+this.displayBuffers[idx].triangleColours[j].length);
                    for(var k=0;k<spline.length;k+=3){
                        var colIndex = parseInt(k/3/thisSplineAccu)*4;
                        var colr = this.displayBuffers[idx].triangleColours[j][colIndex];
                        var colg = this.displayBuffers[idx].triangleColours[j][colIndex+1];
                        var colb = this.displayBuffers[idx].triangleColours[j][colIndex+2];
                        var cola = this.displayBuffers[idx].triangleColours[j][colIndex+3];
                        //console.log(colr+" "+colg+" "+colb+" "+cola+" "+parseInt(k*4/3/thisSplineAccu)+" "+this.displayBuffers[idx].triangleColours[j].length);
                        //console.log(parseInt(colIndex/3)*3+" "+this.displayBuffers[idx].triangleColours[j].length);
                        var thisSize2 = primitiveSizes[j][parseInt(sizeIndex/thisSplineAccu)]* scale[k];
                        var thisSize1 = wormWidth;
                        if(isWorm){
                            thisSize1 = primitiveSizes[j][parseInt(sizeIndex/thisSplineAccu)]* scale[k];
                        }
                        sizeIndex++;
                        var n1p = [ c1 * thisSize2*nspline[k] + s1 * thisSize1*n2spline[k] , c1 * thisSize2*nspline[k+1] + s1 * thisSize1*n2spline[k+1] , c1 * thisSize2*nspline[k+2] + s1 * thisSize1*n2spline[k+2] ];
                        var n2p = [ c2 * thisSize2*nspline[k] + s2 * thisSize1*n2spline[k] , c2 * thisSize2*nspline[k+1] + s2 * thisSize1*n2spline[k+1] , c2 * thisSize2*nspline[k+2] + s2 * thisSize1*n2spline[k+2] ];
                        var n1 = [ c1 * nspline[k] / thisSize2 + s1 * n2spline[k] / thisSize1, c1 * nspline[k+1] / thisSize2 + s1 * n2spline[k+1] / thisSize1 , c1 * nspline[k+2] / thisSize2 + s1 * n2spline[k+2]  / thisSize1];
                        var n2 = [ c2 * nspline[k] / thisSize2 + s2 * n2spline[k]  / thisSize1, c2 * nspline[k+1] / thisSize2 + s2 * n2spline[k+1] / thisSize1 , c2 * nspline[k+2] / thisSize2 + s2 * n2spline[k+2]  / thisSize1];
                        var p1 = [ spline[k] + n1p[0], spline[k+1] + n1p[1], spline[k+2] + n1p[2] ];
                        var p2 = [ spline[k] + n2p[0], spline[k+1] + n2p[1], spline[k+2] + n2p[2] ];
                        triangleNormals.push(n1[0]);
                        triangleNormals.push(n1[1]);
                        triangleNormals.push(n1[2]);
                        triangleNormals.push(n2[0]);
                        triangleNormals.push(n2[1]);
                        triangleNormals.push(n2[2]);
                        triangleVertices.push(p1[0]);
                        triangleVertices.push(p1[1]);
                        triangleVertices.push(p1[2]);
                        triangleVertices.push(p2[0]);
                        triangleVertices.push(p2[1]);
                        triangleVertices.push(p2[2]);
                        triangleColours.push(colr);
                        triangleColours.push(colg);
                        triangleColours.push(colb);
                        triangleColours.push(cola);
                        triangleColours.push(colr);
                        triangleColours.push(colg);
                        triangleColours.push(colb);
                        triangleColours.push(cola);
                        triangleIndexs.push(index); index++;
                        triangleIndexs.push(index); index++;
                    }
                    triangleIndexs.push(index-1);
                    triangleIndexs.push(index-1);
                }

                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexIndexBuffer[j]);
                if(this.ext){
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(triangleIndexs), this.gl.STATIC_DRAW);
                }else{
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(triangleIndexs), this.gl.STATIC_DRAW);
                }

                this.displayBuffers[idx].triangleVertexNormalBuffer[j].numItems   = triangleNormals.length/3;
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].numItems = triangleVertices.length/3;
                this.displayBuffers[idx].triangleColourBuffer[j].numItems         = triangleColours.length/4;
                this.displayBuffers[idx].triangleVertexIndexBuffer[j].numItems    = triangleIndexs.length;
                //console.log("Buffering "+triangleVertices.length/3+" vertices");
                //console.log("Buffering "+triangleNormals.length/3+" normals");
                //console.log("Buffering "+triangleColours.length/4+" colours");
                //console.log("Buffering "+triangleIndexs.length+" indices");

                this.displayBuffers[idx].triangleVertexIndexBuffer[j].itemSize = 1; 
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexNormalBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(triangleNormals), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexNormalBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexPositionBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(triangleVertices), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleColourBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(triangleColours), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleColourBuffer[j].itemSize = 4;
                

            } else if(this.displayBuffers[idx].bufferTypes[j]==="STARS"){
                if(typeof(this.starBuffer) === "undefined"){
                    this.starBuffer = new DisplayBuffer();
                    this.starBuffer.triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.starBuffer.triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.starBuffer.triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.starBuffer.triangleVertexIndexBuffer[0]);
                    this.starBuffer.triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.starBuffer.triangleVertexIndexBuffer[0].numItems = starIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(starIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(starIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.starBuffer.triangleVertexNormalBuffer[0]);
                    this.starBuffer.triangleVertexNormalBuffer[0].itemSize = 3;
                    this.starBuffer.triangleVertexNormalBuffer[0].numItems = starNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(starNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.starBuffer.triangleVertexPositionBuffer[0]);
                    this.starBuffer.triangleVertexPositionBuffer[0].itemSize = 3;
                    this.starBuffer.triangleVertexPositionBuffer[0].numItems = starVertices.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(starVertices), this.gl.STATIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j]==="POINTS_SPHERES"||this.displayBuffers[idx].bufferTypes[j]==="SPHEROIDS"){
                if(typeof(this.sphereBuffer) === "undefined"){
                    this.sphereBuffer = new DisplayBuffer();
                    this.sphereBuffer.triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.sphereBuffer.triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.sphereBuffer.triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.sphereBuffer.triangleVertexIndexBuffer[0]);
                    this.sphereBuffer.triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.sphereBuffer.triangleVertexIndexBuffer[0].numItems = icosaIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(icosaIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(icosaIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.sphereBuffer.triangleVertexNormalBuffer[0]);
                    this.sphereBuffer.triangleVertexNormalBuffer[0].itemSize = 3;
                    this.sphereBuffer.triangleVertexNormalBuffer[0].numItems = icosaNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(icosaNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.sphereBuffer.triangleVertexPositionBuffer[0]);
                    this.sphereBuffer.triangleVertexPositionBuffer[0].itemSize = 3;
                    this.sphereBuffer.triangleVertexPositionBuffer[0].numItems = icosaVertices.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(icosaVertices), this.gl.STATIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j].substring(0,8)==="POLYSTAR"){
                if(typeof(this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]]) === "undefined"){
                    var npoints = parseInt(this.displayBuffers[idx].bufferTypes[j].substring(8));
                    var diskIndices = [];
                    var diskNormals = [];
                    this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]] = [];
                    var accuStep = 360./npoints;
                    var diskIdx = 0;
                    this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                    this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                    this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(-1.0);
                    diskIndices.push(diskIdx++);
                    var theta = 0.0;
                    var p = npoints;
                    var q = 2;
                    var S = 1./(Math.cos(Math.PI/p));
                    var R = Math.sin((p-2*q)/(2*p)*Math.PI)/Math.sin((2*q)/p*Math.PI)/S;
                    for(var istep=0;istep<=npoints;istep++){
                        var theta1 = Math.PI *(theta) / 180.0;
                        var y1 = Math.cos(theta1)*R;
                        var x1 = Math.sin(theta1)*R;
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(x1);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(y1);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        diskIndices.push(diskIdx++);
                        theta += accuStep*.5;
                        theta1 = Math.PI *(theta) / 180.0;
                        // FIXME - Scaling is the tricky bit.
                        y1 = Math.cos(theta1);
                        x1 = Math.sin(theta1);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(x1);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(y1);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        diskIndices.push(diskIdx++);
                        theta += accuStep*.5;
                    }
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]] = new DisplayBuffer();
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer[0]);
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer[0].numItems = diskIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(diskIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(diskIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer[0]);
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer[0].itemSize = 3;
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer[0].numItems = diskNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer[0]);
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer[0].itemSize = 3;
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer[0].numItems = this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]]), this.gl.DYNAMIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j].substring(0,7)==="POLYGON"){
                if(typeof(this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]]) === "undefined"){
                    var npoints = parseInt(this.displayBuffers[idx].bufferTypes[j].substring(7));
                    var diskIndices = [];
                    var diskNormals = [];
                    this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]] = [];
                    var accuStep = 360./npoints;
                    var diskIdx = 0;
                    this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                    this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                    this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(-1.0);
                    diskIndices.push(diskIdx++);
                    var theta = 0.0;
                    for(var istep=0;istep<=npoints;istep++){
                        var theta1 = Math.PI *(theta) / 180.0;
                        var y1 = Math.cos(theta1);
                        var x1 = Math.sin(theta1);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(x1);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(y1);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        diskIndices.push(diskIdx++);
                        theta += accuStep;
                    }
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]] = new DisplayBuffer();
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer[0]);
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer[0].numItems = diskIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(diskIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(diskIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer[0]);
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer[0].itemSize = 3;
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer[0].numItems = diskNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer[0]);
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer[0].itemSize = 3;
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer[0].numItems = this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]]), this.gl.DYNAMIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j]==="DIAMONDS"){
                if(typeof(this.diamondBuffer) === "undefined"){
                    var diskIndices = [];
                    var diskNormals = [];
                    this.diamondVertices = [];
                    var accuStep = 90;
                    var diskIdx = 0;
                    this.diamondVertices.push(0.0);
                    this.diamondVertices.push(0.0);
                    this.diamondVertices.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(-1.0);
                    diskIndices.push(diskIdx++);
                    for(var theta=0;theta<=360;theta+=accuStep){
                        var theta1 = Math.PI *(theta) / 180.0;
                        var y1 = Math.cos(theta1);
                        var x1 = Math.sin(theta1);
                        this.diamondVertices.push(x1);
                        this.diamondVertices.push(y1);
                        this.diamondVertices.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        diskIndices.push(diskIdx++);
                    }
                    this.diamondBuffer = new DisplayBuffer();
                    this.diamondBuffer.triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.diamondBuffer.triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.diamondBuffer.triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.diamondBuffer.triangleVertexIndexBuffer[0]);
                    this.diamondBuffer.triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.diamondBuffer.triangleVertexIndexBuffer[0].numItems = diskIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(diskIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(diskIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.diamondBuffer.triangleVertexNormalBuffer[0]);
                    this.diamondBuffer.triangleVertexNormalBuffer[0].itemSize = 3;
                    this.diamondBuffer.triangleVertexNormalBuffer[0].numItems = diskNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.diamondBuffer.triangleVertexPositionBuffer[0]);
                    this.diamondBuffer.triangleVertexPositionBuffer[0].itemSize = 3;
                    this.diamondBuffer.triangleVertexPositionBuffer[0].numItems = this.diamondVertices.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.diamondVertices), this.gl.DYNAMIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j]==="PERFECT_SPHERES"||this.displayBuffers[idx].bufferTypes[j]==="IMAGES"||this.displayBuffers[idx].bufferTypes[j]==="TEXT"){
                if(typeof(this.imageBuffer) === "undefined"){
                    var diskIndices = [];
                    var diskNormals = [];
                    this.imageVertices = [];
                    var accuStep = 90;
                    var diskIdx = 0;
                    this.imageVertices.push(0.0);
                    this.imageVertices.push(0.0);
                    this.imageVertices.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(-1.0);
                    diskIndices.push(diskIdx++);
                    for(var theta=45;theta<=405;theta+=accuStep){
                        var theta1 = Math.PI *(theta) / 180.0;
                        var x1 = Math.cos(theta1);
                        var y1 = Math.sin(theta1);
                        this.imageVertices.push(x1);
                        this.imageVertices.push(-y1);
                        this.imageVertices.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        diskIndices.push(diskIdx++);
                    }
                    this.imageBuffer = new DisplayBuffer();
                    this.imageBuffer.triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.imageBuffer.triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.imageBuffer.triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.imageBuffer.triangleVertexTextureBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.imageBuffer.triangleVertexIndexBuffer[0]);
                    this.imageBuffer.triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.imageBuffer.triangleVertexIndexBuffer[0].numItems = diskIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(diskIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(diskIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.imageBuffer.triangleVertexNormalBuffer[0]);
                    this.imageBuffer.triangleVertexNormalBuffer[0].itemSize = 3;
                    this.imageBuffer.triangleVertexNormalBuffer[0].numItems = diskNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.imageBuffer.triangleVertexPositionBuffer[0]);
                    this.imageBuffer.triangleVertexPositionBuffer[0].itemSize = 3;
                    this.imageBuffer.triangleVertexPositionBuffer[0].numItems = this.imageVertices.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.imageVertices), this.gl.DYNAMIC_DRAW);

                    var imageTextures = [0.5,0.5, 1.0,1.0, 0.0,1.0, 0.0,0.0, 1.0,0.0, 1.0,1.0];
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.imageBuffer.triangleVertexTextureBuffer[0]);
                    this.imageBuffer.triangleVertexTextureBuffer[0].itemSize = 2;
                    this.imageBuffer.triangleVertexTextureBuffer[0].numItems = imageTextures.length/2;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(imageTextures), this.gl.STATIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j]==="SQUARES"){
                if(typeof(this.squareBuffer) === "undefined"){
                    var diskIndices = [];
                    var diskNormals = [];
                    this.squareVertices = [];
                    var accuStep = 90;
                    var diskIdx = 0;
                    this.squareVertices.push(0.0);
                    this.squareVertices.push(0.0);
                    this.squareVertices.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(-1.0);
                    diskIndices.push(diskIdx++);
                    for(var theta=45;theta<=405;theta+=accuStep){
                        var theta1 = Math.PI *(theta) / 180.0;
                        var y1 = Math.cos(theta1);
                        var x1 = Math.sin(theta1);
                        this.squareVertices.push(x1);
                        this.squareVertices.push(y1);
                        this.squareVertices.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        diskIndices.push(diskIdx++);
                    }
                    this.squareBuffer = new DisplayBuffer();
                    this.squareBuffer.triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.squareBuffer.triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.squareBuffer.triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.squareBuffer.triangleVertexIndexBuffer[0]);
                    this.squareBuffer.triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.squareBuffer.triangleVertexIndexBuffer[0].numItems = diskIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(diskIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(diskIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.squareBuffer.triangleVertexNormalBuffer[0]);
                    this.squareBuffer.triangleVertexNormalBuffer[0].itemSize = 3;
                    this.squareBuffer.triangleVertexNormalBuffer[0].numItems = diskNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.squareBuffer.triangleVertexPositionBuffer[0]);
                    this.squareBuffer.triangleVertexPositionBuffer[0].itemSize = 3;
                    this.squareBuffer.triangleVertexPositionBuffer[0].numItems = this.squareVertices.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.squareVertices), this.gl.DYNAMIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j]==="PENTAGONS"){
                if(typeof(this.pentagonBuffer) === "undefined"){
                    var diskIndices = [];
                    var diskNormals = [];
                    this.pentagonVertices = [];
                    var accuStep = 72;
                    var diskIdx = 0;
                    this.pentagonVertices.push(0.0);
                    this.pentagonVertices.push(0.0);
                    this.pentagonVertices.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(-1.0);
                    diskIndices.push(diskIdx++);
                    for(var theta=0;theta<=360;theta+=accuStep){
                        var theta1 = Math.PI *(theta) / 180.0;
                        var y1 = Math.cos(theta1);
                        var x1 = Math.sin(theta1);
                        this.pentagonVertices.push(x1);
                        this.pentagonVertices.push(y1);
                        this.pentagonVertices.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        diskIndices.push(diskIdx++);
                    }
                    this.pentagonBuffer = new DisplayBuffer();
                    this.pentagonBuffer.triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.pentagonBuffer.triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.pentagonBuffer.triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.pentagonBuffer.triangleVertexIndexBuffer[0]);
                    this.pentagonBuffer.triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.pentagonBuffer.triangleVertexIndexBuffer[0].numItems = diskIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(diskIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(diskIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.pentagonBuffer.triangleVertexNormalBuffer[0]);
                    this.pentagonBuffer.triangleVertexNormalBuffer[0].itemSize = 3;
                    this.pentagonBuffer.triangleVertexNormalBuffer[0].numItems = diskNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.pentagonBuffer.triangleVertexPositionBuffer[0]);
                    this.pentagonBuffer.triangleVertexPositionBuffer[0].itemSize = 3;
                    this.pentagonBuffer.triangleVertexPositionBuffer[0].numItems = this.pentagonVertices.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.pentagonVertices), this.gl.DYNAMIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j]==="HEXAGONS"){
                if(typeof(this.hexagonBuffer) === "undefined"){
                    var diskIndices = [];
                    var diskNormals = [];
                    this.hexagonVertices = [];
                    var accuStep = 60;
                    var diskIdx = 0;
                    this.hexagonVertices.push(0.0);
                    this.hexagonVertices.push(0.0);
                    this.hexagonVertices.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(-1.0);
                    diskIndices.push(diskIdx++);
                    for(var theta=0;theta<=360;theta+=accuStep){
                        var theta1 = Math.PI *(theta) / 180.0;
                        var y1 = Math.cos(theta1);
                        var x1 = Math.sin(theta1);
                        this.hexagonVertices.push(y1);
                        this.hexagonVertices.push(x1);
                        this.hexagonVertices.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        diskIndices.push(diskIdx++);
                    }
                    this.hexagonBuffer = new DisplayBuffer();
                    this.hexagonBuffer.triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.hexagonBuffer.triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.hexagonBuffer.triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.hexagonBuffer.triangleVertexIndexBuffer[0]);
                    this.hexagonBuffer.triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.hexagonBuffer.triangleVertexIndexBuffer[0].numItems = diskIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(diskIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(diskIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.hexagonBuffer.triangleVertexNormalBuffer[0]);
                    this.hexagonBuffer.triangleVertexNormalBuffer[0].itemSize = 3;
                    this.hexagonBuffer.triangleVertexNormalBuffer[0].numItems = diskNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.hexagonBuffer.triangleVertexPositionBuffer[0]);
                    this.hexagonBuffer.triangleVertexPositionBuffer[0].itemSize = 3;
                    this.hexagonBuffer.triangleVertexPositionBuffer[0].numItems = this.hexagonVertices.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.hexagonVertices), this.gl.DYNAMIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j]==="POINTS"){
                if(typeof(this.diskBuffer) === "undefined"){
                    var diskIndices = [];
                    var diskNormals = [];
                    this.diskVertices = [];
                    var accuStep = 10;
                    var diskIdx = 0;
                    this.diskVertices.push(0.0);
                    this.diskVertices.push(0.0);
                    this.diskVertices.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(0.0);
                    diskNormals.push(-1.0);
                    diskIndices.push(diskIdx++);
                    for(var theta=0;theta<=360;theta+=accuStep){
                        var theta1 = Math.PI *(theta) / 180.0;
                        var y1 = Math.cos(theta1);
                        var x1 = Math.sin(theta1);
                        this.diskVertices.push(x1);
                        this.diskVertices.push(y1);
                        this.diskVertices.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        diskIndices.push(diskIdx++);
                    }
                    this.diskBuffer = new DisplayBuffer();
                    this.diskBuffer.triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.diskBuffer.triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.diskBuffer.triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.diskBuffer.triangleVertexIndexBuffer[0]);
                    this.diskBuffer.triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.diskBuffer.triangleVertexIndexBuffer[0].numItems = diskIndices.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(diskIndices), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(diskIndices), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.diskBuffer.triangleVertexNormalBuffer[0]);
                    this.diskBuffer.triangleVertexNormalBuffer[0].itemSize = 3;
                    this.diskBuffer.triangleVertexNormalBuffer[0].numItems = diskNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.diskBuffer.triangleVertexPositionBuffer[0]);
                    this.diskBuffer.triangleVertexPositionBuffer[0].itemSize = 3;
                    this.diskBuffer.triangleVertexPositionBuffer[0].numItems = this.diskVertices.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.diskVertices), this.gl.DYNAMIC_DRAW);
                }
            } else if(this.displayBuffers[idx].bufferTypes[j].substring(0,"CUSTOM_2D_SHAPE_".length)==="CUSTOM_2D_SHAPE_"){
                var customType = this.displayBuffers[idx].bufferTypes[j];
                console.log("A custom 2d shape:"+customType);
                if(typeof(this.shapesBuffers[customType]) === "undefined"){
                    console.log("Make a custom buffer");
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]] = new DisplayBuffer();
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer.push(this.gl.createBuffer());
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer.push(this.gl.createBuffer());
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer.push(this.gl.createBuffer());
                    console.log(this.displayBuffers[idx].supplementary);
                    var vert2d = this.displayBuffers[idx].supplementary["vertices2d"];
                    this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]] = [];
                    var idxs3d = [];
                    var idx3d = 0;
                    var diskNormals = [];
                    for(var i2d=0;i2d<vert2d[0].length;i2d+=2){
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(vert2d[0][i2d]);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(vert2d[0][i2d+1]);
                        this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(0.0);
                        diskNormals.push(-1.0);
                        idxs3d.push(idx3d++);
                    }
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer[0]);
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer[0].itemSize = 1; 
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexIndexBuffer[0].numItems = idxs3d.length;
                    if(this.ext){
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(idxs3d), this.gl.STATIC_DRAW);
                    }else{
                        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(idxs3d), this.gl.STATIC_DRAW);
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer[0]);
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer[0].itemSize = 3;
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexNormalBuffer[0].numItems = diskNormals.length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer[0]);
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer[0].itemSize = 3;
                    this.shapesBuffers[this.displayBuffers[idx].bufferTypes[j]].triangleVertexPositionBuffer[0].numItems = this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]].length/3;
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.shapesVertices[this.displayBuffers[idx].bufferTypes[j]]), this.gl.DYNAMIC_DRAW);
                    /*
                    */
                }
            } else if(this.displayBuffers[idx].bufferTypes[j]==="CIRCLES"){
                var PIBY2 = Math.PI * 2;
                var primitiveSizes = this.displayBuffers[idx].primitiveSizes[j];
                var triangleIndexs = [];
                var triangleNormals = [];
                var triangleVertices = [];
                var triangleColours = [];
                var torusIdx = 0;
                var icol = 0;
                var tor_idx = 0;
                for(var k=0;k<this.displayBuffers[idx].triangleVertices[j].length;k+=3, icol+=4,tor_idx++){
                    var torusOrigin = vec3.create([this.displayBuffers[idx].triangleVertices[j][k],this.displayBuffers[idx].triangleVertices[j][k+1],this.displayBuffers[idx].triangleVertices[j][k+2]]);
                    var torusNormal = vec3.create([this.displayBuffers[idx].triangleNormals[j][k],this.displayBuffers[idx].triangleNormals[j][k+1],this.displayBuffers[idx].triangleNormals[j][k+2]]);
                    var torusColour   = [this.displayBuffers[idx].triangleColours[j][icol],this.displayBuffers[idx].triangleColours[j][icol+1],this.displayBuffers[idx].triangleColours[j][icol+2],this.displayBuffers[idx].triangleColours[j][icol+3]];
                    vec3.normalize(torusNormal);
                    vec3.cross(xaxis,torusNormal,Q);
                    var valid = false;
                    if(vec3.length(Q)>1e-5){
                        vec3.normalize(Q);
                        vec3.cross(torusNormal,Q,R);
                        valid = true;
                    } else {
                        vec3.cross(yaxis,torusNormal,Q);
                        if(vec3.length(Q)>1e-5){
                            vec3.normalize(Q);
                            vec3.cross(torusNormal,Q,R);
                            valid = true;
                        } else {
                            vec3.cross(zaxis,torusNormal,Q);
                            if(vec3.length(Q)>1e-5){
                                vec3.normalize(Q);
                                vec3.cross(torusNormal,Q,R);
                                valid = true;
                            }
                        }
                    }
                    var mat = mat4.create([R[0],R[1],R[2],0.0,   Q[0],Q[1],Q[2],0.0, torusNormal[0],torusNormal[1],torusNormal[2],0.0,  0.0,0.0,0.0,1.0]);
                    var nsectors = 180;
                    var startAngle = 0.0;
                    var sweepAngle = 360.0;
                    var sa;
                    var ea;
                    if(sweepAngle>0){
                        sa = startAngle;
                        ea = startAngle+sweepAngle;
                    }else{
                        sa = startAngle+sweepAngle;
                        ea = startAngle;
                    }
                    var iloop = 0;
                    var radius = this.displayBuffers[idx].supplementary["radii"][j][tor_idx];
                    var majorRadius = radius;
                    for(var jtor=sa;jtor<ea;jtor=jtor+360/nsectors,iloop++){
                        var phi = 1.0*jtor/360.0 * PIBY2;
                        var phi2 = 1.0*(jtor+360.0/nsectors)/360.0 * PIBY2;
                        if(sweepAngle>0&&jtor+360.0/nsectors>startAngle+sweepAngle)  phi2 = 1.0*(startAngle+sweepAngle)/360.0 * PIBY2;
                        if(sweepAngle<0&&jtor+360.0/nsectors>startAngle) phi2 = 1.0*(startAngle)/360.0 * PIBY2;

                        var x = (majorRadius) * Math.cos(phi);
                        var y = (majorRadius) * Math.sin(phi);
                        var z = 0.0;
                        var norm_x = 0.0;
                        var norm_y = 0.0;
                        var norm_z = 1.0;

                        var x2 = (majorRadius) * Math.cos(phi2);
                        var y2 = (majorRadius) * Math.sin(phi2);
                        var z2 = 0.0;

                        var p1 = vec3.create([x,y,z]);
                        var p2 = vec3.create([x2,y2,z2]);

                        mat4.multiplyVec3(mat,p1);
                        mat4.multiplyVec3(mat,p2);

                        x = p1[0]; y = p1[1]; z = p1[2];
                        x2 = p2[0]; y2 = p2[1]; z2 = p2[2];

                        triangleVertices.push(torusOrigin[0]+x); triangleVertices.push(torusOrigin[1]+y); triangleVertices.push(torusOrigin[2]+z);
                        triangleVertices.push(torusOrigin[0]+x2); triangleVertices.push(torusOrigin[1]+y2); triangleVertices.push(torusOrigin[2]+z2);
                        triangleColours.push(torusColour[0]); triangleColours.push(torusColour[1]); triangleColours.push(torusColour[2]); triangleColours.push(torusColour[3]);
                        triangleColours.push(torusColour[0]); triangleColours.push(torusColour[1]); triangleColours.push(torusColour[2]); triangleColours.push(torusColour[3]);
                        triangleIndexs.push(torusIdx++);
                        triangleIndexs.push(torusIdx++);
                    }
                }

                // Try thick lines
                var size = 1.0;
                var thickLines = this.linesToThickLines(triangleVertices,triangleColours,size);
                var Normals_new = thickLines["normals"];
                var Vertices_new = thickLines["vertices"];
                var Colours_new = thickLines["colours"];
                var Indexs_new = thickLines["indices"];
                //console.log("Buffering "+Normals_new.length/3+" normals");
                //console.log("Buffering "+Vertices_new.length/3+" vertices");
                //console.log("Buffering "+Colours_new.length/4+" colours");
                //console.log("Buffering "+Indexs_new.length+" indices");
                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexIndexBuffer[j]);
                if(this.ext){
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(Indexs_new), this.gl.STATIC_DRAW);
                }else{
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(Indexs_new), this.gl.STATIC_DRAW);
                }
                this.displayBuffers[idx].triangleVertexIndexBuffer[j].itemSize = 1; 
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexNormalBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(Normals_new), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexNormalBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexPositionBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(Vertices_new), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleColourBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(Colours_new), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleColourBuffer[j].itemSize = 4;

                this.displayBuffers[idx].triangleVertexIndexBuffer[j].numItems = Indexs_new.length;
                this.displayBuffers[idx].triangleVertexNormalBuffer[j].numItems = Normals_new.length/3;
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].numItems = Vertices_new.length/3;
                this.displayBuffers[idx].triangleColourBuffer[j].numItems = Colours_new.length/4;
                
            } else if(this.displayBuffers[idx].bufferTypes[j]==="TORUSES"){
                var PIBY2 = Math.PI * 2;
                var primitiveSizes = this.displayBuffers[idx].primitiveSizes[j];
                var triangleIndexs = [];
                var triangleNormals = [];
                var triangleVertices = [];
                var triangleColours = [];
                var torusIdx = 0;
                var icol = 0;
                var tor_idx = 0;
                for(var k=0;k<this.displayBuffers[idx].triangleVertices[j].length;k+=3, icol+=4,tor_idx++){
                    var torusOrigin = vec3.create([this.displayBuffers[idx].triangleVertices[j][k],this.displayBuffers[idx].triangleVertices[j][k+1],this.displayBuffers[idx].triangleVertices[j][k+2]]);
                    var torusNormal = vec3.create([this.displayBuffers[idx].triangleNormals[j][k],this.displayBuffers[idx].triangleNormals[j][k+1],this.displayBuffers[idx].triangleNormals[j][k+2]]);
                    var torusColour   = [this.displayBuffers[idx].triangleColours[j][icol],this.displayBuffers[idx].triangleColours[j][icol+1],this.displayBuffers[idx].triangleColours[j][icol+2],this.displayBuffers[idx].triangleColours[j][icol+3]];
                    vec3.normalize(torusNormal);
                    vec3.cross(xaxis,torusNormal,Q);
                    var valid = false;
                    if(vec3.length(Q)>1e-5){
                        vec3.normalize(Q);
                        vec3.cross(torusNormal,Q,R);
                        valid = true;
                    } else {
                        vec3.cross(yaxis,torusNormal,Q);
                        if(vec3.length(Q)>1e-5){
                            vec3.normalize(Q);
                            vec3.cross(torusNormal,Q,R);
                            valid = true;
                        } else {
                            vec3.cross(zaxis,torusNormal,Q);
                            if(vec3.length(Q)>1e-5){
                                vec3.normalize(Q);
                                vec3.cross(torusNormal,Q,R);
                                valid = true;
                            }
                        }
                    }
                    var mat = mat4.create([R[0],R[1],R[2],0.0,   Q[0],Q[1],Q[2],0.0, torusNormal[0],torusNormal[1],torusNormal[2],0.0,  0.0,0.0,0.0,1.0]);
                    var nsectors = 36;
                    var startAngle = 0.0;
                    var sweepAngle = 360.0;
                    var sa;
                    var ea;
                    if(sweepAngle>0){
                        sa = startAngle;
                        ea = startAngle+sweepAngle;
                    }else{
                        sa = startAngle+sweepAngle;
                        ea = startAngle;
                    }
                    var iloop = 0;
                    var radius = this.displayBuffers[idx].supplementary["radii"][j][tor_idx];
                    var majorRadius = radius;
                    var minorRadius = primitiveSizes[tor_idx];
                    for(var jtor=sa;jtor<ea;jtor=jtor+360/nsectors,iloop++){
                        var phi = 1.0*jtor/360.0 * PIBY2;
                        var phi2 = 1.0*(jtor+360.0/nsectors)/360.0 * PIBY2;
                        if(sweepAngle>0&&jtor+360.0/nsectors>startAngle+sweepAngle)  phi2 = 1.0*(startAngle+sweepAngle)/360.0 * PIBY2;
                        if(sweepAngle<0&&jtor+360.0/nsectors>startAngle) phi2 = 1.0*(startAngle)/360.0 * PIBY2;
                        for(var itor=0;itor<=360;itor=itor+360/nsectors){
                            var theta = 1.0*itor/360.0 * PIBY2;
                            var theta2 = (1.0*itor+360.0/nsectors)/360.0 * PIBY2;

                            var x = (majorRadius +  minorRadius * Math.cos(theta)) * Math.cos(phi);
                            var y = (majorRadius +  minorRadius * Math.cos(theta)) * Math.sin(phi);
                            var z = minorRadius * Math.sin(theta);
                            var norm_x = Math.cos(theta) * Math.cos(phi);
                            var norm_y = Math.cos(theta) * Math.sin(phi);
                            var norm_z = Math.sin(theta);

                            var x2 = (majorRadius +  minorRadius * Math.cos(theta)) * Math.cos(phi2);
                            var y2 = (majorRadius +  minorRadius * Math.cos(theta)) * Math.sin(phi2);
                            var z2 = minorRadius * Math.sin(theta);
                            var norm_x2 = Math.cos(theta) * Math.cos(phi2);
                            var norm_y2 = Math.cos(theta) * Math.sin(phi2);
                            var norm_z2 = Math.sin(theta);

                            var x3 = (majorRadius +  minorRadius * Math.cos(theta2)) * Math.cos(phi);
                            var y3 = (majorRadius +  minorRadius * Math.cos(theta2)) * Math.sin(phi);
                            var z3 = minorRadius * Math.sin(theta2);
                            var norm_x3 = Math.cos(theta2) * Math.cos(phi);
                            var norm_y3 = Math.cos(theta2) * Math.sin(phi);
                            var norm_z3 = Math.sin(theta2);

                            var x4 = (majorRadius +  minorRadius * Math.cos(theta2)) * Math.cos(phi2);
                            var y4 = (majorRadius +  minorRadius * Math.cos(theta2)) * Math.sin(phi2);
                            var z4 = minorRadius * Math.sin(theta2);
                            var norm_x4 = Math.cos(theta2) * Math.cos(phi2);
                            var norm_y4 = Math.cos(theta2) * Math.sin(phi2);
                            var norm_z4 = Math.sin(theta2);

                            var p1 = vec3.create([x,y,z]);
                            var p2 = vec3.create([x2,y2,z2]);
                            var p3 = vec3.create([x3,y3,z3]);
                            var p4 = vec3.create([x4,y4,z4]);

                            var n1 = vec3.create([norm_x,norm_y,norm_z]);
                            var n2 = vec3.create([norm_x2,norm_y2,norm_z2]);
                            var n3 = vec3.create([norm_x3,norm_y3,norm_z3]);
                            var n4 = vec3.create([norm_x4,norm_y4,norm_z4]);

                            mat4.multiplyVec3(mat,p1);
                            mat4.multiplyVec3(mat,p2);
                            mat4.multiplyVec3(mat,p3);
                            mat4.multiplyVec3(mat,p4);
                            mat4.multiplyVec3(mat,n1);
                            mat4.multiplyVec3(mat,n2);
                            mat4.multiplyVec3(mat,n3);
                            mat4.multiplyVec3(mat,n4);

                            x = p1[0]; y = p1[1]; z = p1[2];
                            x2 = p2[0]; y2 = p2[1]; z2 = p2[2];
                            x3 = p3[0]; y3 = p3[1]; z3 = p3[2];
                            x4 = p4[0]; y4 = p4[1]; z4 = p4[2];
                            norm_x  = n1[0]; norm_y  = n1[1]; norm_z  = n1[2];
                            norm_x2 = n2[0]; norm_y2 = n2[1]; norm_z2 = n2[2];
                            norm_x3 = n3[0]; norm_y3 = n3[1]; norm_z3 = n3[2];
                            norm_x4 = n4[0]; norm_y4 = n4[1]; norm_z4 = n4[2];

                            triangleVertices.push(torusOrigin[0]+x2); triangleVertices.push(torusOrigin[1]+y2); triangleVertices.push(torusOrigin[2]+z2);
                            triangleVertices.push(torusOrigin[0]+x); triangleVertices.push(torusOrigin[1]+y); triangleVertices.push(torusOrigin[2]+z);
                            triangleVertices.push(torusOrigin[0]+x3); triangleVertices.push(torusOrigin[1]+y3); triangleVertices.push(torusOrigin[2]+z3);
                            triangleNormals.push(norm_x2); triangleNormals.push(norm_y2); triangleNormals.push(norm_z2);
                            triangleNormals.push(norm_x); triangleNormals.push(norm_y); triangleNormals.push(norm_z);
                            triangleNormals.push(norm_x3); triangleNormals.push(norm_y3); triangleNormals.push(norm_z3);
                            triangleColours.push(torusColour[0]); triangleColours.push(torusColour[1]); triangleColours.push(torusColour[2]); triangleColours.push(torusColour[3]);
                            triangleColours.push(torusColour[0]); triangleColours.push(torusColour[1]); triangleColours.push(torusColour[2]); triangleColours.push(torusColour[3]);
                            triangleColours.push(torusColour[0]); triangleColours.push(torusColour[1]); triangleColours.push(torusColour[2]); triangleColours.push(torusColour[3]);
                            triangleIndexs.push(torusIdx++);
                            triangleIndexs.push(torusIdx++);
                            triangleIndexs.push(torusIdx++);

                            triangleVertices.push(torusOrigin[0]+x4); triangleVertices.push(torusOrigin[1]+y4); triangleVertices.push(torusOrigin[2]+z4);
                            triangleVertices.push(torusOrigin[0]+x2); triangleVertices.push(torusOrigin[1]+y2); triangleVertices.push(torusOrigin[2]+z2);
                            triangleVertices.push(torusOrigin[0]+x3); triangleVertices.push(torusOrigin[1]+y3); triangleVertices.push(torusOrigin[2]+z3);
                            triangleNormals.push(norm_x4); triangleNormals.push(norm_y4); triangleNormals.push(norm_z4);
                            triangleNormals.push(norm_x2); triangleNormals.push(norm_y2); triangleNormals.push(norm_z2);
                            triangleNormals.push(norm_x3); triangleNormals.push(norm_y3); triangleNormals.push(norm_z3);
                            triangleColours.push(torusColour[0]); triangleColours.push(torusColour[1]); triangleColours.push(torusColour[2]); triangleColours.push(torusColour[3]);
                            triangleColours.push(torusColour[0]); triangleColours.push(torusColour[1]); triangleColours.push(torusColour[2]); triangleColours.push(torusColour[3]);
                            triangleColours.push(torusColour[0]); triangleColours.push(torusColour[1]); triangleColours.push(torusColour[2]); triangleColours.push(torusColour[3]);
                            triangleIndexs.push(torusIdx++);
                            triangleIndexs.push(torusIdx++);
                            triangleIndexs.push(torusIdx++);
                        }
                    }
                }

                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexIndexBuffer[j]);
                if(this.ext){
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(triangleIndexs), this.gl.STATIC_DRAW);
                }else{
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(triangleIndexs), this.gl.STATIC_DRAW);
                }

                this.displayBuffers[idx].triangleVertexNormalBuffer[j].numItems   = triangleNormals.length/3;
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].numItems = triangleNormals.length/3;
                this.displayBuffers[idx].triangleColourBuffer[j].numItems         = triangleColours.length/4;
                this.displayBuffers[idx].triangleVertexIndexBuffer[j].numItems    = triangleIndexs.length;
                /*
                console.log("Buffering "+triangleVertices.length/3+" vertices");
                console.log("Buffering "+triangleNormals.length/3+" normals");
                console.log("Buffering "+triangleColours.length/4+" colours");
                console.log("Buffering "+triangleIndexs.length+" indices");
                console.log(triangleVertices);
                console.log(triangleNormals);
                console.log(triangleColours);
                console.log(triangleIndexs);
                */

                this.displayBuffers[idx].triangleVertexIndexBuffer[j].itemSize = 1; 
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexNormalBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(triangleNormals), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexNormalBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexPositionBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(triangleVertices), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleColourBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(triangleColours), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleColourBuffer[j].itemSize = 4;
                    
            } else if(this.displayBuffers[idx].bufferTypes[j]==="LINES"){
                console.log("Treating lines specially");
                var size = 1.0;
                var thickLines = this.linesToThickLines(this.displayBuffers[idx].triangleVertices[j],this.displayBuffers[idx].triangleColours[j],size);
                var Normals_new = thickLines["normals"];
                var Vertices_new = thickLines["vertices"];
                var Colours_new = thickLines["colours"];
                var Indexs_new = thickLines["indices"];
                //console.log("Buffering "+Normals_new.length/3+" normals");
                //console.log("Buffering "+Vertices_new.length/3+" vertices");
                //console.log("Buffering "+Colours_new.length/4+" colours");
                //console.log("Buffering "+Indexs_new.length+" indices");
                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexIndexBuffer[j]);
                if(this.ext){
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(Indexs_new), this.gl.STATIC_DRAW);
                }else{
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(Indexs_new), this.gl.STATIC_DRAW);
                }
                this.displayBuffers[idx].triangleVertexIndexBuffer[j].itemSize = 1; 
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexNormalBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(Normals_new), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexNormalBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexPositionBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(Vertices_new), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleColourBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(Colours_new), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleColourBuffer[j].itemSize = 4;

                this.displayBuffers[idx].triangleVertexIndexBuffer[j].numItems = Indexs_new.length;
                this.displayBuffers[idx].triangleVertexNormalBuffer[j].numItems = Normals_new.length/3;
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].numItems = Vertices_new.length/3;
                this.displayBuffers[idx].triangleColourBuffer[j].numItems = Colours_new.length/4;

                //console.log(this.displayBuffers[idx].triangleVertexIndexBuffer[j]);
                //console.log(this.displayBuffers[idx].triangleVertexNormalBuffer[j]);
                //console.log(this.displayBuffers[idx].triangleVertexPositionBuffer[j]);
                //console.log(this.displayBuffers[idx].triangleColourBuffer[j]);
            } else {
                //console.log("This buffer type is "+this.displayBuffers[idx].bufferTypes[j]);
                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexIndexBuffer[j]);
                if(this.ext){
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(this.displayBuffers[idx].triangleIndexs[j]), this.gl.STATIC_DRAW);
                }else{
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.displayBuffers[idx].triangleIndexs[j]), this.gl.STATIC_DRAW);
                }
                this.displayBuffers[idx].triangleVertexIndexBuffer[j].itemSize = 1; 
                if(this.displayBuffers[idx].bufferTypes[j]!="LINES"&&this.displayBuffers[idx].bufferTypes[j]!="LINE_LOOP"&&this.displayBuffers[idx].bufferTypes[j]!="LINE_STRIP"&&this.displayBuffers[idx].bufferTypes[j]!="POINTS"&&this.displayBuffers[idx].bufferTypes[j]!="STARS"&&this.displayBuffers[idx].bufferTypes[j]!="POINTS_SPHERES"&&this.displayBuffers[idx].bufferTypes[j]!="CYLINDERS"&&this.displayBuffers[idx].bufferTypes[j]!="CAPCYLINDERS"&&this.displayBuffers[idx].bufferTypes[j]!="SPLINE"&&this.displayBuffers[idx].bufferTypes[j]!="WORM"&&this.displayBuffers[idx].bufferTypes[j]!="SPHEROIDS"&&this.displayBuffers[idx].bufferTypes[j]!="TORUSES"&&this.displayBuffers[idx].bufferTypes[j]!="CIRCLES"){
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexNormalBuffer[j]);
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[idx].triangleNormals[j]), this.gl.STATIC_DRAW);
                    this.displayBuffers[idx].triangleVertexNormalBuffer[j].itemSize = 3;
                }
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleVertexPositionBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[idx].triangleVertices[j]), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleVertexPositionBuffer[j].itemSize = 3;
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].triangleColourBuffer[j]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[idx].triangleColours[j]), this.gl.STATIC_DRAW);
                this.displayBuffers[idx].triangleColourBuffer[j].itemSize = 4;
            }
        }
    }
    //console.log("Time to build buffers: "+(new Date().getTime()-start));
}


MGWebGL.prototype.drawSymmetryPMV = function(symmetry,buffer,shader,vertexType,bufferIdx,specialDrawBuffer) {
    var triangleVertexIndexBuffer = buffer.triangleVertexIndexBuffer;

    var drawBuffer;
    if(specialDrawBuffer){
        drawBuffer = specialDrawBuffer;
    } else {
        drawBuffer = triangleVertexIndexBuffer[bufferIdx];
    }

    if(buffer.symmetry){
        this.gl.disableVertexAttribArray(shader.vertexColourAttribute);

        var pmvMatrix = mat4.create();
        var screenZ = vec3.create();
        for(isym=0;isym<symmetry["matrices"].length;isym++){
            var symcol = this.symcols[symmetry["symopnums"][isym]%this.symcols.length];
            this.gl.vertexAttrib4f(shader.vertexColourAttribute,symcol[0],symcol[1],symcol[2],symcol[3]);
            var tempMVMatrix = mat4.create();
            var tempMVInvMatrix = mat4.create();
            var symt = mat4.create(symmetry["matrices"][isym]);
            mat4.multiply(this.mvMatrix,symt,tempMVMatrix);
            mat4.multiply(this.pMatrix,tempMVMatrix,pmvMatrix); // Lines
            this.gl.uniformMatrix4fv(shader.pMatrixUniform, false, pmvMatrix); // Lines
            this.gl.uniformMatrix4fv(shader.mvMatrixUniform, false, tempMVMatrix);
            tempMVMatrix[12] = 0.0;
            tempMVMatrix[13] = 0.0;
            tempMVMatrix[14] = 0.0;
            mat4.inverse(tempMVMatrix,tempMVInvMatrix);
            screenZ[0] = 0.0;
            screenZ[1] = 0.0;
            screenZ[2] = 1.0;
            mat4.multiplyVec3(tempMVInvMatrix,screenZ);
            this.gl.uniform3fv(shader.screenZ,screenZ);
            if(this.ext){
                this.gl.drawElements(vertexType, drawBuffer.numItems, this.gl.UNSIGNED_INT, 0);
            }else{
                this.gl.drawElements(vertexType, drawBuffer.numItems, this.gl.UNSIGNED_SHORT, 0);
            }
        }
        this.gl.uniformMatrix4fv(shader.pMatrixUniform, false, this.pmvMatrix); // Lines
        this.gl.uniform3fv(shader.screenZ,this.screenZ); // Lines
        this.gl.enableVertexAttribArray(shader.vertexColourAttribute);
    }
}

MGWebGL.prototype.drawSymmetry = function(symmetry,buffer,shader,vertexType,bufferIdx,specialDrawBuffer) {
    var triangleVertexIndexBuffer = buffer.triangleVertexIndexBuffer;

    var drawBuffer;
    if(specialDrawBuffer){
        drawBuffer = specialDrawBuffer;
    } else {
        drawBuffer = triangleVertexIndexBuffer[bufferIdx];
    }

    if(buffer.symmetry){
        this.gl.disableVertexAttribArray(shader.vertexColourAttribute);

        var pmvMatrix = mat4.create();
        var screenZ = vec3.create();
        for(isym=0;isym<symmetry["matrices"].length;isym++){
            var symcol = this.symcols[symmetry["symopnums"][isym]%this.symcols.length];
            this.gl.vertexAttrib4f(shader.vertexColourAttribute,symcol[0],symcol[1],symcol[2],symcol[3]);
            var tempMVMatrix = mat4.create();
            var tempMVInvMatrix = mat4.create();
            var symt = mat4.create(symmetry["matrices"][isym]);
            mat4.multiply(this.mvMatrix,symt,tempMVMatrix);
            this.gl.uniformMatrix4fv(shader.mvMatrixUniform, false, tempMVMatrix);
            tempMVMatrix[12] = 0.0;
            tempMVMatrix[13] = 0.0;
            tempMVMatrix[14] = 0.0;
            mat4.inverse(tempMVMatrix,tempMVInvMatrix);// All else
            this.gl.uniformMatrix4fv(shader.mvInvMatrixUniform, false, tempMVInvMatrix);// All else
            screenZ[0] = 0.0;
            screenZ[1] = 0.0;
            screenZ[2] = 1.0;
            mat4.multiplyVec3(tempMVInvMatrix,screenZ);
            this.gl.uniform3fv(shader.screenZ,screenZ);
            if(this.ext){
                this.gl.drawElements(vertexType, drawBuffer.numItems, this.gl.UNSIGNED_INT, 0);
            }else{
                this.gl.drawElements(vertexType, drawBuffer.numItems, this.gl.UNSIGNED_SHORT, 0);
            }
        }
        this.gl.uniformMatrix4fv(shader.mvMatrixUniform, false, this.mvMatrix);// All else
        this.gl.uniformMatrix4fv(shader.mvInvMatrixUniform, false, this.mvInvMatrix);// All else
        this.gl.enableVertexAttribArray(shader.vertexColourAttribute);
    }
}

MGWebGL.prototype.render = function(calculatingShadowMap) {
    var width_ratio = this.gl.viewportWidth/this.rttFramebuffer.width;
    var height_ratio = this.gl.viewportHeight/this.rttFramebuffer.height;

    this.mouseDown = false;
    var ratio =  1.0*this.gl.viewportWidth/this.gl.viewportHeight;

    if(calculatingShadowMap){
        this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, this.rttFramebufferDepth);
        this.gl.viewport(0, 0, this.gl.viewportWidth/width_ratio, this.gl.viewportHeight/height_ratio);
    } else {
        this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
        this.gl.viewport(0, 0, this.gl.viewportWidth, this.gl.viewportHeight);
    }
    this.gl.clearColor(this.background_colour[0],this.background_colour[1],this.background_colour[2],this.background_colour[3]);
    this.gl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT);

    //FIXME - This needs to depend on molecule.
    var shadowExtent = 170.0;

    mat4.identity(this.mvMatrix);

    var oldQuat = quat4.create(this.myQuat);
    var newQuat = quat4.create(this.myQuat);

    if(calculatingShadowMap){
        mat4.ortho(-24*ratio/this.zoom,24*ratio/this.zoom,-24/this.zoom,24/this.zoom, 0.000001, shadowExtent, this.pMatrix);
        mat4.translate(this.mvMatrix, [0,0,-shadowExtent*.75]);

        var rotX = quat4.create([0,0,0,-1]);
        var zprime = vec3.create([this.light_positions[0],this.light_positions[1],this.light_positions[2]]);
        vec3.normalize(zprime);
        var zorig = vec3.create([0.0,0.0,1.0]);
        var dp = vec3.dot(zprime,zorig);
        if((1.0-dp)>1e-6){
            var axis = vec3.create();
            vec3.cross(zprime,zorig,axis);
            vec3.normalize(axis);
            var angle = -Math.acos(dp);
            var dval3 = Math.cos(angle/2.0);
            var dval0 = axis[0]*Math.sin(angle/2.0);
            var dval1 = axis[1]*Math.sin(angle/2.0);
            var dval2 = axis[2]*Math.sin(angle/2.0);
            rotX = quat4.create([dval0,dval1,dval2,dval3]);
            quat4.multiply(newQuat,rotX);
        }
        this.gl.enable(this.gl.CULL_FACE);
        this.gl.cullFace(this.gl.FRONT);
    } else {
        this.gl.disable(this.gl.CULL_FACE);
        this.gl.cullFace(this.gl.BACK);
        //mat4.perspective(45, this.gl.viewportWidth / this.gl.viewportHeight, 0.1, 10000.0, this.pMatrix);
        mat4.ortho(-24*ratio,24*ratio,-24,24, 0.000001, 10000.0, this.pMatrix);
        //FIXME - OH hum, this is a problem for fog and clip which assume -500 translation.
        // Hack it for the moment
        if(this.doShadow){
            mat4.translate(this.mvMatrix, [0,0,-shadowExtent*.75]);
        } else {
            mat4.translate(this.mvMatrix, [0,0,-500]);
        }

    }

    this.myQuat = quat4.create(newQuat);
    var theMatrix = quat4.toMat4(this.myQuat);
    mat4.multiply(this.mvMatrix,theMatrix);

    mat4.identity(this.mvInvMatrix);

    var invQuat = quat4.create();
    quat4.inverse(this.myQuat,invQuat);

    //var invMat = mat4.create();
    //mat4.inverse(theMatrix,invMat);
    var invMat =  quat4.toMat4(invQuat);
    this.mvInvMatrix = invMat;

    var right = vec3.create([1.0,0.0,0.0]);
    var up    = vec3.create([0.0,1.0,0.0]);
    mat4.multiplyVec3(invMat,up);
    mat4.multiplyVec3(invMat,right);

    this.screenZ[0] = 0.0;
    this.screenZ[1] = 0.0;
    this.screenZ[2] = 1.0;

    mat4.multiplyVec3(invMat,this.screenZ);

    this.gl.useProgram(this.shaderProgram);
    if(this.backColour==="default"){
        this.gl.uniform1i(this.shaderProgram.defaultColour, true);
    } else {
        this.gl.uniform1i(this.shaderProgram.defaultColour, false);
        this.gl.uniform4fv(this.shaderProgram.backColour, new Float32Array(this.backColour));
    }
    this.gl.uniform1i(this.shaderProgram.shinyBack, this.shinyBack);

    mat4.scale(this.pMatrix, [1./this.zoom,1./this.zoom,1]);
    mat4.translate(this.mvMatrix, this.origin);

    this.pmvMatrix = mat4.create();
    mat4.multiply(this.pMatrix,this.mvMatrix,this.pmvMatrix);

    this.setMatrixUniforms(this.shaderProgram);
    this.setLightUniforms(this.shaderProgram);

    this.gl.useProgram(this.shaderProgramLines);
    this.setMatrixUniforms(this.shaderProgramLines);

    this.gl.useProgram(this.shaderProgramPointSpheres);
    this.setMatrixUniforms(this.shaderProgramPointSpheres);
    this.setLightUniforms(this.shaderProgramPointSpheres);

    this.gl.useProgram(this.shaderProgram);
    this.gl.enableVertexAttribArray(this.shaderProgram.vertexNormalAttribute);

    this.drawTriangles(calculatingShadowMap,invMat);
    this.drawImagesAndText();
    this.drawTransparent(theMatrix);
    this.drawTextLabels(up,right);
    this.drawClickedAtoms(up,right);

    this.myQuat = quat4.create(oldQuat);

    return invMat;

}

MGWebGL.prototype.drawScene = function() {

    var oldMouseDown = this.mouseDown;

    var calculatingShadowMap;

    if(this.doShadow){
        calculatingShadowMap = true;
        this.render(calculatingShadowMap);

        //FIXME - This is all following mgfbo.cc
        this.textureMatrix = mat4.create();
        mat4.identity(this.textureMatrix);
        mat4.translate(this.textureMatrix, [0.5,0.5,0.5]);
        mat4.scale(this.textureMatrix, [0.5,0.5,0.5]);
        mat4.multiply(this.textureMatrix, this.pMatrix);
        mat4.multiply(this.textureMatrix, this.mvMatrix);
    }

    calculatingShadowMap = false;
    var invMat = this.render(calculatingShadowMap);

    if(this.showAxes&&!calculatingShadowMap){
        this.drawAxes(invMat);
    }

    this.mouseDown = oldMouseDown;

    this.div.dispatchEvent(this.viewChangedEvent);

    if(this.save_pixel_data){
        console.log("Saving pixel data");
        var pixels = new Uint8Array(this.canvas.width * this.canvas.height * 4);
        this.gl.readPixels(0, 0, this.canvas.width, this.canvas.height, this.gl.RGBA, this.gl.UNSIGNED_BYTE, pixels);
        this.pixel_data =  pixels;
    }

    if(!this.doShadowDepthDebug||!this.doShadow){
        return
    }

    // This is testing depth buffer.

    this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);

    // This is small example of how we can do more rendering passes
    //this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, this.rttFramebuffer2);

    this.gl.viewport(0, 0, this.gl.viewportWidth, this.gl.viewportHeight);
    this.gl.clearColor(this.background_colour[0],this.background_colour[1],this.background_colour[2],this.background_colour[3]);
    this.gl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT);

    //this.gl.bindTexture(this.gl.TEXTURE_2D, this.rttTexture);
    this.gl.bindTexture(this.gl.TEXTURE_2D, this.rttTextureDepth);

    this.gl.useProgram(this.shaderProgramRenderFrameBuffer);
    this.setMatrixUniforms(this.shaderProgramRenderFrameBuffer);
    this.gl.enableVertexAttribArray(this.shaderProgramRenderFrameBuffer.vertexTextureAttribute);
    mat4.identity(this.mvMatrix);
    mat4.identity(this.mvInvMatrix);
    var right = vec3.create([1.0,0.0,0.0]);
    var up    = vec3.create([0.0,1.0,0.0]);

    mat4.multiplyVec3(invMat,up);
    mat4.multiplyVec3(invMat,right);
    mat4.scale(this.pMatrix, [this.zoom,this.zoom,1]);

    //console.log(up);
    //console.log(right);

    this.gl.uniformMatrix4fv(this.shaderProgramRenderFrameBuffer.pMatrixUniform, false, this.pMatrix);
    this.gl.uniformMatrix4fv(this.shaderProgramRenderFrameBuffer.mvMatrixUniform, false, this.mvMatrix);
    this.gl.uniformMatrix4fv(this.shaderProgramRenderFrameBuffer.mvInvMatrixUniform, false, this.mvInvMatrix);
    this.gl.uniform1f(this.shaderProgramRenderFrameBuffer.fog_start, 0.0);
    this.gl.uniform1f(this.shaderProgramRenderFrameBuffer.fog_end, 1000.0);

    var textNormals   = [   0.0,0.0,1.0,   0.0,0.0,1.0,   0.0,0.0,1.0,   0.0,0.0,1.0,   0.0,0.0,1.0,   0.0,0.0,1.0 ];
    var textColours   = [ 1.0,0.0,0.0,1.0, 1.0,0.0,0.0,1.0, 1.0,0.0,0.0,1.0, 1.0,0.0,0.0,1.0, 1.0,0.0,0.0,1.0, 1.0,0.0,0.0,1.0 ];
    var textTexCoords = [     0.0,0.0,     1.0,0.0,     1.0,1.0,     
                              0.0,0.0,     1.0,1.0,     0.0,1.0 ];
    var textVertices  = [   
        -24.0*right[0]-24.0*up[0],-24.0*right[1]-24.0*up[1],-24.0*right[2]-24.0*up[2],   
        24.0*right[0]-24.0*up[0],24.0*right[1]-24.0*up[1],24.0*right[2]-24.0*up[2],   
        24.0*right[0]+24.0*up[0],24.0*right[1]+24.0*up[1],24.0*right[2]+24.0*up[2],

        -24.0*right[0]-24.0*up[0],-24.0*right[1]-24.0*up[1],-24.0*right[2]-24.0*up[2],
        24.0*up[0]+24.0*right[0],24.0*up[1]+24.0*right[1],24.0*up[2]+24.0*right[2],
        -24.0*right[0]+24.0*up[0],-24.0*right[1]+24.0*up[1],-24.0*right[2]+24.0*up[2] ];
    var textIndexs    = [       0,       1,       2,       3,       4,       5 ];

    var textNormalBuffer = this.gl.createBuffer();
    var textPositionBuffer = this.gl.createBuffer();
    var textColourBuffer = this.gl.createBuffer();
    var textTexCoordBuffer = this.gl.createBuffer();
    var textIndexesBuffer = this.gl.createBuffer();

    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER,textIndexesBuffer);
    if(this.ext){
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(textIndexs), this.gl.STATIC_DRAW);
    }else{
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(textIndexs), this.gl.STATIC_DRAW);
    }

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, textTexCoordBuffer);
    this.gl.vertexAttribPointer(this.shaderProgramRenderFrameBuffer.vertexTextureAttribute, 2, this.gl.FLOAT, false, 0, 0);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(textTexCoords), this.gl.STATIC_DRAW);

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, textNormalBuffer);
    this.gl.vertexAttribPointer(this.shaderProgramRenderFrameBuffer.vertexNormalAttribute, 3, this.gl.FLOAT, false, 0, 0);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(textNormals), this.gl.STATIC_DRAW);

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, textColourBuffer);
    this.gl.vertexAttribPointer(this.shaderProgramRenderFrameBuffer.vertexColourAttribute, 4, this.gl.FLOAT, false, 0, 0);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(textColours), this.gl.STATIC_DRAW);

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, textPositionBuffer);
    this.gl.vertexAttribPointer(this.shaderProgramRenderFrameBuffer.vertexPositionAttribute, 3, this.gl.FLOAT, false, 0, 0);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(textVertices), this.gl.DYNAMIC_DRAW);

    if(this.ext){
        this.gl.drawElements(this.gl.TRIANGLES, textIndexs.length, this.gl.UNSIGNED_INT, 0);
    }else{
        this.gl.drawElements(this.gl.TRIANGLES, textIndexs.length, this.gl.UNSIGNED_SHORT, 0);
    }

    /*
    // This is small example of how we can do more rendering passes
    this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
    this.gl.bindTexture(this.gl.TEXTURE_2D, this.rttTexture2);
    if(this.ext){
        this.gl.drawElements(this.gl.TRIANGLES, textIndexs.length, this.gl.UNSIGNED_INT, 0);
    }else{
        this.gl.drawElements(this.gl.TRIANGLES, textIndexs.length, this.gl.UNSIGNED_SHORT, 0);
    }
    */

    this.gl.disableVertexAttribArray(this.shaderProgramRenderFrameBuffer.vertexTextureAttribute);
}

MGWebGL.prototype.drawTriangles = function(calculatingShadowMap,invMat) {
    var nprims = 0;

    var symmetries = [];
    var symms = [];

    for (var idx = 0; idx < this.displayBuffers.length; idx++){
        if(this.displayBuffers[idx].symmetry){
            var symmetry;
            var thisSym = this.displayBuffers[idx].symmetry;
            if(symms.indexOf(thisSym)==-1){
                symmetry = genSymMats(thisSym["RO"],thisSym["RF"],thisSym["symmats"],[-this.origin[0],-this.origin[1],-this.origin[2]],thisSym["radius"],thisSym["centre"]);
            } else {
                symmetry = symmetries[symms.indexOf(thisSym)];
            }
            symmetries.push(symmetry);
            symms.push(thisSym);
        } else {
            symmetries.push(null);
            symms.push(null);
        }
    }

    for (var idx = 0; idx < this.displayBuffers.length; idx++){

        if(!this.displayBuffers[idx].visible) {
            continue;
        }

        var bufferTypes = this.displayBuffers[idx].bufferTypes;

        var triangleVertexNormalBuffer = this.displayBuffers[idx].triangleVertexNormalBuffer;
        var triangleVertexPositionBuffer = this.displayBuffers[idx].triangleVertexPositionBuffer;
        var triangleVertexIndexBuffer = this.displayBuffers[idx].triangleVertexIndexBuffer;
        var triangleColourBuffer = this.displayBuffers[idx].triangleColourBuffer;

        var triangleIndexs = this.displayBuffers[idx].triangleIndexs;
        var triangleVertices = this.displayBuffers[idx].triangleVertices;
        var triangleColours = this.displayBuffers[idx].triangleColours;
        var triangleNormals = this.displayBuffers[idx].triangleNormals;


        var primitiveSizes = this.displayBuffers[idx].primitiveSizes;

        // FIXME - This is still way too slow, since there can be *lots* of displayBuffers per molecule. 
        //       - recalculating same symmetry for all of them is insane.
        //       - And should only be done when origin changes!
        var symmetry = null;
        if(this.displayBuffers[idx].symmetry){
            symmetry = symmetries[idx];
        }

        //console.log("Drawing object "+idx+" it has "+triangleVertexIndexBuffer.length+" parts");

        for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
            if(bufferTypes[j]=="LINES"||bufferTypes[j]=="LINE_LOOP"||bufferTypes[j]=="LINE_STRIP"||bufferTypes[j]=="DIAMONDS"||bufferTypes[j]=="TEXT"||bufferTypes[j]=="IMAGES"||bufferTypes[j]=="SQUARES"||bufferTypes[j]=="PENTAGONS"||bufferTypes[j]=="HEXAGONS"||bufferTypes[j]=="POINTS"||bufferTypes[j]=="SPHEROIDS"||bufferTypes[j]=="POINTS_SPHERES"||bufferTypes[j]=="STARS"||bufferTypes[j].substring(0,7)=="POLYGON"||bufferTypes[j].substring(0,8)=="POLYSTAR"||bufferTypes[j].substring(0,"CUSTOM_2D_SHAPE_".length)==="CUSTOM_2D_SHAPE_"||bufferTypes[j]=="PERFECT_SPHERES"){
                continue;
            }
            if(this.displayBuffers[idx].transparent) {
                console.log("Not doing normal drawing way ....");
                continue;
            }

            if(this.doShadow){
                if(calculatingShadowMap){
                    // We use a simple shader here.
                    this.gl.useProgram(this.shaderProgramShadow);
                    this.setMatrixUniforms(this.shaderProgramShadow);
                } else {
                    // FIXME We are developing a more complex shader here, one that does shadow texture lookup.
                    // Currently everything is black ....
                    this.gl.useProgram(this.shaderProgramTriangleShadow);
                    this.setMatrixUniforms(this.shaderProgramTriangleShadow);
                    this.setLightUniforms(this.shaderProgramTriangleShadow);
                    //console.log(this.shaderProgramTriangleShadow.textureMatrixUniform);
                    //console.log(this.textureMatrix);
                    var ShadowMapLoc = this.gl.getUniformLocation(this.shaderProgramTriangleShadow,"ShadowMap");
                    this.gl.uniform1i(ShadowMapLoc, 0);
                    this.gl.activeTexture(this.gl.TEXTURE0);
                    this.gl.bindTexture(this.gl.TEXTURE_2D, this.rttTextureDepth);
                    //console.log("Set the active texture!!");
                    this.gl.uniformMatrix4fv(this.shaderProgramTriangleShadow.textureMatrixUniform, false, this.textureMatrix);
                }
            } else {
                this.gl.useProgram(this.shaderProgram);
            }

            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, triangleVertexNormalBuffer[j]);
            this.gl.vertexAttribPointer(this.shaderProgram.vertexNormalAttribute, triangleVertexNormalBuffer[j].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, triangleVertexPositionBuffer[j]);
            this.gl.vertexAttribPointer(this.shaderProgram.vertexPositionAttribute, triangleVertexPositionBuffer[j].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, triangleColourBuffer[j]);
            this.gl.vertexAttribPointer(this.shaderProgram.vertexColourAttribute, triangleColourBuffer[j].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, triangleVertexIndexBuffer[j]);
            if(bufferTypes[j]=="TRIANGLES"||bufferTypes[j]=="CYLINDERS"||bufferTypes[j]=="CAPCYLINDERS"||this.displayBuffers[idx].bufferTypes[j]=="TORUSES"){
                if(this.ext){
                    this.gl.drawElements(this.gl.TRIANGLES, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.drawElements(this.gl.TRIANGLES, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_SHORT, 0);
                }
                if(symmetry) this.drawSymmetry(symmetry,this.displayBuffers[idx],this.shaderProgram,this.gl.TRIANGLES,j);
            }else if(bufferTypes[j]=="TRIANGLE_STRIP"||bufferTypes[j]=="SPLINE"||bufferTypes[j]=="WORM"){
                if(this.ext){
                    this.gl.drawElements(this.gl.TRIANGLE_STRIP, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.drawElements(this.gl.TRIANGLE_STRIP, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_SHORT, 0);
                }
                if(symmetry) this.drawSymmetry(symmetry,this.displayBuffers[idx],this.shaderProgram,this.gl.TRIANGLE_STRIP,j);
            }
            nprims += triangleVertexIndexBuffer[j].numItems;
        }

        /*
           Cylinders here
         */

        var sphereProgram = this.shaderProgramPointSpheres;

        if(this.frag_depth_ext){
            if(this.doShadow && !calculatingShadowMap){
                var sphereProgram = this.shaderProgramPointSpheresShadow;
                this.gl.useProgram(sphereProgram);
                this.setMatrixUniforms(sphereProgram);
                this.setLightUniforms(sphereProgram);
                this.gl.disableVertexAttribArray(sphereProgram.vertexColourAttribute);
                this.gl.enableVertexAttribArray(sphereProgram.vertexTextureAttribute);
                var ShadowMapLoc = this.gl.getUniformLocation(sphereProgram,"ShadowMap");
                this.gl.uniform1i(ShadowMapLoc, 0);
                this.gl.activeTexture(this.gl.TEXTURE0);
                this.gl.bindTexture(this.gl.TEXTURE_2D, this.rttTextureDepth);
                this.gl.uniformMatrix4fv(sphereProgram.textureMatrixUniform, false, this.textureMatrix);
            }
        }

        this.gl.useProgram(sphereProgram);

        var scaleMatrices = this.displayBuffers[idx].supplementary["scale_matrices"];
        this.gl.disableVertexAttribArray(sphereProgram.vertexColourAttribute);

        for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
            var theseScaleMatrices = [];
            if(bufferTypes[j]!="SPHEROIDS"&&bufferTypes[j]!="POINTS_SPHERES"&&bufferTypes[j]!="STARS"){
                continue;
            }
            var buffer;
            var radMult;
            if(bufferTypes[j]=="POINTS_SPHERES"||bufferTypes[j]=="SPHEROIDS"){
                buffer = this.sphereBuffer;
                radMult = 1.0;
                if(bufferTypes[j]=="SPHEROIDS"){
                    theseScaleMatrices = scaleMatrices[j];
                }
            } else {
                buffer = this.starBuffer;
                radMult = 0.32;
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexNormalBuffer[0]);
            this.gl.vertexAttribPointer(sphereProgram.vertexNormalAttribute, buffer.triangleVertexNormalBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexPositionBuffer[0]);
            this.gl.vertexAttribPointer(sphereProgram.vertexPositionAttribute, buffer.triangleVertexPositionBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, buffer.triangleVertexIndexBuffer[0]);
            var isphere;

            // FIXME - The scaling will be a property of each object. e.g. B/U factors.
            //       - Perhaps we should have different shaders for scaled objects?
            var scaleMatrix = mat3.create([1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0]);
            this.gl.uniformMatrix3fv(sphereProgram.scaleMatrix, false, scaleMatrix);

            var theOffSet = new glMatrixArrayType(3);
            if(theseScaleMatrices.length===triangleVertices[j].length/3){
                for(isphere=0;isphere<triangleVertices[j].length/3;isphere++){
                    scaleMatrix = mat3.create(theseScaleMatrices[isphere]);
                    this.gl.uniformMatrix3fv(sphereProgram.scaleMatrix, false, scaleMatrix);
                    theOffSet[0] = triangleVertices[j][isphere*3];
                    theOffSet[1] = triangleVertices[j][isphere*3+1];
                    theOffSet[2] = triangleVertices[j][isphere*3+2];
                    var thisCol = [ triangleColours[j][isphere*4], triangleColours[j][isphere*4+1], triangleColours[j][isphere*4+2], triangleColours[j][isphere*4+3] ];
                    this.gl.vertexAttrib4f(sphereProgram.vertexColourAttribute,triangleColours[j][isphere*4], triangleColours[j][isphere*4+1], triangleColours[j][isphere*4+2], triangleColours[j][isphere*4+3]);
                    this.gl.uniform3fv(sphereProgram.offset, theOffSet);
                    this.gl.uniform1f(sphereProgram.size, primitiveSizes[j][isphere]*radMult);
                    if(this.ext){
                        this.gl.drawElements(this.gl.TRIANGLES, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_INT, 0);
                    }else{
                        this.gl.drawElements(this.gl.TRIANGLES, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_SHORT, 0);
                    }
                    nprims += triangleVertexIndexBuffer[j].numItems;
                }
            } else {
                for(isphere=0;isphere<triangleVertices[j].length/3;isphere++){
                    theOffSet[0] = triangleVertices[j][isphere*3];
                    theOffSet[1] = triangleVertices[j][isphere*3+1];
                    theOffSet[2] = triangleVertices[j][isphere*3+2];
                    var thisCol = [ triangleColours[j][isphere*4], triangleColours[j][isphere*4+1], triangleColours[j][isphere*4+2], triangleColours[j][isphere*4+3] ];
                    this.gl.vertexAttrib4f(sphereProgram.vertexColourAttribute,triangleColours[j][isphere*4], triangleColours[j][isphere*4+1], triangleColours[j][isphere*4+2], triangleColours[j][isphere*4+3]);
                    this.gl.uniform3fv(sphereProgram.offset, theOffSet);
                    this.gl.uniform1f(sphereProgram.size, primitiveSizes[j][isphere]*radMult);
                    if(this.ext){
                        this.gl.drawElements(this.gl.TRIANGLES, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_INT, 0);
                    }else{
                        this.gl.drawElements(this.gl.TRIANGLES, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_SHORT, 0);
                    }
                    if(symmetry){
                        this.drawSymmetry(symmetry,this.displayBuffers[idx],sphereProgram,this.gl.TRIANGLES,0,buffer.triangleVertexIndexBuffer[0]);
                        this.gl.disableVertexAttribArray(sphereProgram.vertexColourAttribute);
                    }
                    nprims += triangleVertexIndexBuffer[j].numItems;
                }
            }
        }
        this.gl.useProgram(this.shaderProgramTwoDShapes);
        this.setMatrixUniforms(this.shaderProgramTwoDShapes);
        this.gl.disableVertexAttribArray(this.shaderProgramTwoDShapes.vertexColourAttribute);
        this.gl.vertexAttrib4f(this.shaderProgramTwoDShapes.vertexColourAttribute,1.0,1.0,0.0,1.0);
        var diskVertices = [];
        if(typeof(this.diskVertices) !== "undefined"){
            for(var iv=0;iv<this.diskVertices.length;iv+=3){
                var vold = vec3.create([this.diskVertices[iv],this.diskVertices[iv+1],this.diskVertices[iv+2]]);
                var vnew = vec3.create();
                mat4.multiplyVec3(invMat,vold,vnew);
                diskVertices[iv]   = vnew[0];
                diskVertices[iv+1] = vnew[1];
                diskVertices[iv+2] = vnew[2];
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.diskBuffer.triangleVertexPositionBuffer[0]);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diskVertices), this.gl.DYNAMIC_DRAW);
        }
        var hexagonVertices = [];
        if(typeof(this.hexagonVertices) !== "undefined"){
            for(var iv=0;iv<this.hexagonVertices.length;iv+=3){
                var vold = vec3.create([this.hexagonVertices[iv],this.hexagonVertices[iv+1],this.hexagonVertices[iv+2]]);
                var vnew = vec3.create();
                mat4.multiplyVec3(invMat,vold,vnew);
                hexagonVertices[iv]   = vnew[0];
                hexagonVertices[iv+1] = vnew[1];
                hexagonVertices[iv+2] = vnew[2];
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.hexagonBuffer.triangleVertexPositionBuffer[0]);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(hexagonVertices), this.gl.DYNAMIC_DRAW);
        }
        var pentagonVertices = [];
        if(typeof(this.pentagonVertices) !== "undefined"){
            for(var iv=0;iv<this.pentagonVertices.length;iv+=3){
                var vold = vec3.create([this.pentagonVertices[iv],this.pentagonVertices[iv+1],this.pentagonVertices[iv+2]]);
                var vnew = vec3.create();
                mat4.multiplyVec3(invMat,vold,vnew);
                pentagonVertices[iv]   = vnew[0];
                pentagonVertices[iv+1] = vnew[1];
                pentagonVertices[iv+2] = vnew[2];
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.pentagonBuffer.triangleVertexPositionBuffer[0]);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(pentagonVertices), this.gl.DYNAMIC_DRAW);
        }
        var imageVertices = [];
        if(typeof(this.imageVertices) !== "undefined"){
            for(var iv=0;iv<this.imageVertices.length;iv+=3){
                var vold = vec3.create([this.imageVertices[iv],this.imageVertices[iv+1],this.imageVertices[iv+2]]);
                var vnew = vec3.create();
                mat4.multiplyVec3(invMat,vold,vnew);
                imageVertices[iv]   = vnew[0];
                imageVertices[iv+1] = vnew[1];
                imageVertices[iv+2] = vnew[2];
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.imageBuffer.triangleVertexPositionBuffer[0]);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(imageVertices), this.gl.DYNAMIC_DRAW);
        }
        var squareVertices = [];
        if(typeof(this.squareVertices) !== "undefined"){
            for(var iv=0;iv<this.squareVertices.length;iv+=3){
                var vold = vec3.create([this.squareVertices[iv],this.squareVertices[iv+1],this.squareVertices[iv+2]]);
                var vnew = vec3.create();
                mat4.multiplyVec3(invMat,vold,vnew);
                squareVertices[iv]   = vnew[0];
                squareVertices[iv+1] = vnew[1];
                squareVertices[iv+2] = vnew[2];
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.squareBuffer.triangleVertexPositionBuffer[0]);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(squareVertices), this.gl.DYNAMIC_DRAW);
        }
        var diamondVertices = [];
        if(typeof(this.diamondVertices) !== "undefined"){
            for(var iv=0;iv<this.diamondVertices.length;iv+=3){
                var vold = vec3.create([this.diamondVertices[iv],this.diamondVertices[iv+1],this.diamondVertices[iv+2]]);
                var vnew = vec3.create();
                mat4.multiplyVec3(invMat,vold,vnew);
                diamondVertices[iv]   = vnew[0];
                diamondVertices[iv+1] = vnew[1];
                diamondVertices[iv+2] = vnew[2];
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.diamondBuffer.triangleVertexPositionBuffer[0]);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(diamondVertices), this.gl.DYNAMIC_DRAW);
        }
        for(shapeVertices in this.shapesVertices){
            var theseVertices = [];
            for(var iv=0;iv<this.shapesVertices[shapeVertices].length;iv+=3){
                var vold = vec3.create([this.shapesVertices[shapeVertices][iv],this.shapesVertices[shapeVertices][iv+1],this.shapesVertices[shapeVertices][iv+2]]);
                var vnew = vec3.create();
                mat4.multiplyVec3(invMat,vold,vnew);
                theseVertices[iv]   = vnew[0];
                theseVertices[iv+1] = vnew[1];
                theseVertices[iv+2] = vnew[2];
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.shapesBuffers[shapeVertices].triangleVertexPositionBuffer[0]);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(theseVertices), this.gl.DYNAMIC_DRAW);
        }
        for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
            if(bufferTypes[j].substring(0,"CUSTOM_2D_SHAPE_".length)==="CUSTOM_2D_SHAPE_"){
                var buffer = this.shapesBuffers[bufferTypes[j]];
                var scaleImage = true;
                if(typeof(this.gl,this.displayBuffers[idx].supplementary["vert_tri_2d"])!=="undefined"){
                    var tempMVMatrix = mat4.create([this.mvMatrix[0],this.mvMatrix[1],this.mvMatrix[2],this.mvMatrix[3],this.mvMatrix[4],this.mvMatrix[5],this.mvMatrix[6],this.mvMatrix[7],this.mvMatrix[8],this.mvMatrix[9],this.mvMatrix[10],this.mvMatrix[11], (-24.0 + this.displayBuffers[idx].supplementary["vert_tri_2d"][0][0] * 48.0)*this.zoom,(-24.0 + this.displayBuffers[idx].supplementary["vert_tri_2d"][0][1] * 48.0)*this.zoom,-500.0,1.0 ]);
                    this.gl.uniformMatrix4fv(this.shaderProgramTwoDShapes.mvMatrixUniform, false, tempMVMatrix);
                    scaleImage = false;
                }

                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexNormalBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramTwoDShapes.vertexNormalAttribute, buffer.triangleVertexNormalBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexPositionBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramTwoDShapes.vertexPositionAttribute, buffer.triangleVertexPositionBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, buffer.triangleVertexIndexBuffer[0]);
                var theOffSet = new glMatrixArrayType(3);
                for(ishape=0;ishape<triangleVertices[j].length/3;ishape++){
                    theOffSet[0] = triangleVertices[j][ishape*3];
                    theOffSet[1] = triangleVertices[j][ishape*3+1];
                    theOffSet[2] = triangleVertices[j][ishape*3+2];
                    this.gl.uniform3fv(this.shaderProgramTwoDShapes.offset, theOffSet);
                    if(scaleImage){
                        this.gl.uniform1f(this.shaderProgramTwoDShapes.size, primitiveSizes[j][ishape]);
                    }else{
                        this.gl.uniform1f(this.shaderProgramTwoDShapes.size, primitiveSizes[j][ishape]*this.zoom);
                    }

                    var thisCol = [ triangleColours[j][ishape*4], triangleColours[j][ishape*4+1], triangleColours[j][ishape*4+2], triangleColours[j][ishape*4+3] ];
                    this.gl.vertexAttrib4f(this.shaderProgramTwoDShapes.vertexColourAttribute,triangleColours[j][ishape*4], triangleColours[j][ishape*4+1], triangleColours[j][ishape*4+2], triangleColours[j][ishape*4+3]);

                    if(this.ext){
                        this.gl.drawElements(this.gl.TRIANGLES, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_INT, 0);
                    }else{
                        this.gl.drawElements(this.gl.TRIANGLES, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_SHORT, 0);
                    }
                }
                if(typeof(this.gl,this.displayBuffers[idx].supplementary["vert_tri_2d"])!=="undefined"){
                    this.setMatrixUniforms(this.shaderProgramTwoDShapes);
                }
            }
        }
        for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
            if(bufferTypes[j]=="POINTS"||bufferTypes[j]=="HEXAGONS"||bufferTypes[j]=="PENTAGONS"||bufferTypes[j]=="SQUARES"||bufferTypes[j]=="DIAMONDS"||bufferTypes[j].substring(0,7)=="POLYGON"||bufferTypes[j].substring(0,8)=="POLYSTAR"){
                var buffer;
                if(bufferTypes[j]=="HEXAGONS"){
                    buffer = this.hexagonBuffer;
                } else if(bufferTypes[j]=="DIAMONDS"){
                    buffer = this.diamondBuffer;
                } else if(bufferTypes[j]=="SQUARES"){
                    buffer = this.squareBuffer;
                } else if(bufferTypes[j]=="PENTAGONS"){
                    buffer = this.pentagonBuffer;
                } else if(bufferTypes[j].substring(0,7)=="POLYGON"){
                    buffer = this.shapesBuffers[bufferTypes[j]];
                } else if(bufferTypes[j].substring(0,8)=="POLYSTAR"){
                    buffer = this.shapesBuffers[bufferTypes[j]];
                } else {
                    buffer = this.diskBuffer;
                }
                var scaleImage = true;
                if(typeof(this.gl,this.displayBuffers[idx].supplementary["vert_tri_2d"])!=="undefined"){
                    var tempMVMatrix = mat4.create([this.mvMatrix[0],this.mvMatrix[1],this.mvMatrix[2],this.mvMatrix[3],this.mvMatrix[4],this.mvMatrix[5],this.mvMatrix[6],this.mvMatrix[7],this.mvMatrix[8],this.mvMatrix[9],this.mvMatrix[10],this.mvMatrix[11], (-24.0 + this.displayBuffers[idx].supplementary["vert_tri_2d"][0][0] * 48.0)*this.zoom,(-24.0 + this.displayBuffers[idx].supplementary["vert_tri_2d"][0][1] * 48.0)*this.zoom,-500.0,1.0 ]);
                    this.gl.uniformMatrix4fv(this.shaderProgramTwoDShapes.mvMatrixUniform, false, tempMVMatrix);
                    scaleImage = false;
                }

                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexNormalBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramTwoDShapes.vertexNormalAttribute, buffer.triangleVertexNormalBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexPositionBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramTwoDShapes.vertexPositionAttribute, buffer.triangleVertexPositionBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, buffer.triangleVertexIndexBuffer[0]);
                // FIXME - And loop here
                var theOffSet = new glMatrixArrayType(3);
                for(ishape=0;ishape<triangleVertices[j].length/3;ishape++){
                    theOffSet[0] = triangleVertices[j][ishape*3];
                    theOffSet[1] = triangleVertices[j][ishape*3+1];
                    theOffSet[2] = triangleVertices[j][ishape*3+2];
                    this.gl.uniform3fv(this.shaderProgramTwoDShapes.offset, theOffSet);
                    if(scaleImage){
                        this.gl.uniform1f(this.shaderProgramTwoDShapes.size, primitiveSizes[j][ishape]);
                    }else{
                        this.gl.uniform1f(this.shaderProgramTwoDShapes.size, primitiveSizes[j][ishape]*this.zoom);
                    }

                    var thisCol = [ triangleColours[j][ishape*4], triangleColours[j][ishape*4+1], triangleColours[j][ishape*4+2], triangleColours[j][ishape*4+3] ];
                    this.gl.vertexAttrib4f(this.shaderProgramTwoDShapes.vertexColourAttribute,triangleColours[j][ishape*4], triangleColours[j][ishape*4+1], triangleColours[j][ishape*4+2], triangleColours[j][ishape*4+3]);

                    if(this.ext){
                        this.gl.drawElements(this.gl.TRIANGLE_FAN, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_INT, 0);
                    }else{
                        this.gl.drawElements(this.gl.TRIANGLE_FAN, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_SHORT, 0);
                    }
                }
                if(typeof(this.gl,this.displayBuffers[idx].supplementary["vert_tri_2d"])!=="undefined"){
                    this.setMatrixUniforms(this.shaderProgramTwoDShapes);
                }
            }
        }
        this.gl.enableVertexAttribArray(this.shaderProgramTwoDShapes.vertexColourAttribute);

        this.gl.useProgram(this.shaderProgramLines);
        for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
            if(bufferTypes[j]!="LINE_LOOP"&&bufferTypes[j]!="LINE_STRIP"){
                continue;
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, triangleVertexPositionBuffer[j]);
            this.gl.vertexAttribPointer(this.shaderProgramLines.vertexPositionAttribute, triangleVertexPositionBuffer[j].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, triangleColourBuffer[j]);
            this.gl.vertexAttribPointer(this.shaderProgramLines.vertexColourAttribute, triangleColourBuffer[j].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, triangleVertexIndexBuffer[j]);
            //this.gl.disableVertexAttribArray(2)
            if(bufferTypes[j]=="LINES"){
                //console.log("Try to draw "+triangleVertexIndexBuffer[j].numItems+" lines, "+triangleColourBuffer[j].numItems+" colours, "+triangleVertexPositionBuffer[j].numItems+" vertices.");
                if(this.ext){
                    this.gl.drawElements(this.gl.LINES, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.drawElements(this.gl.LINES, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_SHORT, 0);
                }
            }
            if(bufferTypes[j]=="LINE_STRIP"){
                //console.log("Try to draw "+triangleVertexIndexBuffer[j].numItems+" lines, "+triangleColourBuffer[j].numItems+" colours, "+triangleVertexPositionBuffer[j].numItems+" vertices.");
                if(this.ext){
                    this.gl.drawElements(this.gl.LINE_STRIP, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.drawElements(this.gl.LINE_STRIP, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_SHORT, 0);
                }
            }
            nprims += triangleVertexIndexBuffer[j].numItems;
        }

        this.gl.useProgram(this.shaderProgramThickLines);
        this.gl.uniform3fv(this.shaderProgramThickLines.screenZ,this.screenZ);
        this.gl.enableVertexAttribArray(this.shaderProgramThickLines.vertexNormalAttribute);
        this.setMatrixUniforms(this.shaderProgramThickLines);
        this.gl.uniformMatrix4fv(this.shaderProgramThickLines.pMatrixUniform, false, this.pmvMatrix);
        
        for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
            if(bufferTypes[j]!="LINES"&&bufferTypes[j]!="CIRCLES"){
                continue;
            }
            // FIXME ? We assume all are the same size. Anything else is a little tricky for now.
            if(typeof(this.displayBuffers[idx].primitiveSizes)!=="undefined" && typeof(this.displayBuffers[idx].primitiveSizes[j])!=="undefined"&&typeof(this.displayBuffers[idx].primitiveSizes[j][0])!=="undefined"){
                this.gl.uniform1f(this.shaderProgramThickLines.pixelZoom, this.displayBuffers[idx].primitiveSizes[j][0]*0.04*this.zoom);
            }else{
                this.gl.uniform1f(this.shaderProgramThickLines.pixelZoom, 1.0*0.04*this.zoom);
            }
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, triangleVertexNormalBuffer[j]);
            this.gl.vertexAttribPointer(this.shaderProgramThickLines.vertexNormalAttribute, triangleVertexNormalBuffer[j].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, triangleVertexPositionBuffer[j]);
            this.gl.vertexAttribPointer(this.shaderProgramThickLines.vertexPositionAttribute, triangleVertexPositionBuffer[j].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, triangleColourBuffer[j]);
            this.gl.vertexAttribPointer(this.shaderProgramThickLines.vertexColourAttribute, triangleColourBuffer[j].itemSize, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, triangleVertexIndexBuffer[j]);
            if(this.ext){
                this.gl.drawElements(this.gl.TRIANGLES, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_INT, 0);
            }else{
                this.gl.drawElements(this.gl.TRIANGLES, triangleVertexIndexBuffer[j].numItems, this.gl.UNSIGNED_SHORT, 0);
            }

            if(symmetry) this.drawSymmetryPMV(symmetry,this.displayBuffers[idx],this.shaderProgramThickLines,this.gl.TRIANGLES,j);

            nprims += triangleVertexIndexBuffer[j].numItems;
        }

        //shaderProgramPerfectSpheres
        if(this.frag_depth_ext){
            var program = this.shaderProgramPerfectSpheres;
            if(this.doShadow && !calculatingShadowMap ){
                program = this.shaderProgramPerfectSpheresShadow;
                this.gl.useProgram(program);
                this.setMatrixUniforms(program);
                this.setLightUniforms(program);
                this.gl.disableVertexAttribArray(program.vertexColourAttribute);
                this.gl.enableVertexAttribArray(program.vertexTextureAttribute);
                var ShadowMapLoc = this.gl.getUniformLocation(program,"ShadowMap");
                this.gl.uniform1i(ShadowMapLoc, 0);
                this.gl.activeTexture(this.gl.TEXTURE0);
                this.gl.bindTexture(this.gl.TEXTURE_2D, this.rttTextureDepth);
                //console.log("Shadowing perfect spheres "+ShadowMapLoc);
                this.gl.uniformMatrix4fv(program.textureMatrixUniform, false, this.textureMatrix);
            } else {
                this.gl.useProgram(program);
                this.setMatrixUniforms(program);
                this.setLightUniforms(program);
                this.gl.disableVertexAttribArray(program.vertexColourAttribute);
                this.gl.enableVertexAttribArray(program.vertexTextureAttribute);
            }

            for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
                if(bufferTypes[j]=="PERFECT_SPHERES"){
                    buffer = this.imageBuffer;

                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexTextureBuffer[0]);
                    this.gl.vertexAttribPointer(program.vertexTextureAttribute, buffer.triangleVertexTextureBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);

                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexNormalBuffer[0]);
                    this.gl.vertexAttribPointer(program.vertexNormalAttribute, buffer.triangleVertexNormalBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexPositionBuffer[0]);
                    this.gl.vertexAttribPointer(program.vertexPositionAttribute, buffer.triangleVertexPositionBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, buffer.triangleVertexIndexBuffer[0]);
                    var isphere;

                    var scaleMatrix = mat3.create([1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0]);
                    //this.gl.uniformMatrix3fv(this.shaderProgramPerfectSpheres.scaleMatrix, false, scaleMatrix);

                    var theOffSet = new glMatrixArrayType(3);
                    for(isphere=0;isphere<triangleVertices[j].length/3;isphere++){
                        theOffSet[0] = triangleVertices[j][isphere*3];
                        theOffSet[1] = triangleVertices[j][isphere*3+1];
                        theOffSet[2] = triangleVertices[j][isphere*3+2];
                        var thisCol = [ triangleColours[j][isphere*4], triangleColours[j][isphere*4+1], triangleColours[j][isphere*4+2], triangleColours[j][isphere*4+3] ];
                        this.gl.vertexAttrib4f(program.vertexColourAttribute,triangleColours[j][isphere*4], triangleColours[j][isphere*4+1], triangleColours[j][isphere*4+2], triangleColours[j][isphere*4+3]);
                        this.gl.uniform3fv(program.offset, theOffSet);
                        this.gl.uniform1f(program.size, primitiveSizes[j][isphere]*Math.sqrt(2));
                        if(this.ext){
                            this.gl.drawElements(this.gl.TRIANGLE_FAN, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_INT, 0);
                        }else{
                            this.gl.drawElements(this.gl.TRIANGLE_FAN, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_SHORT, 0);
                        }
                        nprims += triangleVertexIndexBuffer[j].numItems;
                        /*
                        if(symmetry){
                            this.drawSymmetry(symmetry,this.displayBuffers[idx],this.shaderProgramPerfectSpheres,this.gl.TRIANGLE_FAN,0,buffer.triangleVertexIndexBuffer[0]);
                            this.gl.disableVertexAttribArray(this.shaderProgramPerfectSpheres.vertexColourAttribute);
                        }
                        */
                    }

                }
            }

            this.gl.enableVertexAttribArray(program.vertexColourAttribute);
            this.gl.disableVertexAttribArray(program.vertexTextureAttribute);
        }
    }
}

MGWebGL.prototype.drawTransparent = function(theMatrix) {

    for (var idx = 0; idx < this.displayBuffers.length; idx++){

        if(!this.displayBuffers[idx].visible) {
            continue;
        }

        var triangleVertexNormalBuffer = this.displayBuffers[idx].triangleVertexNormalBuffer;
        var triangleVertexPositionBuffer = this.displayBuffers[idx].triangleVertexPositionBuffer;
        var triangleVertexIndexBuffer = this.displayBuffers[idx].triangleVertexIndexBuffer;
        var triangleColourBuffer = this.displayBuffers[idx].triangleColourBuffer;

        var triangleIndexs = this.displayBuffers[idx].triangleIndexs;
        var triangleVertices = this.displayBuffers[idx].triangleVertices;
        var triangleColours = this.displayBuffers[idx].triangleColours;
        var triangleNormals = this.displayBuffers[idx].triangleNormals;

        var bufferTypes = this.displayBuffers[idx].bufferTypes;

        if(this.displayBuffers[idx].transparent){
            //console.log(idx+" is transparent ;-) "+bufferTypes[0]);

            if(bufferTypes[0]=="TRIANGLES"){
                if(typeof this.displayBuffers[idx].allVertices === "undefined") {
                    this.displayBuffers[idx].allVertices = [];
                    this.displayBuffers[idx].allNormals = [];
                    this.displayBuffers[idx].allColours = [];
                    this.displayBuffers[idx].allIndexs = [];
                    this.displayBuffers[idx].allTriangleVertexNormalBuffer = this.gl.createBuffer();
                    this.displayBuffers[idx].allTriangleVertexPositionBuffer = this.gl.createBuffer();
                    this.displayBuffers[idx].allTriangleVertexColourBuffer = this.gl.createBuffer();
                    this.displayBuffers[idx].allIndexsBuffer = this.gl.createBuffer();
                    var bufferOffset = 0;
                    //console.log("Concatenating "+triangleVertexIndexBuffer.length+" buffers...");
                    for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
                        //console.log("Concatenating "+triangleVertices[j].length/3+" triangles");
                        this.displayBuffers[idx].allVertices = this.displayBuffers[idx].allVertices.concat(triangleVertices[j]);
                        //console.log("Concatenating "+triangleNormals[j].length/3+" normals");
                        this.displayBuffers[idx].allNormals = this.displayBuffers[idx].allNormals.concat(triangleNormals[j]);
                        //console.log("total this.displayBuffers[idx].allNormals.length: "+this.displayBuffers[idx].allNormals.length);
                        //console.log("Concatenating "+triangleColours[j].length/4+" colours");
                        this.displayBuffers[idx].allColours = this.displayBuffers[idx].allColours.concat(triangleColours[j]);
                        //console.log("Pushing "+triangleIndexs[j].length+" indices");
                        for(var i = 0; i < triangleIndexs[j].length; i++){
                            this.displayBuffers[idx].allIndexs.push(triangleIndexs[j][i]+bufferOffset);
                        }
                        bufferOffset += triangleVertices[j].length/3;
                    }
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].allTriangleVertexNormalBuffer);
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[idx].allNormals), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].allTriangleVertexPositionBuffer);
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[idx].allVertices), this.gl.STATIC_DRAW);
                    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].allTriangleVertexColourBuffer);
                    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[idx].allColours), this.gl.STATIC_DRAW);
                }
                var sortThings = [];
                //console.log("Big thing is of size "+this.displayBuffers[idx].allIndexs.length);
                for (var j = 0; j < this.displayBuffers[idx].allIndexs.length; j+=3){
                    var x1 = this.displayBuffers[idx].allVertices[3*this.displayBuffers[idx].allIndexs[j]];
                    var y1 = this.displayBuffers[idx].allVertices[3*this.displayBuffers[idx].allIndexs[j]+1];
                    var z1 = this.displayBuffers[idx].allVertices[3*this.displayBuffers[idx].allIndexs[j]+2];
                    var x2 = this.displayBuffers[idx].allVertices[3*this.displayBuffers[idx].allIndexs[j+1]];
                    var y2 = this.displayBuffers[idx].allVertices[3*this.displayBuffers[idx].allIndexs[j+1]+1];
                    var z2 = this.displayBuffers[idx].allVertices[3*this.displayBuffers[idx].allIndexs[j+1]+2];
                    var x3 = this.displayBuffers[idx].allVertices[3*this.displayBuffers[idx].allIndexs[j+2]];
                    var y3 = this.displayBuffers[idx].allVertices[3*this.displayBuffers[idx].allIndexs[j+2]+1];
                    var z3 = this.displayBuffers[idx].allVertices[3*this.displayBuffers[idx].allIndexs[j+2]+2];
                    var mid_x = (x1+x2+x3)/3.0;
                    var mid_y = (y1+y2+y3)/3.0;
                    var mid_z = (z1+z2+z3)/3.0;
                    var proj = mid_z; // FIXME - Need to calculate a projection along camera z-axis.
                    var projP = vec3.create([mid_x,mid_y,mid_z]);
                    mat4.multiplyVec3(theMatrix,projP);
                    sortThings.push(new SortThing(projP[2],this.displayBuffers[idx].allIndexs[j],this.displayBuffers[idx].allIndexs[j+1],this.displayBuffers[idx].allIndexs[j+2]));
                }
                sortThings.sort(sortIndicesByProj);
                var allIndexs = [];
                var maxInd = 0;
                for(var j=0;j<sortThings.length;j++){
                    allIndexs.push(sortThings[j].id1);
                    allIndexs.push(sortThings[j].id2);
                    allIndexs.push(sortThings[j].id3);
                    if(sortThings[j].id1>maxInd) maxInd = sortThings[j].id1;
                    if(sortThings[j].id2>maxInd) maxInd = sortThings[j].id2;
                    if(sortThings[j].id3>maxInd) maxInd = sortThings[j].id3;
                }
                //console.log("maxInd: "+maxInd);
                //console.log("this.displayBuffers[idx].allNormals.length: "+this.displayBuffers[idx].allNormals.length/3);
                //console.log("this.displayBuffers[idx].allVertices.length: "+this.displayBuffers[idx].allVertices.length/3);
                //console.log("this.displayBuffers[idx].allColours.length: "+this.displayBuffers[idx].allColours.length/3);
                this.gl.useProgram(this.shaderProgram);
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].allTriangleVertexNormalBuffer);
                this.gl.vertexAttribPointer(this.shaderProgram.vertexNormalAttribute, 3, this.gl.FLOAT, false, 0, 0);

                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].allTriangleVertexPositionBuffer);
                this.gl.vertexAttribPointer(this.shaderProgram.vertexPositionAttribute, 3, this.gl.FLOAT, false, 0, 0);

                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[idx].allTriangleVertexColourBuffer);
                this.gl.vertexAttribPointer(this.shaderProgram.vertexColourAttribute, 4, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER,this.displayBuffers[idx].allIndexsBuffer);
                // FIXME - hmm, one is /3, the other is not ....
                if(this.ext){
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(allIndexs), this.gl.STATIC_DRAW);
                    this.gl.drawElements(this.gl.TRIANGLES, allIndexs.length, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(allIndexs), this.gl.STATIC_DRAW);
                    this.gl.drawElements(this.gl.TRIANGLES, allIndexs.length/3, this.gl.UNSIGNED_SHORT, 0);
                }
            }
        }
    }
}

MGWebGL.prototype.drawImagesAndText = function() {
    // Now the "see-through" primitives: images and text.
    for (var idx = 0; idx < this.displayBuffers.length; idx++){

        if(!this.displayBuffers[idx].visible) {
            continue;
        }

        var bufferTypes = this.displayBuffers[idx].bufferTypes;

        var triangleVertexNormalBuffer = this.displayBuffers[idx].triangleVertexNormalBuffer;
        var triangleVertexPositionBuffer = this.displayBuffers[idx].triangleVertexPositionBuffer;
        var triangleVertexIndexBuffer = this.displayBuffers[idx].triangleVertexIndexBuffer;
        var triangleColourBuffer = this.displayBuffers[idx].triangleColourBuffer;

        var triangleIndexs = this.displayBuffers[idx].triangleIndexs;
        var triangleVertices = this.displayBuffers[idx].triangleVertices;
        var triangleColours = this.displayBuffers[idx].triangleColours;
        var triangleNormals = this.displayBuffers[idx].triangleNormals;

        var primitiveSizes = this.displayBuffers[idx].primitiveSizes;

        this.gl.useProgram(this.shaderProgramImages);
        this.gl.depthFunc(this.gl.ALWAYS);
        this.setMatrixUniforms(this.shaderProgramImages);

        this.gl.enableVertexAttribArray(this.shaderProgramImages.vertexTextureAttribute);
        this.gl.disableVertexAttribArray(this.shaderProgramImages.vertexColourAttribute);
        this.gl.vertexAttrib4f(this.shaderProgramImages.vertexColourAttribute,1.0,1.0,0.0,1.0);

        for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
            if(bufferTypes[j]=="IMAGES"){
                var buffer;
                buffer = this.imageBuffer;

                if(!(this.displayBuffers[idx].texture)){
                    this.displayBuffers[idx].texture = initTextures(this.gl,this.displayBuffers[idx].supplementary["imgsrc"][0]);
                }
                var scaleImage = true;
                if(typeof(this.gl,this.displayBuffers[idx].supplementary["vert_tri_2d"])!=="undefined"){
                    var tempMVMatrix = mat4.create([this.mvMatrix[0],this.mvMatrix[1],this.mvMatrix[2],this.mvMatrix[3],this.mvMatrix[4],this.mvMatrix[5],this.mvMatrix[6],this.mvMatrix[7],this.mvMatrix[8],this.mvMatrix[9],this.mvMatrix[10],this.mvMatrix[11], (-24.0 + this.displayBuffers[idx].supplementary["vert_tri_2d"][0][0] * 48.0)*this.zoom,(-24.0 + this.displayBuffers[idx].supplementary["vert_tri_2d"][0][1] * 48.0)*this.zoom,-500.0,1.0 ]);
                    this.gl.uniformMatrix4fv(this.shaderProgramImages.mvMatrixUniform, false, tempMVMatrix);
                    scaleImage = false;
                }

                this.gl.bindTexture(this.gl.TEXTURE_2D, this.displayBuffers[idx].texture);
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexTextureBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramImages.vertexTextureAttribute, buffer.triangleVertexTextureBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexNormalBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramImages.vertexNormalAttribute, buffer.triangleVertexNormalBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexPositionBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramImages.vertexPositionAttribute, buffer.triangleVertexPositionBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, buffer.triangleVertexIndexBuffer[0]);
                // FIXME - And loop here
                var theOffSet = new glMatrixArrayType(3);
                for(ishape=0;ishape<triangleVertices[j].length/3;ishape++){
                    theOffSet[0] = triangleVertices[j][ishape*3];
                    theOffSet[1] = triangleVertices[j][ishape*3+1];
                    theOffSet[2] = triangleVertices[j][ishape*3+2];
                    this.gl.uniform3fv(this.shaderProgramImages.offset, theOffSet);
                    if(scaleImage){
                        this.gl.uniform1f(this.shaderProgramImages.size, primitiveSizes[j][ishape]);
                    }else{
                        this.gl.uniform1f(this.shaderProgramImages.size, primitiveSizes[j][ishape]*this.zoom);
                    }

                    var thisCol = [ triangleColours[j][ishape*4], triangleColours[j][ishape*4+1], triangleColours[j][ishape*4+2], triangleColours[j][ishape*4+3] ];
                    this.gl.vertexAttrib4f(this.shaderProgramImages.vertexColourAttribute,triangleColours[j][ishape*4], triangleColours[j][ishape*4+1], triangleColours[j][ishape*4+2], triangleColours[j][ishape*4+3]);

                    if(this.ext){
                        this.gl.drawElements(this.gl.TRIANGLE_FAN, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_INT, 0);
                    }else{
                        this.gl.drawElements(this.gl.TRIANGLE_FAN, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_SHORT, 0);
                    }
                }
                if(typeof(this.gl,this.displayBuffers[idx].supplementary["vert_tri_2d"])!=="undefined"){
                    this.setMatrixUniforms(this.shaderProgramImages);
                }
            }
        }

        for (var j = 0; j < triangleVertexIndexBuffer.length; j++){
            if(bufferTypes[j]=="TEXT"){
                var buffer;
                buffer = this.imageBuffer;

                var font = this.displayBuffers[idx].supplementary["font"][0][0];
                var fnsize = font.match(/^\d+|\d+\b|\d+(?=\w)/g)[0];
                if(!(this.displayBuffers[idx].texture)){
                    //this.displayBuffers[idx].texture = initStringTextures(this.gl,"Hello World !!!!");
                    console.log(this.displayBuffers[idx].supplementary["imgsrc"][0][0]);
                    var tex_size = {};
                    console.log(font);
                    this.displayBuffers[idx].texture = initStringTextures(this.gl,this.displayBuffers[idx].supplementary["imgsrc"][0][0],tex_size,font);
                    console.log(tex_size);
                    this.displayBuffers[idx].tex_size = tex_size;
                }
                for(var iv=0;iv<this.imageVertices.length;iv+=3){
                    var vold = vec3.create([this.imageVertices[iv]*this.displayBuffers[idx].tex_size["width"]/this.displayBuffers[idx].tex_size["height"],this.imageVertices[iv+1],this.imageVertices[iv+2]]);
                    var vnew = vec3.create();
                    mat4.multiplyVec3(invMat,vold,vnew);
                    imageVertices[iv]   = vnew[0];
                    imageVertices[iv+1] = vnew[1];
                    imageVertices[iv+2] = vnew[2];
                }
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.imageBuffer.triangleVertexPositionBuffer[0]);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(imageVertices), this.gl.DYNAMIC_DRAW);
                var scaleImage = false;
                if(typeof(this.gl,this.displayBuffers[idx].supplementary["vert_tri_2d"])!=="undefined"){
                    var tempMVMatrix = mat4.create([this.mvMatrix[0],this.mvMatrix[1],this.mvMatrix[2],this.mvMatrix[3],this.mvMatrix[4],this.mvMatrix[5],this.mvMatrix[6],this.mvMatrix[7],this.mvMatrix[8],this.mvMatrix[9],this.mvMatrix[10],this.mvMatrix[11], (-24.0 + this.displayBuffers[idx].supplementary["vert_tri_2d"][0][0] * 48.0)*this.zoom,(-24.0 + this.displayBuffers[idx].supplementary["vert_tri_2d"][0][1] * 48.0)*this.zoom,-500.0,1.0 ]);
                    this.gl.uniformMatrix4fv(this.shaderProgramImages.mvMatrixUniform, false, tempMVMatrix);
                    scaleImage = false;
                }

                this.gl.bindTexture(this.gl.TEXTURE_2D, this.displayBuffers[idx].texture);
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexTextureBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramImages.vertexTextureAttribute, buffer.triangleVertexTextureBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexNormalBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramImages.vertexNormalAttribute, buffer.triangleVertexNormalBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer.triangleVertexPositionBuffer[0]);
                this.gl.vertexAttribPointer(this.shaderProgramImages.vertexPositionAttribute, buffer.triangleVertexPositionBuffer[0].itemSize, this.gl.FLOAT, false, 0, 0);
                this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, buffer.triangleVertexIndexBuffer[0]);
                // FIXME - And loop here ?!?
                var theOffSet = new glMatrixArrayType(3);
                for(ishape=0;ishape<triangleVertices[j].length/3;ishape++){
                    theOffSet[0] = triangleVertices[j][ishape*3];
                    theOffSet[1] = triangleVertices[j][ishape*3+1];
                    theOffSet[2] = triangleVertices[j][ishape*3+2];
                    this.gl.uniform3fv(this.shaderProgramImages.offset, theOffSet);
                    if(scaleImage){
                        this.gl.uniform1f(this.shaderProgramImages.size, primitiveSizes[j][ishape]);
                    }else{
                        this.gl.uniform1f(this.shaderProgramImages.size, primitiveSizes[j][ishape]*this.zoom*fnsize/this.canvas.height*48*getDeviceScale());
                    }

                    var thisCol = [ triangleColours[j][ishape*4], triangleColours[j][ishape*4+1], triangleColours[j][ishape*4+2], triangleColours[j][ishape*4+3] ];
                    this.gl.vertexAttrib4f(this.shaderProgramImages.vertexColourAttribute,triangleColours[j][ishape*4], triangleColours[j][ishape*4+1], triangleColours[j][ishape*4+2], triangleColours[j][ishape*4+3]);

                    if(this.ext){
                        this.gl.drawElements(this.gl.TRIANGLE_FAN, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_INT, 0);
                    }else{
                        this.gl.drawElements(this.gl.TRIANGLE_FAN, buffer.triangleVertexIndexBuffer[0].numItems, this.gl.UNSIGNED_SHORT, 0);
                    }
                }
                if(typeof(this.gl,this.displayBuffers[idx].supplementary["vert_tri_2d"])!=="undefined"){
                    this.setMatrixUniforms(this.shaderProgramImages);
                }
            }
        }

        this.gl.disableVertexAttribArray(this.shaderProgramImages.vertexTextureAttribute);
        this.gl.enableVertexAttribArray(this.shaderProgramImages.vertexColourAttribute);
        this.gl.depthFunc(this.gl.LESS);
        
    }
}

MGWebGL.prototype.drawTextLabels = function(up,right) {
    // make sure we can render it even if it's not a power of 2
    this.gl.bindTexture(this.gl.TEXTURE_2D, this.textTex);
    this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_MIN_FILTER, this.gl.LINEAR);
    this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_WRAP_S, this.gl.CLAMP_TO_EDGE);
    this.gl.texParameteri(this.gl.TEXTURE_2D, this.gl.TEXTURE_WRAP_T, this.gl.CLAMP_TO_EDGE);

    if(this.textLabels.length>0){
        if(typeof this.displayBuffers[0].textPositionBuffer === "undefined") {
            this.initTextBuffers();
            this.displayBuffers[0].textIndexs = [];
            this.displayBuffers[0].textTexCoords = [];
            this.displayBuffers[0].textTexCoords = this.displayBuffers[0].textTexCoords.concat([ 0,1, 1,1, 1,0]);
            this.displayBuffers[0].textTexCoords = this.displayBuffers[0].textTexCoords.concat([ 0,1, 1,0, 0,0]);
            this.displayBuffers[0].textIndexs = this.displayBuffers[0].textIndexs.concat([0,1,2]);
            this.displayBuffers[0].textIndexs = this.displayBuffers[0].textIndexs.concat([3,4,5]);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textTexCoordBuffer);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textTexCoords), this.gl.STATIC_DRAW);
            this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexTextureAttribute, 2, this.gl.FLOAT, false, 0, 0);

            this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER,this.displayBuffers[0].textIndexesBuffer);
            if(this.ext){
                this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(this.displayBuffers[0].textIndexs), this.gl.STATIC_DRAW);
            }else{
                this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.displayBuffers[0].textIndexs), this.gl.STATIC_DRAW);
            }
            this.makeTextCanvas("Fluffy", 512, 32, textColour);
            this.gl.texImage2D(this.gl.TEXTURE_2D, 0, this.gl.RGBA, this.gl.RGBA, this.gl.UNSIGNED_BYTE, this.textCtx.canvas);

            this.displayBuffers[0].textNormals = [];
            this.displayBuffers[0].textColours = [];
            this.displayBuffers[0].textNormals = this.displayBuffers[0].textNormals.concat([ 0,0,1, 0,0,1, 0,0,1]);
            this.displayBuffers[0].textNormals = this.displayBuffers[0].textNormals.concat([ 0,0,1, 0,0,1, 0,0,1]);
            this.displayBuffers[0].textColours = this.displayBuffers[0].textColours.concat([ 1,0,0,1, 1,0,0,1, 1,0,0,1]);
            this.displayBuffers[0].textColours = this.displayBuffers[0].textColours.concat([ 1,0,0,1, 1,0,0,1, 1,0,0,1]);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textNormalBuffer);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textNormals), this.gl.STATIC_DRAW);
            this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexNormalAttribute, 3, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textColourBuffer);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textColours), this.gl.STATIC_DRAW);
            this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexColourAttribute, 4, this.gl.FLOAT, false, 0, 0);
        }

        this.gl.useProgram(this.shaderProgramTextBackground);
        this.gl.enableVertexAttribArray(this.shaderProgramTextBackground.vertexTextureAttribute);
        //this.gl.disableVertexAttribArray(this.shaderProgram.vertexColourAttribute);
        //this.gl.disableVertexAttribArray(this.shaderProgram.vertexNormalAttribute);
        this.setMatrixUniforms(this.shaderProgramTextBackground);
        this.gl.depthFunc(this.gl.ALWAYS);

        var textColour = "black";
        var y = this.background_colour[0]*0.299 + this.background_colour[1]*0.587 + this.background_colour[2]*0.114;
        if(y<0.5){
            textColour = "white";
        }
        var idx = 0;
        this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER,this.displayBuffers[0].textIndexesBuffer);
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textPositionBuffer);
        this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexPositionAttribute, 3, this.gl.FLOAT, false, 0, 0);

        for(var itl=0;itl<this.textLabels.length;itl++){
            var thisVis = this.displayBuffers[itl].visible;
            if(!thisVis){
                continue;
            }
            for(var jtl=0;jtl<this.textLabels[itl].length;jtl++){
                this.displayBuffers[0].textVertices = [];

                this.makeTextCanvas(this.textLabels[itl][jtl].label, 512, 32, textColour);
                var textWidth  = this.textCtx.canvas.width;
                var textHeight = this.textCtx.canvas.height;

                // This is slow on FF on SL6.7
                if(typeof(this.textLabels[itl][jtl].imgData)==="undefined"){
                    this.textLabels[itl][jtl].imgData = this.textCtx.getImageData(0,0,512, 32);
                }
                this.gl.texSubImage2D(this.gl.TEXTURE_2D, 0, 0, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, this.textLabels[itl][jtl].imgData);

                var x = this.textLabels[itl][jtl].x;
                var y = this.textLabels[itl][jtl].y;
                var z = this.textLabels[itl][jtl].z;
                var tSizeX = 2.0 * this.textCtx.canvas.width/this.textCtx.canvas.height * this.zoom;
                var tSizeY = 2.0 * this.zoom;
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeX*right[0],y+tSizeX*right[1],z+tSizeX*right[2], x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2]]);
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2], x+tSizeY*up[0],y+tSizeY*up[1],z+tSizeY*up[2]]);

                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textVertices), this.gl.DYNAMIC_DRAW);

                if(this.ext){
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_SHORT, 0);
                }

            }
        }
        //this.gl.enableVertexAttribArray(this.shaderProgram.vertexColourAttribute);
        //this.gl.enableVertexAttribArray(this.shaderProgram.vertexNormalAttribute);
        this.gl.disableVertexAttribArray(this.shaderProgramTextBackground.vertexTextureAttribute);
        this.gl.depthFunc(this.gl.LESS)
    }
}

MGWebGL.prototype.drawClickedAtoms = function(up,right) {
    if(this.clickedAtoms.length>0){

        // FIXME - bugarrific - using first buffer always is wrong!! Should be text labels per object!
        if(typeof this.displayBuffers[0].textPositionBuffer === "undefined") {
            this.initTextBuffers();
            this.displayBuffers[0].textIndexs = [];
            this.displayBuffers[0].textTexCoords = [];
            this.displayBuffers[0].textTexCoords = this.displayBuffers[0].textTexCoords.concat([ 0,1, 1,1, 1,0]);
            this.displayBuffers[0].textTexCoords = this.displayBuffers[0].textTexCoords.concat([ 0,1, 1,0, 0,0]);
            this.displayBuffers[0].textIndexs = this.displayBuffers[0].textIndexs.concat([0,1,2]);
            this.displayBuffers[0].textIndexs = this.displayBuffers[0].textIndexs.concat([3,4,5]);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textTexCoordBuffer);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textTexCoords), this.gl.STATIC_DRAW);
            this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexTextureAttribute, 2, this.gl.FLOAT, false, 0, 0);

            this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER,this.displayBuffers[0].textIndexesBuffer);
            if(this.ext){
                this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(this.displayBuffers[0].textIndexs), this.gl.STATIC_DRAW);
            }else{
                this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.displayBuffers[0].textIndexs), this.gl.STATIC_DRAW);
            }
            this.makeTextCanvas("Fluffy", 512, 32, textColour);
            this.gl.texImage2D(this.gl.TEXTURE_2D, 0, this.gl.RGBA, this.gl.RGBA, this.gl.UNSIGNED_BYTE, this.textCtx.canvas);

            this.displayBuffers[0].textNormals = [];
            this.displayBuffers[0].textColours = [];
            this.displayBuffers[0].textNormals = this.displayBuffers[0].textNormals.concat([ 0,0,1, 0,0,1, 0,0,1]);
            this.displayBuffers[0].textNormals = this.displayBuffers[0].textNormals.concat([ 0,0,1, 0,0,1, 0,0,1]);
            this.displayBuffers[0].textColours = this.displayBuffers[0].textColours.concat([ 1,0,0,1, 1,0,0,1, 1,0,0,1]);
            this.displayBuffers[0].textColours = this.displayBuffers[0].textColours.concat([ 1,0,0,1, 1,0,0,1, 1,0,0,1]);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textNormalBuffer);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textNormals), this.gl.STATIC_DRAW);
            this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexNormalAttribute, 3, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textColourBuffer);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textColours), this.gl.STATIC_DRAW);
            this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexColourAttribute, 4, this.gl.FLOAT, false, 0, 0);
        }


        this.gl.useProgram(this.shaderProgramLines);
        var lineVertices = [];
        var lineColours = [];
        var lineIndexs = [];
        var idx = 0;

        for(var iat=0;iat<this.clickedAtoms.length;iat++){

            for(var jat=1;jat<this.clickedAtoms[iat].length;jat++){
                var x1 = this.clickedAtoms[iat][jat-1].x;
                var y1 = this.clickedAtoms[iat][jat-1].y;
                var z1 = this.clickedAtoms[iat][jat-1].z;
                var x2 = this.clickedAtoms[iat][jat].x;
                var y2 = this.clickedAtoms[iat][jat].y;
                var z2 = this.clickedAtoms[iat][jat].z;
                lineVertices.push(x1);
                lineVertices.push(y1);
                lineVertices.push(z1);
                lineVertices.push(x2);
                lineVertices.push(y2);
                lineVertices.push(z2);
                lineColours.push(1.0);
                lineColours.push(0.0);
                lineColours.push(0.0);
                lineColours.push(1.0);
                lineColours.push(1.0);
                lineColours.push(0.0);
                lineColours.push(0.0);
                lineColours.push(1.0);
                lineIndexs.push(idx++);
                lineIndexs.push(idx++);
            }
        }
        if(lineIndexs.length>0){
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].clickLinePositionBuffer);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(lineVertices), this.gl.DYNAMIC_DRAW);
            this.gl.vertexAttribPointer(this.shaderProgramLines.vertexPositionAttribute, 3, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].clickLineColourBuffer);
            this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(lineColours), this.gl.DYNAMIC_DRAW);
            this.gl.vertexAttribPointer(this.shaderProgramLines.vertexColourAttribute, 4, this.gl.FLOAT, false, 0, 0);
            this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.displayBuffers[0].clickLineIndexesBuffer);
            if(this.ext){
                this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(lineIndexs), this.gl.DYNAMIC_DRAW);
                this.gl.drawElements(this.gl.LINES, lineIndexs.length, this.gl.UNSIGNED_INT, 0);
            }else{
                this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(lineIndexs), this.gl.DYNAMIC_DRAW);
                this.gl.drawElements(this.gl.LINES,lineIndexs.length, this.gl.UNSIGNED_SHORT, 0);
            }
        }

        this.gl.useProgram(this.shaderProgramTextBackground);
        this.gl.enableVertexAttribArray(this.shaderProgramTextBackground.vertexTextureAttribute);
        this.setMatrixUniforms(this.shaderProgramTextBackground);

        this.gl.depthFunc(this.gl.ALWAYS);

        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textNormalBuffer);
        this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexNormalAttribute, 3, this.gl.FLOAT, false, 0, 0);
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textColourBuffer);
        this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexColourAttribute, 4, this.gl.FLOAT, false, 0, 0);

        this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER,this.displayBuffers[0].textIndexesBuffer);
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.displayBuffers[0].textPositionBuffer);
        this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexPositionAttribute, 3, this.gl.FLOAT, false, 0, 0);

        for(var iat=0;iat<this.clickedAtoms.length;iat++){
            var textColour = "black";
            var y = this.background_colour[0]*0.299 + this.background_colour[1]*0.587 + this.background_colour[2]*0.114;
            if(y<0.5){
                textColour = "white";
            }
            for(var jat=0;jat<this.clickedAtoms[iat].length;jat++){
                this.displayBuffers[0].textVertices = [];

                if(typeof(this.clickedAtoms[iat][jat].imgData)==="undefined"){
                    this.makeTextCanvas(this.clickedAtoms[iat][jat].label, 512, 32, textColour);
                    this.clickedAtoms[iat][jat].imgData = this.textCtx.getImageData(0,0,512, 32);
                }
                this.gl.texSubImage2D(this.gl.TEXTURE_2D, 0, 0, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, this.clickedAtoms[iat][jat].imgData);

                var x = this.clickedAtoms[iat][jat].x;
                var y = this.clickedAtoms[iat][jat].y;
                var z = this.clickedAtoms[iat][jat].z;
                var tSizeX = 2.0 * this.textCtx.canvas.width/this.textCtx.canvas.height * this.zoom;
                var tSizeY = 2.0 * this.zoom;
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeX*right[0],y+tSizeX*right[1],z+tSizeX*right[2], x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2]]);
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2], x+tSizeY*up[0],y+tSizeY*up[1],z+tSizeY*up[2]]);

                //console.log(this.displayBuffers[0].textVertices);
                //console.log(this.displayBuffers[0].textIndexs);
                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textVertices), this.gl.DYNAMIC_DRAW);

                if(this.ext){
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_SHORT, 0);
                }
            }
        }
        for(var iat=0;iat<this.clickedAtoms.length;iat++){

            // FIXME - This needs a tweak for at1-at2-at3-at4 dihedrals.
            for(var jat=1;jat<this.clickedAtoms[iat].length;jat++){
                this.displayBuffers[0].textVertices = [];
                var x1 = this.clickedAtoms[iat][jat-1].x;
                var y1 = this.clickedAtoms[iat][jat-1].y;
                var z1 = this.clickedAtoms[iat][jat-1].z;

                var x2 = this.clickedAtoms[iat][jat].x;
                var y2 = this.clickedAtoms[iat][jat].y;
                var z2 = this.clickedAtoms[iat][jat].z;

                var v1 = vec3.create([x1,y1,z1]);
                var v2 = vec3.create([x2,y2,z2]);

                var v1diffv2 = vec3.create();
                vec3.subtract(v1,v2,v1diffv2);
                var linesize = vec3.length(v1diffv2);

                var v1plusv2 = vec3.create();
                vec3.add(v1,v2,v1plusv2);
                var x = v1plusv2[0] * 0.5;
                var y = v1plusv2[1] * 0.5;
                var z = v1plusv2[2] * 0.5;

                if(typeof(this.clickedAtoms[iat][jat].lengthImgData)==="undefined"){
                    this.makeTextCanvas(linesize.toFixed(3), 512, 32, textColour);
                    this.clickedAtoms[iat][jat].lengthImgData = this.textCtx.getImageData(0,0,512, 32);
                }
                this.gl.texSubImage2D(this.gl.TEXTURE_2D, 0, 0, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, this.clickedAtoms[iat][jat].lengthImgData);

                var tSizeX = 2.0 * this.textCtx.canvas.width/this.textCtx.canvas.height * this.zoom;
                var tSizeY = 2.0 * this.zoom;
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeX*right[0],y+tSizeX*right[1],z+tSizeX*right[2], x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2]]);
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2], x+tSizeY*up[0],y+tSizeY*up[1],z+tSizeY*up[2]]);

                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textVertices), this.gl.DYNAMIC_DRAW);

                if(this.ext){
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_SHORT, 0);
                }
            }
            for(var jat=2;jat<this.clickedAtoms[iat].length;jat++){
                this.displayBuffers[0].textVertices = [];
                var x1 = this.clickedAtoms[iat][jat-2].x;
                var y1 = this.clickedAtoms[iat][jat-2].y;
                var z1 = this.clickedAtoms[iat][jat-2].z;

                var x2 = this.clickedAtoms[iat][jat-1].x;
                var y2 = this.clickedAtoms[iat][jat-1].y;
                var z2 = this.clickedAtoms[iat][jat-1].z;

                var x3 = this.clickedAtoms[iat][jat].x;
                var y3 = this.clickedAtoms[iat][jat].y;
                var z3 = this.clickedAtoms[iat][jat].z;

                var v1 = vec3.create([x1,y1,z1]);
                var v2 = vec3.create([x2,y2,z2]);
                var v3 = vec3.create([x3,y3,z3]);

                var v2diffv1 = vec3.create();
                vec3.subtract(v2,v1,v2diffv1);
                vec3.normalize(v2diffv1);

                var v2diffv3 = vec3.create();
                vec3.subtract(v2,v3,v2diffv3);
                vec3.normalize(v2diffv3);

                //console.log(vec3.dot(v2diffv1,v2diffv3));

                var angle = Math.acos(vec3.dot(v2diffv1,v2diffv3))*180.0/Math.PI;

                var x = x2-tSizeY*.5*up[0];
                var y = y2-tSizeY*.5*up[1];
                var z = z2-tSizeY*.5*up[2];

                var textWidth  = this.textCtx.canvas.width;
                var textHeight = this.textCtx.canvas.height;
                if(typeof(this.clickedAtoms[iat][jat].angleImgData)==="undefined"){
                    this.makeTextCanvas(angle.toFixed(1), 512, 32, textColour);
                    this.clickedAtoms[iat][jat].angleImgData = this.textCtx.getImageData(0,0,512, 32);
                }
                this.gl.texSubImage2D(this.gl.TEXTURE_2D, 0, 0, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, this.clickedAtoms[iat][jat].angleImgData);

                var tSizeX = 2.0 * this.textCtx.canvas.width/this.textCtx.canvas.height * this.zoom;
                var tSizeY = 2.0 * this.zoom;
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeX*right[0],y+tSizeX*right[1],z+tSizeX*right[2], x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2]]);
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2], x+tSizeY*up[0],y+tSizeY*up[1],z+tSizeY*up[2]]);

                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textVertices), this.gl.DYNAMIC_DRAW);

                if(this.ext){
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_SHORT, 0);
                }
            }
            for(var jat=3;jat<this.clickedAtoms[iat].length;jat++){
                this.displayBuffers[0].textVertices = [];
                var x1 = this.clickedAtoms[iat][jat-3].x;
                var y1 = this.clickedAtoms[iat][jat-3].y;
                var z1 = this.clickedAtoms[iat][jat-3].z;

                var x2 = this.clickedAtoms[iat][jat-2].x;
                var y2 = this.clickedAtoms[iat][jat-2].y;
                var z2 = this.clickedAtoms[iat][jat-2].z;

                var x3 = this.clickedAtoms[iat][jat-1].x;
                var y3 = this.clickedAtoms[iat][jat-1].y;
                var z3 = this.clickedAtoms[iat][jat-1].z;

                var x4 = this.clickedAtoms[iat][jat].x;
                var y4 = this.clickedAtoms[iat][jat].y;
                var z4 = this.clickedAtoms[iat][jat].z;

                var v1 = vec3.create([x1,y1,z1]);
                var v2 = vec3.create([x2,y2,z2]);
                var v3 = vec3.create([x3,y3,z3]);
                var v4 = vec3.create([x4,y4,z4]);

                var v3plusv2 = vec3.create();
                vec3.add(v3,v2,v3plusv2);
                var x = v3plusv2[0] * 0.5-tSizeY*.5*up[0];
                var y = v3plusv2[1] * 0.5-tSizeY*.5*up[0];
                var z = v3plusv2[2] * 0.5-tSizeY*.5*up[0];

                var angle = DihedralAngle(v1,v2,v3,v4)*180.0/Math.PI;

                var textWidth  = this.textCtx.canvas.width;
                var textHeight = this.textCtx.canvas.height;
                //this.gl.texImage2D(this.gl.TEXTURE_2D, 0, this.gl.RGBA, this.gl.RGBA, this.gl.UNSIGNED_BYTE, this.textCtx.canvas);
                if(typeof(this.clickedAtoms[iat][jat].dihedralImgData)==="undefined"){
                    this.makeTextCanvas(angle.toFixed(1), 512, 32, textColour);
                    this.clickedAtoms[iat][jat].dihedralImgData = this.textCtx.getImageData(0,0,512, 32);
                }
                this.gl.texSubImage2D(this.gl.TEXTURE_2D, 0, 0, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, this.clickedAtoms[iat][jat].dihedralImgData);

                var tSizeX = 2.0 * this.textCtx.canvas.width/this.textCtx.canvas.height * this.zoom;
                var tSizeY = 2.0 * this.zoom;
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeX*right[0],y+tSizeX*right[1],z+tSizeX*right[2], x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2]]);
                this.displayBuffers[0].textVertices = this.displayBuffers[0].textVertices.concat([x,y,z, x+tSizeY*up[0]+tSizeX*right[0],y+tSizeY*up[1]+tSizeX*right[1],z+tSizeY*up[2]+tSizeX*right[2], x+tSizeY*up[0],y+tSizeY*up[1],z+tSizeY*up[2]]);

                this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.displayBuffers[0].textVertices), this.gl.DYNAMIC_DRAW);

                if(this.ext){
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_INT, 0);
                }else{
                    this.gl.drawElements(this.gl.TRIANGLES, this.displayBuffers[0].textIndexs.length, this.gl.UNSIGNED_SHORT, 0);
                }
            }
            //this.gl.enableVertexAttribArray(this.shaderProgram.vertexColourAttribute);
            //this.gl.enableVertexAttribArray(this.shaderProgram.vertexNormalAttribute);

        }

        this.gl.disableVertexAttribArray(this.shaderProgramTextBackground.vertexTextureAttribute);
        this.gl.depthFunc(this.gl.LESS)

    }
}

MGWebGL.prototype.doClick = function(event,self) {
    if(!self.mouseMoved){
        var x;
        var y;
        var e = event;
        if (e.pageX || e.pageY) { 
            x = e.pageX;
            y = e.pageY;
        }
        else { 
            x = e.clientX ; 
            y = e.clientY ; 
        } 

        var offset = getOffsetRect(self.div);
        
        x -= offset.left;
        y -= offset.top;
        x *= getDeviceScale();
        y *= getDeviceScale();


        //document.getElementById("info").innerHTML="Click: "+event.x+" "+event.y;
        var invQuat = quat4.create();
        quat4.inverse(self.myQuat,invQuat);
        var theMatrix = quat4.toMat4(invQuat);
        var ratio =  1.0*self.gl.viewportWidth/self.gl.viewportHeight;
        var minX = (-24.*ratio*self.zoom);
        var maxX =  (24.*ratio*self.zoom);
        var minY = (-24.*self.zoom);
        var maxY =  (24.*self.zoom);
        var fracX = 1.0*x/self.gl.viewportWidth;
        var fracY = 1.0*(y)/self.gl.viewportHeight;
        var theX = minX + fracX*(maxX - minX);
        var theY = maxY - fracY*(maxY - minY);
        var frontPos = vec3.create([theX,theY,0.000001]);
        var backPos  = vec3.create([theX,theY,1000.0]);
        mat4.multiplyVec3(theMatrix,frontPos);
        mat4.multiplyVec3(theMatrix,backPos);
        vec3.subtract(frontPos,self.origin);
        vec3.subtract(backPos,self.origin);

        var mindist = 100000.;
        var mint = -100000.;
        var minidx = -1;
        var minj = -1;
        var clickTol = 0.65*this.zoom;

        /*
           var theForwardMatrix = quat4.toMat4(self.myQuat);

           var npass = 0;
           var npass0 = 0;
           var npass1 = 0;
           var ntest = 0;
         */

        for (var idx = 0; idx < self.displayBuffers.length; idx++){
            if(!self.displayBuffers[idx].visible) {
                continue;
            }
            //console.log( self.displayBuffers[idx].atoms[0]);
            for (var j = 0; j < self.displayBuffers[idx].atoms.length; j++){
                //var p = vec3.create([79.109,59.437,27.316]);
                var atx = self.displayBuffers[idx].atoms[j].x;
                var aty = self.displayBuffers[idx].atoms[j].y;
                var atz = self.displayBuffers[idx].atoms[j].z;
                //console.log(atx+" "+aty+" "+atz);
                var p = vec3.create([atx,aty,atz]);

                /*
                   ntest++;
                   var ptrans = vec3.create([atx,aty,atz]);
                   mat4.multiplyVec3(theForwardMatrix,ptrans);
                   vec3.add(ptrans,self.origin);
                   var dotRes0 = -ptrans[2] -(self.gl_clipPlane0[3]+500);
                   var dotRes1 =  ptrans[2] +(self.gl_clipPlane1[3]-500);
                // Check if clipped, need plane equations without the 500 offset and to take the origin and rotation into account...
                // Hmm, doesn't seem to work. Especially after rotation .... Why am I so dim....
                //res = Nx*x + Ny*y + Nz*z + D

                if(dotRes0<0) {
                npass0++;
                }
                if(dotRes1>0) {
                npass1++;
                }
                if(dotRes1>0&&dotRes0<0) {
                npass++;
                }
                console.log((-(self.gl_clipPlane0[3]+500))+" "+(self.gl_clipPlane1[3]-500)+" "+dotRes0+" "+dotRes1);
                if(dotRes0>0||dotRes1<0) {
                continue;
                }
                 */
                var dpl = DistanceBetweenPointAndLine(frontPos,backPos,p);
                if(dpl[0]<clickTol&&dpl[1]>mint){
                    minidx = idx;
                    minj = j;
                    mindist = dpl[0];
                    mint = dpl[1];
                }
            }
        }
        //console.log(npass+" "+npass0+" "+npass1+" "+ntest);
        if(minidx>-1){
            //console.log("mint: "+mint+", mindist: "+mindist+", minidx: "+minidx+", minj: "+minj);
            var theAtom = {};
            theAtom.x =  self.displayBuffers[minidx].atoms[minj].x;
            theAtom.y =  self.displayBuffers[minidx].atoms[minj].y;
            theAtom.z =  self.displayBuffers[minidx].atoms[minj].z;
            theAtom.charge =  self.displayBuffers[minidx].atoms[minj].charge;
            theAtom.label =  self.displayBuffers[minidx].atoms[minj].label;
            var atx = theAtom.x;
            var aty = theAtom.y;
            var atz = theAtom.z;
            var label = theAtom.label;
            console.log(atx+" "+aty+" "+atz+" "+label);

            if(event.altKey){
                self.origin = [-atx,-aty,-atz];
                self.reContourMaps();
                self.drawScene();
                return;
            }

            if(self.clickedAtoms.length==0||(self.clickedAtoms[self.clickedAtoms.length-1].length>1&&!event.shiftKey)){
                self.clickedAtoms.push([]);
                self.clickedAtoms[self.clickedAtoms.length-1].push(theAtom);
            } else {
                self.clickedAtoms[self.clickedAtoms.length-1].push(theAtom);
            }

        }
        //console.log(dpl);
    }
    self.drawScene();
}

MGWebGL.prototype.doWheel = function(event) {
    var factor;
    if( event.deltaY>0){
        factor = 1. + 5/50.;
    } else {
        factor = 1. - 5/50.;
    }
    this.zoom = this.zoom * factor;
    if(this.zoom < .01){
        this.zoom = 0.01;
    }
    this.drawScene();
    return;
}

MGWebGL.prototype.linesToThickLines = function(axesVertices,axesColours,size) {
    var axesNormals = [];
    var axesVertices_new = [];
    var axesColours_new = [];
    var axesIndexs_new = [];
    var axesIdx_new = 0;

    for(var il=0;il<axesVertices.length;il+=6){
        axesColours_new.push(axesColours[4*il/3]);
        axesColours_new.push(axesColours[4*il/3+1]);
        axesColours_new.push(axesColours[4*il/3+2]);
        axesColours_new.push(axesColours[4*il/3+3]);
        axesColours_new.push(axesColours[4*il/3]);
        axesColours_new.push(axesColours[4*il/3+1]);
        axesColours_new.push(axesColours[4*il/3+2]);
        axesColours_new.push(axesColours[4*il/3+3]);
        axesColours_new.push(axesColours[4*il/3+4]);
        axesColours_new.push(axesColours[4*il/3+5]);
        axesColours_new.push(axesColours[4*il/3+6]);
        axesColours_new.push(axesColours[4*il/3+7]);
        axesColours_new.push(axesColours[4*il/3]);
        axesColours_new.push(axesColours[4*il/3+1]);
        axesColours_new.push(axesColours[4*il/3+2]);
        axesColours_new.push(axesColours[4*il/3+3]);
        axesColours_new.push(axesColours[4*il/3+4]);
        axesColours_new.push(axesColours[4*il/3+5]);
        axesColours_new.push(axesColours[4*il/3+6]);
        axesColours_new.push(axesColours[4*il/3+7]);
        axesColours_new.push(axesColours[4*il/3+4]);
        axesColours_new.push(axesColours[4*il/3+5]);
        axesColours_new.push(axesColours[4*il/3+6]);
        axesColours_new.push(axesColours[4*il/3+7]);

        axesVertices_new.push(axesVertices[il]);
        axesVertices_new.push(axesVertices[il+1]);
        axesVertices_new.push(axesVertices[il+2]);
        axesVertices_new.push(axesVertices[il]);
        axesVertices_new.push(axesVertices[il+1]);
        axesVertices_new.push(axesVertices[il+2]);
        axesVertices_new.push(axesVertices[il+3]);
        axesVertices_new.push(axesVertices[il+4]);
        axesVertices_new.push(axesVertices[il+5]);
        axesNormals.push(axesVertices[il+3]-axesVertices[il]);
        axesNormals.push(axesVertices[il+4]-axesVertices[il+1]);
        axesNormals.push(axesVertices[il+5]-axesVertices[il+2]);
        var d = Math.sqrt(axesNormals[axesNormals.length-1-2]*axesNormals[axesNormals.length-1-2]+axesNormals[axesNormals.length-1-1]*axesNormals[axesNormals.length-1-1]+axesNormals[axesNormals.length-1-0]*axesNormals[axesNormals.length-1-0]);
        if(d>1e-8){
            axesNormals[axesNormals.length-1-2]   *= size/d;
            axesNormals[axesNormals.length-1-1] *= size/d;
            axesNormals[axesNormals.length-1] *= size/d;
        }
        
        axesNormals.push(-(axesVertices[il+3]-axesVertices[il]));
        axesNormals.push(-(axesVertices[il+4]-axesVertices[il+1]));
        axesNormals.push(-(axesVertices[il+5]-axesVertices[il+2]));
        if(d>1e-8){
            axesNormals[axesNormals.length-1-2]   *= size/d;
            axesNormals[axesNormals.length-1-1] *= size/d;
            axesNormals[axesNormals.length-1] *= size/d;
        }
        axesNormals.push(-(axesVertices[il+3]-axesVertices[il]));
        axesNormals.push(-(axesVertices[il+4]-axesVertices[il+1]));
        axesNormals.push(-(axesVertices[il+5]-axesVertices[il+2]));
        if(d>1e-8){
            axesNormals[axesNormals.length-1-2]   *= size/d;
            axesNormals[axesNormals.length-1-1] *= size/d;
            axesNormals[axesNormals.length-1] *= size/d;
        }
        axesVertices_new.push(axesVertices[il]);
        axesVertices_new.push(axesVertices[il+1]);
        axesVertices_new.push(axesVertices[il+2]);
        axesVertices_new.push(axesVertices[il+3]);
        axesVertices_new.push(axesVertices[il+4]);
        axesVertices_new.push(axesVertices[il+5]);
        axesVertices_new.push(axesVertices[il+3]);
        axesVertices_new.push(axesVertices[il+4]);
        axesVertices_new.push(axesVertices[il+5]);
        axesNormals.push(axesVertices[il+3]-axesVertices[il]);
        axesNormals.push(axesVertices[il+4]-axesVertices[il+1]);
        axesNormals.push(axesVertices[il+5]-axesVertices[il+2]);
        if(d>1e-8){
            axesNormals[axesNormals.length-1-2]   *= size/d;
            axesNormals[axesNormals.length-1-1] *= size/d;
            axesNormals[axesNormals.length-1] *= size/d;
        }
        axesNormals.push(axesVertices[il+3]-axesVertices[il]);
        axesNormals.push(axesVertices[il+4]-axesVertices[il+1]);
        axesNormals.push(axesVertices[il+5]-axesVertices[il+2]);
        if(d>1e-8){
            axesNormals[axesNormals.length-1-2]   *= size/d;
            axesNormals[axesNormals.length-1-1] *= size/d;
            axesNormals[axesNormals.length-1] *= size/d;
        }
        axesNormals.push(-(axesVertices[il+3]-axesVertices[il]));
        axesNormals.push(-(axesVertices[il+4]-axesVertices[il+1]));
        axesNormals.push(-(axesVertices[il+5]-axesVertices[il+2]));
        if(d>1e-8){
            axesNormals[axesNormals.length-1-2]   *= size/d;
            axesNormals[axesNormals.length-1-1] *= size/d;
            axesNormals[axesNormals.length-1] *= size/d;
        }
        axesIndexs_new.push(axesIdx_new++);
        axesIndexs_new.push(axesIdx_new++);
        axesIndexs_new.push(axesIdx_new++);
        axesIndexs_new.push(axesIdx_new++);
        axesIndexs_new.push(axesIdx_new++);
        axesIndexs_new.push(axesIdx_new++);
    }

    var ret = {};
    ret["vertices"] = axesVertices_new;
    ret["indices"] = axesIndexs_new;
    ret["normals"] = axesNormals;
    ret["colours"] = axesColours_new;
    return ret;

}

MGWebGL.prototype.drawAxes = function(invMat) {
    var axesOffset = vec3.create([20,20,0]);
    mat4.multiplyVec3(invMat,axesOffset);
    var right = vec3.create([1.0,0.0,0.0]);
    var up    = vec3.create([0.0,1.0,0.0]);
    mat4.multiplyVec3(invMat,up);
    mat4.multiplyVec3(invMat,right);
    var xoff = -this.origin[0]+this.zoom*axesOffset[0];
    var yoff = -this.origin[1]+this.zoom*axesOffset[1];
    var zoff = -this.origin[2]+this.zoom*axesOffset[2];
    this.gl.useProgram(this.shaderProgramThickLines);
    this.setMatrixUniforms(this.shaderProgramThickLines);
    this.gl.uniformMatrix4fv(this.shaderProgramThickLines.pMatrixUniform, false, this.pmvMatrix);
    this.gl.uniform3fv(this.shaderProgramThickLines.screenZ,this.screenZ);
    this.gl.uniform1f(this.shaderProgramThickLines.pixelZoom, 0.04*this.zoom);
    if(typeof(this.axesPositionBuffer) === "undefined"){
        this.axesPositionBuffer = this.gl.createBuffer();
        this.axesColourBuffer = this.gl.createBuffer();
        this.axesIndexBuffer = this.gl.createBuffer();
        this.axesNormalBuffer = this.gl.createBuffer();
        this.axesTextNormalBuffer = this.gl.createBuffer();
        this.axesTextColourBuffer = this.gl.createBuffer();
        this.axesTextPositionBuffer = this.gl.createBuffer();
        this.axesTextTexCoordBuffer = this.gl.createBuffer();
        this.axesTextIndexesBuffer = this.gl.createBuffer();
    }
    var x1 = 0.0;
    var y1 = 0.0;
    var z1 = 0.0;
    var x2 = this.zoom*3.0;
    var y2 = 0.0;
    var z2 = 0.0;
    var axesVertices = [];
    var axesColours = [];
    var axesIndexs = [];
    var axesIdx = 0;

    axesVertices.push(x1+xoff);
    axesVertices.push(y1+yoff);
    axesVertices.push(z1+zoff);
    axesVertices.push(x2+xoff);
    axesVertices.push(y2+yoff);
    axesVertices.push(z2+zoff);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesIndexs.push(axesIdx++);
    axesIndexs.push(axesIdx++);

    x2 = 0.0;
    y2 = this.zoom*3.0;
    z2 = 0.0;
    axesVertices.push(x1+xoff);
    axesVertices.push(y1+yoff);
    axesVertices.push(z1+zoff);
    axesVertices.push(x2+xoff);
    axesVertices.push(y2+yoff);
    axesVertices.push(z2+zoff);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesIndexs.push(axesIdx++);
    axesIndexs.push(axesIdx++);

    x2 = 0.0;
    y2 = 0.0;
    z2 = this.zoom*3.0;
    axesVertices.push(x1+xoff);
    axesVertices.push(y1+yoff);
    axesVertices.push(z1+zoff);
    axesVertices.push(x2+xoff);
    axesVertices.push(y2+yoff);
    axesVertices.push(z2+zoff);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(1.0);
    axesIndexs.push(axesIdx++);
    axesIndexs.push(axesIdx++);

    x1 = 2.0*this.zoom;
    y1 = 0.5*this.zoom;
    z1 = 0.0;
    x2 = this.zoom*3.0;
    y2 = 0.0;
    z2 = 0.0;
    axesVertices.push(x1+xoff);
    axesVertices.push(y1+yoff);
    axesVertices.push(z1+zoff);
    axesVertices.push(x2+xoff);
    axesVertices.push(y2+yoff);
    axesVertices.push(z2+zoff);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesIndexs.push(axesIdx++);
    axesIndexs.push(axesIdx++);

    x1 = 2.0*this.zoom;
    y1 = -0.5*this.zoom;
    z1 = 0.0;
    x2 = this.zoom*3.0;
    y2 = 0.0;
    z2 = 0.0;
    axesVertices.push(x1+xoff);
    axesVertices.push(y1+yoff);
    axesVertices.push(z1+zoff);
    axesVertices.push(x2+xoff);
    axesVertices.push(y2+yoff);
    axesVertices.push(z2+zoff);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesIndexs.push(axesIdx++);
    axesIndexs.push(axesIdx++);

    z1 = 0.5*this.zoom;
    y1 = 2.0*this.zoom;
    x1 = 0.0;
    x2 = 0.0;
    y2 = this.zoom*3.0;
    z2 = 0.0;
    axesVertices.push(x1+xoff);
    axesVertices.push(y1+yoff);
    axesVertices.push(z1+zoff);
    axesVertices.push(x2+xoff);
    axesVertices.push(y2+yoff);
    axesVertices.push(z2+zoff);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesIndexs.push(axesIdx++);
    axesIndexs.push(axesIdx++);

    z1 = -0.5*this.zoom;
    y1 = 2.0*this.zoom;
    x1 = 0.0;
    x2 = 0.0;
    y2 = this.zoom*3.0;
    z2 = 0.0;
    axesVertices.push(x1+xoff);
    axesVertices.push(y1+yoff);
    axesVertices.push(z1+zoff);
    axesVertices.push(x2+xoff);
    axesVertices.push(y2+yoff);
    axesVertices.push(z2+zoff);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesIndexs.push(axesIdx++);
    axesIndexs.push(axesIdx++);

    x1 = 0.5*this.zoom;
    z1 = 2.0*this.zoom;
    y1 = 0.0;
    x2 = 0.0;
    y2 = 0.0;
    z2 = this.zoom*3.0;
    axesVertices.push(x1+xoff);
    axesVertices.push(y1+yoff);
    axesVertices.push(z1+zoff);
    axesVertices.push(x2+xoff);
    axesVertices.push(y2+yoff);
    axesVertices.push(z2+zoff);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(1.0);
    axesIndexs.push(axesIdx++);
    axesIndexs.push(axesIdx++);

    x1 = -0.5*this.zoom;
    z1 = 2.0*this.zoom;
    y1 = 0.0;
    x2 = 0.0;
    y2 = 0.0;
    z2 = this.zoom*3.0;
    axesVertices.push(x1+xoff);
    axesVertices.push(y1+yoff);
    axesVertices.push(z1+zoff);
    axesVertices.push(x2+xoff);
    axesVertices.push(y2+yoff);
    axesVertices.push(z2+zoff);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(1.0);
    axesColours.push(0.0);
    axesColours.push(0.0);
    axesColours.push(1.0);
    axesColours.push(1.0);
    axesIndexs.push(axesIdx++);
    axesIndexs.push(axesIdx++);

    var size = 1.0;
    var thickLines = this.linesToThickLines(axesVertices,axesColours,size);
    var axesNormals = thickLines["normals"];
    var axesVertices_new = thickLines["vertices"];
    var axesColours_new = thickLines["colours"];
    var axesIndexs_new = thickLines["indices"];

    this.gl.depthFunc(this.gl.ALWAYS);

    this.gl.enableVertexAttribArray(this.shaderProgramThickLines.vertexNormalAttribute);

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.axesNormalBuffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(axesNormals), this.gl.DYNAMIC_DRAW);
    this.gl.vertexAttribPointer(this.shaderProgramThickLines.vertexNormalAttribute, 3, this.gl.FLOAT, false, 0, 0);

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.axesPositionBuffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(axesVertices_new), this.gl.DYNAMIC_DRAW);
    this.gl.vertexAttribPointer(this.shaderProgramThickLines.vertexPositionAttribute, 3, this.gl.FLOAT, false, 0, 0);

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.axesColourBuffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(axesColours_new), this.gl.DYNAMIC_DRAW);
    this.gl.vertexAttribPointer(this.shaderProgramThickLines.vertexColourAttribute, 4, this.gl.FLOAT, false, 0, 0);

    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.axesIndexBuffer);
    if(this.ext){
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(axesIndexs_new), this.gl.DYNAMIC_DRAW);
        this.gl.drawElements(this.gl.TRIANGLES, axesIndexs_new.length, this.gl.UNSIGNED_INT, 0);
    }else{
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(axesIndexs_new), this.gl.DYNAMIC_DRAW);
        this.gl.drawElements(this.gl.TRIANGLES,axesIndexs_new.length, this.gl.UNSIGNED_SHORT, 0);
    }

    this.gl.useProgram(this.shaderProgramTextBackground);
    this.gl.enableVertexAttribArray(this.shaderProgramTextBackground.vertexTextureAttribute);
    this.setMatrixUniforms(this.shaderProgramTextBackground);

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.axesTextNormalBuffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array([0,0,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1]), this.gl.STATIC_DRAW);
    this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexNormalAttribute, 3, this.gl.FLOAT, false, 0, 0);

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.axesTextColourBuffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array([1,0,0,1, 1,0,0,1, 1,0,0,1, 1,1,0,1, 1,1,0,1, 1,1,0,1]), this.gl.STATIC_DRAW);
    this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexColourAttribute, 4, this.gl.FLOAT, false, 0, 0);

    var tSizeX = 2.0 * this.textCtx.canvas.width/this.textCtx.canvas.height * this.zoom;
    var tSizeY = 2.0 * this.zoom;


    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.axesTextTexCoordBuffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array([0,1, 1,1, 1,0, 0,1, 1,0, 0,0]), this.gl.STATIC_DRAW);
    this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexTextureAttribute, 2, this.gl.FLOAT, false, 0, 0);

    var textColour = "black";
    var y = this.background_colour[0]*0.299 + this.background_colour[1]*0.587 + this.background_colour[2]*0.114;
    if(y<0.5){
        textColour = "white";
    }

    // Draw an x
    this.makeTextCanvas("x", 512, 32, textColour);
    var data = this.textCtx.getImageData(0,0,512, 32);
    this.gl.texImage2D(this.gl.TEXTURE_2D, 0, this.gl.RGBA, this.gl.RGBA, this.gl.UNSIGNED_BYTE, data);
    this.gl.texSubImage2D(this.gl.TEXTURE_2D, 0, 0, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, data);

    var textPositions = [];
    var base_x = xoff+3.0*this.zoom;
    var base_y = yoff;
    var base_z = zoff;
    textPositions = textPositions.concat([base_x,base_y,base_z, base_x+tSizeX*right[0],base_y+tSizeX*right[1],base_z+tSizeX*right[2], base_x+tSizeY*up[0]+tSizeX*right[0],base_y+tSizeY*up[1]+tSizeX*right[1],base_z+tSizeY*up[2]+tSizeX*right[2]]);
    textPositions = textPositions.concat([base_x,base_y,base_z, base_x+tSizeY*up[0]+tSizeX*right[0],base_y+tSizeY*up[1]+tSizeX*right[1],base_z+tSizeY*up[2]+tSizeX*right[2], base_x+tSizeY*up[0],base_y+tSizeY*up[1],base_z+tSizeY*up[2]]);
    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.axesTextPositionBuffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(textPositions), this.gl.DYNAMIC_DRAW);
    this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexPositionAttribute, 3, this.gl.FLOAT, false, 0, 0);

    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER,this.axesTextIndexesBuffer);
    if(this.ext){
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array([0,1,2,3,4,5]), this.gl.STATIC_DRAW);
        this.gl.drawElements(this.gl.TRIANGLES, 6, this.gl.UNSIGNED_INT, 0);
    }else{
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array([0,1,2,3,4,5]), this.gl.STATIC_DRAW);
        this.gl.drawElements(this.gl.TRIANGLES, 6, this.gl.UNSIGNED_SHORT, 0);
    }

    // Draw an y
    this.makeTextCanvas("y", 512, 32, textColour);
    var data = this.textCtx.getImageData(0,0,512, 32);
    this.gl.texImage2D(this.gl.TEXTURE_2D, 0, this.gl.RGBA, this.gl.RGBA, this.gl.UNSIGNED_BYTE, data);
    this.gl.texSubImage2D(this.gl.TEXTURE_2D, 0, 0, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, data);

    textPositions = [];
    base_x = xoff;
    base_y = yoff+3.0*this.zoom;
    base_z = zoff;
    textPositions = textPositions.concat([base_x,base_y,base_z, base_x+tSizeX*right[0],base_y+tSizeX*right[1],base_z+tSizeX*right[2], base_x+tSizeY*up[0]+tSizeX*right[0],base_y+tSizeY*up[1]+tSizeX*right[1],base_z+tSizeY*up[2]+tSizeX*right[2]]);
    textPositions = textPositions.concat([base_x,base_y,base_z, base_x+tSizeY*up[0]+tSizeX*right[0],base_y+tSizeY*up[1]+tSizeX*right[1],base_z+tSizeY*up[2]+tSizeX*right[2], base_x+tSizeY*up[0],base_y+tSizeY*up[1],base_z+tSizeY*up[2]]);
    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.axesTextPositionBuffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(textPositions), this.gl.DYNAMIC_DRAW);
    this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexPositionAttribute, 3, this.gl.FLOAT, false, 0, 0);

    if(this.ext){
        this.gl.drawElements(this.gl.TRIANGLES, 6, this.gl.UNSIGNED_INT, 0);
    }else{
        this.gl.drawElements(this.gl.TRIANGLES, 6, this.gl.UNSIGNED_SHORT, 0);
    }

    // Draw an z
    this.makeTextCanvas("z", 512, 32, textColour);
    var data = this.textCtx.getImageData(0,0,512, 32);
    this.gl.texImage2D(this.gl.TEXTURE_2D, 0, this.gl.RGBA, this.gl.RGBA, this.gl.UNSIGNED_BYTE, data);
    this.gl.texSubImage2D(this.gl.TEXTURE_2D, 0, 0, 0, this.gl.RGBA, this.gl.UNSIGNED_BYTE, data);

    textPositions = [];
    base_x = xoff;
    base_y = yoff;
    base_z = zoff+3.0*this.zoom;
    textPositions = textPositions.concat([base_x,base_y,base_z, base_x+tSizeX*right[0],base_y+tSizeX*right[1],base_z+tSizeX*right[2], base_x+tSizeY*up[0]+tSizeX*right[0],base_y+tSizeY*up[1]+tSizeX*right[1],base_z+tSizeY*up[2]+tSizeX*right[2]]);
    textPositions = textPositions.concat([base_x,base_y,base_z, base_x+tSizeY*up[0]+tSizeX*right[0],base_y+tSizeY*up[1]+tSizeX*right[1],base_z+tSizeY*up[2]+tSizeX*right[2], base_x+tSizeY*up[0],base_y+tSizeY*up[1],base_z+tSizeY*up[2]]);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(textPositions), this.gl.DYNAMIC_DRAW);
    this.gl.vertexAttribPointer(this.shaderProgramTextBackground.vertexPositionAttribute, 3, this.gl.FLOAT, false, 0, 0);

    if(this.ext){
        this.gl.drawElements(this.gl.TRIANGLES, 6, this.gl.UNSIGNED_INT, 0);
    }else{
        this.gl.drawElements(this.gl.TRIANGLES, 6, this.gl.UNSIGNED_SHORT, 0);
    }

    this.gl.disableVertexAttribArray(this.shaderProgramTextBackground.vertexTextureAttribute);
    this.gl.depthFunc(this.gl.LESS)
    
}

MGWebGL.prototype.doMouseUp = function(event,self) {
    self.mouseDown = false;
}

MGWebGL.prototype.drawSceneDirty = function() {
    this.doRedraw = true;
}

MGWebGL.prototype.drawSceneIfDirty = function() {
    if(this.doRedraw){
       this.doRedraw = false;
       this.drawScene();
    }
}

MGWebGL.prototype.reContourMaps = function() {
    var self = this;
    for(ilm=0;ilm<self.liveUpdatingMaps.length;ilm++){
        var contourLevel = self.liveUpdatingMaps[ilm].contourLevel;
        var mapColour = self.liveUpdatingMaps[ilm].mapColour;
        //console.log(contourLevel);
        var fracMat = self.liveUpdatingMaps[ilm].fracMat;
        var nCells_x = self.liveUpdatingMaps[ilm].nCells_x;
        var nCells_y = self.liveUpdatingMaps[ilm].nCells_y;
        var nCells_z = self.liveUpdatingMaps[ilm].nCells_z;
        var cellLength_x = self.liveUpdatingMaps[ilm].cellLength_x;
        var cellLength_y = self.liveUpdatingMaps[ilm].cellLength_y;
        var cellLength_z = self.liveUpdatingMaps[ilm].cellLength_z;
        var cell = self.liveUpdatingMaps[ilm].cell;
        //var xfrac = -(fracMat[0]*self.origin[0] + fracMat[1]*self.origin[1]+ fracMat[2]*self.origin[2]);
        //var yfrac = -(fracMat[3]*self.origin[0] + fracMat[4]*self.origin[1]+ fracMat[5]*self.origin[2]);
        //var zfrac = -(fracMat[6]*self.origin[0] + fracMat[7]*self.origin[1]+ fracMat[8]*self.origin[2]);

        var xfrac = -(cell.matrix_frac[0]*self.origin[0] + cell.matrix_frac[4]*self.origin[1]+ cell.matrix_frac[8]*self.origin[2]);
        var yfrac = -(cell.matrix_frac[1]*self.origin[0] + cell.matrix_frac[5]*self.origin[1]+ cell.matrix_frac[9]*self.origin[2]);
        var zfrac = -(cell.matrix_frac[2]*self.origin[0] + cell.matrix_frac[6]*self.origin[1]+ cell.matrix_frac[10]*self.origin[2]);
        //console.log(xfrac+" "+xfrac2);
        //console.log(yfrac+" "+yfrac2);
        //console.log(zfrac+" "+zfrac2);
        
        //console.log(self.origin[0]+" "+self.origin[1]+" "+self.origin[2]);
        var xCentre = (xfrac-Math.floor(xfrac))*nCells_x;
        var yCentre = (yfrac-Math.floor(yfrac))*nCells_y;
        var zCentre = (zfrac-Math.floor(zfrac))*nCells_z;
        var xshift_offset = xCentre - (xfrac-Math.floor(xfrac))*nCells_x;
        var yshift_offset = yCentre - (yfrac-Math.floor(yfrac))*nCells_y;
        var zshift_offset = zCentre - (zfrac-Math.floor(zfrac))*nCells_z;
        // A 10 angstrom parallelpiped
        var radCells_x = 10./cellLength_x;
        var radCells_y = 10./cellLength_y;
        var radCells_z = 10./cellLength_z;
        //var xShiftFracOffset = xCentre*cellLength_x-(xfrac-Math.floor(xfrac))*nCells_x*cellLength_x;
        //var yShiftFracOffset = yCentre*cellLength_y-(yfrac-Math.floor(yfrac))*nCells_y*cellLength_y;
        //var zShiftFracOffset = zCentre*cellLength_z-(zfrac-Math.floor(zfrac))*nCells_z*cellLength_z;
        //console.log(radCells_x+" "+radCells_y+" "+radCells_z);
        var minX = parseInt(Math.floor(xCentre - radCells_x));
        var maxX = parseInt(Math.floor(xCentre + radCells_x));
        var minY = parseInt(Math.floor(yCentre - radCells_y));
        var maxY = parseInt(Math.floor(yCentre + radCells_y));
        var minZ = parseInt(Math.floor(zCentre - radCells_z));
        var maxZ = parseInt(Math.floor(zCentre + radCells_z));
        var minXorig = minX;
        var maxXorig = maxX;
        var minYorig = minY;
        var maxYorig = maxY;
        var minZorig = minZ;
        var maxZorig = maxZ;
        if(minX<0) minX = 0;
        if(minY<0) minY = 0;
        if(minZ<0) minZ = 0;
        if(maxX>=nCells_x) maxX = nCells_x;
        if(maxY>=nCells_y) maxY = nCells_y;
        if(maxZ>=nCells_z) maxZ = nCells_z;

        var doChickenWire = true;

        if(minX<maxX && minY<maxY && minZ < maxZ){
            var map = new CIsoSurface();
            map.cell = cell;
            var start = new Date().getTime();
            map.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,minX,minY,minZ,maxX,maxY,maxZ,doChickenWire);
            var end = new Date().getTime();
            var time = end - start;
            //console.log('generate map on the fly: ' + time*0.001+"s");
            var start = new Date().getTime();
            var xShift = xfrac * cellLength_x * nCells_x - 0.5*(cellLength_x * parseFloat(maxX-minX+1)) ;
            var yShift = yfrac * cellLength_y * nCells_y - 0.5*(cellLength_y * parseFloat(maxY-minY+1)) ;
            var zShift = zfrac * cellLength_z * nCells_z - 0.5*(cellLength_z * parseFloat(maxZ-minZ+1)) ;
            var xShift2 = xfrac * cellLength_x * nCells_x - 0.5*(cellLength_x * parseFloat(maxX-minX+1)) ;
            var yShift2 = yfrac * cellLength_y * nCells_y - 0.5*(cellLength_y * parseFloat(maxY-minY+1)) ;
            var zShift2 = zfrac * cellLength_z * nCells_z - 0.5*(cellLength_z * parseFloat(maxZ-minZ+1)) ;

            // I think this correct.
            var xorth = (cell.matrix_orth[0]*xfrac + cell.matrix_orth[4]*yfrac + cell.matrix_orth[8]*zfrac);
            var yorth = (cell.matrix_orth[1]*xfrac + cell.matrix_orth[5]*yfrac + cell.matrix_orth[9]*zfrac);
            var zorth = (cell.matrix_orth[2]*xfrac + cell.matrix_orth[6]*yfrac + cell.matrix_orth[10]*zfrac);

            var xorthhalf = -(cell.matrix_orth[0]*0.5*(maxX-minX+1.0)/nCells_x + cell.matrix_orth[4]*0.5*(maxY-minY+1.0)/nCells_y + cell.matrix_orth[8]*0.5*(maxZ-minZ+1.0)/nCells_z);
            var yorthhalf = -(cell.matrix_orth[1]*0.5*(maxX-minX+1.0)/nCells_x + cell.matrix_orth[5]*0.5*(maxY-minY+1.0)/nCells_y + cell.matrix_orth[9]*0.5*(maxZ-minZ+1.0)/nCells_z);
            var zorthhalf = -(cell.matrix_orth[2]*0.5*(maxX-minX+1.0)/nCells_x + cell.matrix_orth[6]*0.5*(maxY-minY+1.0)/nCells_y + cell.matrix_orth[10]*0.5*(maxZ-minZ+1.0)/nCells_z);
            //xorthhalf *= (maxX-minX+1.0)/nCells_x;
            //yorthhalf *= (maxY-minY+1.0)/nCells_y;
            //zorthhalf *= (maxZ-minZ+1.0)/nCells_z;

            xorth += xorthhalf;
            yorth += yorthhalf;
            zorth += zorthhalf;

            xShift  = xorth;
            xShift2 = xorth;
            yShift  = yorth;
            yShift2 = yorth;
            zShift  = zorth;
            zShift2 = zorth;

            //console.log("xshift: "+xShift+" "+xorth)
            //console.log("yshift: "+yShift+" "+yorth)
            //console.log("zshift: "+zShift+" "+zorth)

            var vertices2 = [];
            var normals2 = [];
            var indices2 = [];
            var colours2 = [];
            var vertices3 = [];
            var normals3 = [];
            var indices3 = [];
            var colours3 = [];
            var vertices4 = [];
            var normals4 = [];
            var indices4 = [];
            var colours4 = [];
            var vertices5 = [];
            var normals5 = [];
            var indices5 = [];
            var colours5 = [];
            var vertices6 = [];
            var normals6 = [];
            var indices6 = [];
            var colours6 = [];
            var vertices7 = [];
            var normals7 = [];
            var indices7 = [];
            var colours7 = [];
            var vertices8 = [];
            var normals8 = [];
            var indices8 = [];
            var colours8 = [];
            //console.log("orig shifts: "+xShift+" "+yShift+" "+zShift);
            var xShift_new  = xShift;
            var xShift2_new = xShift2;
            var yShift_new  = yShift;
            var yShift2_new = yShift2;
            var zShift_new  = zShift;
            var zShift2_new = zShift2;
            if(maxXorig>=nCells_x){
                xShift += 0.5 * (maxX-maxXorig) * cellLength_x;
                xShift2 -= 0.5 * ((maxXorig)-maxX+2*minX-2*nCells_x) * cellLength_x;
            }
            if(minXorig<0){
                xShift += 0.5 * (minX-minXorig) * cellLength_x;
                xShift2 += 0.5 * (minXorig) * cellLength_x;

            }
            if(maxYorig>=nCells_y){
                yShift += 0.5 * (maxY-maxYorig) * cellLength_y;
                yShift2 -= 0.5 * ((maxYorig)-maxY+2*minY-2*nCells_y) * cellLength_y;

            }
            if(minYorig<0){
                yShift += 0.5 * (minY-minYorig) * cellLength_y;
                yShift2 += 0.5 * (minYorig) * cellLength_y;

            }
            if(maxZorig>=nCells_z){
                zShift += 0.5 * (maxZ-maxZorig) * cellLength_z;
                zShift2 -= 0.5 * ((maxZorig)-maxZ+2*minZ-2*nCells_z) * cellLength_z;

            }
            if(minZorig<0){
                zShift += 0.5 * (minZ-minZorig) * cellLength_z;
                zShift2 += 0.5 * (minZorig) * cellLength_z;

            }

            // I think this is "the right way"!
            var xdiff1 = 0.0;
            var xdiff2 = 0.0;
            var ydiff1 = 0.0;
            var ydiff2 = 0.0;
            var zdiff1 = 0.0;
            var zdiff2 = 0.0;
            var xdiff3 = 0.0;
            var xdiff4 = 0.0;
            if(maxXorig>=nCells_x){
                xdiff1 = (maxX-maxXorig)/nCells_x;
                xdiff3 = (maxX-maxXorig)/nCells_x;
                xdiff2 = (((-maxXorig)+maxX-2*minX+2*nCells_x)/nCells_x);
                xdiff4 = (((-maxXorig)+maxX-2*minX+2*nCells_x)/nCells_x);
            }
            if(minXorig<0){
                xdiff1 = (minX-minXorig)/nCells_x;
                xdiff3 = (minX-minXorig)/nCells_x;
                xdiff2 = minXorig/nCells_x;
                xdiff4 = minXorig/nCells_x;
            }
            if(maxYorig>=nCells_y){
                ydiff1 = (maxY-maxYorig)/nCells_y;
                ydiff2 = (((-maxYorig)+maxY-2*minY+2*nCells_y)/nCells_y);
                xdiff2 -= 2.0* (maxY-maxYorig+2.0*radCells_y)/nCells_y * cell.matrix_orth[4] / (nCells_y * cellLength_y);
                xdiff3 += 2.0* (maxY-maxYorig+2.0*radCells_y)/nCells_y * cell.matrix_orth[4] / (nCells_y * cellLength_y);
            }
            if(minYorig<0){
                ydiff1 = (minY-minYorig)/nCells_y;
                ydiff2 = minYorig/nCells_y;
                // Are we to assume x along x, z along z, and y any direction, so that only dx matters?
                xdiff2 -= 2.0* minYorig/nCells_y * cell.matrix_orth[4] / (nCells_y * cellLength_y);
                xdiff3 += 2.0* minYorig/nCells_y * cell.matrix_orth[4] / (nCells_y * cellLength_y);
            }
            if(maxZorig>=nCells_z){
                zdiff1 = (maxZ-maxZorig)/nCells_z;
                zdiff2 = (((-maxZorig)+maxZ-2*minZ+2*nCells_z)/nCells_z);
            }
            if(minZorig<0){
                zdiff1 = (minZ-minZorig)/nCells_z;
                zdiff2 = minZorig/nCells_z;
            }

            var xShift3 = xShift_new + 0.5*(cell.matrix_orth[0]*xdiff3 + cell.matrix_orth[4]*ydiff1 + cell.matrix_orth[8]*zdiff1);
            var xShift4 = xShift2_new + 0.5*(cell.matrix_orth[0]*xdiff4 + cell.matrix_orth[4]*ydiff2 + cell.matrix_orth[8]*zdiff2);

            xShift_new  += 0.5*(cell.matrix_orth[0]*xdiff1 + cell.matrix_orth[4]*ydiff1 + cell.matrix_orth[8]*zdiff1);
            xShift2_new += 0.5*(cell.matrix_orth[0]*xdiff2 + cell.matrix_orth[4]*ydiff2 + cell.matrix_orth[8]*zdiff2);
            yShift_new  += 0.5*(cell.matrix_orth[1]*xdiff1 + cell.matrix_orth[5]*ydiff1 + cell.matrix_orth[9]*zdiff1);
            yShift2_new += 0.5*(cell.matrix_orth[1]*xdiff2 + cell.matrix_orth[5]*ydiff2 + cell.matrix_orth[9]*zdiff2);
            zShift_new  += 0.5*(cell.matrix_orth[2]*xdiff1 + cell.matrix_orth[6]*ydiff1 + cell.matrix_orth[10]*zdiff1);
            zShift2_new += 0.5*(cell.matrix_orth[2]*xdiff2 + cell.matrix_orth[6]*ydiff2 + cell.matrix_orth[10]*zdiff2);


            //console.log(self.origin);
            //console.log(minXorig+" "+maxXorig+", "+minYorig+" "+maxYorig+", "+minZorig+" "+maxZorig+" ,(maxY: "+maxY+")");
            //console.log("new shifts: "+xShift+" "+yShift+" "+zShift);
            //console.log("new shifts2: "+xShift2+" "+yShift2+" "+zShift2);
            //console.log("newnew shifts: "+xShift_new+" "+yShift_new+" "+zShift_new);
            //console.log("newnew shifts2: "+xShift2_new+" "+yShift2_new+" "+zShift2_new);
            //console.log(xfrac+" "+yfrac+" "+zfrac);
            xShift = xShift_new;
            xShift2 = xShift2_new;
            yShift = yShift_new;
            yShift2 = yShift2_new;
            zShift = zShift_new;
            zShift2 = zShift2_new;


            if(maxXorig>=nCells_x){
                var map2 = new CIsoSurface();
                map2.cell = cell;
                map2.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,minY,minZ,maxXorig-nCells_x,maxY,maxZ,doChickenWire);
                vertices2 = map2.returnVertices(xShift2,yShift,zShift,doChickenWire);
                normals2 = map2.returnNormals_new();
                var indices2 = map2.returnIndices();
                for(var i=0;i<vertices2.length/3;i++){
                    colours2.push(mapColour[0]);
                    colours2.push(mapColour[1]);
                    colours2.push(mapColour[2]);
                    colours2.push(mapColour[3]);
                }
            }
            if(minXorig<0){
                var map2 = new CIsoSurface();
                map2.cell = cell;
                map2.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x+minXorig,minY,minZ,nCells_x,maxY,maxZ,doChickenWire);
                vertices2 = map2.returnVertices(xShift2,yShift,zShift,doChickenWire);
                normals2 = map2.returnNormals_new();
                var indices2 = map2.returnIndices();
                for(var i=0;i<vertices2.length/3;i++){
                    colours2.push(mapColour[0]);
                    colours2.push(mapColour[1]);
                    colours2.push(mapColour[2]);
                    colours2.push(mapColour[3]);
                }
            }
            if(maxYorig>=nCells_y){
                var map3 = new CIsoSurface();
                map3.cell = cell;
                map3.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,minX,0,minZ,maxX,maxYorig-nCells_y,maxZ,doChickenWire);
                vertices3 = map3.returnVertices(xShift3,yShift2,zShift,doChickenWire);
                normals3 = map3.returnNormals_new();
                var indices3 = map3.returnIndices();
                for(var i=0;i<vertices3.length/3;i++){
                    colours3.push(mapColour[0]);
                    colours3.push(mapColour[1]);
                    colours3.push(mapColour[2]);
                    colours3.push(mapColour[3]);
                }
            }
            if(minYorig<0){
                var map3 = new CIsoSurface();
                map3.cell = cell;
                map3.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,minX,nCells_y+minYorig,minZ,maxX,nCells_y,maxZ,doChickenWire);
                vertices3 = map3.returnVertices(xShift3,yShift2,zShift,doChickenWire);
                normals3 = map3.returnNormals_new();
                var indices3 = map3.returnIndices();
                for(var i=0;i<vertices3.length/3;i++){
                    colours3.push(mapColour[0]);
                    colours3.push(mapColour[1]);
                    colours3.push(mapColour[2]);
                    colours3.push(mapColour[3]);
                }
            }
            if(maxZorig>=nCells_z){
                var map4 = new CIsoSurface();
                map4.cell = cell;
                map4.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,minX,minY,0,maxX,maxY,maxZorig-nCells_z,doChickenWire);
                vertices4 = map4.returnVertices(xShift,yShift,zShift2,doChickenWire);
                normals4 = map4.returnNormals_new();
                var indices4 = map4.returnIndices();
                for(var i=0;i<vertices4.length/3;i++){
                    colours4.push(mapColour[0]);
                    colours4.push(mapColour[1]);
                    colours4.push(mapColour[2]);
                    colours4.push(mapColour[3]);
                }
            }
            if(minZorig<0){
                var map4 = new CIsoSurface();
                map4.cell = cell;
                map4.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,minX,minY,nCells_z+minZorig,maxX,maxY,nCells_z,doChickenWire);
                vertices4 = map4.returnVertices(xShift,yShift,zShift2,doChickenWire);
                normals4 = map4.returnNormals_new();
                var indices4 = map4.returnIndices();
                for(var i=0;i<vertices4.length/3;i++){
                    colours4.push(mapColour[0]);
                    colours4.push(mapColour[1]);
                    colours4.push(mapColour[2]);
                    colours4.push(mapColour[3]);
                }
            }
            var map5 = new CIsoSurface();
            map5.cell = cell;
            if(maxXorig>=nCells_x&&maxYorig>=nCells_y){ //++
                map5.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,0,minZ,maxXorig-nCells_x,maxYorig-nCells_y,maxZ,doChickenWire);
            }
            if(minXorig<0&&minYorig<0){ //--
                map5.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x+minXorig,nCells_y+minYorig,minZ,nCells_x,nCells_y,maxZ,doChickenWire);
            }
            if(maxXorig>=nCells_x&&minYorig<0){ //+-
                map5.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,nCells_y+minYorig,minZ,maxXorig-nCells_x,nCells_y,maxZ,doChickenWire);
            }
            if(minXorig<0&&maxYorig>=nCells_y){ //-+
                map5.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x+minXorig,0,minZ,nCells_x,maxYorig-nCells_y,maxZ,doChickenWire);
            }
            //Map 5 xShift2, yShift2, zShift
            vertices5 = map5.returnVertices(xShift4,yShift2,zShift,doChickenWire);
            normals5 = map5.returnNormals_new();
            var indices5 = map5.returnIndices();
            for(var i=0;i<vertices5.length/3;i++){
                colours5.push(mapColour[0]);
                colours5.push(mapColour[1]);
                colours5.push(mapColour[2]);
                colours5.push(mapColour[3]);
            }

            var map6 = new CIsoSurface();
            map6.cell = cell;
            if(maxXorig>=nCells_x&&maxZorig>=nCells_z){ //++
                map6.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,minY,0,maxXorig-nCells_x,maxY,maxZorig-nCells_z,doChickenWire);
            }
            if(minXorig<0&&minZorig<0){ //--
                map6.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x+minXorig,minY,nCells_z+minZorig,nCells_x,maxY,nCells_z,doChickenWire);
            }
            if(maxXorig>=nCells_x&&minZorig<0){ //+-
                map6.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,minY,nCells_z+minZorig,maxXorig-nCells_x,maxY,nCells_z,doChickenWire);
            }
            if(minXorig<0&&maxZorig>=nCells_z){ //-+
                map6.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x+minXorig,minY,0,nCells_x,maxY,maxZorig-nCells_z,doChickenWire);
            }
            //Map 6 xShift2, yShift, zShift2
            vertices6 = map6.returnVertices(xShift2,yShift,zShift2,doChickenWire);
            normals6 = map6.returnNormals_new();
            var indices6 = map6.returnIndices();
            for(var i=0;i<vertices6.length/3;i++){
                colours6.push(mapColour[0]);
                colours6.push(mapColour[1]);
                colours6.push(mapColour[2]);
                colours6.push(mapColour[3]);
            }

            var map7 = new CIsoSurface();
            map7.cell = cell;
            if(maxYorig>=nCells_y&&maxZorig>=nCells_z){ //++
                map7.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,minX,0,0,maxX,maxYorig-nCells_y,maxZorig-nCells_z,doChickenWire);
            }
            if(minYorig<0&&minZorig<0){ //--
                map7.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,minX,nCells_y+minYorig,nCells_z+minZorig,maxX,nCells_y,nCells_z,doChickenWire);
            }
            if(maxYorig>=nCells_y&&minZorig<0){ //+-
                map7.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,minX,0,nCells_z+minZorig,maxX,maxYorig-nCells_y,nCells_z,doChickenWire);
            }
            if(minYorig<0&&maxZorig>=nCells_z){ //-+
                map7.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,minX,nCells_y+minYorig,0,maxX,nCells_y,maxZorig-nCells_z,doChickenWire);
            }
            //Map 7 xShift, yShift2, zShift2
            vertices7 = map7.returnVertices(xShift3,yShift2,zShift2,doChickenWire);
            normals7 = map7.returnNormals_new();
            var indices7 = map7.returnIndices();
            for(var i=0;i<vertices7.length/3;i++){
                colours7.push(mapColour[0]);
                colours7.push(mapColour[1]);
                colours7.push(mapColour[2]);
                colours7.push(mapColour[3]);
            }

            var map8 = new CIsoSurface();
            map8.cell = cell;
            if(maxXorig>=nCells_x&&maxYorig>=nCells_y&&maxZorig>=nCells_z){ //+++
                map8.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,0,0,maxXorig-nCells_x,maxYorig-nCells_y,maxZorig-nCells_z,doChickenWire);
            }
            if(minXorig<0&&minYorig<0&&minZorig<0){ //---
                map8.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x+minXorig,nCells_y+minYorig,nCells_z+minZorig,nCells_x,nCells_y,nCells_z,doChickenWire);
            }

            if(maxXorig>=nCells_x&&minYorig<0&&maxZorig>=nCells_z){ //+-+
                map8.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,nCells_y+minYorig,0,maxXorig-nCells_x,nCells_y,maxZorig-nCells_z,doChickenWire);
            }
            if(maxXorig>=nCells_x&&maxYorig>=nCells_y&&minZorig<0){ //++-
                map8.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,0,nCells_z+minZorig,maxXorig-nCells_x,maxYorig-nCells_y,nCells_z,doChickenWire);
            }
            if(minXorig<0&&maxYorig>=nCells_y&&maxZorig>=nCells_z){ //-++
                map8.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x+minXorig,0,0,nCells_x,maxYorig-nCells_y,maxZorig-nCells_z,doChickenWire);
            }

            if(maxXorig>=nCells_x&&minYorig<0&&minZorig<0){ //+--
                map8.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,0,nCells_y+minYorig,nCells_z+minZorig,maxXorig-nCells_x,nCells_y,nCells_z,doChickenWire);
            }
            if(minXorig<0&&minYorig<0&&maxZorig>=nCells_z){ //--+
                map8.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x+minXorig,nCells_y+minYorig,0,nCells_x,nCells_y,maxZorig-nCells_z,doChickenWire);
            }
            if(minXorig<0&&maxYorig>=nCells_y&&minZorig<0){ //-+-
                map8.GenerateSurfacePartial(self.liveUpdatingMaps[ilm].gridData[0],contourLevel,nCells_x,nCells_y,nCells_z,cellLength_x,cellLength_y,cellLength_z,nCells_x+minXorig,0,nCells_z+minZorig,nCells_x,maxYorig-nCells_y,nCells_z,doChickenWire);
            }
            //Map 8 xShift2, yShift2, zShift2
            vertices8 = map8.returnVertices(xShift4,yShift2,zShift2,doChickenWire);
            normals8 = map8.returnNormals_new();
            var indices8 = map8.returnIndices();
            for(var i=0;i<vertices8.length/3;i++){
                colours8.push(mapColour[0]);
                colours8.push(mapColour[1]);
                colours8.push(mapColour[2]);
                colours8.push(mapColour[3]);
            }

            var vertices = map.returnVertices(xShift_new,yShift_new,zShift_new,doChickenWire); // And is this where we add base (integral part of xfrac, etc.)?
            var normals = map.returnNormals_new();
            var indices = map.returnIndices();
            var end = new Date().getTime();
            var time = end - start;
            //console.log('copy buffers on the fly: ' + time*0.001+"s");
            var colours = [];
            for(var i=0;i<vertices.length/3;i++){
                colours.push(mapColour[0]);
                colours.push(mapColour[1]);
                colours.push(mapColour[2]);
                colours.push(mapColour[3]);
            }

            //FIXME - Maps do not need all the colour data. A fixed vertex attribute could be used.

            var doLines = false;
            if(doChickenWire){
                //FIXME - Normal data is not used for lines.

                var lineNormals = [1.0];
                var lineNormals2 = [1.0];
                var lineNormals3 = [1.0];
                var lineNormals4 = [1.0];
                var lineNormals5 = [1.0];
                var lineNormals6 = [1.0];
                var lineNormals7 = [1.0];
                var lineNormals8 = [1.0];
                var lineIndices = [];
                var lineIndices2 = [];
                var lineIndices3 = [];
                var lineIndices4 = [];
                var lineIndices5 = [];
                var lineIndices6 = [];
                var lineIndices7 = [];
                var lineIndices8 = [];

                var mapTriangleData = {"col_tri":[[colours],[colours2],[colours3],[colours4],[colours5],[colours6],[colours7],[colours8]], 
                    "norm_tri":[[lineNormals],[lineNormals2],[lineNormals3],[lineNormals4],[lineNormals5],[lineNormals6],[lineNormals7],[lineNormals8]],
                    "vert_tri":[[vertices],[vertices2],[vertices3],[vertices4],[vertices5],[vertices6],[vertices7],[vertices8]],
                    "idx_tri":[[lineIndices],[lineIndices2],[lineIndices3],[lineIndices4],[lineIndices5],[lineIndices6],[lineIndices7],[lineIndices8]],
                    "prim_types":[["LINES"],["LINES"],["LINES"],["LINES"],["LINES"],["LINES"],["LINES"],["LINES"]] };
            } else {

                var mapTriangleData = {"col_tri":[[colours],[colours2],[colours3],[colours4],[colours5],[colours6],[colours7],[colours8]], 
                    "norm_tri":[[normals],[normals2],[normals3],[normals4],[normals5],[normals6],[normals7],[normals8]],
                    "vert_tri":[[vertices],[vertices2],[vertices3],[vertices4],[vertices5],[vertices6],[vertices7],[vertices8]],
                    "idx_tri":[[indices],[indices2],[indices3],[indices4],[indices5],[indices6],[indices7],[indices8]],
                    "prim_types":[["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"],["TRIANGLES"]] };
            }

            self.updateBuffers(mapTriangleData,self.liveUpdatingMaps[ilm].theseBuffers);
        }
    }
}

MGWebGL.prototype.doMouseMove = function(event,self) {
    self.mouseMoved = true;

    if(true){
        var x;
        var y;
        var e = event;
        if (e.pageX || e.pageY) { 
            x = e.pageX;
            y = e.pageY;
        }
        else { 
            x = e.clientX ; 
            y = e.clientY ; 
        } 

        var offset = getOffsetRect(self.div);
        
        x -= offset.left;
        y -= offset.top;
    
        this.gl_cursorPos[0] = x;
        this.gl_cursorPos[1] = this.canvas.height-y;
        self.drawSceneDirty();
    }

    if(!self.mouseDown){
        self.init_x = event.pageX;
        self.init_y = event.pageY;
        return;
    }
    self.dx = event.pageX - self.init_x;
    self.dy = event.pageY - self.init_y;
    self.init_x = event.pageX;
    self.init_y = event.pageY;

    if((event.altKey&&event.shiftKey)||(self.mouseDownButton==2&&event.shiftKey)){
        var invQuat = quat4.create();
        quat4.inverse(self.myQuat,invQuat);
        var theMatrix = quat4.toMat4(invQuat);
        var xshift = vec3.create([self.dx/getDeviceScale(),0,0]);
        var yshift = vec3.create([0,self.dy/getDeviceScale(),0]);
        mat4.multiplyVec3(theMatrix,xshift);
        mat4.multiplyVec3(theMatrix,yshift);
        self.origin[0] += xshift[0]/8.*self.zoom;
        self.origin[1] += xshift[1]/8.*self.zoom;
        self.origin[2] += xshift[2]/8.*self.zoom;
        self.origin[0] -= yshift[0]/8.*self.zoom;
        self.origin[1] -= yshift[1]/8.*self.zoom;
        self.origin[2] -= yshift[2]/8.*self.zoom;
        self.drawSceneDirty();
        //console.log(self.origin);
        self.reContourMaps();
        return;
    }

    if(event.altKey){
        var factor = 1. - self.dy/50.
            self.zoom = self.zoom * factor
            if(self.zoom < .01){
                self.zoom = 0.01;
            }
        self.drawSceneDirty();
        return;
    }

    if(event.shiftKey){
        var zQ = createYQuatFromDY(0)
            if(event.pageX<self.gl.viewportWidth/2){
                zQ = createZQuatFromDX(-self.dy/5.)
            } else {
                zQ = createZQuatFromDX(self.dy/5.)
            }
        var yQ = createYQuatFromDY(0)
            if(event.pageY<self.gl.viewportHeight/2){
                yQ = createZQuatFromDX(self.dx/5.)
            } else {
                yQ = createZQuatFromDX(-self.dx/5.)
            }
        quat4.multiply(zQ,yQ)
            quat4.multiply(self.myQuat,zQ)
    } else {
        var xQ = createXQuatFromDX(self.dy);
        var yQ = createYQuatFromDY(-self.dx);
        quat4.multiply(xQ,yQ);
        quat4.multiply(self.myQuat,xQ);
    }

    var theMatrix = quat4.toMat4(self.myQuat);
    self.drawSceneDirty();
}

MGWebGL.prototype.doMouseDown = function(event,self) {
    self.init_x = event.pageX;
    self.init_y = event.pageY;
    self.mouseDown = true;
    self.mouseDownButton = event.button;
    self.mouseMoved = false;
}

MGWebGL.prototype.handleKeyDown = function(event,self) {
    if(event.keyCode==115||event.keyCode==83){
        /* Screenshot*/


        var oldOrigin = [this.origin[0],this.origin[1],this.origin[2]];

        // Getting up and right for doing tiling (in future?)
        var invQuat = quat4.create();
        quat4.inverse(this.myQuat,invQuat);

        var invMat =  quat4.toMat4(invQuat);

        var right = vec3.create([1.0,0.0,0.0]);
        var up    = vec3.create([0.0,1.0,0.0]);

        mat4.multiplyVec3(invMat,up);
        mat4.multiplyVec3(invMat,right);

        console.log(right);
        console.log(up);
        console.log(this.zoom);

        var mag = 1;

        var ncells_x = mag;
        var ncells_y = mag;

        var saveCanvas = document.createElement("canvas");
        saveCanvas.width = this.canvas.width*ncells_x;
        saveCanvas.height = this.canvas.height*ncells_y;
        var ctx = saveCanvas.getContext("2d");

        this.zoom /= ncells_x;

        var jj = 0;
        for(var j=Math.floor(-ncells_y/2);j<Math.floor(ncells_y/2);j++){
            var ii = 0;
            for(var i=Math.floor(-ncells_x/2);i<Math.floor(ncells_x/2);i++){
                var x_off = (2.0*i+1+ncells_x%2);
                var y_off = (2.0*j+1+ncells_y%2);

                this.origin = [oldOrigin[0], oldOrigin[1], oldOrigin[2]];
                this.origin[0] += this.zoom*right[0]*24.0*x_off + this.zoom*up[0]*24.0*y_off;
                this.origin[1] += this.zoom*right[1]*24.0*x_off + this.zoom*up[1]*24.0*y_off;
                this.origin[2] += this.zoom*right[2]*24.0*x_off + this.zoom*up[2]*24.0*y_off;


                self.save_pixel_data = true;
                this.drawScene();
                var pixels = self.pixel_data;

                var imgData = ctx.createImageData(this.canvas.width, this.canvas.height);
                var data = imgData.data;

                for (var pixi = 0; pixi < this.canvas.height; pixi++) {
                    for (var pixj = 0; pixj < this.canvas.width * 4; pixj++) {
                        data[(this.canvas.height-pixi-1)*this.canvas.width * 4 + pixj] = pixels[pixi*this.canvas.width * 4 + pixj];
                    }
                }
                ctx.putImageData(imgData, (ncells_x-ii-1)*this.canvas.width, jj*this.canvas.height);
                ii++;
            }
            jj++;
        }
        

        this.zoom *= ncells_x;
        this.origin = [oldOrigin[0], oldOrigin[1], oldOrigin[2]];
        self.save_pixel_data = false;
        this.drawScene();

        /* Now we need to scale down in steps .... */

        /*
        // Bah, useless. Still not as good as browser shrinking to fit (which is excellent generally).
        var can2 = document.createElement('canvas');
        can2.width = this.canvas.width*ncells_x*0.6933612743506347;
        can2.height = this.canvas.width*ncells_y*0.6933612743506347;
        var ctx2 = can2.getContext('2d');

        ctx2.drawImage(saveCanvas, 0, 0, this.canvas.width*ncells_x*0.6933612743506347, this.canvas.width*ncells_y*0.6933612743506347);
        ctx2.drawImage(can2, 0, 0, this.canvas.width*ncells_x*0.6933612743506347*0.6933612743506347, this.canvas.width*ncells_y*0.6933612743506347*0.6933612743506347);
        ctx2.drawImage(can2, 0, 0, this.canvas.width*ncells_x*0.6933612743506347*0.6933612743506347*0.6933612743506347, this.canvas.width*ncells_y*0.6933612743506347*0.6933612743506347*0.6933612743506347);
        var image = can2.toDataURL();
        */

        var image = saveCanvas.toDataURL();


        var newwindow = window.open();
        newwindow.document.write('<img src="'+image+'"/>');
        newwindow.document.close()
    }
    // FIXME, we need an active map, like Coot.
    if(event.keyCode==187){ //+ (actually equals, sigh)
        for(ilm=0;ilm<self.liveUpdatingMaps.length;ilm++){
            self.liveUpdatingMaps[ilm].contourLevel += 0.05;
        }
        self.drawSceneDirty();
        self.reContourMaps();
        self.drawScene();
    }
    if(event.keyCode==189){ //-
        for(ilm=0;ilm<self.liveUpdatingMaps.length;ilm++){
            self.liveUpdatingMaps[ilm].contourLevel -= 0.05;
        }
        self.drawSceneDirty();
        self.reContourMaps();
        self.drawScene();
    }
    if(event.keyCode==114||event.keyCode==82){
        self.myQuat = quat4.create([0,0,0,-1]);
        self.zoom = 1.0;
        self.clickedAtoms = [];
        self.drawScene();
    }
    if(event.keyCode==99||event.keyCode==67){
        self.clickedAtoms = [];
        self.drawScene();
    }
    if(event.keyCode==37){
        var invQuat = quat4.create();
        quat4.inverse(self.myQuat,invQuat);
        var theMatrix = quat4.toMat4(invQuat);
        var xshift = vec3.create([-4./getDeviceScale(),0,0]);
        mat4.multiplyVec3(theMatrix,xshift);
        self.origin[0] += xshift[0]/8.*self.zoom;
        self.origin[1] += xshift[1]/8.*self.zoom;
        self.origin[2] += xshift[2]/8.*self.zoom;
        self.drawSceneDirty();
        //console.log(self.origin);
        self.reContourMaps();
    }
    if(event.keyCode==39){
        var invQuat = quat4.create();
        quat4.inverse(self.myQuat,invQuat);
        var theMatrix = quat4.toMat4(invQuat);
        var xshift = vec3.create([4./getDeviceScale(),0,0]);
        mat4.multiplyVec3(theMatrix,xshift);
        self.origin[0] += xshift[0]/8.*self.zoom;
        self.origin[1] += xshift[1]/8.*self.zoom;
        self.origin[2] += xshift[2]/8.*self.zoom;
        self.drawSceneDirty();
        //console.log(self.origin);
        self.reContourMaps();
    }
    if(event.keyCode==38){
        var invQuat = quat4.create();
        quat4.inverse(self.myQuat,invQuat);
        var theMatrix = quat4.toMat4(invQuat);
        var yshift = vec3.create([0,4./getDeviceScale(),0]);
        mat4.multiplyVec3(theMatrix,yshift);
        self.origin[0] += yshift[0]/8.*self.zoom;
        self.origin[1] += yshift[1]/8.*self.zoom;
        self.origin[2] += yshift[2]/8.*self.zoom;
        self.drawSceneDirty();
        //console.log(self.origin);
        self.reContourMaps();
    }
    if(event.keyCode==40){
        var invQuat = quat4.create();
        quat4.inverse(self.myQuat,invQuat);
        var theMatrix = quat4.toMat4(invQuat);
        var yshift = vec3.create([0,-4./getDeviceScale(),0]);
        mat4.multiplyVec3(theMatrix,yshift);
        self.origin[0] += yshift[0]/8.*self.zoom;
        self.origin[1] += yshift[1]/8.*self.zoom;
        self.origin[2] += yshift[2]/8.*self.zoom;
        self.drawSceneDirty();
        //console.log(self.origin);
        self.reContourMaps();
    }
    
}

// Puts text in center of canvas.
MGWebGL.prototype.makeTextCanvas = function(text, width, height, textColour) {
    this.textCanvasInitialized = false;
    if(!this.textCanvasInitialized){
        this.textCtx.canvas.width  = width;
        this.textCtx.canvas.height = height;
        this.textCtx.font = "20px helvetica";
        this.textCtx.textAlign = "left";
        this.textCtx.textBaseline = "middle";
        this.textCanvasInitialized = true
    }
    this.textCtx.fillStyle = textColour;
    this.textCtx.clearRect(0, 0, this.textCtx.canvas.width, this.textCtx.canvas.height);
    this.textCtx.fillText(text,  0, height / 2);
}

MGWebGL.prototype.createVertexBuffer = function(tri){
    this.displayBuffers[this.currentBufferIdx].triangleVertexPositionBuffer.push(this.gl.createBuffer());
    this.displayBuffers[this.currentBufferIdx].triangleVertexPositionBuffer[this.displayBuffers[this.currentBufferIdx].triangleVertexPositionBuffer.length-1].numItems = 0;
    this.displayBuffers[this.currentBufferIdx].triangleVertices.push([]);
    for (var j=0; j<tri.length; j++){
        this.displayBuffers[this.currentBufferIdx].triangleVertices[this.displayBuffers[this.currentBufferIdx].triangleVertices.length-1].push(parseFloat(tri[j]));
        this.displayBuffers[this.currentBufferIdx].triangleVertexPositionBuffer[this.displayBuffers[this.currentBufferIdx].triangleVertexPositionBuffer.length-1].numItems++;
    }
    this.displayBuffers[this.currentBufferIdx].triangleVertexPositionBuffer[this.displayBuffers[this.currentBufferIdx].triangleVertexPositionBuffer.length-1].numItems /= 3;
}

MGWebGL.prototype.createNormalBuffer = function(norm){
    this.displayBuffers[this.currentBufferIdx].triangleNormals.push([]);
    this.displayBuffers[this.currentBufferIdx].triangleVertexNormalBuffer.push(this.gl.createBuffer());
    this.displayBuffers[this.currentBufferIdx].triangleVertexNormalBuffer[this.displayBuffers[this.currentBufferIdx].triangleVertexNormalBuffer.length-1].numItems = 0;
    for (var j=0; j<norm.length; j++){
        this.displayBuffers[this.currentBufferIdx].triangleNormals[this.displayBuffers[this.currentBufferIdx].triangleNormals.length-1].push(parseFloat(norm[j]));
        this.displayBuffers[this.currentBufferIdx].triangleVertexNormalBuffer[this.displayBuffers[this.currentBufferIdx].triangleVertexNormalBuffer.length-1].numItems++;
    }
    this.displayBuffers[this.currentBufferIdx].triangleVertexNormalBuffer[this.displayBuffers[this.currentBufferIdx].triangleVertexNormalBuffer.length-1].numItems /= 3;
}

MGWebGL.prototype.setSymmetryMatrices = function(symmetry){
    this.displayBuffers[this.currentBufferIdx].setSymmetryMatrices(symmetry);
}

MGWebGL.prototype.createColourBuffer = function(colour){
    this.displayBuffers[this.currentBufferIdx].triangleColourBuffer.push(this.gl.createBuffer());
    this.displayBuffers[this.currentBufferIdx].triangleColourBuffer[this.displayBuffers[this.currentBufferIdx].triangleColourBuffer.length-1].numItems = 0;
    this.displayBuffers[this.currentBufferIdx].triangleColours.push([])
        if(Math.abs(parseFloat(colour[3]))<0.99){
            console.log("This is transparent");
            this.displayBuffers[this.currentBufferIdx].transparent = true;
        }
    for (var j=0; j<colour.length; j++){
        this.displayBuffers[this.currentBufferIdx].triangleColours[this.displayBuffers[this.currentBufferIdx].triangleColours.length-1].push(parseFloat(colour[j]));
        this.displayBuffers[this.currentBufferIdx].triangleColourBuffer[this.displayBuffers[this.currentBufferIdx].triangleColourBuffer.length-1].numItems++;
    }
    this.displayBuffers[this.currentBufferIdx].triangleColourBuffer[this.displayBuffers[this.currentBufferIdx].triangleColourBuffer.length-1].numItems /= 4;
}

MGWebGL.prototype.addSupplementaryInfo = function(info,name){
    if(typeof(this.displayBuffers[this.currentBufferIdx].supplementary[name])==="undefined"){
      this.displayBuffers[this.currentBufferIdx].supplementary[name] = [info];
    } else {
      this.displayBuffers[this.currentBufferIdx].supplementary[name].push(info);
    }
}

MGWebGL.prototype.createIndexBuffer = function(idx){
    this.displayBuffers[this.currentBufferIdx].triangleVertexIndexBuffer.push(this.gl.createBuffer());
    this.displayBuffers[this.currentBufferIdx].triangleVertexIndexBuffer[this.displayBuffers[this.currentBufferIdx].triangleVertexIndexBuffer.length-1].numItems = 0;
    this.displayBuffers[this.currentBufferIdx].triangleIndexs.push([]);
    for (var j=0; j<idx.length; j++){
        this.displayBuffers[this.currentBufferIdx].triangleIndexs[this.displayBuffers[this.currentBufferIdx].triangleIndexs.length-1].push(parseFloat(idx[j]));
        this.displayBuffers[this.currentBufferIdx].triangleVertexIndexBuffer[this.displayBuffers[this.currentBufferIdx].triangleVertexIndexBuffer.length-1].numItems++;
    }
}

MGWebGL.prototype.createSizeBuffer = function(idx){
    this.displayBuffers[this.currentBufferIdx].primitiveSizes.push([]);
    for (var j=0; j<idx.length; j++){
        this.displayBuffers[this.currentBufferIdx].primitiveSizes[this.displayBuffers[this.currentBufferIdx].primitiveSizes.length-1].push(parseFloat(idx[j]));
    }
}
