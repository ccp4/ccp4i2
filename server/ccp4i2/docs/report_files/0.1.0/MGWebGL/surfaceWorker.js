importScripts('edtsurf.js')
importScripts('glMatrix.js')
importScripts('calculateSurface.js')

self.onmessage = function (message) {

    var msg = message.data;
    var contents = msg;

    var EDTSurfInst = EDTSurf();

    EDTSurfInst['onRuntimeInitialized'] = onRuntimeInitialized;

    function onRuntimeInitialized() {
        var surfaceProps = calculateSurface(contents,EDTSurfInst);
        var totalVerts = surfaceProps.totalVerts;
        var normals = surfaceProps.normals;
        var totalCols = surfaceProps.totalCols;
        var totalIdxs = surfaceProps.totalIdxs;
        var surfaceOrigin = surfaceProps.origin;

        console.log("Got stuff back from surface calc");

        //TODO, now need to post back totalCols, normals, totalVerts, totalIdxs and do what is below back in main thread, or perhaps even calc xtot, etc. here.
        var aBufVerts = new ArrayBuffer(totalVerts.length*4);
        var aBufNorms = new ArrayBuffer(normals.length*4);
        var aBufCols = new ArrayBuffer(totalCols.length*4);
        var aBufIdxs = new ArrayBuffer(totalIdxs.length*4);
        var aBufOrigin = new ArrayBuffer(3*4);
        console.log("Made array buffers");

        var bufViewVerts = new Float32Array(aBufVerts);
        var bufViewNorms = new Float32Array(aBufNorms);
        var bufViewCols = new Float32Array(aBufCols);
        var bufViewOrigin = new Float32Array(aBufOrigin);
        var bufViewIdxs = new Int32Array(aBufIdxs);
        console.log("Made arrays");

        for (var i=0, Len=totalVerts.length; i < Len; i++) {
            bufViewVerts[i] = totalVerts[i];
        }
        for (var i=0, Len=normals.length; i < Len; i++) {
            bufViewNorms[i] = normals[i];
        }
        for (var i=0, Len=totalCols.length; i < Len; i++) {
            bufViewCols[i] = totalCols[i];
        }
        for (var i=0, Len=totalIdxs.length; i < Len; i++) {
            bufViewIdxs[i] = totalIdxs[i];
        }
        console.log("Filled array buffers");

        bufViewOrigin[0] = surfaceOrigin[0];
        bufViewOrigin[1] = surfaceOrigin[1];
        bufViewOrigin[2] = surfaceOrigin[2];

        console.log("posting message");
        self.postMessage({aTopic:'do_sendMainArrBuff', aBufVerts:aBufVerts, aBufNorms:aBufNorms, aBufCols:aBufCols, aBufIdxs:aBufIdxs, aBufOrigin:aBufOrigin}, [aBufVerts,aBufNorms,aBufCols,aBufIdxs,aBufOrigin]);
        console.log("posted message");
    }
    
}
