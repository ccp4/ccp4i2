calculateSurface = function (contents,EDTSurfInst) {

    function chunkSubstr(str, size) {
        const numChunks = Math.ceil(str.length / size)
            const chunks = new Array(numChunks)

            for (var i = 0, o = 0; i < numChunks; ++i, o += size) {
                chunks[i] = str.substr(o, size)
            }

        return chunks
    }

    var doSplit = true;
    if(doSplit){
        var clearPDBString = EDTSurfInst.cwrap('clearPDBString', 'void',['number']);
        var EDTSurfAddPartialPDBString = EDTSurfInst.cwrap('EDTSurfAddPartialPDBString', 'void',['string']);
        var EDTSurfMainFromParts = EDTSurfInst.cwrap('EDTSurfMainFromParts', 'void',['void']);
        clearPDBString(contents.length);
        var chunks = chunkSubstr(contents,10000);
        for(var ichunk=0;ichunk<chunks.length;ichunk++){
            EDTSurfAddPartialPDBString(chunks[ichunk]);
        }
        EDTSurfMainFromParts();
    } else {
        var edtsurf = EDTSurfInst.cwrap('EDTSurfMain', 'void',['string']);
        edtsurf(contents);
    }

    console.log("Finished");
    var getVertices = EDTSurfInst.cwrap('getVertices', 'void',['number','number','number']);
    var getColors = EDTSurfInst.cwrap('getColors', 'void',['number','number','number']);
    var getIndices = EDTSurfInst.cwrap('getIndices', 'void',['number','number','number']);
    var numberOfVertices = EDTSurfInst.cwrap('numberOfVertices', 'number',['void']);
    var numberOfTriangles = EDTSurfInst.cwrap('numberOfTriangles', 'number',['void']);

    console.log(numberOfVertices());
    console.log(numberOfTriangles());

    var data = new Float32Array(Math.min(60000,numberOfVertices()*3));
    // Get data byte size, allocate memory on Emscripten heap, and get pointer
    var nDataBytes = data.length * data.BYTES_PER_ELEMENT;
    var dataPtr = EDTSurfInst._malloc(nDataBytes);

    // Copy data to Emscripten heap (directly accessed from EDTSurfInst.HEAPU8)
    var dataHeap = new Uint8Array(EDTSurfInst.HEAPU8.buffer, dataPtr, nDataBytes);
    EDTSurfInst.HEAPU8.set(data,dataPtr);

    var totalVerts = [];
    for(var ib=0;ib<numberOfVertices()*3;ib+=data.length){
        // Call function and get result
        var verts;
        if(ib+data.length<numberOfVertices()*3){
            getVertices(dataPtr, ib, data.length);
            verts = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, data.length);
        } else {
            var nleft = numberOfVertices()*3-ib;
            getVertices(dataPtr, ib, nleft);
            verts = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, nleft);
        }
        for(var ia=0,vertLength=verts.length;ia<vertLength;ia++){
            totalVerts.push(verts[ia]);
        }
    }

    // Free memory
    EDTSurfInst._free(dataHeap.byteOffset);
    console.log(totalVerts.length);

    var dataCol = new Float32Array(Math.min(60000,numberOfVertices()*4));
    // Get dataCol byte size, allocate memory on Emscripten heap, and get pointer
    var nDataColBytes = dataCol.length * dataCol.BYTES_PER_ELEMENT;
    var dataColPtr = EDTSurfInst._malloc(nDataColBytes);

    // Copy dataCol to Emscripten heap (directly accessed from EDTSurfInst.HEAPU8)
    var dataColHeap = new Uint8Array(EDTSurfInst.HEAPU8.buffer, dataColPtr, nDataColBytes);
    EDTSurfInst.HEAPU8.set(dataCol,dataColPtr);

    var totalCols = [];
    for(var ib=0;ib<numberOfVertices()*4;ib+=data.length){
        // Call function and get result
        var cols;
        if(ib+data.length<numberOfVertices()*4){
            getColors(dataColHeap.byteOffset, ib, dataCol.length);
            cols = new Float32Array(dataColHeap.buffer, dataColHeap.byteOffset, dataCol.length);
        } else {
            var nleft = numberOfVertices()*4-ib;
            getColors(dataColHeap.byteOffset, ib, nleft);
            cols = new Float32Array(dataColHeap.buffer, dataColHeap.byteOffset, nleft);
        }
        for(var ia=0,colLength=cols.length;ia<colLength;ia++){
            totalCols.push(cols[ia]);
        }
    }

    // Free memory
    EDTSurfInst._free(dataColHeap.byteOffset);
    console.log(totalCols.length);


    var dataIdx = new Int32Array(Math.min(60000,numberOfTriangles()*3));
    // Get dataIdx byte size, allocate memory on Emscripten heap, and get pointer
    var nDataIdxBytes = dataIdx.length * dataIdx.BYTES_PER_ELEMENT;
    var dataIdxPtr = EDTSurfInst._malloc(nDataIdxBytes);

    // Copy dataIdx to Emscripten heap (directly accessed from EDTSurfInst.HEAPU8)
    var dataIdxHeap = new Uint8Array(EDTSurfInst.HEAPU8.buffer, dataIdxPtr, nDataIdxBytes);
    EDTSurfInst.HEAPU8.set(dataIdx,dataIdxPtr);

    var totalIdxs = [];
    for(var ib=0;ib<numberOfTriangles()*3;ib+=data.length){
        // Call function and get result
        var idxs;
        if(ib+data.length<numberOfTriangles()*3){
            getIndices(dataIdxHeap.byteOffset, ib, dataIdx.length);
            idxs = new Int32Array(dataIdxHeap.buffer, dataIdxHeap.byteOffset, dataIdx.length);
        } else {
            var nleft =numberOfTriangles()*3-ib;
            getIndices(dataIdxHeap.byteOffset, ib, nleft);
            idxs = new Int32Array(dataIdxHeap.buffer, dataIdxHeap.byteOffset, nleft);
        }
        for(var ia=0,idxLength=idxs.length;ia<idxLength;ia++){
            totalIdxs.push(idxs[ia]);
        }
    }

    // Free memory
    EDTSurfInst._free(dataIdxHeap.byteOffset);
    console.log(totalIdxs.length);

    var normals = [];
    var normalsTriCount = [];
    var idxLength = totalIdxs.length;
    var vertLength = totalVerts.length;
    var vertLengthOver3 = totalVerts.length/3;
    for(var ivert=0;ivert<vertLength;ivert++){
        normals[ivert] = 0.0;
    }
    for(var ivert=0;ivert<vertLengthOver3;ivert++){
        normalsTriCount[ivert] = 0;
    }
    for(var idx=0;idx<idxLength;idx+=3){
        var x1 = totalVerts[totalIdxs[idx]*3];
        var y1 = totalVerts[totalIdxs[idx]*3+1];
        var z1 = totalVerts[totalIdxs[idx]*3+2];
        var x2 = totalVerts[totalIdxs[idx+1]*3];
        var y2 = totalVerts[totalIdxs[idx+1]*3+1];
        var z2 = totalVerts[totalIdxs[idx+1]*3+2];
        var x3 = totalVerts[totalIdxs[idx+2]*3];
        var y3 = totalVerts[totalIdxs[idx+2]*3+1];
        var z3 = totalVerts[totalIdxs[idx+2]*3+2];

        var c12 = vec3.create([x2-x1,y2-y1,z2-z1]);
        var c13 = vec3.create([x3-x1,y3-y1,z3-z1]);
        vec3.normalize(c12);
        vec3.normalize(c13);
        var c123 = vec3.create();
        vec3.cross(c12,c13,c123);
        vec3.normalize(c123);

        normals[totalIdxs[idx]*3]     += c123[0];
        normals[totalIdxs[idx]*3+1]   += c123[1];
        normals[totalIdxs[idx]*3+2]   += c123[2];
        normals[totalIdxs[idx+1]*3]   += c123[0];
        normals[totalIdxs[idx+1]*3+1] += c123[1];
        normals[totalIdxs[idx+1]*3+2] += c123[2];
        normals[totalIdxs[idx+2]*3]   += c123[0];
        normals[totalIdxs[idx+2]*3+1] += c123[1];
        normals[totalIdxs[idx+2]*3+2] += c123[2];

        normalsTriCount[totalIdxs[idx]]   += 1;
        normalsTriCount[totalIdxs[idx+1]] += 1;
        normalsTriCount[totalIdxs[idx+2]] += 1;
    }
    for(var ivert=0;ivert<vertLengthOver3;ivert++){
        normals[3*ivert] /= normalsTriCount[ivert];
        normals[3*ivert+1] /= normalsTriCount[ivert];
        normals[3*ivert+2] /= normalsTriCount[ivert];
    }
    var xtot = 0.0;
    var ytot = 0.0;
    var ztot = 0.0;
    for(var il=0;il<vertLengthOver3;il++){
        xtot += totalVerts[3*il];
        ytot += totalVerts[3*il+1];
        ztot += totalVerts[3*il+2];
    }
    xtot /= vertLengthOver3;
    ytot /= vertLengthOver3;
    ztot /= vertLengthOver3;
    var origin = [xtot,ytot,ztot];
    
    var ret = {};
    ret.totalVerts = totalVerts;
    ret.normals = normals;
    ret.totalCols = totalCols;
    ret.totalIdxs = totalIdxs;
    ret.origin = origin;
    console.log(normals.length);
    return ret;
}
