function stacktrace() { 
    var err = new Error();
    return err.stack;
}

function parsePlySurface(contents) {
    var lines = contents.split("\n");
    var linesLength = lines.length;
    var endHeader = -1;
    var nVerts = 0;
    var nTris = 0;
    //TODO - For now, I am assuming vertices x,y,z,r,g,b and indices i,j,k,r,g,b.
    for(var il=0;il<linesLength;il++){
        if(lines[il].substr(0,14)==="element vertex"){
            nVerts = parseInt(lines[il].substr(14));
        } else if(lines[il].substr(0,12)==="element face"){
            nTris = parseInt(lines[il].substr(12));
        } else if(lines[il].substr(0,10)==="end_header"){
            endHeader = il;
            break;
        }
    }
    var verts = [];
    var cols = [];
    var idxs = [];
    var vertStart = endHeader+1;
    var vertEnd = endHeader+nVerts;
    var triStart = endHeader+nVerts+1;
    var triEnd = endHeader+nVerts+nTris;
    var xtot = 0.0;
    var ytot = 0.0;
    var ztot = 0.0;
    for(var il=vertStart;il<=vertEnd;il++){
        var lsplit = lines[il].split(/\s+/);
        verts.push(lsplit[0]);
        verts.push(lsplit[1]);
        verts.push(lsplit[2]);
        xtot += parseFloat(lsplit[0]);
        ytot += parseFloat(lsplit[1]);
        ztot += parseFloat(lsplit[2]);
        cols.push(parseFloat(lsplit[3])/255.);
        cols.push(parseFloat(lsplit[4])/255.);
        cols.push(parseFloat(lsplit[5])/255.);
        cols.push(1.0);
    }
    xtot /= nVerts;
    ytot /= nVerts;
    ztot /= nVerts;
    for(var il=triStart;il<=triEnd;il++){
        var lsplit = lines[il].split(/\s+/);
        idxs.push(parseInt(lsplit[1]));
        idxs.push(parseInt(lsplit[2]));
        idxs.push(parseInt(lsplit[3]));
    }
    var normals = [];
    var normalsTriCount = [];
    var idxLength = idxs.length;
    var vertLength = verts.length;
    var vertLengthOver3 = verts.length/3;
    for(var ivert=0;ivert<vertLength;ivert++){
        normals[ivert] = 0.0;
    }
    for(var ivert=0;ivert<vertLengthOver3;ivert++){
        normalsTriCount[ivert] = 0;
    }
    for(var idx=0;idx<idxLength;idx+=3){
        var x1 = verts[idxs[idx]*3];
        var y1 = verts[idxs[idx]*3+1];
        var z1 = verts[idxs[idx]*3+2];
        var x2 = verts[idxs[idx+1]*3];
        var y2 = verts[idxs[idx+1]*3+1];
        var z2 = verts[idxs[idx+1]*3+2];
        var x3 = verts[idxs[idx+2]*3];
        var y3 = verts[idxs[idx+2]*3+1];
        var z3 = verts[idxs[idx+2]*3+2];

        var c12 = vec3.create([x2-x1,y2-y1,z2-z1]);
        var c13 = vec3.create([x3-x1,y3-y1,z3-z1]);
        vec3.normalize(c12);
        vec3.normalize(c13);
        var c123 = vec3.create();
        vec3.cross(c12,c13,c123);
        vec3.normalize(c123);

        normals[idxs[idx]*3]     += c123[0];
        normals[idxs[idx]*3+1]   += c123[1];
        normals[idxs[idx]*3+2]   += c123[2];
        normals[idxs[idx+1]*3]   += c123[0];
        normals[idxs[idx+1]*3+1] += c123[1];
        normals[idxs[idx+1]*3+2] += c123[2];
        normals[idxs[idx+2]*3]   += c123[0];
        normals[idxs[idx+2]*3+1] += c123[1];
        normals[idxs[idx+2]*3+2] += c123[2];

        normalsTriCount[idxs[idx]]   += 1;
        normalsTriCount[idxs[idx+1]] += 1;
        normalsTriCount[idxs[idx+2]] += 1;
    }
    for(var ivert=0;ivert<vertLengthOver3;ivert++){
        normals[3*ivert] /= normalsTriCount[ivert];
        normals[3*ivert+1] /= normalsTriCount[ivert];
        normals[3*ivert+2] /= normalsTriCount[ivert];
    }
    var surfTriangleData = {"col_tri":[[cols]], "norm_tri":[[normals]], "vert_tri":[[verts]], "idx_tri":[[idxs]], "prim_types":[["TRIANGLES"]] };
    return [surfTriangleData,[xtot,ytot,ztot]];
}

function getSugarAtomNames(){

    var ribAtomNames = [];
    var nagAtomNames = [];
    var ffcAtomNames = [];
    var latAtomNames = [];
    var lmtAtomNames = [];
    var treAtomNames = [];
    var dmuAtomNames = [];
    var siaAtomNames = [];
    var a2gAtomNames = [];

    var ribR1AtomNames = [];
    var nagR1AtomNames = [];
    var ffcR1AtomNames = [];
    var latR1AtomNames = [];
    var lmtR1AtomNames = [];
    var treR1AtomNames = [];
    var dmuR1AtomNames = [];
    var ffcR2AtomNames = [];
    var latR2AtomNames = [];
    var lmtR2AtomNames = [];
    var treR2AtomNames = [];
    var dmuR2AtomNames = [];
    var siaR1AtomNames = [];
    var a2gR1AtomNames = [];

    var sugarAtomNames = [];

    ffcR1AtomNames.push("C1A");
    ffcR1AtomNames.push("C2A");
    ffcR1AtomNames.push("C3A");
    ffcR1AtomNames.push("C4A");
    ffcR1AtomNames.push("C5A");
    ffcR1AtomNames.push("O5A");
    ffcR2AtomNames.push("C1B");
    ffcR2AtomNames.push("C2B");
    ffcR2AtomNames.push("C3B");
    ffcR2AtomNames.push("C4B");
    ffcR2AtomNames.push("C5B");
    ffcR2AtomNames.push("O5B");

    ffcAtomNames.push(ffcR1AtomNames);
    ffcAtomNames.push(ffcR2AtomNames);
    sugarAtomNames.push(ffcAtomNames);

    latR1AtomNames.push("C1");
    latR1AtomNames.push("C2");
    latR1AtomNames.push("C3");
    latR1AtomNames.push("C4");
    latR1AtomNames.push("C5");
    latR1AtomNames.push("O5");
    latR2AtomNames.push("C1'");
    latR2AtomNames.push("C2'");
    latR2AtomNames.push("C3'");
    latR2AtomNames.push("C4'");
    latR2AtomNames.push("C5'");
    latR2AtomNames.push("O5'");

    latAtomNames.push(latR1AtomNames);
    latAtomNames.push(latR2AtomNames);
    sugarAtomNames.push(latAtomNames);

    lmtR1AtomNames.push("C1'");
    lmtR1AtomNames.push("C2'");
    lmtR1AtomNames.push("C3'");
    lmtR1AtomNames.push("C4'");
    lmtR1AtomNames.push("C5'");
    lmtR1AtomNames.push("O5'");
    lmtR2AtomNames.push("C1B");
    lmtR2AtomNames.push("C2B");
    lmtR2AtomNames.push("C3B");
    lmtR2AtomNames.push("C4B");
    lmtR2AtomNames.push("C5B");
    lmtR2AtomNames.push("O5B");

    lmtAtomNames.push(lmtR1AtomNames);
    lmtAtomNames.push(lmtR2AtomNames);
    sugarAtomNames.push(lmtAtomNames);

    treR1AtomNames.push("C1P");
    treR1AtomNames.push("C2P");
    treR1AtomNames.push("C3P");
    treR1AtomNames.push("C4P");
    treR1AtomNames.push("C5P");
    treR1AtomNames.push("O5P");
    treR2AtomNames.push("C1B");
    treR2AtomNames.push("C2B");
    treR2AtomNames.push("C3B");
    treR2AtomNames.push("C4B");
    treR2AtomNames.push("C5B");
    treR2AtomNames.push("O5B");

    treAtomNames.push(treR1AtomNames);
    treAtomNames.push(treR2AtomNames);
    sugarAtomNames.push(treAtomNames);

    ribR1AtomNames.push("C1'");
    ribR1AtomNames.push("C2'");
    ribR1AtomNames.push("C3'");
    ribR1AtomNames.push("C4'");
    ribR1AtomNames.push("O4'");

    ribAtomNames.push(ribR1AtomNames);
    sugarAtomNames.push(ribAtomNames);

    nagR1AtomNames.push("C1");
    nagR1AtomNames.push("C2");
    nagR1AtomNames.push("C3");
    nagR1AtomNames.push("C4");
    nagR1AtomNames.push("C5");
    nagR1AtomNames.push("O5");

    nagAtomNames.push(nagR1AtomNames);
    sugarAtomNames.push(nagAtomNames);

    siaR1AtomNames.push("C2");
    siaR1AtomNames.push("C3");
    siaR1AtomNames.push("C4");
    siaR1AtomNames.push("C5");
    siaR1AtomNames.push("C6");
    siaR1AtomNames.push("O6");

    siaAtomNames.push(siaR1AtomNames);
    sugarAtomNames.push(siaAtomNames);

    a2gR1AtomNames.push("C1");
    a2gR1AtomNames.push("C2");
    a2gR1AtomNames.push("C3");
    a2gR1AtomNames.push("C4");
    a2gR1AtomNames.push("C5");
    a2gR1AtomNames.push("O");

    a2gAtomNames.push(a2gR1AtomNames);
    sugarAtomNames.push(a2gAtomNames);

    return sugarAtomNames;
}

function singletonsToLinesInfo(singletons,size,colourScheme){
    var cylinder_sizes = [];
    var cylinder_col_tri = [];
    var cylinder_vert_tri = [];
    var cylinder_idx_tri = [];
    var sphere_atoms = [];

    for(var ic=0;ic<singletons.length;ic++){
        var at1 = singletons[ic];
        var atom = {};
        atom["x"] = at1.x();
        atom["y"] = at1.y();
        atom["z"] = at1.z();
        atom["tempFactor"] = at1["_atom_site.B_iso_or_equiv"];
        atom["charge"] = at1["_atom_site.pdbx_formal_charge"];
        atom["label"] =  at1.getAtomID();
        sphere_atoms.push(atom);

        cylinder_idx_tri.push(2*ic);
        cylinder_vert_tri.push(at1.x()-.1);
        cylinder_vert_tri.push(at1.y());
        cylinder_vert_tri.push(at1.z());
        cylinder_vert_tri.push(at1.x()+.1);
        cylinder_vert_tri.push(at1.y());
        cylinder_vert_tri.push(at1.z());
        for(ip=0;ip<colourScheme[at1["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at1["_atom_site.id"]][ip]);
        for(ip=0;ip<colourScheme[at1["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at1["_atom_site.id"]][ip]);
        cylinder_sizes.push(size);

        cylinder_idx_tri.push(2*ic+1);
        cylinder_vert_tri.push(at1.x());
        cylinder_vert_tri.push(at1.y()-.1);
        cylinder_vert_tri.push(at1.z());
        cylinder_vert_tri.push(at1.x());
        cylinder_vert_tri.push(at1.y()+.1);
        cylinder_vert_tri.push(at1.z());
        for(ip=0;ip<colourScheme[at1["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at1["_atom_site.id"]][ip]);
        for(ip=0;ip<colourScheme[at1["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at1["_atom_site.id"]][ip]);
        cylinder_sizes.push(size);

        cylinder_idx_tri.push(2*ic+1);
        cylinder_vert_tri.push(at1.x());
        cylinder_vert_tri.push(at1.y());
        cylinder_vert_tri.push(at1.z()-.1);
        cylinder_vert_tri.push(at1.x());
        cylinder_vert_tri.push(at1.y());
        cylinder_vert_tri.push(at1.z()+.1);
        for(ip=0;ip<colourScheme[at1["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at1["_atom_site.id"]][ip]);
        for(ip=0;ip<colourScheme[at1["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at1["_atom_site.id"]][ip]);
        cylinder_sizes.push(size);
 
    }

    var cylinderPrimitiveInfo = {"atoms":[[sphere_atoms]],"sizes": [[cylinder_sizes]], "col_tri":[[cylinder_col_tri]], "norm_tri":[[[]]], "vert_tri":[[cylinder_vert_tri]], "idx_tri":[[cylinder_idx_tri]] , "prim_types":[["LINES"]] };
    return cylinderPrimitiveInfo;
}

function contactsToCylindersLinesInfo(contacts,size,linetype,colourScheme,dashed){
    //console.log(contacts.length);
    //console.log(contacts);
    var cylinder_sizes = [];
    var cylinder_col_tri = [];
    var cylinder_vert_tri = [];
    var cylinder_idx_tri = [];
    var sphere_atoms = [];
    if(typeof(dashed)!=="undefined" && dashed){
        var dashLength = 0.2;
        var idx = 0;
        for(var ic=0;ic<contacts.length;ic++){
            var at1 = contacts[ic][1];
            var at2 = contacts[ic][2];
            if(typeof(at1.residue)!=="undefined"&&typeof(at2.residue)!=="undefined"){
                var atom = {};
                atom["x"] = at1.x();
                atom["y"] = at1.y();
                atom["z"] = at1.z();
                atom["tempFactor"] = at1["_atom_site.B_iso_or_equiv"];
                atom["charge"] = at1["_atom_site.pdbx_formal_charge"];
                atom["label"] =  at1.getAtomID();
                sphere_atoms.push(atom);
                var atom2 = {};
                atom2["x"] = at2.x();
                atom2["y"] = at2.y();
                atom2["z"] = at2.z();
                atom2["tempFactor"] = at2["_atom_site.B_iso_or_equiv"];
                atom2["charge"] = at2["_atom_site.pdbx_formal_charge"];
                atom2["label"] =  at2.getAtomID();
                sphere_atoms.push(atom2);
            }
            var midpoint = [.5*(at1.x()+at2.x()),.5*(at1.y()+at2.y()),.5*(at1.z()+at2.z())];
            var thisLength = Model.prototype.bondLength(at1,at2);
            var nDashes = parseInt(0.5*thisLength/dashLength);
            var thisDashLength = thisLength/nDashes;
            var fracAdd = thisDashLength/thisLength;
            var frac = 1.0-.5*fracAdd;
            var ifrac = 0;
            while(frac>0.0){
                var frac2 = frac+fracAdd;
                if(frac2>1.0) frac2 = 1.0;
                if(frac>frac2) frac = frac2;
                cylinder_vert_tri.push((1.0-frac)*at1.x()+frac*midpoint[0]);
                cylinder_vert_tri.push((1.0-frac)*at1.y()+frac*midpoint[1]);
                cylinder_vert_tri.push((1.0-frac)*at1.z()+frac*midpoint[2]);
                cylinder_vert_tri.push((1.0-frac2)*at1.x()+frac2*midpoint[0]);
                cylinder_vert_tri.push((1.0-frac2)*at1.y()+frac2*midpoint[1]);
                cylinder_vert_tri.push((1.0-frac2)*at1.z()+frac2*midpoint[2]);
                Array.prototype.push.apply(cylinder_col_tri,colourScheme[at1["_atom_site.id"]]);
                Array.prototype.push.apply(cylinder_col_tri,colourScheme[at1["_atom_site.id"]]);
                cylinder_sizes.push(size);
                cylinder_idx_tri.push(idx++);
                cylinder_idx_tri.push(idx++);
                frac -= 2.0*fracAdd;
                ifrac += 1;
            }
            fracAdd = thisDashLength/thisLength;
            frac = 1.0-.5*fracAdd;
            ifrac = 0;
            while(frac>0.0){
                var frac2 = frac+fracAdd;
                if(frac2>1.0) frac2 = 1.0;
                if(frac>frac2) frac = frac2;
                cylinder_vert_tri.push((1.0-frac)*at2.x()+frac*midpoint[0]);
                cylinder_vert_tri.push((1.0-frac)*at2.y()+frac*midpoint[1]);
                cylinder_vert_tri.push((1.0-frac)*at2.z()+frac*midpoint[2]);
                cylinder_vert_tri.push((1.0-frac2)*at2.x()+frac2*midpoint[0]);
                cylinder_vert_tri.push((1.0-frac2)*at2.y()+frac2*midpoint[1]);
                cylinder_vert_tri.push((1.0-frac2)*at2.z()+frac2*midpoint[2]);
                Array.prototype.push.apply(cylinder_col_tri,colourScheme[at2["_atom_site.id"]]);
                Array.prototype.push.apply(cylinder_col_tri,colourScheme[at2["_atom_site.id"]]);
                cylinder_sizes.push(size);
                cylinder_idx_tri.push(idx++);
                cylinder_idx_tri.push(idx++);
                frac -= 2.0*fracAdd;
                ifrac += 1;
            }
        }
    } else {
        for(var ic=0;ic<contacts.length;ic++){
            var at1 = contacts[ic][1];
            var at2 = contacts[ic][2];
            if(typeof(at1.residue)!=="undefined"&&typeof(at2.residue)!=="undefined"){
                var atom = {};
                atom["x"] = at1.x();
                atom["y"] = at1.y();
                atom["z"] = at1.z();
                atom["tempFactor"] = at1["_atom_site.B_iso_or_equiv"];
                atom["charge"] = at1["_atom_site.pdbx_formal_charge"];
                atom["label"] =  at1.getAtomID();
                sphere_atoms.push(atom);
                var atom2 = {};
                atom2["x"] = at2.x();
                atom2["y"] = at2.y();
                atom2["z"] = at2.z();
                atom2["tempFactor"] = at2["_atom_site.B_iso_or_equiv"];
                atom2["charge"] = at2["_atom_site.pdbx_formal_charge"];
                atom2["label"] =  at2.getAtomID();
                sphere_atoms.push(atom2);
            }
            cylinder_idx_tri.push(2*ic);
            cylinder_idx_tri.push(2*ic+1);
            cylinder_vert_tri.push(at2.x());
            cylinder_vert_tri.push(at2.y());
            cylinder_vert_tri.push(at2.z());
            var midpoint = [.5*(at1.x()+at2.x()),.5*(at1.y()+at2.y()),.5*(at1.z()+at2.z())];
            cylinder_vert_tri.push(midpoint[0]);
            cylinder_vert_tri.push(midpoint[1]);
            cylinder_vert_tri.push(midpoint[2]);
            for(ip=0;ip<colourScheme[at2["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at2["_atom_site.id"]][ip]);
            for(ip=0;ip<colourScheme[at2["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at2["_atom_site.id"]][ip]);
            cylinder_sizes.push(size);
            cylinder_vert_tri.push(at1.x());
            cylinder_vert_tri.push(at1.y());
            cylinder_vert_tri.push(at1.z());
            var midpoint = [.5*(at1.x()+at2.x()),.5*(at1.y()+at2.y()),.5*(at1.z()+at2.z())];
            cylinder_vert_tri.push(midpoint[0]);
            cylinder_vert_tri.push(midpoint[1]);
            cylinder_vert_tri.push(midpoint[2]);
            for(ip=0;ip<colourScheme[at1["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at1["_atom_site.id"]][ip]);
            for(ip=0;ip<colourScheme[at1["_atom_site.id"]].length;ip++) cylinder_col_tri.push(colourScheme[at1["_atom_site.id"]][ip]);
            cylinder_sizes.push(size);
        }
    }

    var cylinderPrimitiveInfo = {"atoms":[[sphere_atoms]],"sizes": [[cylinder_sizes]], "col_tri":[[cylinder_col_tri]], "norm_tri":[[[]]], "vert_tri":[[cylinder_vert_tri]], "idx_tri":[[cylinder_idx_tri]] , "prim_types":[[linetype]] };
    return cylinderPrimitiveInfo;
}

function contactsToLinesInfo(contacts,size,colourScheme,dashed){
    return contactsToCylindersLinesInfo(contacts,size,"LINES",colourScheme,dashed);
}

function contactsToCylindersInfo(contacts,size,colourScheme,dashed){
    return contactsToCylindersLinesInfo(contacts,size,"CYLINDERS",colourScheme,dashed);
}

function contactsToCappedCylindersInfo(contacts,size,colourScheme,dashed){
    return contactsToCylindersLinesInfo(contacts,size,"CAPCYLINDERS",colourScheme,dashed);
}

function DrawSugarBlockInt(res, col1, at1,at2,at3,at4,at5,at6, two_colour_in, sugar_block_thickness, sugar_block_scale, block_vert_tri,block_norm_tri,block_col_tri,block_idx_tri ){
    var retval = [];
    if(at1&&at2&&at3&&at4&&at5){
        //console.log("DrawSugarBlockInt have at least 5 atoms");
        var xaxis = vec3.create([1.0,0.0,0.0]);
        var yaxis = vec3.create([0.0,1.0,0.0]);
        var zaxis = vec3.create([0.0,0.0,1.0]);
        var norm_vec = [];
        var mid_vec = [];
        var cat1 = vec3.create([at1.x(),at1.y(),at1.z()]);
        var cat2 = vec3.create([at2.x(),at2.y(),at2.z()]);
        var cat3 = vec3.create([at3.x(),at3.y(),at3.z()]);
        var cat4 = vec3.create([at4.x(),at4.y(),at4.z()]);
        var cat5 = vec3.create([at5.x(),at5.y(),at5.z()]);
        //console.log(cat1);
        //console.log(cat2);
        //console.log(cat3);
        var cat2cat1 = vec3.create();
        var cat3cat2 = vec3.create();
        var cat4cat3 = vec3.create();
        var cat5cat4 = vec3.create();
        vec3.subtract(cat2,cat1,cat2cat1);
        vec3.subtract(cat3,cat2,cat3cat2);
        vec3.subtract(cat4,cat3,cat4cat3);
        vec3.subtract(cat5,cat4,cat5cat4);
        //console.log(cat2cat1);
        //console.log(cat3cat2);
        var d12 = vec3.length(cat2cat1);
        var d23 = vec3.length(cat3cat2);
        var d34 = vec3.length(cat4cat3);
        var d45 = vec3.length(cat5cat4);
        //console.log(d12);
        //console.log(d23);
        //console.log(d34);
        //console.log(d45);
        var c123 = vec3.create();
        var c234 = vec3.create();
        var c345 = vec3.create();
        vec3.cross(cat2cat1,cat3cat2,c123);
        vec3.cross(cat3cat2,cat4cat3,c234);
        vec3.cross(cat4cat3,cat5cat4,c345);
        if((d12>1.0&&d12<3.2)&&(d23>1.0&&d23<3.2)&&(d34>1.0&&d34<3.2)&&(d45>1.0&&d45<3.2)){
            if(at6){
                //console.log("Six membered test.");
                var cat6 = vec3.create([at6.x(),at6.y(),at6.z()]);
                var cat6cat5 = vec3.create();
                var cat1cat6 = vec3.create();
                vec3.subtract(cat6,cat5,cat6cat5);
                vec3.subtract(cat6,cat1,cat1cat6);
                var d56 = vec3.length(cat6cat5);
                var d61 = vec3.length(cat1cat6);
                if((d61>1.0&&d61<3.2)&&(d56>1.0&&d56<3.2)){
                    norm_vec.push(c123);
                    norm_vec.push(c234);
                    norm_vec.push(c345);
                    mid_vec.push(cat1);
                    mid_vec.push(cat2);
                    mid_vec.push(cat3);
                    mid_vec.push(cat4);
                    mid_vec.push(cat5);
                    var c456 = vec3.create();
                    var c561 = vec3.create();
                    vec3.cross(cat5cat4,cat6cat5,c456);
                    vec3.cross(cat6cat5,cat1cat6,c561);
                    norm_vec.push(c456);
                    norm_vec.push(c561);
                    mid_vec.push(cat6);
                    //console.log("We have 6 membered ring");
                } else if(d61>4.0&&d61<7.0&&res.getName()==="XLS") {
                    norm_vec.push(c123);
                    norm_vec.push(c234);
                    norm_vec.push(c345);
                    mid_vec.push(cat1);
                    mid_vec.push(cat2);
                    mid_vec.push(cat3);
                    mid_vec.push(cat4);
                    mid_vec.push(cat5);
                    var c456 = vec3.create();
                    var c561 = vec3.create();
                    vec3.cross(cat5cat4,cat6cat5,c456);
                    vec3.cross(cat6cat5,cat1cat6,c561);
                    norm_vec.push(c456);
                    norm_vec.push(c561);
                    mid_vec.push(cat6);
                    //console.log("Linear?");
                }
            } else {
                //console.log("Five membered test.");
                var cat5cat1 = vec3.create();
                vec3.subtract(cat5,cat1,cat5cat1);
                var d51 = vec3.length(cat5cat1);
                if((d51>1.0&&d51<3.2)){
                    norm_vec.push(c123);
                    norm_vec.push(c234);
                    norm_vec.push(c345);
                    mid_vec.push(cat1);
                    mid_vec.push(cat2);
                    mid_vec.push(cat3);
                    mid_vec.push(cat4);
                    mid_vec.push(cat5);
                    var c451 = vec3.create();
                    vec3.cross(cat5cat4,cat5cat1,c451);
                    norm_vec.push(c451);
                    //console.log("We have 5 membered ring");
                }
            }
        } else {
            console.log("Fail first test");
        }

        var nsectors = 20;
        var sugar_block_radius = 1.4 * sugar_block_scale;

        if(norm_vec.length>0&&mid_vec.length>0){

            var normal = vec3.create([0.0,0.0,0.0]);
            for(var i=0;i<norm_vec.length;i++){
                vec3.add(normal,norm_vec[i]);
            }
            vec3.normalize(normal);

            var centre = vec3.create([0.0,0.0,0.0]);
            for(var i=0;i<mid_vec.length;i++){
                vec3.add(centre,mid_vec[i]);
            }
            centre[0] /= mid_vec.length;
            centre[1] /= mid_vec.length;
            centre[2] /= mid_vec.length;

            retval.push(centre);
            retval.push(normal);
            res["GLYCO_BLOCK_CENTRE"] = [centre[0],centre[1],centre[2]];

            var p = vec3.create([normal[0],normal[1],normal[2]]);
            var up = vec3.create();
            var right = vec3.create();
            if(Math.abs(vec3.dot(xaxis,p))<0.95){
                vec3.cross(xaxis,p,up);
                vec3.normalize(up);
            }else if(Math.abs(vec3.dot(yaxis,p))<0.95){
                vec3.cross(yaxis,p,up);
                vec3.normalize(up);
            } else {
                vec3.cross(zaxis,p,up);
                vec3.normalize(up);
            }
            vec3.cross(up,p,right);

            vec3.normalize(up);
            vec3.normalize(right);

            right[0] *= sugar_block_radius;
            right[1] *= sugar_block_radius;
            right[2] *= sugar_block_radius;
            // Why??
            up[0] *= sugar_block_radius;
            up[1] *= sugar_block_radius;
            up[2] *= sugar_block_radius;

            function setGlycoBlockUDDCircle(udd,atom){
                var cp = vec3.create();
                var cpx = vec3.create();
                var cpp = vec3.create();
                var pp = vec3.create();

                vec3.subtract(centre,atom,cp);
                vec3.normalize(cp);
                vec3.cross(cp,normal,cpx);
                vec3.normalize(cpx);
                vec3.cross(normal,cpx,cpp);
                vec3.normalize(cpp);
                cpp[0] *= vec3.length(right);
                cpp[1] *= vec3.length(right);
                cpp[2] *= vec3.length(right);
                vec3.subtract(centre,cpp,pp);
                res[udd] = [pp[0],pp[1],pp[2]];
                
            }

            var white = [1.0,1.0,1.0,col1[3]];

            if(res.getName()==="NAG"||res.getName()==="NBG"||res.getName()==="NGA"||res.getName()==="NG1"||res.getName()==="NG6"||res.getName()==="GCS"||res.getName()==="6MN"||res.getName()==="GLP"||res.getName()==="GP4"||res.getName()==="A2G"){
                /* Draw a square */
                var diagonal = false;
                if(res.getName()==="GCS") {
                    diagonal = true;
                }
                console.log("Draw a square!");
                var two_colour = two_colour_in && diagonal;

                var c1c4 = vec3.create();
                vec3.subtract(cat1,cat4,c1c4);
                vec3.normalize(c1c4);

                var c1c4Perp = vec3.create();
                vec3.cross(normal,c1c4,c1c4Perp);
                vec3.normalize(c1c4Perp);

                var norm2 = vec3.create();
                vec3.cross(c1c4,c1c4Perp,norm2);

                c1c4[0] *= sugar_block_radius;
                c1c4[1] *= sugar_block_radius;
                c1c4[2] *= sugar_block_radius;
                c1c4Perp[0] *= sugar_block_radius;
                c1c4Perp[1] *= sugar_block_radius;
                c1c4Perp[2] *= sugar_block_radius;

                var ls = vec3.create();
                var le = vec3.create();
                vec3.add(centre,c1c4,ls);
                vec3.add(ls,c1c4Perp);
                vec3.add(centre,c1c4,le);
                vec3.subtract(le,c1c4Perp);

                var t = DistanceBetweenPointAndLine(ls,le,cat1)[1];
                var pp = vec3.create();
                var lels = vec3.create();
                vec3.subtract(le,ls,lels);
                lels[0] *= t;
                lels[1] *= t;
                lels[2] *= t;
                vec3.add(ls,lels,pp);
                res["GLYCO_BLOCK_C1"] = [pp[0],pp[1],pp[2]];

                var ls3 = vec3.create();
                var le3 = vec3.create();
                vec3.add(centre,c1c4,ls3);
                vec3.add(ls3,c1c4Perp);
                vec3.subtract(centre,c1c4,le3);
                vec3.add(le3,c1c4Perp);

                var t3 = DistanceBetweenPointAndLine(ls3,le3,cat3)[1];
                var pp3 = vec3.create();
                var lels3 = vec3.create();
                vec3.subtract(le3,ls3,lels3);
                lels3[0] *= t;
                lels3[1] *= t;
                lels3[2] *= t;
                vec3.add(ls3,lels3,pp3);
                res["GLYCO_BLOCK_C3"] = [pp3[0],pp3[1],pp3[2]];
 
                var ls4 = vec3.create();
                var le4 = vec3.create();
                vec3.subtract(centre,c1c4,ls4);
                vec3.add(ls4,c1c4Perp);
                vec3.subtract(centre,c1c4,le4);
                vec3.subtract(le4,c1c4Perp);

                var t4 = DistanceBetweenPointAndLine(ls4,le4,cat4)[1];
                var pp4 = vec3.create();
                var lels4 = vec3.create();
                vec3.subtract(le4,ls4,lels4);
                lels4[0] *= t;
                lels4[1] *= t;
                lels4[2] *= t;
                vec3.add(ls4,lels4,pp4);
                res["GLYCO_BLOCK_C4"] = [pp4[0],pp4[1],pp4[2]];

                var ls5 = vec3.create();
                var le5 = vec3.create();
                vec3.add(centre,c1c4,ls5);
                vec3.subtract(ls5,c1c4Perp);
                vec3.subtract(centre,c1c4,le5);
                vec3.subtract(le5,c1c4Perp);

                var t5 = DistanceBetweenPointAndLine(ls5,le5,cat5)[1];
                var pp5 = vec3.create();
                var lels5 = vec3.create();
                vec3.subtract(le5,ls5,lels5);
                lels5[0] *= t;
                lels5[1] *= t;
                lels5[2] *= t;
                vec3.add(ls5,lels5,pp5);
                res["GLYCO_BLOCK_C5"] = [pp5[0],pp5[1],pp5[2]];
                
                block_vert_tri.push(centre[0] - c1c4[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                if(two_colour){
                    block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                }else{
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                }

                block_vert_tri.push(centre[0] - c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                if(two_colour){
                    block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                }else{
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                }

                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] - c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] - c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] - c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] - c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_norm_tri.push(-c1c4[0]); block_norm_tri.push(-c1c4[1]); block_norm_tri.push(-c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

            } else if(res.getName()==="BMA"||res.getName()==="MAN"||res.getName()==="GAL"||res.getName()==="GLC"||res.getName()==="BGC"){
                console.log("CIRCLE!!!!!!!!!!!");

                setGlycoBlockUDDCircle("GLYCO_BLOCK_C1",cat1);
                setGlycoBlockUDDCircle("GLYCO_BLOCK_C3",cat3);
                setGlycoBlockUDDCircle("GLYCO_BLOCK_C4",cat4);
                setGlycoBlockUDDCircle("GLYCO_BLOCK_C5",cat5);
                setGlycoBlockUDDCircle("GLYCO_BLOCK_C2",cat2);

                for(var j=0;j<360;j=j+360/nsectors){
                    var theta = (1.0*j)/360.0 * Math.PI*2;
                    var theta2 = 1.0*(j+360/nsectors)/360.0 * Math.PI*2;
                    var x1 = Math.cos(theta);
                    var y1 = Math.sin(theta);
                    var x2 = Math.cos(theta2);
                    var y2 = Math.sin(theta2);

                    // "Cylinder" Normals
                    var n1x = x1*up[0] + y1*right[0];
                    var n1y = x1*up[1] + y1*right[1];
                    var n1z = x1*up[2] + y1*right[2];
                    var n2x = x2*up[0] + y2*right[0];
                    var n2y = x2*up[1] + y2*right[1];
                    var n2z = x2*up[2] + y2*right[2];

                    //block_vert_tri,block_norm_tri,block_col_tri,block_idx_tri

                    block_vert_tri.push(centre[0]+n1x - sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n1y - sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n1z - sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0]+n2x + sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n2y + sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n2z + sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0]+n2x - sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n2y - sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n2z - sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(n1x); block_norm_tri.push(n1y); block_norm_tri.push(n1z);
                    block_norm_tri.push(n2x); block_norm_tri.push(n2y); block_norm_tri.push(n2z);
                    block_norm_tri.push(n2x); block_norm_tri.push(n2y); block_norm_tri.push(n2z);

                    block_vert_tri.push(centre[0]+n1x - sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n1y - sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n1z - sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0]+n1x + sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n1y + sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n1z + sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0]+n2x + sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n2y + sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n2z + sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(n1x); block_norm_tri.push(n1y); block_norm_tri.push(n1z);
                    block_norm_tri.push(n1x); block_norm_tri.push(n1y); block_norm_tri.push(n1z);
                    block_norm_tri.push(n2x); block_norm_tri.push(n2y); block_norm_tri.push(n2z);

                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                    block_vert_tri.push(centre[0] - sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1] - sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2] - sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0]+n1x - sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n1y - sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n1z - sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0]+n2x - sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n2y - sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n2z - sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);

                    block_vert_tri.push(centre[0] + sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1] + sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2] + sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0]+n1x + sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n1y + sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n1z + sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0]+n2x + sugar_block_thickness*normal[0]);
                    block_vert_tri.push(centre[1]+n2y + sugar_block_thickness*normal[1]);
                    block_vert_tri.push(centre[2]+n2z + sugar_block_thickness*normal[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);

                    block_norm_tri.push(-normal[0]); block_norm_tri.push(-normal[1]); block_norm_tri.push(-normal[2]);
                    block_norm_tri.push(-normal[0]); block_norm_tri.push(-normal[1]); block_norm_tri.push(-normal[2]);
                    block_norm_tri.push(-normal[0]); block_norm_tri.push(-normal[1]); block_norm_tri.push(-normal[2]);
                    block_norm_tri.push(-normal[0]); block_norm_tri.push(-normal[1]); block_norm_tri.push(-normal[2]);
                    block_norm_tri.push(-normal[0]); block_norm_tri.push(-normal[1]); block_norm_tri.push(-normal[2]);
                    block_norm_tri.push(-normal[0]); block_norm_tri.push(-normal[1]); block_norm_tri.push(-normal[2]);

                    // FIXME - real colour
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                }

            } else if(res.getName()==="FCA"||res.getName()==="FCA"||res.getName()==="FUC"){
                /* Draw a triangle */
                console.log("Draw a triangle!");
                // Try 4byh

                var c1c4 = vec3.create();
                vec3.subtract(cat1,cat4,c1c4);
                vec3.normalize(c1c4);

                var c1c4Perp = vec3.create();
                vec3.cross(normal,c1c4,c1c4Perp);
                vec3.normalize(c1c4Perp);

                var norm2 = vec3.create();
                vec3.cross(c1c4,c1c4Perp,norm2);

                c1c4[0] *= sugar_block_radius;
                c1c4[1] *= sugar_block_radius;
                c1c4[2] *= sugar_block_radius;
                c1c4Perp[0] *= sugar_block_radius;
                c1c4Perp[1] *= sugar_block_radius;
                c1c4Perp[2] *= sugar_block_radius;

                var ytri = -0.5;
                var xtri = -0.8660254037844387;

                var ptri = vec3.create();
                var mtri = vec3.create();

                ptri[0] = ytri*c1c4Perp[0] + xtri*c1c4[0];
                ptri[1] = ytri*c1c4Perp[1] + xtri*c1c4[1];
                ptri[2] = ytri*c1c4Perp[2] + xtri*c1c4[2];

                mtri[0] = ytri*c1c4Perp[0] - xtri*c1c4[0];
                mtri[1] = ytri*c1c4Perp[1] - xtri*c1c4[1];
                mtri[2] = ytri*c1c4Perp[2] - xtri*c1c4[2];

                var normp = vec3.create();
                var normm = vec3.create();
                normp[0] = 0.5 * c1c4Perp[0] + 0.5 * ptri[0];
                normp[1] = 0.5 * c1c4Perp[1] + 0.5 * ptri[1];
                normp[2] = 0.5 * c1c4Perp[2] + 0.5 * ptri[2];
                normm[0] = 0.5 * c1c4Perp[0] + 0.5 * mtri[0];
                normm[1] = 0.5 * c1c4Perp[1] + 0.5 * mtri[1];
                normm[2] = 0.5 * c1c4Perp[2] + 0.5 * mtri[2];

                var ls = vec3.create();
                var le = vec3.create();

                vec3.add(centre,mtri,ls);
                vec3.add(centre,c1c4Perp,le);

                var t = DistanceBetweenPointAndLine(ls,le,cat1)[1];
                var pp = vec3.create();
                var lels = vec3.create();
                vec3.subtract(le,ls,lels);
                lels[0] *= t;
                lels[1] *= t;
                lels[2] *= t;
                vec3.add(ls,lels,pp);
                res["GLYCO_BLOCK_C1"] = [pp[0],pp[1],pp[2]];

                var t2 = DistanceBetweenPointAndLine(ls,le,cat2)[1];
                var pp2 = vec3.create();
                var lels2 = vec3.create();
                vec3.subtract(le,ls,lels2);
                lels2[0] *= t2;
                lels2[1] *= t2;
                lels2[2] *= t2;
                vec3.add(ls,lels2,pp2);
                res["GLYCO_BLOCK_C2"] = [pp2[0],pp2[1],pp2[2]];

                vec3.add(centre,ptri,ls);
                vec3.add(centre,c1c4Perp,le);

                var t3 = DistanceBetweenPointAndLine(ls,le,cat3)[1];
                var pp3 = vec3.create();
                var lels3 = vec3.create();
                vec3.subtract(le,ls,lels3);
                lels3[0] *= t3;
                lels3[1] *= t3;
                lels3[2] *= t3;
                vec3.add(ls,lels3,pp3);
                res["GLYCO_BLOCK_C3"] = [pp3[0],pp3[1],pp3[2]];

                var t4 = DistanceBetweenPointAndLine(ls,le,cat4)[1];
                var pp4 = vec3.create();
                var lels4 = vec3.create();
                vec3.subtract(le,ls,lels4);
                lels4[0] *= t4;
                lels4[1] *= t4;
                lels4[3] *= t4;
                vec3.add(ls,lels4,pp4);
                res["GLYCO_BLOCK_C4"] = [pp4[0],pp4[1],pp4[2]];
                
                vec3.add(centre,ptri,ls);
                vec3.add(centre,mtri,le);

                var t5 = DistanceBetweenPointAndLine(ls,le,cat5)[1];
                var pp5 = vec3.create();
                var lels5 = vec3.create();
                vec3.subtract(le,ls,lels5);
                lels5[0] *= t5;
                lels5[1] *= t5;
                lels5[3] *= t5;
                vec3.add(ls,lels5,pp5);
                res["GLYCO_BLOCK_C5"] = [pp5[0],pp5[1],pp5[2]];

                block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + mtri[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + mtri[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + mtri[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + ptri[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + ptri[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + ptri[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + mtri[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + mtri[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + mtri[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + ptri[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + ptri[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + ptri[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + mtri[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + mtri[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + mtri[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + mtri[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + mtri[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + mtri[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + ptri[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + ptri[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + ptri[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(c1c4Perp[0]); block_norm_tri.push(c1c4Perp[1]); block_norm_tri.push(c1c4Perp[2]);
                block_norm_tri.push(c1c4Perp[0]); block_norm_tri.push(c1c4Perp[1]); block_norm_tri.push(c1c4Perp[2]);
                block_norm_tri.push(c1c4Perp[0]); block_norm_tri.push(c1c4Perp[1]); block_norm_tri.push(c1c4Perp[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                
                block_vert_tri.push(centre[0] + mtri[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + mtri[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + mtri[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + ptri[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + ptri[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + ptri[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + ptri[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + ptri[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + ptri[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_norm_tri.push(-c1c4Perp[0]); block_norm_tri.push(-c1c4Perp[1]); block_norm_tri.push(-c1c4Perp[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + ptri[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + ptri[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + ptri[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + ptri[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + ptri[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + ptri[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-normp[0]); block_norm_tri.push(-normp[1]); block_norm_tri.push(-normp[2]);
                block_norm_tri.push(-normp[0]); block_norm_tri.push(-normp[1]); block_norm_tri.push(-normp[2]);
                block_norm_tri.push(-normp[0]); block_norm_tri.push(-normp[1]); block_norm_tri.push(-normp[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + ptri[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + ptri[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + ptri[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(normp[0]); block_norm_tri.push(normp[1]); block_norm_tri.push(normp[2]);
                block_norm_tri.push(normp[0]); block_norm_tri.push(normp[1]); block_norm_tri.push(normp[2]);
                block_norm_tri.push(normp[0]); block_norm_tri.push(normp[1]); block_norm_tri.push(normp[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + mtri[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + mtri[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + mtri[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + mtri[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + mtri[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + mtri[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(normm[0]); block_norm_tri.push(normm[1]); block_norm_tri.push(normm[2]);
                block_norm_tri.push(normm[0]); block_norm_tri.push(normm[1]); block_norm_tri.push(normm[2]);
                block_norm_tri.push(normm[0]); block_norm_tri.push(normm[1]); block_norm_tri.push(normm[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + mtri[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + mtri[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + mtri[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-normm[0]); block_norm_tri.push(-normm[1]); block_norm_tri.push(-normm[2]);
                block_norm_tri.push(-normm[0]); block_norm_tri.push(-normm[1]); block_norm_tri.push(-normm[2]);
                block_norm_tri.push(-normm[0]); block_norm_tri.push(-normm[1]); block_norm_tri.push(-normm[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                
            } else if(res.getName()==="BEM"||res.getName()==="GTR"||res.getName()==="ADA"||res.getName()==="DGU"||res.getName()==="KDN"||res.getName()==="SI3"||res.getName()==="NCC"||res.getName()==="IDR"||res.getName()==="GC4"||res.getName()==="GCD"||res.getName()==="GCU"||res.getName()==="GCV"||res.getName()==="GCW"||res.getName()==="IDS"||res.getName()==="REL"||res.getName()==="SIA"){
                /* Draw a diamond */
                // Try 4byh
                var horizontal = false;
                var vertical = false;
                var invert_colour = false;
                if(res.getName()==="IDR"||res.getName()==="GC4"||res.getName()==="GCD"||res.getName()==="GCU"||res.getName()==="GCV"||res.getName()==="GCW"||res.getName()==="IDS"||res.getName()==="REL"){
                    horizontal = true;
                }
                if(res.getName()==="GTR"||res.getName()==="ADA"||res.getName()==="DGU"||res.getName()==="BEM"){
                    vertical = true;
                }
                if(res.getName()==="IDR"||res.getName()==="BEM"){
                    invert_colour = true;
                }
                console.log("Draw a diamond!");
                var two_colour = two_colour_in && (horizontal||vertical);

                var c1c4 = vec3.create();
                vec3.subtract(cat1,cat4,c1c4);
                vec3.normalize(c1c4);

                var c1c4Perp = vec3.create();
                vec3.cross(normal,c1c4,c1c4Perp);
                vec3.normalize(c1c4Perp);

                var norm2 = vec3.create();
                vec3.cross(c1c4,c1c4Perp,norm2);

                c1c4[0] *= sugar_block_radius;
                c1c4[1] *= sugar_block_radius;
                c1c4[2] *= sugar_block_radius;
                c1c4Perp[0] *= sugar_block_radius;
                c1c4Perp[1] *= sugar_block_radius;
                c1c4Perp[2] *= sugar_block_radius;

                var ls = vec3.create();
                var le = vec3.create();
                vec3.add(centre,c1c4,ls);
                vec3.subtract(centre,c1c4Perp,le);

                var t = DistanceBetweenPointAndLine(ls,le,cat1)[1];
                var pp = vec3.create();
                var lels = vec3.create();
                vec3.subtract(le,ls,lels);
                lels[0] *= t;
                lels[1] *= t;
                lels[2] *= t;
                vec3.add(ls,lels,pp);
                res["GLYCO_BLOCK_C1"] = [pp[0],pp[1],pp[2]];
                
                var t2 = DistanceBetweenPointAndLine(ls,le,cat2)[1];
                var pp2 = vec3.create();
                var lels2 = vec3.create();
                vec3.subtract(le,ls,lels2);
                lels2[0] *= t2;
                lels2[1] *= t2;
                lels2[2] *= t2;
                vec3.add(ls,lels2,pp2);
                res["GLYCO_BLOCK_C2"] = [pp2[0],pp2[1],pp2[2]];

                vec3.subtract(centre,c1c4,ls);
                vec3.subtract(centre,c1c4Perp,le);

                var t3 = DistanceBetweenPointAndLine(ls,le,cat3)[1];
                var pp3 = vec3.create();
                var lels3 = vec3.create();
                vec3.subtract(le,ls,lels3);
                lels3[0] *= t3;
                lels3[1] *= t3;
                lels3[2] *= t3;
                vec3.add(ls,lels3,pp3);
                res["GLYCO_BLOCK_C3"] = [pp3[0],pp3[1],pp3[2]];

                var t4 = DistanceBetweenPointAndLine(ls,le,cat4)[1];
                var pp4 = vec3.create();
                var lels4 = vec3.create();
                vec3.subtract(le,ls,lels4);
                lels4[0] *= t4;
                lels4[1] *= t4;
                lels4[2] *= t4;
                vec3.add(ls,lels4,pp4);
                res["GLYCO_BLOCK_C4"] = [pp4[0],pp4[1],pp4[2]];
                
                vec3.subtract(centre,c1c4,ls);
                vec3.add(centre,c1c4Perp,le);
                
                var t5 = DistanceBetweenPointAndLine(ls,le,cat5)[1];
                var pp5 = vec3.create();
                var lels5 = vec3.create();
                vec3.subtract(le,ls,lels5);
                lels5[0] *= t5;
                lels5[1] *= t5;
                lels5[2] *= t5;
                vec3.add(ls,lels5,pp5);
                res["GLYCO_BLOCK_C5"] = [pp5[0],pp5[1],pp5[2]];

                console.log("horizontal");
                console.log(horizontal);
                if(horizontal){
                    block_vert_tri.push(centre[0] - c1c4[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] + c1c4[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    if(two_colour&&invert_colour) {
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    } else {
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    }

                    block_vert_tri.push(centre[0] + c1c4[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] - c1c4[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    if(two_colour&&invert_colour) {
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    } else {
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    }
                    
                    block_vert_tri.push(centre[0] - c1c4[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] + c1c4[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    if(two_colour&&invert_colour) {
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    } else {
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    }

                    block_vert_tri.push(centre[0] + c1c4[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] - c1c4[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    if(two_colour&&invert_colour) {
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    } else {
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    }
                } else {
                    block_vert_tri.push(centre[0] - c1c4[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    if(two_colour&&!invert_colour) {
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    } else {
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    }

                    block_vert_tri.push(centre[0] + c1c4[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    if(two_colour&&!invert_colour) {
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    } else {
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    }

                    block_vert_tri.push(centre[0] - c1c4[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    if(two_colour&&!invert_colour) {
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    } else {
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    }

                    block_vert_tri.push(centre[0] + c1c4[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                    block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                    block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                    block_idx_tri.push(block_vert_tri.length/3-1);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                    if(two_colour&&!invert_colour) {
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                        block_col_tri.push(white[0]); block_col_tri.push(white[1]); block_col_tri.push(white[2]); block_col_tri.push(white[3]);
                    } else {
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                        block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                    }

                }

                block_vert_tri.push(centre[0] + c1c4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(c1c4Perp[0]-c1c4[0]); block_norm_tri.push(c1c4Perp[1]-c1c4[1]); block_norm_tri.push(c1c4Perp[2]-c1c4[2]);
                block_norm_tri.push(c1c4Perp[0]-c1c4[0]); block_norm_tri.push(c1c4Perp[1]-c1c4[1]); block_norm_tri.push(c1c4Perp[2]-c1c4[2]);
                block_norm_tri.push(c1c4Perp[0]-c1c4[0]); block_norm_tri.push(c1c4Perp[1]-c1c4[1]); block_norm_tri.push(c1c4Perp[2]-c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4Perp[0]+c1c4[0]); block_norm_tri.push(-c1c4Perp[1]+c1c4[1]); block_norm_tri.push(-c1c4Perp[2]+c1c4[2]);
                block_norm_tri.push(-c1c4Perp[0]+c1c4[0]); block_norm_tri.push(-c1c4Perp[1]+c1c4[1]); block_norm_tri.push(-c1c4Perp[2]+c1c4[2]);
                block_norm_tri.push(-c1c4Perp[0]+c1c4[0]); block_norm_tri.push(-c1c4Perp[1]+c1c4[1]); block_norm_tri.push(-c1c4Perp[2]+c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] - c1c4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4Perp[0]+c1c4[0]); block_norm_tri.push(-c1c4Perp[1]+c1c4[1]); block_norm_tri.push(-c1c4Perp[2]+c1c4[2]);
                block_norm_tri.push(-c1c4Perp[0]+c1c4[0]); block_norm_tri.push(-c1c4Perp[1]+c1c4[1]); block_norm_tri.push(-c1c4Perp[2]+c1c4[2]);
                block_norm_tri.push(-c1c4Perp[0]+c1c4[0]); block_norm_tri.push(-c1c4Perp[1]+c1c4[1]); block_norm_tri.push(-c1c4Perp[2]+c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] - c1c4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(c1c4Perp[0]-c1c4[0]); block_norm_tri.push(c1c4Perp[1]-c1c4[1]); block_norm_tri.push(c1c4Perp[2]-c1c4[2]);
                block_norm_tri.push(c1c4Perp[0]-c1c4[0]); block_norm_tri.push(c1c4Perp[1]-c1c4[1]); block_norm_tri.push(c1c4Perp[2]-c1c4[2]);
                block_norm_tri.push(c1c4Perp[0]-c1c4[0]); block_norm_tri.push(c1c4Perp[1]-c1c4[1]); block_norm_tri.push(c1c4Perp[2]-c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] - c1c4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4Perp[0]-c1c4[0]); block_norm_tri.push(-c1c4Perp[1]-c1c4[1]); block_norm_tri.push(-c1c4Perp[2]-c1c4[2]);
                block_norm_tri.push(-c1c4Perp[0]-c1c4[0]); block_norm_tri.push(-c1c4Perp[1]-c1c4[1]); block_norm_tri.push(-c1c4Perp[2]-c1c4[2]);
                block_norm_tri.push(-c1c4Perp[0]-c1c4[0]); block_norm_tri.push(-c1c4Perp[1]-c1c4[1]); block_norm_tri.push(-c1c4Perp[2]-c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] - c1c4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(c1c4Perp[0]+c1c4[0]); block_norm_tri.push(c1c4Perp[1]+c1c4[1]); block_norm_tri.push(c1c4Perp[2]+c1c4[2]);
                block_norm_tri.push(c1c4Perp[0]+c1c4[0]); block_norm_tri.push(c1c4Perp[1]+c1c4[1]); block_norm_tri.push(c1c4Perp[2]+c1c4[2]);
                block_norm_tri.push(c1c4Perp[0]+c1c4[0]); block_norm_tri.push(c1c4Perp[1]+c1c4[1]); block_norm_tri.push(c1c4Perp[2]+c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(c1c4Perp[0]+c1c4[0]); block_norm_tri.push(c1c4Perp[1]+c1c4[1]); block_norm_tri.push(c1c4Perp[2]+c1c4[2]);
                block_norm_tri.push(c1c4Perp[0]+c1c4[0]); block_norm_tri.push(c1c4Perp[1]+c1c4[1]); block_norm_tri.push(c1c4Perp[2]+c1c4[2]);
                block_norm_tri.push(c1c4Perp[0]+c1c4[0]); block_norm_tri.push(c1c4Perp[1]+c1c4[1]); block_norm_tri.push(c1c4Perp[2]+c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + c1c4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4Perp[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + c1c4Perp[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + c1c4Perp[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + c1c4Perp[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-c1c4Perp[0]-c1c4[0]); block_norm_tri.push(-c1c4Perp[1]-c1c4[1]); block_norm_tri.push(-c1c4Perp[2]-c1c4[2]);
                block_norm_tri.push(-c1c4Perp[0]-c1c4[0]); block_norm_tri.push(-c1c4Perp[1]-c1c4[1]); block_norm_tri.push(-c1c4Perp[2]-c1c4[2]);
                block_norm_tri.push(-c1c4Perp[0]-c1c4[0]); block_norm_tri.push(-c1c4Perp[1]-c1c4[1]); block_norm_tri.push(-c1c4Perp[2]-c1c4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                
            } else if(res.getName()==="XLS"||res.getName()==="CXY"||res.getName()==="RBY"||res.getName()==="TDX"||res.getName()==="XYL"||res.getName()==="XYS"||res.getName()==="XYP"){
                /* Draw a star! */
                console.log("Draw a star!");
                // Try 4cuo

                var c1c4 = vec3.create();
                vec3.subtract(cat1,cat4,c1c4);
                vec3.normalize(c1c4);

                var c1c4Perp = vec3.create();
                vec3.cross(normal,c1c4,c1c4Perp);
                vec3.normalize(c1c4Perp);

                var norm2 = vec3.create();
                vec3.cross(c1c4,c1c4Perp,norm2);

                c1c4[0] *= sugar_block_radius;
                c1c4[1] *= sugar_block_radius;
                c1c4[2] *= sugar_block_radius;
                c1c4Perp[0] *= sugar_block_radius;
                c1c4Perp[1] *= sugar_block_radius;
                c1c4Perp[2] *= sugar_block_radius;

                // http://mathworld.wolfram.com/Pentagon.html
                var c1 = Math.cos(2.0*Math.PI/5.);
                var c2 = Math.cos(Math.PI/5.);
                var s1 = Math.sin(2.0*Math.PI/5.);
                var s2 = Math.sin(4.0*Math.PI/5.);

                var p1 = vec3.create([c1c4Perp[0],c1c4Perp[1],c1c4Perp[2]]);
                var p2 = vec3.create([s1 * c1c4[0] + c1 * c1c4Perp[0], s1 * c1c4[1] + c1 * c1c4Perp[1], s1 * c1c4[2] + c1 * c1c4Perp[2]]);
                var p3 = vec3.create([s2 * c1c4[0] - c2 * c1c4Perp[0], s2 * c1c4[1] - c2 * c1c4Perp[1], s2 * c1c4[2] - c2 * c1c4Perp[2]]);
                var p4 = vec3.create([-s2 * c1c4[0] - c2 * c1c4Perp[0], -s2 * c1c4[1] - c2 * c1c4Perp[1], -s2 * c1c4[2] - c2 * c1c4Perp[2]]);
                var p5 = vec3.create([-s1 * c1c4[0] + c1 * c1c4Perp[0], -s1 * c1c4[1] + c1 * c1c4Perp[1], -s1 * c1c4[2] + c1 * c1c4Perp[2]]);

                var t6 = DistanceBetweenTwoLines(p1,p3,p5,p2)[1];
                var t7 = DistanceBetweenTwoLines(p1,p3,p2,p4)[1];
                var t8 = DistanceBetweenTwoLines(p2,p4,p3,p5)[1];
                var t9 = DistanceBetweenTwoLines(p4,p1,p3,p5)[1];
                var t10 = DistanceBetweenTwoLines(p5,p2,p4,p1)[1];

                function AddFractionToLine(pp1,pp2,frac){
                    var ret = vec3.create();
                    var diff = vec3.create();
                    vec3.subtract(pp2,pp1,diff);
                    diff[0] *= frac;
                    diff[1] *= frac;
                    diff[2] *= frac;
                    vec3.add(pp1,diff,ret);
                    return ret;
                }

                var p6  = AddFractionToLine(p1,p3,t6);
                var p7  = AddFractionToLine(p1,p3,t7);
                var p8  = AddFractionToLine(p2,p4,t8);
                var p9  = AddFractionToLine(p4,p1,t9);
                var p10 = AddFractionToLine(p5,p2,t10);
                res["GLYCO_BLOCK_C1"] = [centre[0]+p2[0],centre[1]+p2[1],centre[2]+p2[2]];
                res["GLYCO_BLOCK_C2"] = [centre[0]+p10[0],centre[1]+p10[1],centre[2]+p10[2]];
                res["GLYCO_BLOCK_C3"] = [centre[0]+p6[0],centre[1]+p6[1],centre[2]+p6[2]];
                res["GLYCO_BLOCK_C4"] = [centre[0]+p5[0],centre[1]+p5[1],centre[2]+p5[2]];
                res["GLYCO_BLOCK_C5"] = [centre[0]+p3[0],centre[1]+p3[1],centre[2]+p3[2]];

                block_vert_tri.push(centre[0] + p6[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p7[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p8[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p9[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p10[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p1[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p1[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p1[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p2[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p2[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p2[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p3[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p3[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p3[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p5[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p5[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p5[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p6[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p7[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p8[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p9[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p10[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_norm_tri.push(-norm2[0]); block_norm_tri.push(-norm2[1]); block_norm_tri.push(-norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p1[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p1[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p1[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p2[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p2[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p2[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p3[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p3[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p3[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p4[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p4[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p4[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p5[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p5[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p5[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_norm_tri.push(norm2[0]); block_norm_tri.push(norm2[1]); block_norm_tri.push(norm2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p1[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p1[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p1[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p1[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p1[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p1[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p1[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p1[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p1[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p6[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p6[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p2[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p2[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p2[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p6[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p6[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p6[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p2[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p2[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p2[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p2[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p2[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p2[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p2[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p2[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p2[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p2[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p2[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p2[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p2[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p2[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p2[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p7[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p7[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p3[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p3[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p3[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p7[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p7[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p7[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p3[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p3[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p3[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p3[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p3[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p3[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_norm_tri.push(p2[0]); block_norm_tri.push(p2[1]); block_norm_tri.push(p2[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p3[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p3[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p3[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p3[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p3[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p3[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p3[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p3[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p3[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p8[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p8[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p8[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p8[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p8[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p4[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p4[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p4[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_norm_tri.push(p3[0]); block_norm_tri.push(p3[1]); block_norm_tri.push(p3[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p4[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p4[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p4[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p4[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p4[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p4[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p4[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p4[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p4[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p9[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p9[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p5[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p5[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p5[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p9[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p9[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p9[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p5[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p5[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p5[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p5[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p5[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p5[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_norm_tri.push(p4[0]); block_norm_tri.push(p4[1]); block_norm_tri.push(p4[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p5[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p5[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p5[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p5[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p5[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p5[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p5[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p5[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p5[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_norm_tri.push(p1[0]); block_norm_tri.push(p1[1]); block_norm_tri.push(p1[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p10[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p10[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p1[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p1[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p1[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

                block_vert_tri.push(centre[0] + p10[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p10[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p10[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p1[0] + sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p1[1] + sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p1[2] + sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_vert_tri.push(centre[0] + p1[0] - sugar_block_thickness*norm2[0]);
                block_vert_tri.push(centre[1] + p1[1] - sugar_block_thickness*norm2[1]);
                block_vert_tri.push(centre[2] + p1[2] - sugar_block_thickness*norm2[2]);
                block_idx_tri.push(block_vert_tri.length/3-1);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_norm_tri.push(p5[0]); block_norm_tri.push(p5[1]); block_norm_tri.push(p5[2]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);
                block_col_tri.push(col1[0]); block_col_tri.push(col1[1]); block_col_tri.push(col1[2]); block_col_tri.push(col1[3]);

            } else {
                console.log("I don't know what to draw for "+res.getName());
            }
        }
    }

}

function DrawSugarBlock(res1,sugar_block_thickness,sugarAtomNames, block_vert_tri,block_norm_tri,block_col_tri,block_idx_tri,col1){
    //console.log("Consider glyco block "+res1.getName());
    for(var isugartype=0;isugartype<sugarAtomNames.length;isugartype++){
        for(var iring=0;iring<sugarAtomNames[isugartype].length;iring++){
            if(sugarAtomNames[isugartype][iring].length>4){
                var c1name = sugarAtomNames[isugartype][iring][0];
                var c2name = sugarAtomNames[isugartype][iring][1];
                var c3name = sugarAtomNames[isugartype][iring][2];
                var c4name = sugarAtomNames[isugartype][iring][3];
                var c5name = sugarAtomNames[isugartype][iring][4];
                var o5name = null;
                var c1 = res1.getAtom(c1name);
                var c2 = res1.getAtom(c2name);
                var c3 = res1.getAtom(c3name);
                var c4 = res1.getAtom(c4name);
                var c5 = res1.getAtom(c5name);
                var o5 = null;
                if(sugarAtomNames[isugartype][iring].length==6){
                    o5name = sugarAtomNames[isugartype][iring][5];
                    o5 = res1.getAtom(o5name);
                    //std::cout << "Looking for " << c1name << " " << c2name << " " << c3name << " " << c4name << " " << c5name << " " << o5name << "\n";
                } else {
                    //std::cout << "Looking for " << c1name << " " << c2name << " " << c3name << " " << c4name << " " << c5name << "\n";
                }
                if(c1&&c2&&c3&&c4&&c5&&(o5||sugarAtomNames[isugartype][iring].length!=6)){
                    //console.log("Have them all");
                    // FIXME - col1 is colour of this atom. two_colour depends on whether colour by restype or not.
                    var sugar_block_scale = 1.0;
                    var two_colour = true;
                    DrawSugarBlockInt(res1, col1, c1,c2,c3,c4,c5,o5, two_colour, sugar_block_thickness, sugar_block_scale, block_vert_tri,block_norm_tri,block_col_tri,block_idx_tri );
                    //std::vector<Cartesian> res = DrawSugarBlockInt(tris, cyls, res1, col1, params, global_params, selHnd, c1,c2,c3,c4,c5,o5,res1, molHnd, two_colour, sugar_block_thickness, sugar_block_scale );
                    //if(iring==sugarAtomNames[isugartype].length-1&&res.length==2) return res;
                    return;
                }

            }
        }
    }
    //std::vector<Cartesian> retval;
    //return retval;

}

function residuesToGlycoBlocksInfo(glycanResidues,size,colourScheme){
    var block_col_tri = [];
    var block_norm_tri = [];
    var block_vert_tri = [];
    var block_idx_tri = [];
    var glycanAtoms = [];

    var sugarAtomNames = getSugarAtomNames();

    for(var ig=0;ig<glycanResidues.length;ig++){
        for(var iat=0;iat<glycanResidues[ig].atoms.length;iat++){
            glycanAtoms.push(glycanResidues[ig].atoms[iat]);
        }
    }

    for(var ig=0;ig<glycanResidues.length;ig++){
        var col1 = colourScheme[glycanResidues[ig].atoms[0]["_atom_site.id"]];
        var sugarBlockPrims = DrawSugarBlock(glycanResidues[ig],size,sugarAtomNames,block_vert_tri,block_norm_tri,block_col_tri,block_idx_tri,col1);
    }

    var blockPrimitiveInfo = {"atoms":[[glycanAtoms]], "col_tri":[[block_col_tri]], "norm_tri":[[block_norm_tri]], "vert_tri":[[block_vert_tri]], "idx_tri":[[block_idx_tri]] , "prim_types":[["TRIANGLES"]] };
    return blockPrimitiveInfo;
}

function atomsToSpheresInfo(atoms,size,colourScheme,scalebyVDWRadii){
    // FIXME - PERFECT_SPHERES, POINTS_SPHERES option?
    return atomsToCirclesSpheresInfo(atoms,size,"POINTS_SPHERES",colourScheme,scalebyVDWRadii);
}

function atomsToCirclesInfo(atoms,size,colourScheme,scalebyVDWRadii){
    return atomsToCirclesSpheresInfo(atoms,size,"POINTS",colourScheme,scalebyVDWRadii);
}

function ColourScheme(pdbatoms){
    this.pdbatoms = pdbatoms;
    this.hier = this.pdbatoms["atoms"];
    this.ice_blue     = [0.61,0.69,1.00,1.0];
    this.gold         = [0.70,0.69,0.24,1.0];
    this.coral        = [1.00,0.50,0.31,1.0];
    this.grey         = [0.50,0.50,0.50,1.0];
    this.pink         = [1.00,0.57,1.00,1.0];
    this.sea_green    = [0.50,0.73,0.71,1.0];
    this.pale_brown   = [0.66,0.49,0.37,1.0];
    this.lilac        = [0.68,0.53,0.73,1.0];
    this.lemon        = [1.00,1.00,0.50,1.0];
    this.lawn_green   = [0.27,0.61,0.31,1.0];
    this.pale_crimson = [0.82,0.24,0.24,1.0];
    this.light_blue   = [0.25,0.60,0.88,1.0];
    this.tan          = [0.47,0.00,0.00,1.0];
    this.light_green  = [0.60,1.00,0.60,1.0];
    this.yellow       = [1.00,1.00,0.00,1.0];
    this.white        = [1.00,1.00,1.00,1.0];
    this.blue         = [0.00,0.00,1.00,1.0];
    this.red          = [1.00,0.00,0.00,1.0];
    this.green        = [0.00,1.00,0.00,1.0];
    this.magenta      = [1.00,0.00,1.00,1.0];
    this.cyan         = [0.00,1.00,0.88,1.0];
    this.purple       = [0.58,0.00,1.00,1.0];
    this.dark_purple  = [0.57,0.13,0.34,1.0];
    this.dark_cyan    = [0.10,0.60,0.70,1.0];
    this.black        = [0.00,0.00,0.00,1.0];
    this.residueColours = {"PHE":this.magenta,"TRP":this.magenta,"TYR":this.magenta,"PRO":this.coral,"VAL":this.coral,"ALA":this.coral,"ILE":this.coral,"LEU":this.coral,"SER":this.cyan,"THR":this.cyan,"ASN":this.cyan,"GLN":this.cyan,"ARG":this.blue,"LYS":this.blue,"ASP":this.red,"GLU":this.red,"CYS":this.yellow,"MET":this.yellow,"GLY":this.white,"HIS":this.light_blue,"A":this.red,"T":this.yellow,"G":this.green,"C":this.blue,"U":this.magenta,"DA":this.red,"DT":this.yellow,"DG":this.green,"DC":this.blue,"ADE":this.red,"THY":this.yellow,"GUA":this.green,"CYT":this.blue,"URA":this.magenta,"BMA":this.green,"MAN":this.green,"NAG":this.blue,"GLC":this.blue,"BGC":this.blue,"GCS":this.blue,"GAL":this.yellow,"NGA":this.yellow,"MGC":this.yellow,"NG1":this.yellow,"NG6":this.yellow,"A2G":this.yellow,"6MN":this.blue,"GLP":this.blue,"GP4":this.blue,"BEM":this.green,"KDN":this.light_green,"XLS":this.coral,"CXY":this.coral,"RBY":this.coral,"TDX":this.coral,"XYL":this.coral,"XYS":this.coral,"XYP":this.coral,"FCA":this.red,"FCB":this.red,"FUC":this.red,"GTR":this.yellow,"ADA":this.yellow,"DGU":this.yellow,"SI3":this.purple,"NCC":this.purple,"NGF":this.light_blue,"NGE":this.light_blue,"NGC":this.light_blue,"GC4":this.blue,"GCD":this.blue,"GCU":this.blue,"GCV":this.blue,"GCW":this.blue,"IDS":this.blue,"REL":this.blue,"IDR":this.pale_brown,"SIA":this.purple};
    this.order_colours = [this.ice_blue,this.gold,this.coral,this.grey,this.pink,this.sea_green,this.pale_brown,this.lilac,this.lemon,this.lawn_green,this.pale_crimson,this.light_blue,this.tan,this.light_green,this.yellow,this.blue,this.red,this.green,this.magenta,this.cyan,this.purple,this.dark_purple,this.dark_cyan];
}

ColourScheme.prototype.colourByAtomType = function (params){
    //FIXME - Need to consider models.
    var colarray = {};
    var model = this.hier[0];
    var atoms = model.getAllAtoms();
    for(var iat=0;iat<atoms.length;iat++){
        var element = atoms[iat].element();
        var serNum = atoms[iat]["_atom_site.id"];
        if(typeof(params)!=="undefined"&&element in params){
            colarray[serNum] = params[element];
        } else if(element==="O"){
            colarray[serNum] = [1.0,0.0,0.0,1.0];
        }else if(element==="C"){
            colarray[serNum] = [0.0,1.0,0.0,1.0];
        }else if(element==="N"){
            colarray[serNum] = [0.0,0.0,1.0,1.0];
        }else if(element==="S"){
            colarray[serNum] = [1.0,1.0,0.0,1.0];
        }else{
            colarray[serNum] = [0.5,0.5,0.5,1.0];
        }
    }
    return colarray;
}

ColourScheme.prototype.colourBySecondaryStructure = function (params){
    var colarray = {};
    var chains = this.hier[0].chains;
    for(var ic=0;ic<chains.length;ic++){
        var residues = chains[ic].residues;
        for(var ir=0;ir<residues.length;ir++){
            var ca = residues[ir].getAtomTrimmed( "CA" );
            if(ca&&typeof(ca.residue["SSE"])!=="undefined"){
                if(ca.residue["SSE"]==="SSE_Strand"||ca.residue["SSE"]==="SSE_Bulge"){
                    for(var iat=0;iat<residues[ir].atoms.length;iat++){
                        var serNum = residues[ir].atoms[iat]["_atom_site.id"];
                        colarray[serNum] = params["strand"];
                    }
                } else if(ca.residue["SSE"]==="SSE_Helix"){
                    for(var iat=0;iat<residues[ir].atoms.length;iat++){
                        var serNum = residues[ir].atoms[iat]["_atom_site.id"];
                        colarray[serNum] = params["helix"];
                    }
                } else {
                    for(var iat=0;iat<residues[ir].atoms.length;iat++){
                        var serNum = residues[ir].atoms[iat]["_atom_site.id"];
                        colarray[serNum] = [0.5,0.5,0.5,1.0];
                    }
                }
            } else {
                for(var iat=0;iat<residues[ir].atoms.length;iat++){
                    var serNum = residues[ir].atoms[iat]["_atom_site.id"];
                    colarray[serNum] = [0.5,0.5,0.5,1.0];
                }
            }
        }
    }
    return colarray;
}


ColourScheme.prototype.colourByResidueType = function (params){
    //FIXME - Need to consider models.
    var colarray = {};
    var model = this.hier[0];
    var atoms = model.getAllAtoms();
    for(var iat=0;iat<atoms.length;iat++){
        var resName = atoms[iat].residue.getName();
        var serNum = atoms[iat]["_atom_site.id"];
        if(typeof(params)!=="undefined"&&resName in params){
            colarray[serNum] = params[resName];
        } else if(resName in this.residueColours){
            colarray[serNum] = this.residueColours[resName];
        }else{
            colarray[serNum] = [0.5,0.5,0.5,1.0];
        }
    }
    return colarray;
}

ColourScheme.prototype.colourByEntity = function (params){
}

ColourScheme.prototype.colourByModel = function (params){
    var colarray = {};
    for(var im=0;im<this.hier.length;im++){
        var chains = this.hier[im].chains;
        for(var ic=0;ic<chains.length;ic++){
            var residues = chains[ic].residues;
            for(var ir=0;ir<residues.length;ir++){
                for(var iat=0;iat<residues[ir].atoms.length;iat++){
                    var serNum = residues[ir].atoms[iat]["_atom_site.id"];
                    colarray[serNum] = this.order_colours[im%this.order_colours.length];
                    if(typeof(params)!=="undefined"&&typeof(params["nonCByAtomType"]!=="undefined")&&params["nonCByAtomType"]){
                        var element = residues[ir].atoms[iat].element();
                        if(element in params){
                            colarray[serNum] = params[element];
                        } else if(element==="O"){
                            colarray[serNum] = [1.0,0.0,0.0,1.0];
                        }else if(element==="N"){
                            colarray[serNum] = [0.0,0.0,1.0,1.0];
                        }else if(element==="S"){
                            colarray[serNum] = [1.0,1.0,0.0,1.0];
                        }
                    }
                }
            }
        }
    }
    return colarray;
}

ColourScheme.prototype.colourByChain = function (params){
    var colarray = {};
    for(var im=0;im<this.hier.length;im++){
        var chains = this.hier[im].chains;
        for(var ic=0;ic<chains.length;ic++){
            var residues = chains[ic].residues;
            for(var ir=0;ir<residues.length;ir++){
                for(var iat=0;iat<residues[ir].atoms.length;iat++){
                    var serNum = residues[ir].atoms[iat]["_atom_site.id"];
                    colarray[serNum] = this.order_colours[ic%this.order_colours.length];
                    if(typeof(params)!=="undefined"&&typeof(params["nonCByAtomType"]!=="undefined")&&params["nonCByAtomType"]){
                        var element = residues[ir].atoms[iat].element();
                        if(element in params){
                            colarray[serNum] = params[element];
                        } else if(element==="O"){
                            colarray[serNum] = [1.0,0.0,0.0,1.0];
                        }else if(element==="N"){
                            colarray[serNum] = [0.0,0.0,1.0,1.0];
                        }else if(element==="S"){
                            colarray[serNum] = [1.0,1.0,0.0,1.0];
                        }
                    }
                }
            }
        }
    }
    return colarray;
}

ColourScheme.prototype.colourByMainSide = function (params){
    var colarray = {};
    for(var im=0;im<this.hier.length;im++){
        var allAtoms = this.hier[im].getAllAtoms();
        var sideAtoms = this.hier[im].getAtoms("side");
        var mainAtoms = this.hier[im].getAtoms("main");
        for(var iat=0;iat<allAtoms.length;iat++){
            var serNum = allAtoms[iat]["_atom_site.id"];
            if(mainAtoms.indexOf(allAtoms[iat])!==-1){
                colarray[serNum] = this.lemon;
            } else if(sideAtoms.indexOf(allAtoms[iat])!==-1){
                colarray[serNum] = this.ice_blue;
            } else {
                colarray[serNum] = [0.5,0.5,0.5,1.0];
            }
        }
    }
    return colarray;
}

ColourScheme.prototype.colourByProperty = function (params,prop){
    //FIXME - Need to consider models.
    //FIXME - Below bottom and above top ranges need to be done.
    var colarray = {};
    var model = this.hier[0];
    var atoms = model.getAllAtoms();

    for(var ip=0;ip<params.length;ip++){
        for(var iat=0;iat<atoms.length;iat++){
            var Prop = atoms[iat][prop];
            var serNum = atoms[iat]["_atom_site.id"];
            if(isNaN(Prop)){
                    colarray[serNum] = [0.5,0.5,0.5,1.0];
            } else {
                if(Prop<params[0].min){ // FIXME see above
                    colarray[serNum] = params[0].mincolour;
                } else if(Prop>params[params.length-1].max){ // FIXME see above
                    colarray[serNum] = params[params.length-1].maxcolour;
                } else if(Prop>=params[ip].min&&Prop<=params[ip].max){
                    var frac = 1.0*(Prop - params[ip].min) / (params[ip].max- params[ip].min);
                    var r = frac*params[ip].maxcolour[0] + (1.0-frac)*params[ip].mincolour[0];
                    var g = frac*params[ip].maxcolour[1] + (1.0-frac)*params[ip].mincolour[1];
                    var b = frac*params[ip].maxcolour[2] + (1.0-frac)*params[ip].mincolour[2];
                    var a = frac*params[ip].maxcolour[3] + (1.0-frac)*params[ip].mincolour[3];
                    colarray[serNum] = [r,g,b,a];
                }
            }
        }
    }

    return colarray;
}

ColourScheme.prototype.colourByBFactor = function (params){
    return this.colourByProperty(params,"_atom_site.B_iso_or_equiv");
}

ColourScheme.prototype.colourByOccupancy = function (params){
    return this.colourByProperty(params,"_atom_site.occupancy");
}

ColourScheme.prototype.colourByCharge = function (params){
    return this.colourByProperty(params,"_atom_site.pdbx_formal_charge");
}

function hsvtorgb (hsv) {
    var r,g,b;
    var rgb = [];


    rgb[3] = hsv[3];


    if ( Math.abs(hsv[0]) < 1e-6 && Math.abs(hsv[1]) < 1e-6 && Math.abs(hsv[2]) < 1e-6 ) {
        rgb[0] = 1.0; 
        rgb[1] = 1.0; 
        rgb[2] = 1.0; 
        return rgb;
    }
    if ( Math.abs(hsv[1]) < 1e-6 ) {
        rgb[0] = rgb[1] = rgb[2] = hsv[2];
        return rgb;
    }

    var h = hsv[0] / 60.0;
    var i = parseInt(Math.floor(hsv[0]/60.0)); // Need to round up
    var f = h - i;
    var p = hsv[2] * (1.0 - hsv[1]);
    var q = hsv[2] * (1.0 - hsv[1] * f);
    var t = hsv[2] * (1.0 - hsv[1] * (1.0 - f));
    var v = hsv[2];

    switch (i) {
        case 0 :
            r = v;
            g  = t;
            b  = p;
            break;
        case  1 : 
            r  = q;
            g  = v;
            b  = p;
            break;
        case 2:
            r = p;
            g = v;
            b = t;
            break;
        case 3 :
            r  = p;
            g = q;
            b = v;
            break;
        case 4 :
            r = t;
            g = p;
            b = v;
            break;
        default :
            r = v;
            g = p;
            b = q;
            break;
    }
    rgb[0] = r;
    rgb[1] = g;
    rgb[2] = b;

    return rgb;
}

function rgbtohsv (rgb) {
    var s,h,v;
    var hsv = [];

    hsv[3] = rgb[3];

    var maxrgb = Math.max(rgb[0],rgb[1]);
    maxrgb = Math.max(maxrgb,rgb[2]);
    var minrgb = Math.min(rgb[0],rgb[1]);
    minrgb = Math.min(minrgb,rgb[2]);

    v = maxrgb;
    var delta = maxrgb - minrgb;

    if (Math.abs(delta)<1e-6) {
        hsv[0] = 0;
        hsv[1] = 0;
        hsv[2] = v;
        return hsv;
    }

    if ( Math.abs(maxrgb) > 1e-6 ) {
        s = delta / maxrgb;
    } else {
        hsv[0] = 0;
        hsv[1] = -1;
        hsv[2] = v;
        return hsv;
    }

    if ( Math.abs(rgb[0]-maxrgb) < 1e-6 ) {
        h = 1.0*(rgb[1] - rgb[2]) /delta;
    } else if (Math.abs(rgb[1] - maxrgb)< 1e-6 ) {
        h= 2 + (1.0*( rgb[2] - rgb[0] )/delta);
    } else {
        h = 4 + (1.0*( rgb[0] - rgb[1] )  / delta);
    }

    h= h * 60;
    if (h < 0 ) h = h + 360;

    hsv[0] = h;
    hsv[1] = s;
    hsv[2] = v;
    return hsv;

}

function colourRamp(start_colour_in,end_colour_in,frac,mode){
    if(mode=="hsv"){
        /* 
           Now this ought to be more complicated: wheel direction is important. But this simple way works for rainbow.
           We can do the wheel the other way by splitting interpolation into two ranges, one which goes to hue of 359 and other from 0.
           Or other way round depending in which initial colour has highest hue. We'll do this in due course.
         */
        var start_colour = rgbtohsv([start_colour_in[0],start_colour_in[1],start_colour_in[2],start_colour_in[3]]);
        var end_colour = rgbtohsv([end_colour_in[0],end_colour_in[1],end_colour_in[2],end_colour_in[3]]);

        var r = frac*end_colour[0] + (1.0-frac)*start_colour[0];
        var g = frac*end_colour[1] + (1.0-frac)*start_colour[1];
        var b = frac*end_colour[2] + (1.0-frac)*start_colour[2];
        var a = frac*end_colour[3] + (1.0-frac)*start_colour[3];
        return hsvtorgb([r,g,b,a]);
    } else {
        var start_colour = [start_colour_in[0],start_colour_in[1],start_colour_in[2],start_colour_in[3]];
        var end_colour = [end_colour_in[0],end_colour_in[1],end_colour_in[2],end_colour_in[3]];

        var r = frac*end_colour[0] + (1.0-frac)*start_colour[0];
        var g = frac*end_colour[1] + (1.0-frac)*start_colour[1];
        var b = frac*end_colour[2] + (1.0-frac)*start_colour[2];
        var a = frac*end_colour[3] + (1.0-frac)*start_colour[3];
        return [r,g,b,a];
    }
}

ColourScheme.prototype.colourRainbow = function (params){
    var colarray = {};
    var start_colour = [0.0,0.0,1.0,1.0];
    var end_colour   = [1.0,0.0,0.0,1.0];
    for(var im=0;im<this.hier.length;im++){
        var chains = this.hier[im].chains;
        for(var ic=0;ic<chains.length;ic++){
            var residues = chains[ic].residues;
            var nPeptide = 0;
            for(var ir=0;ir<residues.length;ir++){
                if(residues[ir].isAminoAcid()){
                    nPeptide++;
                }
            }
            var iPeptide = 0;
            for(var ir=0;ir<residues.length;ir++){
                if(residues[ir].isAminoAcid()){
                    var frac = 1.0 * iPeptide / nPeptide;
                    //frac = 0.9;
                    var theColour = colourRamp(start_colour,end_colour,frac,"hsv");
                    for(var iat=0;iat<residues[ir].atoms.length;iat++){
                        var serNum = residues[ir].atoms[iat]["_atom_site.id"];
                        colarray[serNum] = theColour;
                    }
                    iPeptide++;
                } else {
                    for(var iat=0;iat<residues[ir].atoms.length;iat++){
                        var serNum = residues[ir].atoms[iat]["_atom_site.id"];
                        colarray[serNum] = [0.5,0.5,0.5,1.0];
                    }
                }
            }
        }
    }
    return colarray;
}

ColourScheme.prototype.colourOneColour = function (col,params){
    //FIXME - Need to consider models.
    var colarray = {};
    var model = this.hier[0];
    var atoms = model.getAllAtoms();

    if(typeof(params)!=="undefined"&&typeof(params["nonCByAtomType"]!=="undefined")&&params["nonCByAtomType"]){
        for(var im=0;im<this.hier.length;im++){
            var chains = this.hier[im].chains;
            for(var ic=0;ic<chains.length;ic++){
                var residues = chains[ic].residues;
                for(var ir=0;ir<residues.length;ir++){
                    for(var iat=0;iat<residues[ir].atoms.length;iat++){
                        var serNum = residues[ir].atoms[iat]["_atom_site.id"];
                        colarray[serNum] = [col[0],col[1],col[2],col[3]];
                        var element = residues[ir].atoms[iat].element();
                        if(element in params){
                            colarray[serNum] = params[element];
                        } else if(element==="O"){
                            colarray[serNum] = [1.0,0.0,0.0,1.0];
                        }else if(element==="N"){
                            colarray[serNum] = [0.0,0.0,1.0,1.0];
                        }else if(element==="S"){
                            colarray[serNum] = [1.0,1.0,0.0,1.0];
                        }
                    }
                }
            }
        }
    } else {
        for(var iat=0;iat<atoms.length;iat++){
            var element = atoms[iat].element();
            var serNum = atoms[iat]["_atom_site.id"];
            colarray[serNum] = [col[0],col[1],col[2],col[3]];
        }
    }

    return colarray;
}

ColourScheme.prototype.colourByHydrophobicity = function (params){
}

ColourScheme.prototype.colourByAtomSAS = function (params){
}

ColourScheme.prototype.colourByResideSAS = function (params){
}

ColourScheme.prototype.colourBySequenceConservation = function (params){
}

function atomsToCirclesSpheresInfo(atoms,size,primType,colourScheme,scalebyVDWRadii){
    var sphere_sizes = [];
    var sphere_col_tri = [];
    var sphere_vert_tri = [];
    var sphere_idx_tri = [];
    var sphere_atoms = [];
    for(var iat=0;iat<atoms.length;iat++){
        sphere_idx_tri.push(iat);
        sphere_vert_tri.push(atoms[iat].x());
        sphere_vert_tri.push(atoms[iat].y());
        sphere_vert_tri.push(atoms[iat].z());
        for(ip=0;ip<colourScheme[atoms[iat]["_atom_site.id"]].length;ip++) sphere_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][ip]);
        //FIXME - Maybe should be outside loop...
        if(scalebyVDWRadii){
            sphere_sizes.push(size*atoms[iat]["_ccp4mg_vdw_radius"]);
        } else {
            sphere_sizes.push(size);
        }
        var atom = {};
        atom["x"] = atoms[iat].x();
        atom["y"] = atoms[iat].y();
        atom["z"] = atoms[iat].z();
        atom["tempFactor"] = atoms[iat]["_atom_site.B_iso_or_equiv"];
        atom["charge"] = atoms[iat]["_atom_site.pdbx_formal_charge"];
        atom["label"] =  atoms[iat].getAtomID();
        sphere_atoms.push(atom);
    }

    var spherePrimitiveInfo = {"atoms":[[sphere_atoms]],"sizes": [[sphere_sizes]], "col_tri":[[sphere_col_tri]], "norm_tri":[[[]]], "vert_tri":[[sphere_vert_tri]], "idx_tri":[[sphere_idx_tri]] , "prim_types":[[primType]] };
    return spherePrimitiveInfo;
}

function atomsToEllipsoidsInfo(atoms,size,colourScheme){
    var sphere_sizes = [];
    var sphere_col_tri = [];
    var sphere_vert_tri = [];
    var sphere_idx_tri = [];
    var scale_matrices = [];
    var sphere_atoms = [];
    var kBA = 1.3806503e-03;
    var TEMP = 298.155;
    var prefactor=8.0*Math.PI*Math.PI*kBA*TEMP;
    for(var iat=0;iat<atoms.length;iat++){
        sphere_idx_tri.push(iat);
        sphere_vert_tri.push(atoms[iat].x());
        sphere_vert_tri.push(atoms[iat].y());
        sphere_vert_tri.push(atoms[iat].z());
        if(typeof(atoms[iat]["_atom_site_anisotrop.U[1][1]"]) !== "undefined"){
            var anisoU = [];
            anisoU.push(atoms[iat]["_atom_site_anisotrop.U[1][1]"]);
            anisoU.push(atoms[iat]["_atom_site_anisotrop.U[1][2]"]);
            anisoU.push(atoms[iat]["_atom_site_anisotrop.U[1][3]"]);
            anisoU.push(atoms[iat]["_atom_site_anisotrop.U[1][2]"]);
            anisoU.push(atoms[iat]["_atom_site_anisotrop.U[2][2]"]);
            anisoU.push(atoms[iat]["_atom_site_anisotrop.U[2][3]"]);
            anisoU.push(atoms[iat]["_atom_site_anisotrop.U[1][3]"]);
            anisoU.push(atoms[iat]["_atom_site_anisotrop.U[2][3]"]);
            anisoU.push(atoms[iat]["_atom_site_anisotrop.U[3][3]"]);
            scale_matrices.push(anisoU);
        } else {
            // FIXME - We should be using B-factor (if that exists).
            var iso = [];
            iso.push(atoms[iat]["_atom_site.B_iso_or_equiv"]/prefactor);
            iso.push(0.0);
            iso.push(0.0);
            iso.push(0.0);
            iso.push(atoms[iat]["_atom_site.B_iso_or_equiv"]/prefactor);
            iso.push(0.0);
            iso.push(0.0);
            iso.push(0.0);
            iso.push(atoms[iat]["_atom_site.B_iso_or_equiv"]/prefactor);
            scale_matrices.push(iso);
        }
        for(ip=0;ip<colourScheme[atoms[iat]["_atom_site.id"]].length;ip++) sphere_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][ip]);
        sphere_sizes.push(size);
        var atom = {};
        atom["x"] = atoms[iat].x();
        atom["y"] = atoms[iat].y();
        atom["z"] = atoms[iat].z();
        atom["tempFactor"] = atoms[iat]["_atom_site.B_iso_or_equiv"];
        atom["charge"] = atoms[iat]["_atom_site.pdbx_formal_charge"];
        atom["label"] =  atoms[iat].getAtomID();
        sphere_atoms.push(atom);
    }

    var spherePrimitiveInfo = {"atoms":[[sphere_atoms]],"scale_matrices": [[scale_matrices]], "sizes": [[sphere_sizes]], "col_tri":[[sphere_col_tri]], "norm_tri":[[[]]], "vert_tri":[[sphere_vert_tri]], "idx_tri":[[sphere_idx_tri]] , "prim_types":[["SPHEROIDS"]] };
    return spherePrimitiveInfo;
}

function ringsToToruses(rings,size,radius,linetype,colourScheme){
    var toruses = [];

    var cylinder_sizes = [];
    var cylinder_col_tri = [];
    var cylinder_vert_tri = [];
    var cylinder_norm_tri = [];
    var cylinder_idx_tri = [];
    var radii = [];
    for(var i=0;i<rings.length;i++){
        var theRing = rings[i];
        var col = colourScheme[theRing[0][0]["_atom_site.id"]];
        // Always try to colour ring by colour of first C atom in the ring.
        for(var j=1;j<theRing[0].length;j++){
            if(theRing[0][1].element()=="C"){
                col = colourScheme[theRing[0][j]["_atom_site.id"]];
                break;
            }
        }
        cylinder_sizes.push(size);
        Array.prototype.push.apply(cylinder_col_tri,col);
        Array.prototype.push.apply(cylinder_vert_tri,theRing[2]);
        Array.prototype.push.apply(cylinder_norm_tri,theRing[3]);
        cylinder_idx_tri.push(i);
        if(theRing[0].length===6){
            radii.push(radius);
        } else {
            radii.push(radius*5.0/6.0);
        }
    }
    var cylinderPrimitiveInfo = {"radii": [[radii]], "sizes": [[cylinder_sizes]], "col_tri":[[cylinder_col_tri]], "norm_tri":[[cylinder_norm_tri]], "vert_tri":[[cylinder_vert_tri]], "idx_tri":[[cylinder_idx_tri]] , "prim_types":[[linetype]] };
    //console.log(cylinderPrimitiveInfo);
    return cylinderPrimitiveInfo;
}

function getMultipleBondsLines(atoms,enerLib,size,colourScheme){
    return getMultipleBonds(atoms,enerLib,size,colourScheme,"LINES");
}

function getMultipleBondsCylinders(atoms,enerLib,size,colourScheme){
    return getMultipleBonds(atoms,enerLib,size,colourScheme,"CYLINDERS");
}

function getMultipleBonds(atoms,enerLib,size,colourScheme,style){

    var thisResId = "";
    var thisChainId = "";
    var thisResAtoms = [];
    var contacts = [];
    var dashedContacts = [];
    var singleContacts = [];
    var tripleContacts = [];
    var drawRings = [];

    var xaxis = vec3.create([1.0,0.0,0.0]);
    var yaxis = vec3.create([0.0,1.0,0.0]);
    var zaxis = vec3.create([0.0,0.0,1.0]);

    // FIXME - This is for CYLINDERS, need different for balls, bonds, etc. And I am making numbers up.
    //var dbOffset = size*.5;

    var tbOffset = size*2./3.;
    var dbOffset = size*.675;
    var tbSize = size*.3;
    var dbSize = size*.4;
    var radius = 0.65;

    if(typeof(style)!=="undefined"&&style==="LINES"){
        tbOffset = size*0.02;
        dbOffset = size*0.02;
        tbSize = size;
        dbSize = size;
        radius = 0.65;
    }

    function GetFuncGroupBonds(name1,name2,bonds) {
        function getAdj(thisBond){
            return ((thisBond["_chem_comp_bond.atom_id_1"]===name1 && thisBond["_chem_comp_bond.atom_id_2"]!==name2) || (thisBond["_chem_comp_bond.atom_id_2"]===name1 && thisBond["_chem_comp_bond.atom_id_1"]!==name2));
        }
        var adjacent_bonds = bonds.filter(getAdj);
        var adjacent_children = [];
        for(var i=0;i<adjacent_bonds.length;i++){
            if(adjacent_bonds[i]["_chem_comp_bond.atom_id_1"]===name1){
                adjacent_children.push([adjacent_bonds[i]["_chem_comp_bond.atom_id_2"],adjacent_bonds[i]["_chem_comp_bond.type"].trim()]);
            } else {
                adjacent_children.push([adjacent_bonds[i]["_chem_comp_bond.atom_id_1"],adjacent_bonds[i]["_chem_comp_bond.type"].trim()]);
            }
        }
        return adjacent_children;
    }

    function GetFuncGroups(name1,name2,bonds) {
        function getAdj(thisBond){
            return ((thisBond["_chem_comp_bond.atom_id_1"]===name1 && thisBond["_chem_comp_bond.atom_id_2"]!==name2) || (thisBond["_chem_comp_bond.atom_id_2"]===name1 && thisBond["_chem_comp_bond.atom_id_1"]!==name2));
        }
        var adjacent_bonds = bonds.filter(getAdj);
        var adjacent_children = [];
        for(var i=0;i<adjacent_bonds.length;i++){
            if(adjacent_bonds[i]["_chem_comp_bond.atom_id_1"]===name1){
                adjacent_children.push(adjacent_bonds[i]["_chem_comp_bond.atom_id_2"]);
            } else {
                adjacent_children.push(adjacent_bonds[i]["_chem_comp_bond.atom_id_1"]);
            }
        }
        return adjacent_children;
    }

    function assignResidueBonding(resatoms){
        var resname = resatoms[0].residue.getName();
        if(resname in enerLib.monLibBonds){
            function bondAtomsInAtoms(thisBond){
                var name1 = thisBond["_chem_comp_bond.atom_id_1"].trim();
                var name2 = thisBond["_chem_comp_bond.atom_id_2"].trim();
                var at1 = null;
                var at2 = null;
                for(var iat=0;iat<resatoms.length;iat++){
                    if(resatoms[iat]["_atom_site.label_atom_id"].trim()===name1){
                        at1 = resatoms[iat];
                    }else if(resatoms[iat]["_atom_site.label_atom_id"].trim()===name2){
                        at2 = resatoms[iat];
                    }
                }
                if(at1&&at2){
                    return true;
                }
                return false;
            }
            var bonds = enerLib.monLibBonds[resname].filter(bondAtomsInAtoms);
            var all_ring_centres = [];

            for(ib=0;ib<bonds.length;ib++){
                var name1 = bonds[ib]["_chem_comp_bond.atom_id_1"].trim();
                var name2 = bonds[ib]["_chem_comp_bond.atom_id_2"].trim();
                var type = bonds[ib]["_chem_comp_bond.type"].trim();
                if(type!=="metal"){
                    var at1 = null;
                    var at2 = null;
                    var iat1 = -1;
                    var iat2 = -1;
                    for(var iat=0;iat<resatoms.length;iat++){
                        if(resatoms[iat]["_atom_site.label_atom_id"].trim()===name1){
                            at1 = resatoms[iat];
                            iat1 = iat;
                        }else if(resatoms[iat]["_atom_site.label_atom_id"].trim()===name2){
                            at2 = resatoms[iat];
                            iat2 = iat;
                        }
                    }
                    if(at1&&at2){
                        if(type==="triple"){
                            //console.log("triple");
                            //console.log(name1+" "+name2);
                            var thisFuncGroups = GetFuncGroups(name1,name2,bonds);
                            var otherFuncGroups = GetFuncGroups(name2,name1,bonds);
                            var adjFuncGroups = [];
                            var planeAtoms = [];
                            //console.log("thisFuncGroups");
                            //console.log(thisFuncGroups);
                            //console.log("otherFuncGroups");
                            //console.log(otherFuncGroups);
                            if(thisFuncGroups.length>0&&otherFuncGroups.length>0){
                                if(iat1>iat2){
                                    //console.log("Look for(adj) func groups on: "+name1+" "+thisFuncGroups[0]);
                                    adjFuncGroups = GetFuncGroups(thisFuncGroups[0],name1,bonds);
                                    planeAtoms.push(name1);
                                    planeAtoms.push(thisFuncGroups[0]);
                                } else {
                                    adjFuncGroups = GetFuncGroups(otherFuncGroups[0],name2,bonds);
                                    //console.log("Look for(adj) func groups on: "+name2+" "+otherFuncGroups[0]);
                                    planeAtoms.push(name2);
                                    planeAtoms.push(otherFuncGroups[0]);
                                }
                                //console.log("adjFuncGroups (case 1)");
                                //console.log(adjFuncGroups);
                            } else if(thisFuncGroups.length>0){
                                adjFuncGroups = GetFuncGroups(thisFuncGroups[0],name1,bonds);
                                planeAtoms.push(name1);
                                planeAtoms.push(thisFuncGroups[0]);
                            } else if(otherFuncGroups.length>0){
                                adjFuncGroups = GetFuncGroups(otherFuncGroups[0],name2,bonds);
                                planeAtoms.push(name2);
                                planeAtoms.push(otherFuncGroups[0]);
                            }
                            for(var iadj=0;iadj<adjFuncGroups.length;iadj++){
                                planeAtoms.push(adjFuncGroups[iadj]);
                            }
                            //console.log("planeAtoms");
                            //console.log(planeAtoms);

                            var atmp = vec3.create();
                            var carts0 = vec3.create([at1.x(),at1.y(),at1.z()]);
                            var carts1 = vec3.create([at2.x(),at2.y(),at2.z()]);
                            var l = vec3.create();
                            vec3.subtract(carts0,carts1,l);
                            vec3.normalize(l);

                            if(planeAtoms.length>2){
                                var l1 = vec3.create();
                                var l2 = vec3.create();
                                var planeAtom0at = at1.residue.getAtomTrimmed(planeAtoms[0]);
                                var planeAtom0cart = vec3.create([planeAtom0at.x(),planeAtom0at.y(),planeAtom0at.z()]);
                                var planeAtom1at = at1.residue.getAtomTrimmed(planeAtoms[1]);
                                var planeAtom1cart = vec3.create([planeAtom1at.x(),planeAtom1at.y(),planeAtom1at.z()]);
                                var planeAtom2at = at1.residue.getAtomTrimmed(planeAtoms[2]);
                                var planeAtom2cart = vec3.create([planeAtom2at.x(),planeAtom2at.y(),planeAtom2at.z()]);
                                vec3.subtract(planeAtom1cart,planeAtom0cart,l1);
                                vec3.subtract(planeAtom2cart,planeAtom1cart,l2);
                                vec3.normalize(l1);
                                vec3.normalize(l2);
                                var cross1 = vec3.create();
                                vec3.cross(l1,l2,cross1);
                                vec3.normalize(cross1);
                                vec3.cross(cross1,l,atmp);
                            } else {
                                vec3.cross(l,zaxis,atmp);
                                if(vec3.length(atmp) < 0.000001){
                                    vec3.cross(l,yaxis,atmp);
                                    if(vec3.length(atmp) < 0.000001){
                                        vec3.cross(l,xaxis,atmp);
                                    }
                                }
                            }

                            if(vec3.length(atmp) > 0.000001){
                                vec3.normalize(atmp);
                                atmp[0] *= tbOffset;
                                atmp[1] *= tbOffset;
                                atmp[2] *= tbOffset;

                                var acontact = [];
                                var thisLength = Model.prototype.bondLength(at1,at2);
                                var newAt1 = new Atom(at1);
                                newAt1["_atom_site.Cartn_x"] = at1.x()+atmp[0];
                                newAt1["_atom_site.Cartn_y"] = at1.y()+atmp[1];
                                newAt1["_atom_site.Cartn_z"] = at1.z()+atmp[2];
                                newAt1["_atom_site.id"] = at1["_atom_site.id"];
                                var newAt2 = new Atom(at2);
                                newAt2["_atom_site.Cartn_x"] = at2.x()+atmp[0];
                                newAt2["_atom_site.Cartn_y"] = at2.y()+atmp[1];
                                newAt2["_atom_site.Cartn_z"] = at2.z()+atmp[2];
                                newAt2["_atom_site.id"] = at2["_atom_site.id"];
                                acontact.push(thisLength);
                                acontact.push(newAt1);
                                acontact.push(newAt2);
                                tripleContacts.push(acontact);

                                var acontact2 = [];
                                var newAt3 = new Atom(at1);
                                newAt3["_atom_site.Cartn_x"] = at1.x()-atmp[0];
                                newAt3["_atom_site.Cartn_y"] = at1.y()-atmp[1];
                                newAt3["_atom_site.Cartn_z"] = at1.z()-atmp[2];
                                newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                var newAt4 = new Atom(at2);
                                newAt4["_atom_site.Cartn_x"] = at2.x()-atmp[0];
                                newAt4["_atom_site.Cartn_y"] = at2.y()-atmp[1];
                                newAt4["_atom_site.Cartn_z"] = at2.z()-atmp[2];
                                newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                acontact2.push(thisLength);
                                acontact2.push(newAt3);
                                acontact2.push(newAt4);
                                tripleContacts.push(acontact2);

                                var acontact3 = [];
                                acontact3.push(thisLength);
                                acontact3.push(at1);
                                acontact3.push(at2);
                                tripleContacts.push(acontact3);
                            }



                        } else if(type==="double"||type==="aromatic"||type==="deloc"||type==="aromat"){
                            function GetDistanceByName(n1,n2){
                                var n1at = at1.residue.getAtomTrimmed(n1);
                                if(!n1at) return 1000.0;
                                var n1cart = vec3.create([n1at.x(),n1at.y(),n1at.z()]);
                                var n2at = at1.residue.getAtomTrimmed(n2);
                                if(!n2at) return 1000.0;
                                var n2cart = vec3.create([n2at.x(),n2at.y(),n2at.z()]);
                                var diff = vec3.create();
                                vec3.subtract(n1cart,n2cart,diff);
                                return vec3.length(diff);
                            }
                            function checkPlanar(planarAtoms){
                                // blah, blah, check normals.
                                if(planarAtoms.length===5) return true; // Really ?!! 8/01/2017 !!
                                var nprev = null;
                                for(var ip=0;ip<planarAtoms.length-2;ip++){
                                    var n1at = at1.residue.getAtomTrimmed(planarAtoms[ip]);
                                    var n2at = at1.residue.getAtomTrimmed(planarAtoms[ip+1]);
                                    var n3at = at1.residue.getAtomTrimmed(planarAtoms[ip+2]);
                                    var b1 = vec3.create([n2at.x()-n1at.x(),n2at.y()-n1at.y(),n2at.z()-n1at.z()]);
                                    var b2 = vec3.create([n3at.x()-n2at.x(),n3at.y()-n2at.y(),n3at.z()-n2at.z()]);
                                    vec3.normalize(b1);
                                    vec3.normalize(b2);
                                    var n = vec3.create();
                                    vec3.cross(b1,b2,n);
                                    if(nprev){
                                        if(vec3.dot(n,nprev)<0.62) return false;
                                    }
                                    nprev = n;
                                }
                                return true;
                            }
                            function CheckRing(theName1,theName2,theBonds) {
                                var thisFuncGroups = GetFuncGroups(theName1,theName2,theBonds);
                                var otherFuncGroups = GetFuncGroups(theName2,theName1,theBonds);
                                var rings = [];

                                //console.log(thisFuncGroups);
                                //console.log(otherFuncGroups);

                                //console.log("other: "+theName2);
                                //console.log("start: "+theName1);
                                for(var it=0;it<thisFuncGroups.length;it++){
                                    var adjFuncGroups = GetFuncGroups(thisFuncGroups[it],theName1,bonds);
                                    //console.log("    neighb0: "+thisFuncGroups[it]);
                                    for(var nb1=0;nb1<adjFuncGroups.length;nb1++){
                                        if(GetDistanceByName(theName1,adjFuncGroups[nb1])<2.6){
                                            //console.log("        neighb1: "+adjFuncGroups[nb1]);
                                            //console.log("        "+adjFuncGroups[nb1]+" "+GetDistanceByName(theName1,adjFuncGroups[nb1]));
                                            var adjFuncGroups2 = GetFuncGroups(adjFuncGroups[nb1],thisFuncGroups[it],bonds);
                                            for(var nb2=0;nb2<adjFuncGroups2.length;nb2++){
                                                if(GetDistanceByName(theName1,adjFuncGroups2[nb2])<3.0){
                                                    //console.log("            neighb2: "+adjFuncGroups2[nb2]);
                                                    //console.log("            "+adjFuncGroups2[nb2]+" "+GetDistanceByName(theName1,adjFuncGroups2[nb2]));
                                                    var adjFuncGroups3 = GetFuncGroups(adjFuncGroups2[nb2],adjFuncGroups[nb1],bonds);
                                                    for(var nb3=0;nb3<adjFuncGroups3.length;nb3++){
                                                        if(GetDistanceByName(theName1,adjFuncGroups3[nb3])<2.6){
                                                            //console.log("                neighb3: "+adjFuncGroups3[nb3]);
                                                            //console.log("                "+adjFuncGroups3[nb3]+" "+GetDistanceByName(theName1,adjFuncGroups3[nb3]));
                                                            if(adjFuncGroups3[nb3]===theName2){
                                                                // 5-membered ring.
                                                                //console.log("                5-membered ring");
                                                                //console.log("                "+theName1+" "+thisFuncGroups[it]+" "+adjFuncGroups[nb1]+" "+adjFuncGroups2[nb2]+" "+theName2);
                                                                if(checkPlanar([theName1,thisFuncGroups[it],adjFuncGroups[nb1],adjFuncGroups2[nb2],adjFuncGroups3[nb3]])){
                                                                    rings.push([theName1,thisFuncGroups[it],adjFuncGroups[nb1],adjFuncGroups2[nb2],adjFuncGroups3[nb3]]);
                                                                }
                                                            } else {
                                                                var adjFuncGroups4 = GetFuncGroups(adjFuncGroups3[nb3],adjFuncGroups2[nb2],bonds);
                                                                for(var nb4=0;nb4<adjFuncGroups4.length;nb4++){
                                                                    if(adjFuncGroups4[nb4]===theName2&&checkPlanar([theName1,thisFuncGroups[it],adjFuncGroups[nb1],adjFuncGroups2[nb2],adjFuncGroups3[nb3],theName2])){
                                                                        //console.log("                    6-membered ring");
                                                                        //console.log("                    "+theName1+" "+thisFuncGroups[it]+" "+adjFuncGroups[nb1]+" "+adjFuncGroups2[nb2]+" "+adjFuncGroups3[nb3]+" "+theName2);
                                                                        rings.push([theName1,thisFuncGroups[it],adjFuncGroups[nb1],adjFuncGroups2[nb2],adjFuncGroups3[nb3],theName2]);
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                return rings;
                            }
                            var skipBond = false;
                            var rings = CheckRing(name1,name2,bonds);
                            if(rings.length>0&&rings[0].length>4){
                                for(var ir=0;ir<rings.length;ir++){
                                    var ring_centre = vec3.create([0.0, 0.0, 0.0]);
                                    for(var ir2=0;ir2<rings[ir].length;ir2++){
                                        var n1at = at1.residue.getAtomTrimmed(rings[ir][ir2]);
                                        ring_centre[0] += n1at.x();
                                        ring_centre[1] += n1at.y();
                                        ring_centre[2] += n1at.z();
                                    }
                                    ring_centre[0] /= rings[ir].length;
                                    ring_centre[1] /= rings[ir].length;
                                    ring_centre[2] /= rings[ir].length;
                                    var doneThisRing = false;
                                    for(var iar=0;iar<all_ring_centres.length;iar++){
                                        var diff = vec3.create();
                                        vec3.subtract(ring_centre,all_ring_centres[iar],diff);
                                        var dist = vec3.length(diff);
                                        if(dist<1e-4){
                                            doneThisRing = true;
                                        }
                                    }
                                    // We now have to do the "single bond" around edge of ring.

                                    var acontact = [];
                                    var thisLength = Model.prototype.bondLength(at1,at2);
                                    acontact.push(thisLength);
                                    acontact.push(at1);
                                    acontact.push(at2);
                                    singleContacts.push(acontact);
                                    
                                    if(!doneThisRing){
                                        all_ring_centres.push(ring_centre);
                                        //console.log("DRAW RING: "+name1+" "+name2);
                                        var drawRing = [];
                                        var drawAt = [];
                                        var ringNorm = vec3.create([0.0,0.0,0.0]);
                                        var nn = 0;
                                        for(var ir2=0;ir2<rings[ir].length;ir2++){
                                            var atDraw = at1.residue.getAtomTrimmed(rings[ir][ir2]);
                                            var cartDraw = vec3.create([atDraw.x(),atDraw.y(),atDraw.z()]);
                                            drawRing.push(cartDraw);
                                            drawAt.push(atDraw);
                                            if(ir2>1){
                                                var d1 = vec3.create();
                                                var d2 = vec3.create();
                                                vec3.subtract(drawRing[drawRing.length-1],drawRing[drawRing.length-2],d1);
                                                vec3.subtract(drawRing[drawRing.length-2],drawRing[drawRing.length-3],d2);
                                                //console.log(d1);
                                                //console.log(d2);
                                                var ncross = vec3.create();
                                                vec3.cross(d1,d2,ncross);
                                                //console.log(ncross);
                                                //console.log();
                                                vec3.add(ringNorm,ncross,ringNorm);
                                                nn += 1;
                                            }
                                        }
                                        vec3.normalize(ringNorm);
                                        drawRings.push([drawAt,drawRing,ring_centre,ringNorm]);
                                    }
                                }
                            } else {
                                var thisFuncGroups = GetFuncGroups(name1,name2,bonds);
                                var otherFuncGroups = GetFuncGroups(name2,name1,bonds);
                                if(thisFuncGroups.length==3||otherFuncGroups.length==3){
                                    console.log("DRAW DOUBLE CASE 1: "+name1+" "+name2);

                                    var planeCarts = [];
                                    var doublePlaneCarts = [];
                                    var singlePlaneCarts = [];
                                    var nDouble = 0;

                                    var conn_lists = [];
                                    if(thisFuncGroups.length==3){
                                        conn_lists = GetFuncGroupBonds(name1,name2,bonds);
                                    } else {
                                        conn_lists = GetFuncGroupBonds(name2,name1,bonds);
                                    }
                                    for(var k=0;k<conn_lists.length;k++){
                                        var planeAtom0at = at1.residue.getAtomTrimmed(conn_lists[k][0]);
                                        var planeAtom0cart = vec3.create([planeAtom0at.x(),planeAtom0at.y(),planeAtom0at.z()]);
                                        if(conn_lists[k][1]=="double"||conn_lists[k][1]=="deloc"||conn_lists[k][1]=="aromatic"||conn_lists[k][1]=="aromat"){
                                            doublePlaneCarts.push(planeAtom0cart);
                                            nDouble++;
                                        } else {
                                            singlePlaneCarts.push(planeAtom0cart);
                                        }
                                    }
                                    if(nDouble===0||nDouble===1){

                                        var mc = vec3.create();
                                        vec3.subtract(singlePlaneCarts[0],singlePlaneCarts[1],mc);
                                        vec3.normalize(mc);
                                        mc[0] *= dbOffset; mc[1] *= dbOffset; mc[2] *= dbOffset;

                                        var acontact2 = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        var newAt3 = new Atom(at1);
                                        newAt3["_atom_site.Cartn_x"] = at1.x()+mc[0];
                                        newAt3["_atom_site.Cartn_y"] = at1.y()+mc[1];
                                        newAt3["_atom_site.Cartn_z"] = at1.z()+mc[2];
                                        newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                        var newAt4 = new Atom(at2);
                                        newAt4["_atom_site.Cartn_x"] = at2.x()+mc[0];
                                        newAt4["_atom_site.Cartn_y"] = at2.y()+mc[1];
                                        newAt4["_atom_site.Cartn_z"] = at2.z()+mc[2];
                                        newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                        acontact2.push(thisLength);
                                        acontact2.push(newAt3);
                                        acontact2.push(newAt4);
                                        contacts.push(acontact2);

                                        var acontact = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        var newAt1 = new Atom(at1);
                                        newAt1["_atom_site.Cartn_x"] = at1.x()-mc[0];
                                        newAt1["_atom_site.Cartn_y"] = at1.y()-mc[1];
                                        newAt1["_atom_site.Cartn_z"] = at1.z()-mc[2];
                                        newAt1["_atom_site.id"] = at1["_atom_site.id"];
                                        var newAt2 = new Atom(at2);
                                        newAt2["_atom_site.Cartn_x"] = at2.x()-mc[0];
                                        newAt2["_atom_site.Cartn_y"] = at2.y()-mc[1];
                                        newAt2["_atom_site.Cartn_z"] = at2.z()-mc[2];
                                        newAt2["_atom_site.id"] = at2["_atom_site.id"];
                                        acontact.push(thisLength);
                                        acontact.push(newAt1);
                                        acontact.push(newAt2);
                                        if(type==="double"){
                                            contacts.push(acontact);
                                        } else {
                                            dashedContacts.push(acontact);
                                        }
                                    }
                                    if(nDouble==2){
                                        // FIXME - There should be length adjustments.
                                        var cross1 = vec3.create();
                                        var cross2 = vec3.create()
                                        var carts0 = vec3.create([at1.x(),at1.y(),at1.z()]);
                                        var carts1 = vec3.create([at2.x(),at2.y(),at2.z()]);
                                        if(thisFuncGroups.length==3){
                                            vec3.subtract(singlePlaneCarts[0],carts0,cross1);
                                            vec3.subtract(carts1,carts0,cross2);
                                        } else {
                                            vec3.subtract(singlePlaneCarts[0],carts1,cross1);
                                            vec3.subtract(carts0,carts1,cross2);
                                        }
                                        vec3.normalize(cross1);
                                        vec3.normalize(cross2);
                                        var mc = vec3.create();
                                        vec3.cross(cross1,cross2,mc);
                                        vec3.normalize(mc);
                                        mc[0] *= dbOffset; mc[1] *= dbOffset; mc[2] *= dbOffset;

                                        var acontact2 = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        var newAt3 = new Atom(at1);
                                        newAt3["_atom_site.Cartn_x"] = at1.x()+mc[0];
                                        newAt3["_atom_site.Cartn_y"] = at1.y()+mc[1];
                                        newAt3["_atom_site.Cartn_z"] = at1.z()+mc[2];
                                        newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                        var newAt4 = new Atom(at2);
                                        newAt4["_atom_site.Cartn_x"] = at2.x()+mc[0];
                                        newAt4["_atom_site.Cartn_y"] = at2.y()+mc[1];
                                        newAt4["_atom_site.Cartn_z"] = at2.z()+mc[2];
                                        newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                        acontact2.push(thisLength);
                                        acontact2.push(newAt3);
                                        acontact2.push(newAt4);
                                        contacts.push(acontact2);

                                        var acontact = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        var newAt1 = new Atom(at1);
                                        newAt1["_atom_site.Cartn_x"] = at1.x()-mc[0];
                                        newAt1["_atom_site.Cartn_y"] = at1.y()-mc[1];
                                        newAt1["_atom_site.Cartn_z"] = at1.z()-mc[2];
                                        newAt1["_atom_site.id"] = at1["_atom_site.id"];
                                        var newAt2 = new Atom(at2);
                                        newAt2["_atom_site.Cartn_x"] = at2.x()-mc[0];
                                        newAt2["_atom_site.Cartn_y"] = at2.y()-mc[1];
                                        newAt2["_atom_site.Cartn_z"] = at2.z()-mc[2];
                                        newAt2["_atom_site.id"] = at2["_atom_site.id"];
                                        acontact.push(thisLength);
                                        acontact.push(newAt1);
                                        acontact.push(newAt2);
                                        if(type==="double"){
                                            contacts.push(acontact);
                                        } else {
                                            dashedContacts.push(acontact);
                                        }
                                    }

                                    if(nDouble==3){
                                        // FIXME - There should be length adjustments.
                                        var carts0 = vec3.create([at1.x(),at1.y(),at1.z()]);
                                        var carts1 = vec3.create([at2.x(),at2.y(),at2.z()]);
                                        var cross1 = vec3.create();
                                        var cross2 = vec3.create();
                                        if(thisFuncGroups.length==3){
                                            // FIXME if(j<2), ??!! see build_tree_primitives.
                                            vec3.subtract(doublePlaneCarts[1],doublePlaneCarts[2],cross1);
                                            vec3.subtract(carts1,carts0,cross2);
                                        } else {
                                            // FIXME if(i<conn_lists[conn_lists[i][j]][2]), ??!! see build_tree_primitives.
                                            vec3.subtract(doublePlaneCarts[1],doublePlaneCarts[2],cross1);
                                            vec3.subtract(carts0,carts1,cross2);
                                        }

                                        vec3.normalize(cross1);
                                        vec3.normalize(cross2);
                                        var mc = vec3.create();
                                        vec3.cross(cross1,cross2,mc);
                                        vec3.normalize(mc);
                                        mc[0] *= dbOffset; mc[1] *= dbOffset; mc[2] *= dbOffset;
                                        
                                        var acontact2 = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        var newAt3 = new Atom(at1);
                                        newAt3["_atom_site.Cartn_x"] = at1.x()+mc[0];
                                        newAt3["_atom_site.Cartn_y"] = at1.y()+mc[1];
                                        newAt3["_atom_site.Cartn_z"] = at1.z()+mc[2];
                                        newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                        var newAt4 = new Atom(at2);
                                        newAt4["_atom_site.Cartn_x"] = at2.x()+mc[0];
                                        newAt4["_atom_site.Cartn_y"] = at2.y()+mc[1];
                                        newAt4["_atom_site.Cartn_z"] = at2.z()+mc[2];
                                        newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                        acontact2.push(thisLength);
                                        acontact2.push(newAt3);
                                        acontact2.push(newAt4);
                                        contacts.push(acontact2);

                                        var acontact = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        var newAt1 = new Atom(at1);
                                        newAt1["_atom_site.Cartn_x"] = at1.x()-mc[0];
                                        newAt1["_atom_site.Cartn_y"] = at1.y()-mc[1];
                                        newAt1["_atom_site.Cartn_z"] = at1.z()-mc[2];
                                        newAt1["_atom_site.id"] = at1["_atom_site.id"];
                                        var newAt2 = new Atom(at2);
                                        newAt2["_atom_site.Cartn_x"] = at2.x()-mc[0];
                                        newAt2["_atom_site.Cartn_y"] = at2.y()-mc[1];
                                        newAt2["_atom_site.Cartn_z"] = at2.z()-mc[2];
                                        newAt2["_atom_site.id"] = at2["_atom_site.id"];
                                        acontact.push(thisLength);
                                        acontact.push(newAt1);
                                        acontact.push(newAt2);
                                        if(type==="double"){
                                            contacts.push(acontact);
                                        } else {
                                            dashedContacts.push(acontact);
                                        }
                                    }
                                }
                                if(thisFuncGroups.length==2&&(otherFuncGroups.length<2||(otherFuncGroups.length==2&&(iat1<iat2)))){
                                    console.log("DRAW DOUBLE CASE 2: "+name1+" "+name2);
                                    // FIXME - "other double" adjustments

                                    var mc = vec3.create();
                                    var thisAt0 = at1.residue.getAtomTrimmed(thisFuncGroups[0]);
                                    var cart0 = vec3.create([thisAt0.x(),thisAt0.y(),thisAt0.z()]);
                                    var thisAt1 = at1.residue.getAtomTrimmed(thisFuncGroups[1]);
                                    var cart1 = vec3.create([thisAt1.x(),thisAt1.y(),thisAt1.z()]);
                                    vec3.subtract(cart0,cart1,mc);
                                    vec3.normalize(mc);
                                    mc[0] *= dbOffset; mc[1] *= dbOffset; mc[2] *= dbOffset;

                                    var acontact2 = [];
                                    var thisLength = Model.prototype.bondLength(at1,at2);
                                    var newAt3 = new Atom(at1);
                                    newAt3["_atom_site.Cartn_x"] = at1.x()+mc[0];
                                    newAt3["_atom_site.Cartn_y"] = at1.y()+mc[1];
                                    newAt3["_atom_site.Cartn_z"] = at1.z()+mc[2];
                                    newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                    var newAt4 = new Atom(at2);
                                    newAt4["_atom_site.Cartn_x"] = at2.x()+mc[0];
                                    newAt4["_atom_site.Cartn_y"] = at2.y()+mc[1];
                                    newAt4["_atom_site.Cartn_z"] = at2.z()+mc[2];
                                    newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                    acontact2.push(thisLength);
                                    acontact2.push(newAt3);
                                    acontact2.push(newAt4);
                                    contacts.push(acontact2);

                                    var acontact = [];
                                    var thisLength = Model.prototype.bondLength(at1,at2);
                                    var newAt1 = new Atom(at1);
                                    newAt1["_atom_site.Cartn_x"] = at1.x()-mc[0];
                                    newAt1["_atom_site.Cartn_y"] = at1.y()-mc[1];
                                    newAt1["_atom_site.Cartn_z"] = at1.z()-mc[2];
                                    newAt1["_atom_site.id"] = at1["_atom_site.id"];
                                    var newAt2 = new Atom(at2);
                                    newAt2["_atom_site.Cartn_x"] = at2.x()-mc[0];
                                    newAt2["_atom_site.Cartn_y"] = at2.y()-mc[1];
                                    newAt2["_atom_site.Cartn_z"] = at2.z()-mc[2];
                                    newAt2["_atom_site.id"] = at2["_atom_site.id"];
                                    acontact.push(thisLength);
                                    acontact.push(newAt1);
                                    acontact.push(newAt2);
                                    if(type==="double"){
                                        contacts.push(acontact);
                                    } else {
                                        dashedContacts.push(acontact);
                                    }
                                }
                                if((thisFuncGroups.length<2||(thisFuncGroups.length==2&&(iat1>iat2)))&&otherFuncGroups.length==2){
                                    console.log("DRAW DOUBLE CASE 3: "+name1+" "+name2);
                                    // FIXME - "other double" adjustments

                                    //console.log(thisFuncGroups);
                                    //console.log(otherFuncGroups);

                                    var mc = vec3.create();
                                    var otherAt0 = at1.residue.getAtomTrimmed(otherFuncGroups[0]);
                                    var cart0 = vec3.create([otherAt0.x(),otherAt0.y(),otherAt0.z()]);
                                    var otherAt1 = at1.residue.getAtomTrimmed(otherFuncGroups[1]);
                                    var cart1 = vec3.create([otherAt1.x(),otherAt1.y(),otherAt1.z()]);
                                    vec3.subtract(cart0,cart1,mc);
                                    vec3.normalize(mc);
                                    mc[0] *= dbOffset; mc[1] *= dbOffset; mc[2] *= dbOffset;

                                    var acontact2 = [];
                                    var thisLength = Model.prototype.bondLength(at1,at2);
                                    var newAt3 = new Atom(at1);
                                    newAt3["_atom_site.Cartn_x"] = at1.x()+mc[0];
                                    newAt3["_atom_site.Cartn_y"] = at1.y()+mc[1];
                                    newAt3["_atom_site.Cartn_z"] = at1.z()+mc[2];
                                    newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                    var newAt4 = new Atom(at2);
                                    newAt4["_atom_site.Cartn_x"] = at2.x()+mc[0];
                                    newAt4["_atom_site.Cartn_y"] = at2.y()+mc[1];
                                    newAt4["_atom_site.Cartn_z"] = at2.z()+mc[2];
                                    newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                    acontact2.push(thisLength);
                                    acontact2.push(newAt3);
                                    acontact2.push(newAt4);
                                    contacts.push(acontact2);

                                    var acontact = [];
                                    var thisLength = Model.prototype.bondLength(at1,at2);
                                    var newAt1 = new Atom(at1);
                                    newAt1["_atom_site.Cartn_x"] = at1.x()-mc[0];
                                    newAt1["_atom_site.Cartn_y"] = at1.y()-mc[1];
                                    newAt1["_atom_site.Cartn_z"] = at1.z()-mc[2];
                                    newAt1["_atom_site.id"] = at1["_atom_site.id"];
                                    var newAt2 = new Atom(at2);
                                    newAt2["_atom_site.Cartn_x"] = at2.x()-mc[0];
                                    newAt2["_atom_site.Cartn_y"] = at2.y()-mc[1];
                                    newAt2["_atom_site.Cartn_z"] = at2.z()-mc[2];
                                    newAt2["_atom_site.id"] = at2["_atom_site.id"];
                                    acontact.push(thisLength);
                                    acontact.push(newAt1);
                                    acontact.push(newAt2);
                                    if(type==="double"){
                                        contacts.push(acontact);
                                    } else {
                                        dashedContacts.push(acontact);
                                    }
                                }
                                if(thisFuncGroups.length==1&&otherFuncGroups.length==1){
                                    var otherAt0 = at1.residue.getAtomTrimmed(otherFuncGroups[0]);
                                    var thisAt0 = at1.residue.getAtomTrimmed(thisFuncGroups[0]);
                                    if(!otherAt0){
                                    }
                                    var c1 = vec3.create([otherAt0.x(),otherAt0.y(),otherAt0.z()]);
                                    var c2 = vec3.create([at1.x(),at1.y(),at1.z()]);
                                    var c3 = vec3.create([at2.x(),at2.y(),at2.z()]);
                                    var c4 = vec3.create([thisAt0.x(),thisAt0.y(),thisAt0.z()]);
                                    if(Math.abs(DihedralAngle(c1,c2,c3,c4))>Math.PI/2.0){
                                        // TRANS
                                        var c4c1 = vec3.create();
                                        var c2c3 = vec3.create();
                                        vec3.subtract(c4,c1,c4c1);
                                        vec3.subtract(c2,c3,c2c3);
                                        var cross1 = vec3.create();
                                        vec3.cross(c4c1,c2c3,cross1);
                                        vec3.normalize(cross1);
                                        var mc = vec3.create();
                                        vec3.normalize(c2c3);
                                        vec3.cross(cross1,c2c3,mc);
                                        vec3.normalize(mc);
                                        mc[0] *= dbOffset; mc[1] *= dbOffset; mc[2] *= dbOffset;
                                        var acontact2 = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        var newAt3 = new Atom(at1);
                                        newAt3["_atom_site.Cartn_x"] = at1.x()+mc[0];
                                        newAt3["_atom_site.Cartn_y"] = at1.y()+mc[1];
                                        newAt3["_atom_site.Cartn_z"] = at1.z()+mc[2];
                                        newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                        var newAt4 = new Atom(at2);
                                        newAt4["_atom_site.Cartn_x"] = at2.x()+mc[0];
                                        newAt4["_atom_site.Cartn_y"] = at2.y()+mc[1];
                                        newAt4["_atom_site.Cartn_z"] = at2.z()+mc[2];
                                        newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                        acontact2.push(thisLength);
                                        acontact2.push(newAt3);
                                        acontact2.push(newAt4);
                                        contacts.push(acontact2);

                                        var acontact = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        var newAt1 = new Atom(at1);
                                        newAt1["_atom_site.Cartn_x"] = at1.x()-mc[0];
                                        newAt1["_atom_site.Cartn_y"] = at1.y()-mc[1];
                                        newAt1["_atom_site.Cartn_z"] = at1.z()-mc[2];
                                        newAt1["_atom_site.id"] = at1["_atom_site.id"];
                                        var newAt2 = new Atom(at2);
                                        newAt2["_atom_site.Cartn_x"] = at2.x()-mc[0];
                                        newAt2["_atom_site.Cartn_y"] = at2.y()-mc[1];
                                        newAt2["_atom_site.Cartn_z"] = at2.z()-mc[2];
                                        newAt2["_atom_site.id"] = at2["_atom_site.id"];
                                        acontact.push(thisLength);
                                        acontact.push(newAt1);
                                        acontact.push(newAt2);
                                        if(type==="double"){
                                            contacts.push(acontact);
                                        } else {
                                            dashedContacts.push(acontact);
                                        }
                                    } else {
                                        // CIS
                                        // FIXME - There should be length adjustments.
                                        var midpoint1 = vec3.create([.5*(c4[0]+c1[0]),.5*(c4[1]+c1[1]),.5*(c4[2]+c1[2])]);
                                        var midpoint2 = vec3.create([.5*(c2[0]+c3[0]),.5*(c2[1]+c3[1]),.5*(c2[2]+c3[2])]);
                                        var mc = vec3.create();
                                        vec3.subtract(midpoint1,midpoint2,mc);
                                        vec3.normalize(mc);
                                        mc[0] *= dbOffset; mc[1] *= dbOffset; mc[2] *= dbOffset;
                                        var acontact2 = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        var newAt3 = new Atom(at1);
                                        newAt3["_atom_site.Cartn_x"] = at1.x()+2.0*mc[0];
                                        newAt3["_atom_site.Cartn_y"] = at1.y()+2.0*mc[1];
                                        newAt3["_atom_site.Cartn_z"] = at1.z()+2.0*mc[2];
                                        newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                        var newAt4 = new Atom(at2);
                                        newAt4["_atom_site.Cartn_x"] = at2.x()+2.0*mc[0];
                                        newAt4["_atom_site.Cartn_y"] = at2.y()+2.0*mc[1];
                                        newAt4["_atom_site.Cartn_z"] = at2.z()+2.0*mc[2];
                                        newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                        acontact2.push(thisLength);
                                        acontact2.push(newAt3);
                                        acontact2.push(newAt4);
                                        contacts.push(acontact2);

                                        var acontact = [];
                                        var thisLength = Model.prototype.bondLength(at1,at2);
                                        acontact.push(thisLength);
                                        acontact.push(at1);
                                        acontact.push(at2);
                                        if(type==="double"){
                                            contacts.push(acontact);
                                        } else {
                                            dashedContacts.push(acontact);
                                        }
                                    }
                                }
                                if((thisFuncGroups.length==1&&otherFuncGroups.length==0)||(thisFuncGroups.length==0&&otherFuncGroups.length==1)){
                                    console.log("DRAW DOUBLE CASE 5: "+name1+" "+name2);
                                    var c2 = vec3.create([at1.x(),at1.y(),at1.z()]);
                                    var c3 = vec3.create([at2.x(),at2.y(),at2.z()]);
                                    var l = vec3.create();

                                    if((thisFuncGroups.length==1&&otherFuncGroups.length==0)){
                                        vec3.subtract(c2,c3,l);
                                    } else {
                                        vec3.subtract(c3,c2,l);
                                    }
                                    vec3.normalize(l);

                                    var mc = vec3.create();
                                    vec3.cross(l,zaxis,mc);
                                    if(vec3.length(mc) < 0.000001){
                                        vec3.cross(l,yaxis,mc);
                                        if(vec3.length(mc) < 0.000001){
                                            vec3.cross(l,xaxis,mc);
                                        }
                                    }
                                    vec3.normalize(mc);
                                    mc[0] *= dbOffset; mc[1] *= dbOffset; mc[2] *= dbOffset;
                                    var acontact2 = [];
                                    var thisLength = Model.prototype.bondLength(at1,at2);
                                    var newAt3 = new Atom(at1);
                                    newAt3["_atom_site.Cartn_x"] = at1.x()+mc[0];
                                    newAt3["_atom_site.Cartn_y"] = at1.y()+mc[1];
                                    newAt3["_atom_site.Cartn_z"] = at1.z()+mc[2];
                                    newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                    var newAt4 = new Atom(at2);
                                    newAt4["_atom_site.Cartn_x"] = at2.x()+mc[0];
                                    newAt4["_atom_site.Cartn_y"] = at2.y()+mc[1];
                                    newAt4["_atom_site.Cartn_z"] = at2.z()+mc[2];
                                    newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                    acontact2.push(thisLength);
                                    acontact2.push(newAt3);
                                    acontact2.push(newAt4);
                                    contacts.push(acontact2);

                                    var acontact = [];
                                    var thisLength = Model.prototype.bondLength(at1,at2);
                                    var newAt1 = new Atom(at1);
                                    newAt1["_atom_site.Cartn_x"] = at1.x()-mc[0];
                                    newAt1["_atom_site.Cartn_y"] = at1.y()-mc[1];
                                    newAt1["_atom_site.Cartn_z"] = at1.z()-mc[2];
                                    newAt1["_atom_site.id"] = at1["_atom_site.id"];
                                    var newAt2 = new Atom(at2);
                                    newAt2["_atom_site.Cartn_x"] = at2.x()-mc[0];
                                    newAt2["_atom_site.Cartn_y"] = at2.y()-mc[1];
                                    newAt2["_atom_site.Cartn_z"] = at2.z()-mc[2];
                                    newAt2["_atom_site.id"] = at2["_atom_site.id"];
                                    acontact.push(thisLength);
                                    acontact.push(newAt1);
                                    acontact.push(newAt2);
                                    if(type==="double"){
                                        contacts.push(acontact);
                                    } else {
                                        dashedContacts.push(acontact);
                                    }
                                }
                                if((thisFuncGroups.length==0&&otherFuncGroups.length==0)) {
                                    console.log("DRAW DOUBLE CASE 6: "+name1+" "+name2);
                                    var c2 = vec3.create([at1.x(),at1.y(),at1.z()]);
                                    var c3 = vec3.create([at2.x(),at2.y(),at2.z()]);
                                    var l = vec3.create();

                                    vec3.subtract(c2,c3,l);
                                    vec3.normalize(l);

                                    var mc = vec3.create();
                                    vec3.cross(l,zaxis,mc);
                                    if(vec3.length(mc) < 0.000001){
                                        vec3.cross(l,yaxis,mc);
                                        if(vec3.length(mc) < 0.000001){
                                            vec3.cross(l,xaxis,mc);
                                        }
                                    }
                                    vec3.normalize(mc);
                                    mc[0] *= dbOffset; mc[1] *= dbOffset; mc[2] *= dbOffset;
                                    var acontact2 = [];
                                    var thisLength = Model.prototype.bondLength(at1,at2);
                                    var newAt3 = new Atom(at1);
                                    newAt3["_atom_site.Cartn_x"] = at1.x()+mc[0];
                                    newAt3["_atom_site.Cartn_y"] = at1.y()+mc[1];
                                    newAt3["_atom_site.Cartn_z"] = at1.z()+mc[2];
                                    newAt3["_atom_site.id"] = at1["_atom_site.id"];
                                    var newAt4 = new Atom(at2);
                                    newAt4["_atom_site.Cartn_x"] = at2.x()+mc[0];
                                    newAt4["_atom_site.Cartn_y"] = at2.y()+mc[1];
                                    newAt4["_atom_site.Cartn_z"] = at2.z()+mc[2];
                                    newAt4["_atom_site.id"] = at2["_atom_site.id"];
                                    acontact2.push(thisLength);
                                    acontact2.push(newAt3);
                                    acontact2.push(newAt4);
                                    contacts.push(acontact2);

                                    var acontact = [];
                                    var thisLength = Model.prototype.bondLength(at1,at2);
                                    var newAt1 = new Atom(at1);
                                    newAt1["_atom_site.Cartn_x"] = at1.x()-mc[0];
                                    newAt1["_atom_site.Cartn_y"] = at1.y()-mc[1];
                                    newAt1["_atom_site.Cartn_z"] = at1.z()-mc[2];
                                    newAt1["_atom_site.id"] = at1["_atom_site.id"];
                                    var newAt2 = new Atom(at2);
                                    newAt2["_atom_site.Cartn_x"] = at2.x()-mc[0];
                                    newAt2["_atom_site.Cartn_y"] = at2.y()-mc[1];
                                    newAt2["_atom_site.Cartn_z"] = at2.z()-mc[2];
                                    newAt2["_atom_site.id"] = at2["_atom_site.id"];
                                    acontact.push(thisLength);
                                    acontact.push(newAt1);
                                    acontact.push(newAt2);
                                    if(type==="double"){
                                        contacts.push(acontact);
                                    } else {
                                        dashedContacts.push(acontact);
                                    }
                                }
                            }
                        } else {
                            var acontact = [];
                            var thisLength = Model.prototype.bondLength(at1,at2);
                            acontact.push(thisLength);
                            acontact.push(at1);
                            acontact.push(at2);
                            singleContacts.push(acontact);
                        }
                    }
                }
            }
        }
    }

    var ligandResidues = [];
    for(var iat=0;iat<atoms.length;iat++){
        if(atoms[iat].residue.getResidueID()!==thisResId||atoms[iat].getChainID()!=thisChainID){
            thisResId = atoms[iat].residue.getResidueID();
            thisChainID = atoms[iat].getChainID();
            if(thisResAtoms.length>0){
                assignResidueBonding(thisResAtoms);
                ligandResidues.push(thisResAtoms[0].residue);
            }
            thisResAtoms = [];
        }
        thisResAtoms.push(atoms[iat]);
    }

    if(thisResAtoms.length>0){
        assignResidueBonding(thisResAtoms);
        ligandResidues.push(thisResAtoms[0].residue);
    }

    var interResCylinders = {};
    if(ligandResidues.length>0){
        var all_lig_lig_contacts = [];
        var model = ligandResidues[0].chain.model;
        for(var ilig=0;ilig<ligandResidues.length;ilig++){
            for(var jlig=0;jlig<ilig;jlig++){
                var lig_lig_contacts = model.SeekContacts(ligandResidues[ilig].atoms,ligandResidues[jlig].atoms,0.6,1.8);
                console.log(lig_lig_contacts);
                 Array.prototype.push.apply(all_lig_lig_contacts,lig_lig_contacts);
            }
        }
        if(typeof(style)=="undefined"){
            interResCylinders = contactsToCylindersLinesInfo(all_lig_lig_contacts,size,"CYLINDERS",colourScheme,false);
        } else {
            interResCylinders = contactsToCylindersLinesInfo(all_lig_lig_contacts,size,style,colourScheme,false);
        }
    }

    var cylinders = [];
    var dashCylinders = [];
    var tripleCylinders = [];
    var singleCylinders = [];

    if(typeof(style)=="undefined"){
        cylinders = contactsToCylindersLinesInfo(contacts,dbSize,"CYLINDERS",colourScheme,false);
        dashCylinders = contactsToCylindersLinesInfo(dashedContacts,dbSize,"CYLINDERS",colourScheme,true);
        tripleCylinders = contactsToCylindersLinesInfo(tripleContacts,tbSize,"CYLINDERS",colourScheme,false);
        singleCylinders = contactsToCylindersLinesInfo(singleContacts,size,"CYLINDERS",colourScheme,false);
        //console.log(cylinders);
        //console.log(dashCylinders);
        //console.log(tripleCylinders);
        //console.log(singleCylinders);

        toruses = ringsToToruses(drawRings,size,radius,"TORUSES",colourScheme);
    } else {
        cylinders = contactsToCylindersLinesInfo(contacts,dbSize,style,colourScheme,false);
        dashCylinders = contactsToCylindersLinesInfo(dashedContacts,dbSize,style,colourScheme,true);
        tripleCylinders = contactsToCylindersLinesInfo(tripleContacts,tbSize,style,colourScheme,false);
        singleCylinders = contactsToCylindersLinesInfo(singleContacts,size,style,colourScheme,false);
        if(style==="LINES"){
            toruses = ringsToToruses(drawRings,size,radius,"CIRCLES",colourScheme);
        } else {
            toruses = ringsToToruses(drawRings,0.1,radius,"TORUSES",colourScheme);
        }
    }

    var ret = [cylinders,dashCylinders,tripleCylinders,singleCylinders,toruses]; 
    if(typeof(interResCylinders["norm_tri"])!=="undefined"){
        ret.push(interResCylinders);
    }
    return ret;
}

function getBaseBlocks(atoms,size,colourScheme){
    var block_col_tri = [];
    var block_norm_tri = [];
    var block_vert_tri = [];
    var block_idx_tri = [];

    var up = vec3.create();
    var n1c2 = vec3.create();
    var c2n3 = vec3.create();

    var n1pcart = vec3.create();
    var c2pcart = vec3.create();
    var n3pcart = vec3.create();
    var c4pcart = vec3.create();
    var c5pcart = vec3.create();
    var c6pcart = vec3.create();
    var n1mcart = vec3.create();
    var c2mcart = vec3.create();
    var n3mcart = vec3.create();
    var c4mcart = vec3.create();
    var c5mcart = vec3.create();
    var c6mcart = vec3.create();

    var n9pcart = vec3.create();
    var c8pcart = vec3.create();
    var n7pcart = vec3.create();
    var n9mcart = vec3.create();
    var c8mcart = vec3.create();
    var n7mcart = vec3.create();

    var midpointp = vec3.create();
    var midpointm = vec3.create();
    var mid1 = vec3.create();
    var mid2 = vec3.create();

    var idx = 0;
    for(var iat=0;iat<atoms.length;iat++){
        if(atoms[iat]["_atom_site.label_atom_id"].trim()==="C5*"||atoms[iat]["_atom_site.label_atom_id"].trim()==="C5'"){
            var res1 = atoms[iat].residue;
            var n1 = res1.getAtomTrimmed("N1");
            var c2 = res1.getAtomTrimmed("C2");
            var n3 = res1.getAtomTrimmed("N3");
            var c4 = res1.getAtomTrimmed("C4");
            var c5 = res1.getAtomTrimmed("C5");
            var c6 = res1.getAtomTrimmed("C6");
            var n7 = res1.getAtomTrimmed("N7");
            var c8 = res1.getAtomTrimmed("C8");
            var n9 = res1.getAtomTrimmed("N9");
            if(n1&&c2&&n3&&c4&&c5&&c6){
                var n1cart = vec3.create([n1.x(),n1.y(),n1.z()]);
                var c2cart = vec3.create([c2.x(),c2.y(),c2.z()]);
                var n3cart = vec3.create([n3.x(),n3.y(),n3.z()]);
                var c4cart = vec3.create([c4.x(),c4.y(),c4.z()]);
                var c5cart = vec3.create([c5.x(),c5.y(),c5.z()]);
                var c6cart = vec3.create([c6.x(),c6.y(),c6.z()]);

                vec3.subtract(n1cart,c2cart,n1c2);
                vec3.subtract(c2cart,n3cart,c2n3);
                vec3.cross(n1c2,c2n3,up);
                up[0] *= size*.5;
                up[1] *= size*.5;
                up[2] *= size*.5;

                vec3.add(n1cart,up,n1pcart);
                vec3.add(c2cart,up,c2pcart);
                vec3.add(n3cart,up,n3pcart);
                vec3.add(c4cart,up,c4pcart);
                vec3.add(c5cart,up,c5pcart);
                vec3.add(c6cart,up,c6pcart);

                vec3.subtract(n1cart,up,n1mcart);
                vec3.subtract(c2cart,up,c2mcart);
                vec3.subtract(n3cart,up,n3mcart);
                vec3.subtract(c4cart,up,c4mcart);
                vec3.subtract(c5cart,up,c5mcart);
                vec3.subtract(c6cart,up,c6mcart);

                var midpoint = vec3.create([0.0,0.0,0.0]);
                vec3.add(midpoint,n1cart);
                vec3.add(midpoint,c2cart);
                vec3.add(midpoint,n3cart);
                vec3.add(midpoint,c4cart);
                vec3.add(midpoint,c5cart);
                vec3.add(midpoint,c6cart);
                midpoint[0] /= 6.0;
                midpoint[1] /= 6.0;
                midpoint[2] /= 6.0;

                vec3.add(midpoint,up,midpointp);
                vec3.subtract(midpoint,up,midpointm);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(n1pcart[0]); block_vert_tri.push(n1pcart[1]); block_vert_tri.push(n1pcart[2]);
                block_vert_tri.push(c2pcart[0]); block_vert_tri.push(c2pcart[1]); block_vert_tri.push(c2pcart[2]);
                block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c2pcart[0]); block_vert_tri.push(c2pcart[1]); block_vert_tri.push(c2pcart[2]);
                block_vert_tri.push(n3pcart[0]); block_vert_tri.push(n3pcart[1]); block_vert_tri.push(n3pcart[2]);
                block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(n3pcart[0]); block_vert_tri.push(n3pcart[1]); block_vert_tri.push(n3pcart[2]);
                block_vert_tri.push(c4pcart[0]); block_vert_tri.push(c4pcart[1]); block_vert_tri.push(c4pcart[2]);
                block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c4pcart[0]); block_vert_tri.push(c4pcart[1]); block_vert_tri.push(c4pcart[2]);
                block_vert_tri.push(c5pcart[0]); block_vert_tri.push(c5pcart[1]); block_vert_tri.push(c5pcart[2]);
                block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c5pcart[0]); block_vert_tri.push(c5pcart[1]); block_vert_tri.push(c5pcart[2]);
                block_vert_tri.push(c6pcart[0]); block_vert_tri.push(c6pcart[1]); block_vert_tri.push(c6pcart[2]);
                block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c6pcart[0]); block_vert_tri.push(c6pcart[1]); block_vert_tri.push(c6pcart[2]);
                block_vert_tri.push(n1pcart[0]); block_vert_tri.push(n1pcart[1]); block_vert_tri.push(n1pcart[2]);
                block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(n1mcart[0]); block_vert_tri.push(n1mcart[1]); block_vert_tri.push(n1mcart[2]);
                block_vert_tri.push(c2mcart[0]); block_vert_tri.push(c2mcart[1]); block_vert_tri.push(c2mcart[2]);
                block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c2mcart[0]); block_vert_tri.push(c2mcart[1]); block_vert_tri.push(c2mcart[2]);
                block_vert_tri.push(n3mcart[0]); block_vert_tri.push(n3mcart[1]); block_vert_tri.push(n3mcart[2]);
                block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(n3mcart[0]); block_vert_tri.push(n3mcart[1]); block_vert_tri.push(n3mcart[2]);
                block_vert_tri.push(c4mcart[0]); block_vert_tri.push(c4mcart[1]); block_vert_tri.push(c4mcart[2]);
                block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c4mcart[0]); block_vert_tri.push(c4mcart[1]); block_vert_tri.push(c4mcart[2]);
                block_vert_tri.push(c5mcart[0]); block_vert_tri.push(c5mcart[1]); block_vert_tri.push(c5mcart[2]);
                block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c5mcart[0]); block_vert_tri.push(c5mcart[1]); block_vert_tri.push(c5mcart[2]);
                block_vert_tri.push(c6mcart[0]); block_vert_tri.push(c6mcart[1]); block_vert_tri.push(c6mcart[2]);
                block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c6mcart[0]); block_vert_tri.push(c6mcart[1]); block_vert_tri.push(c6mcart[2]);
                block_vert_tri.push(n1mcart[0]); block_vert_tri.push(n1mcart[1]); block_vert_tri.push(n1mcart[2]);
                block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                // The edges
                vec3.subtract(n1cart,midpoint,mid1);
                vec3.subtract(c2cart,midpoint,mid2);
                vec3.add(mid1,mid2);
                vec3.normalize(mid1);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(n1mcart[0]); block_vert_tri.push(n1mcart[1]); block_vert_tri.push(n1mcart[2]);
                block_vert_tri.push(c2pcart[0]); block_vert_tri.push(c2pcart[1]); block_vert_tri.push(c2pcart[2]);
                block_vert_tri.push(n1pcart[0]); block_vert_tri.push(n1pcart[1]); block_vert_tri.push(n1pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(n1mcart[0]); block_vert_tri.push(n1mcart[1]); block_vert_tri.push(n1mcart[2]);
                block_vert_tri.push(c2mcart[0]); block_vert_tri.push(c2mcart[1]); block_vert_tri.push(c2mcart[2]);
                block_vert_tri.push(c2pcart[0]); block_vert_tri.push(c2pcart[1]); block_vert_tri.push(c2pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                vec3.subtract(c2cart,midpoint,mid1);
                vec3.subtract(n3cart,midpoint,mid2);
                vec3.add(mid1,mid2);
                vec3.normalize(mid1);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c2mcart[0]); block_vert_tri.push(c2mcart[1]); block_vert_tri.push(c2mcart[2]);
                block_vert_tri.push(n3pcart[0]); block_vert_tri.push(n3pcart[1]); block_vert_tri.push(n3pcart[2]);
                block_vert_tri.push(c2pcart[0]); block_vert_tri.push(c2pcart[1]); block_vert_tri.push(c2pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c2mcart[0]); block_vert_tri.push(c2mcart[1]); block_vert_tri.push(c2mcart[2]);
                block_vert_tri.push(n3mcart[0]); block_vert_tri.push(n3mcart[1]); block_vert_tri.push(n3mcart[2]);
                block_vert_tri.push(n3pcart[0]); block_vert_tri.push(n3pcart[1]); block_vert_tri.push(n3pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                vec3.subtract(n3cart,midpoint,mid1);
                vec3.subtract(c4cart,midpoint,mid2);
                vec3.add(mid1,mid2);
                vec3.normalize(mid1);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(n3mcart[0]); block_vert_tri.push(n3mcart[1]); block_vert_tri.push(n3mcart[2]);
                block_vert_tri.push(c4pcart[0]); block_vert_tri.push(c4pcart[1]); block_vert_tri.push(c4pcart[2]);
                block_vert_tri.push(n3pcart[0]); block_vert_tri.push(n3pcart[1]); block_vert_tri.push(n3pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(n3mcart[0]); block_vert_tri.push(n3mcart[1]); block_vert_tri.push(n3mcart[2]);
                block_vert_tri.push(c4mcart[0]); block_vert_tri.push(c4mcart[1]); block_vert_tri.push(c4mcart[2]);
                block_vert_tri.push(c4pcart[0]); block_vert_tri.push(c4pcart[1]); block_vert_tri.push(c4pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                vec3.subtract(c4cart,midpoint,mid1);
                vec3.subtract(c5cart,midpoint,mid2);
                vec3.add(mid1,mid2);
                vec3.normalize(mid1);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c4mcart[0]); block_vert_tri.push(c4mcart[1]); block_vert_tri.push(c4mcart[2]);
                block_vert_tri.push(c5pcart[0]); block_vert_tri.push(c5pcart[1]); block_vert_tri.push(c5pcart[2]);
                block_vert_tri.push(c4pcart[0]); block_vert_tri.push(c4pcart[1]); block_vert_tri.push(c4pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c4mcart[0]); block_vert_tri.push(c4mcart[1]); block_vert_tri.push(c4mcart[2]);
                block_vert_tri.push(c5mcart[0]); block_vert_tri.push(c5mcart[1]); block_vert_tri.push(c5mcart[2]);
                block_vert_tri.push(c5pcart[0]); block_vert_tri.push(c5pcart[1]); block_vert_tri.push(c5pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                vec3.subtract(c5cart,midpoint,mid1);
                vec3.subtract(c6cart,midpoint,mid2);
                vec3.add(mid1,mid2);
                vec3.normalize(mid1);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c5mcart[0]); block_vert_tri.push(c5mcart[1]); block_vert_tri.push(c5mcart[2]);
                block_vert_tri.push(c6pcart[0]); block_vert_tri.push(c6pcart[1]); block_vert_tri.push(c6pcart[2]);
                block_vert_tri.push(c5pcart[0]); block_vert_tri.push(c5pcart[1]); block_vert_tri.push(c5pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c5mcart[0]); block_vert_tri.push(c5mcart[1]); block_vert_tri.push(c5mcart[2]);
                block_vert_tri.push(c6mcart[0]); block_vert_tri.push(c6mcart[1]); block_vert_tri.push(c6mcart[2]);
                block_vert_tri.push(c6pcart[0]); block_vert_tri.push(c6pcart[1]); block_vert_tri.push(c6pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                vec3.subtract(c6cart,midpoint,mid1);
                vec3.subtract(n1cart,midpoint,mid2);
                vec3.add(mid1,mid2);
                vec3.normalize(mid1);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c6mcart[0]); block_vert_tri.push(c6mcart[1]); block_vert_tri.push(c6mcart[2]);
                block_vert_tri.push(n1pcart[0]); block_vert_tri.push(n1pcart[1]); block_vert_tri.push(n1pcart[2]);
                block_vert_tri.push(c6pcart[0]); block_vert_tri.push(c6pcart[1]); block_vert_tri.push(c6pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                block_vert_tri.push(c6mcart[0]); block_vert_tri.push(c6mcart[1]); block_vert_tri.push(c6mcart[2]);
                block_vert_tri.push(n1mcart[0]); block_vert_tri.push(n1mcart[1]); block_vert_tri.push(n1mcart[2]);
                block_vert_tri.push(n1pcart[0]); block_vert_tri.push(n1pcart[1]); block_vert_tri.push(n1pcart[2]);

                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                if(n9&&c8&&n7){
                    var n9cart = vec3.create([n9.x(),n9.y(),n9.z()]);
                    var c8cart = vec3.create([c8.x(),c8.y(),c8.z()]);
                    var n7cart = vec3.create([n7.x(),n7.y(),n7.z()]);
                    vec3.add(n9cart,up,n9pcart);
                    vec3.add(c8cart,up,c8pcart);
                    vec3.add(n7cart,up,n7pcart);
                    vec3.subtract(n9cart,up,n9mcart);
                    vec3.subtract(c8cart,up,c8mcart);
                    vec3.subtract(n7cart,up,n7mcart);
                    var midpoint = vec3.create([0.0,0.0,0.0]);
                    vec3.add(midpoint,c4cart);
                    vec3.add(midpoint,c5cart);
                    vec3.add(midpoint,n7cart);
                    vec3.add(midpoint,c8cart);
                    vec3.add(midpoint,n9cart);
                    midpoint[0] /= 5.0;
                    midpoint[1] /= 5.0;
                    midpoint[2] /= 5.0;
                    vec3.add(midpoint,up,midpointp);
                    vec3.subtract(midpoint,up,midpointm);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c5pcart[0]); block_vert_tri.push(c5pcart[1]); block_vert_tri.push(c5pcart[2]);
                    block_vert_tri.push(c4pcart[0]); block_vert_tri.push(c4pcart[1]); block_vert_tri.push(c4pcart[2]);
                    block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c4pcart[0]); block_vert_tri.push(c4pcart[1]); block_vert_tri.push(c4pcart[2]);
                    block_vert_tri.push(n9pcart[0]); block_vert_tri.push(n9pcart[1]); block_vert_tri.push(n9pcart[2]);
                    block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(n9pcart[0]); block_vert_tri.push(n9pcart[1]); block_vert_tri.push(n9pcart[2]);
                    block_vert_tri.push(c8pcart[0]); block_vert_tri.push(c8pcart[1]); block_vert_tri.push(c8pcart[2]);
                    block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c8pcart[0]); block_vert_tri.push(c8pcart[1]); block_vert_tri.push(c8pcart[2]);
                    block_vert_tri.push(n7pcart[0]); block_vert_tri.push(n7pcart[1]); block_vert_tri.push(n7pcart[2]);
                    block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(n7pcart[0]); block_vert_tri.push(n7pcart[1]); block_vert_tri.push(n7pcart[2]);
                    block_vert_tri.push(c5pcart[0]); block_vert_tri.push(c5pcart[1]); block_vert_tri.push(c5pcart[2]);
                    block_vert_tri.push(midpointp[0]); block_vert_tri.push(midpointp[1]); block_vert_tri.push(midpointp[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c5mcart[0]); block_vert_tri.push(c5mcart[1]); block_vert_tri.push(c5mcart[2]);
                    block_vert_tri.push(c4mcart[0]); block_vert_tri.push(c4mcart[1]); block_vert_tri.push(c4mcart[2]);
                    block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c4mcart[0]); block_vert_tri.push(c4mcart[1]); block_vert_tri.push(c4mcart[2]);
                    block_vert_tri.push(n9mcart[0]); block_vert_tri.push(n9mcart[1]); block_vert_tri.push(n9mcart[2]);
                    block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(n9mcart[0]); block_vert_tri.push(n9mcart[1]); block_vert_tri.push(n9mcart[2]);
                    block_vert_tri.push(c8mcart[0]); block_vert_tri.push(c8mcart[1]); block_vert_tri.push(c8mcart[2]);
                    block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c8mcart[0]); block_vert_tri.push(c8mcart[1]); block_vert_tri.push(c8mcart[2]);
                    block_vert_tri.push(n7mcart[0]); block_vert_tri.push(n7mcart[1]); block_vert_tri.push(n7mcart[2]);
                    block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(n7mcart[0]); block_vert_tri.push(n7mcart[1]); block_vert_tri.push(n7mcart[2]);
                    block_vert_tri.push(c5mcart[0]); block_vert_tri.push(c5mcart[1]); block_vert_tri.push(c5mcart[2]);
                    block_vert_tri.push(midpointm[0]); block_vert_tri.push(midpointm[1]); block_vert_tri.push(midpointm[2]);

                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);
                    block_norm_tri.push(up[0]); block_norm_tri.push(up[1]); block_norm_tri.push(up[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    // The edges
                    vec3.subtract(c4cart,midpoint,mid1);
                    vec3.subtract(n9cart,midpoint,mid2);
                    vec3.add(mid1,mid2);
                    vec3.normalize(mid1);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c4mcart[0]); block_vert_tri.push(c4mcart[1]); block_vert_tri.push(c4mcart[2]);
                    block_vert_tri.push(n9pcart[0]); block_vert_tri.push(n9pcart[1]); block_vert_tri.push(n9pcart[2]);
                    block_vert_tri.push(c4pcart[0]); block_vert_tri.push(c4pcart[1]); block_vert_tri.push(c4pcart[2]);

                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c4mcart[0]); block_vert_tri.push(c4mcart[1]); block_vert_tri.push(c4mcart[2]);
                    block_vert_tri.push(n9mcart[0]); block_vert_tri.push(n9mcart[1]); block_vert_tri.push(n9mcart[2]);
                    block_vert_tri.push(n9pcart[0]); block_vert_tri.push(n9pcart[1]); block_vert_tri.push(n9pcart[2]);

                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);
                    
                    vec3.subtract(n9cart,midpoint,mid1);
                    vec3.subtract(c8cart,midpoint,mid2);
                    vec3.add(mid1,mid2);
                    vec3.normalize(mid1);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(n9mcart[0]); block_vert_tri.push(n9mcart[1]); block_vert_tri.push(n9mcart[2]);
                    block_vert_tri.push(c8pcart[0]); block_vert_tri.push(c8pcart[1]); block_vert_tri.push(c8pcart[2]);
                    block_vert_tri.push(n9pcart[0]); block_vert_tri.push(n9pcart[1]); block_vert_tri.push(n9pcart[2]);

                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(n9mcart[0]); block_vert_tri.push(n9mcart[1]); block_vert_tri.push(n9mcart[2]);
                    block_vert_tri.push(c8mcart[0]); block_vert_tri.push(c8mcart[1]); block_vert_tri.push(c8mcart[2]);
                    block_vert_tri.push(c8pcart[0]); block_vert_tri.push(c8pcart[1]); block_vert_tri.push(c8pcart[2]);

                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    vec3.subtract(c8cart,midpoint,mid1);
                    vec3.subtract(n7cart,midpoint,mid2);
                    vec3.add(mid1,mid2);
                    vec3.normalize(mid1);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c8mcart[0]); block_vert_tri.push(c8mcart[1]); block_vert_tri.push(c8mcart[2]);
                    block_vert_tri.push(n7pcart[0]); block_vert_tri.push(n7pcart[1]); block_vert_tri.push(n7pcart[2]);
                    block_vert_tri.push(c8pcart[0]); block_vert_tri.push(c8pcart[1]); block_vert_tri.push(c8pcart[2]);

                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(c8mcart[0]); block_vert_tri.push(c8mcart[1]); block_vert_tri.push(c8mcart[2]);
                    block_vert_tri.push(n7mcart[0]); block_vert_tri.push(n7mcart[1]); block_vert_tri.push(n7mcart[2]);
                    block_vert_tri.push(n7pcart[0]); block_vert_tri.push(n7pcart[1]); block_vert_tri.push(n7pcart[2]);

                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    vec3.subtract(n7cart,midpoint,mid1);
                    vec3.subtract(c5cart,midpoint,mid2);
                    vec3.add(mid1,mid2);
                    vec3.normalize(mid1);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(n7mcart[0]); block_vert_tri.push(n7mcart[1]); block_vert_tri.push(n7mcart[2]);
                    block_vert_tri.push(c5pcart[0]); block_vert_tri.push(c5pcart[1]); block_vert_tri.push(c5pcart[2]);
                    block_vert_tri.push(n7pcart[0]); block_vert_tri.push(n7pcart[1]); block_vert_tri.push(n7pcart[2]);

                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);

                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][0]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][1]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][2]);
                    block_col_tri.push(colourScheme[atoms[iat]["_atom_site.id"]][3]);

                    block_vert_tri.push(n7mcart[0]); block_vert_tri.push(n7mcart[1]); block_vert_tri.push(n7mcart[2]);
                    block_vert_tri.push(c5mcart[0]); block_vert_tri.push(c5mcart[1]); block_vert_tri.push(c5mcart[2]);
                    block_vert_tri.push(c5pcart[0]); block_vert_tri.push(c5pcart[1]); block_vert_tri.push(c5pcart[2]);

                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);
                    block_norm_tri.push(mid1[0]); block_norm_tri.push(mid1[1]); block_norm_tri.push(mid1[2]);

                    block_idx_tri.push(idx++); block_idx_tri.push(idx++); block_idx_tri.push(idx++);
                    
                }
            }
        }
    }

    var blockPrimitiveInfo = {"atoms":[[atoms]], "col_tri":[[block_col_tri]], "norm_tri":[[block_norm_tri]], "vert_tri":[[block_vert_tri]], "idx_tri":[[block_idx_tri]] , "prim_types":[["TRIANGLES"]] };
    return blockPrimitiveInfo;
}

function getGlycoBlocks(model,size,colourScheme,filterAtoms){

    //FIXME - This is the one that should take a selection or set of atoms.
    //      - The filtering is totally unneccessary if filterAtoms == all or filterAtoms == nglyc. Hmm.
    var objects = [];

    // Glycoblocks prototyping .... This belongs elsewhere.
    var atomColoursResType = colourScheme.colourByResidueType();
    var blackColours = colourScheme.colourOneColour([0.0,0.0,0.0,1.0]);

    var glycanAtoms = model.getAtoms("nglycosylation");
    if(typeof(filterAtoms)!=="undefined"){
        //FIXME - This probably deserves a method.
        var atoms1 = glycanAtoms;
        var atoms2 = filterAtoms;
        var andAtoms = [];
        var allAtoms = model.getAllAtoms();
        var allLen = allAtoms.length;
        var len1 = atoms1.length;
        var len2 = atoms2.length;
        var uuid = guid();
        var uuid2 = guid();
        for(var iat=0;iat<allLen;iat++){
            allAtoms[iat][uuid] = 0;
            allAtoms[iat][uuid2] = 0;
        }
        for(var iat=0;iat<len1;iat++){
            atoms1[iat][uuid] = 1;
        }
        for(var iat=0;iat<len2;iat++){
            atoms2[iat][uuid2] = 1;
        }
        andAtoms = allAtoms.filter(function(x) { return x[uuid] === 1 && x[uuid2] === 1 });
        for(var iat=0;iat<allLen;iat++){
            delete allAtoms[iat][uuid];
            delete allAtoms[iat][uuid2];
        }
        glycanAtoms = andAtoms;
        
    }

    //var glycanSpheres = atomsToSpheresInfo(glycanAtoms,0.4,atomColours);
    //objects.push(glycanSpheres);

    var glycanResidues = model.getGlycanResidues();
    if(typeof(filterAtoms)!=="undefined"){
        glycanResidues = model.filterResiduesByAtoms(glycanResidues,filterAtoms);
    }
    var glycoBlocks = residuesToGlycoBlocksInfo(glycanResidues,size,atomColoursResType);
    objects.push(glycoBlocks);

    var cylinder_sizes = [];
    var cylinder_col_tri = [];
    var cylinder_vert_tri = [];
    var cylinder_idx_tri = [];
    var sphere_atoms = [];
    var glycConn = model.glycan_cache["glycanGlycanConnections"];
    if(typeof(filterAtoms)!=="undefined"){
        glycConn = model.filterGlycanConnectionsByAtoms(glycConn,filterAtoms);
    }
    console.log(glycConn);
    var gg = vec3.create();
    var uddAtoms = ["GLYCO_BLOCK_C1","GLYCO_BLOCK_C2","GLYCO_BLOCK_C3","GLYCO_BLOCK_C4","GLYCO_BLOCK_C5"];
    for(var ig=0;ig<glycConn.length;ig++){
        var mindist = 1e+8;
        var c1 = null;
        var c2 = null;
        var theT1 = null;
        var theT2 = null;
        for(var it1=0;it1<uddAtoms.length;it1++){
            var t1 = uddAtoms[it1];
            for(var it2=0;it2<uddAtoms.length;it2++){
                var t2 = uddAtoms[it2];
                if(typeof(glycConn[ig][0][t1])!=="undefined"&&typeof(glycConn[ig][1][t2])!=="undefined"){
                    var gg1 = vec3.create(glycConn[ig][0][t1]);
                    var gg2 = vec3.create(glycConn[ig][1][t2]);
                    vec3.subtract(gg1,gg2,gg);
                    var dist = vec3.length(gg);
                    if(dist<mindist){
                        mindist = dist;
                        c1 = glycConn[ig][0][t1];
                        c2 = glycConn[ig][1][t2];
                        theT1 = t1;
                        theT2 = t2;
                    }
                }
            }
        }
        cylinder_vert_tri.push(c1[0]);
        cylinder_vert_tri.push(c1[1]);
        cylinder_vert_tri.push(c1[2]);
        cylinder_vert_tri.push(c2[0]);
        cylinder_vert_tri.push(c2[1]);
        cylinder_vert_tri.push(c2[2]);
        cylinder_col_tri.push(0.0);
        cylinder_col_tri.push(0.0);
        cylinder_col_tri.push(0.0);
        cylinder_col_tri.push(1.0);
        cylinder_col_tri.push(0.0);
        cylinder_col_tri.push(0.0);
        cylinder_col_tri.push(0.0);
        cylinder_col_tri.push(1.0);
        cylinder_sizes.push(0.15);
        cylinder_idx_tri.push(2*ig);

        // FIXME - now make dummy atoms here and do contacts to cylinders.
        /*
        var atom = {};
        atom["x"] = theT1[0];
        atom["y"] = theT1[1];
        atom["z"] = theT1[2];
        atom["tempFactor"] = at1["_atom_site.B_iso_or_equiv"];
        atom["charge"] = at1["_atom_site.pdbx_formal_charge"];
        atom["label"] =  at1.getAtomID();
        sphere_atoms.push(atom);
        */
    }

    var cylinderPrimitiveInfo = {"atoms":[[sphere_atoms]],"sizes": [[cylinder_sizes]], "col_tri":[[cylinder_col_tri]], "norm_tri":[[[]]], "vert_tri":[[cylinder_vert_tri]], "idx_tri":[[cylinder_idx_tri]] , "prim_types":[["CYLINDERS"]] };
    console.log(cylinderPrimitiveInfo);
    objects.push(cylinderPrimitiveInfo);

    var rootGlycConn = model.glycan_cache["rootGlycans"];
    if(typeof(filterAtoms)!=="undefined"){
        rootGlycConn = model.filterConnectionsByAtoms2(rootGlycConn,filterAtoms);
    }
    var rootGlycConnPrimitiveInfo = contactsToCylindersInfo(rootGlycConn,size,atomColoursResType);
    objects.push(rootGlycConnPrimitiveInfo);

    var peptideAtoms = model.getAtoms("peptide");
    var saccharideAtoms = glycanAtoms;
    var glycoPeptideHBonds = model.getHBonds(peptideAtoms,saccharideAtoms);
    for(var ig=0;ig<glycoPeptideHBonds.length;ig++){
        if(glycoPeptideHBonds[ig][1].residue.getAtom("CA")){
            glycoPeptideHBonds[ig][1] =  glycoPeptideHBonds[ig][1].residue.getAtom("CA");
            var atom = new Atom(glycoPeptideHBonds[ig][2]);
            atom["_atom_site.Cartn_x"] = glycoPeptideHBonds[ig][2].residue["GLYCO_BLOCK_CENTRE"][0];
            atom["_atom_site.Cartn_y"] = glycoPeptideHBonds[ig][2].residue["GLYCO_BLOCK_CENTRE"][1];
            atom["_atom_site.Cartn_z"] = glycoPeptideHBonds[ig][2].residue["GLYCO_BLOCK_CENTRE"][2];
            glycoPeptideHBonds[ig][2] = atom;
        } else if(glycoPeptideHBonds[ig][2].residue.getAtom("CA")) {
            glycoPeptideHBonds[ig][2] =  glycoPeptideHBonds[ig][2].residue.getAtom("CA");
            var atom = new Atom(glycoPeptideHBonds[ig][1]);
            atom["_atom_site.Cartn_x"] = glycoPeptideHBonds[ig][1].residue["GLYCO_BLOCK_CENTRE"][0];
            atom["_atom_site.Cartn_y"] = glycoPeptideHBonds[ig][1].residue["GLYCO_BLOCK_CENTRE"][1];
            atom["_atom_site.Cartn_z"] = glycoPeptideHBonds[ig][1].residue["GLYCO_BLOCK_CENTRE"][2];
            glycoPeptideHBonds[ig][1] = atom;
        } else {
            console.log("Confusion in glycan-peptide hbonds ????????????????????????????????");
        }
    }

    var glycoGlycoHBonds = model.getHBonds(saccharideAtoms,saccharideAtoms);
    for(var ig=0;ig<glycoGlycoHBonds.length;ig++){
        var atom = new Atom(glycoGlycoHBonds[ig][2]);
        atom["_atom_site.Cartn_x"] = glycoGlycoHBonds[ig][2].residue["GLYCO_BLOCK_CENTRE"][0];
        atom["_atom_site.Cartn_y"] = glycoGlycoHBonds[ig][2].residue["GLYCO_BLOCK_CENTRE"][1];
        atom["_atom_site.Cartn_z"] = glycoGlycoHBonds[ig][2].residue["GLYCO_BLOCK_CENTRE"][2];
        glycoGlycoHBonds[ig][2] = atom;
        var atom2 = new Atom(glycoGlycoHBonds[ig][1]);
        atom2["_atom_site.Cartn_x"] = glycoGlycoHBonds[ig][1].residue["GLYCO_BLOCK_CENTRE"][0];
        atom2["_atom_site.Cartn_y"] = glycoGlycoHBonds[ig][1].residue["GLYCO_BLOCK_CENTRE"][1];
        atom2["_atom_site.Cartn_z"] = glycoGlycoHBonds[ig][1].residue["GLYCO_BLOCK_CENTRE"][2];
        glycoGlycoHBonds[ig][1] = atom2;
    }

    var glycoPeptideHBondPrimitiveInfo = contactsToCappedCylindersInfo(glycoPeptideHBonds,size*.5,blackColours,true);
    var glycoGlycoHBondPrimitiveInfo = contactsToCappedCylindersInfo(glycoGlycoHBonds,size*.5,blackColours,true);

    objects.push(glycoPeptideHBondPrimitiveInfo);
    objects.push(glycoGlycoHBondPrimitiveInfo);

    return objects;

}
