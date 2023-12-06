function Cell(){
}

Cell.prototype.init = function(a,b,c,alpha,beta,gamma){
    // store cell (allows a/b/c/alpha/beta/gamma methods)
    this.a = a;
    this.b = b;
    this.c = c;
    this.alpha = alpha*Math.PI/180.;
    this.beta  = beta*Math.PI/180.;
    this.gamma = gamma*Math.PI/180.;

    // calculate cell volume (allows a_star../alpha_star.. methods)
    this.vol = this.a * this.b * this.c * Math.sqrt( 2.0*Math.cos(this.alpha)*Math.cos(this.beta)*Math.cos(this.gamma)
            - Math.cos(this.alpha)*Math.cos(this.alpha)
            - Math.cos( this.beta)*Math.cos( this.beta)
            - Math.cos(this.gamma)*Math.cos(this.gamma) + 1.0 );

    this.a_star     = b*c*Math.sin(this.alpha)/this.vol;
    this.b_star     = c*a*Math.sin(this.beta)/this.vol; 
    this.c_star     = a*b*Math.sin(this.gamma)/this.vol; 
    this.alpha_star = Math.acos( (Math.cos(this.gamma)*Math.cos(this.beta)-Math.cos(this.alpha)) / (Math.sin(this.beta)*Math.sin(this.gamma)) ); 
    this.beta_star  = Math.acos( (Math.cos(this.alpha)*Math.cos(this.gamma)-Math.cos(this.beta)) / (Math.sin(this.gamma)*Math.sin(this.alpha)) ); 
    this.gamma_star = Math.acos( (Math.cos(this.beta)*Math.cos(this.alpha)-Math.cos(this.gamma)) / (Math.sin(this.alpha)*Math.sin(this.beta)) ); 

    // deal with null cell
    if(this.vol < 1e-6){
        return;
    }

    // othogonalisation+fractionisation matrices
    var orthmat = mat4.create();
    orthmat[0*4+0] =  this.a;
    orthmat[0*4+1] =  this.b*Math.cos(this.gamma);
    orthmat[0*4+2] =  this.c*Math.cos( this.beta);
    orthmat[1*4+1] =  this.b*Math.sin(this.gamma);
    orthmat[1*4+2] = -this.c*Math.sin( this.beta)*Math.cos(this.alpha_star);
    orthmat[2*4+2] =  this.c*Math.sin( this.beta)*Math.sin(this.alpha_star);
    orthmat[3*4+3] =  1.0;
    var fracmat = mat4.create();
    mat4.inverse(orthmat,fracmat);

    this.fracmat = mat4.create();
    this.orthmat = mat4.create();

    mat4.transpose(orthmat,this.orthmat);
    mat4.transpose(fracmat,this.fracmat);

    this.matrix_frac = this.fracmat;
    this.matrix_orth = this.orthmat;
}

function readMapFromArrayBuffer(contents) {
    var intArray = new Int32Array(contents);
    var floatArray = new Float32Array(contents);

    var header = {};

    header["machst"] = intArray[54];// Machine stamp
    // if machst == 1055298683 little endian? If not, then we have to create arrays "differently"
    // https://www.html5rocks.com/en/tutorials/webgl/typed_arrays/

    header["nc"] = intArray[0];
    header["nr"] = intArray[1];
    header["ns"] = intArray[2];
    header["mode"] = intArray[3];
    header["ncstart"] = intArray[4];
    header["nrstart"] = intArray[5];
    header["nsstart"] = intArray[6];
    header["nx"] = intArray[7];
    header["ny"] = intArray[8];
    header["nz"] = intArray[9];
    header["xlen"] = floatArray[10];
    header["ylen"] = floatArray[11];
    header["zlen"] = floatArray[12];
    header["alpha"] = floatArray[13];
    header["beta"] = floatArray[14];
    header["gamma"] = floatArray[15];
    header["mapc"] = intArray[16];
    header["mapr"] = intArray[17];
    header["maps"] = intArray[18];
    header["amin"] = floatArray[19];
    header["amax"] = floatArray[20];
    header["amean"] = floatArray[21];
    header["ispg"] = intArray[22];
    header["nsymbt"] = intArray[23] / 4;
    header["lskflg"] = intArray[24];
    header["skwmat"] = [];
    if(header["lskflg"]!==0){
        header["skwmat"][0] = intArray[25];
        header["skwmat"][1] = intArray[26];
        header["skwmat"][2] = intArray[27];
        header["skwmat"][3] = intArray[28];
        header["skwmat"][4] = intArray[29];
        header["skwmat"][5] = intArray[30];
        header["skwmat"][6] = intArray[31];
        header["skwmat"][7] = intArray[32];
        header["skwmat"][8] = intArray[33];
    }

    console.log(header);
    console.log(floatArray[256+header["nsymbt"]]);// First element of map.
    console.log(floatArray.length);

    var origin = [0.0,0.0,0.0];

    var gfms0 = origin;
    var gfms1 = [-1,-1,-1];
    var orderxyz = [-1,-1,-1];
    var orderfms = [header["mapc"],header["mapr"],header["maps"]];
    var dim = [header["nc"],header["nr"],header["ns"]];
    console.log(orderfms);

    for(var i=0;i<3;i++){
        gfms1[i] = gfms0[i] + dim[i] - 1;
    }
    for(var i=0;i<3;i++){
        orderxyz[orderfms[i]-1] = i
    }
    console.log(gfms1);
    console.log(orderxyz);

    var n0 = (gfms1[0]-gfms0[0]+1);
    var n1 = n0 * (gfms1[1]-gfms0[1]+1);
    console.log(n0);
    console.log(n1);


    var xmap_0 = [header["nc"],header["nr"],header["ns"]][orderxyz[0]];
    var xmap_1 = [header["nc"],header["nr"],header["ns"]][orderxyz[1]];
    var xmap_2 = [header["nc"],header["nr"],header["ns"]][orderxyz[2]];

    console.log(xmap_0+" "+xmap_1+" "+xmap_2);

    var xmap = [];

    var g = [-1,-1,-1];
    g[2] =  gfms0[2];

    var index = 256+header["nsymbt"];
    var secSkip = header["nc"]* header["nr"];

    while(g[2] <= gfms1[2]){
        g[1] =  gfms0[1];
        while(g[1] <= gfms1[1]){
            g[0] =  g[0] = gfms0[0];
            while(g[0] <= gfms1[0]){
                if(typeof(xmap[parseInt(g[orderxyz[0]])])==="undefined"){
                    xmap[parseInt(g[orderxyz[0]])] = [];
                }
                if(typeof(xmap[parseInt(g[orderxyz[0]])][parseInt(g[orderxyz[1]])])==="undefined"){
                    xmap[parseInt(g[orderxyz[0]])][parseInt(g[orderxyz[1]])] = [];
                }
                xmap[parseInt(g[orderxyz[0]])][parseInt(g[orderxyz[1]])][parseInt(g[orderxyz[2]])] = floatArray[index];
                g[0] += 1;
                index += 1;
            }
            g[1] += 1;
        }
        g[2] += 1;
        //index += secSkip;
    }

    console.log(xmap);

    var symops = getSymOpsFromSpGrpNo(header["ispg"]);
    console.log(symops);

    var symmats = [];
    var invmats = [];
    for(var i=0;i<symops.length;i++){
        var mat = getMatrixFromSymOp(symops[i]);
        var invMat = mat4.create();
        mat4.inverse(mat,invMat);
        symmats.push(mat);
        invmats.push(invMat);
    }
    console.log(symmats);
    console.log(invmats);

    var cell = new Cell();
    cell.init(header["xlen"],header["ylen"],header["zlen"],header["alpha"],header["beta"],header["gamma"])

    var invTMats = [];
    var RF = cell.matrix_frac;
    var RO = cell.matrix_orth;

    function printMat(mat){
        for(var i=0;i<4;i++){
            console.log(mat[i*4]+" "+mat[i*4+1]+" "+mat[i*4+2]+" "+mat[i*4+3]);
        }
        console.log(" ");
    }

    //console.log("RF");
    //printMat(RF);
    for(var i=0;i<symops.length;i++){
        var tm = symmats[i];
        var tmt = mat4.create();
        mat4.transpose(tm,tmt);
        //console.log("tm");
        //printMat(tm);
        var fm = mat4.create();
        var theTMatrix = mat4.create();
        var invTMat = mat4.create();
        mat4.multiply(tmt,RF,fm);
        //console.log("fm");
        //printMat(fm);
        mat4.multiply(RO,fm,theTMatrix);
        //console.log("theTMatrix");
        //printMat(theTMatrix);
        mat4.inverse(theTMatrix,invTMat);
        //console.log("invTMat");
        //printMat(invTMat);
        invTMats.push(invTMat);
        //console.log(" ");
    }

    //console.log(invTMats);

    var ptScalarField = [];
    // And now construct a complete unit cell.
    // FIXME - goes 
    var veco = mat4.create();
    for(var i=0;i<header["nz"];i++){
        for(var j=0;j<header["ny"];j++){
            for(var k=0;k<header["nx"];k++){
                if(i<xmap_2 && j<xmap_1 && k<xmap_0){
                    ptScalarField.push(xmap[k][j][i]);
                    //console.log(k+" "+j+" "+i+" -- "+k+" "+j+" "+i);
                }else{
                    var vec = [ 1.0*k/header["nx"],1.0*j/header["ny"],1.0*i/header["nz"] ];
                    mat4.multiplyVec3(RO,vec,veco);
                    //console.log(vec);
                    //console.log(veco);
                    var doneIJK = false;
                    for(var imat=1;imat<invTMats.length;imat++){
                        var tvec = vec3.create();
                        mat4.multiplyVec3(invTMats[imat],veco,tvec);
                        //console.log(tvec);
                        //console.log(invTMats[imat]);
                        tvec[0] += invTMats[imat][3];
                        tvec[1] += invTMats[imat][7];
                        tvec[2] += invTMats[imat][11];
                        //console.log(tvec);
                        mat4.multiplyVec3(RF,tvec);
                        //console.log(tvec);
                        var ii = parseInt(Math.round(tvec[2]*header["nz"]));
                        var jj = parseInt(Math.round(tvec[1]*header["ny"]));
                        var kk = parseInt(Math.round(tvec[0]*header["nx"]));
                        if(ii<0){
                            ii += header["nz"];
                        }
                        if(jj<0){
                            jj += header["ny"];
                        }
                        if(kk<0){
                            kk += header["nx"];
                        }
                        if(ii>=0 && ii<xmap_2 && jj>=0 && jj<xmap_1 && kk>=0 && kk<xmap_0){
                            ptScalarField.push(xmap[kk][jj][ii]);
                            //console.log(k+" "+j+" "+i+" -- "+kk+" "+jj+" "+ii);
                            doneIJK = true;
                            break;
                        }
                    }
                    if(!doneIJK){
                        console.log("Skipped "+k+" "+j+" "+i);
                        console.log("Giving up ...");
                        return {};
                    }
                }
            }
        }
    }
    var ret = {};
    ret["header"] = header;
    ret["ptScalarField"] = ptScalarField;
    cell.nx = header["nx"];
    cell.ny = header["ny"];
    cell.nz = header["nz"];
    ret["cell"] = cell;
    return ret;
}

function mapToMapGrid(map){
    var header = map["header"];
    var mapo = new CIsoSurface();
    var xCellLen = header["xlen"]/header["nx"];
    var yCellLen = header["ylen"]/header["ny"];
    var zCellLen = header["zlen"]/header["nz"];
    var mapGrid = {};
    mapGrid.nCells_x = header["nx"]-1;
    mapGrid.nCells_y = header["ny"]-1;
    mapGrid.nCells_z = header["nz"]-1;
    mapGrid.cellLength_x = 1.0*header["xlen"]/(header["nx"]-1);
    mapGrid.cellLength_y = 1.0*header["ylen"]/(header["ny"]-1);
    mapGrid.cellLength_z = 1.0*header["zlen"]/(header["nz"]-1);
    mapGrid.grid = [map["ptScalarField"]];
    mapGrid.cell = map["cell"];
    mapGrid.fracToOrthMat = [header["xlen"], 0.0, 0.0, 0.0, header["ylen"], 0.0, 0.0, 0.0, header["zlen"]];
    var tmpMat = mat3.toMat4(mapGrid.fracToOrthMat);
    // Now I think fracMat should be multiplied by skwmat and then we invert for orthMat 
    mapGrid.orthToFracMat = mat4.toInverseMat3(tmpMat);
    return mapGrid;
}
