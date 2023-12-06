function wizardBonds(pdbatoms,enerLib,bruteForceHB){
    var hier = pdbatoms["atoms"];

    var colourScheme = new ColourScheme(pdbatoms);
    var atomColours = colourScheme.colourByAtomType();

    var objects = [];

    for(var imod=0;imod<1;imod++){
        var model = hier[imod];
        /*
        var atoms = model.getAllAtoms();
        contacts = model.SeekContacts(atoms,atoms,0.6,1.6);
        */
        var contactsAndSingletons = model.getBondsContactsAndSingletons();
        var contacts = contactsAndSingletons["contacts"];
        var singletons = contactsAndSingletons["singletons"];
        var linePrimitiveInfo = contactsToLinesInfo(contacts,4,atomColours);
        var singletonPrimitiveInfo = singletonsToLinesInfo(singletons,4,atomColours);
        objects.push(linePrimitiveInfo);
        objects.push(singletonPrimitiveInfo);
    }

    return objects;
}

function wizardRibbons(pdbatoms,enerLib,bruteForceHB){

    var start = new Date().getTime();
    var objects = [];

    var hier = pdbatoms["atoms"];
    var flagBulge = true;
    CalcSecStructure(hier,flagBulge);

    console.log("Time to CalcSecStructure: "+(new Date().getTime()-start));

    var colourScheme = new ColourScheme(pdbatoms);
    var atomColours = colourScheme.colourBySecondaryStructure({"strand":[0.0,0.0,1.0,1.0],"helix":[1.0,0.0,0.0,1.0]});
    var coloured_splines_info = GetSplinesColoured(pdbatoms,atomColours);

    console.log("Time to GetSplinesColoured: "+(new Date().getTime()-start));

    for(itri=0;itri<coloured_splines_info.length;itri++){
        objects.push(coloured_splines_info[itri]);
    }

    console.log("Time to do everything: "+(new Date().getTime()-start));

    return objects;
  
}

function wizardRibbonsByChain(pdbatoms,enerLib,bruteForceHB){

    var start = new Date().getTime();
    var objects = [];

    var hier = pdbatoms["atoms"];
    var flagBulge = true;
    CalcSecStructure(hier,flagBulge);

    console.log("Time to CalcSecStructure: "+(new Date().getTime()-start));

    var colourScheme = new ColourScheme(pdbatoms);
    var atomColours = colourScheme.colourByChain();
    var coloured_splines_info = GetSplinesColoured(pdbatoms,atomColours);

    console.log("Time to GetSplinesColoured: "+(new Date().getTime()-start));

    for(itri=0;itri<coloured_splines_info.length;itri++){
        objects.push(coloured_splines_info[itri]);
    }

    console.log("Time to do everything: "+(new Date().getTime()-start));

    return objects;
  
}

function wizardWorms(pdbatoms,enerLib,bruteForceHB){

    var start = new Date().getTime();
    var objects = [];

    var hier = pdbatoms["atoms"];

    var colourScheme = new ColourScheme(pdbatoms);
    var atomColours = colourScheme.colourByChain();
    var coloured_splines_info = GetWormColoured(pdbatoms,atomColours);

    console.log("Time to GetSplinesColoured: "+(new Date().getTime()-start));

    for(itri=0;itri<coloured_splines_info.length;itri++){
        objects.push(coloured_splines_info[itri]);
    }

    console.log("Time to do everything: "+(new Date().getTime()-start));

    return objects;
  
}

function wizardSiteAndRibbonsByChain(pdbatoms,enerLib,bruteForceHB){

    // TODO 
    //  - Glycoblocks (NEED glycoblocks and glycosylation selection!)
    //  - ligand, neighb, and all-atom surfaces (invis). (NEED surface!)

    var start = new Date().getTime();
    var objects = [];

    if(typeof(pdbatoms["modamino"])!=="undefined"){
        for(imod=0;imod<pdbatoms["modamino"].length;imod++){
            Model.prototype.getPeptideLibraryEntry(pdbatoms["modamino"][imod],enerLib);
        }
    }

    enerLib.AssignHBTypes(pdbatoms,bruteForceHB);
    console.log("Time to AssignHBTypes: "+(new Date().getTime()-start));
    var hier = pdbatoms["atoms"];
    var model = hier[0];
    console.log("Time to 'get model': "+(new Date().getTime()-start));
    model.calculateHBonds();
    //return objects;
    console.log("Time to calculateHBonds: "+(new Date().getTime()-start));

    var flagBulge = true;
    CalcSecStructure(hier,flagBulge);
    console.log("Time to CalcSecStructure: "+(new Date().getTime()-start));

    var colourScheme = new ColourScheme(pdbatoms);
    var atomColours = colourScheme.colourByChain({"nonCByAtomType":true});
    var atomColoursByChainWithC = colourScheme.colourByChain({"nonCByAtomType":false});
    var coloured_splines_info = GetSplinesColoured(pdbatoms,atomColours);
    console.log("Time to GetSplinesColoured: "+(new Date().getTime()-start));
    for(itri=0;itri<coloured_splines_info.length;itri++){
        if(typeof(coloured_splines_info[itri].sizes)!=="undefined"&&coloured_splines_info[itri].sizes.length>0&&coloured_splines_info[itri].sizes[0].length>0&&coloured_splines_info[itri].sizes[0][0].length>0){
            objects.push(coloured_splines_info[itri]);
        }
    }

    var atomColours2 = colourScheme.colourByAtomType();

    // The hidden all atom model, water and solute.
    for(var imod=0;imod<1;imod++){
        var model = hier[imod];
        atoms = model.getAllAtoms();
        var contacts = model.SeekContacts(atoms,atoms,0.6,1.6);
        var linePrimitiveInfo = contactsToLinesInfo(contacts,1,atomColours2);
        linePrimitiveInfo.visibility = [false];
        objects.push(linePrimitiveInfo);
    }
    console.log("Time to all atom objects: "+(new Date().getTime()-start));
    var allWater = model.getAtoms("water");
    var allWaterSpheres = atomsToSpheresInfo(allWater,0.2,atomColours2);
    allWaterSpheres.visibility = [false];
    objects.push(allWaterSpheres);
    console.log("Time to water objects: "+(new Date().getTime()-start));

    var solute = model.getAtoms("solute");
    var soluteSpheres = atomsToSpheresInfo(solute,0.2,atomColours2);
    soluteSpheres.visibility = [false];
    objects.push(soluteSpheres);
    var soluteContacts = model.SeekContacts(solute,solute,0.6,1.6);
    var solutePrimitiveInfo = contactsToCylindersInfo(soluteContacts,0.2,atomColours2,false);
    solutePrimitiveInfo.visibility = [false];
    objects.push(solutePrimitiveInfo);
    console.log("Time to solute objects: "+(new Date().getTime()-start));

    var ligandAtoms = model.getAtoms("ligands");
    console.log("Time to get ligands atoms: "+(new Date().getTime()-start));
    var multipleBonds = getMultipleBonds(ligandAtoms,enerLib,0.15,atomColours2);
    console.log("Time to get ligands multiple bonds: "+(new Date().getTime()-start));
    for(var imbo=0;imbo<multipleBonds.length;imbo++){
        objects.push(multipleBonds[imbo]);
    }
    console.log("Time to get ligands bonds objects: "+(new Date().getTime()-start));
    var spheres = atomsToSpheresInfo(ligandAtoms,0.4,atomColours2);
    console.log("Time to get ligands spheres: "+(new Date().getTime()-start));
    objects.push(spheres);
    console.log("Time to get ligands objects: "+(new Date().getTime()-start));

    var allAtoms = model.getAtoms("all");
    var hbonds = model.getHBonds(allAtoms,ligandAtoms);
    console.log("HBonds: "+hbonds.length);
    var blackColours = colourScheme.colourOneColour([0.0,0.0,0.0,1.0]);
    var hBondPrimitiveInfo = contactsToCappedCylindersInfo(hbonds,0.05,blackColours,true);
    objects.push(hBondPrimitiveInfo);
    console.log("Time to get hbond objects: "+(new Date().getTime()-start));

    var sideChains = model.getAtoms("neighb cid=\"ligands\" maxd=4.0 mind=0.0 group=side excl=water,ligands");
    console.log("Time to get side chain atoms: "+(new Date().getTime()-start));
    var mainChains = model.getAtoms("neighb cid=\"ligands\" maxd=4.0 mind=0.0 group=main excl=water,ligands");
    console.log("Time to get main chain atoms: "+(new Date().getTime()-start));

    var mainChainsHBonded = model.getAtoms("neighb cid=\"ligands\" maxd=4.0 group=main hbonded=1 excl=central,water");
    var mainSpheres = atomsToSpheresInfo(mainChainsHBonded,0.2,atomColours);
    objects.push(mainSpheres);
    var mainHB_contacts = model.SeekContacts(mainChainsHBonded,mainChainsHBonded,0.6,1.6);
    var mainHBPrimitiveInfo = contactsToCylindersInfo(mainHB_contacts,0.2,atomColours);
    objects.push(mainHBPrimitiveInfo);

    var contacts = model.SeekContacts(sideChains,sideChains,0.6,1.6);
    var sideCylinderPrimitiveInfo = contactsToCylindersInfo(contacts,0.2,atomColours);
    objects.push(sideCylinderPrimitiveInfo);
    var sideSpheres = atomsToSpheresInfo(sideChains,0.2,atomColours);
    objects.push(sideSpheres);
    console.log("Time to get neighbourhood objects: "+(new Date().getTime()-start));

    var water = model.getAtoms("water and neighb cid=\"ligands\" maxd=4.0 mind=0.0 group=side excl=ligands");
    console.log("Time to get water atoms: "+(new Date().getTime()-start));
    var waterSpheres = atomsToSpheresInfo(water,0.2,atomColours2);
    objects.push(waterSpheres);

    var sideWaterHBonds = model.getHBonds(water,sideChains);
    var sideWaterHBondPrimitiveInfo = contactsToCappedCylindersInfo(sideWaterHBonds,0.05,blackColours,true);
    objects.push(sideWaterHBondPrimitiveInfo);

    /*
       // FIXME - are these required? Get a funny HBond with 8a3h.
    var mainWaterHBonds = model.getHBonds(water,mainChains);
    var mainWaterHBondPrimitiveInfo = contactsToCappedCylindersInfo(mainWaterHBonds,0.05,blackColours,true);
    objects.push(mainWaterHBondPrimitiveInfo);
    */

    var waterWaterHBonds = model.getHBonds(water,water);
    var waterWaterHBondPrimitiveInfo = contactsToCappedCylindersInfo(waterWaterHBonds,0.05,blackColours,true);
    objects.push(waterWaterHBondPrimitiveInfo);

    var baseAtoms = model.getAtoms("base");
    var base_contacts = model.SeekContacts(baseAtoms,baseAtoms,0.6,1.6);
    var basePrimitiveInfo = contactsToCylindersInfo(base_contacts,0.2,atomColoursByChainWithC);
    var baseSpheresInfo = atomsToSpheresInfo(baseAtoms,0.2,atomColoursByChainWithC);
    objects.push(basePrimitiveInfo);
    objects.push(baseSpheresInfo);

    var nucleic_atoms = model.getAtoms("nucleic");
    var blockPrimitiveInfo = getBaseBlocks(nucleic_atoms,0.2,atomColoursByChainWithC);
    objects.push(blockPrimitiveInfo);

    var metal_atoms = model.getAtoms("metal");
    var metalSpheresInfo = atomsToSpheresInfo(metal_atoms,0.6,atomColours2);
    objects.push(metalSpheresInfo);
    console.log("Time to do everything: "+(new Date().getTime()-start));

    function checkOverlapWithBond(contact){
        // Check with 1kpe.
        var cAt1 = contact[1];
        var cAt2 = contact[2];

        if(cAt1.bonds.indexOf(cAt2)>-1){
            return false;
        }
        // Here we see if this contact is broadly parallel to an existing known bond, and if so reject it.
        // This probably belongs in the CloseContacts method.
        var mc = vec3.create([cAt1.x()-cAt2.x(),cAt1.y()-cAt2.y(),cAt1.z()-cAt2.z()]);
        vec3.normalize(mc);
        for(var ib=0;ib<cAt1.bonds.length;ib++){
            var mb = vec3.create([cAt1.x()-cAt1.bonds[ib].x(),cAt1.y()-cAt1.bonds[ib].y(),cAt1.z()-cAt1.bonds[ib].z()]);
            vec3.normalize(mb);
            if(vec3.dot(mb,mc)>0.7){
                return false;
            }
        }

        return true;
    }

    if(ligandAtoms.length>0){
        console.log("Getting metals");
        var metalsInBindingSite = model.getAtoms("allmetal and neighb cid=\"ligands\" group=atom maxd=4.0 exclude=water");
        console.log("Got metals");
        if(metalsInBindingSite.length>0){
            // FIXME - Should check metal distances.
            console.log("Getting close contacts");
            var metalContacts = model.CloseContacts(metalsInBindingSite,allAtoms,0.0,4.0,true);
            console.log("Got close contacts");
            metalContacts = metalContacts.filter(checkOverlapWithBond);
            console.log("Filtered close contacts");
            var redColours = colourScheme.colourOneColour([1.0,0.0,0.0,1.0]);
            var contactsPrimitiveInfo = contactsToCappedCylindersInfo(metalContacts,0.05,redColours,true);
            console.log("Got close contacts primitives");
            objects.push(contactsPrimitiveInfo);
            console.log("Added close contacts primitives");
        }
    }


    var glycoBlocks = getGlycoBlocks(model,0.2,colourScheme);
    console.log(glycoBlocks);
    Array.prototype.push.apply(objects,glycoBlocks);

    return objects;
 
}

var wizards = {"Bonds":wizardBonds,"Ribbons by chain":wizardRibbonsByChain,"Ribbons by secondary structure":wizardRibbons,"Worms":wizardWorms,"Site and ribbons by chain":wizardSiteAndRibbonsByChain};
