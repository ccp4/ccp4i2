define (['jquery'], function(){

function profileSVG(drawingOffset, drawingSize, profileNode) {
    var maxDataValue = parseInt(profileNode.getAttribute('max_data_value'))
    var nCols = parseInt(profileNode.getAttribute('width'))
    var nRows = parseInt(profileNode.getAttribute('height'))
    
    var cellWidth = parseInt(drawingSize.width)/nCols
    var cellHeight = parseInt(drawingSize.height)/nRows
    var cellWidthStr = cellWidth.toFixed(2)
    var cellHeightStr = cellHeight.toFixed(2)
    /*
     subProfile.set('type','averaged_only')
     subProfile.set('raw_data','')
     subProfile.set('mask','')
     else:
     subProfile.set('type','original_only')
     subProfile.set('alt_raw_data','')
     subProfile.set('alt_mask','')
     */
    var rawData = profileNode.getAttribute('raw_data').replace(/\s+/g, ' ').split(" ")
    var maskData = profileNode.getAttribute('mask').replace(/\s+/g, ' ').split(" ")
    var maskColor = 'cyan';
    
    if (profileNode.hasAttribute('type')){
        switch(profileNode.getAttribute('type')) {
            case 'averaged_only':
                rawData = profileNode.getAttribute('raw_data').replace(/\s+/g, ' ').split(" ")
                maskData = profileNode.getAttribute('mask').replace(/\s+/g, ' ').split(" ")
                maskColor = 'red';
                break;
            case 'dual':
                rawData = profileNode.getAttribute('raw_data').replace(/\s+/g, ' ').split(" ")
                maskData = profileNode.getAttribute('mask').replace(/\s+/g, ' ').split(" ")
                maskColor = 'green';
                break;
            default:
                maskColor = 'cyan';
        }
    }
    
    svgText = ''
    for (var ix=0; ix<nCols; ix++){
        for (var iy=0; iy<nRows; iy++){
            var left = drawingOffset.x + (ix*cellWidth)
            var leftStr = left.toFixed(2)
            var top = drawingOffset.y + (iy*cellHeight)
            var topStr = top.toFixed(2)
            var right = left+cellWidth;
            var rightStr = right.toFixed(2)
            var bottom = top+cellHeight;
            var bottomStr = bottom.toFixed(2)
            
            value = Math.round((255./maxDataValue) * (maxDataValue-rawData[ix*nRows+iy]))
            if (value >255) value = 255;
            var maskClass = 'in'
            if (maskData[ix*nRows+iy] == '-1') maskClass = 'out';
            svgText += '<rect x="'+leftStr+'" y="'+topStr+'" width="'+cellWidthStr+'" height="'+cellHeightStr+'" style="fill:rgb('+value+','+value+','+value+');fill-opacity:1.0;"/>'
            if (ix === 0) {
                svgText += '<line x1="'+leftStr+'" y1="'+topStr+'" x2="'+leftStr+'" y2="'+bottomStr+'" style="stroke:'+maskColor+';stroke-width:1;"/>';
            }
            else {
                var meIn = (maskData[ix*nRows+iy] !== '-1');
                var theeIn = ((maskData[(ix-1)*nRows+iy] !== '-1'))
                if ((meIn && !theeIn) || (theeIn && !meIn)) svgText += '<line x1="'+leftStr+'" y1="'+topStr+'" x2="'+leftStr+'" y2="'+bottomStr+'" style="stroke:'+maskColor+';stroke-width:1;"/>';
            }
            if (ix === nCols-1) svgText += '<line x1="'+rightStr+'" y1="'+topStr+'" x2="'+rightStr+'" y2="'+bottomStr+'" style="stroke:'+maskColor+';stroke-width:1;"/>';
            if (iy === 0) svgText += '<line x1="'+leftStr+'" y1="'+topStr+'" x2="'+rightStr+'" y2="'+topStr+'" style="stroke:'+maskColor+';stroke-width:1;"/>';
            else {
                var meIn = (maskData[(ix*nRows)+iy] !== '-1');
                var theeIn = ((maskData[(ix*nRows)+(iy-1)] !== '-1'))
                if ((meIn && !theeIn) || (theeIn && !meIn)) svgText += '<line x1="'+leftStr+'" y1="'+topStr+'" x2="'+rightStr+'" y2="'+topStr+'" style="stroke:'+maskColor+';stroke-width:1;"/>';
            }
            if (iy === nRows-1) svgText += '<line x1="'+leftStr+'" y1="'+bottomStr+'" x2="'+rightStr+'" y2="'+bottomStr+'" style="stroke:'+maskColor+';stroke-width:1;"/>';
        }
    }
    //Here for outlining masked region
    
    //alert(svgText)
    return svgText;
}

var mosflmModule = {
profileRenderer: function (txt, width, height, div){
    var xmlDoc = $.parseXML( txt ),
    $xml = $( xmlDoc );
    var profileNode = $xml.find( "profile" ).get(0)
    
    drawingSize = {'width':width, 'height':height}
    
    var svgText = '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" style="width:'+drawingSize.width+'; height:'+drawingSize.height+';border-width:0px;margin:0px;padding:0px;float:left;">\n'
    svgText = svgText + '<style type="text/css" >\n'
    svgText = svgText + '<![CDATA[\n'
    svgText = svgText + 'rect.in { stroke:cyan; stroke-width:1;}\n'
    svgText = svgText + 'rect.out { stroke-width:0;}\n'
    svgText = svgText + 'line {  stroke:cyan;stroke-width:4;}\n'
    svgText = svgText + 'text.outline {font-family:Arial, Helvetica, sans-serif; font-size:10px;stroke-width:4px; stroke:white; fill:white; }\n'
    svgText = svgText + 'text.foreground {font-family:Arial, Helvetica, sans-serif; font-size:10px;stroke-width:1; stroke:black; fill:black;}\n'
    svgText = svgText + ']]>\n'
    svgText = svgText+ '</style>\n'
    
    drawingOffset = {"x":0, "y":0}
    
    gridText = profileSVG(drawingOffset, drawingSize, profileNode)
    svgText += gridText
    svgText += '</svg>\n'
    div.innerHTML = svgText
},
    
profileGridRenderer: function (txt, width, height, div){
    var xmlDoc = $.parseXML( txt ),
    $xml = $( xmlDoc );
    var profileGridNode = $xml.find( "profile_grid" ).get(0)
    //alert(new XMLSerializer().serializeToString(profileGridNode))
    
    drawingSize = {'width':width, 'height':height}
    
    var svgText = '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" style="width:'+drawingSize.width+'; height:'+drawingSize.height+';border-width:0px;margin:0px;padding:0px;float:left;">\n'
    svgText = svgText + '<style type="text/css" >\n'
    svgText = svgText + '<![CDATA[\n'
    svgText = svgText + 'rect.in { stroke:cyan; stroke-width:1;}\n'
    svgText = svgText + 'rect.out { stroke-width:0;}\n'
    svgText = svgText + 'line {  stroke:black;stroke-width:1;}\n'
    svgText = svgText + 'text.outline {font-family:Arial, Helvetica, sans-serif; font-size:10px;stroke-width:4px; stroke:white; fill:white; }\n'
    svgText = svgText + 'text.foreground {font-family:Arial, Helvetica, sans-serif; font-size:10px;stroke-width:1; stroke:black; fill:black;}\n'
    svgText = svgText + ']]>\n'
    svgText = svgText+ '</style>\n'
    var num_x = parseInt(profileGridNode.getAttribute('num_x'))
    var num_y = parseInt(profileGridNode.getAttribute('num_y'))
    var max_x = parseInt(profileGridNode.getAttribute('max_width'))
    var max_y = parseInt(profileGridNode.getAttribute('max_height'))
    
    var drawingWidth=parseInt(drawingSize.width)
    var drawingHeight=parseInt(drawingSize.height)
    var blockStepX = drawingWidth / num_x
    var blockSizeX = blockStepX - 1
    var blockStepY = drawingWidth / num_y
    var blockSizeY = blockStepY - 1
    
    var blockSize = Math.min(blockSizeX, blockSizeY)
    var blockStep = Math.min(blockStepX, blockStepY)
    //alert(JSON.stringify(div.style))
    
    var ownerDocument = (profileGridNode.ownerDocument == null ? profileGridNode.documentElement : profileGridNode.ownerDocument.documentElement)
    var nsResolver = document.createNSResolver( ownerDocument );
    $xml.find('regional_profile').each(function(){
                                       var regional_profile = $(this).get(0);
                                       var iBox = parseInt(regional_profile.getAttribute('box')) - 1;
                                       var ixProfile = Math.floor(iBox / num_x);
                                       var iyProfile = iBox % num_x;
                                       svgText += profileSVG({'x':ixProfile * blockStep, 'y':iyProfile * blockStep},
                                                             {'width':blockSize,'height':blockSize}, regional_profile)
                                       })
    
    svgText += '</svg>\n'
    div.innerHTML = svgText
    //alert(svgText)
}
    
}
return mosflmModule

});
