//Some code needed to handle javascrp texpectations of plotly.  The latter could be
//dnagerous !!
if (!Function.prototype.bind) {
    Function.prototype.bind = function(oThis) {
        if (typeof this !== 'function') {
            // closest thing possible to the ECMAScript 5
            // internal IsCallable function
            throw new TypeError('Function.prototype.bind - what is trying to be bound is not callable');
        }
        
        var aArgs   = Array.prototype.slice.call(arguments, 1),
        fToBind = this,
        fNOP    = function() {},
        fBound  = function() {
            return fToBind.apply(this instanceof fNOP
                                 ? this
                                 : oThis,
                                 aArgs.concat(Array.prototype.slice.call(arguments)));
        };
        
        if (this.prototype) {
            // Function.prototype doesn't have a prototype property
            fNOP.prototype = this.prototype;
        }
        fBound.prototype = new fNOP();
        
        return fBound;
    };
}
if (typeof Float64Array === 'undefined') Float64Array = Float32Array;

//Some lingering globals for the report

function fetchAndDecode(url) {
  return fetch(url.innerHTML).then(response => {
    if(!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    } else {
      if(response.headers.get("content-type") === "text/plain") {
        return response.text();
      } else {
        return url.innerHTML+" is not plain text. Ignoring."
      }
    }
  })
  .catch(e => {
    console.log(`There has been a problem with your fetch operation for resource "${url}": ` + e.message);
  })
  .finally(() => {
    //console.log(`fetch attempt for "${url}" finished.`);
  })
}

function toggleview(obj) {
    src = obj
    while ( (' ' + obj.className +' ').indexOf(' hidesection ') == -1){
        if      ( obj.nextSibling != null ) obj = obj.nextSibling;
        else if ( obj.parentNode  != null ) obj = obj.parentNode;
        else                                return;
    }
    if ( (' ' + obj.className +' ').indexOf(' displayed ') > -1 ) {
        $(obj).slideUp(50,function(){
            obj.className = obj.className.replace( 'displayed' , 'hidden' );
        })
        src.childNodes[0].nodeValue = src.childNodes[0].nodeValue.replace( "\u25BC", "\u25B6" );
    } else if ( (' ' + obj.className +' ').indexOf(' hidden ') > -1 )  {
        $(obj).slideDown(50,function(){
            var matches = obj.querySelectorAll(":scope .urlfetchtag");
            if(matches.length===1) {
                fetchAndDecode(matches[0]).then(res => {
                  matches[0].innerHTML = res
                  matches[0].classList.remove("urlfetchtag")
                })
            }
            obj.className = obj.className.replace( 'hidden' , 'displayed' );
        })
        src.childNodes[0].nodeValue = src.childNodes[0].nodeValue.replace( "\u25B6", "\u25BC" );
    }
}

function togglesubjob(obj) {
    
    src = obj
    while ( obj.className != "subjob" ) {
        if  ( obj.nextSibling != null ) obj = obj.nextSibling;
        else if ( obj.parentNode  != null ) obj = obj.parentNode;
        else                                return;
    }
    
    //alert('style.display ', obj.style.display);
    if ( obj.style.display == 'block' ) {
        obj.style.display = 'none';
        //src.childNodes[0].nodeValue = src.childNodes[0].nodeValue.replace( "Hide", "Show" );
        //src.childNodes[0].nodeValue = src.nodeValue.replace( "Hide", "Show" );
    } else {
        obj.style.display = 'block';
        //src.childNodes[0].nodeValue = src.childNodes[0].nodeValue.replace( "Show", "Hide" );
        //src.childNodes[0].nodeValue = src.nodeValue.replace( "Show", "Hide" );
    }
    
    if ( obj.childNodes.length == 0 ) {
        //alert("calling loadSubJobReport"+obj.tagName+obj.id);
        //Invoke pyqt tools to insert subjob report into page - see CCP4WebView
        SubJobReport.load(obj.id);
    }
    
    
}

function ccp4ReportOnLoad(){
    //Define slow scroll function
    (function($){$.fn.goTo = function() {
     $('html, body').animate({
                             scrollTop: $(this).offset().top + 'px'
                             }, 'slow');
     return this; // for chaining...
     }})(jQuery)
    
    drawXMLTexts()
    //Here find all divs that have a  data-widget-type attribute
    $('div[data-initially-drawn="True"]').each(function(iDiv, theDiv){
                                               if ($(theDiv).data("widget-type") == 'CCP4i2DrawnDiv') {
                                               drawDrawnDiv(theDiv);
                                               }
                                               else if ($(theDiv).data("widget-type") == 'CCP4i2LineChooser') {
                                               drawCCP4i2LineChooser(theDiv);
                                               }
                                               });

    galleryListObjClicked = function(obj){
        $(obj).parent().siblings().children('td').each(function(iTd, td){
                                                       $(td).removeClass('Selected');
                                                       $(td).addClass('NotSelected');
                                                       var correspondingDivId = $(td).attr('value') + '_div';
                                                       $("#"+correspondingDivId).removeClass('displayed')
                                                       $("#"+correspondingDivId).addClass('hidden')
                                                       });
        $(obj).removeClass('NotSelected');
        $(obj).addClass('Selected');
        var correspondingDivId = $(obj).attr('value') + '_div';
        $("#"+correspondingDivId).removeClass('hidden')
        $("#"+correspondingDivId).addClass('displayed')
        $("#"+correspondingDivId).find('div[data-initially-drawn]').each(function(iElement, element){drawDrawnDiv(element)});
    }
}

function reloadReportData(){
    drawXMLTexts();
    var xmlTexts =$('span[data-url]');
    var graphElements = $('div[data-renderer="CCP4i2Widgets.CCP4FlotRenderer"]');
    for (var i=0; i<graphElements.length;i++) {
        var graphElement = graphElements[i];
        drawDrawnDiv(graphElement);
    }
    var tableElements = $('div[data-renderer="CCP4i2Widgets.CCP4TableRenderer"]');
    for (var i=0; i<tableElements.length;i++) {
        var tableElement = tableElements[i];
        drawDrawnDiv(tableElement);
    }
    var progressElements = $('div[data-renderer="CCP4i2Widgets.CCP4i2HTMLChunk"]');
    for (var i=0; i<progressElements.length;i++) {
        var progressElement = progressElements[i];
        drawDrawnDiv(progressElement);
    }
    var total = xmlTexts.length + graphElements.length + tableElements.length + progressElements.length;
    window.reloadDataBridge.reloadReportDataResult(total);
    return total;
}

function drawXMLTexts(){
    $("span[data-url]").each(function(){
     var txt = '';
     var xmlhttp = new XMLHttpRequest();
     var self = this;
     xmlhttp.open("GET", $(self).data('url'), false);
     xmlhttp.overrideMimeType('text/plain;');
     xmlhttp.send();
     txt = xmlhttp.responseText;
     $(self).html(txt);
     });
}

var requireGrabbedModules = {};

function drawDrawnDiv(aDiv){
    var myDrawnDiv=aDiv
    /*Important to make the following vars local to this routine, since execution of the draw may be deferred in the event
     *that the drawing involves a requireload of javascript, so that global variables will end up having changed values
     */
    var rendererName = myDrawnDiv.getAttribute('data-renderer')
    var rendererData = myDrawnDiv.getAttribute('data-data')
    var renderWidth = myDrawnDiv.getAttribute('data-width')
    var renderHeight = myDrawnDiv.getAttribute('data-height')
    if (myDrawnDiv.hasAttribute('data-require')){
        var requireKey = myDrawnDiv.getAttribute('data-require')
        if (typeof requireGrabbedModules[requireKey] === 'undefined'){
            require([requireKey],function(grabbedModule){
                    eval('var '+requireKey+'=grabbedModule;')
                    eval(rendererName)(rendererData, renderWidth, renderHeight, myDrawnDiv);
                    requireGrabbedModules[requireKey] = grabbedModule;
                    });
        }
        else {
            eval('var '+requireKey+'=requireGrabbedModules[requireKey];')
            eval(rendererName)(rendererData, renderWidth, renderHeight, myDrawnDiv);
        }
    }
    else{
        var putativeRenderFunction = window[rendererName];
        if (typeof putativeRenderFunction === "function") {
            putativeRenderFunction.apply(null,rendererData, renderWidth, renderHeight, myDrawnDiv);
        }
    }
}

function drawCCP4i2LineChooser(theDiv){
    var plotData = theDiv.getAttribute('data-data');
    var heightStr = theDiv.getAttribute('data-height');
    var tableWidthStr = theDiv.getAttribute('data-table-width');
    var contentWidthStr = theDiv.getAttribute('data-content-width');
    var myDiv = theDiv;
    require(['CCP4i2Widgets'],function(CCP4i2Widgets){
            var newChooser = CCP4i2Widgets.CCP4i2LineChooser(myDiv, tableWidthStr, contentWidthStr, heightStr, {data:plotData});
            });
}

require.config({
               //baseUrl: 'http://127.0.0.1:43434/report_files/0.1.0',
               paths: {
               // the left side is the module ID,
               // the right side is the path to
               // the jQuery file, relative to baseUrl.
               // Also, the path should NOT include
               // the '.js' file extension. This example
               // is using jQuery 1.9.0 located at
               // js/lib/jquery-1.9.0.js, relative to
               // the HTML page.
               jquery: 'jquery.min',
               mosflm: 'mosflm',
               flot: 'flot/jquery.flot',
               CCP4i2Widgets: 'CCP4i2Widgets',
               jspimple: 'jspimple'
               }
               });

require(['jquery','CCP4i2Widgets','jspimple'],function(jquery,CCP4i2Widgets,jspimple){
        $( document ).ready(function(){ccp4ReportOnLoad();})
        })
