function toggleview(obj) {
    src = obj
    while ( obj.className != "hidesection" ) {
      if      ( obj.nextSibling != null ) obj = obj.nextSibling;
      else if ( obj.parentNode  != null ) obj = obj.parentNode;
      else                                return;
    }
    if ( obj.style.display == 'block' ) {
        obj.style.display = 'none';
        src.childNodes[0].nodeValue = src.childNodes[0].nodeValue.replace( "Hide", "Show" );
    } else {
        obj.style.display = 'block';
        src.childNodes[0].nodeValue = src.childNodes[0].nodeValue.replace( "Show", "Hide" );
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
        src.childNodes[0].nodeValue = src.childNodes[0].nodeValue.replace( "Hide", "Show" );
        /rc.childNodes[0].nodeValue = src.nodeValue.replace( "Hide", "Show" );
    } else {
        obj.style.display = 'block';
        src.childNodes[0].nodeValue = src.childNodes[0].nodeValue.replace( "Show", "Hide" );
        src.childNodes[0].nodeValue = src.nodeValue.replace( "Show", "Hide" );
    }

    if ( obj.childNodes.length == 0 ) {
      //alert("calling loadSubJobReport"+obj.tagName+obj.id);
      //Invoke pyqt tools to insert subjob report into page - see CCP4WebView
      SubJobReport.load(obj.id);
    }
    

}

(function($) {
    $.fn.goTo = function() {
        $('html, body').animate({
            scrollTop: $(this).offset().top + 'px'
        }, 'fast');
        return this; // for chaining...
    }
})(jQuery);