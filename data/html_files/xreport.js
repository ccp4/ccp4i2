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

