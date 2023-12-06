function BezierCurve(carts, accu){
  var spline = [];
  var nsteps = accu * (carts.length/3);
  var tstep = 1.0 / parseFloat(nsteps);

  var t = 0.0;
  for (var k = 0; k < nsteps; k++){
   var cartsi = JSON.parse(JSON.stringify(carts));

   for (var j = carts.length-3; j > 0; j-=3){
    for (var i = 0; i < j; i+=3){
     cartsi[i]   = (1-t)*cartsi[i]   + t*cartsi[i+3];
     cartsi[i+1] = (1-t)*cartsi[i+1] + t*cartsi[i+4];
     cartsi[i+2] = (1-t)*cartsi[i+2] + t*cartsi[i+5];
    }
   }

   spline.push(cartsi[0]);
   spline.push(cartsi[1]);
   spline.push(cartsi[2]);
   t += tstep;
  }

  return spline;

}

/****** lincrv.c ******/
/* Ken Shoemake, 1994 */
/* DialASpline(t,a,p,m,n,Cn,interp,val) computes a point val at parameter
    t on a spline with knot values a and control points p. The curve will have
    Cn continuity, and if interp is TRUE it will interpolate the control points.
    Possibilities include Langrange interpolants, Bezier curves, Catmull-Rom
    interpolating splines, and B-spline curves. Points have m coordinates, and
    n+1 of them are provided. The work array must have room for n+1 points.
 */

function DialASpline(t, a,  p, Cn, interp, output, idx, work)
{
    var i, j, k, h, lo, hi;

    var n = p.length/3 - 1;

    if (Cn>n-1) Cn = n-1;       /* Anything greater gives one polynomial */
    for (k=0; t> a[k]&&k<a.length; k++);    /* Find enclosing knot interval */
    for (h=k; t==a[k]; k++);    /* May want to use fewer legs */
    if (k>n) {k = n; if (h>k) h = k;}
    h = 1+Cn - (k-h); k--;
    lo = k-Cn; hi = k+1+Cn;


    if (interp) {               /* Lagrange interpolation steps */
        var drop=0;
        if (lo<0) {lo = 0; drop += Cn-k;
                   if (hi-lo<Cn) {drop += Cn-hi; hi = Cn;}}
        if (hi>n) {hi = n; drop += k+1+Cn-n;
                   if (hi-lo<Cn) {drop += lo-(n-Cn); lo = n-Cn;}}
        for (i=lo; i<=hi; i++){
          work[3*i]   = p[3*i];
          work[3*i+1] = p[3*i+1];
          work[3*i+2] = p[3*i+2];
        }
        for (j=1; j<=Cn; j++) {
            for (i=lo; i<=hi-j; i++) {
                var t0=(a[i+j]-t)/(a[i+j]-a[i]), t1=1-t0;
                work[3*i]   = t0*work[3*i]   + t1*work[3*(i+1)];
                work[3*i+1] = t0*work[3*i+1] + t1*work[3*(i+1)+1];
                work[3*i+2] = t0*work[3*i+2] + t1*work[3*(i+1)+2];
            }
        }
        h = 1+Cn-drop;
    } else {                    /* Prepare for B-spline steps */
        if (lo<0) {h += lo; lo = 0;}
        for (i=lo; i<=lo+h; i++){
          work[3*i]   = p[3*i];
          work[3*i+1] = p[3*i+1];
          work[3*i+2] = p[3*i+2];
        }
        if (h<0) h = 0;
    }
    for (j=0; j<h; j++) {
        var tmp = 1+Cn-j;
        for (i=h-1; i>=j; i--) {
            var t0=(a[lo+i+tmp]-t)/(a[lo+i+tmp]-a[lo+i]), t1=1-t0;
            work[3*(lo+i+1)]   = t0*work[3*(lo+i)]   + t1*work[3*(lo+i+1)];
            work[3*(lo+i+1)+1] = t0*work[3*(lo+i)+1] + t1*work[3*(lo+i+1)+1];
            work[3*(lo+i+1)+2] = t0*work[3*(lo+i)+2] + t1*work[3*(lo+i+1)+2];
        }
    }

    //console.log(t+" "+lo+" "+hi+" "+p.length);
    //console.log("done spline point");
    //console.log(work[3*(lo+h)]+" "+work[3*(lo+h)+1]+" "+work[3*(lo+h)+2]);
    output[3*idx]   = work[3*(lo+h)];
    output[3*idx+1] = work[3*(lo+h)+1];
    output[3*idx+2] = work[3*(lo+h)+2];
}

function SplineCurve(ctlPts, accu, Cn, interp){
   var nsteps = ctlPts.length/3 * accu;
   var  BIG = 1.0e8;
   var knots = [];
   var maxx = -BIG;
   var minx = BIG;
   var maxy = -BIG;
   var miny = BIG;
   var maxz = -BIG;
   var minz = BIG;
   var maxt = -BIG;
   var mint = BIG;
   var tstep;
   var knotstep;

   var outputi; // Cartesian

   for(var i=0;i<ctlPts.length;i+=3){
     minx = Math.min(ctlPts[i],minx);
     maxx = Math.max(ctlPts[i],maxx);
     miny = Math.min(ctlPts[i+1],minx);
     maxy = Math.max(ctlPts[i+1],maxx);
     minz = Math.min(ctlPts[i+2],minx);
     maxz = Math.max(ctlPts[i+2],maxx);
   }

   mint = minx;
   mint = Math.min(mint,miny);
   mint = Math.min(mint,minz);

   maxt = maxx;
   maxt = Math.max(maxt,maxy);
   maxt = Math.max(maxt,maxz);

   tstep = (maxt-mint)/parseFloat(nsteps-1);
   knotstep = (maxt-mint)/(ctlPts.length/3-1);

   for(var i=0;i<ctlPts.length/3;i++){
     knots.push(mint + parseFloat(i)*knotstep);
   }
   knots.push(tstep);

   var output = [];
   var work = [];
   //console.log("mint: " + mint);
   //console.log("maxt: " + maxt);
   
   for (var ii=0;ii<nsteps;ii++){
     var t = mint + ii*tstep;
     DialASpline(t, knots, ctlPts, Cn, interp,output,ii,work);
   }
   //console.log(output.length);
   //console.log(knots.length+" "+nsteps);

   return output;
   
}
