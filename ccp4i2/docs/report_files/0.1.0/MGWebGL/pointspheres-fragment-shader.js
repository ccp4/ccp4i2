function hereDoc(f) {
  return f.toString().
      replace(/^[^\/]+\/\*!?/, '').
      replace(/\*\/[^\/]+$/, '');
}

var pointspheres_fragment_shader_source = hereDoc(function() {/*!
    precision mediump float;
    varying lowp vec4 vColor;
    varying lowp vec3 vNormal;
    varying lowp vec3 vPosition;
    varying lowp vec4 eyePos;
    varying mediump mat4 mvInvMatrix;

    uniform vec4 fogColour;

    uniform float fog_end;
    uniform float fog_start;

    uniform vec4 clipPlane0;
    uniform vec4 clipPlane1;
    uniform vec4 clipPlane2;
    uniform vec4 clipPlane3;
    uniform vec4 clipPlane4;
    uniform vec4 clipPlane5;
    uniform vec4 clipPlane6;
    uniform vec4 clipPlane7;
    uniform int nClipPlanes;

    uniform vec4 light_positions;
    uniform vec4 light_colours_ambient;
    uniform vec4 light_colours_specular;
    uniform vec4 light_colours_diffuse;

    void main(void) {
      if(dot(eyePos, clipPlane0)<0.0){
       discard;
      }
      if(dot(eyePos, clipPlane1)<0.0){
       discard;
      }
      vec3 L;
      vec3 E;
      vec3 R;
      vec4 Iamb =vec4(0.0,0.0,0.0,0.0);
      vec4 Idiff=vec4(0.0,0.0,0.0,0.0);
      vec4 Ispec=vec4(0.0,0.0,0.0,0.0);
      vec3 norm = normalize(vPosition);

      //vec3 fdx = dFdx(normalize(vPosition));
      //vec3 fdy = dFdy(normalize(vPosition));
      //vec3 theNormal = cross(fdx, fdy);

      E = ( vec4(normalize(vNormal),1.0)).xyz;
      //for (i = 0; i<nLights&&i<8; i++) {
       L = normalize((mvInvMatrix *light_positions).xyz);
       R = normalize(-reflect(L,norm));
       //calculate Ambient Term:
       Iamb += light_colours_ambient;
       //calculate Diffuse Term:
       Idiff += light_colours_diffuse * max(dot(E,L), 0.0);
       // calculate Specular Term:
       float y = max(max(light_colours_specular.r,light_colours_specular.g),light_colours_specular.b);
       Ispec += light_colours_specular * pow(max(dot(E,L),0.0),16.);
       Ispec.a *= y;
      //}
      
      float FogFragCoord = abs(eyePos.z/eyePos.w);
      float fogFactor = (fog_end - FogFragCoord)/(fog_end - fog_start);
      fogFactor = 1.0 - clamp(fogFactor,0.0,1.0);

      vec4 theColor = vec4(vColor);

      vec4 color = (1.5*theColor*Iamb + 1.2*theColor* Idiff);
      color.a = vColor.a;
      color += Ispec;
      gl_FragColor = mix(color, fogColour, fogFactor );

 
    }
*/});
