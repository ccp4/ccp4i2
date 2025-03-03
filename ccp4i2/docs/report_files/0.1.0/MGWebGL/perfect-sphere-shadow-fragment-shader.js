function hereDoc(f) {
  return f.toString().
      replace(/^[^\/]+\/\*!?/, '').
      replace(/\*\/[^\/]+$/, '');
}

var perfect_sphere_shadow_fragment_shader_source = hereDoc(function() {/*!
#extension GL_OES_element_index : enable
#extension GL_EXT_frag_depth : enable
    precision mediump float;
    varying lowp vec4 vColor;
    varying lowp vec3 vNormal;
    varying lowp vec4 eyePos;
    varying lowp vec2 vTexture;
    varying mediump mat4 mvMatrix;

    varying lowp vec4 ShadowCoord;

    uniform float fog_end;
    uniform float fog_start;

    uniform sampler2D uSampler;

    uniform vec4 clipPlane0;
    uniform vec4 clipPlane1;
    uniform vec4 clipPlane2;
    uniform vec4 clipPlane3;
    uniform vec4 clipPlane4;
    uniform vec4 clipPlane5;
    uniform vec4 clipPlane6;
    uniform vec4 clipPlane7;
    uniform int nClipPlanes;
    uniform vec4 fogColour;

    uniform vec4 light_positions;
    uniform vec4 light_colours_ambient;
    uniform vec4 light_colours_specular;
    uniform vec4 light_colours_diffuse;

    varying mediump mat4 projMatrix;
    varying float size_v;

    uniform sampler2D ShadowMap;
    //FIXME  - my buffer is currently always 1024 x 1024. This may change.
    //uniform float xPixelOffset;
    //uniform float yPixelOffset;

    float lookup(vec2 offSet){
      float xPixelOffset = 1.0/1024.0;
      float yPixelOffset = 1.0/1024.0;
      vec4 coord = ShadowCoord + vec4(offSet.x * xPixelOffset * ShadowCoord.w, offSet.y * yPixelOffset * ShadowCoord.w, 0.07, 0.0);
      if(coord.s>1.0||coord.s<0.0||coord.t>1.0||coord.t<0.0)
          return 1.0;
      //gl_FragColor = texture2D(ShadowMap, coord.xy );
      float shad2 = texture2D(ShadowMap, coord.xy ).x;
      shad2 = shad2/(coord.p/coord.q);
      shad2 = clamp(shad2,0.0,1.0);
      if(shad2<0.9){
          shad2 = 0.0;
      }
      return shad2;
    }

    void main(void) {
      float x = 2.0*(vTexture.x-.5);
      float y = 2.0*(vTexture.y-.5);
      float zz =  1.0 - x*x - y*y;
      float z = sqrt(zz);

      if(zz <= 0.06)
          discard;
    
      if(dot(eyePos, clipPlane0)<0.0){
       discard;
      }
      if(dot(eyePos, clipPlane1)<0.0){
       discard;
      }
      gl_FragColor = vec4(zz*vColor.r, zz*vColor.g, zz*vColor.b, 1.0);


      vec4 pos = eyePos;
      pos.z += 0.7071*z*size_v;
      pos = projMatrix * pos;

      gl_FragDepthEXT = (pos.z / pos.w + 1.0) / 2.0;

      vec3 L;
      vec3 E;
      vec3 R;
      vec4 Iamb =vec4(0.0,0.0,0.0,0.0);
      vec4 Idiff=vec4(0.0,0.0,0.0,0.0);
      vec4 Ispec=vec4(0.0,0.0,0.0,0.0);
       vec3 norm = vec3(x,y,z);

      // Baffled as to why I have to do this.
      vec3 invertLightPos = vec3(light_positions.x,-light_positions.y,light_positions.z);
      E = ( vec4(normalize(norm),1.0)).xyz;
      //for (i = 0; i<nLights&&i<8; i++) {
       L = normalize(invertLightPos);
       //R = normalize(-reflect(L,norm));
       //calculate Ambient Term:
       Iamb += light_colours_ambient;
       //calculate Diffuse Term:
       Idiff += light_colours_diffuse * max(dot(E,L), 0.0);
       // calculate Specular Term:
       y = max(max(light_colours_specular.r,light_colours_specular.g),light_colours_specular.b);
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
      //gl_FragColor = color;
      
      // FIXME - this should be done before fogging.

      int shadowQuality = 1;
      float shad2 = 0.0;
      if(shadowQuality==0){
        shad2 = lookup(vec2(0,0));
      } else if(shadowQuality==1){
          for(float y = -1.5 ; y <=1.5 ; y+=1.0){
              for (float x = -1.5 ; x <=1.5 ; x+=1.0){
                 shad2 += lookup(vec2(x*0.5,y*0.5));
              }
          }
          shad2 /= 16.0 ;
          shad2 += 0.2;
      } else if(shadowQuality==2){
          for(float y = -3.5 ; y <=3.5 ; y+=1.0){
              for (float x = -3.5 ; x <=3.5 ; x+=1.0){
                 shad2 += lookup(vec2(x*0.5,y*0.5));
              }
          }
          shad2 /= 64.0 ;
          shad2 += 0.2;
      }

      shad2 = clamp(shad2,0.0,1.0);
      

      gl_FragColor = vec4(color.r*shad2, color.g*shad2,color.b*shad2, 1.0);
    }
*/});
