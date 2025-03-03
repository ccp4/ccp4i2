function hereDoc(f) {
  return f.toString().
      replace(/^[^\/]+\/\*!?/, '').
      replace(/\*\/[^\/]+$/, '');
}

var pointspheres_vertex_shader_source = hereDoc(function() {/*!
    // This shader is basic
    attribute vec3 aVertexPosition;
    attribute vec4 aVertexColour;
    attribute vec3 aVertexNormal;

    uniform vec3 offset;
    uniform float size;
    uniform mat3 scaleMatrix;

    uniform mat4 uMVMatrix;
    uniform mat4 uMVINVMatrix;
    uniform mat4 uPMatrix;

    varying lowp vec4 vColor;
    varying lowp vec3 vNormal;
    varying lowp vec3 vPosition;
    varying mediump mat4 mvInvMatrix;

    varying lowp vec4 eyePos;

    void main(void) {

      vec4 theVert = vec4(size*scaleMatrix*aVertexPosition+offset,1.0);

      gl_Position = uPMatrix * uMVMatrix * theVert;
      vColor = aVertexColour;

      eyePos = uMVMatrix * theVert;

      vNormal = aVertexNormal;
      mvInvMatrix = uMVINVMatrix;
      vPosition = aVertexPosition;

    }
*/});
