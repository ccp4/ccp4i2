require.config({
               paths: {
               'jquery': 'jquery-1.11.0.min',
               'fogandclip': 'fogandclip'
               }
               });

define(['jquery'],function(apple,banana,pineapple,orange,plum,pear,strawberry,raspberry,lemon){

function FogClip(theDiv){
  this.myDiv = document.getElementById(theDiv);
  this.canvas = document.createElement("canvas");
  this.myDiv.appendChild(this.canvas);
  var computedStyle = getComputedStyle(this.myDiv, null);
  this.canvas.setAttribute("width",parseInt(computedStyle.width));
  this.canvas.setAttribute("height",parseInt(computedStyle.height));
  this.mouseDown = false;

  this.currentFrontX = parseInt(this.canvas.width/2-this.canvas.height/4);
  this.currentBackX = parseInt(this.canvas.width/2+this.canvas.height/4);

  this.mouseDownButton = -1;

  window.addEventListener("resize",this,false);
  window.addEventListener("mousemove",this,false);
  window.addEventListener("mousewheel",this,false);
  window.addEventListener("DOMMouseScroll",this,false);
  window.addEventListener("mousedown",this,false);
  window.addEventListener("mouseup",this,false);

  this.l1Active = false;
  this.l2Active = false;
  this.bothActive = false;
  this.l1Near = false;
  this.l2Near = false;
  this.bothNear = false;

  this.snapDistance = 15;
  this.barHeight = 14;
  this.barCornerRadius = 6;
  this.barButtonWidth = 10;

  if(this.myDiv.getAttribute("data-barButtonWidth")){
    this.barButtonWidth = parseInt(this.myDiv.getAttribute("data-barButtonWidth"));
  }
  if(this.myDiv.getAttribute("data-barHeight")){
    this.barHeight = parseInt(this.myDiv.getAttribute("data-barHeight"));
  }
  if(this.myDiv.getAttribute("data-barCornerRadius")){
    this.barCornerRadius = parseInt(this.myDiv.getAttribute("data-barCornerRadius"));
  }
  this.draw();

}

FogClip.prototype.handleEvent = function(evt) {
  if(evt.type=="resize"){
    var computedStyle = getComputedStyle(this.myDiv, null);
    this.canvas.setAttribute("width",computedStyle.width);
    this.canvas.setAttribute("height",computedStyle.height);
    console.log(this.canvas.width);
    console.log(this.canvas.height);
    this.draw();
  }
  if(evt.type=="mouseup"){
    this.mouseDown = false;
    this.l1Active = false;
    this.l2Active = false;
    this.bothActive = false;
    var rect = this.canvas.getBoundingClientRect();
    var w = rect["width"];
    var x = evt.clientX-rect["left"];
    var y = evt.clientY-rect["top"];
    this.oldX = x;
    var ctx = this.canvas.getContext("2d");
    if(!this.l1Active&&!this.l2Active&&!this.bothActive){
      var inSB = false;
      ctx.save();
      ctx.beginPath();
      ctx.rect(0,this.canvas.height-this.barHeight,this.canvas.width,this.barHeight-1);
      if(ctx.isPointInPath(x,y)){
        inSB = true;
      }
      ctx.restore();
      ctx.save();
      ctx.beginPath();
      ctx.rect(this.currentFrontX-1,this.canvas.height-this.barHeight,this.currentBackX-this.currentFrontX+2,this.barHeight-1);
      if(ctx.isPointInPath(x,y)){
        inSB = false;
      }
      ctx.restore();
      if(inSB){
        var midpoint = (this.currentBackX+this.currentFrontX)/2;
        var dx = 0;
        if(midpoint>x){
          dx = -this.canvas.width/10;
        } else {
          dx = this.canvas.width/10;
        }
        if(this.currentFrontX+dx>=(this.barButtonWidth+this.barCornerRadius)&&this.currentBackX+dx<w-(this.barButtonWidth+this.barCornerRadius)){
          this.currentBackX += dx;
          this.currentFrontX += dx;
          this.draw();
        } else if(this.currentBackX+dx>=w-(this.barButtonWidth+this.barCornerRadius)){
          dx = w - this.currentBackX-(this.barButtonWidth+this.barCornerRadius);
          this.currentBackX += dx;
          this.currentFrontX += dx;
          this.draw();
        } else if(this.currentFrontX+dx<(this.barButtonWidth+this.barCornerRadius)){
          dx = -this.currentFrontX+(this.barButtonWidth+this.barCornerRadius);
          this.currentBackX += dx;
          this.currentFrontX += dx;
          this.draw();
        }
      }
    }

  }
  if(evt.type=="mousedown"){
    this.mouseDown = true;
    this.mouseDownButton = evt.button;
    var rect = this.canvas.getBoundingClientRect();
    var x = evt.clientX-rect["left"];
    var y = evt.clientY-rect["top"];
    this.oldX = x;
    this.oldY = y;
    this.downX = evt.clientX;
    this.downY = evt.clientY;
  }

  if(evt.type=="DOMMouseScroll"){
    this.downX = evt.clientX;
    this.downY = evt.clientY;
    var rect = this.canvas.getBoundingClientRect();
    var x = evt.clientX-rect["left"];
    var y = evt.clientY-rect["top"];
    var w = rect["width"];
    downX = this.downX-rect["left"];
    downY = this.downY-rect["top"];
    if((downX>=0)&&(downY>=0)&&(downX<this.canvas.width)&&(downY<this.canvas.height)){
      var dy;
      if(evt.detail>0){
        dy = 2;
      } else {
        dy = -2;
      }
      if(this.currentFrontX+dy>=(this.barButtonWidth+this.barCornerRadius)&&this.currentBackX-dy<w-(this.barButtonWidth+this.barCornerRadius)){
        if((this.currentFrontX+dy)<=(this.currentBackX-dy)){
          this.currentBackX -= dy;
          this.currentFrontX += dy;
          this.draw();
        }
      }
    }
  }

  if(evt.type=="mousewheel"){
    this.downX = evt.clientX;
    this.downY = evt.clientY;
    var rect = this.canvas.getBoundingClientRect();
    var x = evt.clientX-rect["left"];
    var y = evt.clientY-rect["top"];
    var w = rect["width"];
    downX = this.downX-rect["left"];
    downY = this.downY-rect["top"];
    if((downX>=0)&&(downY>=0)&&(downX<this.canvas.width)&&(downY<this.canvas.height)){
      var dy;
      if(evt.wheelDelta<0){
        dy = 2;
      } else {
        dy = -2;
      }
      if(this.currentFrontX+dy>=(this.barButtonWidth+this.barCornerRadius)&&this.currentBackX-dy<w-(this.barButtonWidth+this.barCornerRadius)){
        if((this.currentFrontX+dy)<=(this.currentBackX-dy)){
          this.currentBackX -= dy;
          this.currentFrontX += dy;
          this.draw();
        }
      }
    }
  }

  if(evt.type=="mousemove"){
    this.l1Near = false;
    this.l2Near = false;
    this.bothNear = false;
    var rect = this.canvas.getBoundingClientRect();
    var x = evt.clientX-rect["left"];
    var y = evt.clientY-rect["top"];
    var w = rect["width"];
    downX = this.downX-rect["left"];
    downY = this.downY-rect["top"];

    if(((this.mouseDown&&evt.metaKey)||(this.mouseDown&&this.mouseDownButton==1))){
      if((downX>=0)&&(downY>=0)&&(downX<this.canvas.width)&&(downY<this.canvas.height)){
        var dy = (y - this.oldY)/8;
        if(this.currentFrontX+dy>=(this.barButtonWidth+this.barCornerRadius)&&this.currentBackX-dy<w-(this.barButtonWidth+this.barCornerRadius)){
          if((this.currentFrontX+dy)<=(this.currentBackX-dy)){
            this.currentBackX -= dy;
            this.currentFrontX += dy;
            this.draw();
          }
        }
      }
      return;
    }

    if(!this.mouseDown){
      var l1Diff = Number.MAX_VALUE;
      var l2Diff = Number.MAX_VALUE;
      var ctx = this.canvas.getContext("2d");
      ctx.save();
      ctx.beginPath();
      ctx.rect(this.currentFrontX-this.snapDistance,0,this.snapDistance*2,this.canvas.height-this.barHeight);
      if(ctx.isPointInPath(x,y)){
        this.l1Near = true;
        l1Diff = Math.abs(x-this.currentFrontX);
      }
      ctx.restore();
      ctx.save();
      ctx.beginPath();
      ctx.rect(this.currentBackX-this.snapDistance,0,this.snapDistance*2,this.canvas.height-this.barHeight);
      if(ctx.isPointInPath(x,y)){
        this.l2Near = true;
        l2Diff = Math.abs(x-this.currentBackX);
      }
      ctx.restore();
      if(this.l1Near&&this.l2Near){
        if(l1Diff<l2Diff){
          this.l2Near = false;
        } else {
          this.l1Near = false;
        }
      }
      if(!this.l1Near&&!this.l2Near){
        ctx.save();
        ctx.beginPath();
        ctx.rect(this.currentFrontX-1,this.canvas.height-this.barHeight,this.currentBackX-this.currentFrontX+2,this.barHeight-1);
        if(ctx.isPointInPath(x,y)){
          this.bothNear = true;
        }
        ctx.restore();
      }
    }

    if(this.mouseDown&&!this.l1Active&&!this.l2Active&&!this.bothActive){

      var ctx = this.canvas.getContext("2d");
      ctx.save();
      ctx.beginPath();
      ctx.rect(this.currentFrontX-this.snapDistance,0,this.snapDistance*2,this.canvas.height-this.barHeight);
      if(ctx.isPointInPath(x,y)){
        this.l1Active = true;
      }
      ctx.restore();
      ctx.save();
      ctx.beginPath();
      ctx.rect(this.currentBackX-this.snapDistance,0,this.snapDistance*2,this.canvas.height-this.barHeight);
      if(ctx.isPointInPath(x,y)){
        this.l2Active = true;
      }
      ctx.restore();
      if(!this.l1Active&&!this.l2Active){
        ctx.save();
        ctx.beginPath();
        ctx.rect(this.currentFrontX-1,this.canvas.height-this.barHeight,this.currentBackX-this.currentFrontX+2,this.barHeight-1);
        if(ctx.isPointInPath(x,y)){
          this.bothActive = true;
        }
        ctx.restore();
      }
    }

    if(x>=(this.barButtonWidth+this.barCornerRadius)&&x<w-(this.barButtonWidth+this.barCornerRadius)){
      downX = this.downX-rect["left"];
      downY = this.downY-rect["top"];
      if(downX>=0&&downX<this.canvas.width&&downY>=0&&downY<this.canvas.height){
      if(this.l1Active){
        this.currentFrontX = x;
        if(this.currentFrontX>this.currentBackX){
          this.currentFrontX = this.currentBackX;
        }
      } else if(this.l2Active){
        this.currentBackX = x;
        if(this.currentFrontX>this.currentBackX){
          this.currentBackX = this.currentFrontX;
        }
      } else if(this.bothActive){
        var dx = x - this.oldX;
        if(this.currentFrontX+dx>=(this.barButtonWidth+this.barCornerRadius)&&this.currentBackX+dx<w-(this.barButtonWidth+this.barCornerRadius)){
          this.currentBackX += dx;
          this.currentFrontX += dx;
        }
      }
      }
    }

    this.oldX = x;
    this.oldY = y;
    this.draw();

  }

}

Fog.prototype = Object.create(FogClip.prototype);
function Fog(div){
  FogClip.call(this,div);
}

Clip.prototype = Object.create(FogClip.prototype);
function Clip(div){
  FogClip.call(this,div);
}

Clip.prototype.draw = function() {
  var ctx = this.canvas.getContext("2d");
  this.drawMost(ctx);
  this.drawClip(ctx);
  this.drawRest(ctx);
  ctx.font = "9pt Arial";
  ctx.fillText("Front: "+this.currentFrontX, 0, 10)
  ctx.fillText("Back: "+this.currentBackX, 0, 20)
}

Fog.prototype.draw = function() {
  var ctx = this.canvas.getContext("2d");
  this.drawMost(ctx);
  this.drawFog(ctx);
  this.drawRest(ctx);
  ctx.font = "9pt Arial";
  ctx.fillText("Fog start: "+this.currentFrontX, 0, 10)
  ctx.fillText("Fog end: "+this.currentBackX, 0, 20)
}

FogClip.prototype.drawMost = function(ctx) {

  ctx.save();
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.arc(this.canvas.width/2,this.canvas.height/2,this.canvas.height/4,0,Math.PI*2);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  var x1 = 0;
  var y1 = this.canvas.height-this.barHeight;
  var x2 = 0;
  var y2 = this.canvas.height;
  var linearGradient1 = ctx.createLinearGradient(x1, y1, x2, y2);
  linearGradient1.addColorStop(0, 'rgb(134, 134, 134)');
  linearGradient1.addColorStop(0.1, 'rgb(170, 170, 170)');
  linearGradient1.addColorStop(0.25, 'rgb(180, 180, 180)');
  linearGradient1.addColorStop(0.5, 'rgb(198, 198, 198)');
  linearGradient1.addColorStop(0.75, 'rgb(180, 180, 180)');
  linearGradient1.addColorStop(0.9, 'rgb(170, 170, 170)');
  linearGradient1.addColorStop(1, 'rgb(134, 134, 134)');

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.fillStyle = linearGradient1;
  ctx.moveTo(this.barCornerRadius,this.canvas.height-this.barHeight);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barHeight);
  ctx.arc(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barHeight+this.barCornerRadius,this.barCornerRadius,Math.PI*3/2,0);
  ctx.lineTo(this.canvas.width-1,this.canvas.height-this.barCornerRadius-1);
  ctx.arc(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barCornerRadius-1,this.barCornerRadius,0,Math.PI/2);
  ctx.lineTo(this.barCornerRadius,this.canvas.height-1);
  ctx.arc(this.barCornerRadius,this.canvas.height-this.barCornerRadius-1,this.barCornerRadius,Math.PI/2,Math.PI);
  ctx.lineTo(0,this.canvas.height-this.barHeight+this.barCornerRadius);
  ctx.arc(this.barCornerRadius,this.canvas.height-this.barHeight+this.barCornerRadius,this.barCornerRadius,Math.PI,Math.PI*3/2);
  ctx.fill();
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.lineWidth = 1;
  ctx.strokeStyle = "#000000";
  ctx.moveTo(this.barCornerRadius,this.canvas.height-this.barHeight);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barHeight);
  ctx.arc(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barHeight+this.barCornerRadius,this.barCornerRadius,Math.PI*3/2,0);
  ctx.lineTo(this.canvas.width-1,this.canvas.height-this.barCornerRadius-1);
  ctx.arc(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barCornerRadius-1,this.barCornerRadius,0,Math.PI/2);
  ctx.lineTo(this.barCornerRadius,this.canvas.height-1);
  ctx.arc(this.barCornerRadius,this.canvas.height-this.barCornerRadius-1,this.barCornerRadius,Math.PI/2,Math.PI);
  ctx.lineTo(0,this.canvas.height-this.barHeight+this.barCornerRadius);
  ctx.arc(this.barCornerRadius,this.canvas.height-this.barHeight+this.barCornerRadius,this.barCornerRadius,Math.PI,Math.PI*3/2);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  var linearGradientButton = ctx.createLinearGradient(x1, y1, x2, y2);
  linearGradientButton.addColorStop(0, 'rgb(240, 240, 240)');
  linearGradientButton.addColorStop(1, 'rgb(180, 180, 180)');

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.fillStyle = linearGradientButton;
  ctx.moveTo(this.barCornerRadius,this.canvas.height-this.barHeight);
  ctx.lineTo(this.barCornerRadius+this.barButtonWidth,this.canvas.height-this.barHeight);
  ctx.lineTo(this.barCornerRadius+this.barButtonWidth,this.canvas.height-1);
  ctx.lineTo(this.barCornerRadius,this.canvas.height-1);
  ctx.arc(this.barCornerRadius,this.canvas.height-this.barCornerRadius-1,this.barCornerRadius,Math.PI/2,Math.PI);
  ctx.lineTo(0,this.canvas.height-this.barHeight+this.barCornerRadius);
  ctx.arc(this.barCornerRadius,this.canvas.height-this.barHeight+this.barCornerRadius,this.barCornerRadius,Math.PI,Math.PI*3/2);
  ctx.fill();
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.lineWidth = 1;
  ctx.strokeStyle = "#000000";
  ctx.moveTo(this.barCornerRadius,this.canvas.height-this.barHeight);
  ctx.lineTo(this.barCornerRadius+this.barButtonWidth,this.canvas.height-this.barHeight);
  ctx.lineTo(this.barCornerRadius+this.barButtonWidth,this.canvas.height-1);
  ctx.lineTo(this.barCornerRadius,this.canvas.height-1);
  ctx.arc(this.barCornerRadius,this.canvas.height-this.barCornerRadius-1,this.barCornerRadius,Math.PI/2,Math.PI);
  ctx.lineTo(0,this.canvas.height-this.barHeight+this.barCornerRadius);
  ctx.arc(this.barCornerRadius,this.canvas.height-this.barHeight+this.barCornerRadius,this.barCornerRadius,Math.PI,Math.PI*3/2);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.fillStyle = linearGradientButton;
  ctx.moveTo(this.canvas.width-this.barCornerRadius-1-this.barButtonWidth,this.canvas.height-this.barHeight);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barHeight);
  ctx.arc(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barHeight+this.barCornerRadius,this.barCornerRadius,Math.PI*3/2,0);
  ctx.lineTo(this.canvas.width-1,this.canvas.height-this.barCornerRadius-1);
  ctx.arc(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barCornerRadius-1,this.barCornerRadius,0,Math.PI/2);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-1-this.barButtonWidth,this.canvas.height-1);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-1-this.barButtonWidth,this.canvas.height-this.barHeight);
  ctx.fill();
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.lineWidth = 1;
  ctx.strokeStyle = "#000000";
  ctx.moveTo(this.canvas.width-this.barCornerRadius-1-this.barButtonWidth,this.canvas.height-this.barHeight);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barHeight);
  ctx.arc(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barHeight+this.barCornerRadius,this.barCornerRadius,Math.PI*3/2,0);
  ctx.lineTo(this.canvas.width-1,this.canvas.height-this.barCornerRadius-1);
  ctx.arc(this.canvas.width-this.barCornerRadius-1,this.canvas.height-this.barCornerRadius-1,this.barCornerRadius,0,Math.PI/2);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-1-this.barButtonWidth,this.canvas.height-1);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-1-this.barButtonWidth,this.canvas.height-this.barHeight);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.lineWidth = 1;
  ctx.moveTo(this.canvas.width-this.barCornerRadius-1-this.barButtonWidth+4,this.canvas.height-this.barHeight+3);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-2,this.canvas.height-this.barHeight/2-0.5);
  ctx.lineTo(this.canvas.width-this.barCornerRadius-1-this.barButtonWidth+4,this.canvas.height-4);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.lineWidth = 1;
  ctx.moveTo(this.barCornerRadius-1+this.barButtonWidth-3,this.canvas.height-this.barHeight+3);
  ctx.lineTo(this.barCornerRadius,this.canvas.height-this.barHeight/2-0.5);
  ctx.lineTo(this.barCornerRadius-1+this.barButtonWidth-3,this.canvas.height-4);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  var linearGradient2 = ctx.createLinearGradient(x1, y1, x2, y2);
  linearGradient2.addColorStop(0, 'rgb(134, 134, 255)');
  linearGradient2.addColorStop(1, 'rgb(50, 50, 255)');

  var linearGradient3 = ctx.createLinearGradient(x1, y1, x2, y2);
  linearGradient3.addColorStop(0, 'rgb(154, 154, 255)');
  linearGradient3.addColorStop(1, 'rgb(70, 70, 255)');

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  if(this.bothNear||this.bothActive){
    ctx.fillStyle = linearGradient3;
  } else {
    ctx.fillStyle = linearGradient2;
  }
  ctx.fillRect(this.currentFrontX,this.canvas.height-this.barHeight,this.currentBackX-this.currentFrontX,this.barHeight-1);
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.lineWidth = 1;
  ctx.strokeStyle = "#000000";
  ctx.strokeRect(this.currentFrontX,this.canvas.height-this.barHeight,this.currentBackX-this.currentFrontX,this.barHeight-1);
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.lineWidth = 1;
  ctx.moveTo((this.currentFrontX+this.currentBackX)/2,this.canvas.height-this.barHeight+3);
  ctx.lineTo((this.currentFrontX+this.currentBackX)/2,this.canvas.height-3);
  if(this.currentBackX-this.currentFrontX>8){
    ctx.moveTo((this.currentFrontX+this.currentBackX)/2-4,this.canvas.height-this.barHeight+3);
    ctx.lineTo((this.currentFrontX+this.currentBackX)/2-4,this.canvas.height-3);
    ctx.moveTo((this.currentFrontX+this.currentBackX)/2+4,this.canvas.height-this.barHeight+3);
    ctx.lineTo((this.currentFrontX+this.currentBackX)/2+4,this.canvas.height-3);
  }
  ctx.stroke();
  ctx.closePath();
  ctx.restore();
}

FogClip.prototype.drawClip = function(ctx) {
  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0,0,this.currentFrontX,this.canvas.height-this.barHeight);
  ctx.fillRect(this.currentBackX,0,this.canvas.width,this.canvas.height-this.barHeight);
  ctx.closePath();
  ctx.restore();
}

FogClip.prototype.drawFog = function(ctx) {
  x1 = 0;
  y1 = 0;
  x2 = this.canvas.width;
  y2 = 0;
  var linearGradientAlpha = ctx.createLinearGradient(x1, y1, x2, y2);
  var startFrac = parseFloat(this.currentFrontX)/parseFloat(this.canvas.width);
  var endFrac = parseFloat(this.currentBackX)/parseFloat(this.canvas.width);
  linearGradientAlpha.addColorStop(startFrac, 'rgba(255, 255, 255, 0.0)');
  linearGradientAlpha.addColorStop(endFrac, 'rgba(255, 255, 255, 1.0)');

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.fillStyle = linearGradientAlpha;
  ctx.fillRect(0,0,this.canvas.width,this.canvas.height-this.barHeight);
  ctx.closePath();
  ctx.restore();
}

FogClip.prototype.drawRest = function(ctx) {
  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  if(this.l1Near||this.l1Active){
    ctx.lineWidth = 2;
    ctx.strokeStyle = "#77ff77";
  } else {
    ctx.strokeStyle = "#00ff00";
    ctx.lineWidth = 1;
  }
  ctx.moveTo(this.currentFrontX,0);
  ctx.lineTo(this.currentFrontX,this.canvas.height-this.barHeight);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  if(this.l2Near||this.l2Active){
    ctx.lineWidth = 2;
    ctx.strokeStyle = "#ff7777";
  } else {
    ctx.lineWidth = 1;
    ctx.strokeStyle = "#ff0000";
  }
  ctx.moveTo(this.currentBackX,0);
  ctx.lineTo(this.currentBackX,this.canvas.height-this.barHeight);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.lineWidth = 1;
  ctx.moveTo(this.canvas.width/2-3,this.canvas.height/2);
  ctx.lineTo(this.canvas.width/2+3,this.canvas.height/2);
  ctx.moveTo(this.canvas.width/2,this.canvas.height/2-3);
  ctx.lineTo(this.canvas.width/2,this.canvas.height/2+3);
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  this.drawButton(ctx,this.canvas.width-4-this.barHeight,4,this.barHeight,this.barHeight,'+');
  this.drawButton(ctx,this.canvas.width-4-this.barHeight,4+this.barHeight+2,this.barHeight,this.barHeight,'-');
  this.drawButton(ctx,this.canvas.width-4-this.barHeight,4+this.barHeight+2+this.barHeight+2,this.barHeight,this.barHeight,'c');
}

FogClip.prototype.drawButton = function(ctx,x1,y1,w,h,style,active) {
  var linearGradientButton = ctx.createLinearGradient(0, y1, 0, y1+h);
  if(active){
    linearGradientButton.addColorStop(0, 'rgb(250, 250, 250)');
    linearGradientButton.addColorStop(1, 'rgb(220, 220, 220)');
  }else{
    linearGradientButton.addColorStop(0, 'rgb(240, 240, 240)');
    linearGradientButton.addColorStop(1, 'rgb(180, 180, 180)');
  }

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.fillStyle = linearGradientButton;
  ctx.fillRect(x1,y1,w,h);
  ctx.closePath();
  ctx.restore();

  ctx.save();
  ctx.translate(0.5, 0.5);
  ctx.beginPath();
  ctx.strokeStyle = '#000000';
  ctx.strokeRect(x1,y1,w,h);
  if(style=="+"){
    ctx.moveTo(x1+w/2,y1+2);
    ctx.lineTo(x1+w/2,y1+h-2);
  }
  if(style=="+"||style=="-"){
    ctx.moveTo(x1+2,y1+h/2);
    ctx.lineTo(x1+w-2,y1+h/2);
  }
  ctx.stroke();
  ctx.closePath();
  ctx.restore();

  ctx.save();
  //ctx.translate(0.5, 0.5);
  ctx.beginPath();
  if(style=="c"){
    ctx.fillRect(x1+w/2-1,y1+w/2-1,3,3);
  }
  ctx.closePath();
  ctx.restore();
}

return {
    Fog : function(theDiv){
        var fog = new Fog(theDiv);
        return fog;
    },
    Clip : function(theDiv){
        var clip = new Clip(theDiv);
        return clip;
    }
}

});
