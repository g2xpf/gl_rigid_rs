#version 330
in vec2 coord;

uniform vec2 r;
uniform vec2 wh;
uniform float scale;

void main(){
    vec2 p = r + coord * wh;
    gl_Position = vec4(scale * p, 0.0, 1.0);
}
