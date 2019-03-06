#version 330

uniform vec2 resolution;
uniform vec2 r;
uniform float angle;
uniform float radius;
uniform float scale;

in vec2 coord;
out vec2 v_color;

vec2 rot(vec2 v){
    return mat2(cos(angle), sin(angle), -sin(angle), cos(angle)) * v;
}

void main(){
    vec2 p = scale * (r + rot(radius * coord));
    v_color = coord;
    gl_Position = vec4(p, 0.0, 1.0);
}
