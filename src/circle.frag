#version 330

in vec2 v_color;
out vec4 color;

void main(){
    color = vec4(sin(v_color), cos(v_color.x), 0.1);
}
