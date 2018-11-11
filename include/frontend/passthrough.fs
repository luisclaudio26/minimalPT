#version 450

// Texture sampler: this is our actual rendered
// image and should be the same size as the screen
//uniform sampler2D img;
uniform sampler2D color_buffer;

in vec2 uv_frag;
out vec4 pixel;

void main()
{
  pixel = texture2D(color_buffer, uv_frag);
}
