#version 450

in vec2 pos, uv_vertex;
out vec2 uv_frag;

void main()
{
  gl_Position.xy = pos;
  gl_Position.zw = vec2(0.0f, 1.0f);

  uv_frag = uv_vertex;
}
