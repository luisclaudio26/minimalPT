#ifndef GUI_H
#define GUI_H

#include <nanogui/opengl.h>
#include <nanogui/glutil.h>
#include <nanogui/screen.h>
#include <nanogui/window.h>
#include <vector>

#include "../core/integrator.h"
#include "../core/spectrum.h"

class GUI : public nanogui::Screen
{
private:
  nanogui::GLShader shader;
  GLuint color_buffer_gpu;

  Integrator& integrator;
  const Scene& scene;
public:
  GUI(const Scene& scene, Integrator& integrator);
  virtual void drawContents();
  virtual void draw(NVGcontext *ctx);
};

#endif
