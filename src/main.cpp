#include "../include/frontend/gui.h"

int main(int argc, char** args)
{
  srand(0);

  // "load" scene -----------
  Scene scene;
  scene.prims.push_back( Shape(Vec3(0.0f,-2.0f,0.0f), 1.0f) );

  // configure integrator and film settings ----------
  Integrator integrator;

  // --------- invoke renderer ----------
  nanogui::init();

  GUI myGUI(scene, integrator);
  myGUI.drawAll();
  myGUI.setVisible(true);

  // main loop is also responsible for invoking
  // the integrator in order to request samples
  nanogui::mainloop();

  nanogui::shutdown();

  return 0;
}
