#include "../include/frontend/gui.h"

int main(int argc, char** args)
{
  srand(0);

  // "load" scene -----------
  // units are given in meters!!!
  // QUESTION: with 35mm film and 35mm lens an object
  // of radius r should be at a distance of 2r to be
  // completely in scene, but this is not the case for
  // these spheres, which must be at roughly 2.25r for
  // this. why is this? Must test with planes after!!!
  Scene scene;
  scene.prims.push_back( Shape(Vec3(0.0f,0.0f,-0.9f), 0.4f) );

  // configure integrator and film settings ----------
  Integrator integrator;

  // invoke renderer ----------
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
