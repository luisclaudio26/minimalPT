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

  scene.cam.compute_parameters(Vec3(0.0f,0.1f,0.0f), Vec3(0.0f,1.0f,0.0f), Vec3(0.0f,0.0f,-1.0f), 35.0f, 4.0/3.0f, 35.0f);
  //scene.cam.compute_parameters(Vec3(0.0f,1.5f,-1.0f), Vec3(0.0f,0.0f,-1.0f), Vec3(0.0f,-1.0f,-1.0f), 35.0f, 4.0/3.0f, 35.0f);

  Shape ball_floor(Vec3(0.0f,-60.0f,0.0), 60.0f);
  ball_floor.diff_color = RGB(1.0f, 0.0f, 0.0f);
  scene.add_primitive( ball_floor );

  Shape ball_red(Vec3(-0.3f,0.1f,-1.0f), 0.1f);
  ball_red.diff_color = RGB(1.0f, 1.0f, 1.0f);
  ball_red.type = DELTA;
  scene.add_primitive( ball_red );

  // Green ball will emit light.
  // Its radius is 0.1m, so its surface is 0.04pi m².
  // In order to emit 1 W, we need to irradiate
  // 1/(0.04pi) W/m². Then, for an isotropic light field,
  // radiance should be equal in all directions and integrate
  // to 1/(0.04pi) W/m², thus radiance must be
  //
  //    (1/(0.04pi))/2pi = 1/(0.08pi²) ~ 1.2665 W.m⁻².sr
  //
  // TODO: Is it worth creating some routines that automatically
  // compute radiance for uniform, isotropic lightsources for a given power?
  Shape ball_green(Vec3(0.0f,0.4f,-1.0f), 0.1f);
  ball_green.diff_color = RGB(0.0f, 0.0f, 1.0f);
  ball_green.emission = RGB(5.0f, 5.0f, 5.0f); // in W/m²sr!!!
  scene.add_primitive( ball_green );

  Shape ball_blue(Vec3(0.3f,0.1f,-1.0f), 0.1f);
  ball_blue.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive( ball_blue );

  Shape ball_glass(Vec3(0.0f,0.1f,-0.5f), 0.1f);
  //Shape ball_glass(Vec3(0.0f,0.1f,0.0f), 1.0f);
  ball_glass.diff_color = RGB(1.0f, 1.0f, 1.0f);
  ball_glass.type = GLASS;
  ball_glass.eta = 0.45f;
  scene.add_primitive(ball_glass);

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
