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

  scene.cam.compute_parameters(Vec3(0.0f,0.5f,1.5f), Vec3(0.0f,1.0f,0.0f), Vec3(0.0f,0.5f,-1.0f), 35.0f, 4.0/3.0f, 35.0f);
  //scene.cam.compute_parameters(Vec3(0.0f,1.5f,-1.0f), Vec3(0.0f,0.0f,-1.0f), Vec3(0.0f,-1.0f,-1.0f), 35.0f, 4.0/3.0f, 35.0f);

  // <editor-fold> Scene I
  float R = 1000.0f;
  Shape ball_floor(Vec3(0.0f,-R,0.0f), R);
  ball_floor.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive(ball_floor);

  Shape ball_ceil(Vec3(0.0f,R+1.0f,0.0f), R);
  ball_ceil.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive(ball_ceil);

  Shape ball_right(Vec3(R+0.5f, 0.0f, 0.0f), R);
  ball_right.diff_color = RGB(0.0f, 1.0f, 0.0f);
  scene.add_primitive(ball_right);

  Shape ball_left(Vec3(-(R+0.5f), 0.0f, 0.0f), R);
  ball_left.diff_color = RGB(1.0f, 0.0f, 0.0f);
  scene.add_primitive(ball_left);

  Shape ball_back(Vec3(0.0f, 0.0f, -(R+1.0f)), R);
  ball_back.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive(ball_back);

  Shape ball_1(Vec3(0.0f,0.1f,-0.2f), 0.1f);
  ball_1.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive(ball_1);

  Shape ball_2(Vec3(-0.3f,0.1f,-0.2f), 0.1f);
  ball_2.type = DELTA;
  scene.add_primitive(ball_2);

  Shape ball_3(Vec3(0.3f,0.1f,-0.2f), 0.1f);
  ball_3.type = GLASS;
  ball_3.eta = 0.45f;
  scene.add_primitive(ball_3);

  float w = 10.0f;
  Shape ball_light(Vec3(0.0f,0.9f,-0.4f), 0.05f);
  //Shape ball_light(Vec3(0.0f,0.7f,-0.4f), 0.1f);
  ball_light.emission = RGB(w, w, w);
  ball_light.diff_color = RGB(0.0f, 0.0f, 0.0f);
  scene.add_primitive(ball_light);

  Shape ball_occluder(Vec3(0.0f,0.75f,-0.4f), 0.15f);
  ball_occluder.diff_color = RGB(0.0f, 0.0f, 1.0f);
  scene.add_primitive(ball_occluder);

  // </editor-fold>

  // <editor-fold> Scene II
  /*
  Shape ball_light(Vec3(0.0f, 1.0f, 0.0f), 0.2f);
  ball_light.emission = RGB(5.0f, 5.0f, 5.0f);
  scene.add_primitive( ball_light );

  Shape ball1(Vec3(0.2f, 0.1f, -0.3f), 0.1f);
  ball1.diff_color = RGB(1.0f, 0.0f, 0.0f);
  scene.add_primitive( ball1 );

  Shape ball2(Vec3(0.1f, 0.1f, -0.5f), 0.1f);
  ball2.diff_color = RGB(0.0f, 1.0f, 0.0f);
  scene.add_primitive( ball2 );

  Shape ball3(Vec3(0.0f, 0.1f, -0.7f), 0.1f);
  ball3.diff_color = RGB(0.0f, 0.0f, 1.0f);
  scene.add_primitive( ball3 );

  Shape ball4(Vec3(-0.2f, 0.1f, -1.0f), 0.1f);
  ball4.diff_color = RGB(0.0f, 1.0f, 1.0f);
  scene.add_primitive( ball4 );

  Shape ball5(Vec3(-0.4f, 0.1f, -1.3f), 0.1f);
  ball5.diff_color = RGB(1.0f, 0.0f, 1.0f);
  scene.add_primitive( ball5 );

  Shape ball_floor(Vec3(0.0f,-60.0f,0.0f), 60.0f);
  ball_floor.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive( ball_floor );
  */
  // </editor-fold>

  // <editor-fold> Scene III
  /*
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
  ball_glass.diff_color = RGB(1.0f, 1.0f, 1.0f);
  ball_glass.type = GLASS;
  ball_glass.eta = 0.45f;
  scene.add_primitive(ball_glass);
  */
  // </editor-fold>

  // -------------------------------------------------
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
