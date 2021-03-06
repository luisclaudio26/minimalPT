#include "../include/frontend/gui.h"
#include <iostream>

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



  //scene.cam.compute_parameters(Vec3(0.0f,0.0f,0.1f), Vec3(0.0f,1.0f,0.0f), Vec3(0.0f,0.5f,-1.0f), 35.0f, 4.0/3.0f, 35.0f);
  //scene.cam.compute_parameters(Vec3(0.0f,1.5f,-1.0f), Vec3(0.0f,0.0f,-1.0f), Vec3(0.0f,-1.0f,-1.0f), 35.0f, 4.0/3.0f, 35.0f);

  // CORNELL BOX
  scene.cam.compute_parameters(Vec3(0.0f,0.5f,1.5f), Vec3(0.0f,1.0f,0.0f), Vec3(0.0f,0.5f,-1.0f), 35.0f, 4.0/3.0f, 35.0f);

  Shape ball_floor(Vec3(0.0f,-60.0f,0.0f), 60.0f);
  ball_floor.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive(ball_floor);

  Shape ball_ceil(Vec3(0.0f,60.0f+1.0f,0.0f), 60.0f);
  ball_ceil.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive(ball_ceil);

  Shape ball_right(Vec3(60.0f+0.5f, 0.0f, 0.0f), 60.0f);
  ball_right.diff_color = RGB(0.0f, 1.0f, 0.0f);
  scene.add_primitive(ball_right);

  Shape ball_left(Vec3(-(60.0f+0.5f), 0.0f, 0.0f), 60.0f);
  ball_left.diff_color = RGB(1.0f, 0.0f, 0.0f);
  scene.add_primitive(ball_left);

  Shape ball_back(Vec3(0.0f, 0.0f, -(60.0f+1.0f)), 60.0f);
  ball_back.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive(ball_back);

  Shape ball_1(Vec3(0.15f,0.1f,-0.2f), 0.1f);
  ball_1.diff_color = RGB(1.0f, 1.0f, 1.0f);
  scene.add_primitive(ball_1);

  Shape ball_2(Vec3(-0.15f,0.1f,-0.2f), 0.1f);
  ball_2.type = DELTA;
  scene.add_primitive(ball_2);

  Shape ball_light(Vec3(0.0f,1.1f,0.0f), 0.2f);
  ball_light.emission = RGB(4.0f, 4.0f, 2.0f);
  ball_light.diff_color = RGB(0.0f, 0.0f, 0.0f);
  scene.add_primitive(ball_light);


  // SETUP II
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

  //SETUP 1
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



  // -------------------------------------------------
  // configure integrator and film settings ----------
  Integrator integrator;

  integrator.start_rendering(scene);

  std::chrono::high_resolution_clock::time_point tS = std::chrono::high_resolution_clock::now();

  /*
  for(int i = 0; i < 5; ++i)
  {
    integrator.dump_image();
    std::cout<<"\rPreview "<<i;
    std::cout.flush();
  }
  */

  for(auto& t : integrator.render_jobs) t.join();
  std::chrono::high_resolution_clock::time_point tE = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(tE - tS);
  double acc_time = time_span.count();

  printf("\n%f s\n", acc_time);

  // invoke renderer ----------
  /*
  nanogui::init();

  GUI myGUI(scene, integrator);
  myGUI.drawAll();
  myGUI.setVisible(true);

  // main loop is also responsible for invoking
  // the integrator in order to request samples
  nanogui::mainloop();

  nanogui::shutdown();
  */

  return 0;
}
