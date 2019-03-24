#include "../include/frontend/gui.h"
#include <nanogui/button.h>
#include <nanogui/layout.h>
#include <iostream>
#include <chrono>
#include <cstdlib>

GUI::GUI(const Scene& scene, Integrator& integrator, int n_threads)
  : scene(scene),
    integrator(integrator),
    tp(n_threads, integrator, scene),
    nanogui::Screen(Eigen::Vector2i(100, 100), "Andaluz renderer", false)
{
  using namespace nanogui;

  // TODO: initFromFiles is looking for path relative to the
  // first .cpp that is including gui.h. Report to github?!
  shader.initFromFiles("passthrough",
                        "../include/frontend/passthrough.vs",
                        "../include/frontend/passthrough.fs");

  // upload the triangles we'll use to draw onto
  Eigen::MatrixXf quad(2, 6);
  quad.col(0)<<-1.0f, -1.0f; //lower triangle
  quad.col(1)<<+1.0f, -1.0f;
  quad.col(2)<<+1.0f, +1.0f;
  quad.col(3)<<-1.0f, -1.0f; //upper triangle
  quad.col(4)<<+1.0f, +1.0f;
  quad.col(5)<<-1.0f, +1.0f;

  Eigen::MatrixXf uv(2, 6);
  uv.col(0)<<0.0f, 1.0f; //lower triangle
  uv.col(1)<<1.0f, 1.0f;
  uv.col(2)<<1.0f, 0.0f;
  uv.col(3)<<0.0f, 1.0f; //upper triangle
  uv.col(4)<<1.0f, 0.0f;
  uv.col(5)<<0.0f, 0.0f;

  shader.bind();
  shader.uploadAttrib<Eigen::MatrixXf>("pos", quad);
  shader.uploadAttrib<Eigen::MatrixXf>("uv_vertex", uv);

  // --------- set GUI parameters ----------
  //resize according to scene
  this->setSize( Eigen::Vector2i(integrator.hRes, integrator.vRes) );

  // initialize texture which will hold image
  glGenTextures(1, &color_buffer_gpu);
  glBindTexture(GL_TEXTURE_2D, color_buffer_gpu);
  glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, this->width(), this->height());

  // --------- set threadpool ------------
  initialize_job_list();
  tp.spawn();
};

void GUI::initialize_job_list()
{
  const int patch_sz = 64;

  for(int j = 0; j < this->height(); j += patch_sz)
    for(int i = 0; i < this->width(); i += patch_sz)
    {
      Threadpool::Job job;
      job.spp = 0;
      job.i = i; job.j = j;

      int remaining_w = this->width() - i;
      job.w = remaining_w >= patch_sz ? patch_sz : remaining_w;
      int remaining_h = this->height() - j;
      job.h = remaining_h >= patch_sz ? patch_sz : remaining_h;

      tp.jobs.push(job);
    }
}

void GUI::drawContents()
{
  using namespace std::chrono;
  static high_resolution_clock::time_point t_s = high_resolution_clock::now();

  // keep updating buffer
  shader.bind();

  // update buffers after N frames
  // TODO: why is the framerate so low??
  const int refresh_period = 60;
  static int frame_counter = 0;

  frame_counter++;
  if( frame_counter >= refresh_period )
  {
    tp.hold();
    high_resolution_clock::time_point t_e = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t_e - t_s);
    integrator.reconstruct_image(time_span.count());

    // copy CPU color_buffer to GPU
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_buffer_gpu);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, this->width(), this->height(),
                     GL_RGBA, GL_FLOAT, (void*)integrator.frame.data());

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    // resume
    frame_counter = 0;
    t_s = high_resolution_clock::now();
    tp.resume();
  }

  // draw stuff
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, color_buffer_gpu);
  shader.setUniform("color_buffer", 0);
  shader.drawArray(GL_TRIANGLES, 0, 6);
}

void GUI::draw(NVGcontext *ctx) { Screen::draw(ctx); }
