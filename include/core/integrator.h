#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "geometry.h"
#include "scene.h"
#include "spectrum.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

class Integrator
{
private:
  RGB direct_illumination_solidangle(const Scene& scene,
                                      const Ray& primary_ray,
                                      const Isect& isect) const;

  RGB direct_illumination_surfacearea(const Scene& scene,
                                      const Ray& primary_ray,
                                      const Isect& isect) const;

  RGB camera_path(const Scene& scene,
                  const Ray& primary_ray,
                  const Isect& isect,
                  int path_length = 1) const;

  RGB pathtracer(const Scene& scene,
                  const Ray& primary_ray,
                  const Isect& isect) const;

  RGB radiance_measurement(const Scene& scene,
                            const Ray& primary_ray,
                            const Isect& isect) const;

  RGB normal_shading(const Scene& scene,
                      const Ray& primary_ray,
                      const Isect& isect) const;

  // maps a given irradiance measurement to a pixel output.
  // Ideally we should convert ENERGY measurements (irradiance
  // integrated over the pixel surface, integrated over time),
  // but initial tests with power measurements shown that it
  // is not trivial to deal with the extremely low power
  // measurements (as the pixel area is few um²). How realistic
  // it is to work with irradiances instead of actual power
  // measurements I'm not sure.
  RGBA camera_response_curve(const RGB& irradiance) const;

public:
  Integrator();

  // film data
  int vRes, hRes;
  float pixel_area;

  // the samples buffer which accumulate sample contributions,
  // and its respective per pixel weights, which we use for filtering.
  // frame holds the elementwise division of sample_buffer by weight.
  RadianceBuffer samples;
  std::vector<float> weights;
  ColorBuffer frame;

  void render(const Scene& scene);

  std::vector<std::thread> render_jobs;
  bool halt;
  std::mutex mtx;
  std::condition_variable cv;
  std::atomic<int> counter, done;

  void start_rendering(const Scene& scene);
  void dump_image();

  void render_patch(int start_x, int block_sz_x,
                    int start_y, int block_sz_y,
                    const Scene& scene);
};

#endif
