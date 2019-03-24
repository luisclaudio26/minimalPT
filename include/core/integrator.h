#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "geometry.h"
#include "scene.h"
#include "spectrum.h"

class Integrator
{
private:
  RGB direct_illumination_solidangle(const Scene& scene,
                                      const Ray& primary_ray,
                                      const Isect& isect);

  RGB direct_illumination_surfacearea(const Scene& scene,
                                      const Ray& primary_ray,
                                      const Isect& isect);

  RGB camera_path(const Scene& scene,
                  const Ray& primary_ray,
                  const Isect& isect,
                  int path_length = 1);

  RGB pathtracer(const Scene& scene,
                  const Ray& primary_ray,
                  const Isect& isect);

  RGB bd_path(const Scene& scene,
              const Ray& primary_ray,
              const Vec3& lens_normal,
              const Isect& isect,
              int path_length = 2);

  RGB bdpt(const Scene& scene,
            const Ray& primary_ray,
            const Vec3& lens_normal,
            const Isect& isect);

  void light_path(const Scene& scene, int path_length = 1);

  RGB radiance_measurement(const Scene& scene,
                            const Ray& primary_ray,
                            const Isect& isect);

  RGB normal_shading(const Scene& scene,
                      const Ray& primary_ray,
                      const Isect& isect);

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
  void render_patch(const Scene& scene, int I, int J, int w, int h);
  void reconstruct_image();
  void reconstruct_image(double elapsed_time); //instrumented version
};

#endif
