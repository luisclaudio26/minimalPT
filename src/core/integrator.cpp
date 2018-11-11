#include "../../include/core/integrator.h"
#include <cstdio>
#include <chrono>
using namespace std::chrono;

/*
  high_resolution_clock::time_point tS = high_resolution_clock::now();
  high_resolution_clock::time_point tE = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(tE - tS);
*/

Integrator::Integrator()
  : vRes(600), hRes(800)
{
  // clear color buffer
  int n = vRes*hRes;

  /*
  color_buffer.resize(vRes*hRes);
  for(RGBA& p : color_buffer)
    p = RGBA(0.1f, 0.1f, 0.1f, 1.0f);
    */

  samples.insert(samples.begin(), n, RGB(0.0f,0.0f,0.0f));
  weights.insert(weights.begin(), n, 0.0f);
    frame.insert(  frame.begin(), n, RGBA(0.1f,0.1f,0.1f,1.0f));
}

RGBA Integrator::normal_shading(const Scene& scene,
                                const Ray& primary_ray,
                                const Isect& isect)
{
  return (RGBA(isect.normal,1.0f)+1.0f)*0.5f;
}

void Integrator::render(const Scene& scene)
{
  // loop over pixels in film, sampling them
  // on the center. trace primary ray and return
  // black/white if we had an intersection or not
  for(int j = 0; j < vRes; ++j)
    for(int i = 0; i < hRes; ++i)
    {
      RGB sample(0.0f, 0.0f, 0.0f);

      // uniformly sample pixel (i,j) and map sample coordinates
      // to [0,1]² with j-axis flipping (film coordinate system)
      float e1 = (float)rand()/RAND_MAX, e2 = (float)rand()/RAND_MAX;
      Vec2 uv( (i+e1)/hRes, ((vRes-1)-j+e2)/vRes );

      // get primary ray and first intersection,
      // then invoke shader to compute the sample value
      Ray primary_ray = scene.cam.get_primary_ray(uv);

      Isect isect;
      if( scene.cast_ray(primary_ray, isect) )
        sample = normal_shading(scene, primary_ray, isect);

      // sample splatting
      int sample_add = j*hRes+i;
      samples[sample_add] += sample;
      weights[sample_add] += 1.0f;
    }

  // reconstruct image. this is expensive as hell!
  for(int i = 0; i < vRes*hRes; ++i)
    frame[i] = RGBA(glm::min(RGB(1.0f), samples[i]/weights[i]), 1.0f);
}
