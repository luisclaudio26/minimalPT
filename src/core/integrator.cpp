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
  color_buffer.resize(vRes*hRes);
  for(RGBA& p : color_buffer)
    p = RGBA(0.1f, 0.1f, 0.1f, 1.0f);
}

void Integrator::render(const Scene& scene)
{
  // loop over pixels in film, sampling them
  // on the center. trace primary ray and return
  // black/white if we had an intersection or not
  for(int j = 0; j < vRes; ++j)
    for(int i = 0; i < hRes; ++i)
    {
      RGBA sample_color(0.0f, 0.0f, 0.0f, 1.0f);

      // coordinates of the sample in [0,1]² range
      // (film coordinate system), flipping the v-axis (j-axis)
      Vec2 uv( (i+0.5f)/hRes, ((vRes-1)-j+0.5f)/vRes );

      // get primary ray and check intersections
      Ray primary_ray = scene.cam.get_primary_ray(uv);
      Isect isect;
      if( scene.cast_ray(primary_ray, isect) )
      {
        RGBA normal = RGBA(isect.normal.x, isect.normal.y, isect.normal.z, 1.0f);
        sample_color = (normal+1.0f)*0.5f;
      }

      // output to color buffer
      color_buffer[j*hRes+i] = sample_color;
    }
}
