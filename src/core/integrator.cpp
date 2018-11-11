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
  // minimal working test
  static int tmp = 0;
  color_buffer[tmp] = RGBA(1.0f, 1.0f, 1.0f, 1.0f);
  if(tmp > 0) color_buffer[tmp-1] = RGBA(0.0f, 0.0f, 0.0f, 1.0f);
  tmp = (tmp == vRes*hRes-1) ? 0 : tmp+1;
}
