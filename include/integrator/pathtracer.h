#ifndef PATHTRACER_H
#define PATHTRACER_H

#include "../core/integrator.h"

class Pathtracer : public Integrator
{
private:
  RGB integrate(const Vec2& uv, const Scene& scene) const override;
};

#endif
