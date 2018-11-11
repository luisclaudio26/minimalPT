#ifndef DIRECT_H
#define DIRECT_H

#include "../core/integrator.h"

class DirectIllumination : public Integrator
{
private:
  RGB integrate(const Vec2& uv, const Scene& scene) const override;
};

#endif
