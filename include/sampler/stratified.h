#ifndef STRATIFIED_H
#define STRATIFIED_H

#include "../core/sampler.h"

class Stratified : public Sampler
{
private:
  
public:
  float get1D() override;
};

#endif
