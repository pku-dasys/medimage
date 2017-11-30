#pragma once

#include "main.h"

void forward_proj(const Parameter &args,
                  float sx,float sy,float sz,
                  float dx,float dy,float dz,
                  int64 *ind,float *wgt,int &numb);