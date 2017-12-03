#pragma once

#include "ct3d.h"

void forward_proj(const Parameter &args,
                  float sx,float sy,float sz,
                  float dx,float dy,float dz,
                  int64_t *ind,float *wgt,int &numb);

void parallel(const Parameter &args,
              int alpha, int detectorX, int detectorY,
              float &srcX, float &srcY, float &srcZ,
              float &dstX, float &dstY, float &dstZ);

void cone(const Parameter &args,
          int alpha, int detectorX, int detectorY,
          float &srcX, float &srcY, float &srcZ,
          float &dstX, float &dstY, float &dstZ);