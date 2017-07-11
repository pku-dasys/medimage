#ifndef _PROJ_H_
#define _PROJ_H_

#include "config.h"

void A(double proj[nAlpha][projWidth], double orig[imageHight][imageWidth]);
void AStar(double back[imageHight][imageWidth],double proj[nAlpha][projWidth]);

#endif