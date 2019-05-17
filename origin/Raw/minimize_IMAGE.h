#ifndef _MINIMIZE_IMAGE_H_
#define _MINIMIZE_IMAGE_H_

#include "header.h"

void minimize_in_IMAGE_variable(double out[imageHight][imageWidth], double g[nAlpha][imageHight],
                                  double f[imageHight][imageWidth],double v[imageHight][imageWidth],
                                  double alpha,double epsilon,double IM_TOL,int IterImage);

#endif
