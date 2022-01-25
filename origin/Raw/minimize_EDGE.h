#ifndef _MINIMIZE_EDGE_H_
#define _MINIMIZE_EDGE_H_

#include "header.h"

void minimize_in_EDGE_variable(double out[imageHight][imageWidth], double f[imageHight][imageWidth],
                                 double v[imageHight][imageWidth],double alpha,double beta,double epsilon,
                                 double ED_TOL,int IterEdge);

#endif
