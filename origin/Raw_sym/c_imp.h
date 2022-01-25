#ifndef _C_IMP_H_
#define _C_IMP_H_

#include "header.h"

void medfilt2(double f[imageHight][imageWidth]);
double sqr(double x);
void randn(double x[nAlpha][imageHight]);
int iround(double x);
double x2Grid(double x);
void Clip(double *px1, double *py1, double *px2, double *py2);
void RotateImage(double ISrc[imageHight][imageWidth], double IDst[imageHight][imageWidth],
                   double a,double xy[imageHight]);

#endif
