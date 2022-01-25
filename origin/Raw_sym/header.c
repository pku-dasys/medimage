#include "header.h"

double gDelta[nAlpha][projWidth];
double init_f[imageHight][imageWidth];
double f[imageHight][imageWidth], v[imageHight][imageWidth];

/*
variable used by A.c and AStar.c
*/
double alpha[nAlpha];
double xy[imageHight];
double rot[imageHight][imageWidth];
double lines[imageHight][imageWidth];

/* variable used by minimize_*.c */
double update[imageHight][imageWidth], d[imageHight][imageWidth];
double g[imageHight][imageWidth], laplace_v[imageHight][imageWidth];
double d_xf[imageHight][imageWidth], d_yf[imageHight][imageWidth];
double n[imageHight][imageWidth];
double tempZeros[imageHight+2][imageWidth+2], tempOnes[imageHight+2][imageWidth+2];
double padded[imageHight+2][imageWidth+2];

double p[imageHight][imageWidth];
double temp[imageHight][imageWidth];
double Af[nAlpha][projWidth], r[nAlpha][projWidth], q[nAlpha][projWidth], Ad[nAlpha][projWidth];
