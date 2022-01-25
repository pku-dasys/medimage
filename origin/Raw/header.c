#include "header.h"

#define PHANTOM_FILE "phantom.txt"

#define imageHight 512
#define imageWidth 512

#define nAlpha 180
#define BETA 0.0002
#define ALPHA 1000
#define EPSILON 0.0001

#define _IterAlt 10
#define _IterImage 10
#define _IterEdge 10
#define _IM_TOL 0.0001
#define _ED_TOL 0.0001

double gDelta[nAlpha][imageHight];
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
double Af[nAlpha][imageHight], r[nAlpha][imageHight], q[nAlpha][imageHight], Ad[nAlpha][imageHight];
