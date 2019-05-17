#ifndef _HEADER_H_
#define _HEADER_H_

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

extern double gDelta[nAlpha][imageHight];
extern double init_f[imageHight][imageWidth];
extern double f[imageHight][imageWidth], v[imageHight][imageWidth];

/*
variable used by A.c and AStar.c
*/
extern double alpha[nAlpha];
extern double xy[imageHight];
extern double rot[imageHight][imageWidth];
extern double lines[imageHight][imageWidth];

/* variable used by minimize_*.c */
extern double update[imageHight][imageWidth], d[imageHight][imageWidth];
extern double g[imageHight][imageWidth], laplace_v[imageHight][imageWidth];
extern double d_xf[imageHight][imageWidth], d_yf[imageHight][imageWidth];
extern double n[imageHight][imageWidth];
extern double tempZeros[imageHight+2][imageWidth+2], tempOnes[imageHight+2][imageWidth+2];
extern double padded[imageHight+2][imageWidth+2];

extern double p[imageHight][imageWidth];
extern double temp[imageHight][imageWidth];
extern double Af[nAlpha][imageHight], r[nAlpha][imageHight], q[nAlpha][imageHight], Ad[nAlpha][imageHight];

#endif
