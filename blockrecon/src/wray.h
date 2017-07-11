#ifndef _WRAY_H_
#define _WRAY_H_

void wray(int np,int nr,int *line, float *weight, int *numb, float *snorm);

void A(float g[NPROJ][NRAY],float f[IMGSIZE][IMGSIZE]);
void A_blk(float g[NPROJ][NRAY],int _bid,float f[BLKSIZE][BLKSIZE]);

#endif
