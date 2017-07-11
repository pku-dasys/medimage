#ifndef _UTILITY_H_
#define _UTILITY_H_

void window(float &x,float lo,float hi);

float sqr(float x);

void ind2xy(int ind,int &x,int &y);
int xy2bid(int x,int y,int &bx,int &by);

void read_file(float[IMGSIZE][IMGSIZE], const char *);
void write_file(float[IMGSIZE][IMGSIZE], const char *);
void write_data(float g[NPROJ][NRAY], const char *file_name);
void normalize(float img[IMGSIZE][IMGSIZE]);

void block(float img[IMGSIZE][IMGSIZE],float blk[NBLK][BLKSIZE][BLKSIZE]);
void unblock(float blk[NBLK][BLKSIZE][BLKSIZE],float img[IMGSIZE][IMGSIZE]);

#endif
