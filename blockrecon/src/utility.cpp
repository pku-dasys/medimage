#include <fstream>

#include <cassert>

#include "config.h"

using namespace std;

void window(float &x,float lo,float hi) {
    if (x<lo) x = lo;
    if (x>hi) x = hi;
}

float sqr(float x) {
    return x*x;
}

int xy2bid(int x,int y,int &bx,int &by) {
    bx = x%BLKSIZE;
    by = y%BLKSIZE;
    return (x/BLKSIZE)*DIMSIZE+(y/BLKSIZE);
}

void ind2xy(int ind,int &x,int &y) {
    x = ind/IMGSIZE;
    y = ind%IMGSIZE;
}

void block(float img[IMGSIZE][IMGSIZE],float blk[NBLK][BLKSIZE][BLKSIZE]) {
    for (int i = 0; i<IMGSIZE; ++i)
        for (int j = 0; j<IMGSIZE; ++j) {
            int bx,by;
            int bid = xy2bid(i,j,bx,by);
            blk[bid][bx][by] = img[i][j];
        }
}

void unblock(float blk[NBLK][BLKSIZE][BLKSIZE],float img[IMGSIZE][IMGSIZE]) {
    for (int i = 0; i<IMGSIZE; ++i)
        for (int j = 0; j<IMGSIZE; ++j) {
            int bx,by;
            int bid = xy2bid(i,j,bx,by);
            img[i][j] = blk[bid][bx][by];
        }
}

void normalize(float img[IMGSIZE][IMGSIZE]) {
    float minval = 1e30,maxval = -1e30;
    for (int i = 0; i<IMGSIZE; ++i)
        for (int j = 0; j<IMGSIZE; ++j){
            float tmp = img[i][j];
            if (minval>tmp) minval = tmp;
            if (maxval<tmp) maxval = tmp;
        }
    maxval -= minval;
    for (int i = 0; i<IMGSIZE; ++i)
        for (int j = 0; j<IMGSIZE; ++j)
            img[i][j] = (img[i][j]-minval)/maxval*255.;
}

void read_file(float img[IMGSIZE][IMGSIZE], const char *file_name) {
    ifstream fin(file_name);
    assert(fin);
    for (int i = 0; i<IMGSIZE; ++i)
        for (int j = 0; j<IMGSIZE; ++j)
            fin>>img[i][j];
    fin.close();
}

void write_file(float img[IMGSIZE][IMGSIZE], const char *file_name) {
    ofstream fou(file_name);
    assert(fou);
    for (int i = 0; i<IMGSIZE; ++i) {
        for (int j = 0; j<IMGSIZE; ++j)
            fou<<img[i][j]<<' ';
        fou<<endl;
    }
    fou.close();
}

void write_data(float g[NPROJ][NRAY], const char *file_name) {
    ofstream fou(file_name);
    assert(fou);
    for (int i = 0; i<NPROJ; ++i) {
        for (int j = 0; j<NRAY; ++j)
            fou<<g[i][j]<<' ';
        fou<<endl;
    }
    fou.close();
}
