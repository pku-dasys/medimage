#ifndef _CONFIG_H_
#define _CONFIG_H_

#define PHANTOM_FILE "phantom.dat"

const int IMGSIZE = 512;
const int NPROJ = 180;
const int NRAY = (IMGSIZE+IMGSIZE/2);

const int MIDRAY = NRAY/2;
const int MIDPIX = IMGSIZE/2;


const int DIMSIZE = 4;
const int NBLK = DIMSIZE*DIMSIZE;
const int BLKSIZE = IMGSIZE/DIMSIZE;

const int NDIR = 4;
const int DX[NDIR] = {-1,0,1,0};
const int DY[NDIR] = {0,1,0,-1};

#define UP 0
#define RIGHT 1
#define DOWN 2
#define LEFT 3

const float ALPHA = 20;
const float BETA = 1.35;
const float EPSILON = 0.001;

//const float G_RELAX = 1;

const int ALL_ITER = 50;

const float F_TOL = 1e-8;
const float GRID_TOL = 1e-2;

const float PI = 3.14159265359;

#endif
