#include "testcase.h"

#include "utility.h"
#include "tracing.h"
#include "ct3d.h"

#include <cmath>
#include <fstream>

#include <boost/filesystem.hpp>

using namespace std;

void wrapper_shepp_logan(const Parameter &args, float *img, ushort *prj, int a,int x,int y, int *bf,int *br) {
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    if(args.BEAM == "Parallel")
        parallel(args,a,x,y,
                srcX,srcY,srcZ,dstX,dstY,dstZ);
    else if(args.BEAM == "Cone")
        cone(args,a,x,y,
                srcX,srcY,srcZ,dstX,dstY,dstZ);
    else
        assert(false);

    int64_t *ind = new int64_t[args.MAX_RAYLEN];
    float *wgt = new float[args.MAX_RAYLEN];
    int numb;
    
    forward_proj(args,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);

    float sum = 0;

    for (int i = 0; i<numb; ++i) {
        int64_t idx = ind[i];
        sum += img[idx]*wgt[i];
        idx /= args.NX*args.NY;
        if (bf[idx]>x) bf[idx] = x;
        if (br[idx]<x) br[idx] = x;
    }

    prj[x*args.NDY+y] = ushort(sum);
    delete[] ind;
    delete[] wgt;
}

const int nelipse = 10;
const float center[][2] = {
    {0,     0},
    {0,     -0.0184},
    {0.22,  0},
    {-0.22, 0},
    {0,     0.35},
    {0,     0.1},
    {0,     -0.1},
    {-0.08*2, -0.605},
    {0,     -0.605},
    {0.06*2,  -0.605}
};
const float axis[][2] = {
	{0.69,	0.92},
	{0.6624,0.874},
	{0.11,	0.31},
	{0.16,	0.41},
	{0.21,	0.25},
	{0.046,	0.046},
	{0.046,	0.046},
	{0.046*2,	0.023*2},
	{0.023*2,	0.023*2},
	{0.023*2,	0.046*2}
};
const float theta[] = {
    0,
    0,
    PI/10,
    PI/10,
    0,
    0,
    0,
    0,
    0,
    0
};
const float gray[] = {
    4,
    -2.5,
    -1.05,
    -1.05,
    1.1,
    1.1,
    1.1,
    1.1,
    1.1,
    1.1
};

void shepp_logan(const Parameter &args) {

    int *back_front = new int[args.NZ];
    int *back_rear = new int[args.NZ];

    for (int i = 0; i<args.NZ; ++i)
        back_front[i] = args.NDX+1;
    for (int i = 0; i<args.NZ; ++i)
        back_rear[i] = -1;

    boost::filesystem::path dr_file(args.RAW_DATA_FILE);

    if (boost::filesystem::exists(dr_file))
        return;

    float *img = new float[args.NX*args.NY*args.NZ];
    ushort *prj = new ushort[args.NPROJ*args.NDX*args.NDY];

    for (int k = 0; k<args.NZ; ++k) {
        for (int i = 0; i<args.NX; ++i) {
            for (int j = 0; j<args.NY; ++j) {
                float x = ((float)i*2)/args.NX - 1,
                      y = ((float)j*2)/args.NY - 1,
                      z = ((float)k*2)/args.NZ - 1;
                img[(k*args.NZ+i)*args.NX+j] = 0.0;
                for (int l = 0; l < nelipse; l++)
                {
                    const float axis1 = axis[l][1],
                                axis2 = axis[l][0],
                                axis3 = (axis1+axis2)/2;
                    const float xx = (x-center[l][1]) * cos(theta[l]) + (y-center[l][0]) * sin(theta[l]),
                                yy =-(x-center[l][1]) * sin(theta[l]) + (y-center[l][0]) * cos(theta[l]),
                                zz = z;
                    if(sqr(xx/axis1)+sqr(yy/axis2)+sqr(zz/axis3) <= 1)
                        img[(k*args.NZ+i)*args.NX+j] += gray[l];
                }
            }
        }
    }

    for(int a = 0; a < args.NPROJ; a++) {
        for (int ndx = 0; ndx<args.NDX; ++ndx) {
            for (int ndy = 0; ndy<args.NDY; ++ndy) {
                wrapper_shepp_logan(args,img,prj+a*args.NDX*args.NDY,a,ndx,ndy, back_front, back_rear);
            }
        }
    }

    ofstream fou;
    fou.open(args.RAW_DATA_FILE, ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64_t size = args.NDX * args.NDY;
    for (int k = 0; k<args.NPROJ; ++k)
        fou.write((char*)(prj+k*args.NDX*args.NDY), sizeof(ushort) * size);
    fou.close();

    //fou.open(args.PRETRACING_FILE);
    //for (int i = 0; i<args.NZ; ++i)
    //    fou<<back_front[i]<<' '<<back_rear[i]<<endl;
    //fou.close();

    delete [] back_front;
    delete [] back_rear;

    delete [] img;
    delete [] prj;
}