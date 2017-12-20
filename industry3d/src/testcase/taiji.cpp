#include "testcase.h"

#include "utility.h"
#include "tracing.h"
#include "ct3d.h"

#include <cmath>
#include <fstream>

#include <boost/filesystem.hpp>

using namespace std;

void wrapper_taiji(const Parameter &args, float *img, ushort *prj, int a,int x,int y, int *bf,int *br) {
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

void taiji(const Parameter &args) {

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

    float raidus = 0.75;//control the size of taiji, raidus \in [0,1]

    for (int k = 0; k<args.NZ; ++k) {
        for (int i = 0; i<args.NX; ++i) {
            for (int j = 0; j<args.NY; ++j) {
                const float zz = ((float)k*2/args.NZ - 1) / raidus;
                const float sa = sin(zz*PI), ca = cos(zz*PI);
                const float xx = ((float)i*2/args.NX - 1)/raidus,
                            yy = ((float)j*2/args.NY - 1)/raidus;
                const float x = xx * ca - yy * sa,
                            y = yy * ca + xx * sa;
                const float dst1 = sqrt(sqr(x + 0.5)+sqr(y)),
                            dst2 = sqrt(sqr(x - 0.5)+sqr(y)),
                            dst3 = sqrt(sqr(x)+sqr(y));
                if( (dst1 < 0.25) ||
                    (dst2 > 0.25 && dst2 < 0.5) ||
                    (dst1 > 0.5 && dst2 > 0.5 && y >= 0 && dst3 < 1))
                    img[k*args.NX*args.NY+i*args.NY+j] = 10;
                else img[k*args.NX*args.NY+i*args.NY+j] = 0.0;
            }
        }
    }

    for(int a = 0; a < args.NPROJ; a++) {
        for (int ndx = 0; ndx<args.NDX; ++ndx) {
            for (int ndy = 0; ndy<args.NDY; ++ndy) {
                wrapper_taiji(args,img,prj+a*args.NDX*args.NDY,a,ndx,ndy, back_front, back_rear);
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

    //fou.open(args.PRETRACING_FILE, ios::out|ios::trunc);
    //for (int i = 0; i<args.NZ; ++i)
    //    fou<<back_front[i]<<' '<<back_rear[i]<<endl;
    //fou.close();

    delete [] back_front;
    delete [] back_rear;

    delete [] img;
    delete [] prj;
}