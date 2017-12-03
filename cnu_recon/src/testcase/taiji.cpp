#include "testcase.h"

#include "utility.h"
#include "tracing.h"
#include "ct3d.h"

#include <cmath>
#include <fstream>

#include <boost/filesystem.hpp>

using namespace std;

void wrapper_taiji(const Parameter &args, float *img, ushort *prj, int a,int x,int y) {
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
        sum += img[ind[i]]*wgt[i];
    }

    prj[x*args.NDY+y] = ushort(sum);
    delete[] ind;
    delete[] wgt;
}

void taiji(const Parameter &args) {

    boost::filesystem::path dr_file(args.RAW_DATA_FILE);

    if (boost::filesystem::exists(dr_file))
        return;

    float *img = new float[args.NX*args.NY*args.NZ];
    ushort *prj = new ushort[args.NPROJ*args.NDX*args.NDY];

    for (int k = 0; k<args.NZ; ++k) {
        for (int i = 0; i<args.NX; ++i) {
            for (int j = 0; j<args.NY; ++j) {
                const float sa = sin((float)k/args.NZ*2*PI), ca = cos((float)k/args.NZ*2*PI);
                const float x = (i-args.HALFSIZE)*ca-(j-args.HALFSIZE)*sa,
                            y = (j-args.HALFSIZE)*ca+(i-args.HALFSIZE)*sa;
                const float dst1 = sqrt(sqr(x+args.HALFSIZE/2)+sqr(y)),
                            dst2 = sqrt(sqr(x-args.HALFSIZE/2)+sqr(y)),
                            dst3 = sqrt(sqr(x)+sqr(y));
                if( (dst1 < args.HALFSIZE/4) ||
                    (dst2 > args.HALFSIZE/4 && dst2 < args.HALFSIZE/2) ||
                    (dst1 > args.HALFSIZE/2 && dst2 > args.HALFSIZE/2 && y >= 0 && dst3 < args.HALFSIZE))
                    img[k*args.NX*args.NY+i*args.NY+j] = 0.2;
                else img[k*args.NX*args.NY+i*args.NY+j] = 0.0;
            }
        }
    }

    for(int a = 0; a < args.NPROJ; a++) {
        for (int ndx = 0; ndx<args.NDX; ++ndx) {
            for (int ndy = 0; ndy<args.NDY; ++ndy) {
                wrapper_taiji(args,img,prj+a*args.NDX*args.NDY,a,ndx,ndy);
            }
        }
    }

    ofstream fou(args.RAW_DATA_FILE, ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64_t size = args.NDX * args.NDY;
    for (int k = 0; k<args.NPROJ; ++k)
        fou.write((char*)(prj+k*args.NDX*args.NDY), sizeof(ushort) * size);
    fou.close();

    delete [] img;
    delete [] prj;
}