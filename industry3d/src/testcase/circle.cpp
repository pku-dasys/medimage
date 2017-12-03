#include "testcase.h"

#include "utility.h"
#include "tracing.h"
#include "ct3d.h"

#include <cmath>
#include <fstream>

#include <boost/filesystem.hpp>

using namespace std;

void wrapper_parallel_circle(const Parameter &args, float *img, ushort *prj, int a,int x,int y) {
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    parallel(args,a,x,y,
            srcX,srcY,srcZ,dstX,dstY,dstZ);

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

void parallel_circle(const Parameter &args) {

    boost::filesystem::path dr_file(args.RAW_DATA_FILE);

    if (boost::filesystem::exists(dr_file))
        return;

    float *img = new float[args.NX*args.NY*args.NZ];
    ushort *prj = new ushort[args.NDX*args.NDY];

    for (int k = 0; k<args.NZ; ++k) {
        for (int i = 0; i<args.NX; ++i)
            for (int j = 0; j<args.NY; ++j) {
                float dst = sqrt(sqr(i-args.HALFSIZE)+sqr(j-args.HALFSIZE));
                if (equals_draw(dst,24) || equals_draw(dst,48) || equals_draw(dst,36) || equals_draw(dst,12) || equals_draw(dst,60)) {
                    img[k*args.NX*args.NY+i*args.NY+j] = 0.2;
                }
                else img[k*args.NX*args.NY+i*args.NY+j] = 0.0;
            }
    }

    for (int ndx = 0; ndx<args.NDX; ++ndx) {
        for (int ndy = 0; ndy<args.NDY; ++ndy) {
            wrapper_parallel_circle(args,img,prj,0,ndx,ndy);
        }
    }

    ofstream fou(args.RAW_DATA_FILE, ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64_t size = args.NDX * args.NDY;
    for (int k = 0; k<args.NPROJ; ++k)
        fou.write((char*)prj, sizeof(ushort) * size);
    fou.close();

    delete [] img;
    delete [] prj;
}