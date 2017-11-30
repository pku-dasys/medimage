#include <cmath>
#include <fstream>

#include "tracing.h"
#include "main.h"

using namespace std;

typedef unsigned short ushort;
typedef long long int64;

float img[128*128*128];
ushort prj[128*128];

inline float sqr(float x) {
    return x*x;
}

inline bool equals(float x,int y) {
    return fabs(x-y)<1.25;
}

void wrapper(int a,int x,int y) {
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    srcX = y;
    srcZ = 128-x-1;
    srcY = 450;

    
    dstX = y;
    dstZ = 128-x-1;
    dstY = -450;

    int64 *ind = new int64[512];
    float *wgt = new float[512];
    int numb;

    Parameter args;
    args.NX = 128;
    args.NY = 128;
    args.NZ = 128;

    forward_proj(args,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);

    float sum = 0;

    for (int i = 0; i<numb; ++i) {
        sum += img[ind[i]]*wgt[i];
    }

    prj[x*128+y] = ushort(sum);
}

int main(int argc, char** argv) {
    for (int k = 0; k<128; ++k) {
        for (int i = 0; i<128; ++i)
            for (int j = 0; j<128; ++j) {
                float dst = sqrt(sqr(i-64.0)+sqr(j-64.0));
                if (equals(dst,24) || equals(dst,48)) {
                    img[k*128*128+i*128+j] = 0.2;
                }
                else img[k*128*128+i*128+j] = 0.0;
            }
    }

    ofstream pfile("test.txt");
        for (int ndx = 0; ndx<128; ++ndx) {
            for (int ndy = 0; ndy<128; ++ndy) {
                wrapper(0,ndx,ndy);
                pfile<<prj[ndx*128+ndy]<<' ';
            }
            pfile<<endl;
        }
    pfile.close();

    ofstream fou("test.dr", ios::out | ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64 size = (int64)128 * 128;
    for (int k = 0; k<360; ++k)
        fou.write((char*)prj, sizeof(ushort) * size);
    fou.close();

    return 0;
}