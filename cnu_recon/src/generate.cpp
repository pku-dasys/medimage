#include <cmath>
#include <fstream>

#include "tracing.h"
#include "main.h"

using namespace std;

typedef unsigned short ushort;
typedef long long int64;

float img[128*128*128];
ushort prj[360*128*128];

inline float sqr(float x) {
    return x*x;
}

inline bool equals(float x,int y) {
    return fabs(x-y)<2;
}

void wrapper(int a,int x,int y) {
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    srcX = x;
    srcY = y;
    srcZ = 0.0;

    
    dstX = x;
    dstY = y;
    dstZ = 256;

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

    prj[a*128*128+x*128+y] = ushort(sum);
}

int main(int argc, char** argv) {
    for (int k = 0; k<128; ++k) {
        for (int i = 0; i<128; ++i)
            for (int j = 0; j<128; ++j) {
                float dst = sqrt(sqr(i-64.0)+sqr(j-64.0));
                if (equals(dst,8) || equals(dst,32) || equals(dst,64) || equals(dst,96) || equals(dst,16)) {
                    img[k*128*128+i*128+j] = 0.2;
                }
                else img[k*128*128+i*128+j] = 0.0;
            }
    }

    for (int alpha = 0; alpha<360; ++alpha)
        for (int ndx = 0; ndx<128; ++ndx)
            for (int ndy = 0; ndy<128; ++ndy) {
                wrapper(alpha,ndx,ndy);
            }
    
    ofstream fou("test.dr", ios::out | ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64 size = (int64)360 * 128 * 128;
    fou.write((char*)prj, sizeof(ushort) * size);
    fou.close();

    return 0;
}
