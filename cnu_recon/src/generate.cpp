#include <cmath>
#include <fstream>
#include <iostream>

#include "tracing.h"
#include "ct3d.h"
#include "utility.h"

using namespace std;

typedef unsigned short ushort;
typedef long long int64;

inline float sqr(float x) {
    return x*x;
}

inline bool equals(float x,int y) {
    return fabs(x-y)<1.25;
}

void wrapper1(float *img, ushort *prj, int a,int x,int y) {
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    srcX = y+0.5-64.0;
    srcZ = 128-x-1;
    srcY = -450;

    
    dstX = y+0.5-64.0;
    dstZ = 128-x-1;
    dstY = 450;

    float theta = (float)a/360*2*PI;

    rotate_2d(srcX, srcY, theta);
    rotate_2d(dstX, dstY, theta);

    int64_t *ind = new int64_t[512];
    float *wgt = new float[512];
    int numb;

    srcX += 64.0;
    srcY += 64.0;
    dstX += 64.0;
    dstY += 64.0;

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
    delete[] ind;
    delete[] wgt;
}

void test1() {
    float *img = new float[128*128*128];
    ushort *prj = new ushort[128*128];

    for (int k = 0; k<128; ++k) {
        for (int i = 0; i<128; ++i)
            for (int j = 0; j<128; ++j) {
                float dst = sqrt(sqr(i-64.0)+sqr(j-64.0));
                if (equals(dst,24) || equals(dst,48) || equals(dst,36) || equals(dst,12) || equals(dst,60)) {
                    img[k*128*128+i*128+j] = 0.2;
                }
                else img[k*128*128+i*128+j] = 0.0;
            }
    }

    //ofstream pfile("test.txt");
    for (int ndx = 0; ndx<128; ++ndx) {
        for (int ndy = 0; ndy<128; ++ndy) {
            wrapper1(img,prj,0,ndx,ndy);
            //pfile<<prj[ndx*128+ndy]<<' ';
        }
        //pfile<<endl;
    }
    //pfile.close();

    ofstream fou("test_360x128x128_128_parallel_circle.dr", ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64_t size = (int64)128 * 128;
    for (int k = 0; k<360; ++k)
        fou.write((char*)prj, sizeof(ushort) * size);
    fou.close();

    delete [] img;
    delete [] prj;
}

void test2() {
    float *img = new float[128*128*128];
    ushort *prj = new ushort[128*128];

    for (int k = 0; k<128; ++k) {
        for (int i = 0; i<128; ++i)
            for (int j = 0; j<128; ++j) {
                float dst = sqrt(sqr(i-64.0)+sqr(j-64.0)+sqr(k-64.0));
                if (equals(dst,24) || equals(dst,48) || equals(dst,36) || equals(dst,12) || equals(dst,60)) {
                    img[k*128*128+i*128+j] = 0.2;
                }
                else img[k*128*128+i*128+j] = 0.0;
            }
    }

    //ofstream pfile("test.txt");
    for (int ndx = 0; ndx<128; ++ndx) {
        for (int ndy = 0; ndy<128; ++ndy) {
            wrapper1(img,prj,0,ndx,ndy);
            //pfile<<prj[ndx*128+ndy]<<' ';
        }
        //pfile<<endl;
    }
    //pfile.close();

    ofstream fou("test_360x128x128_128_parallel_sphere.dr", ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64_t size = (int64)128 * 128;
    for (int k = 0; k<360; ++k)
        fou.write((char*)prj, sizeof(ushort) * size);
    fou.close();

    delete [] img;
    delete [] prj;
}

void test3() {
    float *img = new float[128*128*128];
    ushort *prj = new ushort[360*128*128];

    for (int k = 0; k<128; ++k) {
        for (int i = 0; i<128; ++i)
            for (int j = 0; j<128; ++j) {
                float dst = sqrt(sqr(1.2*(i-64.0))+sqr(0.8*(j-64.0))+sqr(k-64.0));
                if (equals(dst,24) || equals(dst,36) || equals(dst,12) || equals(dst,48)) {
                    img[k*128*128+i*128+j] = 0.2;
                }
                else img[k*128*128+i*128+j] = 0.0;
            }
    }
    //ofstream pfile("test.txt");
    for(int a = 0; a<360; ++a) {
        for (int ndx = 0; ndx<128; ++ndx) {
            for (int ndy = 0; ndy<128; ++ndy) {
                wrapper1(img,prj+a*128*128,a,ndx,ndy);
            }
            //pfile<<prj[ndx*128+ndy]<<' ';
        }
        //pfile<<endl;
    }
    //pfile.close();

    ofstream fou("test_360x128x128_128_parallel_ellipse.dr", ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64_t size = (int64)128 * 128;
    for (int k = 0; k<360; ++k)
        fou.write((char*)(&prj[k*128*128]), sizeof(ushort) * size);
    fou.close();

/*
    for(int z = 0; z < 180; z++)
    {
        char buf[256];
        snprintf(buf, 256, "origin_img_128x128x128/%d", z);
        fou.open(buf, ios::out);
        fou.precision(6);
        for(int x = 0; x < 128; x++)
        {
            for(int y = 0; y < 128; y++)
            {
                fou << prj[z*128*128+x*128+y] << ' ';
            }
            fou << endl;
        }
        fou.close();
    }
    for(int z = 0; z < 360; z++)
    {
        char buf[256];
        snprintf(buf, 256, "origin_img_128x128x128_prj/%d", z);
        fou.open(buf, ios::out);
        fou.precision(6);
        for(int x = 0; x < 128; x++)
        {
            for(int y = 0; y < 128; y++)
            {
                fou << prj[z*128*128+x*128+y] << ' ';
            }
            fou << endl;
        }
        fou.close();
    }
*/

    delete [] img;
    delete [] prj;
}

void test4() {
    float *img = new float[128*128*128];
    ushort *prj = new ushort[360*128*128];

    for (int k = 0; k<128; ++k) {
        for (int i = 0; i<128; ++i)
            for (int j = 0; j<128; ++j) {
                const float sa = sin((float)k/128*2*PI), ca = cos((float)k/128*2*PI);
                const float x = (i-64.0)*ca-(j-64.0)*sa, y = (j-64.0)*ca+(i-64.0)*sa;
                const float dst1 = sqrt(sqr(x+32.0)+sqr(y)),
                            dst2 = sqrt(sqr(x-32.0)+sqr(y)),
                            dst3 = sqrt(sqr(x)+sqr(y));
                if((dst1 < 16) || (dst2 > 16 && dst2 < 32) || (dst1 > 32 && dst2 > 32 && y >= 0 && dst3 < 64))
                    img[k*128*128+i*128+j] = 0.2;
                else
                    img[k*128*128+i*128+j] = 0.0;

            }
    }
    //ofstream pfile("test.txt");
    for(int a = 0; a<360; ++a) {
        for (int ndx = 0; ndx<128; ++ndx) {
            for (int ndy = 0; ndy<128; ++ndy) {
                wrapper1(img,prj+a*128*128,a,ndx,ndy);
            }
            //pfile<<prj[ndx*128+ndy]<<' ';
        }
        //pfile<<endl;
    }
    //pfile.close();

    ofstream fou("test_360x128x128_128_parallel_taiji.dr", ios::binary);
    char empty[1024] = {};
    fou.write(empty,1024);
    int64_t size = (int64)128 * 128;
    for (int k = 0; k<360; ++k)
        fou.write((char*)(&prj[k*128*128]), sizeof(ushort) * size);
    fou.close();
    for(int z = 0; z < 360; z++)
    {
        char buf[256];
        snprintf(buf, 256, "origin_img_128x128x128/%d", z);
        fou.open(buf, ios::out);
        fou.precision(6);
        for(int x = 0; x < 128; x++)
        {
            for(int y = 0; y < 128; y++)
            {
                fou << prj[z*128*128+x*128+y] << ' ';
            }
            fou << endl;
        }
        fou.close();
    }
    delete [] img;
    delete [] prj;
}

int main(int argc, char** argv) {
    test1();
    test2();
    test3();
    test4();
    return 0;
}
