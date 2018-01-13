#include "testcase.h"

#include "utility.h"
#include "tracing.h"
#include "ct3d.h"

#include <cmath>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

void wrapper_abdomen(const Parameter &args, float *img, ushort *prj, int a,int x,int y, int *bf,int *br) {
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

const float sz = 20;

//C spheres
const float C_delta_x[8] = {
    2.4,
    2,
    1.6,
    1.2,
    1,
    0.8,
    0.6,
    0.4
};
const float C_z[2] = {
    0,
    3
};
const float C_y[8] = {
    1.2,
    -0.6,
    -2.1,
    -3.3,
    -4.2,
    -4.95,
    -5.55,
    -6
};
const float C_r[8] = {
    0.6,
    0.5,
    0.4,
    0.3,
    0.25,
    0.2,
    0.15,
    0.1
};
const float C_gray = 1.05;

void abdomen(const Parameter &args) {

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

    // A
    for (int k = 0; k < args.NZ; ++k)
        for (int i = 0; i < args.NX; ++i)
            for (int j = 0; j < args.NY; ++j) {
                float x = ((float)i*2)/args.NX - 1,
                      y = ((float)j*2)/args.NY - 1;
                      //z = ((float)k*2)/args.NZ - 1;
                img[(k*args.NZ+i)*args.NX+j] = 0.0;
                if (sqr(y/(20/sz))+sqr(x/(12/sz)) <= 1)
                    img[(k*args.NZ+i)*args.NX+j] = 0.8;
            }

    //B
    for (int k = 0; k < args.NZ; ++k)
        for (int i = 0; i < args.NX; ++i)
            for (int j = 0; j < args.NY; ++j) {
                float x = ((float)i*2)/args.NX - 1,
                      y = ((float)j*2)/args.NY - 1,
                      z = ((float)k*2)/args.NZ - 1;
                float xx = cos(PI/3)*(x) - sin(PI/3)*(y+10/sz),
                      yy = sin(PI/3)*(x) + cos(PI/3)*(y+10/sz),
                      zz = z;
                if (sqr(xx/(7/sz))+sqr(yy/(8/sz))+sqr(zz/(10/sz)) <= 1)
                    img[(k*args.NZ+i)*args.NX+j] = 1.0;
            }
    //Wirbelkoerper, unfinished
    for (int lev = 0; lev < 5; ++lev) {
        float xx = 0, yy = -7, zz = -10 + lev*5, l = 2;
        for (int k = 0; k < args.NZ; ++k)
            for (int i = 0; i < args.NX; ++i)
                for (int j = 0; j < args.NY; ++j) {
                    float y = ((float)i*2)/args.NX - 1,
                          x =-(((float)j*2)/args.NY - 1),
                          z = ((float)k*2)/args.NZ - 1;
                    //x= -y; y = x;
                    if (sqr(x-xx/sz)+sqr(y-yy/sz) <= sqr(1.75/sz) && fabs(z-zz/sz) <= l/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                    if (fabs(x-xx/sz) <= sqrt(1.3725)/sz && fabs(y-(yy-1.4)/sz) <= 0.1/sz && fabs(z-zz/sz) <= l/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                    if (sqr(x-xx/sz)+sqr(y-yy/sz) <= sqr(1.5/sz) && fabs(z-zz/sz) <= l/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;

                    if (fabs(x-(xx-2)/sz) <= 2/sz && fabs(y-(yy-1.975)/sz) <= 0.475/sz && fabs(z-zz/sz) <= l/sz &&
                        (x*(-0.95) + y*(4-sqrt(1.3725)) < (-0.95*(xx-4)/sz + (4-sqrt(1.3725)*(yy-2.45)/sz))) &&
                        (-2*x + 32*y) > (-2*(xx-4)/sz + 32*(yy-2.45)/sz))
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                    if (fabs(x-(xx+2)/sz) <= 2/sz && fabs(y-(yy-1.975)/sz) <= 0.475/sz && fabs(z-zz/sz) <= l/sz &&
                        (x*0.95 + y*(4-sqrt(1.3725)) < (0.95*(xx+4)/sz + (4-sqrt(1.3725)*(yy-2.45)/sz))) &&
                        (2*x + 32*y) > (2*(xx+4)/sz + 32*(yy-2.45)/sz))
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                    if (fabs(x-xx/sz) <= 0.8/sz && fabs(y-(yy-2.35)/sz) <= 0.15/sz && fabs(z-zz/sz) <= l/sz &&
                        (-2*x + 4*y) > (-2*(xx+0.4)/sz + 4*(yy-2.45)/sz) &&
                        (-2*x - 4*y) < (-2*(xx-0.4)/sz - 4*(yy-2.45)/sz))
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;

                    if (fabs(x-xx/sz) <= 0.4/sz && fabs(y-(yy-3.25)/sz) <= 0.8/sz && fabs(z-zz/sz) <= l/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                    if (sqr(x-xx/sz)+sqr(y-(yy-4.05)/sz) <= sqr(0.4/sz) && fabs(z-zz/sz) <= l/sz && y < (yy-4.05)/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                }
    }
    // C & D
    for (int id_z = 0; id_z < 1; ++id_z) {
        float zz = C_z[id_z], xx, yy, r;
        for (int lev = 0; lev < 8; ++lev) {
            yy = -14; xx = C_y[lev]; r = C_r[lev];
            for (int id = 0; id < 5; ++id) {
                if (id) yy += C_delta_x[lev];
                for (int k = 0; k < args.NZ; ++k)
                    for (int i = 0; i < args.NX; ++i)
                        for (int j = 0; j < args.NY; ++j) {
                            float x = ((float)i*2)/args.NX - 1,
                                  y = ((float)j*2)/args.NY - 1,
                                  z = ((float)k*2)/args.NZ - 1;
                            if (sqr(x-xx/sz)+sqr(y-yy/sz)+sqr(z-zz/sz) <= sqr(r/sz))
                                img[(k*args.NZ+i)*args.NX+j] = 0.8;
                        }
            }
        }
    }

    //E
    for (int k = 0; k < args.NZ; ++k)
        for (int i = 0; i < args.NX; ++i)
            for (int j = 0; j < args.NY; ++j) {
                float x = ((float)i*2)/args.NX - 1,
                      y = ((float)j*2)/args.NY - 1,
                      z = ((float)k*2)/args.NZ - 1;
                if (sqr(y/(3/sz))+sqr((x-6/sz)/(2/sz)) <= 1 && fabs(z) <= 20/sz/2) {
                    if (x <= 6.5/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 1.0;
                    else
                        img[(k*args.NZ+i)*args.NX+j] = 0.0;
                }
            }

    //F
    for (int k = 0; k < args.NZ; ++k)
        for (int i = 0; i < args.NX; ++i)
            for (int j = 0; j < args.NY; ++j) {
                float x = ((float)i*2)/args.NX - 1,
                      y = ((float)j*2)/args.NY - 1,
                      z = ((float)k*2)/args.NZ - 1;
                float xx = cos(-PI/3)*(x) - sin(-PI/3)*z,
                      yy = y+18/sz,
                      zz = sin(-PI/3)*(x) + cos(-PI/3)*z;
                if (sqr(xx/(0.4/sz))+sqr(yy/(0.4/sz)) <= 1 && fabs(zz) <= 8/sz/2)
                    img[(k*args.NZ+i)*args.NX+j] = 2.0;
                yy = y-18/sz;
                if (sqr(xx/(0.4/sz))+sqr(yy/(0.4/sz)) <= 1 && fabs(zz) <= 8/sz/2)
                    img[(k*args.NZ+i)*args.NX+j] = 2.0;
            }

    //G
    for (int k = 0; k < args.NZ; ++k)
        for (int i = 0; i < args.NX; ++i)
            for (int j = 0; j < args.NY; ++j) {
                float x = ((float)i*2)/args.NX - 1,
                      y = ((float)j*2)/args.NY - 1,
                      z = ((float)k*2)/args.NZ - 1;
                float xx, yy;
                for (int lev = 0; lev < 7; ++lev) {
                    xx = 1;
                    yy = 12-3.0/10+1.0*lev/10;
                    if (fabs(x-xx/sz) <= 0.5/sz && fabs(y-yy/sz) <= 0.25/10/sz && fabs(z) <= 1/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                }
                for (int lev = 0; lev < 7; ++lev) {
                    xx = 1;
                    yy = 14-3.0/12+1.0*lev/12;
                    if (fabs(x-xx/sz) <= 0.5/sz && fabs(y-yy/sz) <= 0.25/12/sz && fabs(z) <= 1/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                }
                for (int lev = 0; lev < 7; ++lev) {
                    yy = 12;
                    xx = -1-3.0/10+1.0*lev/10;
                    if (fabs(x-xx/sz) <= 0.25/10/sz && fabs(y-yy/sz) <= 0.5/sz && fabs(z) <= 1/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                }
                for (int lev = 0; lev < 7; ++lev) {
                    yy = 14;
                    xx = -1-3.0/12+1.0*lev/12;
                    if (fabs(x-xx/sz) <= 0.25/12/sz && fabs(y-yy/sz) <= 0.5/sz && fabs(z) <= 1/sz)
                        img[(k*args.NZ+i)*args.NX+j] = 2.0;
                }
            }

    for(int a = 0; a < args.NPROJ; a++) {
        for (int ndx = 0; ndx<args.NDX; ++ndx) {
            for (int ndy = 0; ndy<args.NDY; ++ndy) {
                wrapper_abdomen(args,img,prj+a*args.NDX*args.NDY,a,ndx,ndy, back_front, back_rear);
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

    // fou.open(args.PRETRACING_FILE);
    // for (int i = 0; i<args.NZ; ++i)
    //    fou<<back_front[i]<<' '<<back_rear[i]<<endl;
    // fou.close();

    delete [] back_front;
    delete [] back_rear;

    delete [] img;
    delete [] prj;
}