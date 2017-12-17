#include "ct3d.h"
#include "utility.h"
#include "tracing.h"

#include <iostream>

#include <omp.h>

#include <cmath>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace std;
using boost::format;

int max_numb;
float ave_numb;

float Afg;

/*
 Z-axis         
 |   /    
 |  / |   
 | /  |  detector  board
 | |  |   
 | |  /      <--   object   <--  source
 | | /   /
 | |/   /  x-axis
 |     /
 |    /
 |   /
 |  /
 | /
 |/               y-axis (negative direction)
 ----------------------->

step 1: the center of the object is at the origin (0,0)
        then we calculate src and dst using rotate
step 2: shift src and dst to right place
*/

float minIMAGE(float Af, int64_t *line, float *weight, int numb, float lambda,
               const Parameter &args, const CTInput &in, CTOutput &out);
void minEDGE(int64_t *line, float *weight, int numb, float lambda,
             const Parameter &args, const CTInput &in, CTOutput &out);

void compute(float lambda,float *sin_table,float *cos_table,
             int alpha, int detectorX, int detectorY,
             const Parameter &args, const CTInput &in, CTOutput &out) {
    //const float sina = sin_table[alpha], cosa = cos_table[alpha];
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    if (args.BEAM=="Parallel")
        parallel(args,alpha,detectorX,detectorY,
                 srcX,srcY,srcZ,dstX,dstY,dstZ);
    else if (args.BEAM=="Cone")
        cone(args,alpha,detectorX,detectorY,
             srcX,srcY,srcZ,dstX,dstY,dstZ);

    //cout << boost::format("%1%  :  %2% %3% %4%     %5% %6% %7%") %alpha%srcX%srcY%srcZ%dstX%dstY%dstZ <<endl;

    int64_t *ind = new int64_t[args.MAX_RAYLEN];
    float *wgt = new float[args.MAX_RAYLEN];
    int numb;

    float Af = -in.sino_data(alpha, detectorX, detectorY);

    forward_proj(args,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);
   
    if (max_numb<numb) max_numb = numb;
    ave_numb += numb;

    minIMAGE(Af, ind, wgt, numb, lambda, args,in,out);
    minEDGE(ind, wgt, numb, lambda, args,in,out);

    delete[] ind;
    delete[] wgt;
}

/*
void Akernel(sino_type &Af, img_type *f,
             int alpha, int detectorX, int detectorY,
             const Parameter &args, const CTInput &in, CTOutput &out) {
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    if (args.BEAM=="Parallel")
        parallel(args,alpha,detectorX,detectorY,
                 srcX,srcY,srcZ,dstX,dstY,dstZ);
    else if (args.BEAM=="Cone")
        cone(args,alpha,detectorX,detectorY,
             srcX,srcY,srcZ,dstX,dstY,dstZ);

    //cout << boost::format("%1%  :  %2% %3% %4%     %5% %6% %7%") %alpha%srcX%srcY%srcZ%dstX%dstY%dstZ <<endl;

    int64_t *ind = new int64_t[args.MAX_RAYLEN];
    float *wgt = new float[args.MAX_RAYLEN];
    int numb;

    float Af = -in.sino_data(alpha, detectorX, detectorY);

    forward_proj(args,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);

    Af = 0;
    for (int i = 0; i<numb; ++i) {
        int64_t idx = ind[i];
        idx %= args.NX*args.NY;
        int64_t nx,ny;
        nx = idx/args.NY;
        ny = idx%args.NY;
        Af += f[(nx+1)*(args.NY+2)+(ny+1)]*wgt[i];
    }

    delete[] ind;
    delete[] wgt;
}

void AStarkernel(sino_type &Af, img_type *f,
                 int alpha, int detectorX, int detectorY,
                 const Parameter &args, const CTInput &in, CTOutput &out) {
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    if (args.BEAM=="Parallel")
        parallel(args,alpha,detectorX,detectorY,
                 srcX,srcY,srcZ,dstX,dstY,dstZ);
    else if (args.BEAM=="Cone")
        cone(args,alpha,detectorX,detectorY,
             srcX,srcY,srcZ,dstX,dstY,dstZ);

    //cout << boost::format("%1%  :  %2% %3% %4%     %5% %6% %7%") %alpha%srcX%srcY%srcZ%dstX%dstY%dstZ <<endl;

    int64_t *ind = new int64_t[args.MAX_RAYLEN];
    float *wgt = new float[args.MAX_RAYLEN];
    int numb;

    float Af = -in.sino_data(alpha, detectorX, detectorY);

    forward_proj(args,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);

    for (int i = 0; i<numb; ++i) {
        int64_t idx = ind[i];
        idx %= args.NX*args.NY;
        int64_t nx,ny;
        nx = idx/args.NY;
        ny = idx%args.NY;
        f[(nx+1)*(args.NY+2)+(ny+1)] += Af*wgt[i];
    }

    delete[] ind;
    delete[] wgt;
}

void A(sino_type *Af, img_type *f,
       int np, int _back_front, int _back_rear,
       const Parameter &args, const CTInput &in, CTOutput &out) {
    for (int ndx = _back_front; ndx<=_back_rear; ++ndx)
        for (int ndy = 0; ndy<args.NDY; ++ndy) {
            Akernel(Af[(ndx-_back_front)*args.NDY+ndy], f, np,ndx,ndy, args,in,out);
        }
}

void AStar(sino_type *Af, img_type *f,
           int np, int _back_front, int _back_rear,
           const Parameter &args, const CTInput &in, CTOutput &out) {
    for (int ndx = _back_front; ndx<=_back_rear; ++ndx)
        for (int ndy = 0; ndy<args.NDY; ++ndy) {
            AStarkernel(Af[(ndx-_back_front)*args.NDY+ndy], f, np,ndx,ndy, args,in,out);
        }
}

void compute(float lambda,float *sin_table,float *cos_table,
             int np, int z,
             const Parameter &args, const CTInput &in, CTOutput &out) {
    //const float sina = sin_table[alpha], cosa = cos_table[alpha];

    _back_front = in.back_front[z];
    _back_rear = in.back_rear[z];

    img_type *f = new img_type[(args.NX+2)*(args.NY+2)];
    edge_type *v = new edge_type[(args.NX+2)*(args.NY+2)];
    sino_type *g = new sino_type[(_back_rear-_back_front+1)*args.NDY];

    for (int i = 0; i<(args.NX+2)*(args.NY+2); ++i) {
        f[i] = 0.;
        v[i] = 1.0;
    }

    img_type *__f = &out.img[z*args.NX*args.NY];
    edge_type *__v = &out.edge[z*args.NX*args.NY];
    for (int i = 1; i<=args.NX; ++i)
        for (int j = 1; j<=args.NY; ++i) {
            f[i*(args.NY+2)+j] = __f[(i-1)*args.NY+(j-1)];
            v[i*(args.NY+2)+j] = __v[(i-1)*args.NY+(j-1)];
        }
    for (int ndx = _back_front; ndx<=_back_rear; ++ndx)
        for (int ndy = 0; ndy<args.NDY; ++ndy)
            g[(ndx-_back_front)*args.NDY+ndy] = in.sino_data(np,ndx,ndy);

    // IMAGE first

    // r = g-Af
    sino_type *Af = new sino_type[(_back_rear-_back_front+1)*args.NDY];
    sino_type *r = new sino_type[(_back_rear-_back_front+1)*args.NDY];
    A(Af, f, np, _back_front, _back_rear, args,in,out);
    for (int i = 0; i<(_back_rear-_back_front+1)*args.NDY; ++i)
        r[i] = g[i]-Af[i];

    // d = A*(r)
    img_type *d = new img_type[(args.NX+2)*(args.NY+2)];
    for (int i = 0; i<(args.NX+2)*(args.NY+2); ++i) d[i] = 0;
    AStar(r, d, np, _back_front, _back_rear, args,in,out);

    // d += alpha*(v^2 \nabla f) + e^2 \laplace f
    for (int i = 1; i<=args.NX; ++i)
        for (int j = 1; j<=args.NY; ++j) {
            d[i*(args.NY+2)] += args.ALPHA*(
                +sqr(v[i    *(args.NY+2)+j    ])*(f[(i+1)*(args.NY+2)+j    ]-f[i    *(args.NY+2)+j    ])
                +sqr(v[i    *(args.NY+2)+j    ])*(f[i    *(args.NY+2)+(j+1)]-f[i    *(args.NY+2)+j    ])
                -sqr(v[(i-1)*(args.NY+2)+j    ])*(f[i    *(args.NY+2)+j    ]-f[(i-1)*(args.NY+2)+j    ])
                -sqr(v[i    *(args.NY+2)+(j-1)])*(f[i    *(args.NY+2)+j    ]-f[i    *(args.NY+2)+(j-1)])
                +0.01*sqr(args.EPSILON)*(
                    +  f[(i+1)*(args.NY+2)+j    ]
                    +  f[i    *(args.NY+2)+(j+1)]
                    -4*f[i    *(args.NY+2)+j    ]
                    +  f[(i-1)*(args.NY+2)+j    ]
                    +  f[i    *(args.NY+2)+(j-1)]
            );
        }

    img_type *p = new img_type[(args.NX+2)*(args.NY+2)];
    sino_type *q = new sino_type[(_back_rear-_back_front+1)*args.NDY];

    // p = d
    for (int i = 0; i<(args.NX+2)*(args.NY+2); ++i) p[i] = d[i];

    A(q, p, np, _back_front, _back_rear, args,in,out);

    float nV_KepsNablap = 0;
    for (int i = 1; i<=args.NX; ++i)
        for (int j = 1; j<=args.NY; ++j) {
            nV_KepsNablap += (sqr(v[i*(args.NY+2)+j])+0.01*sqr(args.EPSILON))*(
                sqr(p[i*(args.NY+2)+j]-p[(i-1)*(args.NY+2)+j])
                +sqr(p[i*(args.NY+2)+j]-p[i*(args.NY+2)+(j-1)])
            );
        }
    float pSkalard = 0;
    for (int i = 1; i<=args.NX; ++i)
        for (int j = 1; j<=args.NY; ++j)
            pSkalard += p[i*(args.NY+2)+j]*d[i*(args.NY+2)+j];
    float normAp = 0;
    for (int i = 0; i<(_back_rear-_back_front+1)*args.NDY; ++i)
        normAp += sqr(q[i]);
    float c_1 = pSkalard/(normAp+args.ALPHA*nV_KepsNablap);

    // update f
    for (int i = 1; i<=args.NX; ++i)
        for (int j = 1; j<=args.NY; ++j) {
            f[i*(args.NY+2)+j] += c_1 * p[i*(args.NY+2)+j];



    // write back
    for (int i = 1; i<=args.NX; ++i)
        for (int j = 1; j<=args.NY; ++i) {
            __f[(i-1)*args.NY+(j-1)] = f[i*(args.NY+2)+j];
            __v[(i-1)*args.NY+(j-1)] = v[i*(args.NY+2)+j];
        }

    delete [] p;
    delete [] q;

    delete [] Af;
    delete [] r;
    delete [] d;

    delete [] f;
    delete [] v;
    delete [] g;
}
*/

void wrapper(float lambda,float *sin_table,float *cos_table,
             const Parameter &args,const CTInput &in,CTOutput &out) {
    for (int k = 0; k < args.NPROJ; ++k)
        for (int i = 0; i < args.NDX; ++i)
            #pragma omp parallel for
            for (int j = 0; j < args.NDY; ++j) {
                compute(lambda, sin_table, cos_table, k, i, j, args,in,out);
                //cout << format("%1% %2% %3%") %k%i%j <<endl;
            }
}

/*
void wrapper(float lambda,float *sin_table,float *cos_table,
             const Parameter &args,const CTInput &in,CTOutput &out) {
    for (int np = 0; np < args.NPROJ; ++np) {
        for (int slice = args.NZ-1; slice>=0; --slice) {
            //for (int i = 0; i < args.NDX; ++i) {
            //#pragma omp parallel for
            compute(lambda, sin_table, cos_table, np, slice, args,in,out);
            //cout << format("%1% %2% %3%") %k%i%j <<endl;
    }
}
*/
void ct3d(const Parameter &args,const CTInput &in,CTOutput &out) {
    float *cos_table = new float[args.NPROJ];
    float *sin_table = new float[args.NPROJ];

    for(int i = 0; i < args.NPROJ; i++) {
        const float alpha = (float)i/args.NPROJ*2*PI;
        sin_table[i] = sin(alpha);
        cos_table[i] = cos(alpha);
    }

    max_numb = 0;

    out.allocate_img(args.NX*args.NY*args.NZ);
    out.allocate_edge(args.NX*args.NY*args.NZ);

    float lambda = 0.01;

    int64_t global_start = timer_s();

    for (int iters = 0; iters < args.ITERATIONS; ++iters) {
        lambda = 1.0/(100.0+iters*2.0);

        cout << format("iter = %1%, lambda = %2%") % iters % lambda <<endl;

        ave_numb = 0;
        Afg = 0;

        int64_t start = timer_s();
        wrapper(lambda,sin_table,cos_table, args,in,out);
        int64_t end = timer_s();
        
        ave_numb /= args.NPROJ*args.NDX*args.NDY;
        cout << format("||Af-g||^2 = %1%") % Afg <<endl;
        cout << format("time used = %1% seconds") % (end-start) <<endl;
        cout << "-------------------------------------------" <<endl;
    }

    int64_t global_end = timer_s();
    
    cout << "===========================================" <<endl;
    cout << format("TOTAL time used = %1% seconds") % (global_end-global_start) <<endl;
    cout<< format("actual max raylen = %1%, average raylen = %2%") % max_numb % ave_numb << endl;

    delete [] sin_table;
    delete [] cos_table;
}

int main(int argc, char** argv) {
    Parameter args;

    args.parse_config(argc, argv);
    args.print_options();

    omp_set_dynamic(0);
    omp_set_num_threads(args.THREAD_NUMB);

    CTInput in = CTInput(args);
    CTOutput out = CTOutput(args);

    in.read_sino(args.RAW_DATA_FILE);
    in.read_pretracing(args.PRETRACING_FILE);
    ct3d(args,in,out);
    out.write_img(args.OUTPUT_DIR);
    out.write_edge(args.OUTPUT_DIR);

    return 0;
}

float minIMAGE(float Af, int64_t *line, float *weight, int numb, float lambda,
               const Parameter &args, const CTInput &in, CTOutput &out) {
    float *d = new float[args.MAX_RAYLEN];
    for (int i = 0; i<numb; ++i) {
        Af += out.img[line[i]] * weight[i];
    }
    for (int i = 0; i<numb; ++i) {
        int64_t ind = line[i];
        int64_t plain = ind%(args.NX*args.NY);
        int x = plain/args.NY, y = plain%args.NY;

        float tmp = 0.;
        float lap = 0.;
        
        if (x+1<args.NX) tmp += sqr(out.edge[ind])*(out.img[ind+args.NY]-out.img[ind]);
        else             tmp += sqr(out.edge[ind])*(          0         -out.img[ind]);
        
        if (y+1<args.NY) tmp += sqr(out.edge[ind])*(out.img[ind+1]-out.img[ind]);
        else             tmp += sqr(out.edge[ind])*(        0     -out.img[ind]);
        
        if (x-1>=0)      tmp -= sqr(out.edge[ind-args.NY])*(out.img[ind]-out.img[ind-args.NY]);
        else             tmp -=                            (out.img[ind]-           0        );
        
        if (y-1>=0)      tmp -= sqr(out.edge[ind-1])*(out.img[ind]-out.img[ind-1]);
        else             tmp -=                      (out.img[ind]-      0       );
        
        if (x+1<args.NX) lap += out.img[ind+args.NY];
        if (y+1<args.NY) lap += out.img[ind+1];
        if (x-1>=0)      lap += out.img[ind-args.NY];
        if (y-1>=0)      lap += out.img[ind-1];
        lap -= 4*out.img[ind];

        d[i] = -Af*weight[i]+args.ALPHA*(tmp+sqr(args.EPSILON)*lap);
    }
    for (int i = 0; i<numb; i++) {
        int64_t ind = line[i];
        img_type tmp = out.img[ind] + lambda * d[i];
        if (tmp<0) tmp = 0;
        out.img[ind] = tmp;
    }
    delete [] d;
    Afg += sqr(Af);
}

void minEDGE(int64_t *line, float *weight, int numb, float lambda,
             const Parameter &args, const CTInput &in, CTOutput &out) {
    float *d = new float[args.MAX_RAYLEN];
    for (int i = 0; i<numb; ++i) {
        int64_t ind = line[i];
        int64_t plain = ind%(args.NX*args.NY);
        int x = plain/args.NY, y = plain%args.NY;

        float a = 0.;
        float b = 0.;
        float c = 0.;

        if (x-1>=0)      a += sqr(out.img[ind]-out.img[ind-args.NY]);
        else             a += sqr(out.img[ind]-           0        );
        
        if (y-1>=0)      a += sqr(out.img[ind]-out.img[ind-1]);
        else             a += sqr(out.img[ind]-      0       );
        
        a *= out.edge[ind];
        
        b = out.edge[ind]-1;

        if (x+1<args.NX) c += out.edge[ind+args.NY];
        if (y+1<args.NY) c += out.edge[ind+1];
        if (x-1>=0)      c += out.edge[ind-args.NY];
        if (y-1>=0)      c += out.edge[ind-1];
        c -= 4*out.edge[ind];
        
        d[i] = -args.ALPHA*a-args.BETA/(4*args.EPSILON)*b+args.BETA*args.EPSILON*c;
    }
    for (int i = 0; i<numb; i++) {
        int64_t ind = line[i];
        edge_type tmp = out.edge[ind] + lambda * d[i];
        if (tmp<0) tmp = 0;
        if (tmp>1) tmp = 1;
        out.edge[ind] = tmp;
    }
    delete [] d;
}
