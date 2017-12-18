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

void gdIMAGE(float lambda, int np, int ndx,
              int64_t *ind, float *wgt, int *numb,
              const Parameter &args, const CTInput &in, CTOutput &out);

void gdEDGE(float lambda,
             int64_t *ind, float *wgt, int *numb,
             const Parameter &args, const CTInput &in, CTOutput &out);

void get_tracing(int alpha,int detectorX,int detectorY,
                 int64_t *ind,float *wgt,int &numb,
                 const Parameter &args) {
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    if (args.BEAM=="Parallel")
        parallel(args,alpha,detectorX,detectorY,
                 srcX,srcY,srcZ,dstX,dstY,dstZ);
    else if (args.BEAM=="Cone")
        cone(args,alpha,detectorX,detectorY,
             srcX,srcY,srcZ,dstX,dstY,dstZ);

    forward_proj(args,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);
}

void compute_single(float lambda, int alpha, int detectorX, int detectorY,
                    const Parameter &args, const CTInput &in, CTOutput &out) {

    int64_t *ind = new int64_t[args.MAX_RAYLEN];
    float *wgt = new float[args.MAX_RAYLEN];
    int numb;

    get_tracing(alpha, detectorX, detectorY,
                ind, wgt, numb, args);

    if (max_numb<numb) max_numb = numb;
    ave_numb += numb;

    float Af = -in.sino_data(alpha, detectorX, detectorY);

    minIMAGE(Af, ind, wgt, numb, lambda, args,in,out);
    minEDGE(ind, wgt, numb, lambda, args,in,out);

    delete[] ind;
    delete[] wgt;
}

void compute(float lambda, int np, int ndx,
             const Parameter &args, const CTInput &in, CTOutput &out) {
    int64_t *ind = new int64_t[args.NDY*args.MAX_RAYLEN];
    float *wgt = new float[args.NDY*args.MAX_RAYLEN];
    int *numb = new int[args.NDY];

    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        get_tracing(np, ndx, ndy,
                    ind+ndy*args.MAX_RAYLEN, wgt+ndy*args.MAX_RAYLEN, numb[ndy],
                    args);
        if (max_numb<numb[ndy]) max_numb = numb[ndy];
        ave_numb += numb[ndy];
    }

    gdIMAGE(lambda, np, ndx,
            ind, wgt, numb,
            args, in, out);

    gdEDGE(lambda, ind, wgt, numb,
           args, in, out);

    delete [] ind;
    delete [] wgt;
    delete [] numb;
}

void wrapper_single(float lambda, const Parameter &args,const CTInput &in,CTOutput &out) {
    for (int k = 0; k < args.NPROJ; ++k)
        for (int i = 0; i < args.NDX; ++i)
            #pragma omp parallel for
            for (int j = 0; j < args.NDY; ++j) {
                compute_single(lambda, k, i, j, args,in,out);
                //cout << format("%1% %2% %3%") %k%i%j <<endl;
            }
}

void wrapper(float lambda, const Parameter &args,const CTInput &in,CTOutput &out) {
    for (int np = 0; np < args.NPROJ; np += 30) {
    //for (int np = 0; np < args.NPROJ; ++np) {
        //for (int slice = args.NZ-1; slice>=0; --slice) {
        //for (int ndx = 0; ndx < 1; ++ndx) {
        for (int ndx = 0; ndx < args.NDX; ++ndx) {
            //#pragma omp parallel for
            compute(lambda, np, ndx, args,in,out);
            //cout << format("%1% %2% %3%") %k%i%j <<endl;
        }
    }
}

void ct3d(const Parameter &args,const CTInput &in,CTOutput &out) {
    max_numb = 0;

    out.allocate();

    float lambda = 1;

    int64_t global_start = timer_s();

    for (int iters = 0; iters < args.ITERATIONS; ++iters) {
        //lambda = 1.0/(100.0+iters*2.0);

        cout << format("iter = %1%, lambda = %2%") % iters % lambda <<endl;

        ave_numb = 0;
        Afg = 0;

        int64_t start = timer_s();
        wrapper(lambda,args,in,out);
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
    //in.read_pretracing(args.PRETRACING_FILE);
    ct3d(args,in,out);
    out.write_img(args.OUTPUT_DIR);
    out.write_edge(args.OUTPUT_DIR);

    return 0;
}

void gdIMAGE(float lambda, int np, int ndx,
              int64_t *ind, float *wgt, int *numb,
              const Parameter &args, const CTInput &in, CTOutput &out) {
    float *r = new float[args.NDY];
    float *g = new float[args.NDY];
    float *Af = new float[args.NDY];

    img_type *d = new img_type[args.NDY*args.MAX_RAYLEN];
    float *q = new float[args.NDY];

    for (int ndy = 0; ndy < args.NDY; ++ndy)
        g[ndy] = in.sino_data(np,ndx,ndy);

    // Af = A(f)
    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        int64_t *__ind = ind + ndy*args.MAX_RAYLEN;
        float *__wgt = wgt + ndy*args.MAX_RAYLEN;
        int __numb = numb[ndy];

        float accum_Af = 0;
        for (int i = 0; i<__numb; ++i) {
            int64_t idx,nz,nx,ny;
            idx = __ind[i];
            nz = idx/(args.NX*args.NY);
            idx %= args.NX*args.NY;
            nx = idx/args.NY;
            ny = idx%args.NY;
            accum_Af += out.img_data(nz,nx,ny)*__wgt[i];
        }
        Af[ndy] = accum_Af;
    }

    // r = g-Af
    for (int i = 0; i<args.NDY; ++i) {
        r[i] = g[i]-Af[i];
        Afg += sqr(r[i]);
    }

    // d = A*(r)
    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        float *__wgt = wgt + ndy*args.MAX_RAYLEN;
        int __numb = numb[ndy];
        for (int i = 0; i<__numb; ++i) {
            d[ndy*args.MAX_RAYLEN+i] += r[ndy]*__wgt[i];
        }
    }

    // d += alpha*(v^2 \nabla f) + e^2 \laplace f
    for (int ndy = 0; ndy < args.NDY; ++ndy) {
        for (int i = 0; i<numb[ndy]; ++i) {
            int64_t idx = ind[ndy*args.MAX_RAYLEN+i],nz,nx,ny;
            nz = idx/(args.NX*args.NY);
            idx %= args.NX*args.NY;
            nx = idx/args.NY;
            ny = idx%args.NY;

            d[ndy*args.MAX_RAYLEN+i] += args.ALPHA*(
                +sqr(out.edge_data(nz  ,nx  ,ny  ))*(out.img_data(nz+1,nx  ,ny  ) - out.img_data(nz  ,nx  ,ny  ))
                -sqr(out.edge_data(nz-1,nx  ,ny  ))*(out.img_data(nz  ,nx  ,ny  ) - out.img_data(nz-1,nx  ,ny  ))
                +sqr(out.edge_data(nz  ,nx  ,ny  ))*(out.img_data(nz  ,nx+1,ny  ) - out.img_data(nz  ,nx  ,ny  ))
                -sqr(out.edge_data(nz  ,nx-1,ny  ))*(out.img_data(nz  ,nx  ,ny  ) - out.img_data(nz  ,nx-1,ny  ))
                +sqr(out.edge_data(nz  ,nx  ,ny  ))*(out.img_data(nz  ,nx  ,ny+1) - out.img_data(nz  ,nx  ,ny  ))
                -sqr(out.edge_data(nz  ,nx  ,ny-1))*(out.img_data(nz  ,nx  ,ny  ) - out.img_data(nz  ,nx  ,ny-1))
                +sqr(args.EPSILON)*(
                    +  out.img_data(nz+1,nx  ,ny  )
                    +  out.img_data(nz  ,nx+1,ny  )
                    +  out.img_data(nz  ,nx  ,ny+1)
                    -6*out.img_data(nz  ,nx  ,ny  )
                    +  out.img_data(nz-1,nx  ,ny  )
                    +  out.img_data(nz  ,nx-1,ny  )
                    +  out.img_data(nz  ,nx  ,ny-1)
                )
            );
        }
    }

    // q = A(d)
    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        float *__wgt = wgt + ndy*args.MAX_RAYLEN;
        int __numb = numb[ndy];

        float accum_q = 0;
        for (int i = 0; i<__numb; ++i) {
            accum_q += d[ndy*args.MAX_RAYLEN+i]*__wgt[i];
        }
        q[ndy] = accum_q;
    }

    float nV_KepsNablap = 0;
    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        int64_t *__ind = ind + ndy*args.MAX_RAYLEN;
        int __numb = numb[ndy];

        for (int i = 0; i<__numb; ++i) {
            int64_t idx,nz,nx,ny;
            idx = __ind[i];
            nz = idx/(args.NX*args.NY);
            idx %= args.NX*args.NY;
            nx = idx/args.NY;
            ny = idx%args.NY;

            nV_KepsNablap += (sqr(out.edge_data(nz,nx,ny))+sqr(args.EPSILON))*(
                +sqr(d[ndy*args.MAX_RAYLEN+i]-out.img_data(nz-1,nx  ,ny  ))
                +sqr(d[ndy*args.MAX_RAYLEN+i]-out.img_data(nz  ,nx-1,ny  ))
                +sqr(d[ndy*args.MAX_RAYLEN+i]-out.img_data(nz  ,nx  ,ny-1))
            );
        }
    }

    
    float pSkalard = 0, normAp = 0;
    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        int __numb = numb[ndy];
        for (int i = 0; i<__numb; ++i) {
            int idx = ndy*args.MAX_RAYLEN+i;
            pSkalard += d[idx]*d[idx];
        }
        normAp += sqr(q[ndy]);
    }

    if (fabs(normAp+args.ALPHA*nV_KepsNablap)<1e-4) return;

    float c_1 = pSkalard/(normAp+args.ALPHA*nV_KepsNablap);

    //cout << format("nV_KepsNablap = %1%, pSkalard = %2%, normAp = %3%, c_1 = %4%") % nV_KepsNablap % pSkalard % normAp % c_1 << endl;

    // update f
    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        int64_t *__ind = ind + ndy*args.MAX_RAYLEN;
        int __numb = numb[ndy];

        for (int i = 0; i<__numb; ++i) {
            int64_t idx,nz,nx,ny;
            idx = __ind[i];
            nz = idx/(args.NX*args.NY);
            idx %= args.NX*args.NY;
            nx = idx/args.NY;
            ny = idx%args.NY;

            out.img_data(nz,nx,ny) += lambda * c_1 * d[ndy*args.MAX_RAYLEN+i];
        }
    }

    delete [] r;
    delete [] g;
    delete [] Af;
    delete [] d;
    delete [] q;
}

void gdEDGE(float lambda,
             int64_t *ind, float *wgt, int *numb,
             const Parameter &args, const CTInput &in, CTOutput &out) {
    float *n = new float[args.NDY*args.MAX_RAYLEN];
    float *d = new float[args.NDY*args.MAX_RAYLEN];

    float nD = 0;

    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        int64_t *__ind = ind + ndy*args.MAX_RAYLEN;
        int __numb = numb[ndy];

        for (int i = 0; i<__numb; ++i) {
            int64_t idx,nz,nx,ny;
            idx = __ind[i];
            nz = idx/(args.NX*args.NY);
            idx %= args.NX*args.NY;
            nx = idx/args.NY;
            ny = idx%args.NY;

            n[ndy*args.MAX_RAYLEN+i] = 0.25*(
                +sqr(out.img_data(nz,nx,ny) - out.img_data(nz-1,nx  ,ny  ))
                +sqr(out.img_data(nz,nx,ny) - out.img_data(nz  ,nx-1,ny  ))
                +sqr(out.img_data(nz,nx,ny) - out.img_data(nz  ,nx  ,ny-1))
            );

            
            d[ndy*args.MAX_RAYLEN+i] = -(
                +args.ALPHA*out.edge_data(nz,nx,ny)*n[ndy*args.MAX_RAYLEN+i]
                +args.BETA/(4*args.EPSILON)*(out.edge_data(nz,nx,ny)-1)
                -args.BETA*args.EPSILON*0.25*(
                    +  out.edge_data(nz+1,nx  ,ny  )
                    +  out.edge_data(nz  ,nx+1,ny  )
                    +  out.edge_data(nz  ,nx  ,ny+1)
                    -6*out.edge_data(nz  ,nx  ,ny  )
                    +  out.edge_data(nz-1,nx  ,ny  )
                    +  out.edge_data(nz  ,nx-1,ny  )
                    +  out.edge_data(nz  ,nx  ,ny-1)
                )
            );

            nD += sqr(d[ndy*args.MAX_RAYLEN+i]);
        }
    }

    float nNablaD = 0;
    float tmp = 0;
    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        int64_t *__ind = ind + ndy*args.MAX_RAYLEN;
        int __numb = numb[ndy];

        for (int i = 0; i<__numb; ++i) {
            int64_t idx,nz,nx,ny;
            idx = __ind[i];
            nz = idx/(args.NX*args.NY);
            idx %= args.NX*args.NY;
            nx = idx/args.NY;
            ny = idx%args.NY;

            nNablaD += 0.25*(
                +sqr(out.edge_data(nz+1,nx  ,ny  ) - out.edge_data(nz-1,nx  ,ny  ))
                +sqr(out.edge_data(nz  ,nx+1,ny  ) - out.edge_data(nz  ,nx-1,ny  ))
                +sqr(out.edge_data(nz  ,nx  ,ny+1) - out.edge_data(nz  ,nx  ,ny-1))
            );

            tmp += n[ndy*args.MAX_RAYLEN+i]*sqr(d[ndy*args.MAX_RAYLEN+i]);
        }
    }

    if (fabs(args.ALPHA*tmp+args.BETA*args.EPSILON*nNablaD+(args.BETA/(4*args.EPSILON))*nD)<1e-4) return;

    float c = nD/(args.ALPHA*tmp+args.BETA*args.EPSILON*nNablaD+(args.BETA/(4*args.EPSILON))*nD);

    //cout << format("nD = %1%, nNablaD = %2%, c = %3%") % nD % nNablaD % c << endl;

    for (int ndy = 0; ndy<args.NDY; ++ndy) {
        int64_t *__ind = ind + ndy*args.MAX_RAYLEN;
        int __numb = numb[ndy];

        for (int i = 0; i<__numb; ++i) {
            int64_t idx,nz,nx,ny;
            idx = __ind[i];
            nz = idx/(args.NX*args.NY);
            idx %= args.NX*args.NY;
            nx = idx/args.NY;
            ny = idx%args.NY;

            edge_type tmp_v = out.edge_data(nz,nx,ny) + c * d[ndy*args.MAX_RAYLEN+i];
            if (tmp_v<0) tmp_v = 0;
            if (tmp_v>1) tmp_v = 1;
            out.edge_data(nz,nx,ny) = tmp_v;
        }
    }

    delete [] n;
    delete [] d;
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
