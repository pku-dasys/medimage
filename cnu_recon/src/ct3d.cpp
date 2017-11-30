#include "utility.h"
#include "main.h"
#include "ct3d.h"
#include "tracing.h"
#include <cstdio>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <cstring>
#include <iostream>
using namespace std;
using namespace boost;

int max_numb;


// alpha is the rotate angle of light source
// (or the object, then result will be a mirror image),
// from the view of +Z to -Z, counterclockwise(+X->+Y->-X->-Y)
// is positive. Initial angle is +Y
// the parameters should be int, as it represents the relative
// coordinate of grid in the detector
// here, detectorY actually represent Z-axis
void compute(float lambda,float *sin_table,float *cos_table, int alpha, int detectorX, int detectorY,
             const Parameter &args, const CTInput &in, CTOutput &out) {
    const float sina = sin_table[alpha], cosa = cos_table[alpha];
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    if (args.BEAM=="Parallel") {
        srcX = detectorX;
        srcY = detectorY;
        srcZ = 0.0;

        
        dstX = detectorX;
        dstY = detectorY;
        dstZ = args.SDD;

        //printf("%f %f %f     %f %f %f\n",srcX, srcY, srcZ, dstX, dstY, dstZ);
    }
    else if (args.BEAM=="Cone") {
        float oridstX,oridstY;
        srcZ = 0.0;

        srcX = (-args.SOD * sina);
        srcY = ( args.SOD * cosa);
        oridstX = detectorX - args.NDX*0.5 + 0.5;
        oridstY = (args.SOD - args.SDD);
        dstX = oridstX * cosa - oridstY * sina;
        dstY = oridstX * sina + oridstY * cosa;
        dstZ = detectorY - args.NDY*0.5 + 0.5;

        srcX += args.NX/2.0;
        srcY += args.NY/2.0;
        //srcZ += NZ/2.0;
        dstX += args.NX/2.0;
        dstY += args.NY/2.0;
        //dstZ += NZ/2.0;
    }

    int64 *ind = new int64[args.MAX_RAYLEN];
    float *wgt = new float[args.MAX_RAYLEN];
    int numb;

    float Af = -in.sino_data(alpha, detectorX, detectorY);

    forward_proj(args,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);

    // if (numb>0) {
    //     printf("%d %d %d\n",alpha,detectorX, detectorY);
    //     printf("%d\n",numb);
    //     for (int i = 0; i<numb; ++i) {
    //         int z,x,y;
    //         int64 tmp = ind[i];
    //         z = tmp/(NX*NY);
    //         tmp %= NX*NY;
    //         x = tmp/NY;
    //         y = tmp%NY;
    //         cout<<ind[i]<<endl;
    //         printf("%lld,%f  -  %d %d %d\n",tmp,wgt[i],z, x, y);
    //     }
    //     printf("\n");
    //     fflush(stdout);
    //     exit(0);
    // }
    
    if (max_numb<numb) max_numb = max_numb;
    
    out.minIMAGE(Af, ind, wgt, numb, lambda);

    delete[] ind;
    delete[] wgt;
}

void wrapper(float lambda,float *sin_table,float *cos_table,
             const Parameter &args,const CTInput &in,CTOutput &out) {
	for (int k = 0; k < args.NPROJ; ++k)
		for (int i = 0; i < args.NDX; ++i)
			for (int j = 0; j < args.NDY; ++j)
				compute(lambda, sin_table, cos_table, k, i, j, args,in,out);
}

void ct3d(const Parameter &args,const CTInput &in,CTOutput &out) {
    float *cos_table = new float[args.NPROJ];
    float *sin_table = new float[args.NPROJ];

    for(int i = 0; i < args.NPROJ; i++) {
        const float alpha = (float)i/args.NPROJ*2*PI;
        sin_table[i] = sin(alpha);
        cos_table[i] = cos(alpha);
    }

    max_numb = 0;

    float lambda = 0.001;
    for (int iters = 0; iters < args.ITERATIONS; ++iters) {
        lambda = 1.0/(1000.0+iters*50.0);

        printf("iter = %d, lambda = %lf\n", iters, lambda);
        fflush(stdout);

        wrapper(lambda,sin_table,cos_table, args,in,out);
    }
    
    printf("max_numb = %d\n",max_numb);

    delete [] sin_table;
    delete [] cos_table;
}
