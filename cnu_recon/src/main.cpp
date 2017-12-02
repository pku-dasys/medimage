#include "ct3d.h"
#include "utility.h"
#include "tracing.h"

#include <iostream>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace std;

int max_numb;

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
void compute(float lambda,float *sin_table,float *cos_table, int alpha, int detectorX, int detectorY,
             const Parameter &args, const CTInput &in, CTOutput &out) {
    //const float sina = sin_table[alpha], cosa = cos_table[alpha];
    float srcX,srcY,srcZ;
    float dstX,dstY,dstZ;

    if (args.BEAM=="Parallel") {
        srcX = detectorY+0.5 - args.HALFSIZE;
        srcX *= args.PIXELSIZE;
        srcY = -args.SOD;

        srcZ = args.NX - detectorX - 1;
        srcZ *= args.PIXELSIZE;

        
        dstX = detectorY+0.5 - args.HALFSIZE;
        dstX *= args.PIXELSIZE;
        dstY = args.SDD-args.SOD;

        dstZ = args.NX-detectorX-1;
        dstZ *= args.PIXELSIZE;

        float theta = (float)alpha/args.NPROJ*2*PI;

        // rotation
        //rotate_axis_3d(srcX, srcY, srcZ, 0, 1, 0, theta);
        //rotate_axis_3d(dstX, dstY, dstZ, 0, 1, 0, theta);
        
        rotate_2d(srcX, srcY, theta);
        rotate_2d(dstX, dstY, theta);

	    srcX /= args.SAMPLESIZE;
	    srcZ /= args.SAMPLESIZE;
	    srcY /= args.SAMPLESIZE;
	    dstX /= args.SAMPLESIZE;
	    dstY /= args.SAMPLESIZE;
	    dstZ /= args.SAMPLESIZE;

        srcX += args.HALFSIZE;
        srcY += args.HALFSIZE;
        dstX += args.HALFSIZE;
        dstY += args.HALFSIZE;
    }
    else if (args.BEAM=="Cone") {
//        float oridstX,oridstY;
        srcZ = 0.0;
        srcX = 0.0;
        srcY = -args.SOD;
        float theta = (float)alpha/args.NPROJ*2*PI;
        rotate_2d(srcX, srcY, theta);

        dstX = detectorX - args.NDX*0.5 + 0.5;
        dstX *= args.PIXELSIZE;
        dstY = (args.SDD - args.SOD);
        dstZ = detectorY - args.NDY*0.5 + 0.5;
        dstZ *= args.PIXELSIZE;
        rotate_2d(dstX, dstY, theta);

	    srcX /= args.SAMPLESIZE;
	    srcZ /= args.SAMPLESIZE;
	    srcY /= args.SAMPLESIZE;
	    dstX /= args.SAMPLESIZE;
	    dstY /= args.SAMPLESIZE;
	    dstZ /= args.SAMPLESIZE;

        srcX += args.NX/2.0;
        srcY += args.NY/2.0;
        srcZ += args.NZ/2.0;
        dstX += args.NX/2.0;
        dstY += args.NY/2.0;
        dstZ += args.NZ/2.0;
    }

    //cout << boost::format("%1%  :  %2% %3% %4%     %5% %6% %7%") %alpha%srcX%srcY%srcZ%dstX%dstY%dstZ <<endl;

    int64_t *ind = new int64_t[args.MAX_RAYLEN];
    float *wgt = new float[args.MAX_RAYLEN];
    int numb;

    float Af = -in.sino_data(alpha, detectorX, detectorY);

    forward_proj(args,
                 srcX, srcY, srcZ,
                 dstX, dstY, dstZ,
                 ind, wgt, numb);

    // if (alpha==270 && detectorX==0 && detectorY==0) {
    //     printf("%d %d %d\n",alpha,detectorX, detectorY);
    //     cout << boost::format("%1%  :  %2% %3% %4%     %5% %6% %7%") %alpha%srcX%srcY%srcZ%dstX%dstY%dstZ <<endl;
    //     printf("%d\n",numb);
    //     for (int i = 0; i<numb; ++i) {
    //         int z,x,y;
    //         int64_t tmp = ind[i];
    //         z = tmp/(args.NX*args.NY);
    //         tmp %= args.NX*args.NY;
    //         x = tmp/args.NY;
    //         y = tmp%args.NY;
    //         cout<<ind[i]<<endl;
    //         printf("%lld,%f  -  %d %d %d\n",tmp,wgt[i],z, x, y);
    //     }
    //     printf("\n");
    //     fflush(stdout);
    //     exit(0);
    // }
    
    if (max_numb<numb) max_numb = numb;
    
    out.minIMAGE(Af, ind, wgt, numb, lambda);

    delete[] ind;
    delete[] wgt;
}

void wrapper(float lambda,float *sin_table,float *cos_table,
             const Parameter &args,const CTInput &in,CTOutput &out) {
    for (int k = 0; k < args.NPROJ; ++k)
    //for (int k = 270; k < 271; ++k)
		for (int i = 0; i < args.NDX; ++i)
			for (int j = 0; j < args.NDY; ++j) {
				compute(lambda, sin_table, cos_table, k, i, j, args,in,out);
                //cout << boost::format("%1% %2% %3%") %k%i%j <<endl;
            }
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

    out.allocate_img(args.NX*args.NY*args.NZ);

    float lambda = 0.01;
    for (int iters = 0; iters < args.ITERATIONS; ++iters) {
        lambda = 1.0/(100.0+iters*2.0);

        cout << boost::format("iter = %1%, lambda = %2%") % iters % lambda <<endl;

        wrapper(lambda,sin_table,cos_table, args,in,out);
    }

    cout<< max_numb << endl;

    delete [] sin_table;
    delete [] cos_table;
}

int main(int argc, char** argv) {
    Parameter args;

    args.parse_config(argc, argv);
    args.print_options();

    CTInput in = CTInput(args);
    CTOutput out = CTOutput(args);

    in.read_sino(args.RAW_DATA_FILE);
    ct3d(args,in,out);
    out.write_img(args.OUTPUT_DIR);

    return 0;
}
