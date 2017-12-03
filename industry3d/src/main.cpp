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

    float tmp_Af = out.minIMAGE(Af, ind, wgt, numb, lambda);
    out.minEDGE(ind, wgt, numb, lambda);
    Afg += sqr(tmp_Af);

    delete[] ind;
    delete[] wgt;
}

void wrapper(float lambda,float *sin_table,float *cos_table,
             const Parameter &args,const CTInput &in,CTOutput &out) {
    for (int k = 0; k < args.NPROJ; ++k)
		for (int i = 0; i < args.NDX; ++i)
            //#pragma omp parallel for
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
    out.allocate_edge(args.NX*args.NY*args.NZ);

    float lambda = 0.01;

    int64_t global_start = timer_s();

    for (int iters = 0; iters < args.ITERATIONS; ++iters) {
        lambda = 1.0/(100.0+iters*2.0);

        cout << boost::format("iter = %1%, lambda = %2%") % iters % lambda <<endl;

        ave_numb = 0;
        Afg = 0;

        int64_t start = timer_s();
        wrapper(lambda,sin_table,cos_table, args,in,out);
        int64_t end = timer_s();
        
        ave_numb /= args.NPROJ*args.NDX*args.NDY;
        cout << boost::format("||Af-g||^2 = %1%") % Afg <<endl;
        cout << boost::format("time used = %1% seconds") % (end-start) <<endl;
        cout << "-------------------------------------------" <<endl;
    }

    int64_t global_end = timer_s();
    
    cout << "===========================================" <<endl;
    cout << boost::format("TOTAL time used = %1% seconds") % (global_end-global_start) <<endl;
    cout<< boost::format("actual max raylen = %1%, average raylen = %2%") % max_numb % ave_numb << endl;

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
    ct3d(args,in,out);
    out.write_img(args.OUTPUT_DIR);
    out.write_edge(args.OUTPUT_DIR);

    return 0;
}
